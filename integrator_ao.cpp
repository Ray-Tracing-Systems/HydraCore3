#include "integrator_pt.h"

// Adapted from HeavyRT/examples/02_rtao. Camera/geometry data and floating-point
// output use Hydra's layout; AO rays retain cosine sampling and closest-hit queries.
void Integrator::AmbientOcclusion(uint tid, uint channels, float* out_color)
{
  for(uint iter = 0; iter < RTVPersistent_Iters(); ++iter)
  {
    RTVPersistent_SetIter(iter);
    float4 rayPos, rayDir, hitPos, hitNorm;
    float visibility;
    kernel_InitAOEyeRay(tid, &rayPos, &rayDir, &visibility);
    kernel_TraceAOEyeRay(tid, &rayPos, &rayDir, &hitPos, &hitNorm);
    // Persistent lanes partition the Sobol samples of a pixel. One invocation
    // computes the requested total spp, including the persistent-rays variant.
    for(uint sample = 0; sample < m_aoSamples / RTVPersistent_Iters(); ++sample)
    {
      kernel_InitAORay(tid, sample, &hitPos, &hitNorm, &rayPos, &rayDir);
      kernel_TraceAORay(tid, &rayPos, &rayDir, &visibility);
    }
    kernel_ContributeAO(tid, channels, &visibility, out_color);
  }
}

void Integrator::kernel_InitAOEyeRay(uint tid, float4* rayPos, float4* rayDir, float* visibility)
{
  const uint xy = m_packedXY[RTVPersistent_ThreadId(tid)];
  const float x = float(xy & 0xFFFF) + 0.5f;
  const float y = float(xy >> 16) + 0.5f;
  float3 origin = float3(0.0f);
  float3 direction = EyeRayDirNormalized((x + float(m_winStartX)) / float(m_fbWidth),
                                        (y + float(m_winStartY)) / float(m_fbHeight), m_projInv);
  transform_ray3f(m_worldViewInv, &origin, &direction);
  *rayPos = to_float4(origin, 0.0f);
  *rayDir = to_float4(direction, FLT_MAX);
  *visibility = 0.0f;
}

void Integrator::kernel_TraceAOEyeRay(uint tid, const float4* rayPos, const float4* rayDir,
                                     float4* hitPos, float4* hitNorm)
{
  const CRT_Hit hit = m_pAccelStruct->RayQuery_NearestHit(*rayPos, *rayDir);
  *hitPos = float4(0.0f);
  *hitNorm = float4(0.0f, 1.0f, 0.0f, 0.0f);
  if(hit.geomId == uint(-1))
    return;
  const uint2 offsets = m_matVertOffset[hit.geomId];
  const uint a = m_triIndices[(offsets.x + hit.primId) * 3 + 0] + offsets.y;
  const uint b = m_triIndices[(offsets.x + hit.primId) * 3 + 1] + offsets.y;
  const uint c = m_triIndices[(offsets.x + hit.primId) * 3 + 2] + offsets.y;
  const float u = hit.coords[0];
  const float v = hit.coords[1];
  float3 normal = to_float3((1.0f - u - v) * m_vData8f[a].normAndTx +
                            v * m_vData8f[b].normAndTx + u * m_vData8f[c].normAndTx);
  normal = normalize(mul3x3(m_normMatrices[hit.instId], normal));
  if(dot(normal, to_float3(*rayDir)) > 0.0f)
    normal = -normal;
  *hitPos = to_float4(to_float3(*rayPos) + (0.9999995f * hit.t) * to_float3(*rayDir), 1.0f);
  *hitNorm = to_float4(normal, 0.0f);
}

void Integrator::kernel_InitAORay(uint tid, uint sample, const float4* hitPos, const float4* hitNorm,
                                  float4* rayPos, float4* rayDir)
{
  *rayPos = float4(0.0f, 0.0f, 0.0f, -1.0f);
  *rayDir = float4(0.0f);
  if(hitPos->w == 0.0f)
    return;
  const uint xy = m_packedXY[RTVPersistent_ThreadId(tid)];
  const uint x = (xy & 0xFFFF) % AO_TILE_SIZE;
  const uint y = (xy >> 16) % AO_TILE_SIZE;
  const uint sampleId = sample * RTVPersistent_Iters() + tid % RTVPersistent_Iters();
  const float2 uv = m_aoRandomsTile[(y * AO_TILE_SIZE + x) * m_aoSamples + sampleId];
  const float3 normal = to_float3(*hitNorm);
  const float3 direction = MapSampleToCosineDistribution(uv.x, uv.y, normal, normal, 1.0f);
  *rayPos = to_float4(OffsRayPos(to_float3(*hitPos), normal, direction), 0.0f);
  *rayDir = to_float4(direction, m_aoRadius);
}

void Integrator::kernel_TraceAORay(uint tid, const float4* rayPos, const float4* rayDir, float* visibility)
{
  if(rayPos->w < 0.0f)
    return;
  const CRT_Hit hit = m_pAccelStruct->RayQuery_NearestHit(*rayPos, *rayDir);
  *visibility += (hit.geomId == uint(-1)) ? 1.0f : 0.0f;
}

void Integrator::kernel_ContributeAO(uint tid, uint channels, const float* visibility, float* out_color)
{
  const float value = *visibility / float(m_aoSamples);
  const float4 color = RTVPersistent_ReduceAdd4f(float4(value, value, value, 0.0f));
  if(RTVPersistent_IsFirst())
  {
    const uint xy = m_packedXY[RTVPersistent_ThreadId(tid)];
    const uint pixel = (xy >> 16) * m_winWidth + (xy & 0xFFFF);
    out_color[pixel * channels] = color.x;
    if(channels == 4)
    {
      out_color[pixel * channels + 1] = color.y;
      out_color[pixel * channels + 2] = color.z;
      out_color[pixel * channels + 3] = 1.0f;
    }
  }
}
