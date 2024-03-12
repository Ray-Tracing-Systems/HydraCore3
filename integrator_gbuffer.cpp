#include "integrator_pt.h"
#include "include/crandom.h"

#include <chrono>
#include <string>

void PlaneHammersley(float *result, int n)
{
  for (int k = 0; k<n; k++)
  {
    float u = 0;
    int kk = k;

    for (float p = 0.5f; kk; p *= 0.5f, kk >>= 1)
      if (kk & 1)                           // kk mod 2 == 1
        u += p;

    float v = (k + 0.5f) / n;

    result[2 * k + 0] = u;
    result[2 * k + 1] = v;
  }
}


static inline float projectedPixelSize(float dist, float FOV, float w, float h)
{
  float ppx = (FOV / w)*dist;
  float ppy = (FOV / h)*dist;

  if (dist > 0.0f)
    return 2.0f*fmax(ppx, ppy);
  else
    return 1000.0f;
}

static inline float surfaceSimilarity(float4 data1, float4 data2, const float MADXDIFF)
{
  const float MANXDIFF = 0.15f;

  float3 n1 = to_float3(data1);
  float3 n2 = to_float3(data2);

  float dist = length(n1 - n2);
  if (dist >= MANXDIFF)
    return 0.0f;

  float d1 = data1.w;
  float d2 = data2.w;

  if (std::abs(d1 - d2) >= MADXDIFF)
    return 0.0f;

  float normalSimilar = std::sqrt(1.0f - (dist / MANXDIFF));
  float depthSimilar  = std::sqrt(1.0f - std::abs(d1 - d2) / MADXDIFF);

  return normalSimilar * depthSimilar;
}

//static inline float gbuffDiff(GBufferAll s1, GBufferAll s2, const float a_fov, float w, float h)
//{
//  const float ppSize         = projectedPixelSize(s1.data1.depth, a_fov, w, h);
//  const float surfaceSimilar = surfaceSimilarity(to_float4(s1.data1.norm, s1.data1.depth),
//                                                 to_float4(s2.data1.norm, s2.data1.depth), ppSize*2.0f);
//
//  const float surfaceDiff    = 1.0f - surfaceSimilar;
//  const float objDiff        = (s1.data2.instId == s2.data2.instId && s1.data2.objId == s2.data2.objId) ? 0.0f : 1.0f;
//  const float matDiff        = (s1.data1.matId  == s2.data1.matId) ? 0.0f : 1.0f;
//  const float alphaDiff      = fabs(s1.data1.rgba.w - s2.data1.rgba.w);
//
//  return surfaceDiff + objDiff + matDiff + alphaDiff;
//}

//static inline float gbuffDiffObj(GBufferAll s1, GBufferAll s2, const float a_fov, int w, int h)
//{
//  const float objDiff = (s1.data2.instId == s2.data2.instId && s1.data2.objId == s2.data2.objId) ? 0.0f : 1.0f;
//  const float matDiff = (s1.data1.matId  == s2.data1.matId) ? 0.0f : 1.0f;
//
//  return objDiff + matDiff;
//}

static constexpr uint GBUFFER_SAMPLES = 16;

void Integrator::EvalGBuffer(uint tidX, uint tidY, GBufferPixel* out_gbuffer)
{
  //float4 rayPosAndNear, rayDirAndFar;
  //kernel_InitEyeRay(tidX, tidY, m_packedXY.data(), &rayPosAndNear, &rayDirAndFar);

  //Lite_Hit hit; 
  //float2   baricentrics; 
  //kernel_RayTrace(tidX, tidY, &rayPosAndNear, &rayDirAndFar, &hit, &baricentrics);
  
  //kernel_GetRayGBuff(tidX, tidY, &hit, &baricentrics, m_packedXY.data(), out_gbuffer);
}

void Integrator::EvalGBufferReduction(uint tidX, uint tidY, GBufferPixel* out_gbuffer)
{
  float minDiff   = 100000000.0f; 
  int   minDiffId = 0;

  //for (int i = 0; i < tidY; i++)
  //{
  //  float diff     = 0.0f;
  //  float coverage = 0.0f;
  //  for (int j = 0; j < tidY; j++)
  //  {
  //    const float thisDiff = gbuffDiff(out_gbuffer[i], out_gbuffer[j], fov, float(m_width), float(m_height));
  //    diff += thisDiff;
  //    if (thisDiff < 1.0f)
  //      coverage += 1.0f;
  //  }
  //
  //  coverage *= (1.0f / (float)tidY);
  //  samples[i].data1.coverage = coverage;
  //
  //  if (diff < minDiff)
  //  {
  //    minDiff   = diff;
  //    minDiffId = i;
  //  }
  //}
  
  out_gbuffer[0] = out_gbuffer[minDiffId];
}

void Integrator::EvalGBufferBlock(uint tidX, GBufferPixel* out_gbuffer)
{
  m_qmcHammersley.resize(GBUFFER_SAMPLES);                 // TODO: move to class constructor sms like that
  PlaneHammersley(&m_qmcHammersley[0].x, GBUFFER_SAMPLES); //

  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for(int i = 0; i < tidX; ++i) 
  {
    GBufferPixel blockData[GBUFFER_SAMPLES]; 

    for(int j = 0; j < GBUFFER_SAMPLES; ++j)
      EvalGBuffer(i,j,blockData);
    
    EvalGBufferReduction(i, GBUFFER_SAMPLES, out_gbuffer);
  }
}
