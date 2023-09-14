#include "integrator_pt.h"
#include "include/crandom.h"

#include <chrono>
#include <string>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

void Integrator::kernel_PackXY(uint tidX, uint tidY, uint* out_pakedXY)
{
  const uint inBlockIdX = tidX % 8; // 8x8 blocks
  const uint inBlockIdY = tidY % 8; // 8x8 blocks
 
  const uint localIndex = inBlockIdY*8 + inBlockIdX;
  const uint wBlocks    = m_winWidth/8;

  const uint blockX     = tidX/8;
  const uint blockY     = tidY/8;
  const uint offset     = (blockX + blockY*wBlocks)*8*8 + localIndex;

  out_pakedXY[offset] = ((tidY << 16) & 0xFFFF0000) | (tidX & 0x0000FFFF);
}

void Integrator::kernel_InitEyeRay(uint tid, const uint* packedXY, float4* rayPosAndNear, float4* rayDirAndFar) // (tid,tidX,tidY,tidZ) are SPECIAL PREDEFINED NAMES!!!
{
  const uint XY = packedXY[tid];

  const uint x = (XY & 0x0000FFFF);
  const uint y = (XY & 0xFFFF0000) >> 16;

  float3 rayDir = EyeRayDirNormalized((float(x)+0.5f)/float(m_winWidth), (float(y)+0.5f)/float(m_winHeight), m_projInv);
  float3 rayPos = float3(0,0,0);

  transform_ray3f(m_worldViewInv, 
                  &rayPos, &rayDir);
  
  *rayPosAndNear = to_float4(rayPos, 0.0f);
  *rayDirAndFar  = to_float4(rayDir, FLT_MAX);
}

void Integrator::kernel_InitEyeRay3(uint tid, const uint* packedXY, 
                                   float4* rayPosAndNear, float4* rayDirAndFar,
                                   float4* accumColor,    float4* accumuThoroughput,
                                   uint* rayFlags) // 
{
  *accumColor        = make_float4(0,0,0,1);
  *accumuThoroughput = make_float4(1,1,1,1);
  //RandomGen genLocal = m_randomGens[tid];
  *rayFlags          = 0;

  const uint XY = packedXY[tid];

  const uint x = (XY & 0x0000FFFF);
  const uint y = (XY & 0xFFFF0000) >> 16;

  float3 rayDir = EyeRayDirNormalized((float(x))/float(m_winWidth), 
                                      (float(y))/float(m_winHeight), m_projInv);
  float3 rayPos = float3(0,0,0);

  transform_ray3f(m_worldViewInv, &rayPos, &rayDir);
  
  *rayPosAndNear = to_float4(rayPos, 0.0f);
  *rayDirAndFar  = to_float4(rayDir, FLT_MAX);
}


bool Integrator::kernel_RayTrace(uint tid, const float4* rayPosAndNear, float4* rayDirAndFar,
                                 Lite_Hit* out_hit, float2* out_bars)
{
  const float4 rayPos = *rayPosAndNear;
  const float4 rayDir = *rayDirAndFar ;

  CRT_Hit hit = m_pAccelStruct->RayQuery_NearestHit(rayPos, rayDir);
  
  Lite_Hit res;
  res.primId = hit.primId;
  res.instId = hit.instId;
  res.geomId = hit.geomId;
  res.t      = hit.t;

  float2 baricentrics = float2(hit.coords[0], hit.coords[1]);
 
  *out_hit  = res;
  *out_bars = baricentrics;
  return (res.primId != -1);
}


void Integrator::kernel_RealColorToUint32(uint tid, float4* a_accumColor, uint* out_color)
{
  out_color[tid] = RealColorToUint32(*a_accumColor);
}

void Integrator::kernel_GetRayColor(uint tid, const Lite_Hit* in_hit, const uint* in_pakedXY, uint* out_color)
{ 
  const Lite_Hit lhit = *in_hit;
  if(lhit.geomId == -1)
  {
    out_color[tid] = 0;
    return;
  }

  const uint32_t matId = m_matIdByPrimId[m_matIdOffsets[lhit.geomId] + lhit.primId];
  const float4 mdata   = m_materials[matId].baseColor;
  const float3 color   = mdata.w > 0.0f ? clamp(float3(mdata.w,mdata.w,mdata.w), 0.0f, 1.0f) : to_float3(mdata);

  const uint XY = in_pakedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;

  out_color[y*m_winWidth+x] = RealColorToUint32_f3(color); 
}


float3 Integrator::MaterialEvalWhitted(int a_materialId, float3 l, float3 v, float3 n, float2 tc)
{
  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const float3 texColor  = to_float3(m_textures[ m_materials[a_materialId].texId[0] ]->sample(texCoordT));
  const float3 color     = to_float3(m_materials[a_materialId].baseColor)*texColor;
  return lambertEvalBSDF(l, v, n)*color;
}

BsdfSample Integrator::MaterialSampleWhitted(int a_materialId, float3 v, float3 n, float2 tc)
{ 
  //uint  type             = m_materials[a_materialId].mtype;
  const float3 specular  = to_float3(m_materials[a_materialId].metalColor);
  const float3 coat      = to_float3(m_materials[a_materialId].coatColor);
  //const float  roughness = 1.0f - m_materials[a_materialId].glosiness;
  float alpha            = m_materials[a_materialId].alpha;
  
  const float3 pefReflDir = reflect((-1.0f)*v, n);
  const float3 reflColor  = alpha*specular + (1.0f - alpha)*coat;

  //if(a_materialId == 4)
  //{
  //  int a = 2;
  //}

  BsdfSample res;
  res.dir = pefReflDir;
  res.val     = reflColor;
  res.pdf       = 1.0f;
  res.flags     = RAY_EVENT_S;
  return res;
}


void Integrator::kernel_RayBounce(uint tid, uint bounce, const float4* in_hitPart1, const float4* in_hitPart2,
                                  float4* rayPosAndNear, float4* rayDirAndFar, float4* accumColor, float4* accumThoroughput,
                                  uint* rayFlags)
{
  const uint currRayFlags = *rayFlags;
  if(isDeadRay(currRayFlags))
    return;

  const uint32_t matId = extractMatId(currRayFlags);

  // process surface hit case
  //
  const float3 ray_dir = to_float3(*rayDirAndFar);
  //const float3 ray_pos = to_float3(*rayPosAndNear);

  const float4 data1 = *in_hitPart1;
  const float4 data2 = *in_hitPart2;

  SurfaceHit hit;
  hit.pos  = to_float3(data1);
  hit.norm = to_float3(data2);
  hit.uv   = float2(data1.w, data2.w);

  // process light hit case
  //
  if(m_materials[matId].mtype == MAT_TYPE_LIGHT_SOURCE)
  {
    const float2 texCoordT = mulRows2x4(m_materials[matId].row0[0], m_materials[matId].row1[0], hit.uv);
    const float3 texColor  = to_float3(m_textures[ m_materials[matId].texId[0] ]->sample(texCoordT));

    const float3 lightIntensity = to_float3(m_materials[matId].baseColor)*texColor;
    const uint lightId          = m_materials[matId].lightId;
    float lightDirectionAtten   = (lightId == 0xFFFFFFFF) ? 1.0f : dot(to_float3(*rayDirAndFar), float3(0,-1,0)) < 0.0f ? 1.0f : 0.0f; // TODO: read light info, gety light direction and e.t.c;


    float4 currAccumColor      = *accumColor;
    float4 currAccumThroughput = *accumThoroughput;

    currAccumColor.x += currAccumThroughput.x * lightIntensity.x * lightDirectionAtten;
    currAccumColor.y += currAccumThroughput.y * lightIntensity.y * lightDirectionAtten;
    currAccumColor.z += currAccumThroughput.z * lightIntensity.z * lightDirectionAtten;

    *accumColor = currAccumColor;
    *rayFlags   = currRayFlags | (RAY_FLAG_IS_DEAD | RAY_FLAG_HIT_LIGHT);
    return;
  }

  float4 shadeColor = float4(0.0f, 0.0f, 0.0f, 1.0f);
  for(uint lightId = 0; lightId < m_lights.size(); ++lightId)
  {
    const float3 lightPos = to_float3(m_lights[lightId].pos);
    const float hitDist   = sqrt(dot(hit.pos - lightPos, hit.pos - lightPos));

    const float3 shadowRayDir = normalize(lightPos - hit.pos);
    const float3 shadowRayPos = hit.pos + hit.norm * std::max(maxcomp(hit.pos), 1.0f) * 5e-6f; // TODO: see Ray Tracing Gems, also use flatNormal for offset

    const bool inShadow = m_pAccelStruct->RayQuery_AnyHit(to_float4(shadowRayPos, 0.0f), to_float4(shadowRayDir, hitDist * 0.9995f));

    if(!inShadow && dot(shadowRayDir, to_float3(m_lights[lightId].norm)) < 0.0f)
    {
      const float3 matSamColor = MaterialEvalWhitted(matId, shadowRayDir, (-1.0f)*ray_dir, hit.norm, hit.uv);
      const float cosThetaOut  = std::max(dot(shadowRayDir, hit.norm), 0.0f);
      shadeColor += to_float4(to_float3(m_lights[lightId].intensity) * matSamColor*cosThetaOut / (hitDist * hitDist), 0.0f);
    }
  }

  const BsdfSample matSam = MaterialSampleWhitted(matId, (-1.0f)*ray_dir, hit.norm, hit.uv);
  const float3 bxdfVal    = matSam.val;
  const float  cosTheta   = dot(matSam.dir, hit.norm);

  const float4 currThoroughput = *accumThoroughput;
  float4 currAccumColor        = *accumColor;

  currAccumColor.x += currThoroughput.x * shadeColor.x;
  currAccumColor.y += currThoroughput.y * shadeColor.y;
  currAccumColor.z += currThoroughput.z * shadeColor.z;

  *accumColor       = currAccumColor;
  *accumThoroughput = currThoroughput * cosTheta * to_float4(bxdfVal, 0.0f);

  *rayPosAndNear = to_float4(OffsRayPos(hit.pos, hit.norm, matSam.dir), 0.0f);
  *rayDirAndFar  = to_float4(matSam.dir, FLT_MAX);
  *rayFlags      = currRayFlags | matSam.flags;
}

void Integrator::kernel_ContributeToImage3(uint tid, const float4* a_accumColor, const uint* in_pakedXY, float4* out_color)
{
  const uint XY = in_pakedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;

  float4 color = *a_accumColor;
  out_color[y*m_winWidth+x] += color;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Integrator::PackXY(uint tidX, uint tidY)
{
  kernel_PackXY(tidX, tidY, m_packedXY.data());
}

void Integrator::CastSingleRay(uint tid, uint* out_color)
{
  float4 rayPosAndNear, rayDirAndFar;
  kernel_InitEyeRay(tid, m_packedXY.data(), &rayPosAndNear, &rayDirAndFar);

  Lite_Hit hit; 
  float2   baricentrics; 
  if(!kernel_RayTrace(tid, &rayPosAndNear, &rayDirAndFar, &hit, &baricentrics))
    return;
  
  kernel_GetRayColor(tid, &hit, m_packedXY.data(), out_color);
}

void Integrator::RayTrace(uint tid, float4* out_color)
{
  float4 accumColor, accumThroughput;
  float4 rayPosAndNear, rayDirAndFar;
  uint      rayFlags = 0;
  kernel_InitEyeRay3(tid, m_packedXY.data(), 
                     &rayPosAndNear, &rayDirAndFar, &accumColor, &accumThroughput, &rayFlags);

  for(uint depth = 0; depth < m_traceDepth; depth++)
  {
    float4 hitPart1, hitPart2;
    uint instId;
    kernel_RayTrace2(tid, &rayPosAndNear, &rayDirAndFar, &hitPart1, &hitPart2, &instId, &rayFlags);
    if(isDeadRay(rayFlags))
      break;

    kernel_RayBounce(tid, depth, &hitPart1, &hitPart2,
                     &rayPosAndNear, &rayDirAndFar, &accumColor, &accumThroughput, &rayFlags);

    if(isDeadRay(rayFlags))
      break;
  }

//  kernel_HitEnvironment(tid, &rayFlags, &rayDirAndFar, &mis, &accumThroughput,
//                        &accumColor);

  kernel_ContributeToImage3(tid, &accumColor, m_packedXY.data(),
                            out_color);
}