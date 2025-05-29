#include "integrator_pt.h"
#include "include/crandom.h"
#include "fourier/fspec.h"

#include <chrono>
#include <string>


#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;


void Integrator::PathTraceF(uint tid, uint channels, float* out_color)
{
  FourierSpec accumColor, accumThroughput;
  float4 rayPosAndNear, rayDirAndFar;
  RandomGen gen; 
  MisData   mis;
  uint      rayFlags;
  float     time;
  kernel_InitEyeRay2F(tid, &rayPosAndNear, &rayDirAndFar, &accumColor, &accumThroughput, &gen, &rayFlags, &mis, &time);

  for(uint depth = 0; depth < m_traceDepth; depth++) 
  {
    FourierSpec shadeColor;
    float4 hitPart1, hitPart2, hitPart3;
    uint instId;
    kernel_RayTrace2(tid, depth, &rayPosAndNear, &rayDirAndFar, &time, 
                     &hitPart1, &hitPart2, &hitPart3, &instId, &rayFlags);
    if(isDeadRay(rayFlags))
      break;
    
    kernel_SampleLightSourceF(tid, &rayPosAndNear, &rayDirAndFar, &hitPart1, &hitPart2, &hitPart3, &rayFlags, &time,
                             depth, &gen, &shadeColor);

    kernel_NextBounce(tid, depth, &hitPart1, &hitPart2, &hitPart3, &instId, &shadeColor,
                      &rayPosAndNear, &rayDirAndFar, &wavelengths, &accumColor, &accumThroughput, &gen, &mis, &rayFlags);

    if(isDeadRay(rayFlags))
      break;
  }

  kernel_HitEnvironment(tid, &rayFlags, &rayDirAndFar, &mis, &accumThroughput,
                        &accumColor);

  kernel_ContributeToImage(tid, &rayFlags, channels, &accumColor, &gen, m_packedXY.data(), &wavelengths, out_color);
}



void kernel_InitEyeRay2F(uint tid, float4* rayPosAndNear, float4* rayDirAndFar,
                         FourierSpec* accumColor, FourierSpec* accumuThoroughput, RandomGen* gen,
                         uint* rayFlags, MisData* misData, float* time)
{
  if(tid >= m_maxThreadId)
    return;

  *accumColor        = FourierSpec();
  *accumuThoroughput = FourierSpec(1.0f);
  RandomGen genLocal = m_randomGens[RandomGenId(tid)];
  *rayFlags          = 0;
  *misData           = makeInitialMisData();

  EyeRayData r = SampleCameraRay(&genLocal, tid);

  *time = r.timeSam;
 
  transform_ray3f(m_worldViewInv, &r.rayPos, &r.rayDir);

  *rayPosAndNear = to_float4(r.rayPos, 0.0f);
  *rayDirAndFar  = to_float4(r.rayDir, FLT_MAX);
  *gen           = genLocal;
}


void Integrator::kernel_SampleLightSourceF(uint tid, const float4* rayPosAndNear, const float4* rayDirAndFar,
                                          const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_hitPart3,
                                          const uint* rayFlags, const float* a_time, uint bounce, RandomGen* a_gen,
                                          FourierSpec* out_shadeColor)
{
  /*
  if(tid >= m_maxThreadId)
    return;
  const uint currRayFlags = *rayFlags;
  if(isDeadRay(currRayFlags))
    return;
    
  const uint32_t matId = extractMatId(currRayFlags);
  const float3 ray_dir = to_float3(*rayDirAndFar);
  
  const float4 data1  = *in_hitPart1;
  const float4 data2  = *in_hitPart2;
  const float4 lambda = *wavelengths;

  SurfaceHit hit;
  hit.pos  = to_float3(data1);
  hit.norm = to_float3(data2);
  hit.tang = to_float3(*in_hitPart3);
  hit.uv   = float2(data1.w, data2.w);
  
  const int bounceTmp = int(bounce); 
  const float4 rands = GetRandomNumbersLgts(tid, a_gen, bounceTmp); 
  const int lightId  = std::min(int(std::floor(rands.w * float(m_lights.size()))), int(m_lights.size() - 1u));
  RecordLightRndIfNeeded(bounce, rands); 

  if(lightId < 0) // no lights or invalid light id
  {
    *out_shadeColor = float4(0.0f, 0.0f, 0.0f, 0.0f);
    return;
  }
  
  const LightSample lSam = LightSampleRev(lightId, to_float3(rands), hit.pos);
  const float  hitDist   = std::sqrt(dot(hit.pos - lSam.pos, hit.pos - lSam.pos));

  const float3 shadowRayDir = normalize(lSam.pos - hit.pos); // explicitSam.direction;
  const float3 shadowRayPos = hit.pos + hit.norm * std::max(maxcomp(hit.pos), 1.0f)*5e-6f; // TODO: see Ray Tracing Gems, also use flatNormal for offset

  float time = *a_time;
  const bool   inIllumArea  = (dot(shadowRayDir, lSam.norm) < 0.0f) || lSam.isOmni || lSam.hasIES;
  const bool   needShade    = inIllumArea && !m_pAccelStruct->RayQuery_AnyHitMotion(to_float4(shadowRayPos, 0.0f), to_float4(shadowRayDir, hitDist*0.9995f), time); /// (!!!) expression-way, RT pipeline bug work around, if change check test_213
  RecordShadowHitIfNeeded(bounce, needShade);

  if(needShade) /// (!!!) expression-way to compute 'needShade', RT pipeline bug work around, if change check test_213
  {
    const BsdfEval bsdfV    = MaterialEval(matId, lambda, shadowRayDir, (-1.0f)*ray_dir, hit.norm, hit.tang, hit.uv);
    float cosThetaOut       = std::max(dot(shadowRayDir, hit.norm), 0.0f);
    
    float      lgtPdfW      = LightPdfSelectRev(lightId) * LightEvalPDF(lightId, shadowRayPos, shadowRayDir, lSam.pos, lSam.norm, lSam.pdf);
    float      misWeight    = (m_intergatorType == INTEGRATOR_MIS_PT) ? misWeightHeuristic(lgtPdfW, bsdfV.pdf) : 1.0f;
    const bool isDirect     = (m_lights[lightId].geomType == LIGHT_GEOM_DIRECT); 
    const bool isPoint      = (m_lights[lightId].geomType == LIGHT_GEOM_POINT); 
    
    if(isDirect)
    {
      misWeight = 1.0f;
      lgtPdfW   = 1.0f;
    }
    else if(isPoint)
      misWeight = 1.0f;

    const bool isDirectLight = !hasNonSpecular(currRayFlags);
    if((m_renderLayer == FB_DIRECT   && !isDirectLight) || 
       (m_renderLayer == FB_INDIRECT && isDirectLight)) // skip some number of bounces if this is set
      misWeight = 0.0f;
      
    
    const float4 lightColor = LightIntensity(lightId, lambda, shadowRayPos, shadowRayDir);
    *out_shadeColor = (lightColor * bsdfV.val / lgtPdfW) * cosThetaOut * misWeight;
  }
  else
    *out_shadeColor = float4(0.0f, 0.0f, 0.0f, 0.0f);
  */
}