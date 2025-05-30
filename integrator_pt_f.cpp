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

    kernel_NextBounceF(tid, depth, &hitPart1, &hitPart2, &hitPart3, &instId, &shadeColor,
                      &rayPosAndNear, &rayDirAndFar, &accumColor, &accumThroughput, &gen, &mis, &rayFlags);

    if(isDeadRay(rayFlags))
      break;
  }

  kernel_HitEnvironment(tid, &rayFlags, &rayDirAndFar, &mis, &accumThroughput,
                        &accumColor);

  kernel_ContributeToImage(tid, &rayFlags, channels, &accumColor, &gen, m_packedXY.data(), &wavelengths, out_color);
}



void Integrator::kernel_InitEyeRay2F(uint tid, float4* rayPosAndNear, float4* rayDirAndFar,
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
  if(tid >= m_maxThreadId)
    return;
  const uint currRayFlags = *rayFlags;
  if(isDeadRay(currRayFlags))
    return;
    
  const uint32_t matId = extractMatId(currRayFlags);
  const float3 ray_dir = to_float3(*rayDirAndFar);
  
  const float4 data1  = *in_hitPart1;
  const float4 data2  = *in_hitPart2;

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
    *out_shadeColor = FourierSpec(0.0f);
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
    const BsdfEvalF bsdfV    = MaterialEvalF(matId, shadowRayDir, (-1.0f)*ray_dir, hit.norm, hit.tang, hit.uv);
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
      
    
    const FourierSpec lightColor = LightIntensityF(lightId, shadowRayPos, shadowRayDir);
    *out_shadeColor = (lightColor * bsdfV.val / lgtPdfW) * cosThetaOut * misWeight;
  }
  else {
    *out_shadeColor = FourierSpec(0.0f);
  }
}

FourierSpec Integrator::LightIntensityF(uint a_lightId, float3 a_rayPos, float3 a_rayDir)
{
  FourierSpec lightColor;// = m_lights[a_lightId].intensity;  
  
  // get spectral data for light source
  //
  const uint specId = m_lights[a_lightId].specId;
  if(KSPEC_SPECTRAL_RENDERING != 0 && specId < 0xFFFFFFFF)
  {
    const uint2 data  = m_spec_offset_sz[specId];
    const uint offset = data.x;
    const uint size   = data.y;
    lightColor = FourierSpec(m_spec_values.data() + offset, size);
  }
  lightColor *= m_lights[a_lightId].mult;
  
  return lightColor;
}



void Integrator::kernel_NextBounceF(uint tid, uint bounce, const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_hitPart3, 
                       const uint* in_instId, const FourierSpec* in_shadeColor, float4* rayPosAndNear, float4* rayDirAndFar,
                       FourierSpec* accumColor, FourierSpec* accumThoroughput, RandomGen* a_gen, MisData* misPrev, uint* rayFlags)
{
  if(tid >= m_maxThreadId)
    return;
  const uint currRayFlags = *rayFlags;
  if(isDeadRay(currRayFlags))
    return;
    
  const uint32_t matId = extractMatId(currRayFlags);

  // process surface hit case
  //
  const float3 ray_dir = to_float3(*rayDirAndFar);
  const float3 ray_pos = to_float3(*rayPosAndNear);
  
  const float4 data1 = *in_hitPart1;
  const float4 data2 = *in_hitPart2;
  
  SurfaceHit hit;
  hit.pos  = to_float3(data1);
  hit.norm = to_float3(data2);
  hit.tang = to_float3(*in_hitPart3);
  hit.uv   = float2(data1.w, data2.w);

  const float hitDist = in_hitPart3->w;
  
  const MisData prevBounce = *misPrev;
  const float   prevPdfW   = prevBounce.matSamplePdf;

  // process light hit case
  //
  if(m_materials[matId].mtype == MAT_TYPE_LIGHT_SOURCE)
  {
    //const uint   texId     = m_materials[matId].texid[0];
    //const float2 texCoordT = mulRows2x4(m_materials[matId].row0[0], m_materials[matId].row1[0], hit.uv);
    //const float4 texColor  = m_textures[texId]->sample(texCoordT);
    const uint   lightId   = m_instIdToLightInstId[*in_instId]; 
    
    FourierSpec lightIntensity;

    if(lightId != 0xFFFFFFFF)
    {
      const float lightCos = dot(to_float3(*rayDirAndFar), to_float3(m_lights[lightId].norm));
      const float lightDirectionAtten = (lightCos < 0.0f || m_lights[lightId].geomType == LIGHT_GEOM_SPHERE) ? 1.0f : 0.0f;
      lightIntensity = LightIntensityF(lightId, ray_pos, to_float3(*rayDirAndFar))*lightDirectionAtten;
    }

    float misWeight = 1.0f;
    if(m_intergatorType == INTEGRATOR_MIS_PT) 
    {
      if(bounce > 0 && lightId != 0xFFFFFFFF)
      {
        const float lgtPdf  = LightPdfSelectRev(lightId) * LightEvalPDF(lightId, ray_pos, ray_dir, hit.pos, hit.norm, 1.0f);
        misWeight           = misWeightHeuristic(prevPdfW, lgtPdf);
        if (prevPdfW <= 0.0f) // specular bounce
          misWeight = 1.0f;
      }
    }
    else if(m_intergatorType == INTEGRATOR_SHADOW_PT && hasNonSpecular(currRayFlags)) {
      misWeight = 0.0f;
    }
    
    const bool isDirectLight  = !hasNonSpecular(currRayFlags);
    const bool isFirstNonSpec = (currRayFlags & RAY_FLAG_FIRST_NON_SPEC) != 0;
    if(m_renderLayer == FB_INDIRECT && (isDirectLight || isFirstNonSpec))
      misWeight = 0.0f;

    FourierSpec currAccumColor      = *accumColor;
    FourierSpec currAccumThroughput = *accumThoroughput;
    
    currAccumColor += currAccumThroughput * lightIntensity * misWeight;
   
    *accumColor = currAccumColor;
    *rayFlags   = currRayFlags | (RAY_FLAG_IS_DEAD | RAY_FLAG_HIT_LIGHT);
    return;
  }
  
  const uint bounceTmp    = bounce;
  const BsdfSampleF matSam = MaterialSampleAndEvalF(matId, tid, bounceTmp, a_gen, (-1.0f)*ray_dir, hit.norm, hit.tang, hit.uv, misPrev, currRayFlags);
  const FourierSpec bxdfVal    = matSam.val * (1.0f / std::max(matSam.pdf, 1e-20f));
  const float  cosTheta   = std::abs(dot(matSam.dir, hit.norm)); 

  MisData nextBounceData      = *misPrev;        // remember current pdfW for next bounce
  nextBounceData.matSamplePdf = (matSam.flags & RAY_EVENT_S) != 0 ? -1.0f : matSam.pdf; 
  nextBounceData.cosTheta     = cosTheta;   
  *misPrev                    = nextBounceData;

  if(m_intergatorType == INTEGRATOR_STUPID_PT)
  {
    *accumThoroughput *= cosTheta * bxdfVal; 
  }
  else if(m_intergatorType == INTEGRATOR_SHADOW_PT || m_intergatorType == INTEGRATOR_MIS_PT)
  {
    const FourierSpec currThoroughput = *accumThoroughput;
    const FourierSpec shadeColor      = *in_shadeColor;
    FourierSpec currAccumColor        = *accumColor;

    currAccumColor += currThoroughput * shadeColor;
    *accumColor       = currAccumColor;
    *accumThoroughput = currThoroughput*cosTheta*bxdfVal; 
  }

  // compute point on the other side of the surface in case of transmission
  if((matSam.flags & RAY_EVENT_T) != 0)
  {
    hit.pos = hit.pos + hitDist * ray_dir * 2 * 1e-6f;
  }  

  *rayPosAndNear = to_float4(OffsRayPos(hit.pos, hit.norm, matSam.dir), 0.0f); // todo: use flatNormal for offset
  *rayDirAndFar  = to_float4(matSam.dir, FLT_MAX);
  
  uint nextFlags = ((currRayFlags & ~RAY_FLAG_FIRST_NON_SPEC) | matSam.flags); // always force reset RAY_FLAG_FIRST_NON_SPEC;
  if(m_renderLayer == FB_DIRECT && hasNonSpecular(currRayFlags))   // NOTE: use currRayFlags for check, not nextFlags because of MIS: a ray may hit light source in next bounce
    nextFlags |= RAY_FLAG_IS_DEAD;                                 //       but if we already have non specular bounce previously, definitely can stop  
  else if(!hasNonSpecular(currRayFlags) && hasNonSpecular(nextFlags))
    nextFlags |= RAY_FLAG_FIRST_NON_SPEC;
  *rayFlags      = nextFlags;                                   
}