#include "integrator_dr.h"
#include "utils.h"

#include "include/cmaterial.h"
#include "include/cmat_gltf.h"
#include "include/cmat_conductor.h"
#include "include/cmat_glass.h"
#include "include/cmat_diffuse.h"
#include "include/cmat_plastic.h"

#include <chrono>
#include <string>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

BsdfEval IntegratorDR::MaterialEval(uint a_materialId, float4 wavelengths, float3 l, float3 v, float3 n, float3 tan, float2 tc, const float* a_data)
{
  BsdfEval res;
  {
    res.val = float4(0,0,0,0);
    res.pdf = 0.0f;
  }

  MatIdWeight currMat = make_id_weight(a_materialId, 1.0f);
  MatIdWeight material_stack[KSPEC_BLEND_STACK_SIZE];
  if(KSPEC_MAT_TYPE_BLEND != 0)
    material_stack[0] = currMat;
  int top = 0;
  bool needPop = false;

  do
  {
    if(KSPEC_MAT_TYPE_BLEND != 0)
    {
      if(needPop)
      {
        top--;
        currMat = material_stack[std::max(top, 0)];
      }
      else
        needPop = true; // if not blend, pop on next iter
    } 
    
    // BSDF is multiplied (outside) by old cosThetaOut.
    // When normal map is enables this becames wrong because normal is changed;
    // First : return cosThetaOut in sam;
    // Second: apply cos(theta2)/cos(theta1) to cos(theta1) to get cos(theta2)
    //
    const float3 geomNormal = n;
          float3 shadeNormal = n;
    float bumpCosMult = 1.0f; 
    const uint normalMapId = as_uint(m_materials[currMat.id].data[UINT_NMAP_ID]);
    if(KSPEC_BUMP_MAPPING != 0 && normalMapId != 0xFFFFFFFF) 
    {
      shadeNormal = BumpMapping(normalMapId, currMat.id, geomNormal, tan, tc);
      const float3 lDir     = l;     
      const float  clampVal = 1e-6f;  
      const float cosThetaOut1 = std::max(dot(lDir, geomNormal),  0.0f);
      const float cosThetaOut2 = std::max(dot(lDir, shadeNormal), 0.0f);
      bumpCosMult              = cosThetaOut2 / std::max(cosThetaOut1, clampVal);
      if (cosThetaOut1 <= 0.0f)
        bumpCosMult = 0.0f;
    }

    const float2 texCoordT = mulRows2x4(m_materials[currMat.id].row0[0], m_materials[currMat.id].row1[0], tc);
    const uint   mtype     = as_uint(m_materials[currMat.id].data[UINT_MTYPE]);

    BsdfEval currVal;
    {
      currVal.val = float4(0,0,0,0);
      currVal.pdf   = 0.0f;
    }
    switch(mtype)
    {
      case MAT_TYPE_GLTF:
      if(KSPEC_MAT_TYPE_GLTF != 0)
      {
        const uint   texId     = as_uint(m_materials[currMat.id].data[GLTF_UINT_TEXID0]);
        const float4 texColor  = HydraTex2DFetch(texId, texCoordT, a_data); //m_textures[texId]->sample(texCoordT);
        const float4 color     = (m_materials[currMat.id].colors[GLTF_COLOR_BASE]) * texColor;
        gltfEval(m_materials.data() + currMat.id, l, v, shadeNormal, tc, color, &currVal);

        res.val += currVal.val * currMat.weight * bumpCosMult;
        res.pdf += currVal.pdf * currMat.weight;

        break;
      }
      case MAT_TYPE_GLASS:
      if(KSPEC_MAT_TYPE_GLASS != 0)
      {
        glassEval(m_materials.data() + currMat.id, l, v, geomNormal, tc, float3(0,0,0), &currVal);

        res.val += currVal.val * currMat.weight * bumpCosMult;
        res.pdf += currVal.pdf * currMat.weight;
        break;
      }
      case MAT_TYPE_CONDUCTOR: 
      if(KSPEC_MAT_TYPE_CONDUCTOR != 0)
      {
        const uint   texId     = as_uint(m_materials[currMat.id].data[CONDUCTOR_TEXID0]);
        const float3 alphaTex  = to_float3(m_textures[texId]->sample(texCoordT));
        const float2 alpha     = float2(m_materials[currMat.id].data[CONDUCTOR_ROUGH_V], m_materials[currMat.id].data[CONDUCTOR_ROUGH_U]);

        if(!trEffectivelySmooth(alpha))
        {
          const float4 etaSpec = SampleMatParamSpectrum(currMat.id, wavelengths, CONDUCTOR_ETA, CONDUCTOR_ETA_SPECID);
          const float4 kSpec   = SampleMatParamSpectrum(currMat.id, wavelengths, CONDUCTOR_K, CONDUCTOR_K_SPECID);
          conductorRoughEval(m_materials.data() + currMat.id, etaSpec, kSpec, l, v, shadeNormal, tc, alphaTex, &currVal);
        }

        res.val += currVal.val * currMat.weight * bumpCosMult;
        res.pdf += currVal.pdf * currMat.weight;
        break;
      }
      case MAT_TYPE_DIFFUSE:
      if(KSPEC_MAT_TYPE_DIFFUSE != 0)
      {
        const uint   texId       = as_uint(m_materials[currMat.id].data[DIFFUSE_TEXID0]);
        const float4 texColor    = HydraTex2DFetch(texId, texCoordT, a_data); //m_textures[texId]->sample(texCoordT))
        const float4 color       = texColor;

        const float4 reflSpec    = SampleMatColorParamSpectrum(currMat.id, wavelengths, DIFFUSE_COLOR, DIFFUSE_SPECID);

        diffuseEval(m_materials.data() + currMat.id, reflSpec, l, v, shadeNormal, tc, color, &currVal);

        res.val += currVal.val * currMat.weight * bumpCosMult;
        res.pdf += currVal.pdf * currMat.weight;
        break;
      }
      case MAT_TYPE_PLASTIC:
      if(KSPEC_MAT_TYPE_PLASTIC != 0)
      {
        const uint   texId       = as_uint(m_materials[currMat.id].data[PLASTIC_COLOR_TEXID]);
        const float4 texColor    = HydraTex2DFetch(texId, texCoordT, a_data); // m_textures[texId]->sample(texCoordT);
        const float4 color       = texColor;

        float4 reflSpec    = SampleMatColorParamSpectrum(currMat.id, wavelengths, PLASTIC_COLOR, PLASTIC_COLOR_SPECID);
        if(m_spectral_mode == 0)
          reflSpec *= color;
        const uint precomp_id = as_uint(m_materials[currMat.id].data[PLASTIC_PRECOMP_ID]);
        plasticEval(m_materials.data() + currMat.id, reflSpec, l, v, shadeNormal, tc, &currVal, 
                    m_precomp_coat_transmittance.data() + precomp_id * MI_ROUGH_TRANSMITTANCE_RES);

        res.val += currVal.val * currMat.weight * bumpCosMult;
        res.pdf += currVal.pdf * currMat.weight;
        break;
        break;
      }
      case MAT_TYPE_BLEND:
      if(KSPEC_MAT_TYPE_BLEND != 0)
      {
        auto childMats = MaterialBlendEval(currMat, wavelengths, l, v, geomNormal, tc);
        currMat = childMats.first;
        needPop = false;                        // we already put 'childMats.first' in 'currMat'
        if(top + 1 <= KSPEC_BLEND_STACK_SIZE)
        {
          material_stack[top] = childMats.second; // remember second mat in stack
          top++;
        }
        break;
      }
      default:
        break;
    }

  } while(KSPEC_MAT_TYPE_BLEND != 0 && top > 0);

  return res;
}

BsdfSample IntegratorDR::MaterialSampleAndEval(uint a_materialId, float4 wavelengths, RandomGen* a_gen, float3 v, float3 n, float3 tan, float2 tc, 
                                               MisData* a_misPrev, const uint a_currRayFlags, const float* a_data)
{
  BsdfSample res;
  {
    res.val   = float4(0, 0, 0, 0);
    res.pdf   = 1.0f;
    res.dir   = float3(0,1,0);
    res.flags = a_currRayFlags;
  }

  uint32_t currMatId = a_materialId;
  uint     mtype     = as_uint(m_materials[currMatId].data[UINT_MTYPE]);
  while(KSPEC_MAT_TYPE_BLEND != 0 && mtype == MAT_TYPE_BLEND)
  {
    currMatId = MaterialBlendSampleAndEval(currMatId, wavelengths, a_gen, v, n, tc, a_misPrev, &res);
    mtype     = as_uint(m_materials[currMatId].data[UINT_MTYPE]);
  }
  
  // BSDF is multiplied (outside) by cosThetaOut1.
  // When normal map is enables this becames wrong because normal is changed;
  // First : return cosThetaOut in sam;
  // Second: apply cos(theta2)/cos(theta1) to cos(theta1) to get cos(theta2)
  //
  const uint normalMapId   = as_uint(m_materials[currMatId].data[UINT_NMAP_ID]);
  const float3 geomNormal  = n;
        float3 shadeNormal = n;

  if(KSPEC_BUMP_MAPPING != 0 && normalMapId != 0xFFFFFFFF)
    shadeNormal = BumpMapping(normalMapId, currMatId, geomNormal, tan, tc);

  const float2 texCoordT = mulRows2x4(m_materials[currMatId].row0[0], m_materials[currMatId].row1[0], tc);
  const float4 rands     = rndFloat4_Pseudo(a_gen);

  switch(mtype)
  {
    case MAT_TYPE_GLTF:
    if(KSPEC_MAT_TYPE_GLTF != 0)
    {
      const uint   texId    = as_uint(m_materials[currMatId].data[GLTF_UINT_TEXID0]);
      const float4 texColor = HydraTex2DFetch(texId, texCoordT, a_data); // m_textures[texId]->sample(texCoordT);
      const float4 color    = m_materials[currMatId].colors[GLTF_COLOR_BASE]*texColor;
      gltfSampleAndEval(m_materials.data() + currMatId, rands, v, shadeNormal, tc, color, &res);
    }
    break;
    case MAT_TYPE_GLASS: 
    if(KSPEC_MAT_TYPE_GLASS != 0)
    {
      glassSampleAndEval(m_materials.data() + currMatId, rands, v, geomNormal, tc, &res, a_misPrev);
    }
    break;
    case MAT_TYPE_CONDUCTOR:
    if(KSPEC_MAT_TYPE_CONDUCTOR != 0)
    {
      const uint   texId     = as_uint(m_materials[currMatId].data[CONDUCTOR_TEXID0]);
      const float3 alphaTex  = to_float3(m_textures[texId]->sample(texCoordT));
      
      const float2 alpha   = float2(m_materials[currMatId].data[CONDUCTOR_ROUGH_V], m_materials[currMatId].data[CONDUCTOR_ROUGH_U]);
      const float4 etaSpec = SampleMatParamSpectrum(currMatId, wavelengths, CONDUCTOR_ETA, CONDUCTOR_ETA_SPECID);
      const float4 kSpec   = SampleMatParamSpectrum(currMatId, wavelengths, CONDUCTOR_K, CONDUCTOR_K_SPECID);
      if(trEffectivelySmooth(alpha))
        conductorSmoothSampleAndEval(m_materials.data() + currMatId, etaSpec, kSpec, rands, v, shadeNormal, tc, &res);
      else
        conductorRoughSampleAndEval(m_materials.data() + currMatId, etaSpec, kSpec, rands, v, shadeNormal, tc, alphaTex, &res);
    }
    break;
    case MAT_TYPE_DIFFUSE:
    if(KSPEC_MAT_TYPE_DIFFUSE != 0)
    {
      const uint   texId       = as_uint(m_materials[currMatId].data[DIFFUSE_TEXID0]);
      const float4 texColor    = HydraTex2DFetch(texId, texCoordT, a_data); //m_textures[texId]->sample(texCoordT);
      const float4 color       = texColor;
      const float4 reflSpec    = SampleMatColorParamSpectrum(currMatId, wavelengths, DIFFUSE_COLOR, DIFFUSE_SPECID);

      diffuseSampleAndEval(m_materials.data() + currMatId, reflSpec, rands, v, shadeNormal, tc, color, &res);
    }
    break;
    case MAT_TYPE_PLASTIC:
    if(KSPEC_MAT_TYPE_PLASTIC != 0)
    {
      const uint   texId       = as_uint(m_materials[currMatId].data[PLASTIC_COLOR_TEXID]);
      const float4 texColor    = HydraTex2DFetch(texId, texCoordT, a_data); //m_textures[texId]->sample(texCoordT);
      const float4 color       = texColor;

      float4 reflSpec    = SampleMatColorParamSpectrum(currMatId, wavelengths, PLASTIC_COLOR, PLASTIC_COLOR_SPECID);
      if(m_spectral_mode == 0)
        reflSpec *= color;

      const uint precomp_id = as_uint(m_materials[currMatId].data[PLASTIC_PRECOMP_ID]);

      plasticSampleAndEval(m_materials.data() + currMatId, reflSpec, rands, v, shadeNormal, tc, &res,
                           m_precomp_coat_transmittance.data() + precomp_id * MI_ROUGH_TRANSMITTANCE_RES);
    }
    break;
    default:
    break;
  }
  
  // BSDF is multiplied (outside) by cosThetaOut1.
  // When normal map is enables this becames wrong because normal is changed;
  // First : return cosThetaOut in sam;
  // Second: apply cos(theta2)/cos(theta1) to cos(theta1) to get cos(theta2)
  //
  if(KSPEC_BUMP_MAPPING != 0 && normalMapId != 0xFFFFFFFF)
  {
    const float cosThetaOut1 = std::abs(dot(res.dir, geomNormal));
    const float cosThetaOut2 = std::abs(dot(res.dir, shadeNormal));
    res.val *= cosThetaOut2 / std::max(cosThetaOut1, 1e-10f);
  }

  return res;
}

void IntegratorDR::kernel_SampleLightSource(uint tid, const float4* rayPosAndNear, const float4* rayDirAndFar, 
                                            const float4* wavelengths, const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_hitPart3,
                                            const uint* rayFlags, uint bounce,
                                            RandomGen* a_gen, float4* out_shadeColor,  const float* a_data)
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
  const float4 lambda = *wavelengths;

  SurfaceHit hit;
  hit.pos  = to_float3(data1);
  hit.norm = to_float3(data2);
  hit.tang = to_float3(*in_hitPart3);
  hit.uv   = float2(data1.w, data2.w);

  const float2 rands = rndFloat2_Pseudo(a_gen); // don't use single rndFloat4 (!!!)
  const float rndId  = rndFloat1_Pseudo(a_gen); // don't use single rndFloat4 (!!!)
  const int lightId  = std::min(int(std::floor(rndId * float(m_lights.size()))), int(m_lights.size() - 1u));

  if(lightId < 0) // no lights or invalid light id
  {
    *out_shadeColor = float4(0.0f, 0.0f, 0.0f, 0.0f);
    return;
  }
  
  const LightSample lSam = LightSampleRev(lightId, rands, hit.pos);
  const float  hitDist   = std::sqrt(dot(hit.pos - lSam.pos, hit.pos - lSam.pos));

  const float3 shadowRayDir = normalize(lSam.pos - hit.pos); // explicitSam.direction;
  const float3 shadowRayPos = hit.pos + hit.norm*std::max(maxcomp(hit.pos), 1.0f)*5e-6f; // TODO: see Ray Tracing Gems, also use flatNormal for offset
  const bool   inShadow     = m_pAccelStruct->RayQuery_AnyHit(to_float4(shadowRayPos, 0.0f), to_float4(shadowRayDir, hitDist*0.9995f));
  const bool   inIllumArea  = (dot(shadowRayDir, lSam.norm) < 0.0f) || lSam.isOmni;

  if(!inShadow && inIllumArea) 
  {
    const BsdfEval bsdfV    = MaterialEval(matId, lambda, shadowRayDir, (-1.0f)*ray_dir, hit.norm, hit.tang, hit.uv, a_data);
    const float cosThetaOut = std::max(dot(shadowRayDir, hit.norm), 0.0f);
    
    float      lgtPdfW      = LightPdfSelectRev(lightId) * LightEvalPDF(lightId, shadowRayPos, shadowRayDir, lSam.pos, lSam.norm);
    float      misWeight    = (m_intergatorType == INTEGRATOR_MIS_PT) ? misWeightHeuristic(lgtPdfW, bsdfV.pdf) : 1.0f;
    const bool isDirect     = (m_lights[lightId].geomType == LIGHT_GEOM_DIRECT); 
    
    if(isDirect)
    {
      misWeight = 1.0f;
      lgtPdfW   = 1.0f;
    }

    if(m_skipBounce >= 1 && int(bounce) < int(m_skipBounce)-1) // skip some number of bounces if this is set
      misWeight = 0.0f;
    
    const float4 lightColor = GetLightSourceIntensity(lightId, wavelengths);
    *out_shadeColor = (lightColor * bsdfV.val / lgtPdfW) * cosThetaOut * misWeight;
  }
  else
    *out_shadeColor = float4(0.0f, 0.0f, 0.0f, 0.0f);
}

void IntegratorDR::kernel_NextBounce(uint tid, uint bounce, const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_hitPart3, const uint* in_instId,
                                     const float4* in_shadeColor, float4* rayPosAndNear, float4* rayDirAndFar, const float4* wavelengths,
                                     float4* accumColor, float4* accumThoroughput, RandomGen* a_gen, MisData* misPrev, uint* rayFlags, const float* a_data)
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
  const float4 lambda  = *wavelengths;
  
  const float4 data1 = *in_hitPart1;
  const float4 data2 = *in_hitPart2;
  
  SurfaceHit hit;
  hit.pos  = to_float3(data1);
  hit.norm = to_float3(data2);
  hit.tang = to_float3(*in_hitPart3);
  hit.uv   = float2(data1.w, data2.w);
  
  const MisData prevBounce = *misPrev;
  const float   prevPdfW   = prevBounce.matSamplePdf;

  // process light hit case
  //
  if(as_uint(m_materials[matId].data[UINT_MTYPE]) == MAT_TYPE_LIGHT_SOURCE)
  {
    const uint texId       = as_uint(m_materials[matId].data[EMISSION_TEXID0]);
    const float2 texCoordT = mulRows2x4(m_materials[matId].row0[0], m_materials[matId].row1[0], hit.uv);
    float4 texColor   = m_textures[texId]->sample(texCoordT);
    float4 lightColor = m_materials[matId].colors[EMISSION_COLOR];
    float  lightMult  = m_materials[matId].data[EMISSION_MULT];

    float4 lightIntensity = lightColor * texColor * lightMult;
    if(KSPEC_SPECTRAL_RENDERING != 0 && m_spectral_mode != 0)
    {
      const uint specId = as_uint(m_materials[matId].data[EMISSION_SPECID0]);
      if(specId < 0xFFFFFFFF)
      {
        const uint2 data  = m_spec_offset_sz[specId];
        const uint offset = data.x;
        const uint size   = data.y;
        lightColor = SampleSpectrum(m_wavelengths.data() + offset, m_spec_values.data() + offset, *wavelengths, size);
      }
      lightIntensity = lightColor * lightMult;
    }

    const uint lightId = m_instIdToLightInstId[*in_instId]; //m_materials[matId].data[UINT_LIGHTID];
    
    float lightCos = 1.0f;
    float lightDirectionAtten = 1.0f;
    if(lightId != 0xFFFFFFFF)
    {
      lightCos = dot(to_float3(*rayDirAndFar), to_float3(m_lights[lightId].norm));
      lightDirectionAtten = (lightCos < 0.0f || m_lights[lightId].geomType == LIGHT_GEOM_SPHERE) ? 1.0f : 0.0f;
    }

    float misWeight = 1.0f;
    if(m_intergatorType == INTEGRATOR_MIS_PT) 
    {
      if(bounce > 0)
      {
        if(lightId != 0xFFFFFFFF)
        {
          const float lgtPdf  = LightPdfSelectRev(lightId) * LightEvalPDF(lightId, ray_pos, ray_dir, hit.pos, hit.norm);
          misWeight           = misWeightHeuristic(prevPdfW, lgtPdf);
          if (prevPdfW <= 0.0f) // specular bounce
            misWeight = 1.0f;
        }
      }
    }
    else if(m_intergatorType == INTEGRATOR_SHADOW_PT && hasNonSpecular(currRayFlags))
      misWeight = 0.0f;
    
    if(m_skipBounce >= 1 && bounce < m_skipBounce) // skip some number of bounces if this is set
      misWeight = 0.0f;

    float4 currAccumColor      = *accumColor;
    float4 currAccumThroughput = *accumThoroughput;
    
    currAccumColor += currAccumThroughput * lightIntensity * misWeight * lightDirectionAtten;
    // currAccumColor.x += currAccumThroughput.x * lightIntensity.x * misWeight * lightDirectionAtten;
    // currAccumColor.y += currAccumThroughput.y * lightIntensity.y * misWeight * lightDirectionAtten;
    // currAccumColor.z += currAccumThroughput.z * lightIntensity.z * misWeight * lightDirectionAtten;
    // if(bounce > 0)
    //   currAccumColor.w *= prevPdfA;
    
    *accumColor = currAccumColor;
    *rayFlags   = currRayFlags | (RAY_FLAG_IS_DEAD | RAY_FLAG_HIT_LIGHT);
    return;
  }
  
  const BsdfSample matSam = MaterialSampleAndEval(matId, lambda, a_gen, (-1.0f)*ray_dir, hit.norm, hit.tang, hit.uv, misPrev, currRayFlags, a_data);
  const float4 bxdfVal    = matSam.val * (1.0f / std::max(matSam.pdf, 1e-20f));
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
    const float4 currThoroughput = *accumThoroughput;
    const float4 shadeColor      = *in_shadeColor;
    float4 currAccumColor        = *accumColor;

    currAccumColor += currThoroughput * shadeColor;
    // currAccumColor.x += currThoroughput.x * shadeColor.x;
    // currAccumColor.y += currThoroughput.y * shadeColor.y;
    // currAccumColor.z += currThoroughput.z * shadeColor.z;
    // if(bounce > 0)
    //   currAccumColor.w *= prevPdfA;

    *accumColor       = currAccumColor;
    *accumThoroughput = currThoroughput*cosTheta*bxdfVal; 
  }

  *rayPosAndNear = to_float4(OffsRayPos(hit.pos, hit.norm, matSam.dir), 0.0f); // todo: use flatNormal for offset
  *rayDirAndFar  = to_float4(matSam.dir, FLT_MAX);
  *rayFlags      = currRayFlags | matSam.flags;
}

float4 IntegratorDR::PathTrace(uint tid, uint channels, float* out_color, const float* a_data)
{
  float4 accumColor, accumThroughput;
  float4 rayPosAndNear, rayDirAndFar;
  float4 wavelengths;
  RandomGen gen; 
  MisData   mis;
  uint      rayFlags;
  kernel_InitEyeRay2(tid, m_packedXY.data(), &rayPosAndNear, &rayDirAndFar, &wavelengths, &accumColor, &accumThroughput, &gen, &rayFlags, &mis);

  for(uint depth = 0; depth < m_traceDepth; depth++) 
  {
    float4   shadeColor, hitPart1, hitPart2, hitPart3;
    uint instId;
    kernel_RayTrace2(tid, &rayPosAndNear, &rayDirAndFar, &hitPart1, &hitPart2, &hitPart3, &instId, &rayFlags);
    if(isDeadRay(rayFlags))
      break;
    
    kernel_SampleLightSource(tid, &rayPosAndNear, &rayDirAndFar, &wavelengths, &hitPart1, &hitPart2, &hitPart3, &rayFlags, depth,
                             &gen, &shadeColor, a_data);

    kernel_NextBounce(tid, depth, &hitPart1, &hitPart2, &hitPart3, &instId, &shadeColor,
                      &rayPosAndNear, &rayDirAndFar, &wavelengths, &accumColor, &accumThroughput, &gen, &mis, &rayFlags, a_data);

    if(isDeadRay(rayFlags))
      break;
  }

  kernel_HitEnvironment(tid, &rayFlags, &rayDirAndFar, &mis, &accumThroughput,
                        &accumColor);

  kernel_ContributeToImage(tid, channels, &accumColor, &gen, m_packedXY.data(), &wavelengths, out_color);

  return accumColor;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void IntegratorDR::kernel_CalcRayColor(uint tid, const Lite_Hit* in_hit, const float2* bars, float4* finalColor, const uint* in_pakedXY, float* out_color, const float* a_data)
{ 
  if(tid >= m_maxThreadId)
    return;

  const Lite_Hit hit = *in_hit;
  if(hit.geomId == -1)
  {
    out_color[tid] = 0;
    return;
  }

  const uint32_t matId  = m_matIdByPrimId[m_matIdOffsets[hit.geomId] + hit.primId];
  const float4 mdata    = m_materials[matId].colors[GLTF_COLOR_BASE];
  const float2 uv       = *bars;

  const uint triOffset  = m_matIdOffsets[hit.geomId];
  const uint vertOffset = m_vertOffset  [hit.geomId];

  const uint A = m_triIndices[(triOffset + hit.primId)*3 + 0];
  const uint B = m_triIndices[(triOffset + hit.primId)*3 + 1];
  const uint C = m_triIndices[(triOffset + hit.primId)*3 + 2];
  const float4 data1 = (1.0f - uv.x - uv.y)*m_vNorm4f[A + vertOffset] + uv.y*m_vNorm4f[B + vertOffset] + uv.x*m_vNorm4f[C + vertOffset];
  const float4 data2 = (1.0f - uv.x - uv.y)*m_vTang4f[A + vertOffset] + uv.y*m_vTang4f[B + vertOffset] + uv.x*m_vTang4f[C + vertOffset];
  float3 hitNorm     = to_float3(data1);
  float3 hitTang     = to_float3(data2);
  float2 hitTexCoord = float2(data1.w, data2.w);

  const uint   texId     = as_uint(m_materials[matId].data[GLTF_UINT_TEXID0]);
  const float2 texCoordT = mulRows2x4(m_materials[matId].row0[0], m_materials[matId].row1[0], hitTexCoord);
  const float4 texColor  = HydraTex2DFetch(texId, texCoordT, a_data); 
  const float3 color     = mdata.w > 0.0f ? clamp(float3(mdata.w,mdata.w,mdata.w), 0.0f, 1.0f) : to_float3(mdata*texColor);

  const uint XY = in_pakedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;

  (*finalColor) = to_float4(color, 0);

  out_color[(y*m_winWidth+x)*4 + 0] = color.x;
  out_color[(y*m_winWidth+x)*4 + 1] = color.y;
  out_color[(y*m_winWidth+x)*4 + 2] = color.z;
  out_color[(y*m_winWidth+x)*4 + 3] = 0.0f;
}


float4 IntegratorDR::CastRayDR(uint tid, uint channels, float* out_color, const float* a_data)
{
  float4 rayPosAndNear, rayDirAndFar;
  kernel_InitEyeRay(tid, m_packedXY.data(), &rayPosAndNear, &rayDirAndFar);

  Lite_Hit hit; 
  float2   baricentrics; 
  if(!kernel_RayTrace(tid, &rayPosAndNear, &rayDirAndFar, &hit, &baricentrics))
    return float4(0,0,0,0);
  
  float4 finalColor;
  kernel_CalcRayColor(tid, &hit, &baricentrics, &finalColor, m_packedXY.data(), out_color, a_data);
  return finalColor;
}


extern double __enzyme_autodiff(void*, ...);
int enzyme_const, enzyme_dup, enzyme_out;

float Loss(IntegratorDR* __restrict__ pIntegrator,
           const float*  __restrict__ a_refImg,
                 float*  __restrict__ out_color,
           const float*  __restrict__ a_data, 
           const uint* in_pakedXY, 
           uint tid, uint channels, uint pitch)
{
  const uint XY = in_pakedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;

  //float4 colorRend = pIntegrator->PathTrace(tid, channels, out_color, a_data);
  float4 colorRend = pIntegrator->CastRayDR(tid, channels, out_color, a_data);

  float4 colorRef  = float4(a_refImg[(y*pitch+x)*channels + 0], 
                            a_refImg[(y*pitch+x)*channels + 1], 
                            a_refImg[(y*pitch+x)*channels + 2], 0.0f);

  float4 diff = colorRend - colorRef;

  return LiteMath::dot3(diff, diff);
}                     


void IntegratorDR::PathTraceDR(uint tid, uint channels, float* out_color, uint a_passNum,
                               const float* a_refImg, const float* a_data, float* a_dataGrad, size_t a_gradSize)
{
  memset(a_dataGrad, 0, sizeof(float)*a_gradSize);

  // init separate gradient for each thread
  //
  //std::vector<float> grads[MAXTHREADS];
  //for(int i=0;i<MAXTHREADS;i++)
  //   std::fill(grads[i].begin(), grads[i].end(), 0.0f);

  double avgLoss = 0.0;

  ConsoleProgressBar progress(tid);
  progress.Start();
  auto start = std::chrono::high_resolution_clock::now();
  //#ifndef _DEBUG
  //#pragma omp parallel for default(shared) // num_threads(MAXTHREADS)
  //#endif
  for (int i = 0; i < int(tid); ++i) {
    float lossVal = 0.0f;
    for (int j = 0; j < int(a_passNum); ++j) {
      //PathTrace(uint(i), channels, out_color);
      lossVal += Loss(this, 
                      a_refImg, 
                      out_color, 
                      a_data,
                      m_packedXY.data(), 
                      uint(i), channels, m_winWidth);

      //__enzyme_autodiff((void*)Loss, 
      //                   enzyme_const, this,
      //                   enzyme_const, a_refImg,
      //                   enzyme_const, out_color,
      //                   enzyme_dup,   a_data, a_dataGrad,
      //                   enzyme_const, m_packedXY.data(),
      //                   enzyme_const, uint(i),
      //                   enzyme_const, channels,
      //                   enzyme_const, m_winWidth);
    }
    avgLoss += double(lossVal)/double(a_passNum);
    progress.Update();
  }
  progress.Done();
  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;

  // accumulate gradient from different threads (parallel reduction/hist)
  //

  //for(int i=0;i<MAXTHREADS;i++) 
  //  for(size_t j=0;j<a_gradSize; j++)
  //    a_dataGrad[j] += grads[i][j];

  avgLoss /= double(m_winWidth*m_winHeight);
  std::cout << "avgLoss = " << avgLoss << std::endl;
  
  //std::ofstream fout("z_grad.txt");
  //for(size_t i=0; i<a_gradSize; i++)
  //  fout << a_dataGrad[i]/float(a_passNum) << std::endl;
  //fout.close();
}

static inline int4 bilinearOffsets(const float ffx, const float ffy, const int w, const int h)
{
	const int sx = (ffx > 0.0f) ? 1 : -1;
	const int sy = (ffy > 0.0f) ? 1 : -1;

	const int px = (int)(ffx);
	const int py = (int)(ffy);

	int px_w0, px_w1, py_w0, py_w1;
  // wrap
	{
		px_w0 = px        % w;
		px_w1 = (px + sx) % w;

		px_w0 = (px_w0 < 0) ? px_w0 + w : px_w0;
		px_w1 = (px_w1 < 0) ? px_w1 + w : px_w1;
	}

  // wrap
	{
		py_w0 = py        % h;
		py_w1 = (py + sy) % h;

		py_w0 = (py_w0 < 0) ? py_w0 + h : py_w0;
		py_w1 = (py_w1 < 0) ? py_w1 + h : py_w1;
	}

	const int offset0 = py_w0*w + px_w0;
	const int offset1 = py_w0*w + px_w1;
	const int offset2 = py_w1*w + px_w0;
	const int offset3 = py_w1*w + px_w1;

	return int4(offset0, offset1, offset2, offset3);
}

float4 IntegratorDR::HydraTex2DFetch(uint texId, float2 a_uv, const float* tex_data)
{
  if(texId == 1 && tex_data != nullptr) 
  {
    const float m_fw = 256.0f;
    const float m_fh = 256.0f;
    const int tex_width  = 256;
    const int tex_height = 256;
    
    float ffx = a_uv.x * m_fw - 0.5f; // a_uv should not be very large, so that the float does not overflow later. 
    float ffy = a_uv.y * m_fh - 0.5f; // This is left to the responsibility of the top level.
    
    if (ffx < 0) ffx = 0.0f;
    if (ffy < 0) ffy = 0.0f;
    
    // Calculate the weights for each pixel
    //
    const int   px = (int)(ffx);
    const int   py = (int)(ffy);
    
    const float fx  = std::abs(ffx - (float)px);
    const float fy  = std::abs(ffy - (float)py);
    const float fx1 = 1.0f - fx;
    const float fy1 = 1.0f - fy;
    
    const float w1 = fx1 * fy1;
    const float w2 = fx  * fy1;
    const float w3 = fx1 * fy;
    const float w4 = fx  * fy;
    
    const int4 offsets = bilinearOffsets(ffx, ffy, tex_width, tex_height);
    
    const float4 f1    = float4(tex_data[offsets.x*4+0], tex_data[offsets.x*4+1], tex_data[offsets.x*4+2], tex_data[offsets.x*4+3]);
    const float4 f2    = float4(tex_data[offsets.y*4+0], tex_data[offsets.y*4+1], tex_data[offsets.y*4+2], tex_data[offsets.y*4+3]);
    const float4 f3    = float4(tex_data[offsets.z*4+0], tex_data[offsets.z*4+1], tex_data[offsets.z*4+2], tex_data[offsets.z*4+3]);
    const float4 f4    = float4(tex_data[offsets.w*4+0], tex_data[offsets.w*4+1], tex_data[offsets.w*4+2], tex_data[offsets.w*4+3]);
    // Calculate the weighted sum of pixels (for each color channel)
    //
    const float outr = f1.x * w1 + f2.x * w2 + f3.x * w3 + f4.x * w4;
    const float outg = f1.y * w1 + f2.y * w2 + f3.y * w3 + f4.y * w4;
    const float outb = f1.z * w1 + f2.z * w2 + f3.z * w3 + f4.z * w4;
    const float outa = f1.w * w1 + f2.w * w2 + f3.w * w3 + f4.w * w4;
    
    return float4(outr, outg, outb, outa);

    //return float4(1.0, 0.1, 0.1,1);
  }
  else
    return m_textures[texId]->sample(a_uv);
}