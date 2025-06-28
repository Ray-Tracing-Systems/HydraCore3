#include "integrator_pt.h"
#include "include/crandom.h"

#include "fourier/fspec.h"
#include "specn.h"

#include "fourier/cmat_diffuse.h"
#include "fourier/cmat_conductor.h"
#include "fourier/cmat_film.h"

#include "include/cmat_diffuse.h"
#include "include/cmat_conductor.h"
#include "include/cmat_film.h"

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;




BsdfEvalF Integrator::MaterialEvalF(uint a_materialId, float3 l, float3 v, float3 n, float3 tan, float2 tc)
{
  BsdfEvalF res;
  {
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
    const uint normalMapId = m_materials[currMat.id].texid[1];
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
    const uint   texId     = m_materials[currMat.id].texid[0];
    const float4 texColor  = m_textures[texId]->sample(texCoordT);
    const uint   mtype     = m_materials[currMat.id].mtype;
    const uint   cflags    = m_materials[currMat.id].cflags;

    float4 fourScalarMatParams = float4(1,1,1,1);
    if(KSPEC_MAT_FOUR_TEXTURES != 0 && (cflags & FLAG_FOUR_TEXTURES) != 0)
    {
      const uint texId2  = m_materials[currMat.id].texid[2];
      const uint texId3  = m_materials[currMat.id].texid[3];

      const float2 texCoord2T = mulRows2x4(m_materials[currMat.id].row0[2], m_materials[currMat.id].row1[2], tc);
      const float2 texCoord3T = mulRows2x4(m_materials[currMat.id].row0[3], m_materials[currMat.id].row1[3], tc);

      const float4 color2 = m_textures[texId2]->sample(texCoord2T);
      const float4 color3 = m_textures[texId3]->sample(texCoord3T);
    
      if((cflags & FLAG_PACK_FOUR_PARAMS_IN_TEXTURE) != 0)
        fourScalarMatParams = color2;
      else
        fourScalarMatParams = float4(color2.x, color3.x, 1, 1);
    }

    BsdfEvalF currVal;
    {
      currVal.pdf = 0.0f;
    }
    switch(mtype)
    {
    case MAT_TYPE_BLEND:
      if(KSPEC_MAT_TYPE_BLEND != 0)
      {
        auto childMats = BlendEval(currMat, float4(), l, v, geomNormal, tc);
        currMat = childMats.first;
        needPop = false;                        // we already put 'childMats.first' in 'currMat'
        if(top + 1 <= KSPEC_BLEND_STACK_SIZE)
        {
          material_stack[top] = childMats.second; // remember second mat in stack
          top++;
        }
      }
      break;
    case MAT_TYPE_DIFFUSE:
      if(KSPEC_MAT_TYPE_DIFFUSE != 0)
      {
        FourierSpec reflSpec;// SampleMatColorSpectrumTexture(currMatId, wavelengths, DIFFUSE_COLOR, 0, tc);
        uint32_t specId = m_materials[currMat.id].spdid[DIFFUSE_COLOR];
        const uint2 data  = m_spec_offset_sz[specId];
        const uint offset = data.x;
        const uint size   = data.y;
        reflSpec = FourierSpec(m_spec_values.data() + offset, size);   
        diffuseEvalF(m_materials.data() + currMat.id, reflSpec, l, v, shadeNormal, tc,  &currVal);

        res.val += currVal.val * currMat.weight * bumpCosMult;
        res.pdf += currVal.pdf * currMat.weight;
      }
      break;
    case MAT_TYPE_CONDUCTOR:
      if(KSPEC_MAT_TYPE_CONDUCTOR != 0)
      {
        const float4 texColor = {1.0, 1.0, 1.0, 1.0}; // placeholder
        const float3 alphaTex  = to_float3(texColor);
        const float2 alpha     = float2(m_materials[currMat.id].data[CONDUCTOR_ROUGH_V], m_materials[currMat.id].data[CONDUCTOR_ROUGH_U]);


        uint precomp_offset = as_uint(m_materials[currMat.id].data[CONDUCTOR_PRECOMP_OFFSET]);

        if(!trEffectivelySmooth(alpha))
        {
          conductorRoughEvalF(m_materials.data() + currMat.id, l, v, shadeNormal, tc, alphaTex, m_precomp_conductor.data() + precomp_offset, &currVal);
        }

        res.val += currVal.val * currMat.weight * bumpCosMult;
        res.pdf += currVal.pdf * currMat.weight;
      }
      /*
      case MAT_TYPE_NEURAL_BRDF:
      if(KSPEC_MAT_TYPE_NEURAL_BRDF != 0) //TODO
      {
        uint weights_offset = m_neural_weights_offsets[currMat.id];


        const uint   ch1texId     = m_materials[currMat.id].texid[1];
        const uint   ch2texId     = m_materials[currMat.id].texid[2];
        const uint   ch3texId     = m_materials[currMat.id].texid[3];
        const float4 ch1 = m_textures[ch1texId]->sample(texCoordT);
        const float4 ch2 = m_textures[ch2texId]->sample(texCoordT);
        const float4 ch3 = m_textures[ch3texId]->sample(texCoordT);

        float buf[16]{texColor.x, texColor.y, texColor.z, texColor.w, ch1.x, ch1.y, ch1.z, ch1.w, ch2.x, ch2.y, ch2.z, ch2.w, ch3.x, ch3.y, ch3.z, ch3.w};

        neuralBrdfEval(m_materials.data() + currMat.id, m_neural_weights.data() + weights_offset, v, l, n, buf, &res);

      }
      break;*/
      default:
        break;
    }

  } while(KSPEC_MAT_TYPE_BLEND != 0 && top > 0);

  return res;
}

BsdfSampleF Integrator::MaterialSampleAndEvalF(uint a_materialId, uint tid, uint bounce, RandomGen* a_gen, float3 v, float3 n, float3 tan, float2 tc, 
                                           MisData* a_misPrev, const uint a_currRayFlags)
{
  BsdfSampleF res;
  {
    res.pdf   = 1.0f;
    res.dir   = float3(0,1,0);
    res.ior   = 1.0f;
    res.flags = a_currRayFlags;
    res.ior   = 1.0f;
  }

  uint32_t currMatId = a_materialId;
  uint     mtype     = m_materials[currMatId].mtype;
  uint     layer     = 0;
  while(KSPEC_MAT_TYPE_BLEND != 0 && mtype == MAT_TYPE_BLEND)
  {
    currMatId = BlendSampleAndEvalF(currMatId, tid, bounce, layer, a_gen, v, n, tc, a_misPrev, &res);
    mtype     = m_materials[currMatId].mtype;
    layer++;
  }
  
  // BSDF is multiplied (outside) by cosThetaOut1.
  // When normal map is enables this becames wrong because normal is changed;
  // First : return cosThetaOut in sam;
  // Second: apply cos(theta2)/cos(theta1) to cos(theta1) to get cos(theta2)
  //
  const uint normalMapId   = m_materials[currMatId].texid[1];
  const float3 geomNormal  = n;
        float3 shadeNormal = n;

  if(KSPEC_BUMP_MAPPING != 0 && normalMapId != 0xFFFFFFFF)
    shadeNormal = BumpMapping(normalMapId, currMatId, geomNormal, tan, tc);

  //const float2 texCoordT = mulRows2x4(m_materials[currMatId].row0[0], m_materials[currMatId].row1[0], tc);
  //const uint   texId     = m_materials[currMatId].texid[0];
  //const float4 texColor  = m_textures[texId]->sample(texCoordT);
  const float4 rands     = GetRandomNumbersMats(tid, a_gen, int(bounce));
  const uint cflags      = m_materials[currMatId].cflags;
  RecordMatRndNeeded(bounce, rands);


  switch(mtype)
  {    
  case MAT_TYPE_DIFFUSE:
    if(KSPEC_MAT_TYPE_DIFFUSE != 0)
    {
      FourierSpec reflSpec;// SampleMatColorSpectrumTexture(currMatId, wavelengths, DIFFUSE_COLOR, 0, tc);
      uint32_t specId = m_materials[currMatId].spdid[DIFFUSE_COLOR];
      const uint2 data  = m_spec_offset_sz[specId];
      const uint offset = data.x;
      const uint size   = data.y;
      reflSpec = FourierSpec(m_spec_values.data() + offset, size);

      diffuseSampleAndEvalF(m_materials.data() + currMatId, reflSpec, rands, v, shadeNormal, tc, &res);
    }    
    break;
  case MAT_TYPE_CONDUCTOR:
    if(KSPEC_MAT_TYPE_CONDUCTOR != 0)
    {
      const float4 texColor = {1.0, 1.0, 1.0, 1.0}; // TO-DO
      const float3 alphaTex = to_float3(texColor);    
      const float2 alpha    = float2(m_materials[currMatId].data[CONDUCTOR_ROUGH_V], m_materials[currMatId].data[CONDUCTOR_ROUGH_U]);

      uint precomp_offset = as_uint(m_materials[currMatId].data[CONDUCTOR_PRECOMP_OFFSET]);
      if(trEffectivelySmooth(alpha))
        conductorSmoothSampleAndEvalF(m_materials.data() + currMatId, rands, v, shadeNormal, tc, m_precomp_conductor.data() + precomp_offset, &res);
      else
        conductorRoughSampleAndEvalF(m_materials.data() + currMatId, rands, v, shadeNormal, tc, alphaTex, m_precomp_conductor.data() + precomp_offset, &res);
    }
    break;
  case MAT_TYPE_THIN_FILM:
    if(KSPEC_MAT_TYPE_THIN_FILM != 0)
    {
      const float4 texColor = {1.0, 1.0, 1.0, 1.0}; // TO-DO
      const float3 alphaTex = to_float3(texColor);  
      const float2 alpha = float2(m_materials[currMatId].data[FILM_ROUGH_V], m_materials[currMatId].data[FILM_ROUGH_U]);

      uint t_offset = as_uint(m_materials[currMatId].data[FILM_THICKNESS_OFFSET]);
      uint layers = as_uint(m_materials[currMatId].data[FILM_LAYERS_COUNT]);

      float extIOR = m_materials[currMatId].data[FILM_ETA_EXT];
      complex intIOR = complex(1, 0); // TO-DO

      uint precomp_offset = as_uint(m_materials[currMatId].data[FILM_PRECOMP_OFFSET]);
      if(trEffectivelySmooth(alpha))
        filmSmoothSampleAndEvalF(m_materials.data() + currMatId, extIOR, intIOR, a_misPrev->ior, rands, v, n, tc, &res,
                          m_precomp_thin_films.data() + precomp_offset);
      //else
      //  filmRoughSampleAndEval(m_materials.data() + currMatId, extIOR, filmIOR, intIOR, thickness, wavelengths_spec, a_misPrev->ior, rands, v, n, tc, alphaTex, &res,
      //                    m_precomp_thin_films.data() + precomp_offset, spectral_mode, precomp_flag);

      //res.flags |= (specId < 0xFFFFFFFF) ? RAY_FLAG_WAVES_DIVERGED : 0;
      res.flags |= RAY_FLAG_WAVES_DIVERGED;

      a_misPrev->ior = res.ior;
    }
    /*case MAT_TYPE_NEURAL_BRDF:
    if(KSPEC_MAT_TYPE_NEURAL_BRDF != 0)
    {
      uint weights_offset = m_neural_weights_offsets[currMatId];

      //const uint   ch1texId     = m_materials[a_materialId].texid[1];
      //const uint   ch2texId     = m_materials[a_materialId].texid[2];
      //const uint   ch3texId     = m_materials[a_materialId].texid[3];
      const float4 ch1 = float4(0.0f, 0.0f, 0.0f, 0.0f);//m_textures[ch1texId]->sample(texCoordT);
      const float4 ch2 = float4(0.0f, 0.0f, 0.0f, 0.0f);//m_textures[ch2texId]->sample(texCoordT);
      const float4 ch3 = float4(0.0f, 0.0f, 0.0f, 0.0f);//m_textures[ch3texId]->sample(texCoordT);

      float buf[16]{texColor.x, texColor.y, texColor.z, texColor.w, ch1.x, ch1.y, ch1.z, ch1.w, ch2.x, ch2.y, ch2.z, ch2.w, ch3.x, ch3.y, ch3.z, ch3.w};

      if(m_spectral_mode == 0)
        neuralBrdfSampleAndEval(m_materials.data() + currMatId, m_neural_weights.data() + weights_offset, rands, v, n, buf, &res);
      else
        neuralSpecSmoothSampleAndEval(m_materials.data() + currMatId, m_neural_weights.data() + weights_offset, v, n, wavelengths, buf, &res);

    }
    break;*/
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
                                    
uint32_t Integrator::BlendSampleAndEvalF(uint a_materialId, uint tid, uint bounce, uint layer, RandomGen* a_gen, float3 v, float3 n, float2 tc, 
                          MisData* a_misPrev, BsdfSampleF* a_pRes)
{
  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const uint   texId     = m_materials[a_materialId].texid[0];
  const float4 weightDat = m_textures[texId]->sample(texCoordT);
  const float  weightTex = weightDat.x;
  const float  weight    = m_materials[a_materialId].data[BLEND_WEIGHT] * weightTex;

  const uint matId1 = m_materials[a_materialId].datai[0];
  const uint matId2 = m_materials[a_materialId].datai[1];

  uint32_t selectedMatId = matId1;
  const float select = GetRandomNumbersMatB(tid, a_gen, int(bounce), int(layer));
  RecordBlendRndNeeded(bounce, layer, select);

  if(select < weight)
  {
    a_pRes->pdf *= weight;
    a_pRes->val *= weight;
    selectedMatId = matId2;
  }
  else
  {
    a_pRes->pdf *= 1.0f - weight;
    a_pRes->val *= 1.0f - weight;
    selectedMatId = matId1;
  }

  return selectedMatId;
}






BsdfSampleN Integrator::MaterialSampleAndEvalN(uint a_materialId, uint tid, uint bounce, const SpecN *wavelengths, RandomGen* a_gen, float3 v, float3 n, float3 tan, float2 tc, 
                                             MisData* a_misPrev, const uint a_currRayFlags)
{
  BsdfSampleN res;
  {
    res.val   = SpecN(0);
    res.pdf   = 1.0f;
    res.dir   = float3(0,1,0);
    res.ior   = 1.0f;
    res.flags = a_currRayFlags;
  }

  const SpecN &lambda = *wavelengths;

  uint32_t currMatId = a_materialId;
  uint     mtype     = m_materials[currMatId].mtype;
  uint     layer     = 0;
  while(KSPEC_MAT_TYPE_BLEND != 0 && mtype == MAT_TYPE_BLEND)
  {
    currMatId = BlendSampleAndEvalN(currMatId, tid, bounce, layer, lambda, a_gen, v, n, tc, a_misPrev, &res);
    mtype     = m_materials[currMatId].mtype;
    layer++;
  }
  
  // BSDF is multiplied (outside) by cosThetaOut1.
  // When normal map is enables this becames wrong because normal is changed;
  // First : return cosThetaOut in sam;
  // Second: apply cos(theta2)/cos(theta1) to cos(theta1) to get cos(theta2)
  //
  const uint normalMapId   = m_materials[currMatId].texid[1];
  const float3 geomNormal  = n;
        float3 shadeNormal = n;

  if(KSPEC_BUMP_MAPPING != 0 && normalMapId != 0xFFFFFFFF)
    shadeNormal = BumpMapping(normalMapId, currMatId, geomNormal, tan, tc);

  const float4 rands     = GetRandomNumbersMats(tid, a_gen, int(bounce));
  const uint cflags      = m_materials[currMatId].cflags;
  RecordMatRndNeeded(bounce, rands);

  SpecN spec = SampleMatParamSpectrumN(currMatId, lambda, 0);

  switch(mtype)
  {
    case MAT_TYPE_CONDUCTOR:
    if(KSPEC_MAT_TYPE_CONDUCTOR != 0)
    {
      const float3 alphaTex = {1.0f, 1.0f, 1.0f}; 
      const float2 alpha    = float2(m_materials[currMatId].data[CONDUCTOR_ROUGH_V], m_materials[currMatId].data[CONDUCTOR_ROUGH_U]);
      const SpecN etaSpec  = SampleMatParamSpectrumN(currMatId, lambda, 0);
      const SpecN kSpec    = SampleMatParamSpectrumN(currMatId, lambda, 1);
      if(trEffectivelySmooth(alpha))
        conductorSmoothSampleAndEvalN(m_materials.data() + currMatId, etaSpec, kSpec, rands, v, shadeNormal, tc, &res);
      else
        conductorRoughSampleAndEvalN(m_materials.data() + currMatId, etaSpec, kSpec, rands, v, shadeNormal, tc, alphaTex, &res);
    }
    break;/*
    case MAT_TYPE_THIN_FILM:
    if(KSPEC_MAT_TYPE_THIN_FILM != 0)
    {
      const float3 alphaTex = to_float3(texColor);  
      const float2 alpha = float2(m_materials[currMatId].data[FILM_ROUGH_V], m_materials[currMatId].data[FILM_ROUGH_U]);

      uint t_offset = as_uint(m_materials[currMatId].data[FILM_THICKNESS_OFFSET]);
      uint layers = as_uint(m_materials[currMatId].data[FILM_LAYERS_COUNT]);
      const bool spectral_mode = wavelengths[0] > 0.0f;
      // sampling 3 wavelengths for naive RGB method
      float4 wavelengths_spec = spectral_mode? float4(wavelengths[0], 0.0f, 0.0f, 0.0f) : float4(645.f, 525.f, 445.f, 0.0f);
      float4 wavelengths_sample = spectral_mode? float4(wavelengths[0], 0.0f, 0.0f, 0.0f) : float4(525.f, 0.f, 0.f, 0.0f);
      float extIOR = m_materials[currMatId].data[FILM_ETA_EXT];
      complex intIOR = complex(
        SampleFilmsSpectrum(currMatId, wavelengths_sample, FILM_ETA_OFFSET, FILM_ETA_SPECID_OFFSET, layers - 1)[0],
        SampleFilmsSpectrum(currMatId, wavelengths_sample, FILM_K_OFFSET, FILM_K_SPECID_OFFSET, layers - 1)[0]
      );
      complex filmIOR = complex(
        SampleFilmsSpectrum(currMatId, wavelengths, FILM_ETA_OFFSET, FILM_ETA_SPECID_OFFSET, 0)[0],
        SampleFilmsSpectrum(currMatId, wavelengths, FILM_K_OFFSET, FILM_K_SPECID_OFFSET, 0)[0]
      );

      float thickness;
      if (as_uint(m_materials[currMatId].data[FILM_THICKNESS_MAP]) > 0u)
      {
        const uint texId  = m_materials[currMatId].texid[2];
        const float2 texCoord = mulRows2x4(m_materials[currMatId].row0[2], m_materials[currMatId].row1[2], tc);
        const float4 thickness_val = m_textures[texId]->sample(texCoord);
        float thickness_max = m_materials[currMatId].data[FILM_THICKNESS_MAX];
        float thickness_min = m_materials[currMatId].data[FILM_THICKNESS_MIN];
        thickness = (thickness_max - thickness_min) * thickness_val.x + thickness_min;
        //std::cout << fourScalarMatParams.x << " " << fourScalarMatParams.y << " " << fourScalarMatParams.z << " " << fourScalarMatParams.w << std::endl;
      }
      else
      {
        thickness = m_materials[currMatId].data[FILM_THICKNESS];
      }

      bool precomp_flag = as_uint(m_materials[currMatId].data[FILM_PRECOMP_FLAG]) > 0u;

      uint precomp_offset = precomp_flag ? as_uint(m_materials[currMatId].data[FILM_PRECOMP_OFFSET]) : 0;
      if(trEffectivelySmooth(alpha))
        filmSmoothSampleAndEval(m_materials.data() + currMatId, extIOR, filmIOR, intIOR, thickness, wavelengths_spec, a_misPrev->ior, rands, v, n, tc, &res,
                          m_precomp_thin_films.data() + precomp_offset, spectral_mode, precomp_flag);
      else
        filmRoughSampleAndEval(m_materials.data() + currMatId, extIOR, filmIOR, intIOR, thickness, wavelengths_spec, a_misPrev->ior, rands, v, n, tc, alphaTex, &res,
                          m_precomp_thin_films.data() + precomp_offset, spectral_mode, precomp_flag);

      //res.flags |= (specId < 0xFFFFFFFF) ? RAY_FLAG_WAVES_DIVERGED : 0;
      res.flags |= RAY_FLAG_WAVES_DIVERGED;

      a_misPrev->ior = res.ior;
    }
    break;*/
    case MAT_TYPE_DIFFUSE:
    if(KSPEC_MAT_TYPE_DIFFUSE != 0)
    {      
      diffuseSampleAndEvalN(m_materials.data() + currMatId, spec, rands, v, shadeNormal, tc, &res);
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

BsdfEvalN Integrator::MaterialEvalN(uint a_materialId, const SpecN *wavelengths, float3 l, float3 v, float3 n, float3 tan, float2 tc)
{
  BsdfEvalN res;
  {
    res.val = SpecN(0);
    res.pdf   = 0.0f;
  }

  const SpecN &lambda = *wavelengths;

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
    const uint normalMapId = m_materials[currMat.id].texid[1];
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

    SpecN spec = SampleMatParamSpectrumN(currMat.id, lambda, 0);

    const SpecN &texColor  = spec;

    const uint   mtype     = m_materials[currMat.id].mtype;
    const uint   cflags    = m_materials[currMat.id].cflags;


    BsdfEvalN currVal;
    {
      currVal.val = SpecN(0);
      currVal.pdf   = 0.0f;
    }


    switch(mtype)
    {
      case MAT_TYPE_CONDUCTOR: 
      if(KSPEC_MAT_TYPE_CONDUCTOR != 0)
      {
        const float3 alphaTex  = {1.0f, 1.0f, 1.0f};
        const float2 alpha     = float2(m_materials[currMat.id].data[CONDUCTOR_ROUGH_V], m_materials[currMat.id].data[CONDUCTOR_ROUGH_U]);

        if(!trEffectivelySmooth(alpha))
        {
          const SpecN etaSpec = SampleMatParamSpectrumN(currMat.id, lambda, 0);
          const SpecN kSpec   = SampleMatParamSpectrumN(currMat.id, lambda, 1);
          conductorRoughEvalN(m_materials.data() + currMat.id, etaSpec, kSpec, l, v, shadeNormal, tc, alphaTex, &currVal);
        }

        res.val += currVal.val * currMat.weight * bumpCosMult;
        res.pdf += currVal.pdf * currMat.weight;
      }
      break;
      /*
      case MAT_TYPE_THIN_FILM: 
      if(KSPEC_MAT_TYPE_THIN_FILM != 0)
      {
        const float3 alphaTex = to_float3(texColor);  
        const float2 alpha = float2(m_materials[currMat.id].data[FILM_ROUGH_V], m_materials[currMat.id].data[FILM_ROUGH_U]);

        if(!trEffectivelySmooth(alpha))
        {
          uint t_offset = as_uint(m_materials[currMat.id].data[FILM_THICKNESS_OFFSET]);
          uint layers = as_uint(m_materials[currMat.id].data[FILM_LAYERS_COUNT]);
          const bool spectral_mode = lambda[0] > 0.0f;
          // sampling 3 wavelengths for naive RGB method
          float4 wavelengths_spec = spectral_mode? float4(lambda[0], 0.0f, 0.0f, 0.0f) : float4(700.f, 525.f, 450.f, 0.0f);
          float4 wavelengths_sample = spectral_mode? float4(lambda[0], 0.0f, 0.0f, 0.0f) : float4(525.f, 0.f, 0.f, 0.0f);
          float extIOR = m_materials[currMat.id].data[FILM_ETA_EXT];
          complex intIOR = complex(
            SampleFilmsSpectrum(currMat.id, wavelengths_sample, FILM_ETA_OFFSET, FILM_ETA_SPECID_OFFSET, layers - 1)[0],
            SampleFilmsSpectrum(currMat.id, wavelengths_sample, FILM_K_OFFSET, FILM_K_SPECID_OFFSET, layers - 1)[0]
          );      
          complex filmIOR = complex(
            SampleFilmsSpectrum(currMat.id, lambda, FILM_ETA_OFFSET, FILM_ETA_SPECID_OFFSET, 0)[0],
            SampleFilmsSpectrum(currMat.id, lambda, FILM_K_OFFSET, FILM_K_SPECID_OFFSET, 0)[0]
          );

          float thickness;
          if (as_uint(m_materials[currMat.id].data[FILM_THICKNESS_MAP]) > 0u)
          {
            const uint texId  = m_materials[currMat.id].texid[2];
            const float2 texCoord = mulRows2x4(m_materials[currMat.id].row0[2], m_materials[currMat.id].row1[2], tc);
            const float4 thickness_val = m_textures[texId]->sample(texCoord);
            float thickness_max = m_materials[currMat.id].data[FILM_THICKNESS_MAX];
            float thickness_min = m_materials[currMat.id].data[FILM_THICKNESS_MIN];
            thickness = (thickness_max - thickness_min) * thickness_val.x + thickness_min;
          }
          else
          {
            thickness = m_materials[currMat.id].data[FILM_THICKNESS];
          }

          bool precomp_flag = as_uint(m_materials[currMat.id].data[FILM_PRECOMP_FLAG]) > 0u;
          uint precomp_offset = precomp_flag ? as_uint(m_materials[currMat.id].data[FILM_PRECOMP_OFFSET]) : 0;
          filmRoughEval(m_materials.data() + currMat.id, extIOR, filmIOR, intIOR, thickness, wavelengths_spec, l, v, n, tc, alphaTex, &currVal,
                          m_precomp_thin_films.data() + precomp_offset, spectral_mode, precomp_flag);
        }

        res.val += currVal.val * currMat.weight * bumpCosMult;
        res.pdf += currVal.pdf * currMat.weight;
      }
      break;*/
      case MAT_TYPE_DIFFUSE:
      if(KSPEC_MAT_TYPE_DIFFUSE != 0)
      {
          
        diffuseEvalN(m_materials.data() + currMat.id, spec, l, v, shadeNormal, tc,  &currVal);

        res.val += currVal.val * currMat.weight * bumpCosMult;
        res.pdf += currVal.pdf * currMat.weight;
      }
      break;
      case MAT_TYPE_BLEND:
      if(KSPEC_MAT_TYPE_BLEND != 0)
      {
        auto childMats = BlendEval(currMat, float4(0.0f), l, v, geomNormal, tc);
        currMat = childMats.first;
        needPop = false;                        // we already put 'childMats.first' in 'currMat'
        if(top + 1 <= KSPEC_BLEND_STACK_SIZE)
        {
          material_stack[top] = childMats.second; // remember second mat in stack
          top++;
        }
      }
      break;
      default:
        break;
    }

  } while(KSPEC_MAT_TYPE_BLEND != 0 && top > 0);

  return res;
}

uint32_t Integrator::BlendSampleAndEvalN(uint a_materialId, uint tid, uint bounce, uint layer, const SpecN &wavelengths, RandomGen* a_gen, float3 v, float3 n, float2 tc, 
                                        MisData* a_misPrev, BsdfSampleN* a_pRes)
{
  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const uint   texId     = m_materials[a_materialId].texid[0];
  const float4 weightDat = m_textures[texId]->sample(texCoordT);
  const float  weightTex = weightDat.x;
  const float  weight    = m_materials[a_materialId].data[BLEND_WEIGHT] * weightTex;

  const uint matId1 = m_materials[a_materialId].datai[0];
  const uint matId2 = m_materials[a_materialId].datai[1];

  uint32_t selectedMatId = matId1;
  const float select = GetRandomNumbersMatB(tid, a_gen, int(bounce), int(layer));
  RecordBlendRndNeeded(bounce, layer, select);

  if(select < weight)
  {
    a_pRes->pdf *= weight;
    a_pRes->val *= weight;
    selectedMatId = matId2;
  }
  else
  {
    a_pRes->pdf *= 1.0f - weight;
    a_pRes->val *= 1.0f - weight;
    selectedMatId = matId1;
  }

  return selectedMatId;
}
