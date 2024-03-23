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
#include <omp.h>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

void IntegratorDR::LoadSceneEnd()
{
  m_texAddressTable.resize(m_textures.size());
  for(auto& texInfo : m_texAddressTable)
    texInfo = {size_t(-1),0,0,0,0};

  m_gradSize = 0;
}

std::pair<size_t, size_t> IntegratorDR::PutDiffTex2D(uint32_t texId, uint32_t width, uint32_t height, uint32_t channels)
{
  if(texId >= m_texAddressTable.size())
  {
    std::cout << "[IntegratorDR::PutDiffTex2D]: bad tex id = " << texId << std::endl;
    return std::make_pair(size_t(-1), 0);
  }

  m_texAddressTable[texId].offset   = m_gradSize;
  m_texAddressTable[texId].width    = width;
  m_texAddressTable[texId].height   = height;
  m_texAddressTable[texId].channels = channels;
  m_texAddressTable[texId].fwidth   = float(width);
  m_texAddressTable[texId].fheight  = float(height);

  size_t oldOffset = m_gradSize;
  size_t currSize  = size_t(width)*size_t(height)*size_t(channels);
  
  m_gradSize += currSize;
  return std::make_pair(oldOffset, currSize);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

uint32_t IntegratorDR::BlendSampleAndEval(uint a_materialId, uint bounce, uint layer, float4 wavelengths, float3 v, float3 n, float2 tc, 
                                          MisData* a_misPrev, BsdfSample* a_pRes, const float* drands, const float* dparams)
{
  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const uint   texId     = m_materials[a_materialId].texid[0];
  const float4 weightDat = Tex2DFetchAD(texId, texCoordT, dparams);
  const float  weightTex = weightDat.x;
  const float  weight    = m_materials[a_materialId].data[BLEND_WEIGHT] * weightTex;

  const uint matId1 = m_materials[a_materialId].datai[0];
  const uint matId2 = m_materials[a_materialId].datai[1];

  uint32_t selectedMatId = matId1;
  const float select     = drands[bounce*RND_PER_BOUNCE + RND_BLD_ID + layer];
  
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

MatIdWeightPair IntegratorDR::BlendEval(MatIdWeight a_mat, float4 wavelengths, float3 l, float3 v, float3 n, float2 tc, const float* dparams)
{
  const float2 texCoordT = mulRows2x4(m_materials[a_mat.id].row0[0], m_materials[a_mat.id].row1[0], tc);
  const uint   texId     = m_materials[a_mat.id].texid[0];
  const float4 weightDat = Tex2DFetchAD(texId, texCoordT, dparams); 
  
  const float  weightTex = weightDat.x;
  const float  weight    = m_materials[a_mat.id].data[BLEND_WEIGHT] * weightTex;

  const uint matId1      = m_materials[a_mat.id].datai[0];
  const uint matId2      = m_materials[a_mat.id].datai[1];

  MatIdWeight p1, p2;
  p1.id     = matId1;
  p1.weight = a_mat.weight * (1.0f - weight);
  p2.id     = matId2;
  p2.weight = a_mat.weight * weight;

  return make_weight_pair(p1, p2);
}

static inline float3 NormalMapTransform(const uint materialFlags, float3 normalFromTex)
{
  float3 normalTS = make_float3(2.0f * normalFromTex.x - 1.0f, 2.0f * normalFromTex.y - 1.0f, normalFromTex.z);

  if((materialFlags & FLAG_NMAP_INVERT_X) != 0)
    normalTS.x *= (-1.0f);

  if((materialFlags & FLAG_NMAP_INVERT_Y) != 0)
    normalTS.y *= (-1.0f);

  if((materialFlags & FLAG_NMAP_SWAP_XY) != 0)
  {
    float temp = normalTS.x;
    normalTS.x = normalTS.y;
    normalTS.y = temp;
  }

  return normalTS; // normalize(normalTS); // do we nedd this normalize here?
}

float3 IntegratorDR::BumpMapping(uint normalMapId, uint currMatId, float3 n, float3 tan, float2 tc, const float* dparams)
{
  const uint   mflags    = m_materials[currMatId].cflags;
  const float2 texCoordT = mulRows2x4(m_materials[currMatId].row0[1], m_materials[currMatId].row1[1], tc);
  const float4 normalTex = Tex2DFetchAD(normalMapId, texCoordT, dparams);
  const float3 normalTS  = NormalMapTransform(mflags, to_float3(normalTex));
  
  const float3   bitan = cross(n, tan);
  const float3x3 tangentTransform = make_float3x3(tan, bitan, n);

  return normalize(inverse3x3(tangentTransform)*normalTS);
}

BsdfSample IntegratorDR::MaterialSampleAndEval(uint a_materialId, uint bounce, float4 wavelengths, float3 v, float3 n, float3 tan, float2 tc, 
                                               MisData* a_misPrev, const uint a_currRayFlags, const float* drands, const float* dparams)
{
  BsdfSample res;
  {
    res.val   = float4(0, 0, 0, 0);
    res.pdf   = 1.0f;
    res.dir   = float3(0,1,0);
    res.flags = a_currRayFlags;
    res.ior   = 1.0f;
  }

  uint32_t currMatId = a_materialId;
  
  //const float2 texCoordT = mulRows2x4(m_materials[currMatId].row0[0], m_materials[currMatId].row1[0], tc);
  //const uint   texId     = m_materials[currMatId].texid[0];
  //const float4 texColor  = Tex2DFetchAD(texId, texCoordT, dparams);
  //const float4 rands = float4(drands[bounce*RND_PER_BOUNCE + RND_MTL_ID + 0], drands[bounce*RND_PER_BOUNCE + RND_MTL_ID + 1],
  //                            drands[bounce*RND_PER_BOUNCE + RND_MTL_ID + 2], drands[bounce*RND_PER_BOUNCE + RND_MTL_ID + 3]);

  //const float4 color = m_materials[currMatId].colors[GLTF_COLOR_BASE]*texColor;
  
  //const float3 lambertDir   = lambertSample(float2(rands.x, rands.y), v, n);
  //const float  lambertPdf   = lambertEvalPDF(lambertDir, v, n);
  //const float  lambertVal   = lambertEvalBSDF(lambertDir, v, n);
  //
  //res.val   = lambertVal*color;
  //res.dir   = lambertDir;
  //res.pdf   = lambertPdf;
  //res.flags = RAY_FLAG_HAS_NON_SPEC;
  
  uint     mtype     = m_materials[currMatId].mtype;
  uint     layer     = 0;
  while(KSPEC_MAT_TYPE_BLEND != 0 && mtype == MAT_TYPE_BLEND)
  {
    currMatId = BlendSampleAndEval(currMatId, bounce, layer, wavelengths, v, n, tc, a_misPrev, &res, drands, dparams);
    mtype     = m_materials[currMatId].mtype;
    layer++;
  }
  
  // BSDF is multiplied (outside) by cosThetaOut1.
  // When normal map is enables this becames wrong because normal is changed;
  // First : return cosThetaOut in sam;
  // Second: apply cos(theta2)/cos(theta1) to cos(theta1) to get cos(theta2)
  //
  const uint cflags = m_materials[currMatId].cflags;
  float4 fourScalarMatParams = float4(1,1,1,1);
  if(KSPEC_MAT_FOUR_TEXTURES != 0 && (cflags & FLAG_FOUR_TEXTURES) != 0)
  {
    const uint texId2  = m_materials[currMatId].texid[2];
    const uint texId3  = m_materials[currMatId].texid[3];
    
    const float2 texCoord2T = mulRows2x4(m_materials[currMatId].row0[2], m_materials[currMatId].row1[2], tc);
    const float2 texCoord3T = mulRows2x4(m_materials[currMatId].row0[3], m_materials[currMatId].row1[3], tc);

    const float4 color2 = m_textures[texId2]->sample(texCoord2T);
    const float4 color3 = m_textures[texId3]->sample(texCoord3T);
    
    if((cflags & FLAG_PACK_FOUR_PARAMS_IN_TEXTURE) != 0)
      fourScalarMatParams = color2;
    else
      fourScalarMatParams = float4(color2.x, color3.x, 1, 1);
  }

  const uint normalMapId   = m_materials[currMatId].texid[1];
  const float3 geomNormal  = n;
        float3 shadeNormal = n;

  if(KSPEC_BUMP_MAPPING != 0 && normalMapId != 0xFFFFFFFF)
    shadeNormal = BumpMapping(normalMapId, currMatId, geomNormal, tan, tc, dparams);

  const float2 texCoordT = mulRows2x4(m_materials[currMatId].row0[0], m_materials[currMatId].row1[0], tc);
  const uint   texId     = m_materials[currMatId].texid[0];
  const float4 texColor  = Tex2DFetchAD(texId, texCoordT, dparams); 
  
  const float4 rands = float4(drands[bounce*RND_PER_BOUNCE + RND_MTL_ID + 0], drands[bounce*RND_PER_BOUNCE + RND_MTL_ID + 1],
                              drands[bounce*RND_PER_BOUNCE + RND_MTL_ID + 2], drands[bounce*RND_PER_BOUNCE + RND_MTL_ID + 3]);

  switch(mtype)
  {
    case MAT_TYPE_GLTF:
    if(KSPEC_MAT_TYPE_GLTF != 0)
    {
      const float4 color = m_materials[currMatId].colors[GLTF_COLOR_BASE]*texColor;
      gltfSampleAndEval(m_materials.data() + currMatId, rands, v, shadeNormal, tc, color, fourScalarMatParams, &res);
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
      const float3 alphaTex = to_float3(texColor);    
      const float2 alpha    = float2(m_materials[currMatId].data[CONDUCTOR_ROUGH_V], m_materials[currMatId].data[CONDUCTOR_ROUGH_U]);
      const float4 etaSpec  = SampleMatParamSpectrum(currMatId, wavelengths, CONDUCTOR_ETA, 0, dparams);
      const float4 kSpec    = SampleMatParamSpectrum(currMatId, wavelengths, CONDUCTOR_K,   1, dparams);
      if(trEffectivelySmooth(alpha))
        conductorSmoothSampleAndEval(m_materials.data() + currMatId, etaSpec, kSpec, rands, v, shadeNormal, tc, &res);
      else
        conductorRoughSampleAndEval(m_materials.data() + currMatId, etaSpec, kSpec, rands, v, shadeNormal, tc, alphaTex, &res);
    }
    break;
    case MAT_TYPE_DIFFUSE:
    if(KSPEC_MAT_TYPE_DIFFUSE != 0)
    {
      float4 reflSpec    = SampleMatColorParamSpectrum(currMatId, wavelengths, DIFFUSE_COLOR, 0, dparams);
      if(m_spectral_mode == 0)
        reflSpec *= texColor;
      diffuseSampleAndEval(m_materials.data() + currMatId, reflSpec, rands, v, shadeNormal, tc, &res);
    }
    break;
    case MAT_TYPE_PLASTIC:
    if(KSPEC_MAT_TYPE_PLASTIC != 0)
    {
      const float4 color = texColor;
      float4 reflSpec    = SampleMatColorParamSpectrum(currMatId, wavelengths, PLASTIC_COLOR, 0, dparams);
      if(m_spectral_mode == 0)
        reflSpec *= color;

      const uint precomp_id = m_materials[currMatId].datai[0];

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

BsdfEval IntegratorDR::MaterialEval(uint a_materialId, float4 wavelengths, float3 l, float3 v, float3 n, float3 tan, float2 tc, const float* dparams)
{
  BsdfEval res;
  {
    res.val = float4(0,0,0,0);
    res.pdf = 0.0f;
  }

  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const uint   texId     = m_materials[a_materialId].texid[0];
  const float4 texColor  = Tex2DFetchAD(texId, texCoordT, dparams); 
  const float4 color     = m_materials[a_materialId].colors[GLTF_COLOR_BASE] * texColor;
  
  //float lambertVal       = lambertEvalBSDF(l, v, n);
  //const float lambertPdf = lambertEvalPDF (l, v, n);    
  //
  //res.val = lambertVal*color;
  //res.pdf = lambertPdf;
  
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
      shadeNormal = BumpMapping(normalMapId, currMat.id, geomNormal, tan, tc, dparams);
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
    const float4 texColor  = Tex2DFetchAD(texId, texCoordT, dparams); 
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
        const float4 color     = (m_materials[currMat.id].colors[GLTF_COLOR_BASE]) * texColor;
        gltfEval(m_materials.data() + currMat.id, l, v, shadeNormal, tc, color, fourScalarMatParams, &currVal);

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
        const float3 alphaTex  = to_float3(texColor);
        const float2 alpha     = float2(m_materials[currMat.id].data[CONDUCTOR_ROUGH_V], m_materials[currMat.id].data[CONDUCTOR_ROUGH_U]);

        if(!trEffectivelySmooth(alpha))
        {
          const float4 etaSpec = SampleMatParamSpectrum(currMat.id, wavelengths, CONDUCTOR_ETA, 0, dparams);
          const float4 kSpec   = SampleMatParamSpectrum(currMat.id, wavelengths, CONDUCTOR_K,   1, dparams);
          conductorRoughEval(m_materials.data() + currMat.id, etaSpec, kSpec, l, v, shadeNormal, tc, alphaTex, &currVal);
        }

        res.val += currVal.val * currMat.weight * bumpCosMult;
        res.pdf += currVal.pdf * currMat.weight;
        break;
      }
      case MAT_TYPE_DIFFUSE:
      if(KSPEC_MAT_TYPE_DIFFUSE != 0)
      {
        float4 reflSpec = SampleMatColorParamSpectrum(currMat.id, wavelengths, DIFFUSE_COLOR, 0, dparams);
        if(m_spectral_mode == 0)
          reflSpec *= texColor; 
        diffuseEval(m_materials.data() + currMat.id, reflSpec, l, v, shadeNormal, tc, &currVal);

        res.val += currVal.val * currMat.weight * bumpCosMult;
        res.pdf += currVal.pdf * currMat.weight;
        break;
      }
      case MAT_TYPE_PLASTIC:
      if(KSPEC_MAT_TYPE_PLASTIC != 0)
      {
        const float4 color = texColor;
        float4 reflSpec    = SampleMatColorParamSpectrum(currMat.id, wavelengths, PLASTIC_COLOR, 0, dparams);
        if(m_spectral_mode == 0)
          reflSpec *= color;
        const uint precomp_id = m_materials[currMat.id].datai[0];
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
        auto childMats = BlendEval(currMat, wavelengths, l, v, geomNormal, tc, dparams);
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


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void IntegratorDR::kernel_InitEyeRay2(uint tid, float4* rayPosAndNear, float4* rayDirAndFar, float4* wavelengths, 
                                      float4* accumColor,    float4* accumuThoroughput,
                                      uint* rayFlags, MisData* misData, 
                                      const float* drands, const float* dparams) // 
{
  if(tid >= m_maxThreadId)
    return;

  *accumColor        = make_float4(0,0,0,0);
  *accumuThoroughput = make_float4(1,1,1,1);
  *rayFlags          = 0;
  *misData           = makeInitialMisData();

  const uint XY = m_packedXY[tid];
  const uint x = (XY & 0x0000FFFF);
  const uint y = (XY & 0xFFFF0000) >> 16;

  const float2 pixelOffsets = float2(drands[RANDS_ARRAY_SIZE - LENS_RANDS + 0], drands[RANDS_ARRAY_SIZE - LENS_RANDS + 1]);

  float3 rayDir = EyeRayDirNormalized((float(x) + pixelOffsets.x)/float(m_winWidth), 
                                      (float(y) + pixelOffsets.y)/float(m_winHeight), m_projInv);
  float3 rayPos = float3(0,0,0);

  transform_ray3f(m_worldViewInv, &rayPos, &rayDir);
  
  if(KSPEC_SPECTRAL_RENDERING !=0 && m_spectral_mode != 0)
  {
    float u = drands[RANDS_ARRAY_SIZE - LENS_RANDS + 3];
    *wavelengths = SampleWavelengths(u, LAMBDA_MIN, LAMBDA_MAX);
  }
  else
  {
    const uint32_t sample_sz = sizeof((*wavelengths).M) / sizeof((*wavelengths).M[0]);
    for (uint32_t i = 0; i < sample_sz; ++i) 
      (*wavelengths)[i] = 0.0f;
  }
 
  *rayPosAndNear = to_float4(rayPos, 0.0f);
  *rayDirAndFar  = to_float4(rayDir, FLT_MAX);
}

void IntegratorDR::kernel_RayTrace2(uint tid, uint bounce, uint cpuThreadId, const float4* rayPosAndNear, const float4* rayDirAndFar,
                                    float4* out_hit1, float4* out_hit2, float4* out_hit3, uint* out_instId, uint* rayFlags,
                                    const float* dparams)
{
  if(tid >= m_maxThreadId)
    return;
  uint currRayFlags = *rayFlags;
  if(isDeadRay(currRayFlags))
    return;

  const float4 rayPos = *rayPosAndNear;
  const float4 rayDir = *rayDirAndFar ;
  const CRT_Hit hit   = m_recorded[cpuThreadId].perBounce[bounce].hit;

  if(hit.geomId != uint32_t(-1))
  {
    const float2 uv     = float2(hit.coords[0], hit.coords[1]);
    const float3 hitPos = to_float3(rayPos) + (hit.t*0.999999f)*to_float3(rayDir); // set hit slightlyt closer to old ray origin to prevent self-interseaction and e.t.c bugs
    
    // alternative, you may consider Johannes Hanika solution from  Ray Tracing Gems2  
    /////////////////////////////////////////////////////////////////////////////////
    // // get distance vectors from triangle vertices
    // vec3 tmpu = P - A, tmpv = P - B, tmpw = P - C
    // // project these onto the tangent planes
    // // defined by the shading normals
    // float dotu = min (0.0, dot(tmpu , nA))
    // float dotv = min (0.0, dot(tmpv , nB))
    // float dotw = min (0.0, dot(tmpw , nC))
    // tmpu -= dotu*nA
    // tmpv -= dotv*nB
    // tmpw -= dotw*nC
    // // finally P' is the barycentric mean of these three
    // vec3 Pp = P + u*tmpu + v*tmpv + w*tmpw
    /////////////////////////////////////////////////////////////////////////////////

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

    // transform surface point with matrix and flip normal if needed
    //
    hitNorm                = normalize(mul3x3(m_normMatrices[hit.instId], hitNorm));
    hitTang                = normalize(mul3x3(m_normMatrices[hit.instId], hitTang));
    const float flipNorm   = dot(to_float3(rayDir), hitNorm) > 0.001f ? -1.0f : 1.0f; // beware of transparent materials which use normal sign to identity "inside/outside" glass for example
    hitNorm                = flipNorm * hitNorm;
    hitTang                = flipNorm * hitTang; // do we need this ??
    
    if (flipNorm < 0.0f) currRayFlags |=  RAY_FLAG_HAS_INV_NORMAL;
    else                 currRayFlags &= ~RAY_FLAG_HAS_INV_NORMAL;

    const uint midOriginal = m_matIdByPrimId[m_matIdOffsets[hit.geomId] + hit.primId];
    const uint midRemaped  = RemapMaterialId(midOriginal, hit.instId, dparams);

    *rayFlags              = packMatId(currRayFlags, midRemaped);
    *out_hit1              = to_float4(hitPos,  hitTexCoord.x); 
    *out_hit2              = to_float4(hitNorm, hitTexCoord.y);
    *out_hit3              = to_float4(hitTang, 0.0f);
    *out_instId            = hit.instId;
  }
  else
    *rayFlags              = currRayFlags | (RAY_FLAG_IS_DEAD | RAY_FLAG_OUT_OF_SCENE);
}

void IntegratorDR::kernel_SampleLightSource(uint tid, uint cpuThreadId, const float4* rayPosAndNear, const float4* rayDirAndFar, 
                                            const float4* wavelengths, const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_hitPart3,
                                            const uint* rayFlags, uint bounce, float4* out_shadeColor, const float* drands, const float* dparams)
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

  const float2 rands = float2(drands[bounce*RND_PER_BOUNCE + RND_LTG_ID + 0], drands[bounce*RND_PER_BOUNCE + RND_LTG_ID + 1]);
  const int lightId  = m_recorded[cpuThreadId].perBounceLightId[bounce];

  if(lightId < 0) // no lights or invalid light id
  {
    *out_shadeColor = float4(0.0f, 0.0f, 0.0f, 0.0f);
    return;
  }
  
  const LightSample lSam = LightSampleRev(lightId, rands, hit.pos, dparams);
  const float  hitDist   = std::sqrt(dot(hit.pos - lSam.pos, hit.pos - lSam.pos));
  
  const float3 shadowRayDir = normalize(lSam.pos - hit.pos);                             //
  const float3 shadowRayPos = hit.pos + hit.norm*std::max(maxcomp(hit.pos), 1.0f)*5e-6f; // 
 
  const bool   inShadow     = (m_recorded[cpuThreadId].perBounce[bounce].inShadow == 1);
  const bool   inIllumArea  = (dot(shadowRayDir, lSam.norm) < 0.0f) || lSam.isOmni;
  
  if(!inShadow && inIllumArea) 
  {
    const BsdfEval bsdfV    = MaterialEval(matId, lambda, shadowRayDir, (-1.0f)*ray_dir, hit.norm, hit.tang, hit.uv, dparams);
    //const BsdfEval bsdfV = {float4(1.0f), 1.0f};
    const float cosThetaOut = std::max(dot(shadowRayDir, hit.norm), 0.0f);
    
    float      lgtPdfW      = LightPdfSelectRev(lightId) * LightEvalPDF(lightId, shadowRayPos, shadowRayDir, lSam.pos, lSam.norm, dparams);
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

    if(m_renderLayer >= 1 && int(bounce) < int(m_renderLayer)-1) // skip some number of bounces if this is set
      misWeight = 0.0f;
    
    const float4 lightColor = LightIntensity(lightId, wavelengths, shadowRayDir, dparams);
    *out_shadeColor = (lightColor * bsdfV.val / lgtPdfW) * cosThetaOut * misWeight;
  }
  else
    *out_shadeColor = float4(0.0f, 0.0f, 0.0f, 0.0f);
}

void IntegratorDR::kernel_NextBounce(uint tid, uint bounce, const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_hitPart3, const uint* in_instId,
                                     const float4* in_shadeColor, float4* rayPosAndNear, float4* rayDirAndFar, const float4* wavelengths,
                                     float4* accumColor, float4* accumThoroughput, MisData* misPrev, uint* rayFlags, 
                                     const float* drands, const float* dparams)
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
  
  if(m_materials[matId].mtype == MAT_TYPE_LIGHT_SOURCE)
  {
    const uint   texId     = m_materials[matId].texid[0];
    const float2 texCoordT = mulRows2x4(m_materials[matId].row0[0], m_materials[matId].row1[0], hit.uv);
    const float4 texColor  = Tex2DFetchAD(texId, texCoordT, dparams);
    const uint   lightId   = m_instIdToLightInstId[*in_instId]; 
    
    const float4 emissColor = m_materials[matId].colors[EMISSION_COLOR];
    float4 lightIntensity   = emissColor * texColor;

    if(lightId != 0xFFFFFFFF)
    {
      const float lightCos = dot(to_float3(*rayDirAndFar), to_float3(m_lights[lightId].norm));
      const float lightDirectionAtten = (lightCos < 0.0f || m_lights[lightId].geomType == LIGHT_GEOM_SPHERE) ? 1.0f : 0.0f;
      lightIntensity = LightIntensity(lightId, wavelengths, to_float3(*rayDirAndFar), dparams)*texColor*lightDirectionAtten;
    }

    float misWeight = 1.0f;
    if(m_intergatorType == INTEGRATOR_MIS_PT) 
    {
      if(bounce > 0)
      {
        if(lightId != 0xFFFFFFFF)
        {
          const float lgtPdf  = LightPdfSelectRev(lightId) * LightEvalPDF(lightId, ray_pos, ray_dir, hit.pos, hit.norm, dparams);
          misWeight           = misWeightHeuristic(prevPdfW, lgtPdf);
          if (prevPdfW <= 0.0f) // specular bounce
            misWeight = 1.0f;
        }
      }
    }
    else if(m_intergatorType == INTEGRATOR_SHADOW_PT && hasNonSpecular(currRayFlags))
      misWeight = 0.0f;
    
    if(m_renderLayer >= 1 && bounce < m_renderLayer) // skip some number of bounces if this is set
      misWeight = 0.0f;

    float4 currAccumColor      = *accumColor;
    float4 currAccumThroughput = *accumThoroughput;
    
    currAccumColor += currAccumThroughput * lightIntensity * misWeight;
    
    *accumColor = currAccumColor;
    *rayFlags   = currRayFlags | (RAY_FLAG_IS_DEAD | RAY_FLAG_HIT_LIGHT);
    return;
  }
  
  const BsdfSample matSam = MaterialSampleAndEval(matId, bounce, lambda, (-1.0f)*ray_dir, hit.norm, hit.tang, hit.uv, misPrev, currRayFlags, drands, dparams);
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
    
    //*accumColor += currThoroughput*float4(1,1,1,1);
    *accumColor      += currThoroughput*shadeColor;
    *accumThoroughput = currThoroughput*cosTheta*bxdfVal; 
  }

  *rayPosAndNear = to_float4(OffsRayPos(hit.pos, hit.norm, matSam.dir), 0.0f); // todo: use flatNormal for offset
  *rayDirAndFar  = to_float4(matSam.dir, FLT_MAX);
  *rayFlags      = currRayFlags | matSam.flags;
}

float4 IntegratorDR::EnvironmentColor(float3 a_dir, const float* dparams)
{
  return m_envColor;
}

void IntegratorDR::kernel_HitEnvironment(uint tid, const uint* rayFlags, const float4* rayDirAndFar, const MisData* a_prevMisData, const float4* accumThoroughput,
                                         float4* accumColor, const float* dparams)
{
  if(tid >= m_maxThreadId)
    return;
  const uint currRayFlags = *rayFlags;
  if(!isOutOfScene(currRayFlags))
    return;
  
  // TODO: HDRI maps
  const float4 envData  = EnvironmentColor(to_float3(*rayDirAndFar), dparams);
  // const float3 envColor = to_float3(envData)/envData.w;    // explicitly account for pdf; when MIS will be enabled, need to deal with MIS weight also!

  const float4 envColor = envData;
  if(m_intergatorType == INTEGRATOR_STUPID_PT)     // todo: when explicit sampling will be added, disable contribution here for 'INTEGRATOR_SHADOW_PT'
    *accumColor = (*accumThoroughput) * envColor;
  else
    *accumColor += (*accumThoroughput) * envColor;
}


void IntegratorDR::kernel_ContributeToImage(uint tid, const uint* rayFlags, uint channels, const float4* a_accumColor, const uint* in_pakedXY,
                                            const float4* wavelengths, float* out_color, const float* dparams)
{
  if(tid >= m_maxThreadId)
    return;

  const uint XY = in_pakedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;

  float4 specSamples = *a_accumColor; 

  float4 tmpVal      = specSamples*m_camRespoceRGB;
  float3 rgb         = to_float3(tmpVal);
  if(KSPEC_SPECTRAL_RENDERING!=0 && m_spectral_mode != 0) 
  {
    float4 waves = *wavelengths;
    
    if(m_camResponseSpectrumId[0] < 0)
    {
      const float3 xyz = SpectrumToXYZ(specSamples, waves, LAMBDA_MIN, LAMBDA_MAX, 
                                       m_cie_x.data(), m_cie_y.data(), m_cie_z.data(), terminateWavelngths(*rayFlags));
      rgb = XYZToRGB(xyz);
    }
    else
    {
      float4 responceX, responceY, responceZ;
      {
        int specId = m_camResponseSpectrumId[0];
        if(specId >= 0)
        {
          const uint2 data  = m_spec_offset_sz[specId];
          const uint offset = data.x;
          const uint size   = data.y;
          responceX = SampleUniformSpectrum(m_spec_values.data() + offset, waves, size);
        }
        else
          responceX = float4(1,1,1,1);

        specId = m_camResponseSpectrumId[1];
        if(specId >= 0)
        {
          const uint2 data  = m_spec_offset_sz[specId];
          const uint offset = data.x;
          const uint size   = data.y;
          responceY = SampleUniformSpectrum(m_spec_values.data() + offset, waves, size);
        }
        else
          responceY = responceX;

        specId = m_camResponseSpectrumId[2];
        if(specId >= 0)
        {
          const uint2 data  = m_spec_offset_sz[specId];
          const uint offset = data.x;
          const uint size   = data.y;
          responceZ = SampleUniformSpectrum(m_spec_values.data() + offset, waves, size);
        }
        else
          responceZ = responceY;
      }

      float3 xyz = float3(0,0,0);
      for (uint32_t i = 0; i < SPECTRUM_SAMPLE_SZ; ++i) {
        xyz.x += specSamples[i]*responceX[i];
        xyz.y += specSamples[i]*responceY[i];
        xyz.z += specSamples[i]*responceZ[i]; 
      } 

      if(m_camResponseType == CAM_RESPONCE_XYZ)
        rgb = XYZToRGB(xyz);
      else
        rgb = xyz;
    }
  }

  float4 colorRes = m_exposureMult * to_float4(rgb, 1.0f);
  //if(x == 415 && (y == 256-130-1))
  //{
  //  int a = 2;
  //  //colorRes = float4(1,0,0,0);
  //}
  
  if(channels == 1)
  {
    const float mono = 0.2126f*colorRes.x + 0.7152f*colorRes.y + 0.0722f*colorRes.z;
    out_color[y*m_winWidth+x] += mono;
  }
  else if(channels <= 4)
  {
    out_color[(y*m_winWidth+x)*channels + 0] += colorRes.x;
    out_color[(y*m_winWidth+x)*channels + 1] += colorRes.y;
    out_color[(y*m_winWidth+x)*channels + 2] += colorRes.z;
  }
  else
  {
    auto waves = (*wavelengths);
    auto color = (*a_accumColor)*m_exposureMult;
    for(int i=0;i<4;i++) {
      const float t         = (waves[i] - LAMBDA_MIN)/(LAMBDA_MAX-LAMBDA_MIN);
      const int channelId   = std::min(int(float(channels)*t), int(channels)-1);
      const int offsetPixel = int(y)*m_winWidth + int(x);
      const int offsetLayer = channelId*m_winWidth*m_winHeight;
      out_color[offsetLayer + offsetPixel] += color[i];
    }
  }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

uint IntegratorDR::RemapMaterialId(uint a_mId, int a_instId, const float* dparams)
{
  const int remapListId  = m_remapInst[a_instId];
  if(remapListId == -1)
    return a_mId;

  const int r_offset     = m_allRemapListsOffsets[remapListId];
  const int r_size       = m_allRemapListsOffsets[remapListId+1] - r_offset;
  const int2 offsAndSize = int2(r_offset, r_size);
  
  uint res = a_mId;
  
  // for (int i = 0; i < offsAndSize.y; i++) // linear search version
  // {
  //   int idRemapFrom = m_allRemapLists[offsAndSize.x + i * 2 + 0];
  //   int idRemapTo   = m_allRemapLists[offsAndSize.x + i * 2 + 1];
  //   if (idRemapFrom == a_mId) {
  //     res = idRemapTo;
  //     break;
  //   }
  // }

  int low  = 0;
  int high = offsAndSize.y - 1;              // binary search version
  
  while (low <= high)
  {
    const int mid         = low + ((high - low) / 2);
    const int idRemapFrom = m_allRemapLists[offsAndSize.x + mid * 2 + 0];
    if (uint(idRemapFrom) >= a_mId)
      high = mid - 1;
    else //if(a[mid]<i)
      low = mid + 1;
  }

  if (high+1 < offsAndSize.y)
  {
    const int idRemapFrom = m_allRemapLists[offsAndSize.x + (high + 1) * 2 + 0];
    const int idRemapTo   = m_allRemapLists[offsAndSize.x + (high + 1) * 2 + 1];
    res                   = (uint(idRemapFrom) == a_mId) ? uint(idRemapTo) : a_mId;
  }

  return res;
} 

float IntegratorDR::LightEvalPDF(int a_lightId, float3 illuminationPoint, float3 ray_dir, const float3 lpos, const float3 lnorm, const float* dparams)
{
  const uint gtype    = m_lights[a_lightId].geomType;
  const float hitDist = length(illuminationPoint - lpos);
  
  float cosVal = 1.0f;
  switch(gtype)
  {
    case LIGHT_GEOM_SPHERE:
    {
      // const float  lradius = m_lights[a_lightId].size.x;
      // const float3 lcenter = to_float3(m_lights[a_lightId].pos);
      //if (DistanceSquared(illuminationPoint, lcenter) - lradius*lradius <= 0.0f)
      //  return 1.0f;
      const float3 dirToV  = normalize(lpos - illuminationPoint);
      cosVal = std::abs(dot(dirToV, lnorm));
    }
    break;

    case LIGHT_GEOM_POINT:
    {
      if(m_lights[a_lightId].distType == LIGHT_DIST_OMNI)
        cosVal = 1.0f;
      else
        cosVal = std::max(dot(ray_dir, -1.0f*lnorm), 0.0f);
    };
    break;

    default:
    cosVal  = std::max(dot(ray_dir, -1.0f*lnorm), 0.0f);
    break;
  };
  
  return PdfAtoW(m_lights[a_lightId].pdfA, hitDist, cosVal);
}

float4 IntegratorDR::SampleMatColorParamSpectrum(uint32_t matId, float4 a_wavelengths, uint32_t paramId, uint32_t paramSpecId, const float* dparams)
{  
  float4 res = m_materials[matId].colors[paramId];
  if(a_wavelengths[0] == 0.0f)
    return res;

  const uint specId = m_materials[matId].spdid[paramSpecId];

  if(specId < 0xFFFFFFFF)
  {
    const uint2 data  = m_spec_offset_sz[specId];
    const uint offset = data.x;
    const uint size   = data.y;
    res = SampleUniformSpectrum(m_spec_values.data() + offset, a_wavelengths, size);
  }

  return res;
}

float4 IntegratorDR::SampleMatParamSpectrum(uint32_t matId, float4 a_wavelengths, uint32_t paramId, uint32_t paramSpecId, const float* dparams)
{  
  float4 res = float4(m_materials[matId].data[paramId]);
  if(a_wavelengths[0] == 0.0f)
    return res;

  const uint specId = m_materials[matId].spdid[paramSpecId];

  if(specId < 0xFFFFFFFF)
  {
    const uint2 data  = m_spec_offset_sz[specId];
    const uint offset = data.x;
    const uint size   = data.y;
    res = SampleUniformSpectrum(m_spec_values.data() + offset, a_wavelengths, size);
  }

  return res;
}

LightSample IntegratorDR::LightSampleRev(int a_lightId, float2 rands, float3 illiminationPoint, const float* dparams)
{
  const uint gtype = m_lights[a_lightId].geomType;
  switch(gtype)
  {
    case LIGHT_GEOM_DIRECT: return directLightSampleRev(m_lights.data() + a_lightId, rands, illiminationPoint);
    case LIGHT_GEOM_SPHERE: return sphereLightSampleRev(m_lights.data() + a_lightId, rands);
    case LIGHT_GEOM_POINT:  return pointLightSampleRev (m_lights.data() + a_lightId);
    default:                return areaLightSampleRev  (m_lights.data() + a_lightId, rands);
  };
}

float4 IntegratorDR::LightIntensity(uint a_lightId, const float4* a_wavelengths, float3 a_rayDir, const float* dparams)
{
  float4 lightColor = m_lights[a_lightId].intensity;  
  if(KSPEC_SPECTRAL_RENDERING !=0 && m_spectral_mode != 0)
  {
    const uint specId = m_lights[a_lightId].specId;
  
    if(specId < 0xFFFFFFFF)
    {
      // lightColor = SampleSpectrum(m_spectra.data() + specId, *a_wavelengths);
      const uint2 data  = m_spec_offset_sz[specId];
      const uint offset = data.x;
      const uint size   = data.y;
      lightColor = SampleUniformSpectrum(m_spec_values.data() + offset, *a_wavelengths, size);
    }
  }
  lightColor *= m_lights[a_lightId].mult;

  uint iesId = m_lights[a_lightId].iesId;
  if(iesId != uint(-1))
  {
    float sintheta        = 0.0f;
    const float2 texCoord = sphereMapTo2DTexCoord((-1.0f)*a_rayDir, &sintheta);
    const float4 texColor = Tex2DFetchAD(iesId, texCoord, dparams);
    lightColor *= texColor;
  }

  return lightColor;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

float4 IntegratorDR::Tex2DFetchAD(uint texId, float2 a_uv, const float* tex_data)
{
  const auto info = m_texAddressTable[texId];

  if(info.offset != size_t(-1) && tex_data != nullptr && m_gradMode != 0) 
  {
    const float m_fw     = info.fwidth;
    const float m_fh     = info.fheight;
    const int tex_width  = info.width;
    const int tex_height = info.height;
    
    float ffx = a_uv.x * m_fw - 0.5f; // a_texCoord should not be very large, so that the float does not overflow later. 
    float ffy = a_uv.y * m_fh - 0.5f; // This is left to the responsibility of the top level.
    
    auto sampler = m_textures[texId]->sampler();

    if ((sampler.addressU == Sampler::AddressMode::CLAMP) != 0 && ffx < 0) ffx = 0.0f;
    if ((sampler.addressV == Sampler::AddressMode::CLAMP) != 0 && ffy < 0) ffy = 0.0f;
    
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

    // Calculate the weighted sum of pixels (for each color channel)
    //
    if(info.channels == 4)
    {
      const float4 f1    = float4(tex_data[info.offset+offsets.x*4+0], tex_data[info.offset+offsets.x*4+1], tex_data[info.offset+offsets.x*4+2], tex_data[info.offset+offsets.x*4+3]);
      const float4 f2    = float4(tex_data[info.offset+offsets.y*4+0], tex_data[info.offset+offsets.y*4+1], tex_data[info.offset+offsets.y*4+2], tex_data[info.offset+offsets.y*4+3]);
      const float4 f3    = float4(tex_data[info.offset+offsets.z*4+0], tex_data[info.offset+offsets.z*4+1], tex_data[info.offset+offsets.z*4+2], tex_data[info.offset+offsets.z*4+3]);
      const float4 f4    = float4(tex_data[info.offset+offsets.w*4+0], tex_data[info.offset+offsets.w*4+1], tex_data[info.offset+offsets.w*4+2], tex_data[info.offset+offsets.w*4+3]);
  
      const float outr = f1.x * w1 + f2.x * w2 + f3.x * w3 + f4.x * w4;
      const float outg = f1.y * w1 + f2.y * w2 + f3.y * w3 + f4.y * w4;
      const float outb = f1.z * w1 + f2.z * w2 + f3.z * w3 + f4.z * w4;
      const float outa = f1.w * w1 + f2.w * w2 + f3.w * w3 + f4.w * w4;
      
      return float4(outr, outg, outb, outa);
    }
    else
    {
      const float f1 = tex_data[info.offset+offsets.x];
      const float f2 = tex_data[info.offset+offsets.y];
      const float f3 = tex_data[info.offset+offsets.z];
      const float f4 = tex_data[info.offset+offsets.w];
  
      const float outVal = f1 * w1 + f2 * w2 + f3 * w3 + f4 * w4;
      
      return float4(outVal, outVal, outVal, outVal);
    }
  }
  else
    return m_textures[texId]->sample(a_uv);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void IntegratorDR::kernel_InitEyeRay(uint tid, const uint* packedXY, float4* rayPosAndNear, float4* rayDirAndFar, const float* a_data) // (tid,tidX,tidY,tidZ) are SPECIAL PREDEFINED NAMES!!!
{
  if(tid >= m_maxThreadId)
    return;
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

bool IntegratorDR::kernel_RayTrace(uint tid, const float4* rayPosAndNear, float4* rayDirAndFar,
                                   Lite_Hit* out_hit, float2* out_bars, const float* a_data)
{
  if(tid >= m_maxThreadId)
    return false;
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

  const uint   texId     = m_materials[matId].texid[0];
  const float2 texCoordT = mulRows2x4(m_materials[matId].row0[0], m_materials[matId].row1[0], hitTexCoord);
  const float4 texColor  = Tex2DFetchAD(texId, texCoordT, a_data); 
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
  kernel_InitEyeRay(tid, m_packedXY.data(), &rayPosAndNear, &rayDirAndFar, a_data);

  Lite_Hit hit; 
  float2   baricentrics; 
  if(!kernel_RayTrace(tid, &rayPosAndNear, &rayDirAndFar, &hit, &baricentrics, a_data))
    return float4(0,0,0,0);
  
  float4 finalColor;
  kernel_CalcRayColor(tid, &hit, &baricentrics, &finalColor, m_packedXY.data(), out_color, a_data);
  return finalColor;
}

extern float __enzyme_autodiff(void*, ...);
int enzyme_const, enzyme_dup, enzyme_out;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double RegLossImage1D(size_t a_size, const float* data)
{
  double summ = 0.0f;

  for(size_t i=1;i<a_size-1;i++) {
    float diffLeft  = data[i] - data[i-1];
    float diffRight = data[i] - data[i+1];
    summ += double(diffLeft*diffLeft + diffRight*diffRight);
  }

  return summ/double(a_size);
}

double RegLossImage2D(int w, int h, const float* data)
{
  double summ = 0.0f;
  std::vector<double> lines(h);

  for(int y=1;y<h-1;y++) {
    lines[y] = 0.0;
    for(int x=1;x<w-1;x++) {
      float diffTop    = data[y*w+x] - data[(y+1)*w+x];
      float diffBottom = data[y*w+x] - data[(y-1)*w+x];
      float diffLeft   = data[y*w+x] - data[y*w+x-1];
      float diffRight  = data[y*w+x] - data[y*w+x+1];
      lines[y] += std::sqrt(double(diffLeft*diffLeft + diffRight*diffRight + diffTop*diffTop + diffBottom*diffBottom));
    }
    summ += lines[y];
  }
 
  return summ;
}

using LiteMath::dot3;

double RegLossImage2D4f(int w, int h, const float* data)
{
  double summ = 0.0f;
  std::vector<double> lines(h);

  for(int y=1;y<h-1;y++) {
    lines[y] = 0.0;
    for(int x=1;x<w-1;x++) {
      float4 p0 = float4(data[(y*w+x)*4+0],     data[(y*w+x)*4+1],     data[(y*w+x)*4+2],     data[(y*w+x)*4+3]);
      float4 p1 = float4(data[((y+1)*w+x)*4+0], data[((y+1)*w+x)*4+1], data[((y+1)*w+x)*4+2], data[((y+1)*w+x)*4+3]);
      float4 p2 = float4(data[((y-1)*w+x)*4+0], data[((y-1)*w+x)*4+1], data[((y-1)*w+x)*4+2], data[((y-1)*w+x)*4+3]); 
      float4 p3 = float4(data[(y*w+x-1)*4+0], data[(y*w+x-1)*4+1], data[(y*w+x-1)*4+2], data[(y*w+x-1)*4+3]);
      float4 p4 = float4(data[(y*w+x+1)*4+0], data[(y*w+x+1)*4+1], data[(y*w+x+1)*4+2], data[(y*w+x+1)*4+3]);

      float4 diffTop    = p0 - p1;
      float4 diffBottom = p0 - p2;
      float4 diffLeft   = p0 - p3;
      float4 diffRight  = p0 - p4;

      lines[y] += std::sqrt(double(dot3(diffLeft,diffLeft) + dot3(diffRight,diffRight) + dot3(diffTop,diffTop) + dot3(diffBottom,diffBottom)));
    }
  }
  
  summ = 0.0;
  for(int sy=0; sy < h/2; sy++) {
    int index1 = h/2 + sy;
    int index2 = h/2 - sy;
    if(index1 < h-1)
      summ += lines[index1];
    if(sy!=0 && index2 >= 1)
      summ += lines[index2];
  }
 
  return summ; //double(w*h);
}

void Image1DRegularizer(size_t a_size, const float* data, float* grad)
{
  __enzyme_autodiff((void*)RegLossImage1D, 
                           enzyme_const, a_size,
                           enzyme_dup,   data, grad);

}

void Image2D4fRegularizer(int w, int h, const float* data, float* grad)
{
  __enzyme_autodiff((void*)RegLossImage2D4f, 
                           enzyme_const, w,
                           enzyme_const, h,
                           enzyme_dup,   data, grad);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float PixelLossRT(IntegratorDR* __restrict__ pIntegrator,
                  const float*  __restrict__ a_refImg,
                        float*  __restrict__ out_color,
                  const float*  __restrict__ a_data, 
                  const uint*   __restrict__ in_pakedXY, 
                  uint tid, uint channels, uint pitch,
                  float*  __restrict__       outLoss)
{
  const uint XY = in_pakedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;

  const uint yRef  = pIntegrator->m_winHeight - y - 1; // in input images and when load data from HDD y has different direction
  float4 colorRend = pIntegrator->CastRayDR(tid, channels, out_color, a_data);
  float4 colorRef  = float4(a_refImg[(yRef*pitch+x)*channels + 0], 
                            a_refImg[(yRef*pitch+x)*channels + 1], 
                            a_refImg[(yRef*pitch+x)*channels + 2], 0.0f);

  float4 diff = colorRend - colorRef;
  float loss = LiteMath::dot3(diff, diff);
  (*outLoss) = loss;
  return loss;
}                      

float IntegratorDR::RayTraceDR(uint tid, uint channels, float* out_color, uint a_passNum,
                                const float* a_refImg, const float* a_data, float* a_dataGrad, size_t a_gradSize)
{
  memset(a_dataGrad, 0, sizeof(float)*a_gradSize);

  // init separate gradient for each thread
  //
  //std::vector<float> grads[MAXTHREADS_CPU];
  //for(int i=0;i<MAXTHREADS_CPU;i++)
  //   std::fill(grads[i].begin(), grads[i].end(), 0.0f);

  //double avgLoss = 0.0;
  auto start = std::chrono::high_resolution_clock::now();
  //#ifndef _DEBUG
  //#pragma omp parallel for default(shared) // num_threads(MAXTHREADS_CPU)
  //#endif

  float avgLoss = 0.0f;
  
  if(m_gradMode != 0)
  {
    for (int i = 0; i < int(tid); ++i) {
      float lossVal = 0.0f;
      __enzyme_autodiff((void*)PixelLossRT, 
                         enzyme_const, this,
                         enzyme_const, a_refImg,
                         enzyme_const, out_color,
                         enzyme_dup,   a_data, a_dataGrad,
                         enzyme_const, m_packedXY.data(),
                         enzyme_const, uint(i),
                         enzyme_const, channels,
                         enzyme_const, m_winWidth,
                         enzyme_const, &lossVal);
      avgLoss += float(lossVal)/float(a_passNum);
    }
  }
  else
  {
    for (int i = 0; i < int(tid); ++i) {
      float lossVal = PixelLossRT(this, a_refImg, out_color, a_data, m_packedXY.data(),
                                uint(i), channels, m_winWidth, &lossVal);
      avgLoss += float(lossVal)/float(a_passNum);
    }
  }

  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;

  // accumulate gradient from different threads (parallel reduction/hist)
  //

  //for(int i=0;i<MAXTHREADS_CPU;i++) 
  //  for(size_t j=0;j<a_gradSize; j++)
  //    a_dataGrad[j] += grads[i][j];

  //avgLoss /= float(m_winWidth*m_winHeight);
  //std::cout << "avgLoss = " << avgLoss << std::endl;
  
  //std::ofstream fout("z_grad.txt");
  //for(size_t i=0; i<a_gradSize; i++)
  //  fout << a_dataGrad[i]/float(a_passNum) << std::endl;
  //fout.close();

  return avgLoss;
}

float4 IntegratorDR::PathTraceReplay(uint tid, uint channels, uint cpuThreadId, float* out_color, 
                                     const float* drands, const float* dparams)
{
  float4 accumColor, accumThroughput;
  float4 rayPosAndNear, rayDirAndFar;
  float4 wavelengths; 
  MisData   mis;
  uint      rayFlags;
  kernel_InitEyeRay2(tid, &rayPosAndNear, &rayDirAndFar, &wavelengths, &accumColor, &accumThroughput, &rayFlags, &mis, drands, dparams);

  for(uint depth = 0; depth < m_traceDepth; depth++) 
  {
    float4 hitPart1, hitPart2, hitPart3, shadeColor;
    uint   instId;
    kernel_RayTrace2(tid, depth, cpuThreadId, &rayPosAndNear, &rayDirAndFar, &hitPart1, &hitPart2, &hitPart3, &instId, &rayFlags, dparams);
    if(isDeadRay(rayFlags))
      break;
    
    kernel_SampleLightSource(tid, cpuThreadId, &rayPosAndNear, &rayDirAndFar, &wavelengths, &hitPart1, &hitPart2, &hitPart3, &rayFlags, depth,
                             &shadeColor, drands, dparams);
    
    kernel_NextBounce(tid, depth, &hitPart1, &hitPart2, &hitPart3, &instId, &shadeColor,
                      &rayPosAndNear, &rayDirAndFar, &wavelengths, &accumColor, &accumThroughput, &mis, &rayFlags, 
                      drands, dparams);
    

    if(isDeadRay(rayFlags))
      break;
  }

  kernel_HitEnvironment(tid, &rayFlags, &rayDirAndFar, &mis, &accumThroughput,
                        &accumColor, dparams);

  kernel_ContributeToImage(tid, &rayFlags, channels, &accumColor, m_packedXY.data(), &wavelengths, out_color, dparams);
  return accumColor;
}

float PixelLossPT(IntegratorDR* __restrict__ pIntegrator,
                  uint tid, uint channels, uint pitch, uint cpuThreadId,
                  const float*  __restrict__ a_refImg,
                        float*  __restrict__ out_color,
                  const uint*   __restrict__ in_pakedXY, 
                  float*        __restrict__ outLoss,
                  const float*  __restrict__ a_drands,
                  const float*  __restrict__ a_dparams)
{
  const uint XY = in_pakedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;

  float4 colorRend = pIntegrator->PathTraceReplay(tid, channels, cpuThreadId, out_color, 
                                                  a_drands, a_dparams);

  const uint yRef  = pIntegrator->m_winHeight - y - 1; // in input images and when load data from HDD y has different direction
  float4 colorRef  = float4(a_refImg[(yRef*pitch+x)*channels + 0], 
                            a_refImg[(yRef*pitch+x)*channels + 1], 
                            a_refImg[(yRef*pitch+x)*channels + 2], 0.0f);

  float4 diff = colorRend - colorRef;
  float loss = LiteMath::dot3(diff, diff);
  (*outLoss) = loss;
  return loss;
}   


float IntegratorDR::PathTraceDR(uint tid, uint channels, float* out_color, uint a_passNum,
                                const float* a_refImg, const float* a_data, float* a_dataGrad, size_t a_gradSize)
{
  m_disableImageContrib = 1;
  memset(a_dataGrad, 0, sizeof(float)*a_gradSize);

  // init separate gradient for each thread
  //
  std::vector<float> grads[MAXTHREADS_CPU];
  for(int i=0;i<MAXTHREADS_CPU;i++) {
    grads[i].resize(a_gradSize);
    std::fill(grads[i].begin(), grads[i].end(), 0.0f);
  }

  //double avgLoss = 0.0;
  auto start = std::chrono::high_resolution_clock::now();

  float avgLoss = 0.0f;
  
  if(m_gradMode != 0)
  {
    #ifndef _DEBUG
    #pragma omp parallel for default(shared) num_threads(MAXTHREADS_CPU)
    #endif
    for (int i = 0; i < int(tid); ++i) {
      float lossVal = 0.0f;
      for(int passId = 0; passId < int(a_passNum); passId++) {

        this->PathTrace(i, channels, out_color); // record non differentiable data during common PT 
        
        auto cpuThreadId = omp_get_thread_num();

        __enzyme_autodiff((void*)PixelLossPT, 
                           enzyme_const, this,
                           enzyme_const, uint(i),
                           enzyme_const, channels,
                           enzyme_const, m_winWidth,
                           enzyme_const, cpuThreadId,
                           enzyme_const, a_refImg,
                           enzyme_const, out_color,
                           enzyme_const, m_packedXY.data(),
                           enzyme_const, &lossVal,
                           enzyme_const, this->m_recorded[cpuThreadId].perBounceRands.data(),
                           enzyme_dup,   a_data, grads[cpuThreadId].data());

        //lossVal = PixelLossPT(this, uint(i), channels, m_winWidth, cpuThreadId, 
        //                      a_refImg, out_color, m_packedXY.data(), &lossVal,
        //                      this->m_recorded[cpuThreadId].perBounceRands.data(), a_data);

        avgLoss += float(lossVal)/float(a_passNum);
      }
    }
  }
  else
  {
    //for (int i = 0; i < int(tid); ++i) {
    //  float lossVal = PixelLossPT(this, a_refImg, out_color, a_data, m_packedXY.data(),
    //                              uint(i), channels, m_winWidth, &lossVal);
    //  avgLoss += float(lossVal)/float(a_passNum);
    //}
  }

  diffPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;

  // accumulate gradient from different threads (parallel reduction/hist)
  //
  for(int i=0;i<MAXTHREADS_CPU;i++) 
    for(size_t j=0;j<a_gradSize; j++)
      a_dataGrad[j] += grads[i][j];

  avgLoss /= float(m_winWidth*m_winHeight);
  //std::cout << "avgLoss = " << avgLoss << std::endl;
  //std::cout.flush();
  
  //std::ofstream fout("z_grad.txt");
  //for(size_t i=0; i<a_gradSize; i++)
  //  fout << a_dataGrad[i]/float(a_passNum) << std::endl;
  //fout.close();

  m_disableImageContrib = 0;
  return avgLoss;
}

