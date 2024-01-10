#include "integrator_pt.h"
#include "include/crandom.h"

#include "include/cmaterial.h"
#include "include/cmat_gltf.h"
#include "include/cmat_conductor.h"
#include "include/cmat_glass.h"
#include "include/cmat_film.h"
#include "include/cmat_diffuse.h"

#include <chrono>
#include <string>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

LightSample Integrator::LightSampleRev(int a_lightId, float2 rands, float3 illiminationPoint)
{
  const uint gtype = m_lights[a_lightId].geomType;
  switch(gtype)
  {
    case LIGHT_GEOM_DIRECT: return directLightSampleRev(m_lights.data() + a_lightId, rands, illiminationPoint);
    case LIGHT_GEOM_SPHERE: return sphereLightSampleRev(m_lights.data() + a_lightId, rands);
    default:                return areaLightSampleRev  (m_lights.data() + a_lightId, rands);
  };
}

float Integrator::LightPdfSelectRev(int a_lightId) 
{ 
  return 1.0f/float(m_lights.size()); // uniform select
}

//static inline float DistanceSquared(float3 a, float3 b)
//{
//  const float3 diff = b - a;
//  return dot(diff, diff);
//}

float Integrator::LightEvalPDF(int a_lightId, float3 illuminationPoint, float3 ray_dir, const float3 lpos, const float3 lnorm)
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

    default:
    cosVal  = std::max(dot(ray_dir, -1.0f*lnorm), 0.0f);
    break;
  };
  
  return PdfAtoW(m_lights[a_lightId].pdfA, hitDist, cosVal);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


uint32_t Integrator::MaterialBlendSampleAndEval(uint a_materialId, float4 wavelengths, RandomGen* a_gen, float3 v, float3 n, float2 tc, 
                                                  MisData* a_misPrev, BsdfSample* a_pRes)
{
  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const uint   texId     = as_uint(m_materials[a_materialId].data[BLEND_TEXID0]);
  const float4 weightDat = m_textures[texId]->sample(texCoordT);
  const float  weightTex = weightDat.x;
  const float  weight    = m_materials[a_materialId].data[BLEND_WEIGHT] * weightTex;

  const uint matId1 = as_uint(m_materials[a_materialId].data[BLEND_MAT_ID_1]);
  const uint matId2 = as_uint(m_materials[a_materialId].data[BLEND_MAT_ID_2]);

  uint32_t selectedMatId = matId1;
  const float select = rndFloat1_Pseudo(a_gen);
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

MatIdWeightPair Integrator::MaterialBlendEval(MatIdWeight a_mat, float4 wavelengths, float3 l, float3 v, float3 n, float2 tc)
{
  const float2 texCoordT = mulRows2x4(m_materials[a_mat.id].row0[0], m_materials[a_mat.id].row1[0], tc);
  const uint   texId     = as_uint(m_materials[a_mat.id].data[BLEND_TEXID0]);
  const float4 weightDat = m_textures[texId]->sample(texCoordT);
  const float  weightTex = weightDat.x;
  const float  weight    = m_materials[a_mat.id].data[BLEND_WEIGHT] * weightTex;

  const uint matId1      = as_uint(m_materials[a_mat.id].data[BLEND_MAT_ID_1]);
  const uint matId2      = as_uint(m_materials[a_mat.id].data[BLEND_MAT_ID_2]);

  MatIdWeight p1, p2;
  p1.id     = matId1;
  p1.weight = a_mat.weight * (1.0f - weight);
  p2.id     = matId2;
  p2.weight = a_mat.weight * weight;

  return make_weight_pair(p1, p2);
}


BsdfSample Integrator::MaterialSampleAndEval(uint a_materialId, float4 wavelengths, RandomGen* a_gen, float3 v, float3 n, float2 tc, 
                                             MisData* a_misPrev, const uint a_currRayFlags)
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
  while(mtype == MAT_TYPE_BLEND)
  {
    currMatId = MaterialBlendSampleAndEval(currMatId, wavelengths, a_gen, v, n, tc, a_misPrev, &res);
    mtype     = as_uint(m_materials[currMatId].data[UINT_MTYPE]);
  }

  //assert(currMatId < m_materials.size());

  const float2 texCoordT = mulRows2x4(m_materials[currMatId].row0[0], m_materials[currMatId].row1[0], tc);
  switch(mtype)
  {
    case MAT_TYPE_GLTF:
    if(KSPEC_MAT_TYPE_GLTF != 0)
    {
      const float4 rands = rndFloat4_Pseudo(a_gen);
      
      const uint   texId     = as_uint(m_materials[currMatId].data[GLTF_UINT_TEXID0]);
      const float4 texColor  = (m_textures[texId]->sample(texCoordT));
      const float4 color     = (m_materials[currMatId].colors[GLTF_COLOR_BASE])*texColor;
      gltfSampleAndEval(m_materials.data() + currMatId, rands, v, n, tc, color, &res);
    }
    break;
    case MAT_TYPE_GLASS: 
    if(KSPEC_MAT_TYPE_GLASS != 0)
    {
      const float4 rands = rndFloat4_Pseudo(a_gen);

      glassSampleAndEval(m_materials.data() + currMatId, rands, v, n, tc, &res, a_misPrev);
    }
    break;
    case MAT_TYPE_CONDUCTOR:
    if(KSPEC_MAT_TYPE_CONDUCTOR != 0)
    {
      const float4 rands = rndFloat4_Pseudo(a_gen);

      const uint   texId     = as_uint(m_materials[currMatId].data[CONDUCTOR_TEXID0]);
      const float3 alphaTex  = to_float3(m_textures[texId]->sample(texCoordT));
      
      const float2 alpha   = float2(m_materials[currMatId].data[CONDUCTOR_ROUGH_V], m_materials[currMatId].data[CONDUCTOR_ROUGH_U]);
      const float4 etaSpec = SampleMatParamSpectrum(currMatId, wavelengths, CONDUCTOR_ETA, CONDUCTOR_ETA_SPECID);
      const float4 kSpec   = SampleMatParamSpectrum(currMatId, wavelengths, CONDUCTOR_K, CONDUCTOR_K_SPECID);
      if(trEffectivelySmooth(alpha))
        conductorSmoothSampleAndEval(m_materials.data() + currMatId, etaSpec, kSpec, rands, v, n, tc, &res);
      else
        conductorRoughSampleAndEval(m_materials.data() + currMatId, etaSpec, kSpec, rands, v, n, tc, alphaTex, &res);
    }
    break;
    case MAT_TYPE_THIN_FILM:
    if(KSPEC_MAT_TYPE_THIN_FILM != 0)
    {
      const float4 rands = rndFloat4_Pseudo(a_gen);

      const uint   texId     = as_uint(m_materials[a_materialId].data[FILM_TEXID0]);
      const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
      const float3 alphaTex  = to_float3(m_textures[texId]->sample(texCoordT));
      
      const float2 alpha = float2(m_materials[a_materialId].data[FILM_ROUGH_V], m_materials[a_materialId].data[FILM_ROUGH_U]);

      std::vector<float4> etaSpec;
      std::vector<float4> kSpec;
      std::vector<float> thickness;
      uint layerId = a_materialId;
      do
      {
        etaSpec.push_back(SampleMatParamSpectrum(layerId, wavelengths, FILM_ETA, FILM_ETA_SPECID));
        kSpec.push_back(SampleMatParamSpectrum(layerId, wavelengths, FILM_K, FILM_K_SPECID));
        thickness.push_back(m_materials[layerId].data[FILM_THICKNESS]);
        layerId = as_uint(m_materials[layerId].data[FILM_NEXT_ID]);
      } while (layerId != 0xFFFFFFFF);

      if(trEffectivelySmooth(alpha))
        filmSmoothSampleAndEval(m_materials.data() + a_materialId, etaSpec.data(), kSpec.data(), thickness.data(), etaSpec.size() - 1, wavelengths, rands, v, n, tc, &res);
      else
        filmRoughSampleAndEval(m_materials.data() + a_materialId, etaSpec.data(), kSpec.data(), thickness.data(), etaSpec.size() - 1, wavelengths, rands, v, n, tc, alphaTex, &res);
    }
    break;
    case MAT_TYPE_DIFFUSE:
    if(KSPEC_MAT_TYPE_DIFFUSE != 0)
    {
      const float4 rands = rndFloat4_Pseudo(a_gen);

      const uint   texId       = as_uint(m_materials[currMatId].data[DIFFUSE_TEXID0]);
      const float4 texColor    = (m_textures[texId]->sample(texCoordT));
      const float4 color       = texColor;

      const float4 reflSpec    = SampleMatColorParamSpectrum(currMatId, wavelengths, DIFFUSE_COLOR, DIFFUSE_SPECID);

      diffuseSampleAndEval(m_materials.data() + currMatId, reflSpec, rands, v, n, tc, color, &res);
    }
    break;
    default:
    break;
  }
  
  return res;
}

BsdfEval Integrator::MaterialEval(uint a_materialId, float4 wavelengths, float3 l, float3 v, float3 n, float2 tc)
{
  BsdfEval res;
  {
    res.val = float4(0,0,0,0);
    res.pdf   = 0.0f;
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
        const float4 texColor  = (m_textures[texId]->sample(texCoordT));
        const float4 color     = (m_materials[currMat.id].colors[GLTF_COLOR_BASE]) * texColor;
        gltfEval(m_materials.data() + currMat.id, l, v, n, tc, color, &currVal);

        res.val += currVal.val * currMat.weight;
        res.pdf += currVal.pdf * currMat.weight;

        break;
      }
      case MAT_TYPE_GLASS:
      if(KSPEC_MAT_TYPE_GLASS != 0)
      {
        glassEval(m_materials.data() + currMat.id, l, v, n, tc, float3(0,0,0), &currVal);

        res.val += currVal.val * currMat.weight;
        res.pdf += currVal.pdf * currMat.weight;
        break;
      }
      case MAT_TYPE_CONDUCTOR: 
      if(KSPEC_MAT_TYPE_CONDUCTOR != 0)
      {
        const uint   texId     = as_uint(m_materials[currMat.id].data[CONDUCTOR_TEXID0]);
        const float3 alphaTex  = to_float3(m_textures[texId]->sample(texCoordT));

        const float2 alpha = float2(m_materials[currMat.id].data[CONDUCTOR_ROUGH_V], m_materials[currMat.id].data[CONDUCTOR_ROUGH_U]);

        if(!trEffectivelySmooth(alpha))
        {
          const float4 etaSpec = SampleMatParamSpectrum(currMat.id, wavelengths, CONDUCTOR_ETA, CONDUCTOR_ETA_SPECID);
          const float4 kSpec   = SampleMatParamSpectrum(currMat.id, wavelengths, CONDUCTOR_K, CONDUCTOR_K_SPECID);
          conductorRoughEval(m_materials.data() + currMat.id, etaSpec, kSpec, l, v, n, tc, alphaTex, &currVal);
        }

        res.val += currVal.val * currMat.weight;
        res.pdf += currVal.pdf * currMat.weight;
        break;
      }
      case MAT_TYPE_THIN_FILM: 
      if(KSPEC_MAT_TYPE_THIN_FILM != 0)
      {
        const uint   texId     = as_uint(m_materials[currMat.id].data[FILM_TEXID0]);
        const float3 alphaTex  = to_float3(m_textures[texId]->sample(texCoordT));

        const float2 alpha = float2(m_materials[currMat.id].data[FILM_ROUGH_V], m_materials[currMat.id].data[FILM_ROUGH_U]);

        if(!trEffectivelySmooth(alpha))
        {
          std::vector<float4> etaSpec;
          std::vector<float4> kSpec;
          std::vector<float> thickness;
          uint layerId = a_materialId;
          do
          {
            etaSpec.push_back(SampleMatParamSpectrum(layerId, wavelengths, FILM_ETA, FILM_ETA_SPECID));
            kSpec.push_back(SampleMatParamSpectrum(layerId, wavelengths, FILM_K, FILM_K_SPECID));
            thickness.push_back(m_materials[layerId].data[FILM_THICKNESS]);
            layerId = as_uint(m_materials[layerId].data[FILM_NEXT_ID]);
          } while (layerId != 0xFFFFFFFF);
          filmRoughEval(m_materials.data() + a_materialId, etaSpec.data(), kSpec.data(), thickness.data(), etaSpec.size() - 1, wavelengths, l, v, n, tc, alphaTex, &currVal);
        }

        res.val += currVal.val * currMat.weight;
        res.pdf += currVal.pdf * currMat.weight;
        break;
      }
      case MAT_TYPE_DIFFUSE:
      if(KSPEC_MAT_TYPE_DIFFUSE != 0)
      {
        const uint   texId       = as_uint(m_materials[currMat.id].data[DIFFUSE_TEXID0]);
        const float4 texColor    = (m_textures[texId]->sample(texCoordT));
        const float4 color       = texColor;

        const float4 reflSpec    = SampleMatColorParamSpectrum(currMat.id, wavelengths, DIFFUSE_COLOR, DIFFUSE_SPECID);

        diffuseEval(m_materials.data() + currMat.id, reflSpec, l, v, n, tc, color, &currVal);

        res.val += currVal.val * currMat.weight;
        res.pdf += currVal.pdf * currMat.weight;
        break;
      }
      case MAT_TYPE_BLEND:
      if(KSPEC_MAT_TYPE_BLEND != 0)
      {
        auto childMats = MaterialBlendEval(currMat, wavelengths, l, v, n, tc);
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

float4 Integrator::GetEnvironmentColorAndPdf(float3 a_dir)
{
  return m_envColor;
}

uint Integrator::RemapMaterialId(uint a_mId, int a_instId)
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

void Integrator::GetExecutionTime(const char* a_funcName, float a_out[4])
{
  if(std::string(a_funcName) == "NaivePathTrace" || std::string(a_funcName) == "NaivePathTraceBlock")
    a_out[0] = naivePtTime;
  else if(std::string(a_funcName) == "PathTrace" || std::string(a_funcName) == "PathTraceBlock")
    a_out[0] = shadowPtTime;
  else if(std::string(a_funcName) == "RayTrace" || std::string(a_funcName) == "RayTraceBlock")
    a_out[0] = raytraceTime;
  else if(std::string(a_funcName) == "PathTraceFromInputRays" || std::string(a_funcName) == "PathTraceFromInputRaysBlock")
    a_out[0] = fromRaysPtTime;
}