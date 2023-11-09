#include "integrator_pt.h"
#include "include/crandom.h"

#include "include/cmaterial.h"
#include "include/cmat_gltf.h"
#include "include/cmat_conductor.h"
#include "include/cmat_glass.h"
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
  const float  weightTex = to_float3(m_textures[texId]->sample(texCoordT)).x;
  const float  weight    = m_materials[a_materialId].data[BLEND_WEIGHT] * weightTex;

  const uint matId1 = as_uint(m_materials[a_materialId].data[BLEND_MAT_ID_1]);
  const uint matId2 = as_uint(m_materials[a_materialId].data[BLEND_MAT_ID_2]);

  // BsdfSample res;
  // {
  //   res.val   = float4(0, 0, 0, 0);
  //   res.pdf   = 1.0f;
  //   res.dir   = float3(0,1,0);
  //   res.flags = a_currRayFlags;
  // }

  uint32_t selectedMatId = matId1;
  const float select = rndFloat1_Pseudo(a_gen);
  if(select < weight)
  {
    // res = MaterialSampleAndEval(matId2, wavelengths, a_gen, v, n, tc, a_misPrev, a_currRayFlags);
    a_pRes->pdf *= weight;
    a_pRes->val *= weight;
    selectedMatId = matId2;
  }
  else
  {
    // res = MaterialSampleAndEval(matId1, wavelengths, a_gen, v, n, tc, a_misPrev, a_currRayFlags);
    a_pRes->pdf *= 1.0f - weight;
    a_pRes->val *= 1.0f - weight;
    selectedMatId = matId1;
  }

  return selectedMatId;
}

std::pair<MatIdWeight, MatIdWeight> Integrator::MaterialBlendEval(MatIdWeight a_mat, float4 wavelengths, float3 l, float3 v, float3 n, float2 tc)
{
  const float2 texCoordT = mulRows2x4(m_materials[a_mat.id].row0[0], m_materials[a_mat.id].row1[0], tc);
  const uint   texId     = as_uint(m_materials[a_mat.id].data[BLEND_TEXID0]);
  const float  weightTex = to_float3(m_textures[texId]->sample(texCoordT)).x;
  const float  weight    = m_materials[a_mat.id].data[BLEND_WEIGHT] * weightTex;

  const uint matId1 = as_uint(m_materials[a_mat.id].data[BLEND_MAT_ID_1]);
  const uint matId2 = as_uint(m_materials[a_mat.id].data[BLEND_MAT_ID_2]);

  // BsdfEval res;
  // {
  //   res.val = float4(0, 0, 0, 0);
  //   res.pdf = 0.0f;
  // }

  // BsdfEval res1 = MaterialEval(matId1, wavelengths, l, v, n, tc);
  // BsdfEval res2 = MaterialEval(matId2, wavelengths, l, v, n, tc);
  // res.pdf = weight * res2.pdf + (1.0f - weight) * res1.pdf;
  // res.val = weight * res2.val + (1.0f - weight) * res1.val;

  return std::make_pair(MatIdWeight{matId1, a_mat.weight * (1.0f - weight)}, MatIdWeight{matId2, a_mat.weight * weight});
}

uint32_t stack_pop(uint32_t* stack, uint32_t stack_sz)
{
  uint32_t val = stack[0];
  int i = 0;
  while(stack[i] != 0xFFFFFFFF || i != stack_sz - 2)
  {
    stack[i] = stack[i + 1];
    ++i;
  }
  return val;
}

void stack_push(uint32_t val, uint32_t* stack, uint32_t stack_sz)
{
  int i = stack_sz - 1;
  assert(stack[i] == 0xFFFFFFFF);

  while(i > 0)
  {
    stack[i] = stack[i - 1];
    --i;
  }
  stack[i] = val;
}

MatIdWeight stack_weight_pop(MatIdWeight* stack, uint32_t stack_sz)
{
  MatIdWeight val = stack[0];
  int i = 0;
  while(stack[i].id != 0xFFFFFFFF || i != stack_sz - 2)
  {
    stack[i] = stack[i + 1];
    ++i;
  }
  return val;
}

void stack_weight_push(MatIdWeight val, MatIdWeight* stack, uint32_t stack_sz)
{
  int i = stack_sz - 1;
  assert(stack[i].id == 0xFFFFFFFF);

  while(i > 0)
  {
    stack[i] = stack[i - 1];
    --i;
  }
  stack[i] = val;
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

  constexpr uint32_t stack_sz = 8;
  uint32_t material_stack[stack_sz] = {a_materialId, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
                                       0xFFFFFFFF,   0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF};

  while(material_stack[0] != 0xFFFFFFFF)
  {
    uint32_t currMatId = stack_pop(material_stack, stack_sz);
    assert(currMatId < m_materials.size());

    const float2 texCoordT = mulRows2x4(m_materials[currMatId].row0[0], m_materials[currMatId].row1[0], tc);
    const uint   mtype     = as_uint(m_materials[currMatId].data[UINT_MTYPE]);

    switch(mtype)
    {
      case MAT_TYPE_GLTF:
      {
        const float4 rands = rndFloat4_Pseudo(a_gen);
        
        const uint   texId     = as_uint(m_materials[currMatId].data[GLTF_UINT_TEXID0]);
        const float4 texColor  = (m_textures[texId]->sample(texCoordT));
        const float4 color     = (m_materials[currMatId].colors[GLTF_COLOR_BASE])*texColor;
        gltfSampleAndEval(m_materials.data() + currMatId, rands, v, n, tc, color, &res);
        break;
      }
      case MAT_TYPE_GLASS:
      {
        const float4 rands = rndFloat4_Pseudo(a_gen);

        glassSampleAndEval(m_materials.data() + currMatId, rands, v, n, tc, &res, a_misPrev);
        break;
      }
      case MAT_TYPE_CONDUCTOR:
      {
        const float4 rands = rndFloat4_Pseudo(a_gen);

        const uint   texId     = as_uint(m_materials[currMatId].data[CONDUCTOR_TEXID0]);
        const float3 alphaTex  = to_float3(m_textures[texId]->sample(texCoordT));
        
        const float2 alpha = float2(m_materials[currMatId].data[CONDUCTOR_ROUGH_V], m_materials[currMatId].data[CONDUCTOR_ROUGH_U]);
        if(trEffectivelySmooth(alpha))
          conductorSmoothSampleAndEval(m_materials.data() + currMatId, m_spectra.data(), wavelengths, rands, v, n, tc, &res);
        else
          conductorRoughSampleAndEval(m_materials.data() + currMatId, m_spectra.data(), wavelengths, rands, v, n, tc, alphaTex, &res);
        
        break;
      }
      case MAT_TYPE_DIFFUSE:
      {
        const float4 rands = rndFloat4_Pseudo(a_gen);

        const uint   texId       = as_uint(m_materials[currMatId].data[DIFFUSE_TEXID0]);
        // const float3 reflectance = to_float3(m_materials[a_materialId].colors[DIFFUSE_COLOR]); 
        const float4 texColor    = (m_textures[texId]->sample(texCoordT));
        const float4 color       = texColor;

        diffuseSampleAndEval(m_materials.data() + currMatId, m_spectra.data(), wavelengths, rands, v, n, tc, color, &res);

        break;
      }
      case MAT_TYPE_BLEND:
      {
        auto id = MaterialBlendSampleAndEval(currMatId, wavelengths, a_gen, v, n, tc, a_misPrev, &res);
        stack_push(id, material_stack, stack_sz);
        break;
      }
      default:
        break;
    }
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

  constexpr uint32_t stack_sz = 8;
  MatIdWeight material_stack[stack_sz] = {{a_materialId, 1.0f}, {0xFFFFFFFF, 0.0f}, {0xFFFFFFFF, 0.0f}, {0xFFFFFFFF, 0.0f},
                                          {0xFFFFFFFF, 0.0f},   {0xFFFFFFFF, 0.0f}, {0xFFFFFFFF, 0.0f}, {0xFFFFFFFF, 0.0f}};

  while(material_stack[0].id != 0xFFFFFFFF)
  {
    MatIdWeight currMat = stack_weight_pop(material_stack, stack_sz);
    assert(currMatId < m_materials.size());

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
      {
        glassEval(m_materials.data() + currMat.id, l, v, n, tc, {}, &currVal);

        res.val += currVal.val * currMat.weight;
        res.pdf += currVal.pdf * currMat.weight;
        break;
      }
      case MAT_TYPE_CONDUCTOR: 
      {
        const uint   texId     = as_uint(m_materials[currMat.id].data[CONDUCTOR_TEXID0]);
        const float3 alphaTex  = to_float3(m_textures[texId]->sample(texCoordT));

        const float2 alpha = float2(m_materials[currMat.id].data[CONDUCTOR_ROUGH_V], m_materials[currMat.id].data[CONDUCTOR_ROUGH_U]);
        if(trEffectivelySmooth(alpha))
          conductorSmoothEval(m_materials.data() + currMat.id, wavelengths, l, v, n, tc, &currVal);
        else
          conductorRoughEval(m_materials.data() + currMat.id, m_spectra.data(), wavelengths, l, v, n, tc, alphaTex, &currVal);

        res.val += currVal.val * currMat.weight;
        res.pdf += currVal.pdf * currMat.weight;
        break;
      }
      case MAT_TYPE_DIFFUSE:
      {
        const uint   texId       = as_uint(m_materials[currMat.id].data[DIFFUSE_TEXID0]);
        // const float3 reflectance = to_float3(m_materials[a_materialId].colors[DIFFUSE_COLOR]); 
        const float4 texColor    = (m_textures[texId]->sample(texCoordT));
        const float4 color       = texColor;

        diffuseEval(m_materials.data() + currMat.id, m_spectra.data(), wavelengths, l, v, n, tc, color, &currVal);

        res.val += currVal.val * currMat.weight;
        res.pdf += currVal.pdf * currMat.weight;
        break;
      }
      case MAT_TYPE_BLEND:
      {
        auto childMats = MaterialBlendEval(currMat, wavelengths, l, v, n, tc);
        stack_weight_push(childMats.second, material_stack, stack_sz);
        stack_weight_push(childMats.first, material_stack, stack_sz);
        break;
      }
      default:
        break;
    }
  }

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
}