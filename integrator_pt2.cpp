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
      const float  lradius = m_lights[a_lightId].size.x;
      const float3 lcenter = to_float3(m_lights[a_lightId].pos);
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


BsdfSample Integrator::MaterialBlendSampleAndEval(uint a_materialId, float3 wavelengths, RandomGen* a_gen, float3 v, float3 n, float2 tc, 
                                                  MisData* a_misPrev, const uint a_currRayFlags)
{
  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const uint   texId     = as_uint(m_materials[a_materialId].data[BLEND_TEXID0]);
  const float  weightTex = to_float3(m_textures[texId]->sample(texCoordT)).x;
  const float  weight    = m_materials[a_materialId].data[BLEND_WEIGHT] * weightTex;

  const uint matId1 = as_uint(m_materials[a_materialId].data[BLEND_MAT_ID_1]);
  const uint matId2 = as_uint(m_materials[a_materialId].data[BLEND_MAT_ID_2]);

  BsdfSample res;
  {
    res.val   = float3(0, 0, 0);
    res.pdf   = 1.0f;
    res.dir   = float3(0,1,0);
    res.flags = a_currRayFlags;
  }

  const float select = rndFloat1_Pseudo(a_gen);
  if(select < weight)
  {
    res = MaterialSampleAndEval(matId2, wavelengths, a_gen, v, n, tc, a_misPrev, a_currRayFlags);
    res.pdf *= weight;
  }
  else
  {
    res = MaterialSampleAndEval(matId1, wavelengths, a_gen, v, n, tc, a_misPrev, a_currRayFlags);
    res.pdf *= 1.0f - weight;
  }


  return res;
}

BsdfEval Integrator::MaterialBlendEval(uint a_materialId, float3 wavelengths, float3 l, float3 v, float3 n, float2 tc)
{
  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const uint   texId     = as_uint(m_materials[a_materialId].data[BLEND_TEXID0]);
  const float  weightTex = to_float3(m_textures[texId]->sample(texCoordT)).x;
  const float  weight    = m_materials[a_materialId].data[BLEND_WEIGHT] * weightTex;

  const uint matId1 = as_uint(m_materials[a_materialId].data[BLEND_MAT_ID_1]);
  const uint matId2 = as_uint(m_materials[a_materialId].data[BLEND_MAT_ID_2]);

  BsdfEval res;
  {
    res.val = float3(0, 0, 0);
    res.pdf = 0.0f;
  }

  BsdfEval res1 = MaterialEval(matId1, wavelengths, l, v, n, tc);
  BsdfEval res2 = MaterialEval(matId2, wavelengths, l, v, n, tc);
  res.pdf = weight * res2.pdf + (1.0f - weight) * res1.pdf;
  res.val = weight * res2.val + (1.0f - weight) * res1.val;

  return res;
}

BsdfSample Integrator::MaterialSampleAndEval(uint a_materialId, float3 wavelengths, RandomGen* a_gen, float3 v, float3 n, float2 tc, 
                                             MisData* a_misPrev, const uint a_currRayFlags)
{
  // implicit strategy
  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const uint   mtype     = as_uint(m_materials[a_materialId].data[UINT_MTYPE]);

  // TODO: read other parameters from texture

  BsdfSample res;
  {
    res.val   = float3(0,0,0);
    res.pdf   = 1.0f;
    res.dir   = float3(0,1,0);
    res.flags = a_currRayFlags;
  }

  switch(mtype)
  {
    case MAT_TYPE_GLTF:
    {
      const float4 rands = rndFloat4_Pseudo(a_gen);
      
      const uint   texId     = as_uint(m_materials[a_materialId].data[GLTF_UINT_TEXID0]);
      const float3 texColor  = to_float3(m_textures[texId]->sample(texCoordT));
      const float3 color     = to_float3(m_materials[a_materialId].colors[GLTF_COLOR_BASE])*texColor;
      gltfSampleAndEval(m_materials.data() + a_materialId, rands, v, n, tc, color, &res);
      break;
    }
    case MAT_TYPE_GLASS:
    {
      const float4 rands = rndFloat4_Pseudo(a_gen);

      glassSampleAndEval(m_materials.data() + a_materialId, rands, v, n, tc, &res, a_misPrev);
      break;
    }
    case MAT_TYPE_CONDUCTOR:
    {
      const float4 rands = rndFloat4_Pseudo(a_gen);

      const uint   texId     = as_uint(m_materials[a_materialId].data[CONDUCTOR_TEXID0]);
      const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
      const float3 alphaTex  = to_float3(m_textures[texId]->sample(texCoordT));
      
      const float2 alpha = float2(m_materials[a_materialId].data[CONDUCTOR_ROUGH_V], m_materials[a_materialId].data[CONDUCTOR_ROUGH_U]);
      if(trEffectivelySmooth(alpha))
        conductorSmoothSampleAndEval(m_materials.data() + a_materialId, m_spectra.data(), wavelengths, rands, v, n, tc, &res);
      else
        conductorRoughSampleAndEval(m_materials.data() + a_materialId, m_spectra.data(), wavelengths, rands, v, n, tc, alphaTex, &res);
      
      break;
    }
    case MAT_TYPE_DIFFUSE:
    {
      const float4 rands = rndFloat4_Pseudo(a_gen);

      const uint   texId       = as_uint(m_materials[a_materialId].data[DIFFUSE_TEXID0]);
      // const float3 reflectance = to_float3(m_materials[a_materialId].colors[DIFFUSE_COLOR]); 
      const float3 texColor    = to_float3(m_textures[texId]->sample(texCoordT));
      const float3 color       = texColor;

      diffuseSampleAndEval(m_materials.data() + a_materialId, m_spectra.data(), wavelengths, rands, v, n, tc, color, &res);

      break;
    }
    case MAT_TYPE_BLEND:
    {
      res = MaterialBlendSampleAndEval(a_materialId, wavelengths, a_gen, v, n, tc, a_misPrev, a_currRayFlags);
      break;
    }
    default:
      break;
  }

  return res;
}

BsdfEval Integrator::MaterialEval(uint a_materialId, float3 wavelengths, float3 l, float3 v, float3 n, float2 tc)
{
  // explicit strategy
  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const uint   mtype     = as_uint(m_materials[a_materialId].data[UINT_MTYPE]);

  // TODO: read other parameters from texture
  BsdfEval res;
  {
    res.val = float3(0,0,0);
    res.pdf   = 0.0f;
  }

  switch(mtype)
  {
    case MAT_TYPE_GLTF:
    {
      const uint   texId     = as_uint(m_materials[a_materialId].data[GLTF_UINT_TEXID0]);
      const float3 texColor  = to_float3(m_textures[texId]->sample(texCoordT));
      const float3 color     = to_float3(m_materials[a_materialId].colors[GLTF_COLOR_BASE])*texColor;
      gltfEval(m_materials.data() + a_materialId, l, v, n, tc, color, &res);
      break;
    }
    case MAT_TYPE_GLASS:
    {
      glassEval(m_materials.data() + a_materialId, l, v, n, tc, {}, &res);
      break;
    }
    case MAT_TYPE_CONDUCTOR: 
    {
      const uint   texId     = as_uint(m_materials[a_materialId].data[CONDUCTOR_TEXID0]);
      const float3 alphaTex  = to_float3(m_textures[texId]->sample(texCoordT));

      const float2 alpha = float2(m_materials[a_materialId].data[CONDUCTOR_ROUGH_V], m_materials[a_materialId].data[CONDUCTOR_ROUGH_U]);
      if(trEffectivelySmooth(alpha))
        conductorSmoothEval(m_materials.data() + a_materialId, wavelengths, l, v, n, tc, &res);
      else
        conductorRoughEval(m_materials.data() + a_materialId, m_spectra.data(), wavelengths, l, v, n, tc, alphaTex, &res);

      break;
    }
    case MAT_TYPE_DIFFUSE:
    {
      const uint   texId       = as_uint(m_materials[a_materialId].data[DIFFUSE_TEXID0]);
      // const float3 reflectance = to_float3(m_materials[a_materialId].colors[DIFFUSE_COLOR]); 
      const float3 texColor    = to_float3(m_textures[texId]->sample(texCoordT));
      const float3 color       = texColor;

      diffuseEval(m_materials.data() + a_materialId, m_spectra.data(), wavelengths, l, v, n, tc, color, &res);

      break;
    }
    case MAT_TYPE_BLEND:
    {
      res = MaterialBlendEval(a_materialId, wavelengths, l, v, n, tc);
      break;
    }
    default:
      break;
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