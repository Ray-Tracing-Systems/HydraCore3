#include "integrator_pt.h"
#include "include/crandom.h"

#include "include/cmaterial.h"
#include "include/cmat_gltf.h"
#include "include/cmat_glass.h"

#include <chrono>
#include <string>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

float Integrator::LightPdfSelectRev(int a_lightId) 
{ 
  return 1.0f;
}

float Integrator::LightEvalPDF(int a_lightId, float3 illuminationPoint, float3 ray_dir, const SurfaceHit* pSurfaceHit)
{
  const float3 lpos   = pSurfaceHit->pos;
  const float3 lnorm  = pSurfaceHit->norm;
  const float hitDist = length(illuminationPoint - lpos);
  const float pdfA    = 1.0f / (4.0f * m_lights[a_lightId].size.x * m_lights[a_lightId].size.y);
  const float cosVal  = std::max(dot(ray_dir, -1.0f*lnorm), 0.0f);
  return PdfAtoW(pdfA, hitDist, cosVal);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BsdfSample Integrator::MaterialSampleAndEval(int a_materialId, float4 rands, float3 v, float3 n, float2 tc)
{
  // implicit strategy

  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const float3 texColor  = to_float3(m_textures[ as_uint(m_materials[a_materialId].data[GLTF_UINT_TEXID0]) ]->sample(texCoordT));
  const float3 color     = to_float3(m_materials[a_materialId].colors[GLTF_COLOR_BASE])*texColor;
  const uint   mtype     = as_uint(m_materials[a_materialId].data[UINT_MTYPE]);

  // TODO: read other parameters from texture

  BsdfSample res;
  res.color = float3(0,0,0);
  res.pdf   = 1.0f;
  
  switch(mtype)
  {
    case MAT_TYPE_GLTF:
      gltfSampleAndEval(m_materials.data() + a_materialId, rands, v, n, tc, color, &res);
      break;
    case MAT_TYPE_GLASS:
      glassSampleAndEval(m_materials.data() + a_materialId, rands, v, n, tc, &res);
      break;
    default:
    break;
  }

  return res;
}

BsdfEval Integrator::MaterialEval(int a_materialId, float3 l, float3 v, float3 n, float2 tc)
{
  // explicit strategy

  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const float3 texColor  = to_float3(m_textures[ as_uint(m_materials[a_materialId].data[GLTF_UINT_TEXID0]) ]->sample(texCoordT));
  const float3 color     = to_float3(m_materials[a_materialId].colors[GLTF_COLOR_BASE])*texColor;
  const uint   mtype     = as_uint(m_materials[a_materialId].data[UINT_MTYPE]);

  // TODO: read other parameters from texture

  BsdfEval res;
  res.color = float3(0,0,0);
  res.pdf   = 0.0f;

  switch(mtype)
  {
    case MAT_TYPE_GLTF:
      gltfEval(m_materials.data() + a_materialId, l, v, n, tc, color, &res);
      break;
    case MAT_TYPE_GLASS:
      glassEval(m_materials.data() + a_materialId, l, v, n, tc, color, &res);
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
    if (idRemapFrom >= a_mId)
      high = mid - 1;
    else //if(a[mid]<i)
      low = mid + 1;
  }

  if (high+1 < offsAndSize.y)
  {
    const int idRemapFrom = m_allRemapLists[offsAndSize.x + (high + 1) * 2 + 0];
    const int idRemapTo   = m_allRemapLists[offsAndSize.x + (high + 1) * 2 + 1];
    res                   = (idRemapFrom == a_mId) ? uint(idRemapTo) : a_mId;
  }

  return res;
} 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Integrator::PackXYBlock(uint tidX, uint tidY, uint a_passNum)
{
  #pragma omp parallel for default(shared)
  for(int y=0;y<tidY;y++)
    for(int x=0;x<tidX;x++)
      PackXY(x, y);
}

void Integrator::CastSingleRayBlock(uint tid, uint* out_color, uint a_passNum)
{ 
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for(int i = 0; i < tid; i++)
    CastSingleRay((uint)i, out_color);
}

void Integrator::NaivePathTraceBlock(uint tid, float4* out_color, uint a_passNum)
{
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for(int i = 0; i < tid; i++)
    for(int j = 0; j < a_passNum; j++)
      NaivePathTrace((uint)i, out_color);
  naivePtTime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count()/1000.f;
}

void Integrator::PathTraceBlock(uint tid, float4* out_color, uint a_passNum)
{
  auto start               = std::chrono::high_resolution_clock::now();
  auto start2              = std::chrono::high_resolution_clock::now();
  
  const int allCountSample = tid * a_passNum;
  int countSample          = 0;

  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for (int i = 0; i < tid; i++)
  {

    for (int j = 0; j < a_passNum; j++)
    {
      PathTrace((uint)i, out_color);
      countSample++;
    }
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count() / 1000.f;

    if (duration > 2)
    {
      std::cout << "progress: " << (float)countSample / (float)allCountSample * 100.0f << "                    \r";
      start2 = std::chrono::high_resolution_clock::now();
    }
  }

  std::cout << std::endl;
  shadowPtTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count()/1000.f;
}

void Integrator::RayTraceBlock(uint tid, float4* out_color, uint a_passNum)
{
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for(int i = 0; i < tid; i++)
    RayTrace((uint)i, out_color);
  raytraceTime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count()/1000.f;
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
