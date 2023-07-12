#include "integrator_pt.h"
#include "include/crandom.h"

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
  const float pdfA    = 1.0f / (4.0f*m_light.size.x*m_light.size.y);
  const float cosVal  = std::max(dot(ray_dir, -1.0f*lnorm), 0.0f);
  return PdfAtoW(pdfA, hitDist, cosVal);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BsdfSample Integrator::MaterialSampleAndEval(int a_materialId, float4 rands, float3 v, float3 n, float2 tc)
{
  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const float3 texColor  = to_float3(m_textures[ m_materials[a_materialId].texId[0] ]->sample(texCoordT));

        uint   mtype     = m_materials[a_materialId].mtype;
  const uint   cflags    = m_materials[a_materialId].cflags;

  const float3 color      = to_float3(m_materials[a_materialId].baseColor)*texColor;
  const float3 specular   = to_float3(m_materials[a_materialId].metalColor);
  const float3 coat       = to_float3(m_materials[a_materialId].coatColor);
  const float  roughness  = 1.0f - m_materials[a_materialId].glosiness;
  float        alpha      = m_materials[a_materialId].alpha;
  const float  fresnelIOR = m_materials[a_materialId].ior;

  // TODO: read color     from texture
  // TODO: read roughness from texture
  // TODO: read alpha     from texture

  BsdfSample res;
  res.color = float3(0,0,0);
  res.pdf   = 1.0f;
  
  switch(mtype)
  {
    case MAT_TYPE_GLTF:
    {
      if(cflags == GLTF_COMPONENT_METAL) // assume only GGX-based metal component set
        alpha = 1.0f;

      float3 ggxDir;
      float  ggxPdf; 
      float  ggxVal;

      if(roughness == 0.0f) // perfect specular reflection in coating or metal layer
      {
        const float3 pefReflDir = reflect((-1.0f)*v, n);
        const float cosThetaOut = dot(pefReflDir, n);
        ggxDir = pefReflDir;
        ggxVal = (cosThetaOut <= 1e-6f) ? 0.0f : (1.0f/std::max(cosThetaOut, 1e-6f));  // BSDF is multiplied (outside) by cosThetaOut. For mirrors this shouldn't be done, so we pre-divide here instead.
        ggxPdf = 1.0f;
      }
      else
      {
        ggxDir = ggxSample(float2(rands.x, rands.y), v, n, roughness);
        ggxPdf = ggxEvalPDF (ggxDir, v, n, roughness); 
        ggxVal = ggxEvalBSDF(ggxDir, v, n, roughness);
      }

      const float3 lambertDir = lambertSample(float2(rands.x, rands.y), v, n);
      const float  lambertPdf = lambertEvalPDF(lambertDir, v, n);
      const float  lambertVal = lambertEvalBSDF(lambertDir, v, n);

      float VdotH = dot(v,normalize(v + ggxDir));

      // (1) select between metal and dielectric via rands.z
      //
      float pdfSelect = 1.0f;
      if(rands.z < alpha) // select metall
      {
        pdfSelect *= alpha;
        res.direction = ggxDir;
        res.color     = ggxVal*alpha*hydraFresnelCond(specular, VdotH, fresnelIOR, roughness); //TODO: disable fresnel here for mirrors
        res.pdf       = ggxPdf;
        res.flags     = (roughness == 0.0f) ? RAY_EVENT_S : RAY_FLAG_HAS_NON_SPEC;
      }
      else                // select dielectric
      {
        pdfSelect *= 1.0f - alpha;
        
        // (2) now select between specular and diffise via rands.w
        //
        const float f_i = FrDielectricPBRT(std::abs(dot(v,n)), 1.0f, fresnelIOR); 
        const float m_specular_sampling_weight = m_materials[a_materialId].data[MI_SSW];
        
        float prob_specular = f_i * m_specular_sampling_weight;
        float prob_diffuse  = (1.f - f_i) * (1.f - m_specular_sampling_weight);
        if(prob_diffuse != 0.0f && prob_diffuse != 0.0f)
        {
          prob_specular = prob_specular / (prob_specular + prob_diffuse);
          prob_diffuse  = 1.f - prob_specular;
        }
        else
        {
          prob_diffuse  = 1.0f;
          prob_specular = 0.0f;
        }

        float choicePdf = ((cflags & GLTF_COMPONENT_COAT) == 0) ? 0.0f : prob_specular; // if don't have coal layer, never select it
        if(rands.w < prob_specular) // specular
        {
          pdfSelect *= choicePdf;
          res.direction = ggxDir;
          res.color     = ggxVal*coat*(1.0f - alpha)*f_i;
          res.pdf       = ggxPdf;
          res.flags     = (roughness == 0.0f) ? RAY_EVENT_S : RAY_FLAG_HAS_NON_SPEC;
        } 
        else
        {
          pdfSelect *= (1.0f-choicePdf); // lambert
          res.direction = lambertDir;
          res.color     = lambertVal*color*(1.0f - alpha);
          res.pdf       = lambertPdf;
          res.flags     = RAY_FLAG_HAS_NON_SPEC;
          
          if((cflags & GLTF_COMPONENT_COAT) != 0 && (cflags & GLTF_COMPONENT_LAMBERT) != 0) // Plastic, account for retroreflection between surface and coating layer
          {
            const float m_fdr_int = m_materials[a_materialId].data[MI_FDR_INT];
            const float f_o = FrDielectricPBRT(std::abs(dot(lambertDir,n)), 1.0f, fresnelIOR); 
            res.color *= (1.f - f_i) * (1.f - f_o) / (fresnelIOR*fresnelIOR*(1.f - m_fdr_int));
          }
        }
      }
      
      res.pdf *= pdfSelect;
    }
    break;

    default:
    break;
  }

  return res;
}

BsdfEval Integrator::MaterialEval(int a_materialId, float3 l, float3 v, float3 n, float2 tc)
{
  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const float3 texColor  = to_float3(m_textures[ m_materials[a_materialId].texId[0] ]->sample(texCoordT));

        uint mtype       = m_materials[a_materialId].mtype;
  const uint cflags      = m_materials[a_materialId].cflags;

  const float3 color     = to_float3(m_materials[a_materialId].baseColor)*texColor;
  const float3 specular  = to_float3(m_materials[a_materialId].metalColor);
  const float3 coat      = to_float3(m_materials[a_materialId].coatColor);
  const float roughness  = 1.0f - m_materials[a_materialId].glosiness;
        float alpha      = m_materials[a_materialId].alpha;
  const float fresnelIOR = m_materials[a_materialId].ior;

  // TODO: read color     from texture
  // TODO: read roughness from texture
  // TODO: read alpha     from texture

  BsdfEval res;
  res.color = float3(0,0,0);
  res.pdf   = 0.0f;

  switch(mtype)
  {
    case MAT_TYPE_GLTF:
    {
      if(cflags == GLTF_COMPONENT_METAL) // assume only GGX-based metal
        alpha = 1.0f;
      
      float ggxVal, ggxPdf, VdotH; 
      if(roughness != 0.0f) // perfect specular reflection in coating layer
      {
        ggxVal = ggxEvalBSDF(l, v, n, roughness);
        ggxPdf = ggxEvalPDF (l, v, n, roughness);
        VdotH  = dot(v,normalize(v + l));
      }
      else
      {
        ggxVal = 0.0f;
        ggxPdf = 0.0f;
        VdotH  = dot(v,n);
      }

      float lambertVal       = lambertEvalBSDF(l, v, n);
      const float lambertPdf = lambertEvalPDF (l, v, n); 

      float f_i = 1.0f;
      float prob_diffuse  = 1.0f;
      float prob_specular = 0.0f;
      float coeffLambertPdf = 1.0f;
      
      if((cflags & GLTF_COMPONENT_COAT) != 0 && (cflags & GLTF_COMPONENT_LAMBERT) != 0) // Plastic, account for retroreflection between surface and coating layer
      {
        f_i = FrDielectricPBRT(std::abs(dot(v,n)), 1.0f, fresnelIOR);
        const float f_o = FrDielectricPBRT(std::abs(dot(l,n)), 1.0f, fresnelIOR);  
        const float m_fdr_int = m_materials[a_materialId].data[MI_FDR_INT];
        const float coeff = (1.f - f_i) * (1.f - f_o) / (fresnelIOR*fresnelIOR*(1.f - m_fdr_int));
        lambertVal *= coeff;
        coeffLambertPdf = coeff; 

        const float m_specular_sampling_weight = m_materials[a_materialId].data[MI_SSW];
        prob_specular = f_i * m_specular_sampling_weight;
        prob_diffuse  = (1.f - f_i) * (1.f - m_specular_sampling_weight);
        if(prob_diffuse != 0.0f && prob_diffuse != 0.0f)
          prob_diffuse = prob_diffuse / (prob_specular + prob_diffuse);
        else
        {
          prob_diffuse  = 1.0f;
          prob_specular = 0.0f;
        }
      }

      const float3 fConductor    = hydraFresnelCond(specular, VdotH, fresnelIOR, roughness); // (1) eval metal component      
      const float3 specularColor = ggxVal*fConductor;                                        // eval metal specular component
      
      const float  dielectricPdf = lambertPdf*coeffLambertPdf + 2.0f*ggxPdf*f_i;
      const float3 dielectricVal = lambertVal*color + ggxVal*coat*f_i;

      res.color = alpha*specularColor + (1.0f - alpha)*dielectricVal; // (3) accumulate final color and pdf
      res.pdf   = alpha*ggxPdf        + (1.0f - alpha)*dielectricPdf; // (3) accumulate final color and pdf
    }
    break;

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
  for(uint i=0;i<tid;i++)
    CastSingleRay(i, out_color);
}

void Integrator::NaivePathTraceBlock(uint tid, float4* out_color, uint a_passNum)
{
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for(uint i=0;i<tid;i++)
    for(int j=0;j<a_passNum;j++)
      NaivePathTrace(i, out_color);
  naivePtTime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count()/1000.f;
}

void Integrator::PathTraceBlock(uint tid, float4* out_color, uint a_passNum)
{
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for(uint i=0;i<tid;i++)
    for(int j=0;j<a_passNum;j++)
      PathTrace(i, out_color);
  shadowPtTime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count()/1000.f;
}

void Integrator::RayTraceBlock(uint tid, float4* out_color, uint a_passNum)
{
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for(uint i=0;i<tid;i++)
    for(int j=0;j<a_passNum;j++)
      RayTrace(i, out_color);
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
