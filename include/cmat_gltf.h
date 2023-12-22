#pragma once
#include "cglobals.h"
#include "crandom.h"
#include "cmaterial.h"


static inline void gltfSampleAndEval(const Material* a_materials, float4 rands, float3 v, 
                                     float3 n, float2 tc, float4 color, BsdfSample* pRes)
{
  // PLEASE! use 'a_materials[0].' for a while ... , not a_materials-> and not *(a_materials).
  const uint   cflags   = as_uint(a_materials[0].data[UINT_CFLAGS]);
  const float4 specular = a_materials[0].colors[GLTF_COLOR_METAL]; 
  const float4 coat     = a_materials[0].colors[GLTF_COLOR_COAT];  
  const float  roughness  = clamp(1.0f - a_materials[0].data[GLTF_FLOAT_GLOSINESS], 0.0f, 1.0f);   
  float        alpha      = a_materials[0].data[GLTF_FLOAT_ALPHA];                 
  const float  fresnelIOR = a_materials[0].data[GLTF_FLOAT_IOR];
  
  if(cflags == GLTF_COMPONENT_METAL) // assume only GGX-based metal component set
    alpha = 1.0f;

  float3 ggxDir;
  float  ggxPdf; 
  float  ggxVal;

  if(roughness == 0.0f) // perfect specular reflection in coating or metal layer
  {
    const float3 pefReflDir = reflect((-1.0f) * v, n);
    const float cosThetaOut = dot(pefReflDir, n);
    ggxDir                  = pefReflDir;
    ggxVal                  = (cosThetaOut <= 1e-6f) ? 0.0f : (1.0f/std::max(cosThetaOut, 1e-6f));  // BSDF is multiplied (outside) by cosThetaOut. For mirrors this shouldn't be done, so we pre-divide here instead.
    ggxPdf                  = 1.0f;
  }
  else
  {
    ggxDir                  = ggxSample(float2(rands.x, rands.y), v, n, roughness);
    ggxPdf                  = ggxEvalPDF (ggxDir, v, n, roughness); 
    ggxVal                  = ggxEvalBSDF(ggxDir, v, n, roughness);
  }

  const float3 lambertDir   = lambertSample(float2(rands.x, rands.y), v, n);
  const float  lambertPdf   = lambertEvalPDF(lambertDir, v, n);
  const float  lambertVal   = lambertEvalBSDF(lambertDir, v, n);

  // (1) select between metal and dielectric via rands.z
  //
  float pdfSelect = 1.0f;
  if(rands.z < alpha) // select metall
  {
    pdfSelect         *= alpha;
    const float  VdotH = dot(v,normalize(v + ggxDir));
    pRes->dir          = ggxDir;
    pRes->val          = ggxVal * alpha * hydraFresnelCond(specular, VdotH, fresnelIOR, roughness); //TODO: disable fresnel here for mirrors
    pRes->pdf          = ggxPdf;
    pRes->flags        = (roughness == 0.0f) ? RAY_EVENT_S : RAY_FLAG_HAS_NON_SPEC;
  }
  else                // select dielectric
  {
    pdfSelect *= 1.0f - alpha;
    
    // (2) now select between specular and diffise via rands.w
    //
    const float f_i = FrDielectricPBRT(std::abs(dot(v,n)), 1.0f, fresnelIOR); 
    const float m_specular_sampling_weight = a_materials[0].data[GLTF_FLOAT_MI_SSW];
    
    float prob_specular = f_i * m_specular_sampling_weight;
    float prob_diffuse  = (1.0f - f_i) * (1.0f - m_specular_sampling_weight);
    if(prob_diffuse != 0.0f && prob_specular != 0.0f)
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
      pdfSelect      *= choicePdf;
      pRes->dir       = ggxDir;
      pRes->val       = ggxVal*coat*(1.0f - alpha)*f_i;
      pRes->pdf       = ggxPdf;
      pRes->flags     = (roughness == 0.0f) ? RAY_EVENT_S : RAY_FLAG_HAS_NON_SPEC;
    } 
    else
    {
      pdfSelect      *= (1.0f-choicePdf); // lambert
      pRes->dir       = lambertDir;
      pRes->val       = lambertVal * color * (1.0f - alpha);
      pRes->pdf       = lambertPdf;
      pRes->flags     = RAY_FLAG_HAS_NON_SPEC;
            
      if ((cflags & GLTF_COMPONENT_ORENNAYAR) != 0)
        pRes->val *= orennayarFunc(lambertDir, (-1.0f) * v, n, a_materials[0].data[GLTF_FLOAT_ROUGH_ORENNAYAR]);
            
      if((cflags & GLTF_COMPONENT_COAT) != 0 && (cflags & GLTF_COMPONENT_LAMBERT) != 0) // Plastic, account for retroreflection between surface and coating layer
      {
        const float m_fdr_int = a_materials[0].data[GLTF_FLOAT_MI_FDR_INT];
        const float f_o       = FrDielectricPBRT(std::abs(dot(lambertDir, n)), 1.0f, fresnelIOR);
        pRes->val          *= (1.0f - f_i) * (1.0f - f_o) / (fresnelIOR * fresnelIOR * (1.0f - m_fdr_int));
      }
    }
  }   
  pRes->pdf *= pdfSelect;
}


static void gltfEval(const Material* a_materials, float3 l, float3 v, float3 n, float2 tc, 
                     float4 color, BsdfEval* res)
{
  const uint   cflags     = as_uint(a_materials[0].data[UINT_CFLAGS]);
  const float4 specular   = a_materials[0].colors[GLTF_COLOR_METAL];
  const float4 coat       = a_materials[0].colors[GLTF_COLOR_COAT];
  const float  roughness  = clamp(1.0f - a_materials[0].data[GLTF_FLOAT_GLOSINESS], 0.0f, 1.0f);
        float  alpha      = a_materials[0].data[GLTF_FLOAT_ALPHA];
  const float  fresnelIOR = a_materials[0].data[GLTF_FLOAT_IOR];

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
  float f_i              = 1.0f;
  float prob_diffuse     = 1.0f;
  float prob_specular    = 0.0f;
  float coeffLambertPdf  = 1.0f;

  if ((cflags & GLTF_COMPONENT_ORENNAYAR) != 0)
    lambertVal *= orennayarFunc(l, v, n, a_materials[0].data[GLTF_FLOAT_ROUGH_ORENNAYAR]);
      
  if((cflags & GLTF_COMPONENT_COAT) != 0 && (cflags & GLTF_COMPONENT_LAMBERT) != 0) // Plastic, account for retroreflection between surface and coating layer
  {
    f_i                                    = FrDielectricPBRT(std::abs(dot(v,n)), 1.0f, fresnelIOR);
    const float f_o                        = FrDielectricPBRT(std::abs(dot(l,n)), 1.0f, fresnelIOR);  
    const float m_fdr_int                  = a_materials[0].data[GLTF_FLOAT_MI_FDR_INT];
    const float coeff                      = (1.f - f_i) * (1.f - f_o) / (fresnelIOR*fresnelIOR*(1.f - m_fdr_int));
    lambertVal                            *= coeff;
    coeffLambertPdf                        = coeff; 
    const float m_specular_sampling_weight = a_materials[0].data[GLTF_FLOAT_MI_SSW];
    prob_specular                          = f_i * m_specular_sampling_weight;
    prob_diffuse                           = (1.f - f_i) * (1.f - m_specular_sampling_weight);
    
    if(prob_diffuse != 0.0f && prob_specular != 0.0f)
      prob_diffuse = prob_diffuse / (prob_specular + prob_diffuse);
    else
    {
      prob_diffuse  = 1.0f;
      prob_specular = 0.0f;
    }
  }

  const float4 fConductor    = hydraFresnelCond(specular, VdotH, fresnelIOR, roughness); // (1) eval metal component      
  const float4 specularColor = ggxVal*fConductor;                                        // eval metal specular component
      
  float  dielectricPdf = lambertPdf * coeffLambertPdf; 
  if((cflags & GLTF_COMPONENT_COAT) != 0)
    dielectricPdf += 2.0f * ggxPdf * f_i; 
                                
  const float4 dielectricVal = lambertVal * color + ggxVal * coat * f_i;

  res->val = alpha * specularColor + (1.0f - alpha) * dielectricVal; // (3) accumulate final color and pdf
  res->pdf = alpha * ggxPdf        + (1.0f - alpha) * dielectricPdf; // (3) accumulate final color and pdf
}