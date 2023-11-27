#pragma once
#include "cglobals.h"
#include "crandom.h"
#include "cmaterial.h"


static float dielectricRoughEvalInternal(float3 wo, float3 wi, float3 wm, float2 alpha, float eta)
{
  if(wo.z * wi.z < 0) // not in the same hemisphere
  {
    return 0.0f;
  }

  float cosTheta_o = AbsCosTheta(wo);
  float cosTheta_i = AbsCosTheta(wi);
  if (cosTheta_i == 0 || cosTheta_o == 0)
    return 0.0f;

  float val = trD(wm, alpha) * trG(wo, wi, alpha) / (4.0f * cosTheta_i * cosTheta_o);

  return val;
}

static inline void plasticSampleAndEval(const Material* a_materials, float4 a_reflSpec, float4 rands_1, float rands_2,
                                        float3 v, float3 n, float2 tc, BsdfSample* pRes)
{
  const uint  cflags    = as_uint(a_materials[0].data[UINT_CFLAGS]);
  const float alpha     = a_materials[0].data[PLASTIC_ROUGHNESS];
  const float eta       = a_materials[0].data[PLASTIC_IOR_RATIO];
  const float fdr_int   = a_materials[0].data[PLASTIC_FDR_INTERIOR];
  const float spec_weight = a_materials[0].data[PLASTIC_SPEC_SAMPLE_WEIGHT];
  const uint  nonlinear   = as_uint(a_materials[0].data[PLASTIC_NONLINEAR]);


  float3 ggxDir;
  float  ggxPdf; 
  float  ggxVal;

  if(alpha == 0.0f) 
  {
    const float3 pefReflDir = reflect((-1.0f) * v, n);
    const float cosThetaOut = dot(pefReflDir, n);
    ggxDir                  = pefReflDir;
    ggxVal                  = (cosThetaOut <= 1e-6f) ? 0.0f : (1.0f / std::max(cosThetaOut, 1e-6f));
    ggxPdf                  = 1.0f;
  }
  else
  {
    const float2 alpha2 = float2(alpha, alpha);

    float3 nx, ny, nz = n;
    CoordinateSystem(nz, &nx, &ny);
    const float3 wo = float3(dot(v, nx), dot(v, ny), dot(v, nz));
    if(wo.z == 0)
      return;

    if(wo.z == 0)
      return;

    float3 wm = trSample(wo, float2(rands_1.x, rands_1.y), alpha2);
    float3 wi = reflect((-1.0f) * wo, wm);

    if(wo.z * wi.z < 0) // not in the same hemisphere
    {
      return;
    }

    ggxDir = normalize(wi.x * nx + wi.y * ny + wi.z * nz);
    ggxPdf = trPDF(wo, wm, alpha2) / (4.0f * std::abs(dot(wo, wm))); 
    ggxVal = dielectricRoughEvalInternal(wo, wi, wm, alpha2, eta);
  }

  const float3 lambertDir   = lambertSample(float2(rands_1.z, rands_1.w), v, n);
  const float  lambertPdf   = lambertEvalPDF(lambertDir, v, n);
  const float  lambertVal   = lambertEvalBSDF(lambertDir, v, n);

  const float f_i = FrDielectric(std::abs(dot(v,n)), eta); 
  
  float prob_specular = f_i * spec_weight;
  float prob_diffuse  = (1.0f - f_i) * (1.0f - spec_weight);
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

  if(rands_2 < prob_specular) // specular
  {
    pRes->dir       = ggxDir;
    pRes->val       = float4(ggxVal) * f_i;
    pRes->pdf       = ggxPdf * prob_specular;
    pRes->flags     = (alpha == 0.0f) ? RAY_EVENT_S : RAY_FLAG_HAS_NON_SPEC;
  } 
  else // lambert
  {
    pRes->dir       = lambertDir;
    pRes->val       = lambertVal * a_reflSpec;
    pRes->pdf       = lambertPdf * (1.0f - prob_specular);
    pRes->flags     = RAY_FLAG_HAS_NON_SPEC;
          

    const float f_o = FrDielectric(std::abs(dot(lambertDir, n)), eta);
    pRes->val      *= (1.0f - f_i) * (1.0f - f_o) / (eta * eta * (1.0f - fdr_int));
  }
            
}


static void plasticEval(const Material* a_materials, float4 a_reflSpec, float3 l, float3 v, float3 n, float2 tc, 
                        BsdfEval* pRes)
{
  const uint  cflags    = as_uint(a_materials[0].data[UINT_CFLAGS]);
  const float alpha     = a_materials[0].data[PLASTIC_ROUGHNESS];
  const float eta       = a_materials[0].data[PLASTIC_IOR_RATIO];
  const float fdr_int   = a_materials[0].data[PLASTIC_FDR_INTERIOR];
  const float spec_weight = a_materials[0].data[PLASTIC_SPEC_SAMPLE_WEIGHT];
  const uint  nonlinear   = as_uint(a_materials[0].data[PLASTIC_NONLINEAR]);

  float ggxVal, ggxPdf;
  if(alpha != 0.0f) 
  {
    const float2 alpha2 = float2(alpha, alpha);

    float3 nx, ny, nz = n;
    CoordinateSystem(nz, &nx, &ny);

    const float3 wo = float3(dot(v, nx), dot(v, ny), dot(v, nz));
    const float3 wi = float3(dot(l, nx), dot(l, ny), dot(l, nz));

    if(wo.z * wi.z < 0.0f)
      return;

    float3 wm = wo + wi;
    if (dot(wm, wm) == 0)
        return;

    wm = normalize(wm);

    ggxVal = dielectricRoughEvalInternal(wo, wi, wm, alpha2, eta);
    wm     = FaceForward(wm, float3(0.0f, 0.0f, 1.0f));
    ggxPdf = trPDF(wo, wm, alpha2) / (4.0f * std::abs(dot(wo, wm)));
  }
  else
  {
    ggxVal = 0.0f;
    ggxPdf = 0.0f;
  }

  float lambertVal       = lambertEvalBSDF(l, v, n);
  const float lambertPdf = lambertEvalPDF (l, v, n);
  float f_i              = 1.0f;
  float prob_diffuse     = 1.0f;
  float prob_specular    = 0.0f;
  float coeffLambertPdf  = 1.0f;

  f_i                 = FrDielectric(std::abs(dot(v, n)), eta);
  const float f_o     = FrDielectric(std::abs(dot(l, n)), eta);  
  const float coeff   = (1.f - f_i) * (1.f - f_o) / (eta * eta * (1.f - fdr_int));
  lambertVal         *= coeff;
  coeffLambertPdf     = coeff; 
  // prob_specular       = f_i * spec_weight;
  // prob_diffuse        = (1.f - f_i) * (1.f - spec_weight);
  
  // if(prob_diffuse != 0.0f && prob_specular != 0.0f)
  //   prob_diffuse = prob_diffuse / (prob_specular + prob_diffuse);
  // else
  // {
  //   prob_diffuse  = 1.0f;
  //   prob_specular = 0.0f;
  // }

  pRes->val = lambertVal * a_reflSpec + ggxVal * f_i; 
  pRes->pdf = lambertPdf * coeffLambertPdf + 2.0f * ggxPdf * f_i;
}