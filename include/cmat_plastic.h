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


static float4 plasticEvalInternal(float3 wo, float3 wi, float3 wm, float2 alpha, float eta, float4 reflSpec, bool nonlinear,
                                 float t_o, float t_i, float internal_refl)
{
  if(wo.z * wi.z < 0) // not in the same hemisphere
  {
    return float4{0.0f};
  }

  float cos_theta_o = wo.z;
  float cos_theta_i = wi.z;
  // if (cos_theta_o == 0 || cos_theta_i == 0)
  //   return 0.0f;

  float3 H = normalize(wo + wi);
  const float F = FrDielectric(std::abs(dot(wi, wm)), eta);

  float val = F * trD(wm, alpha) * trG(wi, wo, alpha) / (4.0f * cos_theta_i * cos_theta_o);

  float4 diff = reflSpec;
  diff /= 1.f - (nonlinear > 0 ? (diff * internal_refl) : float4(internal_refl));
  const float inv_eta_2 = 1.f / (eta * eta);

  return diff * (M_1_PI * inv_eta_2 * cos_theta_o * t_i * t_o) + float4(val);
}

// static inline void plasticSampleAndEval(const Material* a_materials, float4 a_reflSpec, float4 rands_1, float rands_2,
//                                         float3 v, float3 n, float2 tc, BsdfSample* pRes, const float* transmittance)
// {
//   const uint  cflags    = as_uint(a_materials[0].data[UINT_CFLAGS]);
//   const float alpha     = a_materials[0].data[PLASTIC_ROUGHNESS];
//   const float eta       = a_materials[0].data[PLASTIC_IOR_RATIO];
//   const float spec_weight = a_materials[0].data[PLASTIC_SPEC_SAMPLE_WEIGHT];
//   const uint  nonlinear   = as_uint(a_materials[0].data[PLASTIC_NONLINEAR]);
//   const float internal_refl = a_materials[0].data[PLASTIC_PRECOMP_REFLECTANCE];
//   const float2 alpha2 = {alpha, alpha};

//   float3 nx, ny, nz = n;
//   CoordinateSystem(nz, &nx, &ny);
//   const float3 wi = float3(dot(v, nx), dot(v, ny), dot(v, nz));
//   if(wi.z == 0)
//     return;


//   const float cos_theta_i = wi.z;

//   const float t_i = lerp_gather(transmittance, cos_theta_i, MI_ROUGH_TRANSMITTANCE_RES);

//   float prob_specular = (1.f - t_i) * spec_weight;
//   float prob_diffuse  = t_i * (1.f - spec_weight);

//   // if (has_specular != has_diffuse)
//   //   prob_specular = has_specular ? 1.f : 0.f;
//   // else
//   prob_specular = prob_specular / (prob_specular + prob_diffuse);
//   prob_diffuse = 1.f - prob_specular;

//   const bool sample_specular = rands_2 < prob_specular;
//   const bool sample_diffuse  = !sample_specular;

//   float3 wo {0.0f, 0.0f, 0.0f};
//   float3 wm {0.0f, 0.0f, 0.0f};
//   if(sample_specular)
//   {
//     wo = reflect((-1.0f) * wi, wm);
//     wm = trSample(wo, float2(rands_1.x, rands_1.y), alpha2);  
//     if(wo.z * wi.z < 0) // not in the same hemisphere
//     {
//       return;
//     }

//     pRes->dir   = normalize(wi.x * nx + wi.y * ny + wi.z * nz);
//     pRes->flags = (alpha == 0.0f) ? RAY_EVENT_S : RAY_FLAG_HAS_NON_SPEC;
//   }

//   float lambert_pdf = 0.f;
//   if(sample_diffuse)
//   {
//     wo = lambertSample(float2(rands_1.z, rands_1.w), v, n);
//     wm = n;
//     lambert_pdf = lambertEvalPDF(wo, v, n);
//     pRes->dir   = wo;
//     pRes->flags = RAY_FLAG_HAS_NON_SPEC;

//     wo = float3(dot(wo, nx), dot(wo, ny), dot(wo, nz));
//   }

//   float cos_theta_o = wo.z;

//   if(wo.z * wi.z < 0) // not in the same hemisphere
//   {
//     return;
//   }

//   // float3 H = normalize(wo + wi);
//   const float t_o = lerp_gather(transmittance, cos_theta_o, MI_ROUGH_TRANSMITTANCE_RES);

//   float pdf = trPDF(wo, wm, alpha2) / (4.0f * std::abs(dot(wo, wm)));
//   pdf *= prob_specular;
//   pdf += prob_diffuse * lambert_pdf;

//   pRes->pdf = pdf;
//   pRes->val = plasticEvalInternal(wo, wi, wm, alpha2, eta, a_reflSpec, nonlinear,
//                                   t_o, t_i, internal_refl);
// }

// static inline void plasticSampleAndEval(const Material* a_materials, float4 a_reflSpec, float4 rands_1, float rands_2,
//                                         float3 v, float3 n, float2 tc, BsdfSample* pRes, const float* transmittance)
// {
//   const uint  cflags    = as_uint(a_materials[0].data[UINT_CFLAGS]);
//   const float alpha     = a_materials[0].data[PLASTIC_ROUGHNESS];
//   const float eta       = a_materials[0].data[PLASTIC_IOR_RATIO];
//   const float spec_weight = a_materials[0].data[PLASTIC_SPEC_SAMPLE_WEIGHT];
//   const uint  nonlinear   = as_uint(a_materials[0].data[PLASTIC_NONLINEAR]);

//   const float internal_refl = a_materials[0].data[PLASTIC_PRECOMP_REFLECTANCE];

//   float3 ggxDir;
//   float  ggxPdf; 
//   float  ggxVal;

//   float cosTheta_i = 0.0f;
//   float cosTheta_o = 0.0f;
//   if(alpha == 0.0f) 
//   {
//     const float3 pefReflDir = reflect((-1.0f) * v, n);
//     const float cosThetaOut = dot(pefReflDir, n);
//     ggxDir                  = pefReflDir;
//     ggxVal                  = (cosThetaOut <= 1e-6f) ? 0.0f : (1.0f / std::max(cosThetaOut, 1e-6f));
//     ggxPdf                  = 1.0f;
//     cosTheta_i = cosThetaOut;
//     cosTheta_o = dot(v, n);
//   }
//   else
//   {
//     const float2 alpha2 = float2(alpha, alpha);

//     float3 nx, ny, nz = n;
//     CoordinateSystem(nz, &nx, &ny);
//     const float3 wo = float3(dot(v, nx), dot(v, ny), dot(v, nz));
//     if(wo.z == 0)
//       return;
//     if(wo.z == 0)
//       return;

//     float3 wm = trSample(wo, float2(rands_1.x, rands_1.y), alpha2);
//     float3 wi = reflect((-1.0f) * wo, wm);

//     if(wo.z * wi.z < 0) // not in the same hemisphere
//     {
//       return;
//     }

//     ggxDir = normalize(wi.x * nx + wi.y * ny + wi.z * nz);
//     ggxPdf = trPDF(wo, wm, alpha2) / (4.0f * std::abs(dot(wo, wm))); 
//     const float f = FrDielectric(std::abs(dot(wo, wm)), eta); 
//     ggxVal = f * dielectricRoughEvalInternal(wo, wi, wm, alpha2, eta);
//     cosTheta_i = AbsCosTheta(wi);
//     cosTheta_o = AbsCosTheta(wo);
//   }

//   const float3 lambertDir   = lambertSample(float2(rands_1.z, rands_1.w), v, n);
//   const float  lambertPdf   = lambertEvalPDF(lambertDir, v, n);
//   const float  lambertVal   = lambertEvalBSDF(lambertDir, v, n);

//   float t_i = lerp_gather(transmittance, cosTheta_i, MI_ROUGH_TRANSMITTANCE_RES);
  
//   float prob_specular = (1.f - t_i) * spec_weight;
//   float prob_diffuse  = t_i * (1.f - spec_weight);

//   if(prob_diffuse != 0.0f && prob_specular != 0.0f)
//   {
//     prob_specular = prob_specular / (prob_specular + prob_diffuse);
//     prob_diffuse  = 1.f - prob_specular;
//   }
//   else
//   {
//     prob_diffuse  = 1.0f;
//     prob_specular = 0.0f;
//   }

//   if(rands_2 < prob_specular) // specular
//   {
//     pRes->dir = ggxDir;
//     // pRes->pdf = ggxPdf * prob_specular;
//     // pRes->val = float4(ggxVal);
//     // pRes->flags = (alpha == 0.0f) ? RAY_EVENT_S : RAY_FLAG_HAS_NON_SPEC;
//   } 
//   else // lambert
//   {
//     pRes->dir = lambertDir;
//     // pRes->pdf = lambertPdf * (1.0f - prob_specular);
//     // float t_o = lerp_gather(transmittance, cosTheta_o, MI_ROUGH_TRANSMITTANCE_RES);
//     // float4 diff = lambertVal * a_reflSpec;
//     // diff /= 1.f - (nonlinear > 0 ? (diff * internal_refl) : float4(internal_refl));
//     // const float inv_eta_2 = 1.f / (eta * eta);

//     // pRes->val = diff * (inv_eta_2 * cosTheta_o * t_i * t_o);
//     // pRes->flags = RAY_FLAG_HAS_NON_SPEC;
//   }
//   pRes->flags = (alpha == 0.0f) ? RAY_EVENT_S : RAY_FLAG_HAS_NON_SPEC;

//   pRes->pdf   = ggxPdf * prob_specular + lambertPdf * (1.0f - prob_specular);

//   float t_o = lerp_gather(transmittance, cosTheta_o, MI_ROUGH_TRANSMITTANCE_RES);
//   float4 diff = lambertVal * a_reflSpec;
//   diff /= 1.f - (nonlinear > 0 ? (diff * internal_refl) : float4(internal_refl));
//   const float inv_eta_2 = 1.f / (eta * eta);

//   pRes->val = float4(ggxVal) + diff * (inv_eta_2 * cosTheta_o * t_i * t_o);
// }


static inline void plasticSampleAndEval(const Material* a_materials, float4 a_reflSpec, float4 rands,
                                        float3 v, float3 n, float2 tc, BsdfSample* pRes, const float* transmittance)
{
  const uint  cflags    = as_uint(a_materials[0].data[UINT_CFLAGS]);
  const float alpha     = a_materials[0].data[PLASTIC_ROUGHNESS];
  const float eta       = a_materials[0].data[PLASTIC_IOR_RATIO];
  const float spec_weight = a_materials[0].data[PLASTIC_SPEC_SAMPLE_WEIGHT];
  const uint  nonlinear   = as_uint(a_materials[0].data[PLASTIC_NONLINEAR]);
  const float internal_refl = a_materials[0].data[PLASTIC_PRECOMP_REFLECTANCE];
  const float2 alpha2 = {alpha, alpha};

  // n = {0.00234923372, 0.241463557, 0.970407009};
  // v = {0.00075211667, -0.242878497, -0.970056355};

  // float3 nx, ny, nz = n;
  // CoordinateSystem(nz, &nx, &ny);

  float3 s, t;
  CoordinateSystemV2(n, &s, &t);
  t = normalize(cross(n, s));
  
  float3 wi = float3(dot(v, s), dot(v, t), dot(v, n));
  if(wi.z == 0)
    return;

  // wi = float3{-0.000852187979, 0.00306384265, 0.999994874};

  float cos_theta_i = wi.z;
  float t_i = lerp_gather(transmittance, cos_theta_i, MI_ROUGH_TRANSMITTANCE_RES);

  float prob_specular = (1.f - t_i) * spec_weight;
  float prob_diffuse  = t_i * (1.f - spec_weight);

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

  // rands_2 = 0.286642551;

  // rands_1 = {0.385110825, 0.979427397, 0.880511522, 0.258460701};

  const bool sample_specular = rands.z < prob_specular;
  const bool sample_diffuse  = !sample_specular;

  float3 wo {0.0f, 0.0f, 0.0f};
  if(sample_specular)
  {
    const float3 wm = to_float3(sample_visible_normal(wi, {rands.x, rands.y}, alpha2));
    wo = reflect((-1.0f) * wi, wm);
  }

  if(sample_diffuse)
  {
    wo = square_to_cosine_hemisphere({rands.x, rands.y});
  }

  float cos_theta_o = wo.z;
  if(wi.z * wo.z < 0)
  {
    return;
  }

  float3 H = normalize(wo + wi);
  float  D = eval_microfacet(H, alpha2);

  float pdf = D * smith_g1(wi, H, alpha2) / (4.f * cos_theta_i);
  pdf *= prob_specular;
  pdf += prob_diffuse * (wo.z / M_PI);

  const float F = FrDielectric(std::abs(dot(wi, H)), eta); 
  float G = microfacet_G(wi, wo, H, alpha2);
  float val = F * D * G / (4.f * cos_theta_i);

  float t_o = lerp_gather(transmittance, cos_theta_o, MI_ROUGH_TRANSMITTANCE_RES); 
  float4 diffuse = a_reflSpec;
  diffuse /= 1.f - (nonlinear > 0 ? (diffuse * internal_refl) : float4(internal_refl));
  const float inv_eta_2 = 1.f / (eta * eta);

  pRes->dir   = normalize(wo.x * s + wo.y * t + wo.z * n);
  pRes->val   = float4(val) + diffuse * (inv_eta_2 * cos_theta_o * t_i * t_o / M_PI);
  pRes->pdf   = pdf;
  pRes->flags = RAY_FLAG_HAS_NON_SPEC;
}


static void plasticEval(const Material* a_materials, float4 a_reflSpec, float3 l, float3 v, float3 n, float2 tc, 
                        BsdfEval* pRes)
{
  const uint  cflags    = as_uint(a_materials[0].data[UINT_CFLAGS]);
  const float alpha     = a_materials[0].data[PLASTIC_ROUGHNESS];
  const float eta       = a_materials[0].data[PLASTIC_IOR_RATIO];
  const uint  precomp_id  = as_uint(a_materials[0].data[PLASTIC_PRECOMP_ID]);
  const float spec_weight = a_materials[0].data[PLASTIC_SPEC_SAMPLE_WEIGHT];
  const uint  nonlinear   = as_uint(a_materials[0].data[PLASTIC_NONLINEAR]);

  const float internal_refl = a_materials[0].data[PLASTIC_PRECOMP_REFLECTANCE];

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
  // const float coeff   = (1.f - f_i) * (1.f - f_o) / (eta * eta * (1.f - fdr_int));
  // lambertVal         *= coeff;
  // coeffLambertPdf     = coeff; 
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