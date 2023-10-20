#pragma once
#include "cglobals.h"
#include "crandom.h"
#include "cmaterial.h"
#include "../spectrum.h"


static inline float3 conductorSampleEtaSpectrum(const Material* a_materials, const Spectrum* a_spectra, const float3 a_wavelengths)
{
  float3 eta = float3(a_materials[0].data[CONDUCTOR_ETA]);
  if(a_wavelengths.M[0] == 0.0f)
    return eta;

  const uint eta_specId = as_uint(a_materials[0].data[CONDUCTOR_ETA_SPECID]);

  if(eta_specId < 0xFFFFFFFF)
  {
    eta = SampleSpectrum(a_spectra + eta_specId, a_wavelengths);
    // const uint spectralSamples = uint(sizeof(a_wavelengths.M)/sizeof(a_wavelengths.M[0])); 
    // for(uint i = 0; i < spectralSamples; ++i)
    //   eta[i] = a_spectra[eta_specId].Sample(a_wavelengths[i]);
  }
  
  return eta;
}

static inline float3 conductorSampleKSpectrum(const Material* a_materials, const Spectrum* a_spectra, const float3 a_wavelengths)
{
  float3 k = float3(a_materials[0].data[CONDUCTOR_K]);
  if(a_wavelengths.M[0] == 0.0f)
    return k;

  const uint k_specId = as_uint(a_materials[0].data[CONDUCTOR_K_SPECID]);

  if(k_specId < 0xFFFFFFFF)
  {
    k = SampleSpectrum(a_spectra + k_specId, a_wavelengths);
    // const uint spectralSamples = uint(sizeof(a_wavelengths.M)/sizeof(a_wavelengths.M[0])); 
    // for(uint i = 0; i < spectralSamples; ++i)
    //   k[i] = a_spectra[k_specId].Sample(a_wavelengths[i]);
  }
  
  return k;
}

static inline void conductorSmoothSampleAndEval(const Material* a_materials, const Spectrum* a_spectra, const float3 a_wavelengths,
                                                float4 rands, float3 v, float3 n, float2 tc,
                                                BsdfSample* pRes)
{
  const uint cflags = as_uint(a_materials[0].data[UINT_CFLAGS]);

  float3 eta = conductorSampleEtaSpectrum(a_materials, a_spectra, a_wavelengths);
  float3 k   = conductorSampleKSpectrum(a_materials, a_spectra, a_wavelengths);

  const float3 pefReflDir = reflect((-1.0f)*v, n);
  const float cosThetaOut = dot(pefReflDir, n);
  float3 dir              = pefReflDir;
  float  pdf              = 1.0f;
  
  float3 val;
  const uint spectralSamples = sizeof(a_wavelengths.M)/sizeof(a_wavelengths.M[0]); 
  for(int i = 0; i < spectralSamples; ++i)
  {
    val[i] = FrComplexConductor(cosThetaOut, complex{eta[i], k[i]});
    // BSDF is multiplied (outside) by cosThetaOut. For mirrors this shouldn't be done, so we pre-divide here instead
    val[i] = (cosThetaOut <= 1e-6f) ? 0.0f : (val[i] / std::max(cosThetaOut, 1e-6f));  
  }
  
  pRes->val = val; 
  pRes->dir = dir;
  pRes->pdf = pdf;
  pRes->flags = RAY_EVENT_S;
}


static void conductorSmoothEval(const Material* a_materials, float3 wavelengths, float3 l, float3 v, float3 n, float2 tc,
                                BsdfEval* pRes)
{
  pRes->val = {0.0f, 0.0f, 0.0f};
  pRes->pdf = 0.0f;
}


static float conductorRoughEvalInternal(float3 wo, float3 wi, float3 wm, float2 alpha, complex ior)
{
  if(wo.z * wi.z < 0) // not in the same hemisphere
  {
    return 0.0f;
  }

  float cosTheta_o = AbsCosTheta(wo);
  float cosTheta_i = AbsCosTheta(wi);
  if (cosTheta_i == 0 || cosTheta_o == 0)
    return 0.0f;

  float F = FrComplexConductor(std::abs(dot(wo, wm)), ior);
  float val = trD(wm, alpha) * F * trG(wo, wi, alpha) / (4.0f * cosTheta_i * cosTheta_o);

  return val;
}


static inline void conductorRoughSampleAndEval(const Material* a_materials, const Spectrum* a_spectra, const float3 a_wavelengths, 
                                               float4 rands, float3 v, float3 n, float2 tc, float3 alpha_tex, 
                                               BsdfSample* pRes)
{
  if(v.z == 0)
    return;

  const uint cflags = as_uint(a_materials[0].data[UINT_CFLAGS]);

  float3 eta = conductorSampleEtaSpectrum(a_materials, a_spectra, a_wavelengths);
  float3 k   = conductorSampleKSpectrum(a_materials, a_spectra, a_wavelengths);
  // const uint  texId = as_uint(a_materials[0].data[CONDUCTOR_TEXID0]);

  const float2 alpha = float2(min(a_materials[0].data[CONDUCTOR_ROUGH_V], alpha_tex.x), 
                              min(a_materials[0].data[CONDUCTOR_ROUGH_U], alpha_tex.y));
  // if(texId != 0)
  // {
  //   alpha.x = alpha_tex.x;
  //   alpha.y = alpha_tex.y;
  // }

  float3 nx, ny, nz = n;
  CoordinateSystem(nz, &nx, &ny);
  const float3 wo = float3(dot(v, nx), dot(v, ny), dot(v, nz));
  if(wo.z == 0)
    return;

  if(wo.z == 0)
    return;

  float3 wm = trSample(wo, float2(rands.x, rands.y), alpha);
  float3 wi = reflect((-1.0f) * wo, wm);

  if(wo.z * wi.z < 0) // not in the same hemisphere
  {
    return;
  }

  float3 val;
  const uint spectralSamples = sizeof(a_wavelengths.M)/sizeof(a_wavelengths.M[0]); 
  for(int i = 0; i < spectralSamples; ++i)
  {
    val[i] = conductorRoughEvalInternal(wo, wi, wm, alpha, complex{eta[i], k[i]});
  }

  pRes->val   = val; 
  pRes->dir   = normalize(wi.x * nx + wi.y * ny + wi.z * nz);
  pRes->pdf   = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
  pRes->flags = RAY_FLAG_HAS_NON_SPEC;
}


static void conductorRoughEval(const Material* a_materials, const Spectrum* a_spectra, const float3 a_wavelengths,
                               float3 l, float3 v, float3 n, float2 tc, float3 alpha_tex, 
                               BsdfEval* pRes)
{
  const uint cflags = as_uint(a_materials[0].data[UINT_CFLAGS]);

  float3 eta = conductorSampleEtaSpectrum(a_materials, a_spectra, a_wavelengths);
  float3 k   = conductorSampleKSpectrum(a_materials, a_spectra, a_wavelengths);

  // const float2 alpha = float2(a_materials[0].data[CONDUCTOR_ROUGH_V], a_materials[0].data[CONDUCTOR_ROUGH_U]);
  const float2 alpha = float2(min(a_materials[0].data[CONDUCTOR_ROUGH_V], alpha_tex.x), 
                              min(a_materials[0].data[CONDUCTOR_ROUGH_U], alpha_tex.y));

  float3 nx, ny, nz = n;
  CoordinateSystem(nz, &nx, &ny);

  // v = (-1.0f) * v;
  const float3 wo = float3(dot(v, nx), dot(v, ny), dot(v, nz));
  const float3 wi = float3(dot(l, nx), dot(l, ny), dot(l, nz));

  if(wo.z * wi.z < 0.0f)
    return;

  float3 wm = wo + wi;
  if (dot(wm, wm) == 0)
      return;

  wm = normalize(wm);
  float3 val;
  const uint spectralSamples = sizeof(a_wavelengths.M)/sizeof(a_wavelengths.M[0]); 
  for(int i = 0; i < spectralSamples; ++i)
  {
    val[i] = conductorRoughEvalInternal(wo, wi, wm, alpha, complex{eta[i], k[i]});
  }

  pRes->val = val;

  wm        = FaceForward(wm, float3(0.0f, 0.0f, 1.0f));
  pRes->pdf = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
}

// static inline void conductorSampleAndEval(const Material* a_materials, float3 wavelengths, float4 rands, float3 v, float3 n, float2 tc,
//                                           float3 alpha_tex,
//                                           BsdfSample* pRes)
// {
//   const float2 alpha = float2(a_materials[0].data[CONDUCTOR_ROUGH_V], a_materials[0].data[CONDUCTOR_ROUGH_U]);

//   if(trEffectivelySmooth(alpha))
//   {
//     conductorSmoothSampleAndEval(a_materials, a_spectra, wavelengths, rands, v, n, tc, pRes);
//   }
//   else
//   {
//     conductorRoughSampleAndEval(a_materials, a_spectra, wavelengths, rands, v, n, tc, alpha_tex, pRes);
//   }
// }

// static inline void conductorEval(const Material* a_materials, float3 wavelengths, float3 l, float3 v, float3 n, float2 tc,
//                                 float3 alpha_tex,
//                                 BsdfEval* pRes)
// {
//   const float2 alpha = float2(a_materials[0].data[CONDUCTOR_ROUGH_V], a_materials[0].data[CONDUCTOR_ROUGH_U]);

//   if(trEffectivelySmooth(alpha))
//   {
//     conductorSmoothEval(a_materials, wavelengths, l, v, n, tc, pRes);
//   }
//   else
//   {
//     conductorRoughEval(a_materials, wavelengths, l, v, n, tc, alpha_tex, pRes);
//   }
// }