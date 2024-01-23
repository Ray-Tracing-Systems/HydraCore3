#pragma once
#include "cglobals.h"
#include "crandom.h"
#include "cmaterial.h"
#include "../spectrum.h"
#include <iostream>

struct ThinFilmPrecomputed
{
  std::array<float, FILM_ANGLE_RES * FILM_LENGTH_RES> reflectivity;
  std::array<float, FILM_ANGLE_RES * FILM_LENGTH_RES> transmittance;
};

static inline void filmSmoothSampleAndEval(const Material* a_materials, const float4* eta, const float4* k, const float* thickness,
        uint layers, const float4 a_wavelengths, float4 rands, float3 v, float3 n, float2 tc, BsdfSample* pRes, const float* reflectance)
{
  //std::cout << eta_1[0] << k_1[0] << eta_2[0] << k_2[0] << std::endl;
  const uint cflags = as_uint(a_materials[0].data[UINT_CFLAGS]);

  const float3 pefReflDir = reflect((-1.0f)*v, n);
  const float cosThetaOut = dot(pefReflDir, n);
  float3 dir              = pefReflDir;
  float  pdf              = 1.0f;
  
  float4 val;
  const uint spectralSamples = sizeof(a_wavelengths.M)/sizeof(a_wavelengths.M[0]); 

  uint precompFlag = as_uint(a_materials[0].data[FILM_PRECOMP_FLAG]);
  if (!precompFlag)
  {
    for(int i = 0; i < spectralSamples; ++i)
    {
      if (layers == 1)
      {
        val[i] = FrFilmRefl(cosThetaOut, complex(1.0f), complex(eta[0][i], k[0][i]), complex(eta[1][i], k[1][i]), thickness[0], a_wavelengths[i]); 
      }
      else if (layers > 1)
      {
        val[i] = multFrFilmRefl4(cosThetaOut, eta, k, thickness, layers, a_wavelengths[i], i);
      }
      // BSDF is multiplied (outside) by cosThetaOut. For mirrors this shouldn't be done, so we pre-divide here instead
      val[i] = (cosThetaOut <= 1e-6f) ? 0.0f : (val[i] / std::max(cosThetaOut, 1e-6f));  
    }
  }
  else
  {
    for(int i = 0; i < spectralSamples; ++i)
    {
      float angleVal = acosf(cosThetaOut) / M_PI_2;
      float w = (a_wavelengths[i] - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN);
      val[i] = lerp_gather_2d(reflectance, w, angleVal, FILM_LENGTH_RES, FILM_ANGLE_RES);
      // BSDF is multiplied (outside) by cosThetaOut. For mirrors this shouldn't be done, so we pre-divide here instead
      val[i] = (cosThetaOut <= 1e-6f) ? 0.0f : (val[i] / std::max(cosThetaOut, 1e-6f));  
    }
  }
  
  pRes->val = val; 
  pRes->dir = dir;
  pRes->pdf = pdf;
  pRes->flags = RAY_EVENT_S;
}


static void filmSmoothEval(const Material* a_materials, const float4 eta_1, const float4 k_1, const float4 eta_2, const float4 k_2, float4 wavelengths, float3 l, float3 v, float3 n, float2 tc,
                                BsdfEval* pRes)
{
  pRes->val = {0.0f, 0.0f, 0.0f, 0.0f};
  pRes->pdf = 0.0f;
}

static float filmRoughEvalInternalPrecomp(float3 wo, float3 wi, float3 wm, float2 alpha, float lambda, const float* reflectance)
{
  if(wo.z * wi.z < 0) // not in the same hemisphere
  {
    return 0.0f;
  }

  float cosTheta_o = AbsCosTheta(wo);
  float cosTheta_i = AbsCosTheta(wi);
  if (cosTheta_i == 0 || cosTheta_o == 0)
    return 0.0f;
  
  float cosTheta = std::abs(dot(wo, wm));
  float w = (lambda - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN);
  float F = lerp_gather_2d(reflectance, w, cosTheta, FILM_LENGTH_RES, FILM_ANGLE_RES);
  float val = trD(wm, alpha) * F * trG(wo, wi, alpha) / (4.0f * cosTheta_i * cosTheta_o);

  return val;
}

static float filmRoughEvalInternal(float3 wo, float3 wi, float3 wm, float2 alpha, complex ior, complex ior2, float thickness, float lambda)
{
  if(wo.z * wi.z < 0) // not in the same hemisphere
  {
    return 0.0f;
  }

  float cosTheta_o = AbsCosTheta(wo);
  float cosTheta_i = AbsCosTheta(wi);
  if (cosTheta_i == 0 || cosTheta_o == 0)
    return 0.0f;

  float F = FrFilmRefl(std::abs(dot(wo, wm)), complex(1.0f), ior, ior2, thickness, lambda);
  float val = trD(wm, alpha) * F * trG(wo, wi, alpha) / (4.0f * cosTheta_i * cosTheta_o);

  return val;
}

static float filmRoughEvalInternal2(float3 wo, float3 wi, float3 wm, float2 alpha, const float4* eta, const float4* k, const float* thickness, uint layers, float lambda, uint comp)
{
  if(wo.z * wi.z < 0) // not in the same hemisphere
  {
    return 0.0f;
  }

  float cosTheta_o = AbsCosTheta(wo);
  float cosTheta_i = AbsCosTheta(wi);
  if (cosTheta_i == 0 || cosTheta_o == 0)
    return 0.0f;

  float F = multFrFilmRefl4(std::abs(dot(wo, wm)), eta, k, thickness, layers, lambda, comp);
  float val = trD(wm, alpha) * F * trG(wo, wi, alpha) / (4.0f * cosTheta_i * cosTheta_o);

  return val;
}


static inline void filmRoughSampleAndEval(const Material* a_materials, const float4* eta, const float4* k, const float* thickness,
        uint layers, const float4 a_wavelengths, float4 rands, float3 v, float3 n, float2 tc, float3 alpha_tex, BsdfSample* pRes, const float* reflectance)
{
  //std::cout << a_wavelengths[0] << std::endl;
  if(v.z == 0)
    return;

  const uint cflags = as_uint(a_materials[0].data[UINT_CFLAGS]);
  // const uint  texId = as_uint(a_materials[0].data[FILM_TEXID0]);

  const float2 alpha = float2(min(a_materials[0].data[FILM_ROUGH_V], alpha_tex.x), 
                              min(a_materials[0].data[FILM_ROUGH_U], alpha_tex.y));
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

  float4 val;
  const uint spectralSamples = sizeof(a_wavelengths.M)/sizeof(a_wavelengths.M[0]); 

  uint precompFlag = as_uint(a_materials[0].data[FILM_PRECOMP_FLAG]);
  if (!precompFlag)
  {
    for(int i = 0; i < spectralSamples; ++i)
    {
      if (layers == 1)
      {
        val[i] = filmRoughEvalInternal(wo, wi, wm, alpha, complex(eta[0][i], k[0][i]), complex(eta[1][i], k[1][i]), thickness[0], a_wavelengths[i]);
      }
      else if (layers > 1)
      {
        val[i] = filmRoughEvalInternal2(wo, wi, wm, alpha, eta, k, thickness, layers, a_wavelengths[i], i);
      }
    }
  }
  else
  {
    for(int i = 0; i < spectralSamples; ++i)
    {
      val[i] = filmRoughEvalInternalPrecomp(wo, wi, wm, alpha, a_wavelengths[i], reflectance);
    }
  }
  

  pRes->val   = val; 
  pRes->dir   = normalize(wi.x * nx + wi.y * ny + wi.z * nz);
  pRes->pdf   = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
  pRes->flags = RAY_FLAG_HAS_NON_SPEC;
}


static void filmRoughEval(const Material* a_materials, const float4* eta, const float4* k, const float* thickness,
        uint layers, const float4 a_wavelengths, float3 l, float3 v, float3 n, float2 tc, float3 alpha_tex, BsdfEval* pRes, const float* reflectance)
{
  const uint cflags = as_uint(a_materials[0].data[UINT_CFLAGS]);

  // const float2 alpha = float2(a_materials[0].data[FILM_ROUGH_V], a_materials[0].data[FILM_ROUGH_U]);
  const float2 alpha = float2(min(a_materials[0].data[FILM_ROUGH_V], alpha_tex.x), 
                              min(a_materials[0].data[FILM_ROUGH_U], alpha_tex.y));

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
  float4 val;
  const uint spectralSamples = sizeof(a_wavelengths.M)/sizeof(a_wavelengths.M[0]); 

  uint precompFlag = as_uint(a_materials[0].data[FILM_PRECOMP_FLAG]);
  if (!precompFlag)
  {
    for(int i = 0; i < spectralSamples; ++i)
    {
      if (layers == 1)
      {
        val[i] = filmRoughEvalInternal(wo, wi, wm, alpha, complex(eta[0][i], k[0][i]), complex(eta[1][i], k[1][i]), thickness[0], a_wavelengths[i]);
      }
      else if (layers > 1)
      {
        val[i] = filmRoughEvalInternal2(wo, wi, wm, alpha, eta, k, thickness, layers, a_wavelengths[i], i);
      }
    }
  }
  else
  {
    for(int i = 0; i < spectralSamples; ++i)
    {
      val[i] = filmRoughEvalInternalPrecomp(wo, wi, wm, alpha, a_wavelengths[i], reflectance);
    }
  }

  pRes->val = val;
  wm        = FaceForward(wm, float3(0.0f, 0.0f, 1.0f));
  pRes->pdf = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
}