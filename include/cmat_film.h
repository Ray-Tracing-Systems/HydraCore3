#pragma once
#include "cglobals.h"
#include "crandom.h"
#include "cmaterial.h"
#include "../spectrum.h"
#include <iostream>

struct ThinFilmPrecomputed
{
  std::array<float, FILM_ANGLE_RES * FILM_LENGTH_RES * 2> reflectivity;
  std::array<float, FILM_ANGLE_RES * FILM_LENGTH_RES * 2> transmittivity;
};

static inline void filmSmoothSampleAndEval(const Material* a_materials, const float4* eta, const float4* k, const float* thickness,
        uint layers, const float4 a_wavelengths, const float _extIOR, float4 rands, float3 v, float3 n, float2 tc, BsdfSample* pRes, const float* reflectance)
{
  const float extIOR = 1.f;
  bool reversed = false;
  bool opaque = false;
  if ((pRes->flags & RAY_FLAG_HAS_INV_NORMAL) != 0) 
  {
    n = -1 * n;
    reversed = true;
    reflectance += FILM_ANGLE_RES * FILM_LENGTH_RES;
  }
  else
  {
    if (length(k[layers - 1]) > 1e-3f)
    {
      opaque = true;
    }
  }

  float3 s, t = n;
  CoordinateSystemV2(n, &s, &t);
  float3 wi = float3(dot(v, s), dot(v, t), dot(v, n));
  float cosThetaI = clamp(fabs(wi.z), 0.0001, 1.0f);

  float ior = eta[layers - 1].x;

  float4 fr = FrDielectricDetailedV2(wi.z, ior);

  const float cosThetaT = fr.y;
  const float eta_it = fr.z;
  const float eta_ti = fr.w;  
  
  float4 R = {0.f, 0.f, 0.f, 0.f};
  const uint spectralSamples = opaque ? sizeof(a_wavelengths.M)/sizeof(a_wavelengths.M[0]) : 1; 

  uint precompFlag = as_uint(a_materials[0].data[FILM_PRECOMP_FLAG]);
  if (!precompFlag)
  {
    for(int i = 0; i < spectralSamples; ++i)
    {
      if (layers == 2)
      {
        if (!reversed)
          R[i] = FrFilmRefl(cosThetaI, complex(1.0f), complex(eta[0][i], k[0][i]), complex(eta[1][i], k[1][i]), thickness[0], a_wavelengths[i]); 
        else
          R[i] = FrFilmRefl(cosThetaI, complex(eta[1][i], k[1][i]), complex(eta[0][i], k[0][i]), complex(1.0f), thickness[0], a_wavelengths[i]); 
      }
      else if (layers > 2)
      {
        if (!reversed)
          R[i] = multFrFilmRefl4(cosThetaI, eta, k, thickness, layers, a_wavelengths[i], i);
        else
          R[i] = multFrFilmRefl4_r(cosThetaI, eta, k, thickness, layers, a_wavelengths[i], i);
      } 
    }
  }
  else
  {
    for(int i = 0; i < spectralSamples; ++i)
    {
      float angleVal = acosf(cosThetaI) / M_PI_2;
      float w = (a_wavelengths[i] - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN);
      R[i] = lerp_gather_2d(reflectance, w, angleVal, FILM_LENGTH_RES, FILM_ANGLE_RES);
    }
  }
  if (opaque)
  {
    float3 wo = float3(-wi.x, -wi.y, wi.z);
    pRes->val = R;
    pRes->pdf = 1.f;
    pRes->dir = normalize(wo.x * s + wo.y * t + wo.z * n);
    pRes->flags |= RAY_EVENT_S;
    pRes->ior = _extIOR;
  }
  else
  {
    float T;
    if (layers == 2)
    {
      if (!reversed)
        T = FrFilmRefr(cosThetaI, complex(1.0f), complex(eta[0][0], k[0][0]), complex(eta[1][0], k[1][0]), thickness[0], a_wavelengths[0]); 
      else
        T = FrFilmRefr(cosThetaI, complex(eta[1][0], k[1][0]), complex(eta[0][0], k[0][0]), complex(1.0f), thickness[0], a_wavelengths[0]); 
    }
    T *= eta_it * clamp(fabs(cosThetaT), 0.0f, 1.0f) / cosThetaI;
    if (rands.x / (R.x + T) < R.x)
    {
      float3 wo = float3(-wi.x, -wi.y, wi.z);
      pRes->val = R;
      pRes->pdf = R.x / (R.x + T);
      pRes->dir = normalize(wo.x * s + wo.y * t + wo.z * n);
      pRes->flags |= RAY_EVENT_S;
      pRes->ior = _extIOR;
    }
    else
    {
      float3 wo = refract(wi, cosThetaT, eta_ti);
      pRes->val = float4(T);
      pRes->pdf = T / (R.x + T);
      pRes->dir = normalize(wo.x * s + wo.y * t + wo.z * n);
      pRes->flags |= (RAY_EVENT_S | RAY_EVENT_T);
      pRes->ior = (_extIOR == ior) ? extIOR : ior;
    }
  }

  if (!opaque)
  {
    pRes->flags |= RAY_FLAG_WAVES_DIVERGED;
  }

  pRes->val /= std::max(std::abs(dot(pRes->dir, n)), 1e-6f);

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
  
  float angle = acosf(std::abs(dot(wo, wm))) / M_PI_2;
  float w = (lambda - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN);
  float F = lerp_gather_2d(reflectance, w, angle, FILM_LENGTH_RES, FILM_ANGLE_RES);
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
      if (layers == 2)
      {
        val[i] = filmRoughEvalInternal(wo, wi, wm, alpha, complex(eta[0][i], k[0][i]), complex(eta[1][i], k[1][i]), thickness[0], a_wavelengths[i]);
      }
      else if (layers > 2)
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
      if (layers == 2)
      {
        val[i] = filmRoughEvalInternal(wo, wi, wm, alpha, complex(eta[0][i], k[0][i]), complex(eta[1][i], k[1][i]), thickness[0], a_wavelengths[i]);
      }
      else if (layers > 2)
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