#pragma once
#include "cglobals.h"
#include "crandom.h"
#include "cmaterial.h"
#include "../spectrum.h"
#include <iostream>

static inline void filmSmoothSampleAndEval(const Material* a_materials, const complex* a_ior, const float* thickness,
        uint layers, const float4 a_wavelengths, const float _extIOR, float4 rands, float3 v, float3 n, float2 tc, BsdfSample* pRes,
        const float* precomputed_data)
{
  const float extIOR = a_materials[0].data[FILM_ETA_EXT];

  bool reversed = false;
  uint32_t refl_offset;
  uint32_t refr_offset;
  if ((pRes->flags & RAY_FLAG_HAS_INV_NORMAL) != 0) // inside of object
  {
    n = -1 * n;
    reversed = true;
    refl_offset = FILM_ANGLE_RES * FILM_LENGTH_RES * 2;
    refr_offset = FILM_ANGLE_RES * FILM_LENGTH_RES * 3;
  }
  else // outside of object
  {
    refl_offset = 0;
    refr_offset = FILM_ANGLE_RES * FILM_LENGTH_RES;
  }

  float3 s, t = n;
  CoordinateSystemV2(n, &s, &t);
  float3 wi = float3(dot(v, s), dot(v, t), dot(v, n));

  float cosThetaI = clamp(fabs(wi.z), 0.0001, 1.0f);

  float ior = a_ior[layers].re / extIOR;

  float4 fr = FrDielectricDetailedV2(wi.z, ior);

  const float cosThetaT = fr.y;
  const float eta_it = fr.z;
  const float eta_ti = fr.w;  
  
  float R, T;
  //FrReflRefr result;
  /*
  uint precompFlag = as_uint(a_materials[0].data[FILM_PRECOMP_FLAG]);
  if (!precompFlag)
  {
    if (layers == 2)
    {
      if (!reversed)
      {
        result.refl = FrFilmRefl(cosThetaI, a_ior[0], a_ior[1], a_ior[2], thickness[0], a_wavelengths[0]); 
        result.refr = FrFilmRefr(cosThetaI, a_ior[0], a_ior[1], a_ior[2], thickness[0], a_wavelengths[0]); 
      }
      else
      {
        result.refl = FrFilmRefl(cosThetaI, a_ior[2], a_ior[1], a_ior[0], thickness[0], a_wavelengths[0]); 
        result.refr = FrFilmRefr(cosThetaI, a_ior[2], a_ior[1], a_ior[0], thickness[0], a_wavelengths[0]); 
      }
    }
    else if (layers > 2)
    {
      if (!reversed)
        result = multFrFilm(cosThetaI, a_ior, thickness, layers, a_wavelengths[0], a_cosTheta, a_phaseDiff);
      else
        result = multFrFilm_r(cosThetaI, a_ior, thickness, layers, a_wavelengths[0], a_cosTheta, a_phaseDiff);
    } 
  }
  else
  */
  {
    float w = (a_wavelengths[0] - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN);
    float theta = acos(cosThetaI) * 2.f / M_PI;
    //result.refl = lerp_gather_2d(reflectance, w, theta, FILM_LENGTH_RES, FILM_ANGLE_RES);
    //result.refr = lerp_gather_2d(transmittance, w, theta, FILM_LENGTH_RES, FILM_ANGLE_RES);
    w *= FILM_LENGTH_RES - 1;
    theta *= FILM_ANGLE_RES - 1;
    uint32_t index1 = std::min(uint32_t(w), uint32_t(FILM_LENGTH_RES - 2));
    uint32_t index2 = std::min(uint32_t(theta), uint32_t(FILM_LENGTH_RES - 2));

    float alpha = w - float(index1);
    float beta = theta - float(index2);

    float v0 = lerp(precomputed_data[refl_offset + index1 * FILM_LENGTH_RES + index2], precomputed_data[refl_offset + (index1 + 1) * FILM_LENGTH_RES + index2], alpha);
    float v1 = lerp(precomputed_data[refl_offset + index1 * FILM_LENGTH_RES + index2 + 1], precomputed_data[refl_offset + (index1 + 1) * FILM_LENGTH_RES + index2 + 1], alpha);
    R = lerp(v0, v1, beta);

    v0 = lerp(precomputed_data[refr_offset + index1 * FILM_LENGTH_RES + index2], precomputed_data[refr_offset + (index1 + 1) * FILM_LENGTH_RES + index2], alpha);
    v1 = lerp(precomputed_data[refr_offset + index1 * FILM_LENGTH_RES + index2 + 1], precomputed_data[refr_offset + (index1 + 1) * FILM_LENGTH_RES + index2 + 1], alpha);
    T = lerp(v0, v1, beta);
  }
  //R = result.refl;
  //T = result.refr;

  if (a_ior[layers].im > 0.001)
  {
    float3 wo = float3(-wi.x, -wi.y, wi.z);
    pRes->val = float4(R);
    pRes->pdf = 1.f;
    pRes->dir = normalize(wo.x * s + wo.y * t + wo.z * n);
    pRes->flags |= RAY_EVENT_S;
    pRes->ior = _extIOR;
  }
  else
  {
    if (rands.x * (R + T) < R)
    {
      float3 wo = float3(-wi.x, -wi.y, wi.z);
      pRes->val = float4(R);
      pRes->pdf = R / (R + T);
      pRes->dir = normalize(wo.x * s + wo.y * t + wo.z * n);
      pRes->flags |= RAY_EVENT_S;
      pRes->ior = _extIOR;
    }
    else
    {
      float3 wo = refract(wi, cosThetaT, eta_ti);
      pRes->val = float4(T);
      pRes->pdf = T / (R + T);
      pRes->dir = normalize(wo.x * s + wo.y * t + wo.z * n);
      pRes->flags |= (RAY_EVENT_S | RAY_EVENT_T);
      pRes->ior = (_extIOR == a_ior[layers].re) ? extIOR : a_ior[layers].re;
    }
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
  return 0.0f; //TO DO
  /*
  if(wo.z * wi.z < 0) // not in the same hemisphere
  {
    return 0.0f;
  }

  float cosTheta_o = AbsCosTheta(wo);
  float cosTheta_i = AbsCosTheta(wi);
  if (cosTheta_i == 0 || cosTheta_o == 0)
    return 0.0f;
  float w = (lambda - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN);
  float theta = acos(std::abs(dot(wo, wm))) * 2.f / M_PI;

  //float F = lerp_gather_2d(reflectance, w, theta, FILM_LENGTH_RES, FILM_ANGLE_RES);
  w *= FILM_LENGTH_RES - 1;
  theta *= FILM_ANGLE_RES - 1;
  uint32_t index1 = std::min(uint32_t(w), uint32_t(FILM_LENGTH_RES - 2));
  uint32_t index2 = std::min(uint32_t(theta), uint32_t(FILM_LENGTH_RES - 2));

  float a = w - float(index1);
  float b = theta - float(index2);

  float v0 = lerp(reflectance[index1 * FILM_LENGTH_RES + index2], reflectance[(index1 + 1) * FILM_LENGTH_RES + index2], a);
  float v1 = lerp(reflectance[index1 * FILM_LENGTH_RES + index2 + 1], reflectance[(index1 + 1) * FILM_LENGTH_RES + index2 + 1], a);
  float R = lerp(v0, v1, b);

  //v0 = lerp(transmittance[index1 * FILM_LENGTH_RES + index2], transmittance[(index1 + 1) * FILM_LENGTH_RES + index2], alpha);
  //v1 = lerp(transmittance[index1 * FILM_LENGTH_RES + index2 + 1], transmittance[(index1 + 1) * FILM_LENGTH_RES + index2 + 1], alpha);
  //T = lerp(v0, v1, beta);


  float val = trD(wm, alpha) * R * trG(wo, wi, alpha) / (4.0f * cosTheta_i * cosTheta_o);

  return val;
  */
}

static float filmRoughEvalInternal(float3 wo, float3 wi, float3 wm, float2 alpha, complex ior1, complex ior2, complex ior3, float thickness, float lambda)
{
  if(wo.z * wi.z < 0) // not in the same hemisphere
  {
    return 0.0f;
  }

  float cosTheta_o = AbsCosTheta(wo);
  float cosTheta_i = AbsCosTheta(wi);
  if (cosTheta_i == 0 || cosTheta_o == 0)
    return 0.0f;

  float F = FrFilmRefl(std::abs(dot(wo, wm)), ior1, ior2, ior3, thickness, lambda);
  float val = trD(wm, alpha) * F * trG(wo, wi, alpha) / (4.0f * cosTheta_i * cosTheta_o);

  return val;
}

/*
static float filmRoughEvalInternal2(float3 wo, float3 wi, float3 wm, float2 alpha, const complex* a_ior, const float* thickness, uint layers, float lambda, uint comp)
{
  if(wo.z * wi.z < 0) // not in the same hemisphere
  {
    return 0.0f;
  }

  float cosTheta_o = AbsCosTheta(wo);
  float cosTheta_i = AbsCosTheta(wi);
  if (cosTheta_i == 0 || cosTheta_o == 0)
    return 0.0f;

  float F = multFrFilm(std::abs(dot(wo, wm)), a_ior, thickness, layers, lambda).refl;
  float val = trD(wm, alpha) * F * trG(wo, wi, alpha) / (4.0f * cosTheta_i * cosTheta_o);

  return val;
}
*/

static inline void filmRoughSampleAndEval(const Material* a_materials, const complex* a_ior, const float* thickness,
        uint layers, const float4 a_wavelengths, float4 rands, float3 v, float3 n, float2 tc, float3 alpha_tex, BsdfSample* pRes, const float* reflectance)
{
  //TO DO
  /*
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
*/
  /*
  uint precompFlag = as_uint(a_materials[0].data[FILM_PRECOMP_FLAG]);
  if (!precompFlag)
  {
    for(int i = 0; i < spectralSamples; ++i)
    {
      if (layers == 2)
      {
        val[i] = filmRoughEvalInternal(wo, wi, wm, alpha, a_ior[0], a_ior[1], a_ior[2], thickness[0], a_wavelengths[i]);
      }
      else if (layers > 2)
      {
        val[i] = filmRoughEvalInternal2(wo, wi, wm, alpha, a_ior, thickness, layers, a_wavelengths[i], i, a_cosTheta, a_phaseDiff);
      }
    }
  }
  else
  */
 /*
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
*/
}


static void filmRoughEval(const Material* a_materials, const complex* a_ior, const float* thickness,
        uint layers, const float4 a_wavelengths, float3 l, float3 v, float3 n, float2 tc, float3 alpha_tex, BsdfEval* pRes, const float* reflectance)
{
  // TO DO
  /*
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
*/
/*
  uint precompFlag = as_uint(a_materials[0].data[FILM_PRECOMP_FLAG]);
  if (!precompFlag)
  {
    for(int i = 0; i < spectralSamples; ++i)
    {
      if (layers == 2)
      {
        val[i] = filmRoughEvalInternal(wo, wi, wm, alpha, a_ior[0], a_ior[1], a_ior[2], thickness[0], a_wavelengths[i]);
      }
      else if (layers > 2)
      {
        val[i] = filmRoughEvalInternal2(wo, wi, wm, alpha, a_ior, thickness, layers, a_wavelengths[i], i, a_cosTheta, a_phaseDiff);
      }
    }
  }
  else
  */
 /*
  {
    for(int i = 0; i < spectralSamples; ++i)
    {
      val[i] = filmRoughEvalInternalPrecomp(wo, wi, wm, alpha, a_wavelengths[i], reflectance);
    }
  }

  pRes->val = val;
  wm        = FaceForward(wm, float3(0.0f, 0.0f, 1.0f));
  pRes->pdf = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
  */
}