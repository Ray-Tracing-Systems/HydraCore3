#pragma once
#include "cglobals.h"
#include "crandom.h"
#include "cmaterial.h"
#include "../spectrum.h"
#include <iostream>

static inline void filmSmoothSampleAndEvalPrecomputed(const Material* a_materials, 
        const float extIOR, const complex intIOR, const float4 a_wavelengths, const float _extIOR,
        float4 rands, float3 v, float3 n, float2 tc, BsdfSample* pRes, const float* precomputed_data, const bool spectral_mode)
{
  bool reversed = false;
  uint32_t refl_offset;
  uint32_t refr_offset;
  if ((pRes->flags & RAY_FLAG_HAS_INV_NORMAL) != 0) // inside of object
  {
    n = -1 * n;
  }
  if (dot(n, v) < 0.f && intIOR.im < 0.001)
  {
    reversed = true;
    refl_offset = FILM_ANGLE_RES * 2;
    refr_offset = FILM_ANGLE_RES * 3;
  }
  else
  {
    refl_offset = 0;
    refr_offset = FILM_ANGLE_RES;
  }

  float3 s, t = n;
  CoordinateSystemV2(n, &s, &t);
  float3 wi = float3(dot(v, s), dot(v, t), dot(v, n));

  float cosThetaI = clamp(fabs(wi.z), 0.0001, 1.0f);
  float ior = intIOR.re / extIOR;
  float4 R = float4(0.0f), T = float4(0.0f);

  if (spectral_mode)
  {
    float w = clamp((a_wavelengths[0] - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN), 0.f, 1.f);
    float theta = clamp(acos(cosThetaI) * 2.f / M_PI, 0.f, 1.f);
    //result.refl = lerp_gather_2d(reflectance, w, theta, FILM_LENGTH_RES, FILM_ANGLE_RES);
    //result.refr = lerp_gather_2d(transmittance, w, theta, FILM_LENGTH_RES, FILM_ANGLE_RES);
    w *= FILM_LENGTH_RES - 1;
    theta *= FILM_ANGLE_RES - 1;
    uint32_t index1 = std::min(uint32_t(w), uint32_t(FILM_LENGTH_RES - 2));
    uint32_t index2 = std::min(uint32_t(theta), uint32_t(FILM_ANGLE_RES - 2));

    float alpha = w - float(index1);
    float beta = theta - float(index2);

    float v0 = lerp(precomputed_data[refl_offset * FILM_LENGTH_RES + index1 * FILM_ANGLE_RES + index2], precomputed_data[refl_offset * FILM_LENGTH_RES + (index1 + 1) * FILM_ANGLE_RES + index2], alpha);
    float v1 = lerp(precomputed_data[refl_offset * FILM_LENGTH_RES + index1 * FILM_ANGLE_RES + index2 + 1], precomputed_data[refl_offset * FILM_LENGTH_RES + (index1 + 1) * FILM_ANGLE_RES + index2 + 1], alpha);
    R[0] = lerp(v0, v1, beta);

    v0 = lerp(precomputed_data[refr_offset * FILM_LENGTH_RES + index1 * FILM_ANGLE_RES + index2], precomputed_data[refr_offset * FILM_LENGTH_RES + (index1 + 1) * FILM_ANGLE_RES + index2], alpha);
    v1 = lerp(precomputed_data[refr_offset * FILM_LENGTH_RES + index1 * FILM_ANGLE_RES + index2 + 1], precomputed_data[refr_offset * FILM_LENGTH_RES + (index1 + 1) * FILM_ANGLE_RES + index2 + 1], alpha);
    T[0] = lerp(v0, v1, beta);
  }
  else
  {
    float theta = clamp(acos(cosThetaI) * 2.f / M_PI, 0.f, 1.f);
    theta *= FILM_ANGLE_RES - 1;
    uint32_t index = std::min(uint32_t(theta), uint32_t(FILM_ANGLE_RES - 2));

    float alpha = theta - float(index);

    R[0] = lerp(precomputed_data[(refl_offset + index) * 3], precomputed_data[(refl_offset + index + 1) * 3], alpha);
    R[1] = lerp(precomputed_data[(refl_offset + index) * 3 + 1], precomputed_data[(refl_offset + index + 1) * 3 + 1], alpha);
    R[2] = lerp(precomputed_data[(refl_offset + index) * 3 + 2], precomputed_data[(refl_offset + index + 1) * 3 + 2], alpha);

    T[0] = lerp(precomputed_data[(refr_offset + index) * 3], precomputed_data[(refr_offset + index + 1) * 3], alpha);
    T[1] = lerp(precomputed_data[(refr_offset + index) * 3 + 1], precomputed_data[(refr_offset + index + 1) * 3 + 1], alpha);
    T[2] = lerp(precomputed_data[(refr_offset + index) * 3 + 2], precomputed_data[(refr_offset + index + 1) * 3 + 2], alpha);
  }

  if (intIOR.im > 0.001)
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
    if (rands.x * (sum(R) + sum(T)) < sum(R))
    {
      float3 wo = float3(-wi.x, -wi.y, wi.z);
      pRes->val = R;
      pRes->pdf = sum(R) / (sum(R) + sum(T));
      pRes->dir = normalize(wo.x * s + wo.y * t + wo.z * n);
      pRes->flags |= RAY_EVENT_S;
      pRes->ior = _extIOR;
    }
    else
    {
      float4 fr = FrDielectricDetailedV2(wi.z, ior);
      const float cosThetaT = fr.y;
      const float eta_ti = fr.w;  

      float3 wo = refract(wi, cosThetaT, eta_ti);
      pRes->val = T;
      pRes->pdf = sum(T) / (sum(R) + sum(T));
      pRes->dir = normalize(wo.x * s + wo.y * t + wo.z * n);
      pRes->flags |= (RAY_EVENT_S | RAY_EVENT_T);
      pRes->ior = (_extIOR == intIOR.re) ? extIOR : intIOR.re;
    }
  }

  pRes->val /= std::max(std::abs(dot(pRes->dir, n)), 1e-6f);
}

static inline void filmRoughSampleAndEvalPrecomputed(const Material* a_materials, 
        const float extIOR, const complex intIOR, const float4 a_wavelengths, const float _extIOR,
        float4 rands, float3 v, float3 n, float2 tc, float3 alpha_tex, BsdfSample* pRes, const float* precomputed_data, const bool spectral_mode)
{
  bool reversed = false;
  uint32_t refl_offset;
  uint32_t refr_offset;
  if ((pRes->flags & RAY_FLAG_HAS_INV_NORMAL) != 0) // inside of object
  {
    n = -1 * n;
  }

  if (dot(v, n) < 0.f && intIOR.im < 0.001)
  {
    reversed = true;
    refl_offset = FILM_ANGLE_RES * 2;
    refr_offset = FILM_ANGLE_RES * 3;
  }
  else
  {
    refl_offset = 0;
    refr_offset = FILM_ANGLE_RES;
  }

  const float2 alpha = float2(min(a_materials[0].data[FILM_ROUGH_V], alpha_tex.x), 
                              min(a_materials[0].data[FILM_ROUGH_U], alpha_tex.y));

  float3 s, t = n;
  CoordinateSystemV2(n, &s, &t);
  float3 wi = float3(dot(v, s), dot(v, t), dot(v, n));

  float ior = intIOR.re / extIOR;
  if (reversed)
  {
    wi = -1 * wi;
    ior = 1.f / ior;
  }

  const float4 wm_pdf = sample_visible_normal(wi, {rands.x, rands.y}, alpha);
  const float3 wm = to_float3(wm_pdf);
  if(wm_pdf.w == 0.0f) // not in the same hemisphere
  {
    return;
  }

  float cosThetaI = clamp(fabs(dot(wi, wm)), 0.00001, 1.0f);
  float4 R = float4(0.0f), T = float4(0.0f);

  if (spectral_mode)
  {
    float w = clamp((a_wavelengths[0] - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN), 0.f, 1.f);
    float theta = clamp(acos(cosThetaI) * 2.f / M_PI, 0.f, 1.f);
    //result.refl = lerp_gather_2d(reflectance, w, theta, FILM_LENGTH_RES, FILM_ANGLE_RES);
    //result.refr = lerp_gather_2d(transmittance, w, theta, FILM_LENGTH_RES, FILM_ANGLE_RES);
    w *= FILM_LENGTH_RES - 1;
    theta *= FILM_ANGLE_RES - 1;
    uint32_t index1 = std::min(uint32_t(w), uint32_t(FILM_LENGTH_RES - 2));
    uint32_t index2 = std::min(uint32_t(theta), uint32_t(FILM_ANGLE_RES - 2));

    float alpha = w - float(index1);
    float beta = theta - float(index2);

    float v0 = lerp(precomputed_data[refl_offset * FILM_LENGTH_RES + index1 * FILM_ANGLE_RES + index2], precomputed_data[refl_offset * FILM_LENGTH_RES + (index1 + 1) * FILM_ANGLE_RES + index2], alpha);
    float v1 = lerp(precomputed_data[refl_offset * FILM_LENGTH_RES + index1 * FILM_ANGLE_RES + index2 + 1], precomputed_data[refl_offset * FILM_LENGTH_RES + (index1 + 1) * FILM_ANGLE_RES + index2 + 1], alpha);
    R[0] = lerp(v0, v1, beta);

    v0 = lerp(precomputed_data[refr_offset * FILM_LENGTH_RES + index1 * FILM_ANGLE_RES + index2], precomputed_data[refr_offset * FILM_LENGTH_RES + (index1 + 1) * FILM_ANGLE_RES + index2], alpha);
    v1 = lerp(precomputed_data[refr_offset * FILM_LENGTH_RES + index1 * FILM_ANGLE_RES + index2 + 1], precomputed_data[refr_offset * FILM_LENGTH_RES + (index1 + 1) * FILM_ANGLE_RES + index2 + 1], alpha);
    T[0] = lerp(v0, v1, beta);
  }
  else
  {
    float theta = clamp(acos(cosThetaI) * 2.f / M_PI, 0.f, 1.f);
    theta *= FILM_ANGLE_RES - 1;
    uint32_t index = std::min(uint32_t(theta), uint32_t(FILM_ANGLE_RES - 2));

    float alpha = theta - float(index);

    R[0] = lerp(precomputed_data[(refl_offset + index) * 3], precomputed_data[(refl_offset + index + 1) * 3], alpha);
    R[1] = lerp(precomputed_data[(refl_offset + index) * 3 + 1], precomputed_data[(refl_offset + index + 1) * 3 + 1], alpha);
    R[2] = lerp(precomputed_data[(refl_offset + index) * 3 + 2], precomputed_data[(refl_offset + index + 1) * 3 + 2], alpha);

    T[0] = lerp(precomputed_data[(refr_offset + index) * 3], precomputed_data[(refr_offset + index + 1) * 3], alpha);
    T[1] = lerp(precomputed_data[(refr_offset + index) * 3 + 1], precomputed_data[(refr_offset + index + 1) * 3 + 1], alpha);
    T[2] = lerp(precomputed_data[(refr_offset + index) * 3 + 2], precomputed_data[(refr_offset + index + 1) * 3 + 2], alpha);
  }
  
  if (intIOR.im > 0.001)
  {
    float3 wo = reflect((-1.0f) * wi, wm);
    if (wi.z < 0.f || wo.z <= 0.f)
    {
      return;
    }
    float G = microfacet_G(wi, wo, wm, alpha);
    pRes->pdf = wm_pdf.w / (4.0f * std::abs(dot(wo, wm)));
    pRes->val = wm_pdf.w / (4.0f * wo.z * wm.z) * smith_g1(wo, wm, alpha) * R / std::max(wo.z, EPSILON_32);
    if (reversed)
    {
      wo = -1 * wo;
    }
    pRes->dir = normalize(wo.x * s + wo.y * t + wo.z * n);
    pRes->flags = RAY_FLAG_HAS_NON_SPEC;
    pRes->ior = _extIOR;
  }
  else
  {
    if (rands.w * (sum(R) + sum(T)) < sum(R))
    {
      float3 wo = reflect((-1.0f) * wi, wm);
      if (wi.z < 0.f || wo.z <= 0.f)
      {
        return;
      }
      pRes->pdf = wm_pdf.w / (4.0f * std::abs(dot(wo, wm))) * sum(R) / (sum(R) + sum(T));
      pRes->val = wm_pdf.w / (4.0f * wo.z * wm.z) * smith_g1(wo, wm, alpha) * R / std::max(wo.z, EPSILON_32);
      if (reversed)
      {
        wo = -1 * wo;
      }
      pRes->dir = normalize(wo.x * s + wo.y * t + wo.z * n);
      pRes->flags = RAY_FLAG_HAS_NON_SPEC;
      pRes->ior = _extIOR;
    }
    else
    {
      float4 fr = FrDielectricDetailedV2(dot(wi, wm), ior);
      const float cosThetaT = fr.y;
      const float eta_it = fr.z;
      const float eta_ti = fr.w;  

      float3 ws, wt;
      CoordinateSystemV2(wm, &ws, &wt);
      const float3 local_wi = {dot(ws, wi), dot(wt, wi), dot(wm, wi)};
      const float3 local_wo = refract(local_wi, cosThetaT, eta_ti);
      float3 wo = local_wo.x * ws + local_wo.y * wt + local_wo.z * wm;
      if (wo.z > 0.f)
      {
        return;
      }
      float G = microfacet_G(wi, wo, wm, alpha);
      float denom = sqr(dot(wo, wm) + dot(wi, wm) / eta_it);
      float dwm_dwi = fabs(dot(wo, wm)) / denom;
      pRes->val = wm_pdf.w * G * T * fabs(dot(wi, wm) * dot(wo, wm) / (wi.z * wo.z * denom));
      pRes->pdf = wm_pdf.w * dwm_dwi * sum(T) / (sum(R) + sum(T));
      if (reversed)
      {
        wo = -1 * wo;
      }
      pRes->dir = normalize(wo.x * s + wo.y * t + wo.z * n);
      pRes->flags = RAY_FLAG_HAS_NON_SPEC;;
      pRes->ior = (_extIOR == intIOR.re) ? extIOR : intIOR.re;
    }
  }
}


static void filmRoughEvalPrecomputed(const Material* a_materials, 
        const float extIOR, const complex intIOR, const float4 a_wavelengths, float3 l, float3 v, float3 n, float2 tc,
        float3 alpha_tex, BsdfEval* pRes, const float* precomputed_data, const bool spectral_mode)
{
  if (intIOR.im < 0.001)
  {
    return;
  }

  uint32_t refl_offset;
  uint32_t refr_offset;

  bool reversed = false;
  if (dot(v, n) < 0.f && intIOR.im < 0.001)
  {
    reversed = true;
    refl_offset = FILM_ANGLE_RES * 2;
    refr_offset = FILM_ANGLE_RES * 3;
  }
  else
  {
    refl_offset = 0;
    refr_offset = FILM_ANGLE_RES;
  }

  const float2 alpha = float2(min(a_materials[0].data[FILM_ROUGH_V], alpha_tex.x), 
                              min(a_materials[0].data[FILM_ROUGH_U], alpha_tex.y));

  float3 s, t = n;
  CoordinateSystemV2(n, &s, &t);
  const float3 wo = float3(dot(l, s), dot(l, t), dot(l, n));
  const float3 wi = float3(dot(v, s), dot(v, t), dot(v, n));
  const float3 wm = normalize(wo + wi);

  if (wi.z * wo.z < 0.f)
  {
    return;
  }

  float ior = intIOR.re / extIOR;
  if (reversed)
  {
    ior = 1.f / ior;
  }

  float cosThetaI = clamp(fabs(dot(wo, wm)), 0.00001, 1.0f);
  
  float4 R = float4(0.0f);
  if (spectral_mode)
  {
    float w = clamp((a_wavelengths[0] - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN), 0.f, 1.f);
    float theta = clamp(acos(cosThetaI) * 2.f / M_PI, 0.f, 1.f);
    //result.refl = lerp_gather_2d(reflectance, w, theta, FILM_LENGTH_RES, FILM_ANGLE_RES);
    w *= FILM_LENGTH_RES - 1;
    theta *= FILM_ANGLE_RES - 1;
    uint32_t index1 = std::min(uint32_t(w), uint32_t(FILM_LENGTH_RES - 2));
    uint32_t index2 = std::min(uint32_t(theta), uint32_t(FILM_ANGLE_RES - 2));

    float alpha = w - float(index1);
    float beta = theta - float(index2);

    float v0 = lerp(precomputed_data[refl_offset * FILM_LENGTH_RES + index1 * FILM_ANGLE_RES + index2], precomputed_data[refl_offset * FILM_LENGTH_RES + (index1 + 1) * FILM_ANGLE_RES + index2], alpha);
    float v1 = lerp(precomputed_data[refl_offset * FILM_LENGTH_RES + index1 * FILM_ANGLE_RES + index2 + 1], precomputed_data[refl_offset * FILM_LENGTH_RES + (index1 + 1) * FILM_ANGLE_RES + index2 + 1], alpha);
    R[0] = lerp(v0, v1, beta);
  }
  else
  {
    float theta = clamp(acos(cosThetaI) * 2.f / M_PI, 0.f, 1.f);
    theta *= FILM_ANGLE_RES - 1;
    uint32_t index = std::min(uint32_t(theta), uint32_t(FILM_ANGLE_RES - 2));

    float alpha = theta - float(index);

    R[0] = lerp(precomputed_data[(refl_offset + index) * 3], precomputed_data[(refl_offset + index + 1) * 3], alpha);
    R[1] = lerp(precomputed_data[(refl_offset + index) * 3 + 1], precomputed_data[(refl_offset + index + 1) * 3 + 1], alpha);
    R[2] = lerp(precomputed_data[(refl_offset + index) * 3 + 2], precomputed_data[(refl_offset + index + 1) * 3 + 2], alpha);
  }

  const float cos_theta_i = std::max(wi.z, EPSILON_32);
  const float cos_theta_o = std::max(wo.z, EPSILON_32);

  float D = eval_microfacet(wm, alpha, 1);
  float G = microfacet_G(wi, wo, wm, alpha);
  pRes->val = D * G * R / (4.0f * cos_theta_i * cos_theta_o);
  pRes->pdf = D * smith_g1(wi, wm, alpha) / (4.0f * cos_theta_i);
}