#pragma once

#include "../include/cglobals.h"
#include "../include/crandom.h"
#include "../include/cmaterial.h"

#include "fspec.h"


static inline void conductorSmoothSampleAndEvalF(const Material* a_materials,
                                                float4 rands, float3 v, float3 n, float2 tc, 
                                                const float* precomputed_data, BsdfSampleF* pRes)
{
  const float3 pefReflDir = reflect((-1.0f)*v, n);
  const float cosThetaOut = dot(pefReflDir, n);
  float3 dir              = pefReflDir;
  float  pdf              = 1.0f;

  float theta = clamp(std::acos(cosThetaOut) * 2.f / M_PI, 0.f, 1.f);
  theta *= FOURIER_ANGLE_RES - 1;
  uint32_t index = std::min(uint32_t(theta), uint32_t(FILM_ANGLE_RES - 2));

  float alpha = theta - float(index);

  FourierSpec val;
  for (uint i = 0; i < FourierSpec::SIZE; ++i)
  {
    val[i] = lerp(precomputed_data[FourierSpec::SIZE * index + i], precomputed_data[FourierSpec::SIZE * (index + 1) + i], alpha);
    val[i] = (cosThetaOut <= 1e-6f) ? 0.0f : (val[i] / std::max(cosThetaOut, 1e-6f)); 
  }
  
  pRes->val = val; 
  pRes->dir = dir;
  pRes->pdf = pdf;
  pRes->flags = RAY_EVENT_S;
}


static FourierSpec conductorRoughEvalInternalF(float3 wo, float3 wi, float3 wm, float2 alpha, const float* precomputed_data)
{
  if(wo.z * wi.z < 0) // not in the same hemisphere
    return FourierSpec();

  float cosTheta_o = AbsCosTheta(wo);
  float cosTheta_i = AbsCosTheta(wi);
  if (cosTheta_i == 0 || cosTheta_o == 0)
    return FourierSpec();

  auto cosThetaI = std::abs(dot(wo, wm));

  float theta = clamp(std::acos(cosThetaI) * 2.f / M_PI, 0.f, 1.f);
  theta *= FOURIER_ANGLE_RES - 1;
  uint32_t index = std::min(uint32_t(theta), uint32_t(FILM_ANGLE_RES - 2));

  float beta = theta - float(index);

  FourierSpec F;
  for (uint i = 0; i < FourierSpec::SIZE; ++i)
  {
    F[i] = lerp(precomputed_data[FourierSpec::SIZE * index + i], precomputed_data[FourierSpec::SIZE * (index + 1) + i], beta);
  }

  FourierSpec val = trD(wm, alpha) * F * trG(wo, wi, alpha) / (4.0f * cosTheta_i * cosTheta_o);

  return val;
}


static inline void conductorRoughSampleAndEvalF(const Material* a_materials, float4 rands, float3 v, float3 n, float2 tc, 
                                                float3 alpha_tex, const float* precomputed_data, BsdfSampleF* pRes)
{
  if(v.z == 0)
    return;

  const float2 alpha = float2(min(a_materials[0].data[CONDUCTOR_ROUGH_U], alpha_tex.x), 
                              min(a_materials[0].data[CONDUCTOR_ROUGH_V], alpha_tex.y));

  float3 nx, ny, nz = n;
  CoordinateSystemV2(nz, &nx, &ny);
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

  FourierSpec val = conductorRoughEvalInternalF(wo, wi, wm, alpha, precomputed_data);

  pRes->val   = val; 
  pRes->dir   = normalize(wi.x * nx + wi.y * ny + wi.z * nz);
  pRes->pdf   = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
  pRes->flags = RAY_FLAG_HAS_NON_SPEC;
}


static void conductorRoughEvalF(const Material* a_materials, float3 l, float3 v, float3 n, float2 tc, 
                                float3 alpha_tex, const float* precomputed_data, BsdfEvalF* pRes)
{
  const float2 alpha = float2(min(a_materials[0].data[CONDUCTOR_ROUGH_U], alpha_tex.x), 
                              min(a_materials[0].data[CONDUCTOR_ROUGH_V], alpha_tex.y));

  float3 nx, ny, nz = n;
  CoordinateSystemV2(nz, &nx, &ny);

  // v = (-1.0f) * v;
  const float3 wo = float3(dot(v, nx), dot(v, ny), dot(v, nz));
  const float3 wi = float3(dot(l, nx), dot(l, ny), dot(l, nz));

  if(wo.z * wi.z < 0.0f)
    return;

  float3 wm = wo + wi;
  if (dot(wm, wm) == 0)
      return;

  wm = normalize(wm);
  FourierSpec val = conductorRoughEvalInternalF(wo, wi, wm, alpha, precomputed_data);

  pRes->val = val;

  wm        = FaceForward(wm, float3(0.0f, 0.0f, 1.0f));
  pRes->pdf = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
}