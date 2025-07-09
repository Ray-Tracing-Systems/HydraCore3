#pragma once
#include "../include/cglobals.h"
#include "../include/crandom.h"
#include "../include/cmaterial.h"
#include "../include/airy_reflectance.h"
#include "../include/cmat_film.h"

#include "fspec.h"

#include <cassert>

#include <iostream>

static inline void filmSmoothSampleAndEvalF(const Material* a_materials, 
        const float extIOR, const complex intIOR, const float _extIOR,
        float4 rands, float3 v, float3 n, float2 tc, BsdfSampleF* pRes, const float* precomputed_data)
{
  assert(!FourierSpec::unpack_on_multiply);

  const uint transparFlag = as_uint((a_materials[0].data[FILM_TRANSPARENT]));
  bool reversed = false;
  uint32_t refl_offset;
  uint32_t refr_offset;
  if ((pRes->flags & RAY_FLAG_HAS_INV_NORMAL) != 0) // inside of object
  {
    n = -1 * n;
  }
  if (dot(n, v) < 0.f && intIOR.im < 0.001f)
  {
    reversed = true;
    refl_offset = FILM_ANGLE_RES * FourierSpec::SIZE * 2;
    refr_offset = FILM_ANGLE_RES * FourierSpec::SIZE * 3;
  }
  else
  {
    refl_offset = 0;
    refr_offset = FILM_ANGLE_RES * FourierSpec::SIZE;
  }

  float3 s, t = n;
  CoordinateSystemV2(n, &s, &t);
  float3 wi = float3(dot(v, s), dot(v, t), dot(v, n));

  float cosThetaI = clamp(std::abs(wi.z), 0.0001f, 1.0f);
  float ior = intIOR.re / extIOR;
  
  //float4 R = float4(0.0f), T = float4(0.0f);
  float theta = clamp(std::acos(cosThetaI) * 2.f / M_PI, 0.f, 1.f);
  theta *= FILM_ANGLE_RES - 1;
  uint32_t index = std::min(uint32_t(theta), uint32_t(FILM_ANGLE_RES - 2));

  float alpha = theta - float(index);

  FourierSpec R, T;
  for (uint i = 0; i < FourierSpec::SIZE; ++i)
  {
    R[i] = lerp(precomputed_data[refl_offset + FourierSpec::SIZE * index + i], precomputed_data[refl_offset + FourierSpec::SIZE * (index + 1) + i], alpha);
    T[i] = lerp(precomputed_data[refr_offset + FourierSpec::SIZE * index + i], precomputed_data[refr_offset + FourierSpec::SIZE * (index + 1) + i], alpha);
  }

  if (intIOR.im > 0.001f || transparFlag == 0)
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

    if (rands.x < 0.1) // rough estimation, needs better heuristic
    {
      float3 wo = float3(-wi.x, -wi.y, wi.z);
      pRes->val = R;
      pRes->pdf = 0.1;
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
      pRes->pdf = 0.9;
      pRes->dir = normalize(wo.x * s + wo.y * t + wo.z * n);
      pRes->flags |= (RAY_EVENT_S | RAY_EVENT_T);
      pRes->ior = (_extIOR == intIOR.re) ? extIOR : intIOR.re;
    }
  }

  pRes->val /= std::max(std::abs(dot(pRes->dir, n)), 1e-6f);
}
