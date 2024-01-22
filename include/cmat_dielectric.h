#pragma once
#include "cglobals.h"
#include "crandom.h"
#include "cmaterial.h"
#include "../spectrum.h"


static inline void dielectricSmoothSampleAndEval(const Material* a_materials, const float4 etaSpec,
                                                 const float _extIOR, float4 rands, float3 v, float3 n, float2 tc, 
                                                 BsdfSample* pRes)
{
  const float extIOR = a_materials[0].data[DIELECTRIC_ETA_EXT];

  if ((pRes->flags & RAY_FLAG_HAS_INV_NORMAL) != 0) // hit the reverse side of the polygon from the volume
  {
    n = -1 * n;
  }

  float3 s, t = n;
  CoordinateSystemV2(n, &s, &t);
  float3 wi = float3(dot(v, s), dot(v, t), dot(v, n));

  float eta = etaSpec.x / extIOR; // TODO: spectral eta - kill all other wavelengths

  // if ((pRes->flags & RAY_FLAG_HAS_INV_NORMAL) != 0) // hit the reverse side of the polygon from the volume
  // {
  //   if (_extIOR == etaSpec.x) // TODO: spectral eta
  //     eta = 1.0f / etaSpec.x;
  // }

// 0.745328009  0.666666687
  float4 fr = FrDielectricDetailedV2(wi.z, eta); 
  // float4 fr1 = FrDielectricDetailed(wi.z, eta); 
  // float4 fr2 = FrDielectricDetailedV2(wi.z, eta); 
  // if(std::abs(fr1.x - fr2.x) > 1e6f || std::abs(fr1.y - fr2.y) > 1e6f || std::abs(fr1.z - fr2.z) > 1e6f || std::abs(fr1.w - fr2.w) > 1e6f)
  //   int a = 1;


  const float R = fr.x;
  const float cos_theta_t = fr.y;
  const float eta_it = fr.z;
  const float eta_ti = fr.w;  
  const float T = 1 - R;

  if(rands.x < R) // perfect specular reflection
  {
    float3 wo = float3(-wi.x, -wi.y, wi.z);
    pRes->val = float4(R);
    pRes->pdf = R;
    pRes->dir = normalize(wo.x * s + wo.y * t + wo.z * n);
    pRes->flags |= RAY_EVENT_S;
    pRes->ior = _extIOR;
  }
  else // perfect specular transmission
  {
    float3 wo = refract(wi, cos_theta_t, eta_ti);
    pRes->val = float4((eta_ti * eta_ti) * T);
    pRes->pdf = T;
    pRes->dir = normalize(wo.x * s + wo.y * t + wo.z * n);
    pRes->flags |= (RAY_EVENT_S | RAY_EVENT_T);
    pRes->ior = (_extIOR == etaSpec.x) ? extIOR : etaSpec.x;
  }

  pRes->val /= std::max(std::abs(dot(pRes->dir, n)), 1e-6f);
}


static void dielectricSmoothEval(BsdfEval* pRes)
{
  pRes->val = float4(0.0f);
  pRes->pdf = 0.0f;
}