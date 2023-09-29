#pragma once
#include "cglobals.h"
#include "crandom.h"
#include "cmaterial.h"


static inline void conductorSmoothSampleAndEval(const Material* a_materials, float4 rands, float3 v, float3 n, float2 tc,
                                                float3 color,
                                                BsdfSample* pRes)
{
  const uint  cflags = as_uint(a_materials[0].data[UINT_CFLAGS]);
  const float eta    = a_materials[0].data[CONDUCTOR_ETA];
  const float k      = a_materials[0].data[CONDUCTOR_K];
  
  const float3 pefReflDir = reflect((-1.0f)*v, n);
  const float cosThetaOut = dot(pefReflDir, n);
  float3 dir              = pefReflDir;
  float  pdf              = 1.0f;
  float  val              = FrComplexConductor(cosThetaOut, complex{eta, k});
  
  val = (cosThetaOut <= 1e-6f) ? 0.0f : (val / std::max(cosThetaOut, 1e-6f));  // BSDF is multiplied (outside) by cosThetaOut. For mirrors this shouldn't be done, so we pre-divide here instead.

  pRes->val = float3(val, val, val); 
  pRes->dir = dir;
  pRes->pdf = pdf;
  pRes->flags = RAY_EVENT_S;
}


static void conductorSmoothEval(const Material* a_materials, float3 l, float3 v, float3 n, float2 tc,
                                float3 color, 
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


static inline void conductorRoughSampleAndEval(const Material* a_materials, float4 rands, float3 v, float3 n, float2 tc,
                                                float3 color,
                                                BsdfSample* pRes)
{
  if(v.z == 0)
    return;

  const uint  cflags = as_uint(a_materials[0].data[UINT_CFLAGS]);
  const float eta    = a_materials[0].data[CONDUCTOR_ETA];
  const float k      = a_materials[0].data[CONDUCTOR_K];
  const float2 alpha = float2(a_materials[0].data[CONDUCTOR_ROUGH_V], a_materials[0].data[CONDUCTOR_ROUGH_U]);

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

  float val = conductorRoughEvalInternal(wo, wi, wm, alpha, complex{eta, k});

  pRes->val   = float3(val, val, val); 
  pRes->dir   = normalize(wi.x * nx + wi.y * ny + wi.z * nz);
  pRes->pdf   = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
  pRes->flags = RAY_FLAG_HAS_NON_SPEC;
}


static void conductorRoughEval(const Material* a_materials, float3 l, float3 v, float3 n, float2 tc,
                                float3 color, 
                                BsdfEval* pRes)
{
  const uint  cflags = as_uint(a_materials[0].data[UINT_CFLAGS]);
  const float eta    = a_materials[0].data[CONDUCTOR_ETA];
  const float k      = a_materials[0].data[CONDUCTOR_K];
  const float2 alpha = float2(a_materials[0].data[CONDUCTOR_ROUGH_V], a_materials[0].data[CONDUCTOR_ROUGH_U]);

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

  float val = conductorRoughEvalInternal(wo, wi, wm, alpha, complex{eta, k});


  pRes->val = float3(val, val, val);

  wm        = FaceForward(wm, float3(0.0f, 0.0f, 1.0f));
  pRes->pdf = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
}

static inline void conductorSampleAndEval(const Material* a_materials, float4 rands, float3 v, float3 n, float2 tc,
                                          float3 color,
                                          BsdfSample* pRes)
{
  const float2 alpha = float2(a_materials[0].data[CONDUCTOR_ROUGH_V], a_materials[0].data[CONDUCTOR_ROUGH_U]);

  if(trEffectivelySmooth(alpha))
  {
    conductorSmoothSampleAndEval(a_materials, rands, v, n, tc, color, pRes);
  }
  else
  {
    conductorRoughSampleAndEval(a_materials, rands, v, n, tc, color, pRes);
  }
}

static inline void conductorEval(const Material* a_materials, float3 l, float3 v, float3 n, float2 tc,
                                float3 color,
                                BsdfEval* pRes)
{
  const float2 alpha = float2(a_materials[0].data[CONDUCTOR_ROUGH_V], a_materials[0].data[CONDUCTOR_ROUGH_U]);

  if(trEffectivelySmooth(alpha))
  {
    conductorSmoothEval(a_materials, l, v, n, tc, color, pRes);
  }
  else
  {
    conductorRoughEval(a_materials, l, v, n, tc, color, pRes);
  }
}