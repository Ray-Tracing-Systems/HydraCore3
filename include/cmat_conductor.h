#pragma once
#include "cglobals.h"
#include "crandom.h"
#include "cmaterial.h"

static inline void conductorSmoothSampleAndEval(const GLTFMaterial* a_materials, float4 rands, float3 v, float3 n, float2 tc, 
                                                float3 color,
                                                BsdfSample* pRes)
{
  const uint  cflags  = a_materials[0].cflags;
  const float eta     = a_materials[0].metalColor[2];
  const float k       = a_materials[0].metalColor[3];
  
  const float3 pefReflDir = reflect((-1.0f)*v, n);
  const float cosThetaOut = dot(pefReflDir, n);
  float3 dir = pefReflDir;
  float  pdf = 1.0f;
  float  val = FrComplexConductor(cosThetaOut, complex{eta, k});
  
  val = (cosThetaOut <= 1e-6f) ? 0.0f : (val / std::max(cosThetaOut, 1e-6f));  // BSDF is multiplied (outside) by cosThetaOut. For mirrors this shouldn't be done, so we pre-divide here instead.

  pRes->val = float3(val, val, val); 
  pRes->dir = dir;
  pRes->pdf = pdf;
  pRes->flags = RAY_EVENT_S;
}


static void conductorSmoothEval(const GLTFMaterial* a_materials, float3 l, float3 v, float3 n, float2 tc, 
                                float3 color, 
                                BsdfEval* res)
{
  res->color = {0.0f, 0.0f, 0.0f};
  res->pdf = 0.0f;
}