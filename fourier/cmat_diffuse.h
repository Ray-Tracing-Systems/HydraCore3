#ifndef FOURIER_CMAT__DIFFUSE_H_
#define FOURIER_CMAT__DIFFUSE_H_
#include "../include/cglobals.h"
#include "../include/cmaterial.h"

#include "fspec.h"

static inline void diffuseSampleAndEvalF(const Material* a_materials, FourierSpec a_reflSpec, float4 rands, float3 v, 
                                        float3 n, float2 tc, BsdfSampleF* pRes)
{
  const uint   cflags     = a_materials[0].cflags;
  const float3 lambertDir = lambertSample(float2(rands.x, rands.y), v, n);
  const float  lambertPdf = lambertEvalPDF(lambertDir, v, n);
  const float  lambertVal = lambertEvalBSDF(lambertDir, v, n);

  pRes->dir   = lambertDir;
  pRes->val   = lambertVal * a_reflSpec;
  pRes->pdf   = lambertPdf;
  pRes->flags = RAY_FLAG_HAS_NON_SPEC;

        
  if ((cflags & GLTF_COMPONENT_ORENNAYAR) != 0)
    pRes->val *= orennayarFunc(lambertDir, (-1.0f) * v, n, a_materials[0].data[DIFFUSE_ROUGHNESS]);
            
}


static inline void diffuseEvalF(const Material* a_materials, FourierSpec a_reflSpec, float3 l, float3 v, float3 n, float2 tc, BsdfEvalF* res)
{
  const uint cflags = a_materials[0].cflags;
 
  float lambertVal       = lambertEvalBSDF(l, v, n);
  const float lambertPdf = lambertEvalPDF (l, v, n);

  if ((cflags & GLTF_COMPONENT_ORENNAYAR) != 0)
    lambertVal *= orennayarFunc(l, v, n, a_materials[0].data[DIFFUSE_ROUGHNESS]);

  res->val = lambertVal * a_reflSpec; 
  res->pdf = lambertPdf; 
}

#endif