#pragma once
#include "cglobals.h"
#include "crandom.h"
#include "cmaterial.h"


static inline float3 sampleDiffuseSpectrum(const Material* a_materials, const Spectrum* a_spectra, float3 a_wavelengths)
{  
  float3 reflSpec = float3(a_materials[0].data[DIFFUSE_COLOR]);
  if(a_wavelengths.M[0] == 0.0f)
    return reflSpec;

  const uint specId = as_uint(a_materials[0].data[DIFFUSE_SPECID]);

  if(specId < 0xFFFFFFFF)
  {
    const uint spectralSamples = uint(sizeof(a_wavelengths.M) / sizeof(a_wavelengths.M[0])); 
    for(uint i = 0; i < spectralSamples; ++i)
      reflSpec[i] = a_spectra[specId].Sample(a_wavelengths[i]);
  }

  return reflSpec;
}

static inline void diffuseSampleAndEval(const Material* a_materials, const Spectrum* a_spectra, float3 a_wavelengths, float4 rands, float3 v, 
                                        float3 n, float2 tc, float3 color, BsdfSample* pRes)
{
  const uint   cflags       = as_uint(a_materials[0].data[UINT_CFLAGS]);
  const float3 lambertDir   = lambertSample(float2(rands.x, rands.y), v, n);
  const float  lambertPdf   = lambertEvalPDF(lambertDir, v, n);
  const float  lambertVal   = lambertEvalBSDF(lambertDir, v, n);

  const float3 reflSpec = sampleDiffuseSpectrum(a_materials, a_spectra, a_wavelengths);

  pRes->dir       = lambertDir;
  pRes->val       = lambertVal * reflSpec;
  pRes->pdf       = lambertPdf;
  pRes->flags     = RAY_FLAG_HAS_NON_SPEC;
        
  if ((cflags & GLTF_COMPONENT_ORENNAYAR) != 0)
    pRes->val *= orennayarFunc(lambertDir, (-1.0f) * v, n, a_materials[0].data[DIFFUSE_ROUGHNESS]);
            
}


static void diffuseEval(const Material* a_materials, const Spectrum* a_spectra, float3 a_wavelengths, float3 l, float3 v, float3 n, float2 tc, 
                        float3 color, BsdfEval* res)
{
  const uint   cflags    = as_uint(a_materials[0].data[UINT_CFLAGS]);
  const float3 reflSpec  = sampleDiffuseSpectrum(a_materials, a_spectra, a_wavelengths);
 
  float lambertVal       = lambertEvalBSDF(l, v, n);
  const float lambertPdf = lambertEvalPDF (l, v, n);

  if ((cflags & GLTF_COMPONENT_ORENNAYAR) != 0)
    lambertVal *= orennayarFunc(l, v, n, a_materials[0].data[DIFFUSE_ROUGHNESS]);

  res->val = lambertVal * reflSpec; 
  res->pdf = lambertPdf; 
}