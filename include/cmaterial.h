#ifndef RTC_MATERIAL
#define RTC_MATERIAL

#include "cglobals.h"

struct BsdfSample
{
  float3 color;
  float3 direction;
  float  pdf; 
  int    flags;
};

struct BsdfEval
{
  float3 color;
  float  pdf; 
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////

enum BRDF_TYPES { BRDF_TYPE_LAMBERT         = 1, 
                  BRDF_TYPE_GGX             = 2, 
                  BRDF_TYPE_GLTF            = 5,
                  BRDF_TYPE_GLASS           = 6,
                  BRDF_TYPE_MIRROR          = 7,
                  BRDF_TYPE_LIGHT_SOURCE = 0xEFFFFFFF };

// The BRDF of the metallic-roughness material is a linear interpolation of a metallic BRDF and a dielectric BRDF. 
// The BRDFs **share** the parameters for roughness and base color.

struct GLTFMaterial
{
  float4 row0[1];     ///< texture matrix
  float4 row1[1];     ///< texture matrix
  uint   texId[4];    ///< texture id

  float4 baseColor;   ///< color for both lambert and emissive lights; baseColor.w store emission
  float4 metalColor;  ///< in our implementation we allow different color for metals and diffuse
  float4 coatColor;   ///< in our implementation we allow different color for coating (fresnel) and diffuse

  uint  brdfType;     ///<
  uint  lightId;      ///< identifier of light if this material is light  
  float alpha;        ///< blend factor between lambert and reflection : alpha*baseColor + (1.0f-alpha)*baseColor
  float glosiness;    ///< material glosiness or intensity for lights, take color from baseColor

  float ior;
  float dummy1;
  float dummy2;
  float dummy3;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline float3 lambertSample(float2 rands, float3 v, float3 n)
{
   return MapSampleToCosineDistribution(rands.x, rands.y, n, n, 1.0f);
}

static inline float lambertEvalPDF(float3 l, float3 v, float3 n) 
{ 
  return std::abs(dot(l, n)) * INV_PI;
}

static inline float lambertEvalBSDF(float3 l, float3 v, float3 n)
{
  return INV_PI;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline float3 SphericalDirectionPBRT(const float sintheta, const float costheta, const float phi) 
{ 
  return float3(sintheta * cos(phi), sintheta * sin(phi), costheta); 
}

static inline float GGX_Distribution(const float cosThetaNH, const float alpha)
{
  const float alpha2  = alpha * alpha;
  const float NH_sqr  = clamp(cosThetaNH * cosThetaNH, 0.0f, 1.0f);
  const float den     = NH_sqr * alpha2 + (1.0f - NH_sqr);
  return alpha2 / std::max((float)(M_PI) * den * den, 1e-6f);
}

static inline float GGX_GeomShadMask(const float cosThetaN, const float alpha)
{
  // Height - Correlated G.
  //const float tanNV      = sqrt(1.0f - cosThetaN * cosThetaN) / cosThetaN;
  //const float a          = 1.0f / (alpha * tanNV);
  //const float lambda     = (-1.0f + sqrt(1.0f + 1.0f / (a*a))) / 2.0f;
  //const float G          = 1.0f / (1.0f + lambda);

  // Optimized and equal to the commented-out formulas on top.
  const float cosTheta_sqr = clamp(cosThetaN*cosThetaN, 0.0f, 1.0f);
  const float tan2         = (1.0f - cosTheta_sqr) / std::max(cosTheta_sqr, 1e-6f);
  const float GP           = 2.0f / (1.0f + std::sqrt(1.0f + alpha * alpha * tan2));
  return GP;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline float3 ggxSample(float2 rands, float3 v, float3 n, float roughness)
{
  const float roughSqr = roughness * roughness;
    
  float3 nx, ny, nz = n;
  CoordinateSystem(nz, &nx, &ny);
    
  const float3 wo       = float3(dot(v, nx), dot(v, ny), dot(v, nz));
  const float phi       = rands.x * M_TWOPI;
  const float cosTheta  = clamp(std::sqrt((1.0f - rands.y) / (1.0f + roughSqr * roughSqr * rands.y - rands.y)), 0.0f, 1.0f);
  const float sinTheta  = std::sqrt(1.0f - cosTheta * cosTheta);
  const float3 wh       = SphericalDirectionPBRT(sinTheta, cosTheta, phi);
    
  const float3 wi = 2.0f * dot(wo, wh) * wh - wo;      // Compute incident direction by reflecting about wm  
  return normalize(wi.x * nx + wi.y * ny + wi.z * nz); // back to normal coordinate system
}

static inline float ggxEvalPDF(float3 l, float3 v, float3 n, float roughness) 
{ 
  const float dotNV = dot(n, v);
  const float dotNL = dot(n, l);
  if (dotNV < 1e-6f || dotNL < 1e-6f)
    return 1.0f;

  const float  roughSqr  = roughness * roughness;
    
  const float3 h    = normalize(v + l); // half vector.
  const float dotNH = dot(n, h);
  const float dotHV = dot(h, v);
  const float D     = GGX_Distribution(dotNH, roughSqr);
  return  D * dotNH / (4.0f * std::max(dotHV,1e-6f));
}

static inline float ggxEvalBSDF(float3 l, float3 v, float3 n, float roughness)
{
  if(std::abs(dot(l, n)) < 1e-5f)
    return 0.0f; 
 
  const float dotNV = dot(n, v);  
  const float dotNL = dot(n, l);
  if (dotNV < 1e-6f || dotNL < 1e-6f)
    return 0.0f; 

  const float  roughSqr = roughness * roughness;
  const float3 h    = normalize(v + l); // half vector.
  const float dotNH = dot(n, h);
  const float D     = GGX_Distribution(dotNH, roughSqr);
  const float G     = GGX_GeomShadMask(dotNV, roughSqr)*GGX_GeomShadMask(dotNL, roughSqr);      

  return (D * G / std::max(4.0f * dotNV * dotNL, 1e-6f));  // Pass single-scattering
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline float3 gltfConductorFresnel(float3 f0, float VdotH) 
{
  //const float tmp = 1.0f - std::abs(VdotH);
  const float tmp = std::abs(VdotH);
  return f0 + (float3(1.0f,1.0f,1.0f) - f0) * (tmp*tmp*tmp*tmp*tmp);
}

static inline float fresnelSlick(float VdotH)
{
  const float tmp = std::abs(VdotH);
  return (tmp*tmp)*(tmp*tmp)*tmp;
}

static inline float fresnelConductorHydra(float cosTheta, float eta, float roughness)
{
  float tmp = (eta*eta + roughness*roughness) * (cosTheta * cosTheta);
  float rParl2 = (tmp - (eta * (2.0f * cosTheta)) + 1.0f) / (tmp + (eta * (2.0f * cosTheta)) + 1.0f);
  float tmpF = eta*eta + roughness*roughness;
  float rPerp2 = (tmpF - (eta * (2.0f * cosTheta)) + (cosTheta*cosTheta)) / (tmpF + (eta * (2.0f * cosTheta)) + (cosTheta*cosTheta));
  return (rParl2 + rPerp2) / 2.0f;
}

//  The following functions calculate the reflected and refracted 
//	directions in addition to the fresnel coefficients. Based on PBRT
//	and the paper "Derivation of Refraction Formulas" by Paul S. Heckbert.
//

static inline float fresnelDielectric(float cosTheta1, float cosTheta2, float etaExt, float etaInt)
{
  float Rs = (etaExt * cosTheta1 - etaInt * cosTheta2) / (etaExt * cosTheta1 + etaInt * cosTheta2);
  float Rp = (etaInt * cosTheta1 - etaExt * cosTheta2) / (etaInt * cosTheta1 + etaExt * cosTheta2);

  return (Rs * Rs + Rp * Rp) / 2.0f;
}

static inline float fresnelReflectionCoeff(float cosTheta1, float etaExt, float etaInt)
{
  // Swap the indices of refraction if the interaction starts
  // at the inside of the object
  //
  if (cosTheta1 < 0.0f)
  {
    float temp = etaInt;
    etaInt = etaExt;
    etaExt = temp;
  }

  // Using Snell's law, calculate the sine of the angle
  // between the transmitted ray and the surface normal 
  //
  float sinTheta2 = etaExt / etaInt * std::sqrt(std::max(0.0f, 1.0f - cosTheta1*cosTheta1));

  if (sinTheta2 > 1.0f)
    return 1.0f;  // Total internal reflection!

  // Use the sin^2+cos^2=1 identity - max() guards against
  //	numerical imprecision
  //
  float cosTheta2 = std::sqrt(std::max(0.0f, 1.0f - sinTheta2*sinTheta2));

  // Finally compute the reflection coefficient
  //
  return fresnelDielectric(std::abs(cosTheta1), cosTheta2, etaInt, etaExt);
}

static inline float3 colorExtrusionStrong(float3 f0)
{
  const float maxVal = max(max(f0.x,f0.y), max(f0.z, 5e-6f));
  return f0*(1.0f/maxVal);
}

static inline float3 conductorFresnel(float3 f0, float VdotH, float ior, float roughness) 
{
  if(ior == 0.0f) // fresnel reflactance is disabled
    return f0;

  //const float fresnelCoeff = fresnelSlick(VdotH);
  const float fresnelCoeff = fresnelConductorHydra(VdotH, ior, roughness);
  //const float fresnelCoeff = clamp(fresnelReflectionCoeff(std::abs(VdotH), 1.0f, ior), 0.0f, 1.0f);
  return (f0 + (float3(1.0f,1.0f,1.0f) - f0) * fresnelCoeff)*colorExtrusionStrong(f0);
}

static inline float gltfFresnelMix(float VdotH) 
{
  //return 0.25f;
  //const float f1  = (1.0f-ior)/(1+ior);
  //const float f0  = f1*f1;
  // Note that the dielectric index of refraction ior = 1.5 is now f0 = 0.04
  const float tmp = 1.0f - std::abs(VdotH);
  return 0.04f + 0.96f*(tmp*tmp*tmp*tmp*tmp);
}

#endif