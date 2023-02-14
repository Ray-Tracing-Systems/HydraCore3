#ifndef RTC_MATERIAL
#define RTC_MATERIAL

#include "cglobals.h"

struct BsdfSample
{
  float3 color;
  float3 direction;
  float  pdf; 
  uint   flags;
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

enum MATERIAL_EVENT {
  RAY_EVENT_S         = 1,  ///< Indicates Specular reflection or refraction (check for RAY_EVENT_T)
  RAY_EVENT_D         = 2,  ///< Indicates Diffuse  reflection or translucent (check for RAY_EVENT_T)
  RAY_EVENT_G         = 4,  ///< Indicates GLossy   reflection or refraction (check for RAY_EVENT_T)
  RAY_EVENT_T         = 8,  ///< Indicates Transparensy or reftacrion. 
  RAY_EVENT_V         = 16, ///< Indicates Volume scattering, not used for a while
  RAY_EVENT_TOUT      = 32, ///< Indicates Transparensy Outside of water or glass or e.t.c. (old RAY_IS_INSIDE_TRANSPARENT_OBJECT = 128)
  RAY_EVENT_TNINGLASS = 64,
};

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

static inline float FrDielectricPBRT(float cosThetaI, float etaI, float etaT) 
{
  cosThetaI = clamp(cosThetaI, -1.0f, 1.0f);
  // Potentially swap indices of refraction
  bool entering = cosThetaI > 0.f;
  if (!entering) 
  {
    const float tmp = etaI;
    etaI = etaT;
    etaT = tmp;
    cosThetaI = std::abs(cosThetaI);
  }

  // Compute _cosThetaT_ using Snell's law
  float sinThetaI = std::sqrt(std::max(0.0f, 1.0f - cosThetaI * cosThetaI));
  float sinThetaT = etaI / etaT * sinThetaI;

  // Handle total internal reflection
  if (sinThetaT >= 1.0f) 
    return 1.0f;

  const float cosThetaT = std::sqrt(std::max(0.0f, 1.0f - sinThetaT * sinThetaT));
  const float Rparl     = ((etaT * cosThetaI) - (etaI * cosThetaT)) / ((etaT * cosThetaI) + (etaI * cosThetaT));
  const float Rperp     = ((etaI * cosThetaI) - (etaT * cosThetaT)) / ((etaI * cosThetaI) + (etaT * cosThetaT));
  return 0.5f*(Rparl * Rparl + Rperp * Rperp);
}

//static inline float fresnelConductor(float cosTheta, float eta, float roughness)
//{
//  float tmp = (eta*eta + roughness*roughness) * (cosTheta * cosTheta);
//  float rParl2 = (tmp - (eta * (2.0f * cosTheta)) + 1.0f) / (tmp + (eta * (2.0f * cosTheta)) + 1.0f);
//  float tmpF = eta*eta + roughness*roughness;
//  float rPerp2 = (tmpF - (eta * (2.0f * cosTheta)) + (cosTheta*cosTheta)) / (tmpF + (eta * (2.0f * cosTheta)) + (cosTheta*cosTheta));
//  return 0.5f*(rParl2 + rPerp2);
//}


static inline float fresnelSlick(float VdotH)
{
  const float tmp = 1.0f - std::abs(VdotH);
  return (tmp*tmp)*(tmp*tmp)*tmp;
}

static inline float3 hydraFresnelCond(float3 f0, float VdotH, float ior, float roughness) 
{
  if(ior == 0.0f) // fresnel reflactance is disabled
    return f0;

  return f0 + (float3(1.0f,1.0f,1.0f) - f0) * fresnelSlick(VdotH); // return bsdf * (f0 + (1 - f0) * (1 - abs(VdotH))^5)
}


static inline float pow5(float x) { return (x*x)*(x*x)*x; }


//// BSDF Inline Functions
//inline float CosTheta   (float3 w, float3 n) { return dot(w,n); }
//inline float Cos2Theta  (float3 w, float3 n) { float z = dot(w,n); return z * z; }
//inline float AbsCosTheta(float3 w, float3 n) { return std::abs(CosTheta(w,n)); }
//inline float Sin2Theta  (float3 w, float3 n) { return std::max(0.0f, 1.0f - Cos2Theta(w,n)); }
//
//inline float SinTheta    (float3 w, float3 n) { return std::sqrt(Sin2Theta(w,n)); }
//inline float TanTheta    (float3 w, float3 n) { return SinTheta(w,n) / CosTheta(w,n); }
//inline float Tan2Theta   (float3 w, float3 n) { return Sin2Theta(w,n) / Cos2Theta(w,n); }
//
//static inline float CosPhi(float3 w, float3 n) 
//{
//  float sinTheta = SinTheta(w,n);
//  return (sinTheta == 0) ? 1 : clamp(w.x / sinTheta, -1.0f, 1.0f);
//}
//
//static inline float SinPhi(float3 w, float3 n) 
//{
//  float sinTheta = SinTheta(w,n);
//  return (sinTheta == 0) ? 0 : clamp(w.y / sinTheta, -1.0f, 1.0f);
//}
//
//static inline float Cos2Phi(float3 w, float3 n) { return CosPhi(w,n) * CosPhi(w,n); }
//static inline float Sin2Phi(float3 w, float3 n) { return SinPhi(w,n) * SinPhi(w,n); }
//
//static inline float PBRT_GGX_D(const float3 wh, const float3 n, float roughness) 
//{
//  float tan2Theta = Tan2Theta(wh,n);
//  if (std::isinf(tan2Theta)) 
//    return 0.0f;
//  const float cos4Theta = Cos2Theta(wh,n) * Cos2Theta(wh,n);
//  const float e         = (Cos2Phi(wh,n) / (roughness * roughness) + Sin2Phi(wh,n) / (roughness * roughness)) * tan2Theta;
//  return 1.0f / (M_PI * roughness * roughness * cos4Theta * (1.0f + e) * (1.0f + e));
//}

static inline float3 pbrtFresnelBlendBRDF(float3 Rd, float3 Rs, float3 l, float3 v, float3 n, float roughness, float a_ggxVal) 
{
  const float  cosThetaL = std::abs(dot(l,n)); 
  const float  cosThetaV = std::abs(dot(v,n));
  const float  diffMult  = (28.f/(23.f *(float)(M_PI)))*(1 - pow5(1 - .5f*cosThetaV))*(1 - pow5(1 - .5f * cosThetaL)); 

  // const float3 wh = l + v;
  // if (wh.x == 0 && wh.y == 0 && wh.z == 0) 
  //   return float3(0,0,0);
  // //const float3 schlickFresnel = Rs + pow5(1 - dot(v, wh)) * (float3(1.0f) - Rs);
  // //(D * G / std::max(4.0f * dotNV * dotNL, 1e-6f));
  // const float3 specular = PBRT_GGX_D(wh,n,roughness) /(4.0f * abs(dot(v, wh)) * std::max(cosThetaV, cosThetaL))*float3(1,1,1); // * schlickFresnel;
  return Rd*diffMult; // + specular;
}

static inline float pbrtFresnelDiffuseMult(float3 l, float3 v, float3 n) 
{
  const float  cosThetaL = std::abs(dot(l,n)); 
  const float  cosThetaV = std::abs(dot(v,n));
  const float  diffMult  = (28.f/(23.f *(float)(M_PI)))*(1 - pow5(1 - .5f*cosThetaV))*(1 - pow5(1 - .5f * cosThetaL)); 
  return diffMult;
}

static inline float hydraFresnelDiel(float VdotH, float ior, float roughness) 
{
  return FrDielectricPBRT(std::abs(VdotH), 1.0f, ior);  
}

static inline float hydraFresnelDiel2(float3 l, float3 v, float3 n, float ior, float roughness) 
{
  const float cosThetaL = std::abs(dot(l,n)); 
  const float cosThetaV = std::abs(dot(v,n));
  return std::sqrt(FrDielectricPBRT(cosThetaL, 1.0f, ior)*FrDielectricPBRT(cosThetaV, 1.0f, ior));  
}

#endif