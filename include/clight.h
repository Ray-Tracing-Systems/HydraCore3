#pragma once

#include "cglobals.h"

static constexpr uint LIGHT_GEOM_RECT   = 1; 
static constexpr uint LIGHT_GEOM_DISC   = 2;
static constexpr uint LIGHT_GEOM_SPHERE = 3;

struct LightSource
{
  float4x4 matrix;    ///<! translation in matrix is always (0,0,0,1)
  float4   pos;       ///<! translation aclually stored here
  float4   intensity; ///<! brightress, i.e. screen value if light is visable directly
  float4   norm;
  float2   size;
  uint     geomType;  ///<! LIGHT_GEOM_RECT, LIGHT_GEOM_DISC, LIGHT_GEOM_SPHERE
  float    dummy;
};

struct LightSample
{
  float3 pos;
  float  pdfA;
};

static inline LightSample areaLightSampleRev(const LightSource* a_pLight, float2 rands)
{
  const float2 sampleOff = 2.0f * (float2(-0.5f,-0.5f) + rands) * a_pLight[0].size;
  float3 samplePos = float3(sampleOff.x, 0.0f, sampleOff.y); // -1e-5f*std::max(m_light.size.x, m_light.size.y)

  samplePos = mul3x3(a_pLight[0].matrix, samplePos) +
              epsilonOfPos(to_float3(a_pLight[0].pos)) * to_float3(a_pLight[0].norm) +
              to_float3(a_pLight[0].pos);
  
  LightSample res;
  res.pos  = samplePos;
  res.pdfA = 1.0f / (4.0f * a_pLight[0].size.x * a_pLight[0].size.y);
  return res;
}


static inline LightSample sphereLightSampleRev(const LightSource* a_pLight, float2 rands)
{
  LightSample res;
  return res;
}



//hydra2:
//static inline void SphereLightSampleRev(__global const PlainLight* pLight, float3 rands, float3 illuminatingPoint,
//                                        __private ShadowSample* a_out)
//{
//  const float theta = 2.0f * M_PI * rands.x;
//  const float phi   = acos(1.0f - 2.0f * rands.y);
//  const float x     = sin(phi) * cos(theta);
//  const float y     = sin(phi) * sin(theta);
//  const float z     = cos(phi);
//
//  const float3 lcenter = lightPos(pLight);
//  const float  lradius = pLight->data[SPHERE_LIGHT_RADIUS];
//
//  const float3 samplePos = lcenter + lradius*make_float3(x, y, z);
//  const float3 lightNorm = normalize(samplePos - lcenter);
//  const float3 dirToV    = normalize(samplePos - illuminatingPoint);
//
//  a_out->isPoint    = false;
//  a_out->pos        = samplePos;
//  a_out->color      = sphereLightGetIntensity(pLight);
//  a_out->pdf        = sphereLightEvalPDF(pLight, illuminatingPoint, samplePos, lightNorm);
//  a_out->maxDist    = length(samplePos - illuminatingPoint);
//  a_out->cosAtLight = fabs(dot(lightNorm, dirToV));
//}

//hydra2:
//static inline float sphereLightEvalPDF(__global const PlainLight* pLight, float3 illuminatingPoint, float3 lpos, float3 lnorm)
//{
//  float  lradius = pLight->data[SPHERE_LIGHT_RADIUS];
//  float3 lcenter = lightPos(pLight);
//
//  if (DistanceSquared(illuminatingPoint, lcenter) - lradius*lradius <= 0.0f)
//    return 1.0f; // 
//
//  const float  pdfA   = 1.0f / pLight->data[PLIGHT_SURFACE_AREA];
//  const float  dist   = length(lpos - illuminatingPoint);
//
//  const float3 dirToV    = normalize(lpos - illuminatingPoint);
//  const float cosAtLight = fabs(dot(dirToV, lnorm));
//
//  return PdfAtoW(pdfA, dist, cosAtLight);
//}
