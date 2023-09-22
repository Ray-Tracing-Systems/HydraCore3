#pragma once

#include "cglobals.h"

static constexpr uint LIGHT_GEOM_RECT   = 1; 
static constexpr uint LIGHT_GEOM_DISC   = 2;
static constexpr uint LIGHT_GEOM_SPHERE = 3;
static constexpr uint LIGHT_GEOM_DIRECT = 4;

struct LightSource
{
  float4x4 matrix;    ///<! translation in matrix is always (0,0,0,1)
  float4   pos;       ///<! translation aclually stored here
  float4   intensity; ///<! brightress, i.e. screen value if light is visable directly
  float4   norm;
  float2   size;
  float    pdfA;
  uint     geomType;  ///<! LIGHT_GEOM_RECT, LIGHT_GEOM_DISC, LIGHT_GEOM_SPHERE
  //float4 dummy;     ///<! don't forget to manually align stucture to 16 bytes when add new fields
};

struct LightSample
{
  float3 pos;
  float3 norm;
};

static inline LightSample areaLightSampleRev(const LightSource* a_pLight, float2 rands)
{
  float2 sampleOff = 2.0f * (float2(-0.5f,-0.5f) + rands) * a_pLight[0].size;  // PLEASE! use 'a_pLight[0].' for a while ... , not a_pLight-> and not *(a_pLight[0])
  if(a_pLight[0].geomType == LIGHT_GEOM_DISC)
  {
    const float offsetX = rands.x * 2.0f - 1.0f;
    const float offsetY = rands.y * 2.0f - 1.0f;
    sampleOff = MapSamplesToDisc(float2(offsetX, offsetY))*a_pLight[0].size.x; 
  }
  const float3 samplePos = mul3x3(a_pLight[0].matrix, float3(sampleOff.x, 0.0f, sampleOff.y)) + to_float3(a_pLight[0].pos) + epsilonOfPos(to_float3(a_pLight[0].pos)) * to_float3(a_pLight[0].norm);
  LightSample res;
  res.pos  = samplePos;
  res.norm = to_float3(a_pLight[0].norm);
  return res;
}

static inline LightSample sphereLightSampleRev(const LightSource* a_pLight, float2 rands)
{
  const float theta = 2.0f * M_PI * rands.x;
  const float phi   = std::acos(1.0f - 2.0f * rands.y);
  const float x     = std::sin(phi) * std::cos(theta);
  const float y     = std::sin(phi) * std::sin(theta);
  const float z     = std::cos(phi);
  const float3 lcenter   = to_float3(a_pLight[0].pos);
  const float  lradius   = a_pLight[0].size.x;
  const float3 samplePos = lcenter + (lradius*1.000001f)*make_float3(x, y, z);
  LightSample res;
  res.pos  = samplePos;
  res.norm = normalize(samplePos - lcenter);
  return res;
}

static inline LightSample directLightSampleRev(const LightSource* a_pLight, float2 rands, float3 illuminationPoint)
{
  const float3 norm = to_float3(a_pLight[0].norm);
  LightSample res;
  res.pos  = illuminationPoint - norm*100000.0f;
  res.norm = norm;
  return res;
}

//static inline void DirectLightSampleRev(__global const PlainLight* pLight, float3 rands, float3 illuminatingPoint,
//                                        __private ShadowSample* a_out)
//{
//  float3 lpos = lightPos(pLight);
//  float3 norm = lightNorm(pLight);
//
//  float pdfW = 1.0f;
//  if (pLight->data[DIRECT_LIGHT_SSOFTNESS] > 1e-5f)
//  {
//    const float cosAlpha = pLight->data[DIRECT_LIGHT_ALPHA_COS];
//    norm = MapSamplesToCone(cosAlpha, make_float2(rands.x, rands.y), norm);
//    //pdfW = 1.0f / (2.0f * M_PI * (1.0f - cosAlpha));
//  }
//
//  float3 AC    = illuminatingPoint - lpos;
//  float  CBLen = dot(normalize(AC), norm)*length(AC);
//
//  float3 samplePos = illuminatingPoint - norm*CBLen;
//  float hitDist    = CBLen; // length(samplePos - illuminatingPoint);
//
//  const float3 color = lightBaseColor(pLight)*directLightAttenuation(pLight, illuminatingPoint);
//
//  a_out->isPoint    = true; // (pdfW == 1.0f);
//  a_out->pos        = samplePos;
//  a_out->color      = color*pdfW;
//  a_out->pdf        = pdfW;
//  a_out->maxDist    = hitDist;
//  a_out->cosAtLight = 1.0f;
//}