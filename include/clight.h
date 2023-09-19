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
  float    pdfA;
  uint     geomType;  ///<! LIGHT_GEOM_RECT, LIGHT_GEOM_DISC, LIGHT_GEOM_SPHERE
  //float4 dummy;     ///<! don't forget to manually align stucture to 16 bytes when add new fields
};

static inline float3 areaLightSampleRev(const LightSource* a_pLight, float2 rands)
{
  const float2 sampleOff = 2.0f * (float2(-0.5f,-0.5f) + rands) * a_pLight[0].size;  // PLEASE! use 'a_pLight[0].' for a while ... , not a_pLight-> and not *(a_pLight[0])
  return mul3x3(a_pLight[0].matrix, float3(sampleOff.x, 0.0f, sampleOff.y)) + to_float3(a_pLight[0].pos) + epsilonOfPos(to_float3(a_pLight[0].pos)) * to_float3(a_pLight[0].norm);
}

static inline float3 sphereLightSampleRev(const LightSource* a_pLight, float2 rands)
{
  const float theta = 2.0f * M_PI * rands.x;
  const float phi   = std::acos(1.0f - 2.0f * rands.y);
  const float x     = std::sin(phi) * std::cos(theta);
  const float y     = std::sin(phi) * std::sin(theta);
  const float z     = std::cos(phi);
  const float3 lcenter = to_float3(a_pLight[0].pos);
  const float  lradius = a_pLight[0].size.x;
  return lcenter + (lradius*1.000001f)*make_float3(x, y, z);
}
