#version 460
#extension GL_GOOGLE_include_directive : require
//#extension GL_EXT_ray_query : require
#extension GL_EXT_ray_tracing : require
#extension GL_EXT_nonuniform_qualifier : require

#include "common_generated.h"
layout (constant_id = 1) const int KSPEC_MAT_TYPE_GLTF = 1;
layout (constant_id = 10) const int KSPEC_MAT_TYPE_DIELECTRIC = 10;
layout (constant_id = 11) const int KSPEC_MAT_FOUR_TEXTURES = 11;
layout (constant_id = 12) const int KSPEC_LIGHT_IES = 12;
layout (constant_id = 13) const int KSPEC_LIGHT_ENV = 13;
layout (constant_id = 2) const int KSPEC_MAT_TYPE_GLASS = 2;
layout (constant_id = 3) const int KSPEC_MAT_TYPE_CONDUCTOR = 3;
layout (constant_id = 4) const int KSPEC_MAT_TYPE_DIFFUSE = 4;
layout (constant_id = 5) const int KSPEC_MAT_TYPE_PLASTIC = 5;
layout (constant_id = 6) const int KSPEC_SPECTRAL_RENDERING = 6;
layout (constant_id = 7) const int KSPEC_MAT_TYPE_BLEND = 7;
layout (constant_id = 8) const int KSPEC_BLEND_STACK_SIZE = 8;
layout (constant_id = 9) const int KSPEC_BUMP_MAPPING = 9;

layout(binding = 0, set = 0) buffer data0 { float out_color[]; }; //
layout(binding = 1, set = 0) buffer data1 { vec4 m_vNorm4f[]; }; //
layout(binding = 2, set = 0) buffer data2 { vec4 m_vTang4f[]; }; //
layout(binding = 3, set = 0) buffer data3 { uint m_triIndices[]; }; //
layout(binding = 4, set = 0) buffer data4 { int m_allRemapListsOffsets[]; }; //
layout(binding = 5, set = 0) buffer data5 { int m_allRemapLists[]; }; //
layout(binding = 6, set = 0) buffer data6 { float m_precomp_coat_transmittance[]; }; //
layout(binding = 7, set = 0) buffer data7 { LightSource m_lights[]; }; //
layout(binding = 8, set = 0) buffer data8 { Material m_materials[]; }; //
layout(binding = 9, set = 0) buffer data9 { int m_remapInst[]; }; //
layout(binding = 10, set = 0) buffer data10 { uint m_vertOffset[]; }; //
layout(binding = 11, set = 0) buffer data11 { float m_spec_values[]; }; //
layout(binding = 12, set = 0) buffer data12 { uint m_packedXY[]; }; //
layout(binding = 13, set = 0) uniform sampler2D m_textures[]; //
layout(binding = 14, set = 0) buffer data14 { float m_pdfLightData[]; }; //
layout(binding = 15, set = 0) buffer data15 { uvec2 m_spec_offset_sz[]; }; //
layout(binding = 16, set = 0) uniform accelerationStructureEXT m_pAccelStruct;
layout(binding = 17, set = 0) buffer data17 { uint m_matIdByPrimId[]; }; //
layout(binding = 18, set = 0) buffer data18 { uint m_matIdOffsets[]; }; //
layout(binding = 19, set = 0) buffer dataUBO { Integrator_Generated_UBO_Data ubo; };

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
layout(location = 0) rayPayloadEXT CRT_Hit kgen_hitValue;
layout(location = 1) rayPayloadEXT bool    kgen_inShadow;

// RayScene intersection with 'm_pAccelStruct'
//
CRT_Hit m_pAccelStruct_RayQuery_NearestHit(const vec4 rayPos, const vec4 rayDir)
{
	traceRayEXT(m_pAccelStruct, gl_RayFlagsOpaqueEXT, 0xff, 0, 0, 0, rayPos.xyz, rayPos.w, rayDir.xyz, rayDir.w, 0);
  return kgen_hitValue;
  
  //rayQueryEXT rayQuery;
  //rayQueryInitializeEXT(rayQuery, m_pAccelStruct, gl_RayFlagsOpaqueEXT, 0xff, rayPos.xyz, rayPos.w, rayDir.xyz, rayDir.w);
  //
  //while(rayQueryProceedEXT(rayQuery)) { } // actually may omit 'while' when 'gl_RayFlagsOpaqueEXT' is used
  //
  //CRT_Hit res;
  //res.primId = -1;
  //res.instId = -1;
  //res.geomId = -1;
  //res.t      = rayDir.w;
  //
  //if(rayQueryGetIntersectionTypeEXT(rayQuery, true) == gl_RayQueryCommittedIntersectionTriangleEXT)
  //{    
	//  res.primId    = rayQueryGetIntersectionPrimitiveIndexEXT(rayQuery, true);
	//  res.geomId    = rayQueryGetIntersectionInstanceCustomIndexEXT(rayQuery, true);
  //  res.instId    = rayQueryGetIntersectionInstanceIdEXT    (rayQuery, true);
	//  res.t         = rayQueryGetIntersectionTEXT(rayQuery, true);
  //  vec2 bars     = rayQueryGetIntersectionBarycentricsEXT(rayQuery, true);
  //  
  //  res.coords[0] = bars.y;
  //  res.coords[1] = bars.x;
  //  res.coords[2] = 1.0f - bars.y - bars.x;
  //}
  //return res;
}

bool m_pAccelStruct_RayQuery_AnyHit(const vec4 rayPos, const vec4 rayDir)
{
  traceRayEXT(m_pAccelStruct, gl_RayFlagsOpaqueEXT | gl_RayFlagsTerminateOnFirstHitEXT | gl_RayFlagsSkipClosestHitShaderEXT,
              0xff, 0, 0, 1, rayPos.xyz, rayPos.w, rayDir.xyz, rayDir.w, 1);
  return kgen_inShadow;

  //rayQueryEXT rayQuery;
  //rayQueryInitializeEXT(rayQuery, m_pAccelStruct, gl_RayFlagsTerminateOnFirstHitEXT, 0xff, rayPos.xyz, rayPos.w, rayDir.xyz, rayDir.w);
  //rayQueryProceedEXT(rayQuery);
  //return (rayQueryGetIntersectionTypeEXT(rayQuery, true) == gl_RayQueryCommittedIntersectionTriangleEXT);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

layout( push_constant ) uniform kernelArgs
{
  uint iNumElementsX; 
  uint iNumElementsY; 
  uint iNumElementsZ; 
  uint tFlagsMask;    
} kgenArgs;

///////////////////////////////////////////////////////////////// subkernels here
void kernel_GetRayColor_m_packedXY_out_color(uint tid, in Lite_Hit in_hit, in vec2 bars, uint in_pakedXYOffset, uint out_colorOffset) 
{
   
  if(tid >= ubo.m_maxThreadId)
    return;

  const Lite_Hit hit = in_hit;
  if(hit.geomId == -1)
  {
    out_color[tid + out_colorOffset] = 0;
    return;
  }

  const uint matId = m_matIdByPrimId[m_matIdOffsets[hit.geomId] + uint(hit.primId)];
  const vec4 mdata = m_materials[matId].colors[GLTF_COLOR_BASE];
  const vec2 uv = bars;

  const uint triOffset  = m_matIdOffsets[hit.geomId];
  const uint vertOffset = m_vertOffset  [hit.geomId];

  const uint A = m_triIndices[(triOffset + uint(hit.primId))*3 + 0];
  const uint B = m_triIndices[(triOffset + uint(hit.primId))*3 + 1];
  const uint C = m_triIndices[(triOffset + uint(hit.primId))*3 + 2];
  const vec4 data1 = (1.0f - uv.x - uv.y)*m_vNorm4f[A + vertOffset] + uv.y*m_vNorm4f[B + vertOffset] + uv.x*m_vNorm4f[C + vertOffset];
  const vec4 data2 = (1.0f - uv.x - uv.y)*m_vTang4f[A + vertOffset] + uv.y*m_vTang4f[B + vertOffset] + uv.x*m_vTang4f[C + vertOffset];
  //float3 hitNorm     = to_float3(data1);
  //float3 hitTang     = to_float3(data2);
  vec2 hitTexCoord = vec2(data1.w,data2.w);

  const uint   texId     = m_materials[matId].texid[0];
  const vec2 texCoordT = mulRows2x4(m_materials[matId].row0[0], m_materials[matId].row1[0], hitTexCoord);
  const vec4 texColor = texture(m_textures[texId], texCoordT);

  const vec3 color = mdata.w > 0.0f ? clamp(vec3(mdata.w,mdata.w,mdata.w), 0.0f, 1.0f) : (mdata*texColor).xyz;

  const uint XY = m_packedXY[tid + in_pakedXYOffset];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;

  out_color[(y*uint(ubo.m_winWidth)+x)*4 + 0 + out_colorOffset] = color.x;
  out_color[(y*uint(ubo.m_winWidth)+x)*4 + 1 + out_colorOffset] = color.y;
  out_color[(y*uint(ubo.m_winWidth)+x)*4 + 2 + out_colorOffset] = color.z;
  out_color[(y*uint(ubo.m_winWidth)+x)*4 + 3 + out_colorOffset] = 0.0f;

}

void kernel_InitEyeRay_m_packedXY(uint tid, uint packedXYOffset, inout vec4 rayPosAndNear, inout vec4 rayDirAndFar) 
{
  
  if(tid >= ubo.m_maxThreadId)
    return;
  const uint XY = m_packedXY[tid + packedXYOffset];

  const uint x = (XY & 0x0000FFFF);
  const uint y = (XY & 0xFFFF0000) >> 16;

  vec3 rayDir = EyeRayDirNormalized((float(x)+0.5f)/float(ubo.m_winWidth), (float(y)+0.5f)/float(ubo.m_winHeight), ubo.m_projInv);
  vec3 rayPos = vec3(0,0,0);

  transform_ray3f(ubo.m_worldViewInv, 
                  rayPos, rayDir);
  
  rayPosAndNear = vec4(rayPos, 0.0f);
  rayDirAndFar  = vec4(rayDir, FLT_MAX);

}

bool kernel_RayTrace(uint tid, in vec4 rayPosAndNear, inout vec4 rayDirAndFar, inout Lite_Hit out_hit, inout vec2 out_bars) 
{
  
  if(tid >= ubo.m_maxThreadId)
    return false;
  const vec4 rayPos = rayPosAndNear;
  const vec4 rayDir = rayDirAndFar ;

  CRT_Hit hit = m_pAccelStruct_RayQuery_NearestHit(rayPos, rayDir);
  
  Lite_Hit res;
  res.primId = int(hit.primId);
  res.instId = int(hit.instId);
  res.geomId = int(hit.geomId);
  res.t      = hit.t;

  vec2 baricentrics = vec2(hit.coords[0],hit.coords[1]);
 
  out_hit  = res;
  out_bars = baricentrics;
  return (res.primId != -1);

}

///////////////////////////////////////////////////////////////// subkernels here

void main()
{
  ///////////////////////////////////////////////////////////////// prolog
  //const uint tid = uint(gl_GlobalInvocationID[0]); 
  const uint tid = uint(gl_LaunchIDEXT.x); 
  ///////////////////////////////////////////////////////////////// prolog

  
  vec4 rayPosAndNear,  rayDirAndFar;
  kernel_InitEyeRay_m_packedXY(tid, 0, rayPosAndNear, rayDirAndFar);

  Lite_Hit hit; 
  vec2 baricentrics; 
  if(!kernel_RayTrace(tid, rayPosAndNear, rayDirAndFar, hit, baricentrics))
    return;
  
  kernel_GetRayColor_m_packedXY_out_color(tid, hit, baricentrics, 0, 0);

}

