#ifndef Integrator_UBO_H
#define Integrator_UBO_H

#ifndef GLSL
#define LAYOUT_STD140
#include "LiteMath.h"
using LiteMath::uint;
typedef LiteMath::float4x4 mat4;
typedef LiteMath::float2   vec2;
typedef LiteMath::float3   vec3;
typedef LiteMath::float4   vec4;
typedef LiteMath::int2     ivec2;
typedef LiteMath::int3     ivec3;
typedef LiteMath::int4     ivec4;
typedef LiteMath::uint2    uvec2;
typedef LiteMath::uint3    uvec3;
typedef LiteMath::uint4    uvec4;
#else
#define M_PI          3.14159265358979323846f
#define M_TWOPI       6.28318530717958647692f
#define INV_PI        0.31830988618379067154f
#define INV_TWOPI     0.15915494309189533577f
#endif

struct Integrator_Generated_UBO_Data
{
  mat4 m_projInv; 
  mat4 m_worldViewInv; 
  vec4 m_camRespoceRGB; 
  vec4 m_envColor; 
  int m_camResponseSpectrumId[3]; 
  int m_camResponseType; 
  uint m_disableImageContrib; 
  float m_exposureMult; 
  uint m_intergatorType; 
  uint m_maxThreadId; 
  uint m_skipBounce; 
  int m_spectral_mode; 
  uint m_tileSize; 
  uint m_traceDepth; 
  int m_winHeight; 
  int m_winWidth; 
  uint m_allRemapListsOffsets_capacity; 
  uint m_allRemapListsOffsets_size; 
  uint m_allRemapLists_capacity; 
  uint m_allRemapLists_size; 
  uint m_cie_x_capacity; 
  uint m_cie_x_size; 
  uint m_cie_y_capacity; 
  uint m_cie_y_size; 
  uint m_cie_z_capacity; 
  uint m_cie_z_size; 
  uint m_instIdToLightInstId_capacity; 
  uint m_instIdToLightInstId_size; 
  uint m_lights_capacity; 
  uint m_lights_size; 
  uint m_matIdByPrimId_capacity; 
  uint m_matIdByPrimId_size; 
  uint m_matIdOffsets_capacity; 
  uint m_matIdOffsets_size; 
  uint m_materials_capacity; 
  uint m_materials_size; 
  uint m_normMatrices_capacity; 
  uint m_normMatrices_size; 
  uint m_pAccelStruct_capacity; 
  uint m_pAccelStruct_size; 
  uint m_packedXY_capacity; 
  uint m_packedXY_size; 
  uint m_precomp_coat_transmittance_capacity; 
  uint m_precomp_coat_transmittance_size; 
  uint m_randomGens_capacity; 
  uint m_randomGens_size; 
  uint m_remapInst_capacity; 
  uint m_remapInst_size; 
  uint m_spec_offset_sz_capacity; 
  uint m_spec_offset_sz_size; 
  uint m_spec_values_capacity; 
  uint m_spec_values_size; 
  uint m_textures_capacity; 
  uint m_textures_size; 
  uint m_triIndices_capacity; 
  uint m_triIndices_size; 
  uint m_vNorm4f_capacity; 
  uint m_vNorm4f_size; 
  uint m_vTang4f_capacity; 
  uint m_vTang4f_size; 
  uint m_vertOffset_capacity; 
  uint m_vertOffset_size; 
  uint dummy_last;
};

#endif

