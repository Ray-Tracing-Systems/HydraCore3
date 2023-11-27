#pragma once

#include "include/cmaterial.h"

void SetMiPlastic(Material* material, float int_ior, float ext_ior, float4 diffuse_reflectance, float4 specular_reflectance = float4(1));

namespace mi
{
  float fresnel_diffuse_reflectance(float eta);
}