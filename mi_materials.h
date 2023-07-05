#pragma once

#include "include/cmaterial.h"

void SetMiPlastic(GLTFMaterial* material, float int_ior, float ext_ior, float3 diffuse_reflectance, float3 specular_reflectance = float3(1,1,1));
