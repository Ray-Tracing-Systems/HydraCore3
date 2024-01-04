#pragma once

#include "include/cglobals.h" // We assume that all code that should pe passed to kernels will be just included both for CPU and OpenCL
#include "include/crandom.h"
#include "include/clight.h"
#include "include/cmaterial.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <cfloat>

#include "CrossRT.h" // special include for ray tracing
#include "Image2d.h" // special include for textures

#include "spectrum.h"
#include "integrator_pt.h"

using LiteImage::ICombinedImageSampler;


struct MaterialNonDiff
{
  uint lambertTexId;
};

//struct MaterialDiff
//{
//  float4 color;
//};


class IntegratorDR : public Integrator
{
public:

  IntegratorDR(int a_maxThreads = 1, int a_spectral_mode = 0, std::vector<uint32_t> a_features = {}) : Integrator(a_maxThreads, a_spectral_mode, a_features){}
  virtual ~IntegratorDR() {}

  void LoadSceneEnd() override;

  float4 CastRayDR(uint tid, uint channels, float* out_color, const float* a_data);
  
  bool kernel_RayTrace(uint tid, const float4* rayPosAndNear, float4* rayDirAndFar,
                       Lite_Hit* out_hit, float2* out_bars, const float* a_data);
  void kernel_InitEyeRay(uint tid, const uint* packedXY, float4* rayPosAndNear, float4* rayDirAndFar, const float* a_data);

  void kernel_CalcRayColor(uint tid, const Lite_Hit* in_hit, const float2* bars, float4* finalColor, const uint* in_pakedXY, float* out_color, const float* a_data);

  void PathTraceDR(uint tid, uint channels, float* out_color, uint a_passNum,
                   const float* a_refImg, const float* a_data, float* a_dataGrad, size_t a_gradSize);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //float4 PathTrace(uint tid, uint channels, float* out_color, const float* a_texData);
  //
  //void kernel_SampleLightSource(uint tid, const float4* rayPosAndNear, const float4* rayDirAndFar, 
  //                              const float4* wavelengths, const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_hitPart3,
  //                              const uint* rayFlags, uint bounce,
  //                              RandomGen* a_gen, float4* out_shadeColor,  const float* a_data);
  //
  //void kernel_NextBounce(uint tid, uint bounce, const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_hitPart3, const uint* in_instId,
  //                       const float4* in_shadeColor, float4* rayPosAndNear, float4* rayDirAndFar, const float4* wavelengths,
  //                       float4* accumColor, float4* accumThoroughput, RandomGen* a_gen, MisData* misPrev, uint* rayFlags, const float* a_data);
  //
  //BsdfEval MaterialEval(uint a_materialId, float4 wavelengths, float3 l, float3 v, float3 n, float3 tan, float2 tc, const float* a_data);
  //BsdfSample MaterialSampleAndEval(uint a_materialId, float4 wavelengths, RandomGen* a_gen, float3 v, float3 n, float3 tan, float2 tc, 
  //                                 MisData* a_misPrev, const uint a_currRayFlags, const float* a_data);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//protected:
  float4 HydraTex2DFetch(uint texId, float2 texCoord, const float* tex_data);
  
  std::vector<MaterialNonDiff> m_matNonDiff;
  //std::vector<MaterialDiff>    m_matDiff;
};

