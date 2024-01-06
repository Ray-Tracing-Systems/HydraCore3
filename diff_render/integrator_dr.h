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

class IntegratorDR : public Integrator
{
public:

  IntegratorDR(int a_maxThreads = 1, int a_spectral_mode = 0, int a_gradMode = 1, std::vector<uint32_t> a_features = {}) : Integrator(a_maxThreads, a_spectral_mode, a_features), m_gradMode(a_gradMode) {}
  virtual ~IntegratorDR() {}

  void LoadSceneEnd() override;

  float4 CastRayDR(uint tid, uint channels, float* out_color, const float* a_data);
  
  bool kernel_RayTrace(uint tid, const float4* rayPosAndNear, float4* rayDirAndFar,
                       Lite_Hit* out_hit, float2* out_bars, const float* a_data);
  void kernel_InitEyeRay(uint tid, const uint* packedXY, float4* rayPosAndNear, float4* rayDirAndFar, const float* a_data);

  void kernel_CalcRayColor(uint tid, const Lite_Hit* in_hit, const float2* bars, float4* finalColor, const uint* in_pakedXY, float* out_color, const float* a_data);

  float RayTraceDR(uint tid, uint channels, float* out_color, uint a_passNum,
                   const float* a_refImg, const float* a_data, float* a_dataGrad, size_t a_gradSize); ///<! return loss for printing

  float PathTraceDR(uint tid, uint channels, float* out_color, uint a_passNum,
                    const float* a_refImg, const float* a_data, float* a_dataGrad, size_t a_gradSize); ///<! return loss for printing                   
  
  // interaction with diff rendering
  //
  std::pair<size_t, size_t> PutDiffTex2D(uint32_t texId, uint32_t width, uint32_t height, uint32_t channels);

protected:

  float4 Tex2DFetchAD(uint texId, float2 texCoord, const float* tex_data);

  struct TexInfo
  {
    size_t   offset;
    int32_t  width;
    int32_t  height;
    int32_t  channels;
    float    fwidth;
    float    fheight;
  };

  std::vector<TexInfo> m_texAddressTable;
  int m_gradMode;
  size_t m_gradSize = 0;

public:

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  float4 PathTrace(uint tid, uint channels, float* out_color, const float* dparams);

  void kernel_InitEyeRay2(uint tid, const uint* packedXY, float4* rayPosAndNear, float4* rayDirAndFar, float4* wavelengths,
                          float4* accumColor, float4* accumuThoroughput, RandomGen* gen, uint* rayFlags, MisData* misData,
                          const float* dparams);

  void kernel_RayTrace2(uint tid, const float4* rayPosAndNear, const float4* rayDirAndFar,
                        float4* out_hit1, float4* out_hit2, float4* out_hit3, uint* out_instId, uint* rayFlags,
                        const float* dparams);

  void kernel_SampleLightSource(uint tid, const float4* rayPosAndNear, const float4* rayDirAndFar, const float4* wavelengths, 
                                const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_hitPart3,
                                const uint* rayFlags, uint bounce,
                                RandomGen* a_gen, float4* out_shadeColor, const float* dparams);
  
  void kernel_NextBounce(uint tid, uint bounce, const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_hitPart3, const uint* in_instId,
                         const float4* in_shadeColor, float4* rayPosAndNear, float4* rayDirAndFar, const float4* wavelengths,
                         float4* accumColor, float4* accumThoroughput, RandomGen* a_gen, MisData* a_prevMisData, uint* rayFlags, const float* dparams);

  void kernel_HitEnvironment(uint tid, const uint* rayFlags, const float4* rayDirAndFar, const MisData* a_prevMisData, const float4* accumThoroughput,
                             float4* accumColor, const float* dparams);                              

  void kernel_ContributeToImage(uint tid, uint channels, const float4* a_accumColor, const RandomGen* gen, const uint* in_pakedXY, 
                                const float4* wavelengths, float* out_color, const float* dparams);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  float3 BumpMapping(uint normalMapId, uint currMatId, float3 n, float3 tan, float2 tc, const float* dparams);


  BsdfSample MaterialSampleAndEval(uint a_materialId, float4 wavelengths, RandomGen* a_gen, float3 v, float3 n, float3 tan, float2 tc, 
                                   MisData* a_misPrev, const uint a_currRayFlags, const float* dparams);
                                    
  BsdfEval   MaterialEval         (uint a_materialId, float4 wavelengths, float3 l, float3 v, float3 n, float3 tan, float2 tc, const float* dparams);

  uint32_t BlendSampleAndEval(uint a_materialId, float4 wavelengths, RandomGen* a_gen, float3 v, float3 n, float2 tc, 
                              MisData* a_misPrev, BsdfSample* a_pRes, const float* dparams);

  MatIdWeightPair BlendEval(MatIdWeight a_mat, float4 wavelengths, float3 l, float3 v, float3 n, float2 tc, const float* dparams);

  uint  RemapMaterialId(uint a_mId, int a_instId, const float* dparams);
  float LightEvalPDF(int a_lightId, float3 illuminationPoint, float3 ray_dir, const float3 lpos, const float3 lnorm, const float* dparams);

  float4 SampleMatColorParamSpectrum(uint32_t matId, float4 a_wavelengths, uint32_t paramId, uint32_t paramSpecId, const float* dparams);
  float4 SampleMatParamSpectrum(uint32_t matId, float4 a_wavelengths, uint32_t paramId, uint32_t paramSpecId, const float* dparams);
  LightSample LightSampleRev(int a_lightId, float2 rands, float3 illiminationPoint, const float* dparams);
  float4 GetEnvironmentColorAndPdf(float3 a_dir, const float* dparams);
  float4 GetLightSourceIntensity(uint a_lightId, const float4* a_wavelengths, const float* dparams);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

};

