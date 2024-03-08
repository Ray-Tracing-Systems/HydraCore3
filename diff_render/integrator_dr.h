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

  float RayTraceDR(uint tid, uint channels, float* out_color, uint a_passNum,
                   const float* a_refImg, const float* a_data, float* a_dataGrad, size_t a_gradSize); ///<! return loss for printing

  float PathTraceDR(uint tid, uint channels, float* out_color, uint a_passNum,
                    const float* a_refImg, const float* a_data, float* a_dataGrad, size_t a_gradSize); ///<! return loss for printing                   
  
  // interaction with diff rendering
  //
  std::pair<size_t, size_t> PutDiffTex2D(uint32_t texId, uint32_t width, uint32_t height, uint32_t channels);
  
  void GetExecutionTime(const char* a_funcName, float a_out[4]) override;

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

  float4 PathTraceReplay(uint tid, uint channels, uint cpuThreadId, float* out_color, 
                         const float* drands, const float* dparams);

  void kernel_InitEyeRay2(uint tid, float4* rayPosAndNear, float4* rayDirAndFar, float4* wavelengths,
                          float4* accumColor, float4* accumuThoroughput, uint* rayFlags, MisData* misData,
                          const float* drands, const float* dparams);

  void kernel_RayTrace2(uint tid, uint bounce, uint cpuThreadId, const float4* rayPosAndNear, const float4* rayDirAndFar,
                        float4* out_hit1, float4* out_hit2, float4* out_hit3, uint* out_instId, uint* rayFlags,
                        const float* dparams);

  void kernel_SampleLightSource(uint tid, uint cpuThreadId, const float4* rayPosAndNear, const float4* rayDirAndFar, const float4* wavelengths, 
                                const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_hitPart3,
                                const uint* rayFlags, uint bounce, float4* out_shadeColor, 
                                const float* drands, const float* dparams);
  
  void kernel_NextBounce(uint tid, uint bounce, const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_hitPart3, const uint* in_instId,
                         const float4* in_shadeColor, float4* rayPosAndNear, float4* rayDirAndFar, const float4* wavelengths,
                         float4* accumColor, float4* accumThoroughput, MisData* a_prevMisData, uint* rayFlags, const float* drands, const float* dparams);

  void kernel_HitEnvironment(uint tid, const uint* rayFlags, const float4* rayDirAndFar, const MisData* a_prevMisData, const float4* accumThoroughput,
                             float4* accumColor, const float* dparams);                              

  void kernel_ContributeToImage(uint tid, const uint* rayFlags, uint channels, const float4* a_accumColor, const uint* in_pakedXY, 
                                const float4* wavelengths, float* out_color, const float* dparams);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  float3 BumpMapping(uint normalMapId, uint currMatId, float3 n, float3 tan, float2 tc, const float* dparams);


  BsdfSample MaterialSampleAndEval(uint a_materialId, uint bounce, float4 wavelengths, float3 v, float3 n, float3 tan, float2 tc, 
                                   MisData* a_misPrev, const uint a_currRayFlags, const float* drands, const float* dparams);
                                    
  BsdfEval   MaterialEval         (uint a_materialId, float4 wavelengths, float3 l, float3 v, float3 n, float3 tan, float2 tc, const float* dparams);

  uint32_t BlendSampleAndEval(uint a_materialId, uint bounce, uint layer, float4 wavelengths, float3 v, float3 n, float2 tc, 
                              MisData* a_misPrev, BsdfSample* a_pRes, const float* drands, const float* dparams);

  MatIdWeightPair BlendEval(MatIdWeight a_mat, float4 wavelengths, float3 l, float3 v, float3 n, float2 tc, const float* dparams);

  uint  RemapMaterialId(uint a_mId, int a_instId, const float* dparams);
  float LightEvalPDF(int a_lightId, float3 illuminationPoint, float3 ray_dir, const float3 lpos, const float3 lnorm, const float* dparams);

  float4 SampleMatColorParamSpectrum(uint32_t matId, float4 a_wavelengths, uint32_t paramId, uint32_t paramSpecId, const float* dparams);
  float4 SampleMatParamSpectrum(uint32_t matId, float4 a_wavelengths, uint32_t paramId, uint32_t paramSpecId, const float* dparams);
  LightSample LightSampleRev(int a_lightId, float2 rands, float3 illiminationPoint, const float* dparams);
  float4 EnvironmentColor(float3 a_dir, const float* dparams);
  float4 LightIntensity(uint a_lightId, const float4* a_wavelengths, float3 a_rayDir, const float* dparams);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  static constexpr int MAX_REC_BLEND  = 2;
  static constexpr int RND_PER_BOUNCE = 2 + 4 + MAX_REC_BLEND;
  static constexpr int RND_LTG_ID     = 0; // float2  lgtRands;
  static constexpr int RND_MTL_ID     = 2; // float4  matRands;
  static constexpr int RND_BLD_ID     = 6; // float   blendRnd[MAX_REC_BLEND];
  static constexpr int LENS_RANDS     = 4; // pixels offsets, wave selector, lens offset in future 

  int MAXTHREADS_CPU   = 1;
  int RANDS_ARRAY_SIZE = 0;

  virtual void SetMaxThreadsAndBounces(int a_maxThreads, int a_maxBounce)
  {
    m_recorded.resize(a_maxThreads);
    for(auto& rec : m_recorded) {
      rec.perBounce.resize(a_maxBounce);
      rec.perBounceLightId.resize(a_maxBounce);
      rec.perBounceRands.resize(RND_PER_BOUNCE*a_maxBounce + LENS_RANDS); // 4 for pixel offsets and waves selector
    }

    MAXTHREADS_CPU   = a_maxThreads;
    m_traceDepth     = a_maxBounce;
    RANDS_ARRAY_SIZE = int(m_recorded[0].perBounceRands.size());
  }
 
  void RecordPixelRndIfNeeded(float2 offsets, float u) override;
  void RecordRayHitIfNeeded(uint32_t bounceId, CRT_Hit hit) override;
  void RecordShadowHitIfNeeded(uint32_t bounceId, bool inShadow) override;
  void RecordLightRndIfNeeded(uint32_t bounceId, int lightId, float2 rands) override;
  void RecordMatRndNeeded(uint32_t bounceId, float4 rands) override;
  void RecordBlendRndNeeded(uint32_t bounceId, uint layer, float rand) override;

  struct PerBounce
  {
    CRT_Hit hit;
    int     inShadow;
  };

  struct PerThreadData
  {
    std::vector<PerBounce> perBounce;
    std::vector<int>       perBounceLightId;
    std::vector<float>     perBounceRands;
  };

  std::vector<PerThreadData>  m_recorded;
  float diffPtTime = 0.0f;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  bool kernel_RayTrace(uint tid, const float4* rayPosAndNear, float4* rayDirAndFar,
                       Lite_Hit* out_hit, float2* out_bars, const float* a_data);
  void kernel_InitEyeRay(uint tid, const uint* packedXY, float4* rayPosAndNear, float4* rayDirAndFar, const float* a_data);
  void kernel_CalcRayColor(uint tid, const Lite_Hit* in_hit, const float2* bars, float4* finalColor, const uint* in_pakedXY, float* out_color, const float* a_data);
};

