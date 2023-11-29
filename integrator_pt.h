#ifndef TEST_CLASS_H
#define TEST_CLASS_H

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
#include "cam_plugin/CamPluginAPI.h"

using LiteImage::ICombinedImageSampler;


class Integrator // : public DataClass, IRenderer
{
public:

  Integrator(int a_maxThreads = 1, int a_spectral_mode = 0, std::vector<uint32_t> a_features = {}) : m_spectral_mode(a_spectral_mode), m_enabledFeatures(a_features)
  {
    InitRandomGens(a_maxThreads);
    m_pAccelStruct = std::shared_ptr<ISceneObject>(CreateSceneRT(""), [](ISceneObject *p) { DeleteSceneRT(p); } );
  }

  virtual ~Integrator() { m_pAccelStruct = nullptr; }

  virtual void SceneRestrictions(uint32_t a_restrictions[4]) const
  {
    uint32_t maxMeshes            = 4096;
    uint32_t maxTotalVertices     = 16'000'000;
    uint32_t maxTotalPrimitives   = 16'000'000;
    uint32_t maxPrimitivesPerMesh = 4'000'000;

    a_restrictions[0] = maxMeshes;
    a_restrictions[1] = maxTotalVertices;
    a_restrictions[2] = maxTotalPrimitives;
    a_restrictions[3] = maxPrimitivesPerMesh;
  }

  void InitRandomGens(int a_maxThreads);

  void SetAccelStruct(std::shared_ptr<ISceneObject> a_customAccelStruct) { m_pAccelStruct = a_customAccelStruct; };
  virtual bool LoadScene(const char* a_scehePath, const char* a_sncDir);

  void PackXY         (uint tidX, uint tidY);

#ifdef KERNEL_SLICER
  void CastSingleRay  (uint tid, uint* out_color   __attribute__((size("tid"))) ); ///<! ray casting, draw diffuse or emisive color
  void RayTrace       (uint tid, float4* out_color __attribute__((size("tid"))) ); ///<! whitted ray tracing

  void NaivePathTrace (uint tid, float4* out_color __attribute__((size("tid"))) ); ///<! NaivePT
  void PathTrace      (uint tid, float4* out_color __attribute__((size("tid"))) ); ///<! MISPT and ShadowPT
  void PathTraceFromInputRays(uint tid, const RayPart1*  in_rayPosAndNear __attribute__((size("tid"))), 
                                        const RayPart2*  in_rayDirAndFar  __attribute__((size("tid"))),
                                        float4*        out_color        __attribute__((size("tid"))));
#else
  void CastSingleRay  (uint tid, uint* out_color); ///<! ray casting, draw diffuse or emisive color
  void RayTrace       (uint tid, float4* out_color); ///<! whitted ray tracing

  void NaivePathTrace (uint tid, float4* out_color); ///<! NaivePT
  void PathTrace      (uint tid, float4* out_color); ///<! MISPT and ShadowPT
  void PathTraceFromInputRays(uint tid, const RayPart1* in_rayPosAndNear, const RayPart2* in_rayDirAndFar,
                              float4* out_color);
#endif

  virtual void PackXYBlock(uint tidX, uint tidY, uint a_passNum);
  virtual void CastSingleRayBlock(uint tid, uint* out_color, uint a_passNum);
  virtual void NaivePathTraceBlock(uint tid, float4* out_color, uint a_passNum);
  virtual void PathTraceBlock(uint tid, float4* out_color, uint a_passNum);
  virtual void PathTraceFromInputRaysBlock(uint tid, const RayPart1* in_rayPosAndNear, const RayPart2* in_rayDirAndFar, 
                                           float4* out_color, uint a_passNum);
  virtual void RayTraceBlock(uint tid, float4* out_color, uint a_passNum);

  virtual void CommitDeviceData() {}                                     // will be overriden in generated class
  virtual void GetExecutionTime(const char* a_funcName, float a_out[4]); // will be overriden in generated class

  virtual void UpdateMembersPlainData() {}                               // will be overriden in generated class, optional function
  //virtual void UpdateMembersVectorData() {}                              // will be overriden in generated class, optional function
  //virtual void UpdateMembersTexureData() {}                              // will be overriden in generated class, optional function

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void kernel_PackXY(uint tidX, uint tidY, uint* out_pakedXY);

  void kernel_InitEyeRay(uint tid, const uint* packedXY, float4* rayPosAndNear, float4* rayDirAndFar);        // (tid,tidX,tidY,tidZ) are SPECIAL PREDEFINED NAMES!!!
  void kernel_InitEyeRay2(uint tid, const uint* packedXY, float4* rayPosAndNear, float4* rayDirAndFar, float4* wavelengths,
                          float4* accumColor, float4* accumuThoroughput, RandomGen* gen, uint* rayFlags, MisData* misData);
  void kernel_InitEyeRay3(uint tid, const uint* packedXY, float4* rayPosAndNear, float4* rayDirAndFar, float4* accumColor,
                          float4* accumuThoroughput, uint* rayFlags);        

  void kernel_InitEyeRayFromInput(uint tid, const RayPart1* in_rayPosAndNear, const RayPart2* in_rayDirAndFar,
                                  float4* rayPosAndNear, float4* rayDirAndFar, float4* accumColor, float4* accumuThoroughput, 
                                  RandomGen* gen, uint* rayFlags, MisData* misData, float4* wavelengths);

  bool kernel_RayTrace(uint tid, const float4* rayPosAndNear, float4* rayDirAndFar,
                       Lite_Hit* out_hit, float2* out_bars);

  void kernel_RayTrace2(uint tid, const float4* rayPosAndNear, const float4* rayDirAndFar,
                        float4* out_hit1, float4* out_hit2, uint* out_instId, uint* rayFlags);

  void kernel_GetRayColor(uint tid, const Lite_Hit* in_hit, const uint* in_pakedXY, uint* out_color);

  void kernel_NextBounce(uint tid, uint bounce, const float4* in_hitPart1, const float4* in_hitPart2, const uint* in_instId,
                         const float4* in_shadeColor, float4* rayPosAndNear, float4* rayDirAndFar, const float4* wavelengths,
                         float4* accumColor, float4* accumThoroughput, RandomGen* a_gen, MisData* a_prevMisData, uint* rayFlags);

  void kernel_RayBounce(uint tid, uint bounce, const float4* in_hitPart1, const float4* in_hitPart2,
                        float4* rayPosAndNear, float4* rayDirAndFar, float4* accumColor, float4* accumThoroughput, uint* rayFlags);

  void kernel_SampleLightSource(uint tid, const float4* rayPosAndNear, const float4* rayDirAndFar, const float4* wavelengths, 
                                const float4* in_hitPart1, const float4* in_hitPart2, 
                                const uint* rayFlags, uint bounce,
                                RandomGen* a_gen, float4* out_shadeColor);

  void kernel_HitEnvironment(uint tid, const uint* rayFlags, const float4* rayDirAndFar, const MisData* a_prevMisData, const float4* accumThoroughput,
                             float4* accumColor);

  void kernel_RealColorToUint32(uint tid, float4* a_accumColor, uint* out_color);

  void kernel_ContributeToImage(uint tid, const float4* a_accumColor, const RandomGen* gen, const uint* in_pakedXY, 
                                const float4* wavelengths, float4* out_color);

  void kernel_CopyColorToOutput(uint tid, const float4* a_accumColor, const RandomGen* gen, 
                                float4* out_color);

  void kernel_ContributeToImage3(uint tid, const float4* a_accumColor, const uint* in_pakedXY, float4* out_color);                               
  void kernel_ContributePathRayToImage3(float4* out_color, const std::vector<float4>& a_rayColor, std::vector<float3>& a_rayPos);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  static constexpr uint INTEGRATOR_STUPID_PT = 0;
  static constexpr uint INTEGRATOR_SHADOW_PT = 1;
  static constexpr uint INTEGRATOR_MIS_PT    = 2;

  static inline bool isDeadRay     (uint a_flags)  { return (a_flags & RAY_FLAG_IS_DEAD)        != 0; }
  static inline bool hasNonSpecular(uint a_flags)  { return (a_flags & RAY_FLAG_HAS_NON_SPEC)   != 0; }
  static inline bool hasInvNormal  (uint a_flags)  { return (a_flags & RAY_FLAG_HAS_INV_NORMAL) != 0; }
  static inline bool isOutOfScene  (uint a_flags)  { return (a_flags & RAY_FLAG_OUT_OF_SCENE)   != 0; }

  static inline uint extractMatId(uint a_flags)    { return (a_flags & 0x00FFFFFF); }       
  static inline uint packMatId(uint a_flags, uint a_matId) { return (a_flags & 0xFF000000) | (a_matId & 0x00FFFFFF); }       
  static inline uint maxMaterials()             { return 0x00FFFFFF+1; }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////// CPU API

  void SetIntegratorType(const uint a_type) { m_intergatorType = a_type; }
  void SetViewport(int a_xStart, int a_yStart, int a_width, int a_height) 
  { 
    m_winStartX = a_xStart; 
    m_winStartY = a_yStart;
    m_winWidth  = a_width;  // todo: remember a_width for first call as pitch and dont change pitch anymore?
    m_winHeight = a_height;
    m_packedXY.resize(m_winWidth*m_winHeight); // todo: use a_xStart,a_yStart

    const auto sizeX = a_width  - a_xStart;
    const auto sizeY = a_height - a_yStart;

    if(sizeX % 8 == 0 && sizeY % 8 == 0)
      m_tileSize = 8;
    else if(sizeX % 4 == 0 && sizeY % 4 == 0)
      m_tileSize = 4;
    else if(sizeX % 2 == 0 && sizeY % 2 == 0)
      m_tileSize = 2;
    else 
      m_tileSize = 1;
    
    m_maxThreadId = sizeX*sizeY;
  }
  
  void SetWorldView(const float4x4& a_mat)
  {
    m_worldView = a_mat;
    m_worldViewInv = LiteMath::inverse4x4(m_worldView);
  }

  uint GetSPP() const { return m_spp; } 

protected:
  int m_winStartX   = 0;
  int m_winStartY   = 0;
  int m_winWidth    = 512;
  int m_winHeight   = 512;
  uint m_traceDepth = 10;
  uint m_skipBounce = 0; ///!< when greater than 1, skip all bounce before this one: 2 for secondary light, 3 for thertiary and e.t.c. 
                         ///!< TODO: don't account specular bounces(!)
  uint m_spp         = 1024;
  uint m_tileSize    = 8; ///!< screen mini tile, 2x2, 4x4 or 8x8 pixels.
  uint m_maxThreadId = m_winWidth*m_winHeight;

  LightSample LightSampleRev(int a_lightId, float2 rands, float3 illiminationPoint);
  float LightPdfSelectRev(int a_lightId);
  float4 GetLightSourceIntensity(uint a_lightId, const float4* a_wavelengths);

  /**
  \brief offset reflected ray position by epsilon;
  \param  a_lightId   - light id
  \param  ray_pos     - surface point from which we shoot shadow ray (i.e. ShadowRayPos)
  \param  ray_dir     - direction of the shadow ray                  (i.e. shadowRayDir)
  \param  lpos        - position on light surface
  \param  lnorm       - normal   on light surface
  \return PdfW (solid-angle probability density) for sampling target light from point 'ray_pos' with direction 'ray_dir' to surface point on light (lpos, lnorm)
  */
  float  LightEvalPDF(int a_lightId, float3 ray_pos, float3 ray_dir, const float3 lpos, const float3 lnorm);

  float4 GetEnvironmentColorAndPdf(float3 a_dir);

  BsdfSample MaterialSampleWhitted(uint a_materialId, float3 v, float3 n, float2 tc);
  float3     MaterialEvalWhitted  (uint a_materialId, float3 l, float3 v, float3 n, float2 tc);

  BsdfSample MaterialSampleAndEval(uint a_materialId, float4 wavelengths, RandomGen* a_gen, float3 v, float3 n, float2 tc,
                                   MisData* a_misPrev, const uint a_currRayFlags); 
  BsdfEval   MaterialEval         (uint a_materialId, float4 wavelengths, float3 l, float3 v, float3 n, float2 tc);

  uint32_t MaterialBlendSampleAndEval(uint a_materialId, float4 wavelengths, RandomGen* a_gen, float3 v, float3 n, float2 tc, 
                                      MisData* a_misPrev, BsdfSample* a_pRes);
  MatIdWeightPair MaterialBlendEval(MatIdWeight a_mat, float4 wavelengths, float3 l, float3 v, float3 n, float2 tc);

  uint RemapMaterialId(uint a_mId, int a_instId); 
  
  ////////////////////////////////////////////////////////////////////////////////////////////////

  float3 m_camPos = float3(0.0f, 0.85f, 4.5f);
  void InitSceneMaterials(int a_numSpheres, int a_seed = 0);

  std::vector<Material>         m_materials;
  std::vector<uint32_t>         m_matIdOffsets;  ///< offset = m_matIdOffsets[geomId]
  std::vector<uint32_t>         m_matIdByPrimId; ///< matId  = m_matIdByPrimId[offset + primId]
  std::vector<uint32_t>         m_triIndices;    ///< (A,B,C) = m_triIndices[(offset + primId)*3 + 0/1/2]
  std::vector<uint32_t>         m_packedXY;
                                
  std::vector<uint32_t>         m_vertOffset;    ///< vertOffs = m_vertOffset[geomId]
  std::vector<float4>           m_vNorm4f;       ///< vertNorm = m_vNorm4f[vertOffs + vertId]
  std::vector<float2>           m_vTexc2f;       ///< vertTexc = m_vTexc2f[vertOffs + vertId]
                                
  std::vector<int>              m_remapInst;
  std::vector<int>              m_allRemapLists;
  std::vector<int>              m_allRemapListsOffsets;
  std::vector<uint32_t>         m_instIdToLightInstId;
                                
  float4x4                      m_proj;
  float4x4                      m_worldView;
  float4x4                      m_projInv;
  float4x4                      m_worldViewInv;

  std::vector<RandomGen>        m_randomGens;
  std::vector<float4x4>         m_normMatrices; ///< per instance normal matrix, local to world

  std::shared_ptr<ISceneObject> m_pAccelStruct = nullptr;

  std::vector<LightSource> m_lights;
  float4 m_envColor = float4{0.0f};

  uint m_intergatorType = INTEGRATOR_STUPID_PT;
  int  m_spectral_mode  = 0;
  float m_exposureMult = 1.0f;

  float naivePtTime  = 0.0f;
  float shadowPtTime = 0.0f;
  float raytraceTime = 0.0f;
  float fromRaysPtTime = 0.0f;

  //// textures
  //
  std::vector< std::shared_ptr<ICombinedImageSampler> > m_textures; ///< all textures, right now represented via combined image/sampler

  // std::vector<Spectrum> m_spectra;
  std::vector<float> m_wavelengths; 
  std::vector<float> m_spec_values;
  std::vector<uint2> m_spec_offset_sz;
  std::vector<float> m_cie_x;
  std::vector<float> m_cie_y;
  std::vector<float> m_cie_z;

  float4 SampleMatColorParamSpectrum(uint32_t matId, float4 a_wavelengths, uint32_t paramId, uint32_t paramSpecId);
  float4 SampleMatParamSpectrum(uint32_t matId, float4 a_wavelengths, uint32_t paramId, uint32_t paramSpecId);

  static constexpr uint32_t KSPEC_MAT_TYPE_GLTF      = 1;
  static constexpr uint32_t KSPEC_MAT_TYPE_GLASS     = 2;
  static constexpr uint32_t KSPEC_MAT_TYPE_CONDUCTOR = 3;
  static constexpr uint32_t KSPEC_MAT_TYPE_THIN_FILM = 4;
  static constexpr uint32_t KSPEC_MAT_TYPE_DIFFUSE   = 5;
  static constexpr uint32_t KSPEC_SOME_FEATURE_DUMMY = 6;

  static constexpr uint32_t KSPEC_SPECTRAL_RENDERING = 7;
  static constexpr uint32_t KSPEC_MAT_TYPE_BLEND     = 8;
  static constexpr uint32_t KSPEC_BLEND_STACK_SIZE   = 9;
  
  static constexpr uint32_t TOTAL_FEATURES_NUM       = 10; // (!!!) DON'T rename it to KSPEC_TOTAL_FEATURES_NUM.

  //virtual std::vector<uint32_t> ListRequiredFeatures()  { return {1,1,1,1,1,1,1,1,4,1}; } 
  virtual std::vector<uint32_t> ListRequiredFeatures()  { return m_enabledFeatures; } 

  std::vector<uint32_t>         m_enabledFeatures;
  std::vector<uint32_t>         m_actualFeatures;
  std::string                   GetFeatureName(uint32_t a_featureId);

  static std::string g_lastScenePath;
  static std::string g_lastSceneDir;
public:
  static std::vector<uint32_t> PreliminarySceneAnalysis(const char* a_scenePath, const char* a_sncDir,
                                                        int& width, int& height, int& spectral_mode);
};

#endif
