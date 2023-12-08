#pragma once

#include <vector>
#include <memory>
#include <string>
#include <unordered_map>
#include <array>

#include "vk_pipeline.h"
#include "vk_buffers.h"
#include "vk_utils.h"
#include "vk_copy.h"
#include "vk_context.h"

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using namespace LiteMath;

#include "integrator_pt.h"

#include "include/Integrator_generated_ubo.h"
class Integrator_Generated : public Integrator
{
public:

  Integrator_Generated(int a_maxThreads, int a_spectral_mode, std::vector<uint32_t> a_features) : Integrator(a_maxThreads, a_spectral_mode, a_features) 
  {
  }
  virtual void InitVulkanObjects(VkDevice a_device, VkPhysicalDevice a_physicalDevice, size_t a_maxThreadsCount);
  
  virtual void SetVulkanContext(vk_utils::VulkanContext a_ctx) { m_ctx = a_ctx; }
  virtual void SetVulkanInOutFor_RayTrace(
    VkBuffer out_colorBuffer,
    size_t   out_colorOffset,
    uint32_t dummyArgument = 0)
  {
    RayTrace_local.out_colorBuffer = out_colorBuffer;
    RayTrace_local.out_colorOffset = out_colorOffset;
    InitAllGeneratedDescriptorSets_RayTrace();
  }

  virtual void SetVulkanInOutFor_CastSingleRay(
    VkBuffer out_colorBuffer,
    size_t   out_colorOffset,
    uint32_t dummyArgument = 0)
  {
    CastSingleRay_local.out_colorBuffer = out_colorBuffer;
    CastSingleRay_local.out_colorOffset = out_colorOffset;
    InitAllGeneratedDescriptorSets_CastSingleRay();
  }

  virtual void SetVulkanInOutFor_PackXY(
    uint32_t dummyArgument = 0)
  {
    InitAllGeneratedDescriptorSets_PackXY();
  }

  virtual void SetVulkanInOutFor_PathTraceFromInputRays(
    VkBuffer in_rayPosAndNearBuffer,
    size_t   in_rayPosAndNearOffset,
    VkBuffer in_rayDirAndFarBuffer,
    size_t   in_rayDirAndFarOffset,
    VkBuffer out_colorBuffer,
    size_t   out_colorOffset,
    uint32_t dummyArgument = 0)
  {
    PathTraceFromInputRays_local.in_rayPosAndNearBuffer = in_rayPosAndNearBuffer;
    PathTraceFromInputRays_local.in_rayPosAndNearOffset = in_rayPosAndNearOffset;
    PathTraceFromInputRays_local.in_rayDirAndFarBuffer = in_rayDirAndFarBuffer;
    PathTraceFromInputRays_local.in_rayDirAndFarOffset = in_rayDirAndFarOffset;
    PathTraceFromInputRays_local.out_colorBuffer = out_colorBuffer;
    PathTraceFromInputRays_local.out_colorOffset = out_colorOffset;
    InitAllGeneratedDescriptorSets_PathTraceFromInputRays();
  }

  virtual void SetVulkanInOutFor_PathTrace(
    VkBuffer out_colorBuffer,
    size_t   out_colorOffset,
    uint32_t dummyArgument = 0)
  {
    PathTrace_local.out_colorBuffer = out_colorBuffer;
    PathTrace_local.out_colorOffset = out_colorOffset;
    InitAllGeneratedDescriptorSets_PathTrace();
  }

  virtual void SetVulkanInOutFor_NaivePathTrace(
    VkBuffer out_colorBuffer,
    size_t   out_colorOffset,
    uint32_t dummyArgument = 0)
  {
    NaivePathTrace_local.out_colorBuffer = out_colorBuffer;
    NaivePathTrace_local.out_colorOffset = out_colorOffset;
    InitAllGeneratedDescriptorSets_NaivePathTrace();
  }

  virtual ~Integrator_Generated();

  
  virtual void InitMemberBuffers();
  virtual void UpdateAll(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
  {
    UpdatePlainMembers(a_pCopyEngine);
    UpdateVectorMembers(a_pCopyEngine);
    UpdateTextureMembers(a_pCopyEngine);
  }
  
  virtual void CommitDeviceData() override // you have to define this virtual function in the original imput class
  {
    InitMemberBuffers();
    UpdateAll(m_ctx.pCopyHelper);
  }  
  void GetExecutionTime(const char* a_funcName, float a_out[4]) override; 
  void UpdateMembersPlainData() override { UpdatePlainMembers(m_ctx.pCopyHelper); } 
  
  virtual void UpdatePlainMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine);
  virtual void UpdateVectorMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine);
  virtual void UpdateTextureMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine);
  virtual void ReadPlainMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine);
  static VkPhysicalDeviceFeatures2 ListRequiredDeviceFeatures(std::vector<const char*>& deviceExtensions);
  
  virtual void RayTraceCmd(VkCommandBuffer a_commandBuffer, uint tid, float4* out_color);
  virtual void CastSingleRayCmd(VkCommandBuffer a_commandBuffer, uint tid, uint* out_color);
  virtual void PackXYCmd(VkCommandBuffer a_commandBuffer, uint tidX, uint tidY);
  virtual void PathTraceFromInputRaysCmd(VkCommandBuffer a_commandBuffer, uint tid, const RayPosAndW* in_rayPosAndNear, const RayDirAndT* in_rayDirAndFar, float4* out_color);
  virtual void PathTraceCmd(VkCommandBuffer a_commandBuffer, uint tid, float4* out_color);
  virtual void NaivePathTraceCmd(VkCommandBuffer a_commandBuffer, uint tid, float4* out_color);

  void RayTraceBlock(uint tid, float4* out_color, uint32_t a_numPasses) override;
  void CastSingleRayBlock(uint tid, uint* out_color, uint32_t a_numPasses) override;
  void PackXYBlock(uint tidX, uint tidY, uint32_t a_numPasses) override;
  void PathTraceFromInputRaysBlock(uint tid, const RayPosAndW* in_rayPosAndNear, const RayDirAndT* in_rayDirAndFar, float4* out_color, uint32_t a_numPasses) override;
  void PathTraceBlock(uint tid, float4* out_color, uint32_t a_numPasses) override;
  void NaivePathTraceBlock(uint tid, float4* out_color, uint32_t a_numPasses) override;

  inline vk_utils::ExecTime GetRayTraceExecutionTime() const { return m_exTimeRayTrace; }
  inline vk_utils::ExecTime GetCastSingleRayExecutionTime() const { return m_exTimeCastSingleRay; }
  inline vk_utils::ExecTime GetPackXYExecutionTime() const { return m_exTimePackXY; }
  inline vk_utils::ExecTime GetPathTraceFromInputRaysExecutionTime() const { return m_exTimePathTraceFromInputRays; }
  inline vk_utils::ExecTime GetPathTraceExecutionTime() const { return m_exTimePathTrace; }
  inline vk_utils::ExecTime GetNaivePathTraceExecutionTime() const { return m_exTimeNaivePathTrace; }

  vk_utils::ExecTime m_exTimeRayTrace;
  vk_utils::ExecTime m_exTimeCastSingleRay;
  vk_utils::ExecTime m_exTimePackXY;
  vk_utils::ExecTime m_exTimePathTraceFromInputRays;
  vk_utils::ExecTime m_exTimePathTrace;
  vk_utils::ExecTime m_exTimeNaivePathTrace;

  virtual void copyKernelFloatCmd(uint32_t length);
  
  virtual void RayTraceMegaCmd(uint tid, float4* out_color);
  virtual void CastSingleRayMegaCmd(uint tid, uint* out_color);
  virtual void PackXYMegaCmd(uint tidX, uint tidY);
  virtual void PathTraceFromInputRaysMegaCmd(uint tid, const RayPosAndW* in_rayPosAndNear, const RayDirAndT* in_rayDirAndFar, float4* out_color);
  virtual void PathTraceMegaCmd(uint tid, float4* out_color);
  virtual void NaivePathTraceMegaCmd(uint tid, float4* out_color);
  
  struct MemLoc
  {
    VkDeviceMemory memObject = VK_NULL_HANDLE;
    size_t         memOffset = 0;
    size_t         allocId   = 0;
  };

  virtual MemLoc AllocAndBind(const std::vector<VkBuffer>& a_buffers); ///< replace this function to apply custom allocator
  virtual MemLoc AllocAndBind(const std::vector<VkImage>& a_image);    ///< replace this function to apply custom allocator
  virtual void   FreeAllAllocations(std::vector<MemLoc>& a_memLoc);    ///< replace this function to apply custom allocator

protected:

  VkPhysicalDevice           physicalDevice = VK_NULL_HANDLE;
  VkDevice                   device         = VK_NULL_HANDLE;
  vk_utils::VulkanContext    m_ctx          = {};
  VkCommandBuffer            m_currCmdBuffer   = VK_NULL_HANDLE;
  uint32_t                   m_currThreadFlags = 0;
  std::vector<MemLoc>        m_allMems;
  VkPhysicalDeviceProperties m_devProps;

  VkBufferMemoryBarrier BarrierForClearFlags(VkBuffer a_buffer);
  VkBufferMemoryBarrier BarrierForSingleBuffer(VkBuffer a_buffer);
  void BarriersForSeveralBuffers(VkBuffer* a_inBuffers, VkBufferMemoryBarrier* a_outBarriers, uint32_t a_buffersNum);

  virtual void InitHelpers();
  virtual void InitBuffers(size_t a_maxThreadsCount, bool a_tempBuffersOverlay = true);
  virtual void InitKernels(const char* a_filePath);
  virtual void AllocateAllDescriptorSets();

  virtual void InitAllGeneratedDescriptorSets_RayTrace();
  virtual void InitAllGeneratedDescriptorSets_CastSingleRay();
  virtual void InitAllGeneratedDescriptorSets_PackXY();
  virtual void InitAllGeneratedDescriptorSets_PathTraceFromInputRays();
  virtual void InitAllGeneratedDescriptorSets_PathTrace();
  virtual void InitAllGeneratedDescriptorSets_NaivePathTrace();

  virtual void AssignBuffersToMemory(const std::vector<VkBuffer>& a_buffers, VkDeviceMemory a_mem);

  virtual void AllocMemoryForMemberBuffersAndImages(const std::vector<VkBuffer>& a_buffers, const std::vector<VkImage>& a_image);
  virtual std::string AlterShaderPath(const char* in_shaderPath) { return std::string("") + std::string(in_shaderPath); }

  
  

  struct RayTrace_Data
  {
    VkBuffer out_colorBuffer = VK_NULL_HANDLE;
    size_t   out_colorOffset = 0;
    bool needToClearOutput = true;
  } RayTrace_local;

  struct CastSingleRay_Data
  {
    VkBuffer out_colorBuffer = VK_NULL_HANDLE;
    size_t   out_colorOffset = 0;
    bool needToClearOutput = true;
  } CastSingleRay_local;

  struct PackXY_Data
  {
    bool needToClearOutput = true;
  } PackXY_local;

  struct PathTraceFromInputRays_Data
  {
    VkBuffer in_rayPosAndNearBuffer = VK_NULL_HANDLE;
    size_t   in_rayPosAndNearOffset = 0;
    VkBuffer in_rayDirAndFarBuffer = VK_NULL_HANDLE;
    size_t   in_rayDirAndFarOffset = 0;
    VkBuffer out_colorBuffer = VK_NULL_HANDLE;
    size_t   out_colorOffset = 0;
    bool needToClearOutput = true;
  } PathTraceFromInputRays_local;

  struct PathTrace_Data
  {
    VkBuffer out_colorBuffer = VK_NULL_HANDLE;
    size_t   out_colorOffset = 0;
    bool needToClearOutput = true;
  } PathTrace_local;

  struct NaivePathTrace_Data
  {
    VkBuffer out_colorBuffer = VK_NULL_HANDLE;
    size_t   out_colorOffset = 0;
    bool needToClearOutput = true;
  } NaivePathTrace_local;



  struct MembersDataGPU
  {
    VkBuffer m_allRemapListsBuffer = VK_NULL_HANDLE;
    size_t   m_allRemapListsOffset = 0;
    VkBuffer m_allRemapListsOffsetsBuffer = VK_NULL_HANDLE;
    size_t   m_allRemapListsOffsetsOffset = 0;
    VkBuffer m_cie_xBuffer = VK_NULL_HANDLE;
    size_t   m_cie_xOffset = 0;
    VkBuffer m_cie_yBuffer = VK_NULL_HANDLE;
    size_t   m_cie_yOffset = 0;
    VkBuffer m_cie_zBuffer = VK_NULL_HANDLE;
    size_t   m_cie_zOffset = 0;
    VkBuffer m_instIdToLightInstIdBuffer = VK_NULL_HANDLE;
    size_t   m_instIdToLightInstIdOffset = 0;
    VkBuffer m_lightsBuffer = VK_NULL_HANDLE;
    size_t   m_lightsOffset = 0;
    VkBuffer m_matIdByPrimIdBuffer = VK_NULL_HANDLE;
    size_t   m_matIdByPrimIdOffset = 0;
    VkBuffer m_matIdOffsetsBuffer = VK_NULL_HANDLE;
    size_t   m_matIdOffsetsOffset = 0;
    VkBuffer m_materialsBuffer = VK_NULL_HANDLE;
    size_t   m_materialsOffset = 0;
    VkBuffer m_normMatricesBuffer = VK_NULL_HANDLE;
    size_t   m_normMatricesOffset = 0;
    VkBuffer m_packedXYBuffer = VK_NULL_HANDLE;
    size_t   m_packedXYOffset = 0;
    VkBuffer m_randomGensBuffer = VK_NULL_HANDLE;
    size_t   m_randomGensOffset = 0;
    VkBuffer m_remapInstBuffer = VK_NULL_HANDLE;
    size_t   m_remapInstOffset = 0;
    VkBuffer m_spec_offset_szBuffer = VK_NULL_HANDLE;
    size_t   m_spec_offset_szOffset = 0;
    VkBuffer m_spec_valuesBuffer = VK_NULL_HANDLE;
    size_t   m_spec_valuesOffset = 0;
    VkBuffer m_triIndicesBuffer = VK_NULL_HANDLE;
    size_t   m_triIndicesOffset = 0;
    VkBuffer m_vNorm4fBuffer = VK_NULL_HANDLE;
    size_t   m_vNorm4fOffset = 0;
    VkBuffer m_vTang4fBuffer = VK_NULL_HANDLE;
    size_t   m_vTang4fOffset = 0;
    VkBuffer m_vertOffsetBuffer = VK_NULL_HANDLE;
    size_t   m_vertOffsetOffset = 0;
    VkBuffer m_wavelengthsBuffer = VK_NULL_HANDLE;
    size_t   m_wavelengthsOffset = 0;
    std::vector<VkImage>     m_texturesArrayTexture;
    std::vector<VkImageView> m_texturesArrayView   ;
    std::vector<VkSampler>   m_texturesArraySampler; ///<! samplers for texture arrays are always used
    size_t                   m_texturesArrayMaxSize; ///<! used when texture array size is not known after constructor of base class is finished
  } m_vdata;
  
  
  VkImage   CreateTexture2D(const int a_width, const int a_height, VkFormat a_format, VkImageUsageFlags a_usage);
  VkSampler CreateSampler(const Sampler& a_sampler);
  struct TexAccessPair
  {
    TexAccessPair() : image(VK_NULL_HANDLE), access(0) {}
    TexAccessPair(VkImage a_image, VkAccessFlags a_access) : image(a_image), access(a_access) {}
    VkImage image;
    VkAccessFlags access;  
  };
  void TrackTextureAccess(const std::vector<TexAccessPair>& a_pairs, std::unordered_map<uint64_t, VkAccessFlags>& a_currImageFlags);
  size_t m_maxThreadCount = 0;
  VkBuffer m_classDataBuffer = VK_NULL_HANDLE;

  VkPipelineLayout      RayTraceMegaLayout   = VK_NULL_HANDLE;
  VkPipeline            RayTraceMegaPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout RayTraceMegaDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreateRayTraceMegaDSLayout();
  virtual void InitKernel_RayTraceMega(const char* a_filePath);
  VkPipelineLayout      CastSingleRayMegaLayout   = VK_NULL_HANDLE;
  VkPipeline            CastSingleRayMegaPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout CastSingleRayMegaDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreateCastSingleRayMegaDSLayout();
  virtual void InitKernel_CastSingleRayMega(const char* a_filePath);
  VkPipelineLayout      PackXYMegaLayout   = VK_NULL_HANDLE;
  VkPipeline            PackXYMegaPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout PackXYMegaDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatePackXYMegaDSLayout();
  virtual void InitKernel_PackXYMega(const char* a_filePath);
  VkPipelineLayout      PathTraceFromInputRaysMegaLayout   = VK_NULL_HANDLE;
  VkPipeline            PathTraceFromInputRaysMegaPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout PathTraceFromInputRaysMegaDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatePathTraceFromInputRaysMegaDSLayout();
  virtual void InitKernel_PathTraceFromInputRaysMega(const char* a_filePath);
  VkPipelineLayout      PathTraceMegaLayout   = VK_NULL_HANDLE;
  VkPipeline            PathTraceMegaPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout PathTraceMegaDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatePathTraceMegaDSLayout();
  virtual void InitKernel_PathTraceMega(const char* a_filePath);
  VkPipelineLayout      NaivePathTraceMegaLayout   = VK_NULL_HANDLE;
  VkPipeline            NaivePathTraceMegaPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout NaivePathTraceMegaDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreateNaivePathTraceMegaDSLayout();
  virtual void InitKernel_NaivePathTraceMega(const char* a_filePath);


  virtual VkBufferUsageFlags GetAdditionalFlagsForUBO() const;
  virtual uint32_t           GetDefaultMaxTextures() const;

  VkPipelineLayout      copyKernelFloatLayout   = VK_NULL_HANDLE;
  VkPipeline            copyKernelFloatPipeline = VK_NULL_HANDLE;
  VkDescriptorSetLayout copyKernelFloatDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatecopyKernelFloatDSLayout();

  VkDescriptorPool m_dsPool = VK_NULL_HANDLE;
  VkDescriptorSet  m_allGeneratedDS[6];

  Integrator_Generated_UBO_Data m_uboData;
  
  constexpr static uint32_t MEMCPY_BLOCK_SIZE = 256;
  constexpr static uint32_t REDUCTION_BLOCK_SIZE = 256;

  virtual void MakeComputePipelineAndLayout(const char* a_shaderPath, const char* a_mainName, const VkSpecializationInfo *a_specInfo, const VkDescriptorSetLayout a_dsLayout, 
                                            VkPipelineLayout* pPipelineLayout, VkPipeline* pPipeline);
  virtual void MakeComputePipelineOnly(const char* a_shaderPath, const char* a_mainName, const VkSpecializationInfo *a_specInfo, const VkDescriptorSetLayout a_dsLayout, VkPipelineLayout pipelineLayout, 
                                       VkPipeline* pPipeline);

  std::vector<VkPipelineLayout> m_allCreatedPipelineLayouts; ///<! remenber them here to delete later
  std::vector<VkPipeline>       m_allCreatedPipelines;       ///<! remenber them here to delete later
  std::vector<uint32_t>                  m_allSpecConstVals; ///<! require user to define "ListRequiredFeatures" func.    
  std::vector<VkSpecializationMapEntry>  m_allSpecConstInfo;
  VkSpecializationInfo                   m_allSpecInfo;
  const VkSpecializationInfo*            GetAllSpecInfo();
public:

  struct MegaKernelIsEnabled
  {
    bool enableRayTraceMega = true;
    bool enableCastSingleRayMega = true;
    bool enablePackXYMega = true;
    bool enablePathTraceFromInputRaysMega = true;
    bool enablePathTraceMega = true;
    bool enableNaivePathTraceMega = true;
    bool dummy = 0;
  };

  static MegaKernelIsEnabled  m_megaKernelFlags;
  static MegaKernelIsEnabled& EnabledPipelines() { return m_megaKernelFlags; }

};


