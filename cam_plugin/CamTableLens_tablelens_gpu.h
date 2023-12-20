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


#include "CamTableLens.h"

#include "include/CamTableLens_tablelens_gpu_ubo.h"
class CamTableLens_TABLELENS_GPU : public CamTableLens
{
public:

  CamTableLens_TABLELENS_GPU() 
  {
  }
  virtual void InitVulkanObjects(VkDevice a_device, VkPhysicalDevice a_physicalDevice, size_t a_maxThreadsCount);
  
  virtual void SetVulkanContext(vk_utils::VulkanContext a_ctx) { m_ctx = a_ctx; }
  virtual void SetVulkanInOutFor_AddSamplesContributionBlock(
    VkBuffer out_color4fBuffer,
    size_t   out_color4fOffset,
    VkBuffer colors4fBuffer,
    size_t   colors4fOffset,
    uint32_t dummyArgument = 0)
  {
    AddSamplesContributionBlock_local.out_color4fBuffer = out_color4fBuffer;
    AddSamplesContributionBlock_local.out_color4fOffset = out_color4fOffset;
    AddSamplesContributionBlock_local.colors4fBuffer = colors4fBuffer;
    AddSamplesContributionBlock_local.colors4fOffset = colors4fOffset;
    InitAllGeneratedDescriptorSets_AddSamplesContributionBlock();
  }

  virtual void SetVulkanInOutFor_MakeRaysBlock(
    VkBuffer out_rayPosAndNear4fBuffer,
    size_t   out_rayPosAndNear4fOffset,
    VkBuffer out_rayDirAndFar4fBuffer,
    size_t   out_rayDirAndFar4fOffset,
    uint32_t dummyArgument = 0)
  {
    MakeRaysBlock_local.out_rayPosAndNear4fBuffer = out_rayPosAndNear4fBuffer;
    MakeRaysBlock_local.out_rayPosAndNear4fOffset = out_rayPosAndNear4fOffset;
    MakeRaysBlock_local.out_rayDirAndFar4fBuffer = out_rayDirAndFar4fBuffer;
    MakeRaysBlock_local.out_rayDirAndFar4fOffset = out_rayDirAndFar4fOffset;
    InitAllGeneratedDescriptorSets_MakeRaysBlock();
  }

  virtual ~CamTableLens_TABLELENS_GPU();

  
  virtual void InitMemberBuffers();
  virtual void UpdateAll(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
  {
    UpdatePlainMembers(a_pCopyEngine);
    UpdateVectorMembers(a_pCopyEngine);
    UpdateTextureMembers(a_pCopyEngine);
  }
  
  virtual void CommitDeviceData(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyHelper) // you have to define this virtual function in the original imput class
  {
    InitMemberBuffers();
    UpdateAll(a_pCopyHelper);
  }  
  void CommitDeviceData() override { CommitDeviceData(m_ctx.pCopyHelper); }  
  void GetExecutionTime(const char* a_funcName, float a_out[4]) override; 
  

  virtual void ReserveEmptyVectors();
  virtual void UpdatePlainMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine);
  virtual void UpdateVectorMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine);
  virtual void UpdateTextureMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine);
  virtual void ReadPlainMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine);
  static VkPhysicalDeviceFeatures2 ListRequiredDeviceFeatures(std::vector<const char*>& deviceExtensions);
  
  virtual void AddSamplesContributionBlockCmd(VkCommandBuffer a_commandBuffer, float* out_color4f, const float* colors4f, uint32_t in_blockSize, uint32_t a_width, uint32_t a_height, int subPassId);
  virtual void MakeRaysBlockCmd(VkCommandBuffer a_commandBuffer, RayPosAndW* out_rayPosAndNear4f, RayDirAndT* out_rayDirAndFar4f, uint32_t in_blockSize, int subPassId);

  void AddSamplesContributionBlock(float* out_color4f, const float* colors4f, uint32_t in_blockSize, uint32_t a_width, uint32_t a_height, int subPassId) override;
  void MakeRaysBlock(RayPosAndW* out_rayPosAndNear4f, RayDirAndT* out_rayDirAndFar4f, uint32_t in_blockSize, int subPassId) override;

  inline vk_utils::ExecTime GetAddSamplesContributionBlockExecutionTime() const { return m_exTimeAddSamplesContributionBlock; }
  inline vk_utils::ExecTime GetMakeRaysBlockExecutionTime() const { return m_exTimeMakeRaysBlock; }

  vk_utils::ExecTime m_exTimeAddSamplesContributionBlock;
  vk_utils::ExecTime m_exTimeMakeRaysBlock;

  virtual void copyKernelFloatCmd(uint32_t length);
  
  virtual void MakeEyeRayCmd(int in_blockSize, RayPosAndW* out_rayPosAndNear4f, RayDirAndT* out_rayDirAndFar4f, int subPassId);
  virtual void ContribSampleCmd(int in_blockSize, const float* in_color, float* out_color, int subPassId);
  
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

  virtual void InitAllGeneratedDescriptorSets_AddSamplesContributionBlock();
  virtual void InitAllGeneratedDescriptorSets_MakeRaysBlock();

  virtual void AssignBuffersToMemory(const std::vector<VkBuffer>& a_buffers, VkDeviceMemory a_mem);

  virtual void AllocMemoryForMemberBuffersAndImages(const std::vector<VkBuffer>& a_buffers, const std::vector<VkImage>& a_image);
  virtual std::string AlterShaderPath(const char* in_shaderPath) { return std::string("cam_plugin/") + std::string(in_shaderPath); }

  
  

  struct AddSamplesContributionBlock_Data
  {
    VkBuffer out_color4fBuffer = VK_NULL_HANDLE;
    size_t   out_color4fOffset = 0;
    VkBuffer colors4fBuffer = VK_NULL_HANDLE;
    size_t   colors4fOffset = 0;
    bool needToClearOutput = false;
  } AddSamplesContributionBlock_local;

  struct MakeRaysBlock_Data
  {
    VkBuffer out_rayPosAndNear4fBuffer = VK_NULL_HANDLE;
    size_t   out_rayPosAndNear4fOffset = 0;
    VkBuffer out_rayDirAndFar4fBuffer = VK_NULL_HANDLE;
    size_t   out_rayDirAndFar4fOffset = 0;
    bool needToClearOutput = false;
  } MakeRaysBlock_local;



  struct MembersDataGPU
  {
    VkBuffer linesBuffer = VK_NULL_HANDLE;
    size_t   linesOffset = 0;
    VkBuffer m_cie_xBuffer = VK_NULL_HANDLE;
    size_t   m_cie_xOffset = 0;
    VkBuffer m_cie_yBuffer = VK_NULL_HANDLE;
    size_t   m_cie_yOffset = 0;
    VkBuffer m_cie_zBuffer = VK_NULL_HANDLE;
    size_t   m_cie_zOffset = 0;
    VkBuffer m_randomGensBuffer = VK_NULL_HANDLE;
    size_t   m_randomGensOffset = 0;
    VkBuffer m_storedCos4Buffer = VK_NULL_HANDLE;
    size_t   m_storedCos4Offset = 0;
    VkBuffer m_storedWavesBuffer = VK_NULL_HANDLE;
    size_t   m_storedWavesOffset = 0;
  } m_vdata;
  
  
  size_t m_maxThreadCount = 0;
  VkBuffer m_classDataBuffer = VK_NULL_HANDLE;

  VkPipelineLayout      MakeEyeRayLayout   = VK_NULL_HANDLE;
  VkPipeline            MakeEyeRayPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout MakeEyeRayDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreateMakeEyeRayDSLayout();
  virtual void InitKernel_MakeEyeRay(const char* a_filePath);
  VkPipelineLayout      ContribSampleLayout   = VK_NULL_HANDLE;
  VkPipeline            ContribSamplePipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout ContribSampleDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreateContribSampleDSLayout();
  virtual void InitKernel_ContribSample(const char* a_filePath);


  virtual VkBufferUsageFlags GetAdditionalFlagsForUBO() const;
  virtual uint32_t           GetDefaultMaxTextures() const;

  VkPipelineLayout      copyKernelFloatLayout   = VK_NULL_HANDLE;
  VkPipeline            copyKernelFloatPipeline = VK_NULL_HANDLE;
  VkDescriptorSetLayout copyKernelFloatDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatecopyKernelFloatDSLayout();

  VkDescriptorPool m_dsPool = VK_NULL_HANDLE;
  VkDescriptorSet  m_allGeneratedDS[2];

  CamTableLens_TABLELENS_GPU_UBO_Data m_uboData;
  
  constexpr static uint32_t MEMCPY_BLOCK_SIZE = 256;
  constexpr static uint32_t REDUCTION_BLOCK_SIZE = 256;

  virtual void MakeComputePipelineAndLayout(const char* a_shaderPath, const char* a_mainName, const VkSpecializationInfo *a_specInfo, const VkDescriptorSetLayout a_dsLayout, 
                                            VkPipelineLayout* pPipelineLayout, VkPipeline* pPipeline);
  virtual void MakeComputePipelineOnly(const char* a_shaderPath, const char* a_mainName, const VkSpecializationInfo *a_specInfo, const VkDescriptorSetLayout a_dsLayout, VkPipelineLayout pipelineLayout, 
                                       VkPipeline* pPipeline);

  std::vector<VkPipelineLayout> m_allCreatedPipelineLayouts; ///<! remenber them here to delete later
  std::vector<VkPipeline>       m_allCreatedPipelines;       ///<! remenber them here to delete later
public:

  struct MegaKernelIsEnabled
  {
    bool dummy = 0;
  };

  static MegaKernelIsEnabled  m_megaKernelFlags;
  static MegaKernelIsEnabled& EnabledPipelines() { return m_megaKernelFlags; }

};


