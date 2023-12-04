#include <vector>
#include <array>
#include <memory>
#include <limits>
#include <cassert>
#include "vk_copy.h"
#include "vk_context.h"

#include "integrator_pt1_generated.h"
#include "include/Integrator_generated_ubo.h"

static uint32_t ComputeReductionAuxBufferElements(uint32_t whole_size, uint32_t wg_size)
{
  uint32_t sizeTotal = 0;
  while (whole_size > 1)
  {
    whole_size  = (whole_size + wg_size - 1) / wg_size;
    sizeTotal  += std::max<uint32_t>(whole_size, 1);
  }
  return sizeTotal;
}

VkBufferUsageFlags Integrator_Generated::GetAdditionalFlagsForUBO() const
{
  return VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
}

uint32_t Integrator_Generated::GetDefaultMaxTextures() const { return 256; }

void Integrator_Generated::MakeComputePipelineAndLayout(const char* a_shaderPath, const char* a_mainName, const VkSpecializationInfo *a_specInfo, const VkDescriptorSetLayout a_dsLayout, VkPipelineLayout* pPipelineLayout, VkPipeline* pPipeline)
{
  VkPipelineShaderStageCreateInfo shaderStageInfo = {};
  shaderStageInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
  shaderStageInfo.stage = VK_SHADER_STAGE_COMPUTE_BIT;

  auto shaderCode   = vk_utils::readSPVFile(a_shaderPath);
  auto shaderModule = vk_utils::createShaderModule(device, shaderCode);

  shaderStageInfo.module              = shaderModule;
  shaderStageInfo.pName               = a_mainName;
  shaderStageInfo.pSpecializationInfo = a_specInfo;

  VkPushConstantRange pcRange = {};
  pcRange.stageFlags = shaderStageInfo.stage;
  pcRange.offset     = 0;
  pcRange.size       = 128; // at least 128 bytes for push constants for all Vulkan implementations

  VkPipelineLayoutCreateInfo pipelineLayoutInfo = {};
  pipelineLayoutInfo.sType                  = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
  pipelineLayoutInfo.pushConstantRangeCount = 1;
  pipelineLayoutInfo.pPushConstantRanges    = &pcRange;
  pipelineLayoutInfo.pSetLayouts            = &a_dsLayout;
  pipelineLayoutInfo.setLayoutCount         = 1;
   
  VkResult res = vkCreatePipelineLayout(device, &pipelineLayoutInfo, nullptr, pPipelineLayout);
  if(res != VK_SUCCESS)
  {
    std::string errMsg = vk_utils::errorString(res);
    std::cout << "[ShaderError]: vkCreatePipelineLayout have failed for '" << a_shaderPath << "' with '" << errMsg.c_str() << "'" << std::endl;
  }
  else
    m_allCreatedPipelineLayouts.push_back(*pPipelineLayout);

  VkComputePipelineCreateInfo pipelineInfo = {};
  pipelineInfo.sType              = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
  pipelineInfo.flags              = 0;
  pipelineInfo.stage              = shaderStageInfo;
  pipelineInfo.layout             = (*pPipelineLayout);
  pipelineInfo.basePipelineHandle = VK_NULL_HANDLE;
  res = vkCreateComputePipelines(device, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, pPipeline);
  if(res != VK_SUCCESS)
  {
    std::string errMsg = vk_utils::errorString(res);
    std::cout << "[ShaderError]: vkCreateComputePipelines have failed for '" << a_shaderPath << "' with '" << errMsg.c_str() << "'" << std::endl;
  }
  else
    m_allCreatedPipelines.push_back(*pPipeline);

  if (shaderModule != VK_NULL_HANDLE)
    vkDestroyShaderModule(device, shaderModule, VK_NULL_HANDLE);
}

void Integrator_Generated::MakeComputePipelineOnly(const char* a_shaderPath, const char* a_mainName, const VkSpecializationInfo *a_specInfo, const VkDescriptorSetLayout a_dsLayout, VkPipelineLayout pipelineLayout, VkPipeline* pPipeline)
{
  VkPipelineShaderStageCreateInfo shaderStageInfo = {};
  shaderStageInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
  shaderStageInfo.stage = VK_SHADER_STAGE_COMPUTE_BIT;

  auto shaderCode   = vk_utils::readSPVFile(a_shaderPath);
  auto shaderModule = vk_utils::createShaderModule(device, shaderCode);

  shaderStageInfo.module              = shaderModule;
  shaderStageInfo.pName               = a_mainName;
  shaderStageInfo.pSpecializationInfo = a_specInfo;

  VkComputePipelineCreateInfo pipelineInfo = {};
  pipelineInfo.sType              = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
  pipelineInfo.flags              = 0;
  pipelineInfo.stage              = shaderStageInfo;
  pipelineInfo.layout             = pipelineLayout;
  pipelineInfo.basePipelineHandle = VK_NULL_HANDLE;
  VkResult res = vkCreateComputePipelines(device, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, pPipeline);
  if(res != VK_SUCCESS)
  {
    std::string errMsg = vk_utils::errorString(res);
    std::cout << "[ShaderError]: vkCreateComputePipelines have failed for '" << a_shaderPath << "' with '" << errMsg.c_str() << "'" << std::endl;
  }
  else
    m_allCreatedPipelines.push_back(*pPipeline);

  if (shaderModule != VK_NULL_HANDLE)
    vkDestroyShaderModule(device, shaderModule, VK_NULL_HANDLE);
}


Integrator_Generated::~Integrator_Generated()
{
  for(size_t i=0;i<m_allCreatedPipelines.size();i++)
    vkDestroyPipeline(device, m_allCreatedPipelines[i], nullptr);
  for(size_t i=0;i<m_allCreatedPipelineLayouts.size();i++)
    vkDestroyPipelineLayout(device, m_allCreatedPipelineLayouts[i], nullptr);

  vkDestroyDescriptorSetLayout(device, RayTraceMegaDSLayout, nullptr);
  RayTraceMegaDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, CastSingleRayMegaDSLayout, nullptr);
  CastSingleRayMegaDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, PackXYMegaDSLayout, nullptr);
  PackXYMegaDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, PathTraceFromInputRaysMegaDSLayout, nullptr);
  PathTraceFromInputRaysMegaDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, PathTraceMegaDSLayout, nullptr);
  PathTraceMegaDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, NaivePathTraceMegaDSLayout, nullptr);
  NaivePathTraceMegaDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorPool(device, m_dsPool, NULL); m_dsPool = VK_NULL_HANDLE;

 
  vkDestroyBuffer(device, m_classDataBuffer, nullptr);

  vkDestroyBuffer(device, m_vdata.m_allRemapListsBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_allRemapListsOffsetsBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_cie_xBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_cie_yBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_cie_zBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_instIdToLightInstIdBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_lightsBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_matIdByPrimIdBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_matIdOffsetsBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_materialsBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_normMatricesBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_packedXYBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_randomGensBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_remapInstBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_spec_offset_szBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_spec_valuesBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_triIndicesBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_vNorm4fBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_vTexc2fBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_vertOffsetBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_wavelengthsBuffer, nullptr);
  for(auto obj : m_vdata.m_texturesArrayTexture)
    vkDestroyImage(device, obj, nullptr);
  for(auto obj : m_vdata.m_texturesArrayView)
    vkDestroyImageView(device, obj, nullptr);
  for(auto obj : m_vdata.m_texturesArraySampler)
  vkDestroySampler(device, obj, nullptr);
  FreeAllAllocations(m_allMems);
}

void Integrator_Generated::InitHelpers()
{
  vkGetPhysicalDeviceProperties(physicalDevice, &m_devProps);
}

const VkSpecializationInfo* Integrator_Generated::GetAllSpecInfo()
{
  if(m_allSpecConstInfo.size() == m_allSpecConstVals.size()) // already processed
    return &m_allSpecInfo;
  m_allSpecConstInfo.resize(m_allSpecConstVals.size());
  m_allSpecConstInfo[0].constantID = 0;
  m_allSpecConstInfo[0].size       = sizeof(uint32_t);
  m_allSpecConstInfo[0].offset     = 0;
  m_allSpecConstInfo[1].constantID = 1;
  m_allSpecConstInfo[1].size       = sizeof(uint32_t);
  m_allSpecConstInfo[1].offset     = 1*sizeof(uint32_t);
  m_allSpecConstInfo[2].constantID = 2;
  m_allSpecConstInfo[2].size       = sizeof(uint32_t);
  m_allSpecConstInfo[2].offset     = 2*sizeof(uint32_t);
  m_allSpecConstInfo[3].constantID = 3;
  m_allSpecConstInfo[3].size       = sizeof(uint32_t);
  m_allSpecConstInfo[3].offset     = 3*sizeof(uint32_t);
  m_allSpecConstInfo[4].constantID = 4;
  m_allSpecConstInfo[4].size       = sizeof(uint32_t);
  m_allSpecConstInfo[4].offset     = 4*sizeof(uint32_t);
  m_allSpecConstInfo[5].constantID = 5;
  m_allSpecConstInfo[5].size       = sizeof(uint32_t);
  m_allSpecConstInfo[5].offset     = 5*sizeof(uint32_t);
  m_allSpecConstInfo[6].constantID = 6;
  m_allSpecConstInfo[6].size       = sizeof(uint32_t);
  m_allSpecConstInfo[6].offset     = 6*sizeof(uint32_t);
  m_allSpecConstInfo[7].constantID = 7;
  m_allSpecConstInfo[7].size       = sizeof(uint32_t);
  m_allSpecConstInfo[7].offset     = 7*sizeof(uint32_t);
  m_allSpecConstInfo[8].constantID = 8;
  m_allSpecConstInfo[8].size       = sizeof(uint32_t);
  m_allSpecConstInfo[8].offset     = 8*sizeof(uint32_t);
  m_allSpecInfo.dataSize      = m_allSpecConstVals.size()*sizeof(uint32_t);
  m_allSpecInfo.mapEntryCount = static_cast<uint32_t>(m_allSpecConstInfo.size());
  m_allSpecInfo.pMapEntries   = m_allSpecConstInfo.data();
  m_allSpecInfo.pData         = m_allSpecConstVals.data();
  return &m_allSpecInfo;  
}

void Integrator_Generated::InitKernel_RayTraceMega(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_generated/RayTraceMega.comp.spv"); 
  const VkSpecializationInfo* kspec = GetAllSpecInfo();
  RayTraceMegaDSLayout = CreateRayTraceMegaDSLayout();
  if(m_megaKernelFlags.enableRayTraceMega)
    MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, RayTraceMegaDSLayout, &RayTraceMegaLayout, &RayTraceMegaPipeline);
  else
  {
    RayTraceMegaLayout   = nullptr;
    RayTraceMegaPipeline = nullptr;
  }
}

void Integrator_Generated::InitKernel_CastSingleRayMega(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_generated/CastSingleRayMega.comp.spv"); 
  const VkSpecializationInfo* kspec = GetAllSpecInfo();
  CastSingleRayMegaDSLayout = CreateCastSingleRayMegaDSLayout();
  if(m_megaKernelFlags.enableCastSingleRayMega)
    MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, CastSingleRayMegaDSLayout, &CastSingleRayMegaLayout, &CastSingleRayMegaPipeline);
  else
  {
    CastSingleRayMegaLayout   = nullptr;
    CastSingleRayMegaPipeline = nullptr;
  }
}

void Integrator_Generated::InitKernel_PackXYMega(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_generated/PackXYMega.comp.spv"); 
  const VkSpecializationInfo* kspec = GetAllSpecInfo();
  PackXYMegaDSLayout = CreatePackXYMegaDSLayout();
  if(m_megaKernelFlags.enablePackXYMega)
    MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, PackXYMegaDSLayout, &PackXYMegaLayout, &PackXYMegaPipeline);
  else
  {
    PackXYMegaLayout   = nullptr;
    PackXYMegaPipeline = nullptr;
  }
}

void Integrator_Generated::InitKernel_PathTraceFromInputRaysMega(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_generated/PathTraceFromInputRaysMega.comp.spv"); 
  const VkSpecializationInfo* kspec = GetAllSpecInfo();
  PathTraceFromInputRaysMegaDSLayout = CreatePathTraceFromInputRaysMegaDSLayout();
  if(m_megaKernelFlags.enablePathTraceFromInputRaysMega)
    MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, PathTraceFromInputRaysMegaDSLayout, &PathTraceFromInputRaysMegaLayout, &PathTraceFromInputRaysMegaPipeline);
  else
  {
    PathTraceFromInputRaysMegaLayout   = nullptr;
    PathTraceFromInputRaysMegaPipeline = nullptr;
  }
}

void Integrator_Generated::InitKernel_PathTraceMega(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_generated/PathTraceMega.comp.spv"); 
  const VkSpecializationInfo* kspec = GetAllSpecInfo();
  PathTraceMegaDSLayout = CreatePathTraceMegaDSLayout();
  if(m_megaKernelFlags.enablePathTraceMega)
    MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, PathTraceMegaDSLayout, &PathTraceMegaLayout, &PathTraceMegaPipeline);
  else
  {
    PathTraceMegaLayout   = nullptr;
    PathTraceMegaPipeline = nullptr;
  }
}

void Integrator_Generated::InitKernel_NaivePathTraceMega(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_generated/NaivePathTraceMega.comp.spv"); 
  const VkSpecializationInfo* kspec = GetAllSpecInfo();
  NaivePathTraceMegaDSLayout = CreateNaivePathTraceMegaDSLayout();
  if(m_megaKernelFlags.enableNaivePathTraceMega)
    MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, NaivePathTraceMegaDSLayout, &NaivePathTraceMegaLayout, &NaivePathTraceMegaPipeline);
  else
  {
    NaivePathTraceMegaLayout   = nullptr;
    NaivePathTraceMegaPipeline = nullptr;
  }
}


void Integrator_Generated::InitKernels(const char* a_filePath)
{
  InitKernel_RayTraceMega(a_filePath);
  InitKernel_CastSingleRayMega(a_filePath);
  InitKernel_PackXYMega(a_filePath);
  InitKernel_PathTraceFromInputRaysMega(a_filePath);
  InitKernel_PathTraceMega(a_filePath);
  InitKernel_NaivePathTraceMega(a_filePath);
}

void Integrator_Generated::InitBuffers(size_t a_maxThreadsCount, bool a_tempBuffersOverlay)
{
  m_maxThreadCount = a_maxThreadsCount;
  std::vector<VkBuffer> allBuffers;
  allBuffers.reserve(64);

  struct BufferReqPair
  {
    BufferReqPair() {  }
    BufferReqPair(VkBuffer a_buff, VkDevice a_dev) : buf(a_buff) { vkGetBufferMemoryRequirements(a_dev, a_buff, &req); }
    VkBuffer             buf = VK_NULL_HANDLE;
    VkMemoryRequirements req = {};
  };

  struct LocalBuffers
  {
    std::vector<BufferReqPair> bufs;
    size_t                     size = 0;
    std::vector<VkBuffer>      bufsClean;
  };

  std::vector<LocalBuffers> groups;
  groups.reserve(16);


  size_t largestIndex = 0;
  size_t largestSize  = 0;
  for(size_t i=0;i<groups.size();i++)
  {
    if(groups[i].size > largestSize)
    {
      largestIndex = i;
      largestSize  = groups[i].size;
    }
    groups[i].bufsClean.resize(groups[i].bufs.size());
    for(size_t j=0;j<groups[i].bufsClean.size();j++)
      groups[i].bufsClean[j] = groups[i].bufs[j].buf;
  }
  
  auto& allBuffersRef = allBuffers;

  m_classDataBuffer = vk_utils::createBuffer(device, sizeof(m_uboData),  VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT | GetAdditionalFlagsForUBO());
  allBuffersRef.push_back(m_classDataBuffer);


  auto internalBuffersMem = AllocAndBind(allBuffersRef);
  if(a_tempBuffersOverlay)
  {
    for(size_t i=0;i<groups.size();i++)
      if(i != largestIndex)
        AssignBuffersToMemory(groups[i].bufsClean, internalBuffersMem.memObject);
  }
}

void Integrator_Generated::InitMemberBuffers()
{
  std::vector<VkBuffer> memberVectors;
  std::vector<VkImage>  memberTextures;

  m_vdata.m_allRemapListsBuffer = vk_utils::createBuffer(device, m_allRemapLists.capacity()*sizeof(int), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_allRemapListsBuffer);
  m_vdata.m_allRemapListsOffsetsBuffer = vk_utils::createBuffer(device, m_allRemapListsOffsets.capacity()*sizeof(int), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_allRemapListsOffsetsBuffer);
  m_vdata.m_cie_xBuffer = vk_utils::createBuffer(device, m_cie_x.capacity()*sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_cie_xBuffer);
  m_vdata.m_cie_yBuffer = vk_utils::createBuffer(device, m_cie_y.capacity()*sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_cie_yBuffer);
  m_vdata.m_cie_zBuffer = vk_utils::createBuffer(device, m_cie_z.capacity()*sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_cie_zBuffer);
  m_vdata.m_instIdToLightInstIdBuffer = vk_utils::createBuffer(device, m_instIdToLightInstId.capacity()*sizeof(unsigned int), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_instIdToLightInstIdBuffer);
  m_vdata.m_lightsBuffer = vk_utils::createBuffer(device, m_lights.capacity()*sizeof(struct LightSource), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_lightsBuffer);
  m_vdata.m_matIdByPrimIdBuffer = vk_utils::createBuffer(device, m_matIdByPrimId.capacity()*sizeof(unsigned int), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_matIdByPrimIdBuffer);
  m_vdata.m_matIdOffsetsBuffer = vk_utils::createBuffer(device, m_matIdOffsets.capacity()*sizeof(unsigned int), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_matIdOffsetsBuffer);
  m_vdata.m_materialsBuffer = vk_utils::createBuffer(device, m_materials.capacity()*sizeof(struct Material), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_materialsBuffer);
  m_vdata.m_normMatricesBuffer = vk_utils::createBuffer(device, m_normMatrices.capacity()*sizeof(struct LiteMath::float4x4), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_normMatricesBuffer);
  m_vdata.m_packedXYBuffer = vk_utils::createBuffer(device, m_packedXY.capacity()*sizeof(unsigned int), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_packedXYBuffer);
  m_vdata.m_randomGensBuffer = vk_utils::createBuffer(device, m_randomGens.capacity()*sizeof(struct RandomGenT), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_randomGensBuffer);
  m_vdata.m_remapInstBuffer = vk_utils::createBuffer(device, m_remapInst.capacity()*sizeof(int), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_remapInstBuffer);
  m_vdata.m_spec_offset_szBuffer = vk_utils::createBuffer(device, m_spec_offset_sz.capacity()*sizeof(struct LiteMath::uint2), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_spec_offset_szBuffer);
  m_vdata.m_spec_valuesBuffer = vk_utils::createBuffer(device, m_spec_values.capacity()*sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_spec_valuesBuffer);
  m_vdata.m_triIndicesBuffer = vk_utils::createBuffer(device, m_triIndices.capacity()*sizeof(unsigned int), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_triIndicesBuffer);
  m_vdata.m_vNorm4fBuffer = vk_utils::createBuffer(device, m_vNorm4f.capacity()*sizeof(struct LiteMath::float4), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_vNorm4fBuffer);
  m_vdata.m_vTexc2fBuffer = vk_utils::createBuffer(device, m_vTexc2f.capacity()*sizeof(struct LiteMath::float2), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_vTexc2fBuffer);
  m_vdata.m_vertOffsetBuffer = vk_utils::createBuffer(device, m_vertOffset.capacity()*sizeof(unsigned int), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_vertOffsetBuffer);
  m_vdata.m_wavelengthsBuffer = vk_utils::createBuffer(device, m_wavelengths.capacity()*sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_wavelengthsBuffer);

  m_vdata.m_texturesArrayTexture.resize(0);
  m_vdata.m_texturesArrayView.resize(0);
  m_vdata.m_texturesArraySampler.resize(0);
  m_vdata.m_texturesArrayTexture.reserve(64);
  m_vdata.m_texturesArrayView.reserve(64);
  m_vdata.m_texturesArraySampler.reserve(64);
  for(auto imageObj : m_textures) 
  {
    auto tex = CreateTexture2D(imageObj->width(), imageObj->height(), VkFormat(imageObj->format()), VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_SAMPLED_BIT);
    auto sam = CreateSampler(imageObj->sampler());
    m_vdata.m_texturesArrayTexture.push_back(tex);
    m_vdata.m_texturesArrayView.push_back(VK_NULL_HANDLE);
    m_vdata.m_texturesArraySampler.push_back(sam);
    memberTextures.push_back(tex);
  }

  AllocMemoryForMemberBuffersAndImages(memberVectors, memberTextures);
}



VkImage Integrator_Generated::CreateTexture2D(const int a_width, const int a_height, VkFormat a_format, VkImageUsageFlags a_usage)
{
  VkImage result = VK_NULL_HANDLE;
  VkImageCreateInfo imgCreateInfo = {};
  imgCreateInfo.sType         = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
  imgCreateInfo.pNext         = nullptr;
  imgCreateInfo.flags         = 0; // not sure about this ...
  imgCreateInfo.imageType     = VK_IMAGE_TYPE_2D;
  imgCreateInfo.format        = a_format;
  imgCreateInfo.extent        = VkExtent3D{uint32_t(a_width), uint32_t(a_height), 1};
  imgCreateInfo.mipLevels     = 1;
  imgCreateInfo.samples       = VK_SAMPLE_COUNT_1_BIT;
  imgCreateInfo.tiling        = VK_IMAGE_TILING_OPTIMAL;
  imgCreateInfo.usage         = a_usage; 
  imgCreateInfo.sharingMode   = VK_SHARING_MODE_EXCLUSIVE;
  imgCreateInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
  imgCreateInfo.arrayLayers   = 1;
  VK_CHECK_RESULT(vkCreateImage(device, &imgCreateInfo, nullptr, &result));
  return result;
}

VkSampler Integrator_Generated::CreateSampler(const Sampler& a_sampler) // TODO: implement this function correctly
{
  VkSampler result = VK_NULL_HANDLE;
  VkSamplerCreateInfo samplerInfo = {};
  samplerInfo.sType        = VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO;
  samplerInfo.pNext        = nullptr;
  samplerInfo.flags        = 0;
  samplerInfo.magFilter    = VkFilter(int(a_sampler.filter));
  samplerInfo.minFilter    = VkFilter(int(a_sampler.filter));
  samplerInfo.mipmapMode   = (samplerInfo.magFilter == VK_FILTER_LINEAR ) ? VK_SAMPLER_MIPMAP_MODE_LINEAR : VK_SAMPLER_MIPMAP_MODE_NEAREST;
  samplerInfo.addressModeU = VkSamplerAddressMode(int(a_sampler.addressU));
  samplerInfo.addressModeV = VkSamplerAddressMode(int(a_sampler.addressV));
  samplerInfo.addressModeW = VkSamplerAddressMode(int(a_sampler.addressW));
  samplerInfo.mipLodBias   = a_sampler.mipLODBias;
  samplerInfo.compareOp    = VK_COMPARE_OP_NEVER;
  samplerInfo.minLod           = a_sampler.minLOD;
  samplerInfo.maxLod           = a_sampler.maxLOD;
  samplerInfo.maxAnisotropy    = a_sampler.maxAnisotropy;
  samplerInfo.anisotropyEnable = (a_sampler.maxAnisotropy > 1) ? VK_TRUE : VK_FALSE;
  samplerInfo.borderColor      = VK_BORDER_COLOR_FLOAT_OPAQUE_BLACK;
  samplerInfo.unnormalizedCoordinates = VK_FALSE;
  VK_CHECK_RESULT(vkCreateSampler(device, &samplerInfo, nullptr, &result));
  return result;
}


void Integrator_Generated::AssignBuffersToMemory(const std::vector<VkBuffer>& a_buffers, VkDeviceMemory a_mem)
{
  if(a_buffers.size() == 0 || a_mem == VK_NULL_HANDLE)
    return;

  std::vector<VkMemoryRequirements> memInfos(a_buffers.size());
  for(size_t i=0;i<memInfos.size();i++)
  {
    if(a_buffers[i] != VK_NULL_HANDLE)
      vkGetBufferMemoryRequirements(device, a_buffers[i], &memInfos[i]);
    else
    {
      memInfos[i] = memInfos[0];
      memInfos[i].size = 0;
    }
  }
  
  for(size_t i=1;i<memInfos.size();i++)
  {
    if(memInfos[i].memoryTypeBits != memInfos[0].memoryTypeBits)
    {
      std::cout << "[Integrator_Generated::AssignBuffersToMemory]: error, input buffers has different 'memReq.memoryTypeBits'" << std::endl;
      return;
    }
  }

  auto offsets = vk_utils::calculateMemOffsets(memInfos);
  for (size_t i = 0; i < memInfos.size(); i++)
  {
    if(a_buffers[i] != VK_NULL_HANDLE)
      vkBindBufferMemory(device, a_buffers[i], a_mem, offsets[i]);
  }
}

Integrator_Generated::MemLoc Integrator_Generated::AllocAndBind(const std::vector<VkBuffer>& a_buffers)
{
  MemLoc currLoc;
  if(a_buffers.size() > 0)
  {
    currLoc.memObject = vk_utils::allocateAndBindWithPadding(device, physicalDevice, a_buffers);
    currLoc.allocId   = m_allMems.size();
    m_allMems.push_back(currLoc);
  }
  return currLoc;
}

Integrator_Generated::MemLoc Integrator_Generated::AllocAndBind(const std::vector<VkImage>& a_images)
{
  MemLoc currLoc;
  if(a_images.size() > 0)
  {
    std::vector<VkMemoryRequirements> reqs(a_images.size()); 
    for(size_t i=0; i<reqs.size(); i++)
      vkGetImageMemoryRequirements(device, a_images[i], &reqs[i]);

    for(size_t i=0; i<reqs.size(); i++)
    {
      if(reqs[i].memoryTypeBits != reqs[0].memoryTypeBits)
      {
        std::cout << "Integrator_Generated::AllocAndBind(textures): memoryTypeBits warning, need to split mem allocation (override me)" << std::endl;
        break;
      }
    } 

    auto offsets  = vk_utils::calculateMemOffsets(reqs);
    auto memTotal = offsets[offsets.size() - 1];

    VkMemoryAllocateInfo allocateInfo = {};
    allocateInfo.sType           = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    allocateInfo.pNext           = nullptr;
    allocateInfo.allocationSize  = memTotal;
    allocateInfo.memoryTypeIndex = vk_utils::findMemoryType(reqs[0].memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, physicalDevice);
    VK_CHECK_RESULT(vkAllocateMemory(device, &allocateInfo, NULL, &currLoc.memObject));
    
    for(size_t i=0;i<a_images.size();i++) {
      VK_CHECK_RESULT(vkBindImageMemory(device, a_images[i], currLoc.memObject, offsets[i]));
    }

    currLoc.allocId = m_allMems.size();
    m_allMems.push_back(currLoc);
  }
  return currLoc;
}

void Integrator_Generated::FreeAllAllocations(std::vector<MemLoc>& a_memLoc)
{
  // in general you may check 'mem.allocId' for unique to be sure you dont free mem twice
  // for default implementation this is not needed
  for(auto mem : a_memLoc)
    vkFreeMemory(device, mem.memObject, nullptr);
  a_memLoc.resize(0);
}     

void Integrator_Generated::AllocMemoryForMemberBuffersAndImages(const std::vector<VkBuffer>& a_buffers, const std::vector<VkImage>& a_images)
{
  std::vector<VkMemoryRequirements> bufMemReqs(a_buffers.size()); // we must check that all buffers have same memoryTypeBits;
  for(size_t i = 0; i < a_buffers.size(); ++i)                    // if not, split to multiple allocations
  {
    if(a_buffers[i] != VK_NULL_HANDLE)
      vkGetBufferMemoryRequirements(device, a_buffers[i], &bufMemReqs[i]);
    else
    {
      bufMemReqs[i] = bufMemReqs[0];
      bufMemReqs[i].size = 0;
    }
  }

  bool needSplit = false;
  for(size_t i = 1; i < bufMemReqs.size(); ++i)
  {
    if(bufMemReqs[i].memoryTypeBits != bufMemReqs[0].memoryTypeBits)
    {
      needSplit = true;
      break;
    }
  }

  if(needSplit)
  {
    std::unordered_map<uint32_t, std::vector<uint32_t> > bufferSets;
    for(uint32_t j = 0; j < uint32_t(bufMemReqs.size()); ++j)
    {
      uint32_t key = uint32_t(bufMemReqs[j].memoryTypeBits);
      bufferSets[key].push_back(j);
    }

    for(const auto& buffGroup : bufferSets)
    {
      std::vector<VkBuffer> currGroup;
      for(auto id : buffGroup.second)
        currGroup.push_back(a_buffers[id]);
      AllocAndBind(currGroup);
    }
  }
  else
    AllocAndBind(a_buffers);

  std::vector<VkFormat>             formats;  formats.reserve(0);
  std::vector<VkImageView*>         views;    views.reserve(0);
  std::vector<VkImage>              textures; textures.reserve(0);
  VkMemoryRequirements memoryRequirements;

  for(size_t i=0;i< m_vdata.m_texturesArrayTexture.size(); i++) 
  {
    formats.push_back (VkFormat(m_textures[i]->format()));
    views.push_back   (&m_vdata.m_texturesArrayView[i]);
    textures.push_back(m_vdata.m_texturesArrayTexture[i]);
  }

  AllocAndBind(textures);
  for(size_t i=0;i<textures.size();i++)
  {
    VkImageViewCreateInfo imageViewInfo = {};
    imageViewInfo.sType                           = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
    imageViewInfo.flags                           = 0;
    imageViewInfo.viewType                        = VK_IMAGE_VIEW_TYPE_2D;
    imageViewInfo.format                          = formats[i];
    imageViewInfo.components                      = { VK_COMPONENT_SWIZZLE_R, VK_COMPONENT_SWIZZLE_G, VK_COMPONENT_SWIZZLE_B, VK_COMPONENT_SWIZZLE_A };
    imageViewInfo.subresourceRange.aspectMask     = VK_IMAGE_ASPECT_COLOR_BIT;
    imageViewInfo.subresourceRange.baseMipLevel   = 0;
    imageViewInfo.subresourceRange.baseArrayLayer = 0;
    imageViewInfo.subresourceRange.layerCount     = 1;
    imageViewInfo.subresourceRange.levelCount     = 1;
    imageViewInfo.image                           = textures[i];     // The view will be based on the texture's image
    VK_CHECK_RESULT(vkCreateImageView(device, &imageViewInfo, nullptr, views[i]));
  }
}


VkPhysicalDeviceFeatures2 Integrator_Generated::ListRequiredDeviceFeatures(std::vector<const char*>& deviceExtensions)
{
  static VkPhysicalDeviceFeatures2 features2 = {};
  features2.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2;
  features2.pNext = nullptr; 
  features2.features.shaderInt64   = false;
  features2.features.shaderFloat64 = false;  
  features2.features.shaderInt16   = false;
  void** ppNext = &features2.pNext;
  {
    static VkPhysicalDeviceAccelerationStructureFeaturesKHR enabledAccelStructFeatures = {};
    static VkPhysicalDeviceBufferDeviceAddressFeatures      enabledDeviceAddressFeatures = {};
    static VkPhysicalDeviceRayQueryFeaturesKHR              enabledRayQueryFeatures =  {};
    static VkPhysicalDeviceDescriptorIndexingFeatures       indexingFeatures = {};

    indexingFeatures.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_DESCRIPTOR_INDEXING_FEATURES;
    indexingFeatures.pNext = nullptr;
    indexingFeatures.shaderSampledImageArrayNonUniformIndexing = VK_TRUE; // TODO: move bindless texture to seperate feature!
    indexingFeatures.runtimeDescriptorArray                    = VK_TRUE; // TODO: move bindless texture to seperate feature!

    enabledRayQueryFeatures.sType    = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_RAY_QUERY_FEATURES_KHR;
    enabledRayQueryFeatures.rayQuery = VK_TRUE;
    enabledRayQueryFeatures.pNext    = &indexingFeatures;
  
    enabledDeviceAddressFeatures.sType               = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_BUFFER_DEVICE_ADDRESS_FEATURES;
    enabledDeviceAddressFeatures.bufferDeviceAddress = VK_TRUE;
    enabledDeviceAddressFeatures.pNext               = &enabledRayQueryFeatures;
  
    enabledAccelStructFeatures.sType                 = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_ACCELERATION_STRUCTURE_FEATURES_KHR;
    enabledAccelStructFeatures.accelerationStructure = VK_TRUE;
    enabledAccelStructFeatures.pNext                 = &enabledDeviceAddressFeatures;

    (*ppNext) = &enabledAccelStructFeatures; ppNext = &indexingFeatures.pNext;
    
    // Required by VK_KHR_RAY_QUERY
    deviceExtensions.push_back(VK_KHR_ACCELERATION_STRUCTURE_EXTENSION_NAME);
    deviceExtensions.push_back(VK_KHR_RAY_QUERY_EXTENSION_NAME);
    deviceExtensions.push_back("VK_KHR_spirv_1_4");
    deviceExtensions.push_back("VK_KHR_shader_float_controls");  
    // Required by VK_KHR_acceleration_structure
    deviceExtensions.push_back(VK_KHR_BUFFER_DEVICE_ADDRESS_EXTENSION_NAME);
    deviceExtensions.push_back(VK_KHR_DEFERRED_HOST_OPERATIONS_EXTENSION_NAME);
    deviceExtensions.push_back(VK_EXT_DESCRIPTOR_INDEXING_EXTENSION_NAME);
    // // Required by VK_KHR_ray_tracing_pipeline
    // m_deviceExtensions.push_back(VK_KHR_SPIRV_1_4_EXTENSION_NAME);
    // // Required by VK_KHR_spirv_1_4
    // m_deviceExtensions.push_back(VK_KHR_SHADER_FLOAT_CONTROLS_EXTENSION_NAME);
    deviceExtensions.push_back("VK_EXT_descriptor_indexing"); // TODO: move bindless texture it to seperate feature!
  }
  return features2;
}

Integrator_Generated::MegaKernelIsEnabled Integrator_Generated::m_megaKernelFlags;
