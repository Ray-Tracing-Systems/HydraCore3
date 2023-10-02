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

Integrator_Generated::~Integrator_Generated()
{
  m_pMaker = nullptr;
  vkDestroyDescriptorSetLayout(device, RayTraceMegaDSLayout, nullptr);
  RayTraceMegaDSLayout = VK_NULL_HANDLE;

  vkDestroyPipeline(device, RayTraceMegaPipeline, nullptr);
  vkDestroyPipelineLayout(device, RayTraceMegaLayout, nullptr);
  RayTraceMegaLayout   = VK_NULL_HANDLE;
  RayTraceMegaPipeline = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, CastSingleRayMegaDSLayout, nullptr);
  CastSingleRayMegaDSLayout = VK_NULL_HANDLE;

  vkDestroyPipeline(device, CastSingleRayMegaPipeline, nullptr);
  vkDestroyPipelineLayout(device, CastSingleRayMegaLayout, nullptr);
  CastSingleRayMegaLayout   = VK_NULL_HANDLE;
  CastSingleRayMegaPipeline = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, PackXYMegaDSLayout, nullptr);
  PackXYMegaDSLayout = VK_NULL_HANDLE;

  vkDestroyPipeline(device, PackXYMegaPipeline, nullptr);
  vkDestroyPipelineLayout(device, PackXYMegaLayout, nullptr);
  PackXYMegaLayout   = VK_NULL_HANDLE;
  PackXYMegaPipeline = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, PathTraceMegaDSLayout, nullptr);
  PathTraceMegaDSLayout = VK_NULL_HANDLE;

  vkDestroyPipeline(device, PathTraceMegaPipeline, nullptr);
  vkDestroyPipelineLayout(device, PathTraceMegaLayout, nullptr);
  PathTraceMegaLayout   = VK_NULL_HANDLE;
  PathTraceMegaPipeline = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, NaivePathTraceMegaDSLayout, nullptr);
  NaivePathTraceMegaDSLayout = VK_NULL_HANDLE;

  vkDestroyPipeline(device, NaivePathTraceMegaPipeline, nullptr);
  vkDestroyPipelineLayout(device, NaivePathTraceMegaLayout, nullptr);
  NaivePathTraceMegaLayout   = VK_NULL_HANDLE;
  NaivePathTraceMegaPipeline = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, copyKernelFloatDSLayout, nullptr);
  vkDestroyDescriptorPool(device, m_dsPool, NULL); m_dsPool = VK_NULL_HANDLE;

 
  vkDestroyBuffer(device, m_classDataBuffer, nullptr);

  vkDestroyBuffer(device, m_vdata.m_allRemapListsBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_allRemapListsOffsetsBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_instIdToLightInstIdBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_lightsBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_matIdByPrimIdBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_matIdOffsetsBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_materialsBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_normMatricesBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_packedXYBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_randomGensBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_remapInstBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_triIndicesBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_vNorm4fBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_vTexc2fBuffer, nullptr);
  vkDestroyBuffer(device, m_vdata.m_vertOffsetBuffer, nullptr);
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
  m_pMaker = std::make_unique<vk_utils::ComputePipelineMaker>();
}

VkDescriptorSetLayout Integrator_Generated::CreateRayTraceMegaDSLayout()
{
  std::array<VkDescriptorSetLayoutBinding, 16+1> dsBindings;

  // binding for out_color
  dsBindings[0].binding            = 0;
  dsBindings[0].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[0].descriptorCount    = 1;
  dsBindings[0].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[0].pImmutableSamplers = nullptr;

  // binding for m_normMatrices
  dsBindings[1].binding            = 1;
  dsBindings[1].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[1].descriptorCount    = 1;
  dsBindings[1].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[1].pImmutableSamplers = nullptr;

  // binding for m_allRemapListsOffsets
  dsBindings[2].binding            = 2;
  dsBindings[2].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[2].descriptorCount    = 1;
  dsBindings[2].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[2].pImmutableSamplers = nullptr;

  // binding for m_allRemapLists
  dsBindings[3].binding            = 3;
  dsBindings[3].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[3].descriptorCount    = 1;
  dsBindings[3].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[3].pImmutableSamplers = nullptr;

  // binding for m_vNorm4f
  dsBindings[4].binding            = 4;
  dsBindings[4].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[4].descriptorCount    = 1;
  dsBindings[4].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[4].pImmutableSamplers = nullptr;

  // binding for m_pAccelStruct
  dsBindings[5].binding            = 5;
  dsBindings[5].descriptorType     = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
  dsBindings[5].descriptorCount    = 1;
  dsBindings[5].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[5].pImmutableSamplers = nullptr;

  // binding for m_lights
  dsBindings[6].binding            = 6;
  dsBindings[6].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[6].descriptorCount    = 1;
  dsBindings[6].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[6].pImmutableSamplers = nullptr;

  // binding for m_remapInst
  dsBindings[7].binding            = 7;
  dsBindings[7].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[7].descriptorCount    = 1;
  dsBindings[7].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[7].pImmutableSamplers = nullptr;

  // binding for m_packedXY
  dsBindings[8].binding            = 8;
  dsBindings[8].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[8].descriptorCount    = 1;
  dsBindings[8].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[8].pImmutableSamplers = nullptr;

  // binding for m_textures
  dsBindings[9].binding            = 9;
  dsBindings[9].descriptorType     = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
  m_vdata.m_texturesArrayMaxSize = m_textures.size();
  if(m_vdata.m_texturesArrayMaxSize == 0)
    m_vdata.m_texturesArrayMaxSize = GetDefaultMaxTextures();
  dsBindings[9].descriptorCount    = m_vdata.m_texturesArrayMaxSize;
  dsBindings[9].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[9].pImmutableSamplers = nullptr;

  // binding for m_triIndices
  dsBindings[10].binding            = 10;
  dsBindings[10].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[10].descriptorCount    = 1;
  dsBindings[10].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[10].pImmutableSamplers = nullptr;

  // binding for m_vTexc2f
  dsBindings[11].binding            = 11;
  dsBindings[11].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[11].descriptorCount    = 1;
  dsBindings[11].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[11].pImmutableSamplers = nullptr;

  // binding for m_materials
  dsBindings[12].binding            = 12;
  dsBindings[12].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[12].descriptorCount    = 1;
  dsBindings[12].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[12].pImmutableSamplers = nullptr;

  // binding for m_matIdOffsets
  dsBindings[13].binding            = 13;
  dsBindings[13].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[13].descriptorCount    = 1;
  dsBindings[13].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[13].pImmutableSamplers = nullptr;

  // binding for m_vertOffset
  dsBindings[14].binding            = 14;
  dsBindings[14].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[14].descriptorCount    = 1;
  dsBindings[14].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[14].pImmutableSamplers = nullptr;

  // binding for m_matIdByPrimId
  dsBindings[15].binding            = 15;
  dsBindings[15].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[15].descriptorCount    = 1;
  dsBindings[15].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[15].pImmutableSamplers = nullptr;

  // binding for POD members stored in m_classDataBuffer
  dsBindings[16].binding            = 16;
  dsBindings[16].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[16].descriptorCount    = 1;
  dsBindings[16].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[16].pImmutableSamplers = nullptr;
  
  VkDescriptorSetLayoutCreateInfo descriptorSetLayoutCreateInfo = {};
  descriptorSetLayoutCreateInfo.sType        = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
  descriptorSetLayoutCreateInfo.bindingCount = uint32_t(dsBindings.size());
  descriptorSetLayoutCreateInfo.pBindings    = dsBindings.data();
  
  VkDescriptorSetLayout layout = nullptr;
  VK_CHECK_RESULT(vkCreateDescriptorSetLayout(device, &descriptorSetLayoutCreateInfo, NULL, &layout));
  return layout;
}
VkDescriptorSetLayout Integrator_Generated::CreateCastSingleRayMegaDSLayout()
{
  std::array<VkDescriptorSetLayoutBinding, 11+1> dsBindings;

  // binding for out_color
  dsBindings[0].binding            = 0;
  dsBindings[0].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[0].descriptorCount    = 1;
  dsBindings[0].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[0].pImmutableSamplers = nullptr;

  // binding for m_matIdByPrimId
  dsBindings[1].binding            = 1;
  dsBindings[1].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[1].descriptorCount    = 1;
  dsBindings[1].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[1].pImmutableSamplers = nullptr;

  // binding for m_matIdOffsets
  dsBindings[2].binding            = 2;
  dsBindings[2].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[2].descriptorCount    = 1;
  dsBindings[2].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[2].pImmutableSamplers = nullptr;

  // binding for m_pAccelStruct
  dsBindings[3].binding            = 3;
  dsBindings[3].descriptorType     = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
  dsBindings[3].descriptorCount    = 1;
  dsBindings[3].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[3].pImmutableSamplers = nullptr;

  // binding for m_materials
  dsBindings[4].binding            = 4;
  dsBindings[4].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[4].descriptorCount    = 1;
  dsBindings[4].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[4].pImmutableSamplers = nullptr;

  // binding for m_textures
  dsBindings[5].binding            = 5;
  dsBindings[5].descriptorType     = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
  m_vdata.m_texturesArrayMaxSize = m_textures.size();
  if(m_vdata.m_texturesArrayMaxSize == 0)
    m_vdata.m_texturesArrayMaxSize = GetDefaultMaxTextures();
  dsBindings[5].descriptorCount    = m_vdata.m_texturesArrayMaxSize;
  dsBindings[5].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[5].pImmutableSamplers = nullptr;

  // binding for m_remapInst
  dsBindings[6].binding            = 6;
  dsBindings[6].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[6].descriptorCount    = 1;
  dsBindings[6].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[6].pImmutableSamplers = nullptr;

  // binding for m_allRemapLists
  dsBindings[7].binding            = 7;
  dsBindings[7].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[7].descriptorCount    = 1;
  dsBindings[7].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[7].pImmutableSamplers = nullptr;

  // binding for m_allRemapListsOffsets
  dsBindings[8].binding            = 8;
  dsBindings[8].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[8].descriptorCount    = 1;
  dsBindings[8].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[8].pImmutableSamplers = nullptr;

  // binding for m_packedXY
  dsBindings[9].binding            = 9;
  dsBindings[9].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[9].descriptorCount    = 1;
  dsBindings[9].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[9].pImmutableSamplers = nullptr;

  // binding for m_lights
  dsBindings[10].binding            = 10;
  dsBindings[10].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[10].descriptorCount    = 1;
  dsBindings[10].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[10].pImmutableSamplers = nullptr;

  // binding for POD members stored in m_classDataBuffer
  dsBindings[11].binding            = 11;
  dsBindings[11].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[11].descriptorCount    = 1;
  dsBindings[11].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[11].pImmutableSamplers = nullptr;
  
  VkDescriptorSetLayoutCreateInfo descriptorSetLayoutCreateInfo = {};
  descriptorSetLayoutCreateInfo.sType        = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
  descriptorSetLayoutCreateInfo.bindingCount = uint32_t(dsBindings.size());
  descriptorSetLayoutCreateInfo.pBindings    = dsBindings.data();
  
  VkDescriptorSetLayout layout = nullptr;
  VK_CHECK_RESULT(vkCreateDescriptorSetLayout(device, &descriptorSetLayoutCreateInfo, NULL, &layout));
  return layout;
}
VkDescriptorSetLayout Integrator_Generated::CreatePackXYMegaDSLayout()
{
  std::array<VkDescriptorSetLayoutBinding, 7+1> dsBindings;

  // binding for m_materials
  dsBindings[0].binding            = 0;
  dsBindings[0].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[0].descriptorCount    = 1;
  dsBindings[0].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[0].pImmutableSamplers = nullptr;

  // binding for m_textures
  dsBindings[1].binding            = 1;
  dsBindings[1].descriptorType     = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
  m_vdata.m_texturesArrayMaxSize = m_textures.size();
  if(m_vdata.m_texturesArrayMaxSize == 0)
    m_vdata.m_texturesArrayMaxSize = GetDefaultMaxTextures();
  dsBindings[1].descriptorCount    = m_vdata.m_texturesArrayMaxSize;
  dsBindings[1].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[1].pImmutableSamplers = nullptr;

  // binding for m_remapInst
  dsBindings[2].binding            = 2;
  dsBindings[2].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[2].descriptorCount    = 1;
  dsBindings[2].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[2].pImmutableSamplers = nullptr;

  // binding for m_allRemapLists
  dsBindings[3].binding            = 3;
  dsBindings[3].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[3].descriptorCount    = 1;
  dsBindings[3].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[3].pImmutableSamplers = nullptr;

  // binding for m_allRemapListsOffsets
  dsBindings[4].binding            = 4;
  dsBindings[4].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[4].descriptorCount    = 1;
  dsBindings[4].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[4].pImmutableSamplers = nullptr;

  // binding for m_packedXY
  dsBindings[5].binding            = 5;
  dsBindings[5].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[5].descriptorCount    = 1;
  dsBindings[5].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[5].pImmutableSamplers = nullptr;

  // binding for m_lights
  dsBindings[6].binding            = 6;
  dsBindings[6].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[6].descriptorCount    = 1;
  dsBindings[6].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[6].pImmutableSamplers = nullptr;

  // binding for POD members stored in m_classDataBuffer
  dsBindings[7].binding            = 7;
  dsBindings[7].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[7].descriptorCount    = 1;
  dsBindings[7].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[7].pImmutableSamplers = nullptr;
  
  VkDescriptorSetLayoutCreateInfo descriptorSetLayoutCreateInfo = {};
  descriptorSetLayoutCreateInfo.sType        = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
  descriptorSetLayoutCreateInfo.bindingCount = uint32_t(dsBindings.size());
  descriptorSetLayoutCreateInfo.pBindings    = dsBindings.data();
  
  VkDescriptorSetLayout layout = nullptr;
  VK_CHECK_RESULT(vkCreateDescriptorSetLayout(device, &descriptorSetLayoutCreateInfo, NULL, &layout));
  return layout;
}
VkDescriptorSetLayout Integrator_Generated::CreatePathTraceMegaDSLayout()
{
  std::array<VkDescriptorSetLayoutBinding, 18+1> dsBindings;

  // binding for out_color
  dsBindings[0].binding            = 0;
  dsBindings[0].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[0].descriptorCount    = 1;
  dsBindings[0].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[0].pImmutableSamplers = nullptr;

  // binding for m_matIdByPrimId
  dsBindings[1].binding            = 1;
  dsBindings[1].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[1].descriptorCount    = 1;
  dsBindings[1].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[1].pImmutableSamplers = nullptr;

  // binding for m_normMatrices
  dsBindings[2].binding            = 2;
  dsBindings[2].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[2].descriptorCount    = 1;
  dsBindings[2].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[2].pImmutableSamplers = nullptr;

  // binding for m_allRemapListsOffsets
  dsBindings[3].binding            = 3;
  dsBindings[3].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[3].descriptorCount    = 1;
  dsBindings[3].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[3].pImmutableSamplers = nullptr;

  // binding for m_allRemapLists
  dsBindings[4].binding            = 4;
  dsBindings[4].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[4].descriptorCount    = 1;
  dsBindings[4].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[4].pImmutableSamplers = nullptr;

  // binding for m_randomGens
  dsBindings[5].binding            = 5;
  dsBindings[5].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[5].descriptorCount    = 1;
  dsBindings[5].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[5].pImmutableSamplers = nullptr;

  // binding for m_vNorm4f
  dsBindings[6].binding            = 6;
  dsBindings[6].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[6].descriptorCount    = 1;
  dsBindings[6].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[6].pImmutableSamplers = nullptr;

  // binding for m_remapInst
  dsBindings[7].binding            = 7;
  dsBindings[7].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[7].descriptorCount    = 1;
  dsBindings[7].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[7].pImmutableSamplers = nullptr;

  // binding for m_packedXY
  dsBindings[8].binding            = 8;
  dsBindings[8].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[8].descriptorCount    = 1;
  dsBindings[8].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[8].pImmutableSamplers = nullptr;

  // binding for m_textures
  dsBindings[9].binding            = 9;
  dsBindings[9].descriptorType     = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
  m_vdata.m_texturesArrayMaxSize = m_textures.size();
  if(m_vdata.m_texturesArrayMaxSize == 0)
    m_vdata.m_texturesArrayMaxSize = GetDefaultMaxTextures();
  dsBindings[9].descriptorCount    = m_vdata.m_texturesArrayMaxSize;
  dsBindings[9].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[9].pImmutableSamplers = nullptr;

  // binding for m_instIdToLightInstId
  dsBindings[10].binding            = 10;
  dsBindings[10].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[10].descriptorCount    = 1;
  dsBindings[10].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[10].pImmutableSamplers = nullptr;

  // binding for m_triIndices
  dsBindings[11].binding            = 11;
  dsBindings[11].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[11].descriptorCount    = 1;
  dsBindings[11].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[11].pImmutableSamplers = nullptr;

  // binding for m_vTexc2f
  dsBindings[12].binding            = 12;
  dsBindings[12].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[12].descriptorCount    = 1;
  dsBindings[12].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[12].pImmutableSamplers = nullptr;

  // binding for m_materials
  dsBindings[13].binding            = 13;
  dsBindings[13].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[13].descriptorCount    = 1;
  dsBindings[13].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[13].pImmutableSamplers = nullptr;

  // binding for m_pAccelStruct
  dsBindings[14].binding            = 14;
  dsBindings[14].descriptorType     = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
  dsBindings[14].descriptorCount    = 1;
  dsBindings[14].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[14].pImmutableSamplers = nullptr;

  // binding for m_lights
  dsBindings[15].binding            = 15;
  dsBindings[15].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[15].descriptorCount    = 1;
  dsBindings[15].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[15].pImmutableSamplers = nullptr;

  // binding for m_matIdOffsets
  dsBindings[16].binding            = 16;
  dsBindings[16].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[16].descriptorCount    = 1;
  dsBindings[16].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[16].pImmutableSamplers = nullptr;

  // binding for m_vertOffset
  dsBindings[17].binding            = 17;
  dsBindings[17].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[17].descriptorCount    = 1;
  dsBindings[17].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[17].pImmutableSamplers = nullptr;

  // binding for POD members stored in m_classDataBuffer
  dsBindings[18].binding            = 18;
  dsBindings[18].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[18].descriptorCount    = 1;
  dsBindings[18].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[18].pImmutableSamplers = nullptr;
  
  VkDescriptorSetLayoutCreateInfo descriptorSetLayoutCreateInfo = {};
  descriptorSetLayoutCreateInfo.sType        = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
  descriptorSetLayoutCreateInfo.bindingCount = uint32_t(dsBindings.size());
  descriptorSetLayoutCreateInfo.pBindings    = dsBindings.data();
  
  VkDescriptorSetLayout layout = nullptr;
  VK_CHECK_RESULT(vkCreateDescriptorSetLayout(device, &descriptorSetLayoutCreateInfo, NULL, &layout));
  return layout;
}
VkDescriptorSetLayout Integrator_Generated::CreateNaivePathTraceMegaDSLayout()
{
  std::array<VkDescriptorSetLayoutBinding, 18+1> dsBindings;

  // binding for out_color
  dsBindings[0].binding            = 0;
  dsBindings[0].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[0].descriptorCount    = 1;
  dsBindings[0].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[0].pImmutableSamplers = nullptr;

  // binding for m_matIdByPrimId
  dsBindings[1].binding            = 1;
  dsBindings[1].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[1].descriptorCount    = 1;
  dsBindings[1].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[1].pImmutableSamplers = nullptr;

  // binding for m_normMatrices
  dsBindings[2].binding            = 2;
  dsBindings[2].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[2].descriptorCount    = 1;
  dsBindings[2].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[2].pImmutableSamplers = nullptr;

  // binding for m_allRemapListsOffsets
  dsBindings[3].binding            = 3;
  dsBindings[3].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[3].descriptorCount    = 1;
  dsBindings[3].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[3].pImmutableSamplers = nullptr;

  // binding for m_allRemapLists
  dsBindings[4].binding            = 4;
  dsBindings[4].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[4].descriptorCount    = 1;
  dsBindings[4].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[4].pImmutableSamplers = nullptr;

  // binding for m_randomGens
  dsBindings[5].binding            = 5;
  dsBindings[5].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[5].descriptorCount    = 1;
  dsBindings[5].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[5].pImmutableSamplers = nullptr;

  // binding for m_vNorm4f
  dsBindings[6].binding            = 6;
  dsBindings[6].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[6].descriptorCount    = 1;
  dsBindings[6].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[6].pImmutableSamplers = nullptr;

  // binding for m_lights
  dsBindings[7].binding            = 7;
  dsBindings[7].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[7].descriptorCount    = 1;
  dsBindings[7].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[7].pImmutableSamplers = nullptr;

  // binding for m_remapInst
  dsBindings[8].binding            = 8;
  dsBindings[8].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[8].descriptorCount    = 1;
  dsBindings[8].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[8].pImmutableSamplers = nullptr;

  // binding for m_instIdToLightInstId
  dsBindings[9].binding            = 9;
  dsBindings[9].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[9].descriptorCount    = 1;
  dsBindings[9].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[9].pImmutableSamplers = nullptr;

  // binding for m_packedXY
  dsBindings[10].binding            = 10;
  dsBindings[10].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[10].descriptorCount    = 1;
  dsBindings[10].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[10].pImmutableSamplers = nullptr;

  // binding for m_textures
  dsBindings[11].binding            = 11;
  dsBindings[11].descriptorType     = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
  m_vdata.m_texturesArrayMaxSize = m_textures.size();
  if(m_vdata.m_texturesArrayMaxSize == 0)
    m_vdata.m_texturesArrayMaxSize = GetDefaultMaxTextures();
  dsBindings[11].descriptorCount    = m_vdata.m_texturesArrayMaxSize;
  dsBindings[11].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[11].pImmutableSamplers = nullptr;

  // binding for m_triIndices
  dsBindings[12].binding            = 12;
  dsBindings[12].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[12].descriptorCount    = 1;
  dsBindings[12].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[12].pImmutableSamplers = nullptr;

  // binding for m_vTexc2f
  dsBindings[13].binding            = 13;
  dsBindings[13].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[13].descriptorCount    = 1;
  dsBindings[13].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[13].pImmutableSamplers = nullptr;

  // binding for m_materials
  dsBindings[14].binding            = 14;
  dsBindings[14].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[14].descriptorCount    = 1;
  dsBindings[14].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[14].pImmutableSamplers = nullptr;

  // binding for m_pAccelStruct
  dsBindings[15].binding            = 15;
  dsBindings[15].descriptorType     = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
  dsBindings[15].descriptorCount    = 1;
  dsBindings[15].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[15].pImmutableSamplers = nullptr;

  // binding for m_matIdOffsets
  dsBindings[16].binding            = 16;
  dsBindings[16].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[16].descriptorCount    = 1;
  dsBindings[16].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[16].pImmutableSamplers = nullptr;

  // binding for m_vertOffset
  dsBindings[17].binding            = 17;
  dsBindings[17].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[17].descriptorCount    = 1;
  dsBindings[17].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[17].pImmutableSamplers = nullptr;

  // binding for POD members stored in m_classDataBuffer
  dsBindings[18].binding            = 18;
  dsBindings[18].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[18].descriptorCount    = 1;
  dsBindings[18].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[18].pImmutableSamplers = nullptr;
  
  VkDescriptorSetLayoutCreateInfo descriptorSetLayoutCreateInfo = {};
  descriptorSetLayoutCreateInfo.sType        = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
  descriptorSetLayoutCreateInfo.bindingCount = uint32_t(dsBindings.size());
  descriptorSetLayoutCreateInfo.pBindings    = dsBindings.data();
  
  VkDescriptorSetLayout layout = nullptr;
  VK_CHECK_RESULT(vkCreateDescriptorSetLayout(device, &descriptorSetLayoutCreateInfo, NULL, &layout));
  return layout;
}

VkDescriptorSetLayout Integrator_Generated::CreatecopyKernelFloatDSLayout()
{
  std::array<VkDescriptorSetLayoutBinding, 2> dsBindings;

  dsBindings[0].binding            = 0;
  dsBindings[0].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[0].descriptorCount    = 1;
  dsBindings[0].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[0].pImmutableSamplers = nullptr;

  dsBindings[1].binding            = 1;
  dsBindings[1].descriptorType     = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  dsBindings[1].descriptorCount    = 1;
  dsBindings[1].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
  dsBindings[1].pImmutableSamplers = nullptr;

  VkDescriptorSetLayoutCreateInfo descriptorSetLayoutCreateInfo = {};
  descriptorSetLayoutCreateInfo.sType        = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
  descriptorSetLayoutCreateInfo.bindingCount = dsBindings.size();
  descriptorSetLayoutCreateInfo.pBindings    = dsBindings.data();

  VkDescriptorSetLayout layout = nullptr;
  VK_CHECK_RESULT(vkCreateDescriptorSetLayout(device, &descriptorSetLayoutCreateInfo, NULL, &layout));
  return layout;
}


void Integrator_Generated::InitKernel_RayTraceMega(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_generated/RayTraceMega.comp.spv"); 
  
  m_pMaker->LoadShader(device, shaderPath.c_str(), nullptr, "main");
  RayTraceMegaDSLayout = CreateRayTraceMegaDSLayout();
  RayTraceMegaLayout   = m_pMaker->MakeLayout(device, { RayTraceMegaDSLayout }, 128); // at least 128 bytes for push constants
  RayTraceMegaPipeline = m_pMaker->MakePipeline(device);  
}

void Integrator_Generated::InitKernel_CastSingleRayMega(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_generated/CastSingleRayMega.comp.spv"); 
  
  m_pMaker->LoadShader(device, shaderPath.c_str(), nullptr, "main");
  CastSingleRayMegaDSLayout = CreateCastSingleRayMegaDSLayout();
  CastSingleRayMegaLayout   = m_pMaker->MakeLayout(device, { CastSingleRayMegaDSLayout }, 128); // at least 128 bytes for push constants
  CastSingleRayMegaPipeline = m_pMaker->MakePipeline(device);  
}

void Integrator_Generated::InitKernel_PackXYMega(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_generated/PackXYMega.comp.spv"); 
  
  m_pMaker->LoadShader(device, shaderPath.c_str(), nullptr, "main");
  PackXYMegaDSLayout = CreatePackXYMegaDSLayout();
  PackXYMegaLayout   = m_pMaker->MakeLayout(device, { PackXYMegaDSLayout }, 128); // at least 128 bytes for push constants
  PackXYMegaPipeline = m_pMaker->MakePipeline(device);  
}

void Integrator_Generated::InitKernel_PathTraceMega(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_generated/PathTraceMega.comp.spv"); 
  
  m_pMaker->LoadShader(device, shaderPath.c_str(), nullptr, "main");
  PathTraceMegaDSLayout = CreatePathTraceMegaDSLayout();
  PathTraceMegaLayout   = m_pMaker->MakeLayout(device, { PathTraceMegaDSLayout }, 128); // at least 128 bytes for push constants
  PathTraceMegaPipeline = m_pMaker->MakePipeline(device);  
}

void Integrator_Generated::InitKernel_NaivePathTraceMega(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_generated/NaivePathTraceMega.comp.spv"); 
  
  m_pMaker->LoadShader(device, shaderPath.c_str(), nullptr, "main");
  NaivePathTraceMegaDSLayout = CreateNaivePathTraceMegaDSLayout();
  NaivePathTraceMegaLayout   = m_pMaker->MakeLayout(device, { NaivePathTraceMegaDSLayout }, 128); // at least 128 bytes for push constants
  NaivePathTraceMegaPipeline = m_pMaker->MakePipeline(device);  
}


void Integrator_Generated::InitKernels(const char* a_filePath)
{
  InitKernel_RayTraceMega(a_filePath);
  InitKernel_CastSingleRayMega(a_filePath);
  InitKernel_PackXYMega(a_filePath);
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
  m_vdata.m_triIndicesBuffer = vk_utils::createBuffer(device, m_triIndices.capacity()*sizeof(unsigned int), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_triIndicesBuffer);
  m_vdata.m_vNorm4fBuffer = vk_utils::createBuffer(device, m_vNorm4f.capacity()*sizeof(struct LiteMath::float4), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_vNorm4fBuffer);
  m_vdata.m_vTexc2fBuffer = vk_utils::createBuffer(device, m_vTexc2f.capacity()*sizeof(struct LiteMath::float2), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_vTexc2fBuffer);
  m_vdata.m_vertOffsetBuffer = vk_utils::createBuffer(device, m_vertOffset.capacity()*sizeof(unsigned int), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.m_vertOffsetBuffer);

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


