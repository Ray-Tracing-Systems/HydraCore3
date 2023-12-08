#include <vector>
#include <memory>
#include <limits>
#include <cassert>
#include <chrono>
#include <array>

#include "vk_copy.h"
#include "vk_context.h"
#include "vk_images.h"

#include "integrator_pt1_generated.h"
#include "include/Integrator_generated_ubo.h"


#include "CrossRT.h"
ISceneObject* CreateVulkanRTX(VkDevice a_device, VkPhysicalDevice a_physDevice, uint32_t a_graphicsQId, std::shared_ptr<vk_utils::ICopyEngine> a_pCopyHelper,
                              uint32_t a_maxMeshes, uint32_t a_maxTotalVertices, uint32_t a_maxTotalPrimitives, uint32_t a_maxPrimitivesPerMesh,
                              bool build_as_add);

std::shared_ptr<Integrator> CreateIntegrator_Generated(int a_maxThreads, int a_spectral_mode, std::vector<uint32_t> a_features, vk_utils::VulkanContext a_ctx, size_t a_maxThreadsGenerated) 
{ 
  auto pObj = std::make_shared<Integrator_Generated>(a_maxThreads, a_spectral_mode, a_features); 
  pObj->SetVulkanContext(a_ctx);
  pObj->InitVulkanObjects(a_ctx.device, a_ctx.physicalDevice, a_maxThreadsGenerated); 
  return pObj;
}

static uint32_t ComputeReductionSteps(uint32_t whole_size, uint32_t wg_size)
{
  uint32_t steps = 0;
  while (whole_size > 1)
  {
    steps++;
    whole_size = (whole_size + wg_size - 1) / wg_size;
  }
  return steps;
}

constexpr uint32_t KGEN_FLAG_RETURN            = 1;
constexpr uint32_t KGEN_FLAG_BREAK             = 2;
constexpr uint32_t KGEN_FLAG_DONT_SET_EXIT     = 4;
constexpr uint32_t KGEN_FLAG_SET_EXIT_NEGATIVE = 8;
constexpr uint32_t KGEN_REDUCTION_LAST_STEP    = 16;

void Integrator_Generated::InitVulkanObjects(VkDevice a_device, VkPhysicalDevice a_physicalDevice, size_t a_maxThreadsCount) 
{
  physicalDevice = a_physicalDevice;
  device         = a_device;
  m_allCreatedPipelineLayouts.reserve(256);
  m_allCreatedPipelines.reserve(256);
  m_allSpecConstVals = ListRequiredFeatures();
  InitHelpers();
  InitBuffers(a_maxThreadsCount, true);
  InitKernels(".spv");
  AllocateAllDescriptorSets();

  auto queueAllFID = vk_utils::getQueueFamilyIndex(physicalDevice, VK_QUEUE_GRAPHICS_BIT | VK_QUEUE_TRANSFER_BIT);
  uint32_t userRestrictions[4];
  this->SceneRestrictions(userRestrictions);
  uint32_t maxMeshes            = userRestrictions[0];
  uint32_t maxTotalVertices     = userRestrictions[1];
  uint32_t maxTotalPrimitives   = userRestrictions[2];
  uint32_t maxPrimitivesPerMesh = userRestrictions[3];
  m_pAccelStruct = std::shared_ptr<ISceneObject>(CreateVulkanRTX(a_device, a_physicalDevice, queueAllFID, m_ctx.pCopyHelper,
                                                             maxMeshes, maxTotalVertices, maxTotalPrimitives, maxPrimitivesPerMesh, true),
                                                            [](ISceneObject *p) { DeleteSceneRT(p); } );
}

void Integrator_Generated::UpdatePlainMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{
  const size_t maxAllowedSize = std::numeric_limits<uint32_t>::max();
  m_uboData.m_projInv = m_projInv;
  m_uboData.m_worldViewInv = m_worldViewInv;
  m_uboData.m_envColor = m_envColor;
  m_uboData.m_exposureMult = m_exposureMult;
  m_uboData.m_intergatorType = m_intergatorType;
  m_uboData.m_maxThreadId = m_maxThreadId;
  m_uboData.m_skipBounce = m_skipBounce;
  m_uboData.m_spectral_mode = m_spectral_mode;
  m_uboData.m_tileSize = m_tileSize;
  m_uboData.m_traceDepth = m_traceDepth;
  m_uboData.m_winHeight = m_winHeight;
  m_uboData.m_winWidth = m_winWidth;
  m_uboData.m_allRemapLists_size     = uint32_t( m_allRemapLists.size() );     assert( m_allRemapLists.size() < maxAllowedSize );
  m_uboData.m_allRemapLists_capacity = uint32_t( m_allRemapLists.capacity() ); assert( m_allRemapLists.capacity() < maxAllowedSize );
  m_uboData.m_allRemapListsOffsets_size     = uint32_t( m_allRemapListsOffsets.size() );     assert( m_allRemapListsOffsets.size() < maxAllowedSize );
  m_uboData.m_allRemapListsOffsets_capacity = uint32_t( m_allRemapListsOffsets.capacity() ); assert( m_allRemapListsOffsets.capacity() < maxAllowedSize );
  m_uboData.m_cie_x_size     = uint32_t( m_cie_x.size() );     assert( m_cie_x.size() < maxAllowedSize );
  m_uboData.m_cie_x_capacity = uint32_t( m_cie_x.capacity() ); assert( m_cie_x.capacity() < maxAllowedSize );
  m_uboData.m_cie_y_size     = uint32_t( m_cie_y.size() );     assert( m_cie_y.size() < maxAllowedSize );
  m_uboData.m_cie_y_capacity = uint32_t( m_cie_y.capacity() ); assert( m_cie_y.capacity() < maxAllowedSize );
  m_uboData.m_cie_z_size     = uint32_t( m_cie_z.size() );     assert( m_cie_z.size() < maxAllowedSize );
  m_uboData.m_cie_z_capacity = uint32_t( m_cie_z.capacity() ); assert( m_cie_z.capacity() < maxAllowedSize );
  m_uboData.m_instIdToLightInstId_size     = uint32_t( m_instIdToLightInstId.size() );     assert( m_instIdToLightInstId.size() < maxAllowedSize );
  m_uboData.m_instIdToLightInstId_capacity = uint32_t( m_instIdToLightInstId.capacity() ); assert( m_instIdToLightInstId.capacity() < maxAllowedSize );
  m_uboData.m_lights_size     = uint32_t( m_lights.size() );     assert( m_lights.size() < maxAllowedSize );
  m_uboData.m_lights_capacity = uint32_t( m_lights.capacity() ); assert( m_lights.capacity() < maxAllowedSize );
  m_uboData.m_matIdByPrimId_size     = uint32_t( m_matIdByPrimId.size() );     assert( m_matIdByPrimId.size() < maxAllowedSize );
  m_uboData.m_matIdByPrimId_capacity = uint32_t( m_matIdByPrimId.capacity() ); assert( m_matIdByPrimId.capacity() < maxAllowedSize );
  m_uboData.m_matIdOffsets_size     = uint32_t( m_matIdOffsets.size() );     assert( m_matIdOffsets.size() < maxAllowedSize );
  m_uboData.m_matIdOffsets_capacity = uint32_t( m_matIdOffsets.capacity() ); assert( m_matIdOffsets.capacity() < maxAllowedSize );
  m_uboData.m_materials_size     = uint32_t( m_materials.size() );     assert( m_materials.size() < maxAllowedSize );
  m_uboData.m_materials_capacity = uint32_t( m_materials.capacity() ); assert( m_materials.capacity() < maxAllowedSize );
  m_uboData.m_normMatrices_size     = uint32_t( m_normMatrices.size() );     assert( m_normMatrices.size() < maxAllowedSize );
  m_uboData.m_normMatrices_capacity = uint32_t( m_normMatrices.capacity() ); assert( m_normMatrices.capacity() < maxAllowedSize );
  m_uboData.m_packedXY_size     = uint32_t( m_packedXY.size() );     assert( m_packedXY.size() < maxAllowedSize );
  m_uboData.m_packedXY_capacity = uint32_t( m_packedXY.capacity() ); assert( m_packedXY.capacity() < maxAllowedSize );
  m_uboData.m_randomGens_size     = uint32_t( m_randomGens.size() );     assert( m_randomGens.size() < maxAllowedSize );
  m_uboData.m_randomGens_capacity = uint32_t( m_randomGens.capacity() ); assert( m_randomGens.capacity() < maxAllowedSize );
  m_uboData.m_remapInst_size     = uint32_t( m_remapInst.size() );     assert( m_remapInst.size() < maxAllowedSize );
  m_uboData.m_remapInst_capacity = uint32_t( m_remapInst.capacity() ); assert( m_remapInst.capacity() < maxAllowedSize );
  m_uboData.m_spec_offset_sz_size     = uint32_t( m_spec_offset_sz.size() );     assert( m_spec_offset_sz.size() < maxAllowedSize );
  m_uboData.m_spec_offset_sz_capacity = uint32_t( m_spec_offset_sz.capacity() ); assert( m_spec_offset_sz.capacity() < maxAllowedSize );
  m_uboData.m_spec_values_size     = uint32_t( m_spec_values.size() );     assert( m_spec_values.size() < maxAllowedSize );
  m_uboData.m_spec_values_capacity = uint32_t( m_spec_values.capacity() ); assert( m_spec_values.capacity() < maxAllowedSize );
  m_uboData.m_triIndices_size     = uint32_t( m_triIndices.size() );     assert( m_triIndices.size() < maxAllowedSize );
  m_uboData.m_triIndices_capacity = uint32_t( m_triIndices.capacity() ); assert( m_triIndices.capacity() < maxAllowedSize );
  m_uboData.m_vNorm4f_size     = uint32_t( m_vNorm4f.size() );     assert( m_vNorm4f.size() < maxAllowedSize );
  m_uboData.m_vNorm4f_capacity = uint32_t( m_vNorm4f.capacity() ); assert( m_vNorm4f.capacity() < maxAllowedSize );
  m_uboData.m_vTang4f_size     = uint32_t( m_vTang4f.size() );     assert( m_vTang4f.size() < maxAllowedSize );
  m_uboData.m_vTang4f_capacity = uint32_t( m_vTang4f.capacity() ); assert( m_vTang4f.capacity() < maxAllowedSize );
  m_uboData.m_vertOffset_size     = uint32_t( m_vertOffset.size() );     assert( m_vertOffset.size() < maxAllowedSize );
  m_uboData.m_vertOffset_capacity = uint32_t( m_vertOffset.capacity() ); assert( m_vertOffset.capacity() < maxAllowedSize );
  m_uboData.m_wavelengths_size     = uint32_t( m_wavelengths.size() );     assert( m_wavelengths.size() < maxAllowedSize );
  m_uboData.m_wavelengths_capacity = uint32_t( m_wavelengths.capacity() ); assert( m_wavelengths.capacity() < maxAllowedSize );
  a_pCopyEngine->UpdateBuffer(m_classDataBuffer, 0, &m_uboData, sizeof(m_uboData));
}

void Integrator_Generated::ReadPlainMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{
  a_pCopyEngine->ReadBuffer(m_classDataBuffer, 0, &m_uboData, sizeof(m_uboData));
  m_projInv = m_uboData.m_projInv;
  m_worldViewInv = m_uboData.m_worldViewInv;
  m_envColor = m_uboData.m_envColor;
  m_exposureMult = m_uboData.m_exposureMult;
  m_intergatorType = m_uboData.m_intergatorType;
  m_maxThreadId = m_uboData.m_maxThreadId;
  m_skipBounce = m_uboData.m_skipBounce;
  m_spectral_mode = m_uboData.m_spectral_mode;
  m_tileSize = m_uboData.m_tileSize;
  m_traceDepth = m_uboData.m_traceDepth;
  m_winHeight = m_uboData.m_winHeight;
  m_winWidth = m_uboData.m_winWidth;
  m_allRemapLists.resize(m_uboData.m_allRemapLists_size);
  m_allRemapListsOffsets.resize(m_uboData.m_allRemapListsOffsets_size);
  m_cie_x.resize(m_uboData.m_cie_x_size);
  m_cie_y.resize(m_uboData.m_cie_y_size);
  m_cie_z.resize(m_uboData.m_cie_z_size);
  m_instIdToLightInstId.resize(m_uboData.m_instIdToLightInstId_size);
  m_lights.resize(m_uboData.m_lights_size);
  m_matIdByPrimId.resize(m_uboData.m_matIdByPrimId_size);
  m_matIdOffsets.resize(m_uboData.m_matIdOffsets_size);
  m_materials.resize(m_uboData.m_materials_size);
  m_normMatrices.resize(m_uboData.m_normMatrices_size);
  m_packedXY.resize(m_uboData.m_packedXY_size);
  m_randomGens.resize(m_uboData.m_randomGens_size);
  m_remapInst.resize(m_uboData.m_remapInst_size);
  m_spec_offset_sz.resize(m_uboData.m_spec_offset_sz_size);
  m_spec_values.resize(m_uboData.m_spec_values_size);
  m_triIndices.resize(m_uboData.m_triIndices_size);
  m_vNorm4f.resize(m_uboData.m_vNorm4f_size);
  m_vTang4f.resize(m_uboData.m_vTang4f_size);
  m_vertOffset.resize(m_uboData.m_vertOffset_size);
  m_wavelengths.resize(m_uboData.m_wavelengths_size);
}

void Integrator_Generated::UpdateVectorMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{
  if(m_allRemapLists.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_allRemapListsBuffer, 0, m_allRemapLists.data(), m_allRemapLists.size()*sizeof(int) );
  if(m_allRemapListsOffsets.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_allRemapListsOffsetsBuffer, 0, m_allRemapListsOffsets.data(), m_allRemapListsOffsets.size()*sizeof(int) );
  if(m_cie_x.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_cie_xBuffer, 0, m_cie_x.data(), m_cie_x.size()*sizeof(float) );
  if(m_cie_y.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_cie_yBuffer, 0, m_cie_y.data(), m_cie_y.size()*sizeof(float) );
  if(m_cie_z.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_cie_zBuffer, 0, m_cie_z.data(), m_cie_z.size()*sizeof(float) );
  if(m_instIdToLightInstId.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_instIdToLightInstIdBuffer, 0, m_instIdToLightInstId.data(), m_instIdToLightInstId.size()*sizeof(unsigned int) );
  if(m_lights.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_lightsBuffer, 0, m_lights.data(), m_lights.size()*sizeof(struct LightSource) );
  if(m_matIdByPrimId.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_matIdByPrimIdBuffer, 0, m_matIdByPrimId.data(), m_matIdByPrimId.size()*sizeof(unsigned int) );
  if(m_matIdOffsets.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_matIdOffsetsBuffer, 0, m_matIdOffsets.data(), m_matIdOffsets.size()*sizeof(unsigned int) );
  if(m_materials.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_materialsBuffer, 0, m_materials.data(), m_materials.size()*sizeof(struct Material) );
  if(m_normMatrices.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_normMatricesBuffer, 0, m_normMatrices.data(), m_normMatrices.size()*sizeof(struct LiteMath::float4x4) );
  if(m_packedXY.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_packedXYBuffer, 0, m_packedXY.data(), m_packedXY.size()*sizeof(unsigned int) );
  if(m_randomGens.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_randomGensBuffer, 0, m_randomGens.data(), m_randomGens.size()*sizeof(struct RandomGenT) );
  if(m_remapInst.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_remapInstBuffer, 0, m_remapInst.data(), m_remapInst.size()*sizeof(int) );
  if(m_spec_offset_sz.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_spec_offset_szBuffer, 0, m_spec_offset_sz.data(), m_spec_offset_sz.size()*sizeof(struct LiteMath::uint2) );
  if(m_spec_values.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_spec_valuesBuffer, 0, m_spec_values.data(), m_spec_values.size()*sizeof(float) );
  if(m_triIndices.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_triIndicesBuffer, 0, m_triIndices.data(), m_triIndices.size()*sizeof(unsigned int) );
  if(m_vNorm4f.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_vNorm4fBuffer, 0, m_vNorm4f.data(), m_vNorm4f.size()*sizeof(struct LiteMath::float4) );
  if(m_vTang4f.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_vTang4fBuffer, 0, m_vTang4f.data(), m_vTang4f.size()*sizeof(struct LiteMath::float4) );
  if(m_vertOffset.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_vertOffsetBuffer, 0, m_vertOffset.data(), m_vertOffset.size()*sizeof(unsigned int) );
  if(m_wavelengths.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_wavelengthsBuffer, 0, m_wavelengths.data(), m_wavelengths.size()*sizeof(float) );
}

void Integrator_Generated::UpdateTextureMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{ 
  for(int i=0;i<m_vdata.m_texturesArrayTexture.size();i++)
    a_pCopyEngine->UpdateImage(m_vdata.m_texturesArrayTexture[i], m_textures[i]->data(), m_textures[i]->width(), m_textures[i]->height(), m_textures[i]->bpp(), VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL); 
  
  std::array<VkImageMemoryBarrier, 0> barriers;

  
  VkCommandBuffer cmdBuff       = a_pCopyEngine->CmdBuffer();
  VkQueue         transferQueue = a_pCopyEngine->TransferQueue();

  vkResetCommandBuffer(cmdBuff, 0);
  VkCommandBufferBeginInfo beginInfo = {};
  beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
  beginInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
  if (vkBeginCommandBuffer(cmdBuff, &beginInfo) != VK_SUCCESS)
    throw std::runtime_error("Integrator_Generated::UpdateTextureMembers: failed to begin command buffer!");
  vkCmdPipelineBarrier(cmdBuff,VK_PIPELINE_STAGE_TRANSFER_BIT,VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,0,0,nullptr,0,nullptr,uint32_t(barriers.size()),barriers.data());
  vkEndCommandBuffer(cmdBuff);
  
  vk_utils::executeCommandBufferNow(cmdBuff, transferQueue, device);
}

void Integrator_Generated::RayTraceMegaCmd(uint tid, float4* out_color)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(tid);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = tid;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;

  vkCmdPushConstants(m_currCmdBuffer, RayTraceMegaLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, RayTraceMegaPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void Integrator_Generated::CastSingleRayMegaCmd(uint tid, uint* out_color)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(tid);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = tid;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;

  vkCmdPushConstants(m_currCmdBuffer, CastSingleRayMegaLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, CastSingleRayMegaPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void Integrator_Generated::PackXYMegaCmd(uint tidX, uint tidY)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(tidX);
  uint32_t sizeY  = uint32_t(tidY);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = tidX;
  pcData.m_sizeY  = tidY;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;

  vkCmdPushConstants(m_currCmdBuffer, PackXYMegaLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, PackXYMegaPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void Integrator_Generated::PathTraceFromInputRaysMegaCmd(uint tid, const RayPosAndW* in_rayPosAndNear, const RayDirAndT* in_rayDirAndFar, float4* out_color)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(tid);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = tid;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;

  vkCmdPushConstants(m_currCmdBuffer, PathTraceFromInputRaysMegaLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, PathTraceFromInputRaysMegaPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void Integrator_Generated::PathTraceMegaCmd(uint tid, float4* out_color)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(tid);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = tid;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;

  vkCmdPushConstants(m_currCmdBuffer, PathTraceMegaLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, PathTraceMegaPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void Integrator_Generated::NaivePathTraceMegaCmd(uint tid, float4* out_color)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(tid);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = tid;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;

  vkCmdPushConstants(m_currCmdBuffer, NaivePathTraceMegaLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, NaivePathTraceMegaPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}


void Integrator_Generated::copyKernelFloatCmd(uint32_t length)
{
  uint32_t blockSizeX = MEMCPY_BLOCK_SIZE;

  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, copyKernelFloatPipeline);
  vkCmdPushConstants(m_currCmdBuffer, copyKernelFloatLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(uint32_t), &length);
  vkCmdDispatch(m_currCmdBuffer, (length + blockSizeX - 1) / blockSizeX, 1, 1);
}

VkBufferMemoryBarrier Integrator_Generated::BarrierForClearFlags(VkBuffer a_buffer)
{
  VkBufferMemoryBarrier bar = {};
  bar.sType               = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
  bar.pNext               = NULL;
  bar.srcAccessMask       = VK_ACCESS_TRANSFER_WRITE_BIT;
  bar.dstAccessMask       = VK_ACCESS_SHADER_READ_BIT;
  bar.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
  bar.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
  bar.buffer              = a_buffer;
  bar.offset              = 0;
  bar.size                = VK_WHOLE_SIZE;
  return bar;
}

VkBufferMemoryBarrier Integrator_Generated::BarrierForSingleBuffer(VkBuffer a_buffer)
{
  VkBufferMemoryBarrier bar = {};
  bar.sType               = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
  bar.pNext               = NULL;
  bar.srcAccessMask       = VK_ACCESS_SHADER_WRITE_BIT;
  bar.dstAccessMask       = VK_ACCESS_SHADER_READ_BIT;
  bar.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
  bar.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
  bar.buffer              = a_buffer;
  bar.offset              = 0;
  bar.size                = VK_WHOLE_SIZE;
  return bar;
}

void Integrator_Generated::BarriersForSeveralBuffers(VkBuffer* a_inBuffers, VkBufferMemoryBarrier* a_outBarriers, uint32_t a_buffersNum)
{
  for(uint32_t i=0; i<a_buffersNum;i++)
  {
    a_outBarriers[i].sType               = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
    a_outBarriers[i].pNext               = NULL;
    a_outBarriers[i].srcAccessMask       = VK_ACCESS_SHADER_WRITE_BIT;
    a_outBarriers[i].dstAccessMask       = VK_ACCESS_SHADER_READ_BIT;
    a_outBarriers[i].srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    a_outBarriers[i].dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    a_outBarriers[i].buffer              = a_inBuffers[i];
    a_outBarriers[i].offset              = 0;
    a_outBarriers[i].size                = VK_WHOLE_SIZE;
  }
}

void Integrator_Generated::RayTraceCmd(VkCommandBuffer a_commandBuffer, uint tid, float4* out_color)
{
  m_currCmdBuffer = a_commandBuffer;
  VkMemoryBarrier memoryBarrier = { VK_STRUCTURE_TYPE_MEMORY_BARRIER, nullptr, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT }; 
  vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, RayTraceMegaLayout, 0, 1, &m_allGeneratedDS[0], 0, nullptr);
  RayTraceMegaCmd(tid, out_color);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr); 
}

void Integrator_Generated::CastSingleRayCmd(VkCommandBuffer a_commandBuffer, uint tid, uint* out_color)
{
  m_currCmdBuffer = a_commandBuffer;
  VkMemoryBarrier memoryBarrier = { VK_STRUCTURE_TYPE_MEMORY_BARRIER, nullptr, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT }; 
  vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, CastSingleRayMegaLayout, 0, 1, &m_allGeneratedDS[1], 0, nullptr);
  CastSingleRayMegaCmd(tid, out_color);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr); 
}

void Integrator_Generated::PackXYCmd(VkCommandBuffer a_commandBuffer, uint tidX, uint tidY)
{
  m_currCmdBuffer = a_commandBuffer;
  VkMemoryBarrier memoryBarrier = { VK_STRUCTURE_TYPE_MEMORY_BARRIER, nullptr, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT }; 
  vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, PackXYMegaLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  PackXYMegaCmd(tidX, tidY);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr); 
}

void Integrator_Generated::PathTraceFromInputRaysCmd(VkCommandBuffer a_commandBuffer, uint tid, const RayPosAndW* in_rayPosAndNear, const RayDirAndT* in_rayDirAndFar, float4* out_color)
{
  m_currCmdBuffer = a_commandBuffer;
  VkMemoryBarrier memoryBarrier = { VK_STRUCTURE_TYPE_MEMORY_BARRIER, nullptr, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT }; 
  vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, PathTraceFromInputRaysMegaLayout, 0, 1, &m_allGeneratedDS[3], 0, nullptr);
  PathTraceFromInputRaysMegaCmd(tid, in_rayPosAndNear, in_rayDirAndFar, out_color);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr); 
}

void Integrator_Generated::PathTraceCmd(VkCommandBuffer a_commandBuffer, uint tid, float4* out_color)
{
  m_currCmdBuffer = a_commandBuffer;
  VkMemoryBarrier memoryBarrier = { VK_STRUCTURE_TYPE_MEMORY_BARRIER, nullptr, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT }; 
  vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, PathTraceMegaLayout, 0, 1, &m_allGeneratedDS[4], 0, nullptr);
  PathTraceMegaCmd(tid, out_color);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr); 
}

void Integrator_Generated::NaivePathTraceCmd(VkCommandBuffer a_commandBuffer, uint tid, float4* out_color)
{
  m_currCmdBuffer = a_commandBuffer;
  VkMemoryBarrier memoryBarrier = { VK_STRUCTURE_TYPE_MEMORY_BARRIER, nullptr, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT }; 
  vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, NaivePathTraceMegaLayout, 0, 1, &m_allGeneratedDS[5], 0, nullptr);
  NaivePathTraceMegaCmd(tid, out_color);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr); 
}



void Integrator_Generated::RayTraceBlock(uint tid, float4* out_color, uint32_t a_numPasses)
{
  // (1) get global Vulkan context objects
  //
  VkInstance       instance       = m_ctx.instance;
  VkPhysicalDevice physicalDevice = m_ctx.physicalDevice;
  VkDevice         device         = m_ctx.device;
  VkCommandPool    commandPool    = m_ctx.commandPool; 
  VkQueue          computeQueue   = m_ctx.computeQueue; 
  VkQueue          transferQueue  = m_ctx.transferQueue;
  auto             pCopyHelper    = m_ctx.pCopyHelper;
  auto             pAllocatorSpec = m_ctx.pAllocatorSpecial;

  // (2) create GPU objects
  //
  auto outFlags = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
  if(RayTrace_local.needToClearOutput)
    outFlags |= VK_IMAGE_USAGE_TRANSFER_DST_BIT;
  std::vector<VkBuffer> buffers;
  std::vector<VkImage>  images2;
  std::vector<vk_utils::VulkanImageMem*> images;
  auto beforeCreateObjects = std::chrono::high_resolution_clock::now();
  VkBuffer out_colorGPU = vk_utils::createBuffer(device, tid*sizeof(float4 ), outFlags);
  buffers.push_back(out_colorGPU);
  

  VkDeviceMemory buffersMem = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, buffers);
  VkDeviceMemory imagesMem  = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, std::vector<VkBuffer>(), images);
  
  vk_utils::MemAllocInfo tempMemoryAllocInfo;
  tempMemoryAllocInfo.memUsage = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT; // TODO, select depending on device and sample/application (???)
  if(buffers.size() != 0)
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, buffers);
  if(images.size() != 0)
  {
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, images2);
    for(auto imgMem : images)
    {
      VkImageViewCreateInfo imageView{};
      imageView.sType    = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
      imageView.viewType = VK_IMAGE_VIEW_TYPE_2D;
      imageView.image    = imgMem->image;
      imageView.format   = imgMem->format;
      imageView.subresourceRange = {};
      imageView.subresourceRange.aspectMask     = imgMem->aspectMask;
      imageView.subresourceRange.baseMipLevel   = 0;
      imageView.subresourceRange.levelCount     = imgMem->mipLvls;
      imageView.subresourceRange.baseArrayLayer = 0;
      imageView.subresourceRange.layerCount     = 1;
      VK_CHECK_RESULT(vkCreateImageView(device, &imageView, nullptr, &imgMem->view));
    }
  }
  
  auto afterCreateObjects = std::chrono::high_resolution_clock::now();
  m_exTimeRayTrace.msAPIOverhead = std::chrono::duration_cast<std::chrono::microseconds>(afterCreateObjects - beforeCreateObjects).count()/1000.f;
  
  auto afterCopy2 = std::chrono::high_resolution_clock::now(); // just declare it here, replace value later
  
  auto afterInitBuffers = std::chrono::high_resolution_clock::now();
  m_exTimeRayTrace.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(afterInitBuffers - afterCreateObjects).count()/1000.f;
  
  auto beforeSetInOut = std::chrono::high_resolution_clock::now();
  this->SetVulkanInOutFor_RayTrace(out_colorGPU, 0); 

  // (3) copy input data to GPU
  //
  auto beforeCopy = std::chrono::high_resolution_clock::now();
  m_exTimeRayTrace.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(beforeCopy - beforeSetInOut).count()/1000.f;
  auto afterCopy = std::chrono::high_resolution_clock::now();
  m_exTimeRayTrace.msCopyToGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy - beforeCopy).count()/1000.f;
  //  
  m_exTimeRayTrace.msExecuteOnGPU = 0;
  //// (3.1) clear all outputs if we are in RTV pattern
  //
  if(RayTrace_local.needToClearOutput)
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    vkCmdFillBuffer(commandBuffer, out_colorGPU, 0, VK_WHOLE_SIZE, 0); // zero output buffer out_colorGPU
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimeRayTrace.msExecuteOnGPU  += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }

  // (4) now execute algorithm on GPU
  //
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    RayTraceCmd(commandBuffer, tid, out_color);      
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    for(uint32_t pass = 0; pass < a_numPasses; pass++)
      vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimeRayTrace.msExecuteOnGPU += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }
  
  // (5) copy output data to CPU
  //
  auto beforeCopy2 = std::chrono::high_resolution_clock::now();
  pCopyHelper->ReadBuffer(out_colorGPU, 0, out_color, tid*sizeof(float4 ));
  this->ReadPlainMembers(pCopyHelper);
  afterCopy2 = std::chrono::high_resolution_clock::now();
  m_exTimeRayTrace.msCopyFromGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy2 - beforeCopy2).count()/1000.f;

  // (6) free resources 
  //
  //vkDestroyBuffer(device, out_colorGPU, nullptr);
  if(buffersMem != VK_NULL_HANDLE)
    vkFreeMemory(device, buffersMem, nullptr);
  if(imagesMem != VK_NULL_HANDLE)
    vkFreeMemory(device, imagesMem, nullptr);
  
  m_exTimeRayTrace.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - afterCopy2).count()/1000.f;
}

void Integrator_Generated::CastSingleRayBlock(uint tid, uint* out_color, uint32_t a_numPasses)
{
  // (1) get global Vulkan context objects
  //
  VkInstance       instance       = m_ctx.instance;
  VkPhysicalDevice physicalDevice = m_ctx.physicalDevice;
  VkDevice         device         = m_ctx.device;
  VkCommandPool    commandPool    = m_ctx.commandPool; 
  VkQueue          computeQueue   = m_ctx.computeQueue; 
  VkQueue          transferQueue  = m_ctx.transferQueue;
  auto             pCopyHelper    = m_ctx.pCopyHelper;
  auto             pAllocatorSpec = m_ctx.pAllocatorSpecial;

  // (2) create GPU objects
  //
  auto outFlags = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
  if(CastSingleRay_local.needToClearOutput)
    outFlags |= VK_IMAGE_USAGE_TRANSFER_DST_BIT;
  std::vector<VkBuffer> buffers;
  std::vector<VkImage>  images2;
  std::vector<vk_utils::VulkanImageMem*> images;
  auto beforeCreateObjects = std::chrono::high_resolution_clock::now();
  VkBuffer out_colorGPU = vk_utils::createBuffer(device, tid*sizeof(uint ), outFlags);
  buffers.push_back(out_colorGPU);
  

  VkDeviceMemory buffersMem = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, buffers);
  VkDeviceMemory imagesMem  = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, std::vector<VkBuffer>(), images);
  
  vk_utils::MemAllocInfo tempMemoryAllocInfo;
  tempMemoryAllocInfo.memUsage = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT; // TODO, select depending on device and sample/application (???)
  if(buffers.size() != 0)
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, buffers);
  if(images.size() != 0)
  {
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, images2);
    for(auto imgMem : images)
    {
      VkImageViewCreateInfo imageView{};
      imageView.sType    = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
      imageView.viewType = VK_IMAGE_VIEW_TYPE_2D;
      imageView.image    = imgMem->image;
      imageView.format   = imgMem->format;
      imageView.subresourceRange = {};
      imageView.subresourceRange.aspectMask     = imgMem->aspectMask;
      imageView.subresourceRange.baseMipLevel   = 0;
      imageView.subresourceRange.levelCount     = imgMem->mipLvls;
      imageView.subresourceRange.baseArrayLayer = 0;
      imageView.subresourceRange.layerCount     = 1;
      VK_CHECK_RESULT(vkCreateImageView(device, &imageView, nullptr, &imgMem->view));
    }
  }
  
  auto afterCreateObjects = std::chrono::high_resolution_clock::now();
  m_exTimeCastSingleRay.msAPIOverhead = std::chrono::duration_cast<std::chrono::microseconds>(afterCreateObjects - beforeCreateObjects).count()/1000.f;
  
  auto afterCopy2 = std::chrono::high_resolution_clock::now(); // just declare it here, replace value later
  
  auto afterInitBuffers = std::chrono::high_resolution_clock::now();
  m_exTimeCastSingleRay.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(afterInitBuffers - afterCreateObjects).count()/1000.f;
  
  auto beforeSetInOut = std::chrono::high_resolution_clock::now();
  this->SetVulkanInOutFor_CastSingleRay(out_colorGPU, 0); 

  // (3) copy input data to GPU
  //
  auto beforeCopy = std::chrono::high_resolution_clock::now();
  m_exTimeCastSingleRay.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(beforeCopy - beforeSetInOut).count()/1000.f;
  auto afterCopy = std::chrono::high_resolution_clock::now();
  m_exTimeCastSingleRay.msCopyToGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy - beforeCopy).count()/1000.f;
  //  
  m_exTimeCastSingleRay.msExecuteOnGPU = 0;
  //// (3.1) clear all outputs if we are in RTV pattern
  //
  if(CastSingleRay_local.needToClearOutput)
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    vkCmdFillBuffer(commandBuffer, out_colorGPU, 0, VK_WHOLE_SIZE, 0); // zero output buffer out_colorGPU
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimeCastSingleRay.msExecuteOnGPU  += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }

  // (4) now execute algorithm on GPU
  //
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    CastSingleRayCmd(commandBuffer, tid, out_color);      
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    for(uint32_t pass = 0; pass < a_numPasses; pass++)
      vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimeCastSingleRay.msExecuteOnGPU += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }
  
  // (5) copy output data to CPU
  //
  auto beforeCopy2 = std::chrono::high_resolution_clock::now();
  pCopyHelper->ReadBuffer(out_colorGPU, 0, out_color, tid*sizeof(uint ));
  this->ReadPlainMembers(pCopyHelper);
  afterCopy2 = std::chrono::high_resolution_clock::now();
  m_exTimeCastSingleRay.msCopyFromGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy2 - beforeCopy2).count()/1000.f;

  // (6) free resources 
  //
  //vkDestroyBuffer(device, out_colorGPU, nullptr);
  if(buffersMem != VK_NULL_HANDLE)
    vkFreeMemory(device, buffersMem, nullptr);
  if(imagesMem != VK_NULL_HANDLE)
    vkFreeMemory(device, imagesMem, nullptr);
  
  m_exTimeCastSingleRay.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - afterCopy2).count()/1000.f;
}

void Integrator_Generated::PackXYBlock(uint tidX, uint tidY, uint32_t a_numPasses)
{
  // (1) get global Vulkan context objects
  //
  VkInstance       instance       = m_ctx.instance;
  VkPhysicalDevice physicalDevice = m_ctx.physicalDevice;
  VkDevice         device         = m_ctx.device;
  VkCommandPool    commandPool    = m_ctx.commandPool; 
  VkQueue          computeQueue   = m_ctx.computeQueue; 
  VkQueue          transferQueue  = m_ctx.transferQueue;
  auto             pCopyHelper    = m_ctx.pCopyHelper;
  auto             pAllocatorSpec = m_ctx.pAllocatorSpecial;

  // (2) create GPU objects
  //
  auto outFlags = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
  if(PackXY_local.needToClearOutput)
    outFlags |= VK_IMAGE_USAGE_TRANSFER_DST_BIT;
  std::vector<VkBuffer> buffers;
  std::vector<VkImage>  images2;
  std::vector<vk_utils::VulkanImageMem*> images;
  auto beforeCreateObjects = std::chrono::high_resolution_clock::now();
  

  VkDeviceMemory buffersMem = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, buffers);
  VkDeviceMemory imagesMem  = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, std::vector<VkBuffer>(), images);
  
  vk_utils::MemAllocInfo tempMemoryAllocInfo;
  tempMemoryAllocInfo.memUsage = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT; // TODO, select depending on device and sample/application (???)
  if(buffers.size() != 0)
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, buffers);
  if(images.size() != 0)
  {
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, images2);
    for(auto imgMem : images)
    {
      VkImageViewCreateInfo imageView{};
      imageView.sType    = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
      imageView.viewType = VK_IMAGE_VIEW_TYPE_2D;
      imageView.image    = imgMem->image;
      imageView.format   = imgMem->format;
      imageView.subresourceRange = {};
      imageView.subresourceRange.aspectMask     = imgMem->aspectMask;
      imageView.subresourceRange.baseMipLevel   = 0;
      imageView.subresourceRange.levelCount     = imgMem->mipLvls;
      imageView.subresourceRange.baseArrayLayer = 0;
      imageView.subresourceRange.layerCount     = 1;
      VK_CHECK_RESULT(vkCreateImageView(device, &imageView, nullptr, &imgMem->view));
    }
  }
  
  auto afterCreateObjects = std::chrono::high_resolution_clock::now();
  m_exTimePackXY.msAPIOverhead = std::chrono::duration_cast<std::chrono::microseconds>(afterCreateObjects - beforeCreateObjects).count()/1000.f;
  
  auto afterCopy2 = std::chrono::high_resolution_clock::now(); // just declare it here, replace value later
  
  auto afterInitBuffers = std::chrono::high_resolution_clock::now();
  m_exTimePackXY.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(afterInitBuffers - afterCreateObjects).count()/1000.f;
  
  auto beforeSetInOut = std::chrono::high_resolution_clock::now();
  this->SetVulkanInOutFor_PackXY(); 

  // (3) copy input data to GPU
  //
  auto beforeCopy = std::chrono::high_resolution_clock::now();
  m_exTimePackXY.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(beforeCopy - beforeSetInOut).count()/1000.f;
  auto afterCopy = std::chrono::high_resolution_clock::now();
  m_exTimePackXY.msCopyToGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy - beforeCopy).count()/1000.f;
  //  
  m_exTimePackXY.msExecuteOnGPU = 0;
  //// (3.1) clear all outputs if we are in RTV pattern
  //
  if(PackXY_local.needToClearOutput)
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimePackXY.msExecuteOnGPU  += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }

  // (4) now execute algorithm on GPU
  //
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    PackXYCmd(commandBuffer, tidX, tidY);      
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    for(uint32_t pass = 0; pass < a_numPasses; pass++)
      vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimePackXY.msExecuteOnGPU += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }
  
  // (5) copy output data to CPU
  //
  auto beforeCopy2 = std::chrono::high_resolution_clock::now();
  this->ReadPlainMembers(pCopyHelper);
  afterCopy2 = std::chrono::high_resolution_clock::now();
  m_exTimePackXY.msCopyFromGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy2 - beforeCopy2).count()/1000.f;

  // (6) free resources 
  //
  if(buffersMem != VK_NULL_HANDLE)
    vkFreeMemory(device, buffersMem, nullptr);
  if(imagesMem != VK_NULL_HANDLE)
    vkFreeMemory(device, imagesMem, nullptr);
  
  m_exTimePackXY.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - afterCopy2).count()/1000.f;
}

void Integrator_Generated::PathTraceFromInputRaysBlock(uint tid, const RayPosAndW* in_rayPosAndNear, const RayDirAndT* in_rayDirAndFar, float4* out_color, uint32_t a_numPasses)
{
  // (1) get global Vulkan context objects
  //
  VkInstance       instance       = m_ctx.instance;
  VkPhysicalDevice physicalDevice = m_ctx.physicalDevice;
  VkDevice         device         = m_ctx.device;
  VkCommandPool    commandPool    = m_ctx.commandPool; 
  VkQueue          computeQueue   = m_ctx.computeQueue; 
  VkQueue          transferQueue  = m_ctx.transferQueue;
  auto             pCopyHelper    = m_ctx.pCopyHelper;
  auto             pAllocatorSpec = m_ctx.pAllocatorSpecial;

  // (2) create GPU objects
  //
  auto outFlags = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
  if(PathTraceFromInputRays_local.needToClearOutput)
    outFlags |= VK_IMAGE_USAGE_TRANSFER_DST_BIT;
  std::vector<VkBuffer> buffers;
  std::vector<VkImage>  images2;
  std::vector<vk_utils::VulkanImageMem*> images;
  auto beforeCreateObjects = std::chrono::high_resolution_clock::now();
  VkBuffer in_rayPosAndNearGPU = vk_utils::createBuffer(device, tid*sizeof(const RayPosAndW ), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  buffers.push_back(in_rayPosAndNearGPU);
  VkBuffer in_rayDirAndFarGPU = vk_utils::createBuffer(device, tid*sizeof(const RayDirAndT ), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  buffers.push_back(in_rayDirAndFarGPU);
  VkBuffer out_colorGPU = vk_utils::createBuffer(device, tid*sizeof(float4 ), outFlags);
  buffers.push_back(out_colorGPU);
  

  VkDeviceMemory buffersMem = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, buffers);
  VkDeviceMemory imagesMem  = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, std::vector<VkBuffer>(), images);
  
  vk_utils::MemAllocInfo tempMemoryAllocInfo;
  tempMemoryAllocInfo.memUsage = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT; // TODO, select depending on device and sample/application (???)
  if(buffers.size() != 0)
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, buffers);
  if(images.size() != 0)
  {
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, images2);
    for(auto imgMem : images)
    {
      VkImageViewCreateInfo imageView{};
      imageView.sType    = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
      imageView.viewType = VK_IMAGE_VIEW_TYPE_2D;
      imageView.image    = imgMem->image;
      imageView.format   = imgMem->format;
      imageView.subresourceRange = {};
      imageView.subresourceRange.aspectMask     = imgMem->aspectMask;
      imageView.subresourceRange.baseMipLevel   = 0;
      imageView.subresourceRange.levelCount     = imgMem->mipLvls;
      imageView.subresourceRange.baseArrayLayer = 0;
      imageView.subresourceRange.layerCount     = 1;
      VK_CHECK_RESULT(vkCreateImageView(device, &imageView, nullptr, &imgMem->view));
    }
  }
  
  auto afterCreateObjects = std::chrono::high_resolution_clock::now();
  m_exTimePathTraceFromInputRays.msAPIOverhead = std::chrono::duration_cast<std::chrono::microseconds>(afterCreateObjects - beforeCreateObjects).count()/1000.f;
  
  auto afterCopy2 = std::chrono::high_resolution_clock::now(); // just declare it here, replace value later
  
  auto afterInitBuffers = std::chrono::high_resolution_clock::now();
  m_exTimePathTraceFromInputRays.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(afterInitBuffers - afterCreateObjects).count()/1000.f;
  
  auto beforeSetInOut = std::chrono::high_resolution_clock::now();
  this->SetVulkanInOutFor_PathTraceFromInputRays(in_rayPosAndNearGPU, 0, in_rayDirAndFarGPU, 0, out_colorGPU, 0); 

  // (3) copy input data to GPU
  //
  auto beforeCopy = std::chrono::high_resolution_clock::now();
  m_exTimePathTraceFromInputRays.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(beforeCopy - beforeSetInOut).count()/1000.f;
  pCopyHelper->UpdateBuffer(in_rayPosAndNearGPU, 0, in_rayPosAndNear, tid*sizeof(const RayPosAndW )); 
  pCopyHelper->UpdateBuffer(in_rayDirAndFarGPU, 0, in_rayDirAndFar, tid*sizeof(const RayDirAndT )); 
  auto afterCopy = std::chrono::high_resolution_clock::now();
  m_exTimePathTraceFromInputRays.msCopyToGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy - beforeCopy).count()/1000.f;
  //  
  m_exTimePathTraceFromInputRays.msExecuteOnGPU = 0;
  //// (3.1) clear all outputs if we are in RTV pattern
  //
  if(PathTraceFromInputRays_local.needToClearOutput)
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    vkCmdFillBuffer(commandBuffer, out_colorGPU, 0, VK_WHOLE_SIZE, 0); // zero output buffer out_colorGPU
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimePathTraceFromInputRays.msExecuteOnGPU  += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }

  // (4) now execute algorithm on GPU
  //
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    PathTraceFromInputRaysCmd(commandBuffer, tid, in_rayPosAndNear, in_rayDirAndFar, out_color);      
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    for(uint32_t pass = 0; pass < a_numPasses; pass++)
      vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimePathTraceFromInputRays.msExecuteOnGPU += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }
  
  // (5) copy output data to CPU
  //
  auto beforeCopy2 = std::chrono::high_resolution_clock::now();
  pCopyHelper->ReadBuffer(out_colorGPU, 0, out_color, tid*sizeof(float4 ));
  this->ReadPlainMembers(pCopyHelper);
  afterCopy2 = std::chrono::high_resolution_clock::now();
  m_exTimePathTraceFromInputRays.msCopyFromGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy2 - beforeCopy2).count()/1000.f;

  // (6) free resources 
  //
  //vkDestroyBuffer(device, in_rayPosAndNearGPU, nullptr);
  //vkDestroyBuffer(device, in_rayDirAndFarGPU, nullptr);
  //vkDestroyBuffer(device, out_colorGPU, nullptr);
  if(buffersMem != VK_NULL_HANDLE)
    vkFreeMemory(device, buffersMem, nullptr);
  if(imagesMem != VK_NULL_HANDLE)
    vkFreeMemory(device, imagesMem, nullptr);
  
  m_exTimePathTraceFromInputRays.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - afterCopy2).count()/1000.f;
}

void Integrator_Generated::PathTraceBlock(uint tid, float4* out_color, uint32_t a_numPasses)
{
  // (1) get global Vulkan context objects
  //
  VkInstance       instance       = m_ctx.instance;
  VkPhysicalDevice physicalDevice = m_ctx.physicalDevice;
  VkDevice         device         = m_ctx.device;
  VkCommandPool    commandPool    = m_ctx.commandPool; 
  VkQueue          computeQueue   = m_ctx.computeQueue; 
  VkQueue          transferQueue  = m_ctx.transferQueue;
  auto             pCopyHelper    = m_ctx.pCopyHelper;
  auto             pAllocatorSpec = m_ctx.pAllocatorSpecial;

  // (2) create GPU objects
  //
  auto outFlags = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
  if(PathTrace_local.needToClearOutput)
    outFlags |= VK_IMAGE_USAGE_TRANSFER_DST_BIT;
  std::vector<VkBuffer> buffers;
  std::vector<VkImage>  images2;
  std::vector<vk_utils::VulkanImageMem*> images;
  auto beforeCreateObjects = std::chrono::high_resolution_clock::now();
  VkBuffer out_colorGPU = vk_utils::createBuffer(device, tid*sizeof(float4 ), outFlags);
  buffers.push_back(out_colorGPU);
  

  VkDeviceMemory buffersMem = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, buffers);
  VkDeviceMemory imagesMem  = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, std::vector<VkBuffer>(), images);
  
  vk_utils::MemAllocInfo tempMemoryAllocInfo;
  tempMemoryAllocInfo.memUsage = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT; // TODO, select depending on device and sample/application (???)
  if(buffers.size() != 0)
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, buffers);
  if(images.size() != 0)
  {
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, images2);
    for(auto imgMem : images)
    {
      VkImageViewCreateInfo imageView{};
      imageView.sType    = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
      imageView.viewType = VK_IMAGE_VIEW_TYPE_2D;
      imageView.image    = imgMem->image;
      imageView.format   = imgMem->format;
      imageView.subresourceRange = {};
      imageView.subresourceRange.aspectMask     = imgMem->aspectMask;
      imageView.subresourceRange.baseMipLevel   = 0;
      imageView.subresourceRange.levelCount     = imgMem->mipLvls;
      imageView.subresourceRange.baseArrayLayer = 0;
      imageView.subresourceRange.layerCount     = 1;
      VK_CHECK_RESULT(vkCreateImageView(device, &imageView, nullptr, &imgMem->view));
    }
  }
  
  auto afterCreateObjects = std::chrono::high_resolution_clock::now();
  m_exTimePathTrace.msAPIOverhead = std::chrono::duration_cast<std::chrono::microseconds>(afterCreateObjects - beforeCreateObjects).count()/1000.f;
  
  auto afterCopy2 = std::chrono::high_resolution_clock::now(); // just declare it here, replace value later
  
  auto afterInitBuffers = std::chrono::high_resolution_clock::now();
  m_exTimePathTrace.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(afterInitBuffers - afterCreateObjects).count()/1000.f;
  
  auto beforeSetInOut = std::chrono::high_resolution_clock::now();
  this->SetVulkanInOutFor_PathTrace(out_colorGPU, 0); 

  // (3) copy input data to GPU
  //
  auto beforeCopy = std::chrono::high_resolution_clock::now();
  m_exTimePathTrace.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(beforeCopy - beforeSetInOut).count()/1000.f;
  auto afterCopy = std::chrono::high_resolution_clock::now();
  m_exTimePathTrace.msCopyToGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy - beforeCopy).count()/1000.f;
  //  
  m_exTimePathTrace.msExecuteOnGPU = 0;
  //// (3.1) clear all outputs if we are in RTV pattern
  //
  if(PathTrace_local.needToClearOutput)
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    vkCmdFillBuffer(commandBuffer, out_colorGPU, 0, VK_WHOLE_SIZE, 0); // zero output buffer out_colorGPU
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimePathTrace.msExecuteOnGPU  += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }

  // (4) now execute algorithm on GPU
  //
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    PathTraceCmd(commandBuffer, tid, out_color);      
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    for(uint32_t pass = 0; pass < a_numPasses; pass++)
      vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimePathTrace.msExecuteOnGPU += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }
  
  // (5) copy output data to CPU
  //
  auto beforeCopy2 = std::chrono::high_resolution_clock::now();
  pCopyHelper->ReadBuffer(out_colorGPU, 0, out_color, tid*sizeof(float4 ));
  this->ReadPlainMembers(pCopyHelper);
  afterCopy2 = std::chrono::high_resolution_clock::now();
  m_exTimePathTrace.msCopyFromGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy2 - beforeCopy2).count()/1000.f;

  // (6) free resources 
  //
  //vkDestroyBuffer(device, out_colorGPU, nullptr);
  if(buffersMem != VK_NULL_HANDLE)
    vkFreeMemory(device, buffersMem, nullptr);
  if(imagesMem != VK_NULL_HANDLE)
    vkFreeMemory(device, imagesMem, nullptr);
  
  m_exTimePathTrace.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - afterCopy2).count()/1000.f;
}

void Integrator_Generated::NaivePathTraceBlock(uint tid, float4* out_color, uint32_t a_numPasses)
{
  // (1) get global Vulkan context objects
  //
  VkInstance       instance       = m_ctx.instance;
  VkPhysicalDevice physicalDevice = m_ctx.physicalDevice;
  VkDevice         device         = m_ctx.device;
  VkCommandPool    commandPool    = m_ctx.commandPool; 
  VkQueue          computeQueue   = m_ctx.computeQueue; 
  VkQueue          transferQueue  = m_ctx.transferQueue;
  auto             pCopyHelper    = m_ctx.pCopyHelper;
  auto             pAllocatorSpec = m_ctx.pAllocatorSpecial;

  // (2) create GPU objects
  //
  auto outFlags = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
  if(NaivePathTrace_local.needToClearOutput)
    outFlags |= VK_IMAGE_USAGE_TRANSFER_DST_BIT;
  std::vector<VkBuffer> buffers;
  std::vector<VkImage>  images2;
  std::vector<vk_utils::VulkanImageMem*> images;
  auto beforeCreateObjects = std::chrono::high_resolution_clock::now();
  VkBuffer out_colorGPU = vk_utils::createBuffer(device, tid*sizeof(float4 ), outFlags);
  buffers.push_back(out_colorGPU);
  

  VkDeviceMemory buffersMem = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, buffers);
  VkDeviceMemory imagesMem  = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, std::vector<VkBuffer>(), images);
  
  vk_utils::MemAllocInfo tempMemoryAllocInfo;
  tempMemoryAllocInfo.memUsage = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT; // TODO, select depending on device and sample/application (???)
  if(buffers.size() != 0)
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, buffers);
  if(images.size() != 0)
  {
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, images2);
    for(auto imgMem : images)
    {
      VkImageViewCreateInfo imageView{};
      imageView.sType    = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
      imageView.viewType = VK_IMAGE_VIEW_TYPE_2D;
      imageView.image    = imgMem->image;
      imageView.format   = imgMem->format;
      imageView.subresourceRange = {};
      imageView.subresourceRange.aspectMask     = imgMem->aspectMask;
      imageView.subresourceRange.baseMipLevel   = 0;
      imageView.subresourceRange.levelCount     = imgMem->mipLvls;
      imageView.subresourceRange.baseArrayLayer = 0;
      imageView.subresourceRange.layerCount     = 1;
      VK_CHECK_RESULT(vkCreateImageView(device, &imageView, nullptr, &imgMem->view));
    }
  }
  
  auto afterCreateObjects = std::chrono::high_resolution_clock::now();
  m_exTimeNaivePathTrace.msAPIOverhead = std::chrono::duration_cast<std::chrono::microseconds>(afterCreateObjects - beforeCreateObjects).count()/1000.f;
  
  auto afterCopy2 = std::chrono::high_resolution_clock::now(); // just declare it here, replace value later
  
  auto afterInitBuffers = std::chrono::high_resolution_clock::now();
  m_exTimeNaivePathTrace.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(afterInitBuffers - afterCreateObjects).count()/1000.f;
  
  auto beforeSetInOut = std::chrono::high_resolution_clock::now();
  this->SetVulkanInOutFor_NaivePathTrace(out_colorGPU, 0); 

  // (3) copy input data to GPU
  //
  auto beforeCopy = std::chrono::high_resolution_clock::now();
  m_exTimeNaivePathTrace.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(beforeCopy - beforeSetInOut).count()/1000.f;
  auto afterCopy = std::chrono::high_resolution_clock::now();
  m_exTimeNaivePathTrace.msCopyToGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy - beforeCopy).count()/1000.f;
  //  
  m_exTimeNaivePathTrace.msExecuteOnGPU = 0;
  //// (3.1) clear all outputs if we are in RTV pattern
  //
  if(NaivePathTrace_local.needToClearOutput)
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    vkCmdFillBuffer(commandBuffer, out_colorGPU, 0, VK_WHOLE_SIZE, 0); // zero output buffer out_colorGPU
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimeNaivePathTrace.msExecuteOnGPU  += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }

  // (4) now execute algorithm on GPU
  //
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    NaivePathTraceCmd(commandBuffer, tid, out_color);      
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    for(uint32_t pass = 0; pass < a_numPasses; pass++)
      vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimeNaivePathTrace.msExecuteOnGPU += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }
  
  // (5) copy output data to CPU
  //
  auto beforeCopy2 = std::chrono::high_resolution_clock::now();
  pCopyHelper->ReadBuffer(out_colorGPU, 0, out_color, tid*sizeof(float4 ));
  this->ReadPlainMembers(pCopyHelper);
  afterCopy2 = std::chrono::high_resolution_clock::now();
  m_exTimeNaivePathTrace.msCopyFromGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy2 - beforeCopy2).count()/1000.f;

  // (6) free resources 
  //
  //vkDestroyBuffer(device, out_colorGPU, nullptr);
  if(buffersMem != VK_NULL_HANDLE)
    vkFreeMemory(device, buffersMem, nullptr);
  if(imagesMem != VK_NULL_HANDLE)
    vkFreeMemory(device, imagesMem, nullptr);
  
  m_exTimeNaivePathTrace.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - afterCopy2).count()/1000.f;
}

void Integrator_Generated::GetExecutionTime(const char* a_funcName, float a_out[4])
{
  vk_utils::ExecTime res = {};
  if(std::string(a_funcName) == "RayTrace" || std::string(a_funcName) == "RayTraceBlock")
    res = m_exTimeRayTrace;
  if(std::string(a_funcName) == "CastSingleRay" || std::string(a_funcName) == "CastSingleRayBlock")
    res = m_exTimeCastSingleRay;
  if(std::string(a_funcName) == "PackXY" || std::string(a_funcName) == "PackXYBlock")
    res = m_exTimePackXY;
  if(std::string(a_funcName) == "PathTraceFromInputRays" || std::string(a_funcName) == "PathTraceFromInputRaysBlock")
    res = m_exTimePathTraceFromInputRays;
  if(std::string(a_funcName) == "PathTrace" || std::string(a_funcName) == "PathTraceBlock")
    res = m_exTimePathTrace;
  if(std::string(a_funcName) == "NaivePathTrace" || std::string(a_funcName) == "NaivePathTraceBlock")
    res = m_exTimeNaivePathTrace;
  a_out[0] = res.msExecuteOnGPU;
  a_out[1] = res.msCopyToGPU;
  a_out[2] = res.msCopyFromGPU;
  a_out[3] = res.msAPIOverhead;             
}

