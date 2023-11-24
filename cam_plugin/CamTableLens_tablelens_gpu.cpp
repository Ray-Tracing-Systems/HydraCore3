#include <vector>
#include <memory>
#include <limits>
#include <cassert>
#include <chrono>

#include "vk_copy.h"
#include "vk_context.h"
#include "vk_images.h"

#include "CamTableLens_tablelens_gpu.h"
#include "include/CamTableLens_tablelens_gpu_ubo.h"



std::shared_ptr<CamTableLens> CreateCamTableLens_TABLELENS_GPU(vk_utils::VulkanContext a_ctx, size_t a_maxThreadsGenerated) 
{ 
  auto pObj = std::make_shared<CamTableLens_TABLELENS_GPU>(); 
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

constexpr uint32_t KGEN_REDUCTION_LAST_STEP    = 16;

void CamTableLens_TABLELENS_GPU::InitVulkanObjects(VkDevice a_device, VkPhysicalDevice a_physicalDevice, size_t a_maxThreadsCount) 
{
  physicalDevice = a_physicalDevice;
  device         = a_device;
  m_allCreatedPipelineLayouts.reserve(256);
  m_allCreatedPipelines.reserve(256);
  InitHelpers();
  InitBuffers(a_maxThreadsCount, true);
  InitKernels(".spv");
  AllocateAllDescriptorSets();

}

void CamTableLens_TABLELENS_GPU::UpdatePlainMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{
  const size_t maxAllowedSize = std::numeric_limits<uint32_t>::max();
  m_uboData.m_physSize = m_physSize;
  m_uboData.m_height = m_height;
  m_uboData.m_spectral_mode = m_spectral_mode;
  m_uboData.m_width = m_width;
  m_uboData.lines_size     = uint32_t( lines.size() );     assert( lines.size() < maxAllowedSize );
  m_uboData.lines_capacity = uint32_t( lines.capacity() ); assert( lines.capacity() < maxAllowedSize );
  m_uboData.m_cie_x_size     = uint32_t( m_cie_x.size() );     assert( m_cie_x.size() < maxAllowedSize );
  m_uboData.m_cie_x_capacity = uint32_t( m_cie_x.capacity() ); assert( m_cie_x.capacity() < maxAllowedSize );
  m_uboData.m_cie_y_size     = uint32_t( m_cie_y.size() );     assert( m_cie_y.size() < maxAllowedSize );
  m_uboData.m_cie_y_capacity = uint32_t( m_cie_y.capacity() ); assert( m_cie_y.capacity() < maxAllowedSize );
  m_uboData.m_cie_z_size     = uint32_t( m_cie_z.size() );     assert( m_cie_z.size() < maxAllowedSize );
  m_uboData.m_cie_z_capacity = uint32_t( m_cie_z.capacity() ); assert( m_cie_z.capacity() < maxAllowedSize );
  m_uboData.m_randomGens_size     = uint32_t( m_randomGens.size() );     assert( m_randomGens.size() < maxAllowedSize );
  m_uboData.m_randomGens_capacity = uint32_t( m_randomGens.capacity() ); assert( m_randomGens.capacity() < maxAllowedSize );
  m_uboData.m_storedCos4_size     = uint32_t( m_storedCos4.size() );     assert( m_storedCos4.size() < maxAllowedSize );
  m_uboData.m_storedCos4_capacity = uint32_t( m_storedCos4.capacity() ); assert( m_storedCos4.capacity() < maxAllowedSize );
  m_uboData.m_storedWaves_size     = uint32_t( m_storedWaves.size() );     assert( m_storedWaves.size() < maxAllowedSize );
  m_uboData.m_storedWaves_capacity = uint32_t( m_storedWaves.capacity() ); assert( m_storedWaves.capacity() < maxAllowedSize );
  a_pCopyEngine->UpdateBuffer(m_classDataBuffer, 0, &m_uboData, sizeof(m_uboData));
}

void CamTableLens_TABLELENS_GPU::ReadPlainMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{
  a_pCopyEngine->ReadBuffer(m_classDataBuffer, 0, &m_uboData, sizeof(m_uboData));
  m_physSize = m_uboData.m_physSize;
  m_height = m_uboData.m_height;
  m_spectral_mode = m_uboData.m_spectral_mode;
  m_width = m_uboData.m_width;
  lines.resize(m_uboData.lines_size);
  m_cie_x.resize(m_uboData.m_cie_x_size);
  m_cie_y.resize(m_uboData.m_cie_y_size);
  m_cie_z.resize(m_uboData.m_cie_z_size);
  m_randomGens.resize(m_uboData.m_randomGens_size);
  m_storedCos4.resize(m_uboData.m_storedCos4_size);
  m_storedWaves.resize(m_uboData.m_storedWaves_size);
}

void CamTableLens_TABLELENS_GPU::UpdateVectorMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{
  if(lines.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.linesBuffer, 0, lines.data(), lines.size()*sizeof(struct LensElementInterface) );
  if(m_cie_x.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_cie_xBuffer, 0, m_cie_x.data(), m_cie_x.size()*sizeof(float) );
  if(m_cie_y.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_cie_yBuffer, 0, m_cie_y.data(), m_cie_y.size()*sizeof(float) );
  if(m_cie_z.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_cie_zBuffer, 0, m_cie_z.data(), m_cie_z.size()*sizeof(float) );
  if(m_randomGens.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_randomGensBuffer, 0, m_randomGens.data(), m_randomGens.size()*sizeof(struct RandomGenT) );
  if(m_storedCos4.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_storedCos4Buffer, 0, m_storedCos4.data(), m_storedCos4.size()*sizeof(float) );
  if(m_storedWaves.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.m_storedWavesBuffer, 0, m_storedWaves.data(), m_storedWaves.size()*sizeof(struct LiteMath::uint2) );
}

void CamTableLens_TABLELENS_GPU::UpdateTextureMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{ 
}

void CamTableLens_TABLELENS_GPU::MakeEyeRayCmd(int in_blockSize, RayPart1* out_rayPosAndNear4f, RayPart2* out_rayDirAndFar4f, int subPassId)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    int m_subPassId; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(in_blockSize);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = in_blockSize;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_subPassId = subPassId; 

  vkCmdPushConstants(m_currCmdBuffer, MakeEyeRayLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, MakeEyeRayPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void CamTableLens_TABLELENS_GPU::ContribSampleCmd(int in_blockSize, const float* in_color, float* out_color, int subPassId)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    int m_subPassId; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(in_blockSize);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = in_blockSize;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_subPassId = subPassId; 

  vkCmdPushConstants(m_currCmdBuffer, ContribSampleLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, ContribSamplePipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}


void CamTableLens_TABLELENS_GPU::copyKernelFloatCmd(uint32_t length)
{
  uint32_t blockSizeX = MEMCPY_BLOCK_SIZE;

  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, copyKernelFloatPipeline);
  vkCmdPushConstants(m_currCmdBuffer, copyKernelFloatLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(uint32_t), &length);
  vkCmdDispatch(m_currCmdBuffer, (length + blockSizeX - 1) / blockSizeX, 1, 1);
}

VkBufferMemoryBarrier CamTableLens_TABLELENS_GPU::BarrierForClearFlags(VkBuffer a_buffer)
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

VkBufferMemoryBarrier CamTableLens_TABLELENS_GPU::BarrierForSingleBuffer(VkBuffer a_buffer)
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

void CamTableLens_TABLELENS_GPU::BarriersForSeveralBuffers(VkBuffer* a_inBuffers, VkBufferMemoryBarrier* a_outBarriers, uint32_t a_buffersNum)
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

void CamTableLens_TABLELENS_GPU::AddSamplesContributionBlockCmd(VkCommandBuffer a_commandBuffer, float* out_color4f, const float* colors4f, uint32_t in_blockSize, uint32_t a_width, uint32_t a_height, int subPassId)
{
  m_currCmdBuffer = a_commandBuffer;
  VkMemoryBarrier memoryBarrier = { VK_STRUCTURE_TYPE_MEMORY_BARRIER, nullptr, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT }; 
    vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, ContribSampleLayout, 0, 1, &m_allGeneratedDS[0], 0, nullptr);
  ContribSampleCmd(int(in_blockSize), colors4f, out_color4f, subPassId);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr); 

}

void CamTableLens_TABLELENS_GPU::MakeRaysBlockCmd(VkCommandBuffer a_commandBuffer, RayPart1* out_rayPosAndNear4f, RayPart2* out_rayDirAndFar4f, uint32_t in_blockSize, int subPassId)
{
  m_currCmdBuffer = a_commandBuffer;
  VkMemoryBarrier memoryBarrier = { VK_STRUCTURE_TYPE_MEMORY_BARRIER, nullptr, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT }; 
    vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, MakeEyeRayLayout, 0, 1, &m_allGeneratedDS[1], 0, nullptr);
  MakeEyeRayCmd(int(in_blockSize), out_rayPosAndNear4f, out_rayDirAndFar4f, subPassId);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);

}



void CamTableLens_TABLELENS_GPU::AddSamplesContributionBlock(float* out_color4f, const float* colors4f, uint32_t in_blockSize, uint32_t a_width, uint32_t a_height, int subPassId)
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
  if(AddSamplesContributionBlock_local.needToClearOutput)
    outFlags |= VK_IMAGE_USAGE_TRANSFER_DST_BIT;
  std::vector<VkBuffer> buffers;
  std::vector<VkImage>  images2;
  std::vector<vk_utils::VulkanImageMem*> images;
  auto beforeCreateObjects = std::chrono::high_resolution_clock::now();
  VkBuffer colors4fGPU = vk_utils::createBuffer(device, in_blockSize*4*sizeof(const float ), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  buffers.push_back(colors4fGPU);
  VkBuffer out_color4fGPU = vk_utils::createBuffer(device, a_width*a_height*4*sizeof(float ), outFlags);
  buffers.push_back(out_color4fGPU);
  

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
  m_exTimeAddSamplesContributionBlock.msAPIOverhead = std::chrono::duration_cast<std::chrono::microseconds>(afterCreateObjects - beforeCreateObjects).count()/1000.f;
  
  auto afterCopy2 = std::chrono::high_resolution_clock::now(); // just declare it here, replace value later
  
  auto afterInitBuffers = std::chrono::high_resolution_clock::now();
  m_exTimeAddSamplesContributionBlock.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(afterInitBuffers - afterCreateObjects).count()/1000.f;
  
  auto beforeSetInOut = std::chrono::high_resolution_clock::now();
  this->SetVulkanInOutFor_AddSamplesContributionBlock(out_color4fGPU, 0, colors4fGPU, 0, 0); 

  // (3) copy input data to GPU
  //
  auto beforeCopy = std::chrono::high_resolution_clock::now();
  m_exTimeAddSamplesContributionBlock.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(beforeCopy - beforeSetInOut).count()/1000.f;
  pCopyHelper->UpdateBuffer(colors4fGPU, 0, colors4f, in_blockSize*4*sizeof(const float )); 
  auto afterCopy = std::chrono::high_resolution_clock::now();
  m_exTimeAddSamplesContributionBlock.msCopyToGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy - beforeCopy).count()/1000.f;
  //  
  m_exTimeAddSamplesContributionBlock.msExecuteOnGPU = 0;

  // (4) now execute algorithm on GPU
  //
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    AddSamplesContributionBlockCmd(commandBuffer, out_color4f, colors4f, in_blockSize, a_width, a_height, subPassId);      
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimeAddSamplesContributionBlock.msExecuteOnGPU += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }
  
  // (5) copy output data to CPU
  //
  auto beforeCopy2 = std::chrono::high_resolution_clock::now();
  pCopyHelper->ReadBuffer(out_color4fGPU, 0, out_color4f, a_width*a_height*4*sizeof(float ));
  this->ReadPlainMembers(pCopyHelper);
  afterCopy2 = std::chrono::high_resolution_clock::now();
  m_exTimeAddSamplesContributionBlock.msCopyFromGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy2 - beforeCopy2).count()/1000.f;

  // (6) free resources 
  //
  //vkDestroyBuffer(device, colors4fGPU, nullptr);
  //vkDestroyBuffer(device, out_color4fGPU, nullptr);
  if(buffersMem != VK_NULL_HANDLE)
    vkFreeMemory(device, buffersMem, nullptr);
  if(imagesMem != VK_NULL_HANDLE)
    vkFreeMemory(device, imagesMem, nullptr);
  
  m_exTimeAddSamplesContributionBlock.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - afterCopy2).count()/1000.f;
}

void CamTableLens_TABLELENS_GPU::MakeRaysBlock(RayPart1* out_rayPosAndNear4f, RayPart2* out_rayDirAndFar4f, uint32_t in_blockSize, int subPassId)
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
  if(MakeRaysBlock_local.needToClearOutput)
    outFlags |= VK_IMAGE_USAGE_TRANSFER_DST_BIT;
  std::vector<VkBuffer> buffers;
  std::vector<VkImage>  images2;
  std::vector<vk_utils::VulkanImageMem*> images;
  auto beforeCreateObjects = std::chrono::high_resolution_clock::now();
  VkBuffer out_rayPosAndNear4fGPU = vk_utils::createBuffer(device, in_blockSize*sizeof(RayPart1 ), outFlags);
  buffers.push_back(out_rayPosAndNear4fGPU);
  VkBuffer out_rayDirAndFar4fGPU = vk_utils::createBuffer(device, in_blockSize*sizeof(RayPart2 ), outFlags);
  buffers.push_back(out_rayDirAndFar4fGPU);
  

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
  m_exTimeMakeRaysBlock.msAPIOverhead = std::chrono::duration_cast<std::chrono::microseconds>(afterCreateObjects - beforeCreateObjects).count()/1000.f;
  
  auto afterCopy2 = std::chrono::high_resolution_clock::now(); // just declare it here, replace value later
  
  auto afterInitBuffers = std::chrono::high_resolution_clock::now();
  m_exTimeMakeRaysBlock.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(afterInitBuffers - afterCreateObjects).count()/1000.f;
  
  auto beforeSetInOut = std::chrono::high_resolution_clock::now();
  this->SetVulkanInOutFor_MakeRaysBlock(out_rayPosAndNear4fGPU, 0, out_rayDirAndFar4fGPU, 0, 0); 

  // (3) copy input data to GPU
  //
  auto beforeCopy = std::chrono::high_resolution_clock::now();
  m_exTimeMakeRaysBlock.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(beforeCopy - beforeSetInOut).count()/1000.f;
  auto afterCopy = std::chrono::high_resolution_clock::now();
  m_exTimeMakeRaysBlock.msCopyToGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy - beforeCopy).count()/1000.f;
  //  
  m_exTimeMakeRaysBlock.msExecuteOnGPU = 0;

  // (4) now execute algorithm on GPU
  //
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    MakeRaysBlockCmd(commandBuffer, out_rayPosAndNear4f, out_rayDirAndFar4f, in_blockSize, subPassId);      
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimeMakeRaysBlock.msExecuteOnGPU += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }
  
  // (5) copy output data to CPU
  //
  auto beforeCopy2 = std::chrono::high_resolution_clock::now();
  pCopyHelper->ReadBuffer(out_rayPosAndNear4fGPU, 0, out_rayPosAndNear4f, in_blockSize*sizeof(RayPart1 ));
  pCopyHelper->ReadBuffer(out_rayDirAndFar4fGPU, 0, out_rayDirAndFar4f, in_blockSize*sizeof(RayPart2 ));
  this->ReadPlainMembers(pCopyHelper);
  afterCopy2 = std::chrono::high_resolution_clock::now();
  m_exTimeMakeRaysBlock.msCopyFromGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy2 - beforeCopy2).count()/1000.f;

  // (6) free resources 
  //
  //vkDestroyBuffer(device, out_rayPosAndNear4fGPU, nullptr);
  //vkDestroyBuffer(device, out_rayDirAndFar4fGPU, nullptr);
  if(buffersMem != VK_NULL_HANDLE)
    vkFreeMemory(device, buffersMem, nullptr);
  if(imagesMem != VK_NULL_HANDLE)
    vkFreeMemory(device, imagesMem, nullptr);
  
  m_exTimeMakeRaysBlock.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - afterCopy2).count()/1000.f;
}

void CamTableLens_TABLELENS_GPU::GetExecutionTime(const char* a_funcName, float a_out[4])
{
  vk_utils::ExecTime res = {};
  if(std::string(a_funcName) == "AddSamplesContributionBlock" || std::string(a_funcName) == "AddSamplesContributionBlockBlock")
    res = m_exTimeAddSamplesContributionBlock;
  if(std::string(a_funcName) == "MakeRaysBlock" || std::string(a_funcName) == "MakeRaysBlockBlock")
    res = m_exTimeMakeRaysBlock;
  a_out[0] = res.msExecuteOnGPU;
  a_out[1] = res.msCopyToGPU;
  a_out[2] = res.msCopyFromGPU;
  a_out[3] = res.msAPIOverhead;             
}

