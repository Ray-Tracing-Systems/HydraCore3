#pragma once

#include "vk_context.h"
#include <string>

// Optional, single-device profiler for this application's single-threaded Vulkan host.
// No hooks are installed unless -profile-json is passed to hydra.
namespace hydra_profile
{
  void Configure(vk_utils::VulkanDeviceFeatures& features, unsigned deviceIndex);
  void Install(const vk_utils::VulkanContext& context);
  void BeginRender();
  void Write(const std::string& path);
  void Uninstall();
  void SetShaderPath(const char* path);
}
