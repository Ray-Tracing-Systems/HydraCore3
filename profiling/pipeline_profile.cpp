#include "pipeline_profile.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <vector>

namespace hydra_profile
{
namespace
{
  bool supported = false;
  bool recording = false;
  std::string shaderPath;
  vk_utils::VulkanContext context;
  VkPhysicalDevicePipelineExecutablePropertiesFeaturesKHR executableFeature{};
  PFN_vkCreateComputePipelines originalCreate = nullptr;
  PFN_vkCmdBindPipeline originalBind = nullptr;
  PFN_vkCmdDispatch originalDispatch = nullptr;
  PFN_vkCmdDispatchIndirect originalDispatchIndirect = nullptr;

  struct Pipeline
  {
    VkPipeline handle;
    std::string path;
    uint64_t recordedDispatches = 0;
  };
  std::vector<Pipeline> pipelines;
  std::map<VkCommandBuffer, VkPipeline> bound;

  void Check(VkResult result, const char* operation)
  {
    if(result != VK_SUCCESS)
      throw std::runtime_error(std::string(operation) + ": VkResult " + std::to_string(result));
  }

  // Escape arbitrary driver descriptions, including embedded newlines, as JSON.
  void String(std::ostream& out, const std::string& value)
  {
    out << '"';
    for(unsigned char c : value)
    {
      if(c == '"' || c == '\\') out << '\\' << c;
      else if(c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << unsigned(c) << std::dec;
      else out << c;
    }
    out << '"';
  }

  VKAPI_ATTR VkResult VKAPI_CALL Create(VkDevice device, VkPipelineCache cache,
      uint32_t count, const VkComputePipelineCreateInfo* input,
      const VkAllocationCallbacks* allocator, VkPipeline* output)
  {
    std::vector<VkComputePipelineCreateInfo> infos(input, input + count);
    if(supported)
      for(auto& info : infos)
        info.flags |= VK_PIPELINE_CREATE_CAPTURE_STATISTICS_BIT_KHR;
    VkResult result = originalCreate(device, cache, count, infos.data(), allocator, output);
    for(uint32_t i = 0; i < count; ++i)
      if(output[i] != VK_NULL_HANDLE)
        pipelines.push_back({output[i], shaderPath.empty() ? input[i].stage.pName : shaderPath, 0});
    return result;
  }

  VKAPI_ATTR void VKAPI_CALL Bind(VkCommandBuffer cmd, VkPipelineBindPoint point, VkPipeline pipeline)
  {
    if(point == VK_PIPELINE_BIND_POINT_COMPUTE)
      bound[cmd] = pipeline;
    originalBind(cmd, point, pipeline);
  }

  void Mark(VkCommandBuffer cmd)
  {
    if(!recording) return;
    auto it = bound.find(cmd);
    if(it == bound.end()) return;
    for(auto& pipeline : pipelines)
      if(pipeline.handle == it->second)
      {
        ++pipeline.recordedDispatches;
        break;
      }
  }

  VKAPI_ATTR void VKAPI_CALL Dispatch(VkCommandBuffer cmd, uint32_t x, uint32_t y, uint32_t z)
  {
    if(x && y && z) Mark(cmd);
    originalDispatch(cmd, x, y, z);
  }

  VKAPI_ATTR void VKAPI_CALL DispatchIndirect(VkCommandBuffer cmd, VkBuffer buffer, VkDeviceSize offset)
  {
    Mark(cmd);
    originalDispatchIndirect(cmd, buffer, offset);
  }

  void Executables(std::ostream& out, VkPipeline pipeline)
  {
    if(!supported) { out << "[]"; return; }
    VkPipelineInfoKHR info{VK_STRUCTURE_TYPE_PIPELINE_INFO_KHR};
    info.pipeline = pipeline;
    uint32_t count = 0;
    Check(vkGetPipelineExecutablePropertiesKHR(context.device, &info, &count, nullptr), "executable count");
    std::vector<VkPipelineExecutablePropertiesKHR> props;
    VkResult result;
    do
    {
      props.assign(count, {VK_STRUCTURE_TYPE_PIPELINE_EXECUTABLE_PROPERTIES_KHR});
      result = vkGetPipelineExecutablePropertiesKHR(context.device, &info, &count, props.data());
      if(result == VK_INCOMPLETE)
        Check(vkGetPipelineExecutablePropertiesKHR(context.device, &info, &count, nullptr), "executable count");
    } while(result == VK_INCOMPLETE);
    Check(result, "executable properties");
    props.resize(count);
    out << '[';
    for(uint32_t i = 0; i < count; ++i)
    {
      if(i) out << ',';
      out << "{\"name\":"; String(out, props[i].name);
      out << ",\"description\":"; String(out, props[i].description);
      out << ",\"subgroup_size\":" << props[i].subgroupSize << ",\"statistics\":[";
      VkPipelineExecutableInfoKHR exec{VK_STRUCTURE_TYPE_PIPELINE_EXECUTABLE_INFO_KHR};
      exec.pipeline = pipeline;
      exec.executableIndex = i;
      uint32_t n = 0;
      Check(vkGetPipelineExecutableStatisticsKHR(context.device, &exec, &n, nullptr), "statistics count");
      std::vector<VkPipelineExecutableStatisticKHR> stats;
      do
      {
        stats.assign(n, {VK_STRUCTURE_TYPE_PIPELINE_EXECUTABLE_STATISTIC_KHR});
        result = vkGetPipelineExecutableStatisticsKHR(context.device, &exec, &n, stats.data());
        if(result == VK_INCOMPLETE)
          Check(vkGetPipelineExecutableStatisticsKHR(context.device, &exec, &n, nullptr), "statistics count");
      } while(result == VK_INCOMPLETE);
      Check(result, "executable statistics");
      for(uint32_t j = 0; j < n; ++j)
      {
        if(j) out << ',';
        out << "{\"name\":"; String(out, stats[j].name);
        out << ",\"description\":"; String(out, stats[j].description);
        out << ",\"format\":" << stats[j].format << ",\"value\":";
        switch(stats[j].format)
        {
          case VK_PIPELINE_EXECUTABLE_STATISTIC_FORMAT_BOOL32_KHR: out << (stats[j].value.b32 ? "true" : "false"); break;
          case VK_PIPELINE_EXECUTABLE_STATISTIC_FORMAT_INT64_KHR: out << stats[j].value.i64; break;
          case VK_PIPELINE_EXECUTABLE_STATISTIC_FORMAT_UINT64_KHR: out << stats[j].value.u64; break;
          case VK_PIPELINE_EXECUTABLE_STATISTIC_FORMAT_FLOAT64_KHR: out << std::setprecision(17) << stats[j].value.f64; break;
          default: out << "null";
        }
        out << '}';
      }
      out << "]}";
    }
    out << ']';
  }
}

void Configure(vk_utils::VulkanDeviceFeatures& features, unsigned deviceIndex)
{
  // Probe without modifying vkutils or constructing a device with unsupported extensions.
  Check(volkInitialize(), "volkInitialize");
  VkApplicationInfo app{VK_STRUCTURE_TYPE_APPLICATION_INFO};
  app.apiVersion = VK_API_VERSION_1_1;
  VkInstanceCreateInfo ci{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
  ci.pApplicationInfo = &app;
  VkInstance instance = VK_NULL_HANDLE;
  Check(vkCreateInstance(&ci, nullptr, &instance), "profiling probe instance");
  volkLoadInstance(instance);
  uint32_t n = 0;
  Check(vkEnumeratePhysicalDevices(instance, &n, nullptr), "device count");
  std::vector<VkPhysicalDevice> devices(n);
  Check(vkEnumeratePhysicalDevices(instance, &n, devices.data()), "devices");
  if(deviceIndex >= n)
  {
    vkDestroyInstance(instance, nullptr);
    throw std::runtime_error("-gpu_id is outside the Vulkan device list");
  }
  VkPhysicalDevice device = devices[deviceIndex];
  Check(vkEnumerateDeviceExtensionProperties(device, nullptr, &n, nullptr), "extension count");
  std::vector<VkExtensionProperties> extensions(n);
  Check(vkEnumerateDeviceExtensionProperties(device, nullptr, &n, extensions.data()), "extensions");
  supported = std::any_of(extensions.begin(), extensions.end(), [](const auto& e) {
    return std::strcmp(e.extensionName, VK_KHR_PIPELINE_EXECUTABLE_PROPERTIES_EXTENSION_NAME) == 0;
  });
  executableFeature = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PIPELINE_EXECUTABLE_PROPERTIES_FEATURES_KHR};
  if(supported)
  {
    VkPhysicalDeviceFeatures2 query{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2};
    query.pNext = &executableFeature;
    vkGetPhysicalDeviceFeatures2(device, &query);
    supported = executableFeature.pipelineExecutableInfo == VK_TRUE;
  }
  vkDestroyInstance(instance, nullptr);
  if(supported)
  {
    features.extensionNames.push_back(VK_KHR_PIPELINE_EXECUTABLE_PROPERTIES_EXTENSION_NAME);
    executableFeature.pNext = features.features2.pNext;
    features.features2.pNext = &executableFeature;
  }
}

void Install(const vk_utils::VulkanContext& ctx)
{
  context = ctx;
  pipelines.clear();
  bound.clear();
  recording = false;
  if(supported && (!vkGetPipelineExecutablePropertiesKHR || !vkGetPipelineExecutableStatisticsKHR))
    throw std::runtime_error("pipeline executable entry points were not loaded");
  originalCreate = vkCreateComputePipelines; vkCreateComputePipelines = Create;
  originalBind = vkCmdBindPipeline; vkCmdBindPipeline = Bind;
  originalDispatch = vkCmdDispatch; vkCmdDispatch = Dispatch;
  originalDispatchIndirect = vkCmdDispatchIndirect; vkCmdDispatchIndirect = DispatchIndirect;
}

void SetShaderPath(const char* path) { shaderPath = path ? path : ""; }

void BeginRender()
{
  for(auto& pipeline : pipelines) pipeline.recordedDispatches = 0;
  recording = true;
}

void Write(const std::string& path)
{
  recording = false;
  Check(vkDeviceWaitIdle(context.device), "profiling wait");
  VkPhysicalDeviceDriverProperties driver{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_DRIVER_PROPERTIES};
  VkPhysicalDeviceProperties2 props{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2};
  props.pNext = &driver;
  vkGetPhysicalDeviceProperties2(context.physicalDevice, &props);
  std::ofstream out(path);
  if(!out) throw std::runtime_error("cannot write profile: " + path);
  out << "{\"schema_version\":1,\"source\":\"VK_KHR_pipeline_executable_properties\",\"supported\":"
      << (supported ? "true" : "false") << ",\"device\":{\"name\":";
  String(out, props.properties.deviceName);
  out << ",\"vendor_id\":" << props.properties.vendorID << ",\"device_id\":" << props.properties.deviceID
      << ",\"api_version\":" << props.properties.apiVersion << ",\"driver_version\":" << props.properties.driverVersion
      << ",\"driver_name\":"; String(out, driver.driverName);
  out << ",\"driver_info\":"; String(out, driver.driverInfo);
  out << "},\"kernels\":[";
  bool first = true;
  for(const auto& pipeline : pipelines)
  {
    if(!pipeline.recordedDispatches) continue;
    if(!first) out << ',';
    first = false;
    out << "{\"shader\":"; String(out, pipeline.path);
    out << ",\"recorded_dispatches\":" << pipeline.recordedDispatches << ",\"executables\":";
    Executables(out, pipeline.handle);
    out << '}';
  }
  out << "]}\n";
  out.close();
  if(!out) throw std::runtime_error("failed writing profile: " + path);
}

void Uninstall()
{
  vkCreateComputePipelines = originalCreate;
  vkCmdBindPipeline = originalBind;
  vkCmdDispatch = originalDispatch;
  vkCmdDispatchIndirect = originalDispatchIndirect;
  recording = false;
  context = vk_utils::VulkanContext{};
}
}
