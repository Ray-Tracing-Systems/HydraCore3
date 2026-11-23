#pragma once

#include "integrator_pt_generated.h"
#include "pipeline_profile.h"

// Keep shader names at the generated class's public extension point. Generated
// sources and vkutils do not need to be patched after each translation.
class Integrator_Profiled : public Integrator_Generated
{
public:
  using Integrator_Generated::Integrator_Generated;

protected:
  void MakeComputePipelineAndLayout(const char* path, const char* entry,
      const VkSpecializationInfo* spec, VkDescriptorSetLayout ds,
      VkPipelineLayout* layout, VkPipeline* pipeline) override
  {
    hydra_profile::SetShaderPath(path);
    Integrator_Generated::MakeComputePipelineAndLayout(path, entry, spec, ds, layout, pipeline);
    hydra_profile::SetShaderPath(nullptr);
  }

  void MakeComputePipelineOnly(const char* path, const char* entry,
      const VkSpecializationInfo* spec, VkDescriptorSetLayout ds,
      VkPipelineLayout layout, VkPipeline* pipeline) override
  {
    hydra_profile::SetShaderPath(path);
    Integrator_Generated::MakeComputePipelineOnly(path, entry, spec, ds, layout, pipeline);
    hydra_profile::SetShaderPath(nullptr);
  }
};
