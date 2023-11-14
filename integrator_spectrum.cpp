#include "integrator_pt.h"


float4 Integrator::SampleMatColorParamSpectrum(uint32_t matId, float4 a_wavelengths, uint32_t paramId, uint32_t paramSpecId)
{  
  float4 res = m_materials[matId].colors[paramId];
  if(a_wavelengths.M[0] == 0.0f)
    return res;

  const uint specId = as_uint(m_materials[matId].data[paramSpecId]);

  if(specId < 0xFFFFFFFF)
  {
    auto &[offset, size] = m_spec_offset_sz[specId];
    res = SampleSpectrum(m_wavelengths.data() + offset, m_spec_values.data() + offset, a_wavelengths, size);
  }

  return res;
}

float4 Integrator::SampleMatParamSpectrum(uint32_t matId, float4 a_wavelengths, uint32_t paramId, uint32_t paramSpecId)
{  
  float4 res = float4(m_materials[matId].data[paramId]);
  if(a_wavelengths.M[0] == 0.0f)
    return res;

  const uint specId = as_uint(m_materials[matId].data[paramSpecId]);

  if(specId < 0xFFFFFFFF)
  {
    auto &[offset, size] = m_spec_offset_sz[specId];
    res = SampleSpectrum(m_wavelengths.data() + offset, m_spec_values.data() + offset, a_wavelengths, size);
  }

  return res;
}
