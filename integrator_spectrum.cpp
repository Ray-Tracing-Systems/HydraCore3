#include "integrator_pt.h"


float4 Integrator::SampleMatColorParamSpectrum(uint32_t matId, float4 a_wavelengths, uint32_t paramId, uint32_t paramSpecId)
{  
  float4 res = m_materials[matId].colors[paramId];
  if(a_wavelengths[0] == 0.0f)
    return res;

  const uint specId = m_materials[matId].spdid[paramSpecId];
  if(specId < 0xFFFFFFFF)
  {
    const uint2 data  = m_spec_offset_sz[specId];
    const uint offset = data.x;
    const uint size   = data.y;
    //res = SampleSpectrum(m_wavelengths.data() + offset, m_spec_values.data() + offset, a_wavelengths, size);
    res = SampleUniformSpectrum(m_spec_values.data() + offset, a_wavelengths, size);
  }

  return res;
}

float4 Integrator::SampleMatParamSpectrum(uint32_t matId, float4 a_wavelengths, uint32_t paramId, uint32_t paramSpecId)
{  
  float4 res = float4(m_materials[matId].data[paramId]);
  if(a_wavelengths[0] == 0.0f)
    return res;

  const uint specId = m_materials[matId].spdid[paramSpecId];
  if(specId < 0xFFFFFFFF)
  {
    const uint2 data  = m_spec_offset_sz[specId];
    const uint offset = data.x;
    const uint size   = data.y;
    //res = SampleSpectrum(m_wavelengths.data() + offset, m_spec_values.data() + offset, a_wavelengths, size);
    res = SampleUniformSpectrum(m_spec_values.data() + offset, a_wavelengths, size);
  }

  return res;
}

float4 Integrator::SampleFilmsSpectrum(uint32_t matId, float4 a_wavelengths, uint32_t paramId, uint32_t paramSpecId, uint32_t layer)
{  
  float4 res = float4(m_films_eta_k_vec[as_uint(m_materials[matId].data[paramId]) + layer]);
  if(a_wavelengths[0] == 0.0f)
    return res;

  //const uint specId = m_films_eta_id_vec[as_uint(m_materials[a_materialId].data[FILM_ETA_SPECID_OFFSET]) + layer];
  const uint specId = m_films_spec_id_vec[as_uint(m_materials[matId].data[paramSpecId]) + layer];

  if(specId < 0xFFFFFFFF)
  {
    const uint2 data  = m_spec_offset_sz[specId];
    const uint offset = data.x;
    const uint size   = data.y;
    res = SampleUniformSpectrum(m_spec_values.data() + offset, a_wavelengths, size);
  }

  return res;
}