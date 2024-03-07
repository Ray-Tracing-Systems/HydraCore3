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

float3 Integrator::SpectralCamRespoceToRGB(float4 specSamples, float4 waves, uint32_t rayFlags)
{
  float3 rgb = to_float3(specSamples);

  if(m_camResponseSpectrumId[0] < 0)
  {
    const float3 xyz = SpectrumToXYZ(specSamples, waves, LAMBDA_MIN, LAMBDA_MAX, m_cie_x.data(), m_cie_y.data(), m_cie_z.data(), terminateWavelngths(rayFlags));
    rgb = XYZToRGB(xyz);
  }
  else
  {
    float4 responceX, responceY, responceZ;
    {
      int specId = m_camResponseSpectrumId[0];
      if(specId >= 0)
      {
        const uint2 data  = m_spec_offset_sz[specId];
        const uint offset = data.x;
        const uint size   = data.y;
        responceX = SampleUniformSpectrum(m_spec_values.data() + offset, waves, size);
      }
      else
        responceX = float4(1,1,1,1);
      specId = m_camResponseSpectrumId[1];
      if(specId >= 0)
      {
        const uint2 data  = m_spec_offset_sz[specId];
        const uint offset = data.x;
        const uint size   = data.y;
        responceY = SampleUniformSpectrum(m_spec_values.data() + offset, waves, size);
      }
      else
        responceY = responceX;
      specId = m_camResponseSpectrumId[2];
      if(specId >= 0)
      {
        const uint2 data  = m_spec_offset_sz[specId];
        const uint offset = data.x;
        const uint size   = data.y;
        responceZ = SampleUniformSpectrum(m_spec_values.data() + offset, waves, size);
      }
      else
        responceZ = responceY;
    }
    float3 xyz = float3(0,0,0);
    for (uint32_t i = 0; i < SPECTRUM_SAMPLE_SZ; ++i) {
      xyz.x += specSamples[i]*responceX[i];
      xyz.y += specSamples[i]*responceY[i];
      xyz.z += specSamples[i]*responceZ[i]; 
    } 
    if(m_camResponseType == CAM_RESPONCE_XYZ)
      rgb = XYZToRGB(xyz);
    else
      rgb = xyz;
  }
  
  return rgb;
}
