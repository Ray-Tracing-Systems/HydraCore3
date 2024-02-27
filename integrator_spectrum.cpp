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

float4 Integrator::SampleMatParamSpectrumTexture(uint32_t matId, float4 a_wavelengths, uint32_t paramId, uint32_t paramSpecId, float2 texCoords)
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
    
    if(size == 0)
    {
      const uint2 tex_data  = m_spec_tex_offset_sz[specId];
      const uint tex_offset = tex_data.x;
      const uint tex_size   = tex_data.y;
      if(tex_size > 0)
      {
        for(int i = 0; i < 4; ++i)
        {
          if (a_wavelengths[i] < m_spec_tex_ids_wavelengths[tex_offset].y ||
              a_wavelengths[i] > m_spec_tex_ids_wavelengths[tex_offset + tex_size - 1].y)
          {
            res[i] = 0.0f;
            continue;
          }

          int o = BinarySearchU2(m_spec_tex_ids_wavelengths.data() + tex_offset, tex_size, a_wavelengths[i]);

          int texID1 = m_spec_tex_ids_wavelengths[tex_offset + o].x;
          int texID2 = m_spec_tex_ids_wavelengths[tex_offset + o + 1].x;

          const float2 texCoordT = mulRows2x4(m_materials[matId].row0[0], m_materials[matId].row1[0], texCoords);
          const float4 texColor1 = m_textures[texID1]->sample(texCoordT);
          const float4 texColor2 = m_textures[texID2]->sample(texCoordT);
  
          float t = (a_wavelengths[i] - m_spec_tex_ids_wavelengths[tex_offset + o].y) / 
                    (m_spec_tex_ids_wavelengths[tex_offset + o + 1].y - m_spec_tex_ids_wavelengths[tex_offset + o].y );
                    
          float4 outColor = lerp(texColor1, texColor2, t);

          if(std::isinf(outColor.x) || std::isnan(outColor.x) || outColor.x < 0)
            std::cout << t << " " << outColor.x << std::endl;
          
          res[i] = outColor.x;
        }
      }
    }
    else
    {
      //res = SampleSpectrum(m_wavelengths.data() + offset, m_spec_values.data() + offset, a_wavelengths, size);
      res = SampleUniformSpectrum(m_spec_values.data() + offset, a_wavelengths, size);
    }
  }

  return res;
}
