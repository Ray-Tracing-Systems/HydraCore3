#include "spherical_integrator_pt.h"

Integrator::EyeRayData SphericalIntegrator::SampleCameraRay(RandomGen* pGen, uint tid)
{
  const uint XY = m_packedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;

  const float4 pixelOffsets = GetRandomNumbersLens(tid, pGen);

  const float fx = float(x) + pixelOffsets.x;
  const float fy = float(y) + pixelOffsets.y;

  float xCoordNormalized = (fx + float(m_winStartX))/float(m_fbWidth);
  const float yCoordNormalized = (fy + float(m_winStartY))/float(m_fbHeight);

  float3 rayDir = float3(0, -1, 0);
  float3 rayPos = float3(0,-10000000,0);
  if (m_mode == Mode::Hemisphere || m_mode == Mode::Sphere)
  {
    bool isBackside = false;
    if (m_mode == Mode::Sphere)
    {
      if (xCoordNormalized > 0.5f)
      {
        xCoordNormalized -= 0.5f;
        isBackside = true;
      }
      xCoordNormalized *= 2.0f;
    }

    float xWeight = xCoordNormalized * 2.0f - 1.0f;
    float yWeight = yCoordNormalized * 2.0f - 1.0f;

    float3 forward = to_float3(float4(0, 0, isBackside ? 1.0f : -1.0f, 0));
    float3 right   = to_float3(float4(1, 0, 0, 0));
    float3 up      = to_float3(float4(0, 1, 0, 0));
    const bool validPixel = xWeight * xWeight + yWeight * yWeight < 1.0f;

    if (validPixel)
    {
      const float zWeight = std::sqrt(1.0f - xWeight * xWeight - yWeight * yWeight);
      rayDir = normalize(forward * zWeight + right * xWeight + up * yWeight);
      rayPos = float3(0, 0, 0);
    }
  }
  else if (m_mode == Mode::PanoramaHemisphere || m_mode == Mode::PanoramaSphere)
  {
    const float theta = (1 - yCoordNormalized) * M_PI;
    const float phi = (m_mode == Mode::PanoramaHemisphere ? (xCoordNormalized + 1) : xCoordNormalized * 2.0f) * M_PI;
    float x_coord = std::sin(theta) * std::cos(phi);
    float y_coord = std::sin(theta) * std::sin(phi);
    float z_coord = std::cos(theta);
    rayDir = normalize(float3(x_coord, y_coord, z_coord));
    rayPos = float3(0, 0, 0);
  }

  EyeRayData res;
  {
    res.rayPos = rayPos;
    res.rayDir = rayDir;
    res.x      = x;
    res.y      = y;
    res.timeSam = 0.0f;
    res.waveSam = 1.0f;
    if(m_normMatrices2.size() != 0)
      res.timeSam = GetRandomNumbersTime(tid, pGen);
    if(KSPEC_SPECTRAL_RENDERING !=0 && m_spectral_mode != 0)
      res.waveSam = GetRandomNumbersSpec(tid, pGen);
    res.cosTheta = 1.0f;
  }

  RecordPixelRndIfNeeded(pixelOffsets, float2(res.waveSam, res.timeSam));

  return res;
}