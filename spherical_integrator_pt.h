#pragma once

#include "integrator_pt.h"

class SphericalIntegrator : public Integrator
{
public:
  enum class Mode
  {
    Hemisphere,
    Sphere,
    PanoramaHemisphere,
    PanoramaSphere,
  };

  SphericalIntegrator(int a_maxThreads = 1, std::vector<uint32_t> a_features = {}, Mode mode = Mode::Hemisphere) :
    Integrator(a_maxThreads, a_features),
    m_mode(mode)
    {}
  EyeRayData SampleCameraRay(RandomGen* pGen, uint tid) override;

  uint32_t getOptimialHeight() const
  {
    return ((m_mode == Mode::PanoramaSphere || m_mode == Mode::Sphere) ? 2 : 1) * m_fbWidth;
  }

  Mode m_mode = Mode::Hemisphere;
};
