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

  Mode m_mode = Mode::Hemisphere;
};
