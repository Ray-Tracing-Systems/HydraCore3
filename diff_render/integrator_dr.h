#pragma once

#include "include/cglobals.h" // We assume that all code that should pe passed to kernels will be just included both for CPU and OpenCL
#include "include/crandom.h"
#include "include/clight.h"
#include "include/cmaterial.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <cfloat>

#include "CrossRT.h" // special include for ray tracing
#include "Image2d.h" // special include for textures

#include "spectrum.h"
#include "integrator_pt.h"

using LiteImage::ICombinedImageSampler;

class IntegratorDR : public Integrator
{
public:

  IntegratorDR(int a_maxThreads = 1, int a_spectral_mode = 0, std::vector<uint32_t> a_features = {}) : Integrator(a_maxThreads, a_spectral_mode, a_features){}
  virtual ~IntegratorDR() {}

  void PathTraceDR(uint tid, uint channels, float* out_color, uint a_passNum,
                   const float* a_data, float* a_dataGrad);

protected:
  
};

