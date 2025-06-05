#include "integrator_pt.h"
#include "fourier/fspec.h"

#include "include/cmaterial.h"
#include "include/cmat_gltf.h"
#include "include/cmat_conductor.h"
#include "utils.h"

#include <chrono>
#include <string>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;


void Integrator::PathTraceBlockF(uint tid, uint channels, float* fourierBuf, uint a_passNum)
{
  ConsoleProgressBar progress(tid);
  progress.Start();
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for (int i = 0; i < tid; ++i) {
    for (int j = 0; j < a_passNum; ++j) {
      PathTraceF(uint(i), fourierBuf);
    }
    progress.Update();
  }
  progress.Done();
  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}


