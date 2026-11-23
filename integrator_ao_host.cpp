#include "integrator_pt.h"
#include "mlt/rnd_qmc.h"
#include <chrono>
#include <stdexcept>

void Integrator::InitAmbientOcclusion(uint samples, float radius)
{
  if(samples == 0 || (samples & (samples - 1)) != 0 || radius <= 0.0f)
    throw std::runtime_error("RTAO requires a positive power-of-two spp and a positive -ao-radius");
  m_aoSamples = samples;
  m_aoRadius = radius;
  // The same four-dimensional Sobol bucketing as HeavyRT's InitQMCTable:
  // dimensions 0/1 select a pixel in a 128x128 tile, dimensions 2/3 select a ray.
  m_aoRandomsTile.resize(size_t(AO_TILE_SIZE) * AO_TILE_SIZE * samples);
  std::vector<uint> counters(AO_TILE_SIZE * AO_TILE_SIZE, 0);
  unsigned table[qmc::QRNG_DIMENSIONS][qmc::QRNG_RESOLUTION];
  qmc::init(table);
  for(size_t i = 0; i < m_aoRandomsTile.size(); ++i)
  {
    const float x = qmc::rndFloat(uint(i), 0, table[0]);
    const float y = qmc::rndFloat(uint(i), 1, table[0]);
    const float z = qmc::rndFloat(uint(i), 2, table[0]);
    const float w = qmc::rndFloat(uint(i), 3, table[0]);
    const uint ix = std::min(uint(x * float(AO_TILE_SIZE)), AO_TILE_SIZE - 1);
    const uint iy = std::min(uint(y * float(AO_TILE_SIZE)), AO_TILE_SIZE - 1);
    const uint pixel = iy * AO_TILE_SIZE + ix;
    if(counters[pixel] < samples)
      m_aoRandomsTile[size_t(pixel) * samples + counters[pixel]++] = float2(z, w);
  }
  if(std::any_of(counters.begin(), counters.end(), [samples](uint n) { return n != samples; }))
    throw std::runtime_error("Sobol AO tile did not contain the requested samples for every pixel");
}

void Integrator::AmbientOcclusionBlock(uint tid, uint channels, float* out_color, uint a_passNum)
{
  const auto start = std::chrono::high_resolution_clock::now();
  // HeavyRT's optional traversal metrics mutate shared counters and sets.
  #ifndef ENABLE_METRICS
  #pragma omp parallel for
  #endif
  for(int i = 0; i < int(tid); ++i)
    for(uint pass = 0; pass < a_passNum; ++pass)
      AmbientOcclusion(uint(i), channels, out_color);
  aoTime = float(std::chrono::duration_cast<std::chrono::microseconds>(
      std::chrono::high_resolution_clock::now() - start).count()) / 1000.0f;
}
