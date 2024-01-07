#include "integrator_dr.h"
#include "utils.h"

#include "include/cmaterial.h"
#include "include/cmat_gltf.h"
#include "include/cmat_conductor.h"
#include "include/cmat_glass.h"
#include "include/cmat_diffuse.h"
#include "include/cmat_plastic.h"

#include <chrono>
#include <string>

#include <omp.h>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

void IntegratorDR::RecordPixelRndIfNeeded(float2 offsets, float u)
{
  auto cpuThreadId = omp_get_thread_num();
  m_recorded[cpuThreadId].pixelOffsets = offsets;
  m_recorded[cpuThreadId].waveSelector = u;
}

void IntegratorDR::RecordRayHitIfNeeded(uint32_t bounceId, CRT_Hit hit)
{
  auto cpuThreadId = omp_get_thread_num();
  m_recorded[cpuThreadId].perBounce[bounceId].hit = hit;
}

void IntegratorDR::RecordShadowHitIfNeeded(uint32_t bounceId, bool inShadow)
{
  auto cpuThreadId = omp_get_thread_num();
  m_recorded[cpuThreadId].perBounce[bounceId].inShadow = inShadow ? 1 : 0;
}

void IntegratorDR::RecordLightRndIfNeeded(uint32_t bounceId, int lightId, float2 rands)
{
  auto cpuThreadId = omp_get_thread_num();
  m_recorded[cpuThreadId].perBounce[bounceId].lightId  = lightId;
  m_recorded[cpuThreadId].perBounce[bounceId].lgtRands = rands;
}

void IntegratorDR::RecordMatRndNeeded(uint32_t bounceId, float4 rands)
{
  auto cpuThreadId = omp_get_thread_num();
  m_recorded[cpuThreadId].perBounce[bounceId].matRands = rands;
}

void IntegratorDR::RecordBlendRndNeeded(uint32_t bounceId, uint layer, float rand)
{
  auto cpuThreadId = omp_get_thread_num();
  m_recorded[cpuThreadId].perBounce[bounceId].blendRnd[layer] = rand;
}
