#pragma once

#include "integrator_pt.h"
#include "include/crandom.h"

#include "include/cmaterial.h"
#include "include/cmat_gltf.h"
#include "include/cmat_conductor.h"

#include <chrono>
#include <string>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

#include "utils.h"

void Integrator::PackXYBlock(uint tidX, uint tidY, uint a_passNum)
{
  #pragma omp parallel for default(shared)
  for(int y = 0; y < tidY; ++y)
    for(int x = 0; x < tidX; ++x)
      PackXY(uint(x), (uint)(y));
}

void Integrator::CastSingleRayBlock(uint tid, uint* out_color, uint a_passNum)
{ 
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for(int i = 0; i < tid; ++i)
    CastSingleRay(uint(i), out_color);
}

void Integrator::NaivePathTraceBlock(uint tid, float4* out_color, uint a_passNum)
{
  ConsoleProgressBar progress(tid);
  progress.Start();

  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for(int i = 0; i < tid; ++i) {
    for(int j = 0; j < a_passNum; ++j) {
      NaivePathTrace(uint(i), out_color);
    }
    progress.Update();
  }
  progress.Done();
  naivePtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}

void Integrator::PathTraceBlock(uint tid, float4* out_color, uint a_passNum)
{
  ConsoleProgressBar progress(tid);
  progress.Start();
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for (int i = 0; i < tid; ++i) {
    for (int j = 0; j < a_passNum; ++j) {
      PathTrace(uint(i), out_color);
    }
    progress.Update();
  }
  progress.Done();
  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}

void Integrator::RayTraceBlock(uint tid, float4* out_color, uint a_passNum)
{
  ConsoleProgressBar progress(tid);
  progress.Start();
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for(int i = 0; i < tid; ++i)
  {
    RayTrace(uint(i), out_color);
    progress.Update();
  }
  progress.Done();
  raytraceTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}

void Integrator::PathTraceFromInputRaysBlock(uint tid, const float4* in_rayPosAndNear, const float4* in_rayDirAndFar, const ushort4* in_wavesPacked, float4* out_color, uint a_passNum)
{
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for (int i = 0; i < tid; ++i) {
    for (int j = 0; j < a_passNum; ++j)
      PathTraceFromInputRays(uint(i), in_rayPosAndNear, in_rayDirAndFar, in_wavesPacked, out_color);
  }
  fromRaysPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}

