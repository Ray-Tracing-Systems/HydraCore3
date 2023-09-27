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
  for(uint y=0;y<tidY;y++)
    for(uint x=0;x<tidX;x++)
      PackXY(x, y);
}

void Integrator::CastSingleRayBlock(uint tid, uint* out_color, uint a_passNum)
{ 
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for(uint i=0;i<tid;i++)
    CastSingleRay(i, out_color);
}

void Integrator::NaivePathTraceBlock(uint tid, float4* out_color, uint a_passNum)
{
  ConsoleProgressBar progress(tid);
  progress.Start();

  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for(uint i=0;i<tid;i++)
  {
    for(uint j=0;j<a_passNum;j++)
    {
      NaivePathTrace(i, out_color);
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
  for(uint i=0;i<tid;i++)
  {
    for(uint j=0;j<a_passNum;j++)
    {
      PathTrace(i, out_color);
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
  for(uint i=0;i<tid;i++)
  {
    RayTrace(i, out_color);
    progress.Update();
  }
  progress.Done();
  raytraceTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}

