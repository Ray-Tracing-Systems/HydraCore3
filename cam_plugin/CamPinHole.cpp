#pragma once
#include "CamPinHole.h"
#include "../include/cglobals.h"

#include <cassert>
#include <cmath>
#include <cfloat>

using LiteMath::perspectiveMatrix;
using LiteMath::lookAt;
using LiteMath::inverse4x4;

CamPinHole::CamPinHole(){}
CamPinHole::~CamPinHole() {}

void CamPinHole::SetParameters(int a_width, int a_height, const CamParameters& a_params)
{
  m_width  = a_width;
  m_height = a_height;
      
  m_proj    = perspectiveMatrix(a_params.fov, a_params.aspect, a_params.nearPlane, a_params.farPlane);
  m_projInv = inverse4x4(m_proj);
}

void CamPinHole::MakeRaysBlock(float* out_rayPosAndNear4f, float* out_rayDirAndFar4f, size_t in_blockSize, int passId)
{
  assert(in_blockSize == m_width*m_height);
  #pragma omp parallel for default(shared)
  for(int pixelId = 0; pixelId < int(in_blockSize); ++pixelId)
    MakeEyeRay(pixelId, (float4*)out_rayPosAndNear4f, (float4*)out_rayDirAndFar4f);
}

void CamPinHole::AddSamplesContributionBlock(float* out_color4f, const float* colors4f, size_t in_blockSize, 
                                             uint32_t a_width, uint32_t a_height, int passId)
{
  #pragma omp parallel for default(shared)
  for(int pixelId = 0; pixelId < int(in_blockSize); ++pixelId)
    ContribSample(pixelId, (const float4*)colors4f, (float4*)out_color4f);
}


void CamPinHole::MakeEyeRay(int tid, float4* out_rayPosAndNear4f, float4* out_rayDirAndFar4f)
{
  kernel_MakeEyeRay(tid, out_rayPosAndNear4f, out_rayDirAndFar4f);
}

void CamPinHole::ContribSample(int tid, const float4* in_color, float4* out_color)
{
  kernel_ContribSample(tid, in_color, out_color);  
}

void CamPinHole::kernel_MakeEyeRay(int tid, float4* out_rayPosAndNear4f, float4* out_rayDirAndFar4f)
{
  const int x = tid % m_width;
  const int y = tid / m_height;

  float3 rayDir = EyeRayDirNormalized(float(x+0.5f)/float(m_width), float(y+0.5f)/float(m_height), m_projInv);
  float3 rayPos = float3(0,0,0);
  
  out_rayPosAndNear4f[tid] = to_float4(rayPos, 0.0f);
  out_rayDirAndFar4f [tid] = to_float4(rayDir, FLT_MAX);
}

void CamPinHole::kernel_ContribSample(int tid, const float4* in_color, float4* out_color)
{
  const int x = tid % m_width;
  const int y = tid / m_height;
  out_color[y*m_width+x] += in_color[tid];
}
