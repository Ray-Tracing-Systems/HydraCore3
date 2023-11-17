#pragma once
#include "CamPinHole.h"

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
  m_spectral_mode = a_params.spectralMode;
      
  m_proj    = perspectiveMatrix(a_params.fov, a_params.aspect, a_params.nearPlane, a_params.farPlane);
  m_projInv = inverse4x4(m_proj);
}

void CamPinHole::Init(int a_maxThreads)
{
  m_storedWaves.resize(a_maxThreads);
  m_randomGens.resize(a_maxThreads);
  #pragma omp parallel for default(shared)
  for(int i=0;i<a_maxThreads;i++)
    m_randomGens[i] = RandomGenInit(i + 12345*i);

  // init spectral curves
  m_cie_x      = Get_CIE_X();
  m_cie_y      = Get_CIE_Y();
  m_cie_z      = Get_CIE_Z();
}

void CamPinHole::MakeRaysBlock(RayPart1* out_rayPosAndNear4f, RayPart2* out_rayDirAndFar4f, size_t in_blockSize, int subPassId)
{
  kernel1D_MakeEyeRay(int(in_blockSize), out_rayPosAndNear4f, out_rayDirAndFar4f, subPassId);
}

void CamPinHole::AddSamplesContributionBlock(float* out_color4f, const float* colors4f, size_t in_blockSize, 
                                             uint32_t a_width, uint32_t a_height, int subPassId)
{
  kernel1D_ContribSample(int(in_blockSize), (const float4*)colors4f, (float4*)out_color4f, subPassId); 
}

void CamPinHole::kernel1D_MakeEyeRay(int in_blockSize, RayPart1* out_rayPosAndNear4f, RayPart2* out_rayDirAndFar4f, int subPassId)
{
  #pragma omp parallel for default(shared)
  for(int tid = 0; tid < int(in_blockSize); tid++)
  {
    const int x = (tid + subPassId*in_blockSize) % m_width;  // pitch-linear layout
    const int y = (tid + subPassId*in_blockSize) / m_height; // subPas is just a uniform slitting of image along the lines
    
    //if(x == 512 && y == 1023-100) // to debug target pixel
    //{
    //  int a = 2;
    //}

    float3 rayDir = EyeRayDirNormalized(float(x+0.5f)/float(m_width), float(y+0.5f)/float(m_height), m_projInv);
    float3 rayPos = float3(0,0,0);
  
    float4 wavelengths = float4(0,0,0,0);
    if(m_spectral_mode != 0)
    {
      auto genLocal     = m_randomGens[tid];
      float u           = rndFloat1_Pseudo(&genLocal);
      wavelengths       = SampleWavelengths(u, CAM_LAMBDA_MIN, CAM_LAMBDA_MAX);
      m_randomGens[tid] = genLocal;
    }

    RayPart1 p1;
    RayPart2 p2;

    const uint32_t wavesX = uint32_t(65535.0f*((wavelengths.x - CAM_LAMBDA_MIN) / (CAM_LAMBDA_MAX - CAM_LAMBDA_MIN)));
    const uint32_t wavesY = uint32_t(65535.0f*((wavelengths.y - CAM_LAMBDA_MIN) / (CAM_LAMBDA_MAX - CAM_LAMBDA_MIN)));
    const uint32_t wavesZ = uint32_t(65535.0f*((wavelengths.z - CAM_LAMBDA_MIN) / (CAM_LAMBDA_MAX - CAM_LAMBDA_MIN)));
    const uint32_t wavesW = uint32_t(65535.0f*((wavelengths.w - CAM_LAMBDA_MIN) / (CAM_LAMBDA_MAX - CAM_LAMBDA_MIN)));

    p1.origin[0] = rayPos[0];
    p1.origin[1] = rayPos[1];
    p1.origin[2] = rayPos[2];
    p1.pwaves01  = packXY1616(wavesX,wavesY); // 

    p2.direction[0] = rayDir[0];
    p2.direction[1] = rayDir[1];
    p2.direction[2] = rayDir[2];
    p2.pwaves23     = packXY1616(wavesZ,wavesW); // 
  
    out_rayPosAndNear4f[tid] = p1;
    out_rayDirAndFar4f [tid] = p2;
    m_storedWaves      [tid] = uint2(p1.pwaves01,  p2.pwaves23); // just remember waves in our buffer for camera
  }
}

void CamPinHole::kernel1D_ContribSample(int in_blockSize, const float4* in_color, float4* out_color, int subPassId)
{
  for(int tid = 0; tid < int(in_blockSize); tid++)
  {
    const int x = (tid + subPassId*in_blockSize) % m_width;  // pitch-linear layout
    const int y = (tid + subPassId*in_blockSize) / m_height; // subPas is just a uniform slitting of image along the lines

    float4 color = in_color[tid];
    if(m_spectral_mode != 0) // TODO: spectral framebuffer
    {
      const float scale = (1.0f/65535.0f)*(CAM_LAMBDA_MAX - CAM_LAMBDA_MIN);
      const uint2 wavesPacked = m_storedWaves[tid];
      const uint2 wavesXY = unpackXY1616(wavesPacked.x);
      const uint2 wavesZW = unpackXY1616(wavesPacked.y);
      float4 wavelengths = float4(float(wavesXY[0])*scale + CAM_LAMBDA_MIN,
                                  float(wavesXY[1])*scale + CAM_LAMBDA_MIN,
                                  float(wavesZW[0])*scale + CAM_LAMBDA_MIN,
                                  float(wavesZW[1])*scale + CAM_LAMBDA_MIN);
                                  
      const float3 xyz = SpectrumToXYZ(color, wavelengths, CAM_LAMBDA_MIN, CAM_LAMBDA_MAX, m_cie_x.data(), m_cie_y.data(), m_cie_z.data());
      color = to_float4(XYZToRGB(xyz), 1.0f);
    }

    out_color[y*m_width+x] += color;
  }
}
