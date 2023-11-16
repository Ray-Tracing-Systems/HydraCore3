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
  m_randomGens.resize(a_maxThreads);
  #pragma omp parallel for default(shared)
  for(int i=0;i<a_maxThreads;i++)
    m_randomGens[i] = RandomGenInit(i);

  // init spectral curves
  m_cie_lambda = Get_CIE_lambda();
  m_cie_x      = Get_CIE_X();
  m_cie_y      = Get_CIE_Y();
  m_cie_z      = Get_CIE_Z();
}

void CamPinHole::MakeRaysBlock(float* out_rayPosAndNear4f, float* out_rayDirAndFar4f, AuxRayData* out_auxData, size_t in_blockSize, int passId)
{
  kernel1D_MakeEyeRay(int(in_blockSize), (float4*)out_rayPosAndNear4f, (float4*)out_rayDirAndFar4f, out_auxData);
}

void CamPinHole::AddSamplesContributionBlock(float* out_color4f, const float* colors4f, const AuxRayData* in_auxData, size_t in_blockSize, 
                                             uint32_t a_width, uint32_t a_height, int passId)
{
  kernel1D_ContribSample(int(in_blockSize), (const float4*)colors4f, in_auxData, (float4*)out_color4f); 
}

void CamPinHole::kernel1D_MakeEyeRay(int in_blockSize, float4* out_rayPosAndNear4f, float4* out_rayDirAndFar4f, AuxRayData* out_auxData)
{
  #pragma omp parallel for default(shared)
  for(int tid = 0; tid < int(in_blockSize); tid++)
  {
    const int x = tid % m_width;
    const int y = tid / m_height;
  
    float3 rayDir = EyeRayDirNormalized(float(x+0.5f)/float(m_width), float(y+0.5f)/float(m_height), m_projInv);
    float3 rayPos = float3(0,0,0);
  
    float4 wavelengths = float4(0,0,0,0);
    if(m_spectral_mode != 0)
    {
      auto genLocal = m_randomGens[tid];
      float u = rndFloat1_Pseudo(&genLocal);
      wavelengths = SampleWavelengths(u, LAMBDA_MIN, LAMBDA_MAX);
      m_randomGens[tid] = genLocal;
    }
    
    AuxRayData auxDat;
    auxDat.wavelengthsFixp[0] = uint16_t(65535.0f*((wavelengths.x - CAM_LAMBDA_MIN) / (CAM_LAMBDA_MAX - CAM_LAMBDA_MIN)));
    auxDat.wavelengthsFixp[1] = uint16_t(65535.0f*((wavelengths.y - CAM_LAMBDA_MIN) / (CAM_LAMBDA_MAX - CAM_LAMBDA_MIN)));
    auxDat.wavelengthsFixp[2] = uint16_t(65535.0f*((wavelengths.z - CAM_LAMBDA_MIN) / (CAM_LAMBDA_MAX - CAM_LAMBDA_MIN)));
    auxDat.wavelengthsFixp[3] = uint16_t(65535.0f*((wavelengths.w - CAM_LAMBDA_MIN) / (CAM_LAMBDA_MAX - CAM_LAMBDA_MIN)));
  
    out_rayPosAndNear4f[tid] = to_float4(rayPos, 0.0f);
    out_rayDirAndFar4f [tid] = to_float4(rayDir, FLT_MAX);
    out_auxData        [tid] = auxDat;
  }
}

void CamPinHole::kernel1D_ContribSample(int in_blockSize, const float4* in_color, const AuxRayData* in_auxData, float4* out_color)
{
  for(int tid = 0; tid < int(in_blockSize); tid++)
  {
    const int x = tid % m_width;
    const int y = tid / m_height;
    
    float4 color = in_color[tid];
    if(m_spectral_mode != 0) // TODO: spectral framebuffer
    {
      const float scale = (1.0f/65535.0f)*(CAM_LAMBDA_MAX - CAM_LAMBDA_MIN);
      AuxRayData auxDat  = in_auxData[tid];
      float4 wavelengths = float4(float(auxDat.wavelengthsFixp[0])*scale + CAM_LAMBDA_MIN,
                                  float(auxDat.wavelengthsFixp[1])*scale + CAM_LAMBDA_MIN,
                                  float(auxDat.wavelengthsFixp[2])*scale + CAM_LAMBDA_MIN,
                                  float(auxDat.wavelengthsFixp[3])*scale + CAM_LAMBDA_MIN);
      const float3 xyz = SpectrumToXYZ(color, wavelengths, CAM_LAMBDA_MIN, CAM_LAMBDA_MAX, m_cie_x.data(), m_cie_y.data(), m_cie_z.data());
      color = to_float4(XYZToRGB(xyz), 0.0f);
    }
    out_color[y*m_width+x] += color;
  }
}
