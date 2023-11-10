#pragma once
#include "CamPinHole.h"

using LiteMath::perspectiveMatrix;
using LiteMath::lookAt;
using LiteMath::inverse4x4;

CamPinHole::CamPinHole(){}
CamPinHole::~CamPinHole() {}

void CamPinHole::SetParameters(int a_width, int a_height, const CamParameters& a_params)
{
  m_width  = a_width;
  m_height = a_height;
      
  m_proj         = perspectiveMatrix(a_params.fov, a_params.aspect, a_params.nearPlane, a_params.farPlane);;
  m_worldView    = lookAt(float3(a_params.pos), float3(a_params.lookAt), float3(a_params.up));;
  m_projInv      = inverse4x4(m_proj);
  m_worldViewInv = inverse4x4(m_worldViewInv);
}

void CamPinHole::MakeRaysBlock(float* out_rayPosAndNear4f, float* out_rayDirAndFar4f, size_t in_blockSize, int passId)
{

}

void CamPinHole::AddSamplesContributionBlock(float* out_color4f, const float* colors4f, size_t in_blockSize, 
                                             uint32_t a_width, uint32_t a_height, int passId)
{
 
}
