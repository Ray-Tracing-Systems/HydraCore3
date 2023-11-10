#pragma once
#include "CamPinHole.h"

CamPinHole::CamPinHole()
{

}

CamPinHole::~CamPinHole()
{

}

void CamPinHole::SetParameters(int a_width, int a_height, const float a_projInvMatrix[16], const CamParameters& a_params)
{
  
}

void CamPinHole::MakeRaysBlock(float* out_rayPosAndNear4f, float* out_rayDirAndFar4f, size_t in_blockSize, int passId)
{

}

void CamPinHole::AddSamplesContributionBlock(float* out_color4f, const float* colors4f, size_t in_blockSize, 
                                             uint32_t a_width, uint32_t a_height, int passId)
{
 
}
