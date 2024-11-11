#include "integrator_pt.h"
#include "include/crandom.h"

#include "include/cmaterial.h"
#include "include/cmat_gltf.h"
#include "include/cmat_conductor.h"
#include "include/cmat_glass.h"
#include "include/cmat_diffuse.h"
#include "include/cmat_plastic.h"
#include "include/cmat_dielectric.h"

#include <chrono>
#include <string>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

LightSample Integrator::LightSampleRev(int a_lightId, float3 rands, float3 illiminationPoint)
{
  const float2 rands2 = float2(rands.x, rands.y);

  return areaLightSampleRev  (m_lights.data() + a_lightId, rands2);
}

float Integrator::LightPdfSelectRev(int a_lightId) 
{ 
  return 1.0f/float(m_lights.size()); // uniform select
}

//static inline float DistanceSquared(float3 a, float3 b)
//{
//  const float3 diff = b - a;
//  return dot(diff, diff);
//}

float Integrator::LightEvalPDF(int a_lightId, float3 illuminationPoint, float3 ray_dir, const float3 lpos, const float3 lnorm, float a_envPdf)
{
  const uint gtype = m_lights[a_lightId].geomType;
  if(gtype == LIGHT_GEOM_ENV)
    return a_envPdf;

  const float hitDist   = length(illuminationPoint - lpos);
  const float cosValTmp = dot(ray_dir, -1.0f*lnorm);
  float cosVal = std::max(cosValTmp, 0.0f);
  
  return PdfAtoW(m_lights[a_lightId].pdfA, hitDist, cosVal);
}

float4 Integrator::LightIntensity(uint a_lightId, float4 a_wavelengths, float3 a_rayPos, float3 a_rayDir)
{
  float4 lightColor = m_lights[a_lightId].intensity;  
  
  lightColor *= m_lights[a_lightId].mult;

  return lightColor;
}

float4 Integrator::EnvironmentColor(float3 a_dir, float& outPdf)
{
  return m_envColor;
}

Integrator::Map2DPiecewiseSample Integrator::SampleMap2D(float3 rands, uint32_t a_tableOffset, int sizeX, int sizeY)
{
  const float fw = (float)sizeX;
  const float fh = (float)sizeY;
  const float fN = fw*fh;

  float pdf = 1.0f;
  int pixelOffset = SelectIndexPropToOpt(rands.z, m_pdfLightData.data() + a_tableOffset, sizeX*sizeY+1, &pdf);

  if (pixelOffset >= sizeX*sizeY)
    pixelOffset = sizeX*sizeY - 1;

  const int yPos = pixelOffset / sizeX;
  const int xPos = pixelOffset - yPos*sizeX;

  const float texX = (1.0f / fw)*(((float)(xPos) + 0.5f) + (rands.x*2.0f - 1.0f)*0.5f);
  const float texY = (1.0f / fh)*(((float)(yPos) + 0.5f) + (rands.y*2.0f - 1.0f)*0.5f);

  Map2DPiecewiseSample result;
  result.mapPdf   = pdf*fN; 
  result.texCoord = make_float2(texX, texY);
  return result;
}

void Integrator::GetExecutionTime(const char* a_funcName, float a_out[4])
{
  if(std::string(a_funcName) == "NaivePathTrace" || std::string(a_funcName) == "NaivePathTraceBlock")
    a_out[0] = naivePtTime;
  else if(std::string(a_funcName) == "PathTrace" || std::string(a_funcName) == "PathTraceBlock")
    a_out[0] = shadowPtTime;
  else if(std::string(a_funcName) == "RayTrace" || std::string(a_funcName) == "RayTraceBlock")
    a_out[0] = raytraceTime;
  else if(std::string(a_funcName) == "PathTraceFromInputRays" || std::string(a_funcName) == "PathTraceFromInputRaysBlock")
    a_out[0] = fromRaysPtTime;
}

void Integrator::ProgressBarStart()                  { _ProgressBarStart(); }
void Integrator::ProgressBarAccum(float a_progress)  { _ProgressBarAccum(a_progress); }
void Integrator::ProgressBarDone()                   { _ProgressBarDone(); }
