#include "integrator_pt.h"
#include "include/crandom.h"

#include <chrono>
#include <string>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

void Integrator::InitRandomGens(int a_maxThreads)
{
  m_randomGens.resize(a_maxThreads);
  #pragma omp parallel for default(shared)
  for(int i=0;i<a_maxThreads;i++)
    m_randomGens[i] = RandomGenInit(i);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Integrator::kernel_InitEyeRay2(uint tid, const uint* packedXY, 
                                   float4* rayPosAndNear, float4* rayDirAndFar, float4* wavelengths, 
                                   float4* accumColor,    float4* accumuThoroughput,
                                   RandomGen* gen, uint* rayFlags, MisData* misData) // 
{
  if(tid >= m_maxThreadId)
    return;

  *accumColor        = make_float4(0,0,0,0);
  *accumuThoroughput = make_float4(1,1,1,1);
  RandomGen genLocal = m_randomGens[tid];
  *rayFlags          = 0;
  *misData           = makeInitialMisData();

  const uint XY = packedXY[tid];

  const uint x = (XY & 0x0000FFFF);
  const uint y = (XY & 0xFFFF0000) >> 16;
  const float2 pixelOffsets = rndFloat2_Pseudo(&genLocal);

  if(x == 70 && y == 512 - 250 - 1)
  {
    int a = 2;
  }

  float3 rayDir = EyeRayDirNormalized((float(x) + pixelOffsets.x)/float(m_winWidth), 
                                      (float(y) + pixelOffsets.y)/float(m_winHeight), m_projInv);
  float3 rayPos = float3(0,0,0);

  transform_ray3f(m_worldViewInv, &rayPos, &rayDir);
  
  float tmp = 0.0f;
  if(KSPEC_SPECTRAL_RENDERING !=0 && m_spectral_mode != 0)
  {
    float u = rndFloat1_Pseudo(&genLocal);
    *wavelengths = SampleWavelengths(u, LAMBDA_MIN, LAMBDA_MAX);
    tmp = u;
  }
  else
  {
    const uint32_t sample_sz = sizeof((*wavelengths).M) / sizeof((*wavelengths).M[0]);
    for (uint32_t i = 0; i < sample_sz; ++i) 
      (*wavelengths)[i] = 0.0f;
  }

  RecordPixelRndIfNeeded(pixelOffsets, tmp);
 
  *rayPosAndNear = to_float4(rayPos, 0.0f);
  *rayDirAndFar  = to_float4(rayDir, FLT_MAX);
  *gen           = genLocal;
}

void Integrator::kernel_InitEyeRayFromInput(uint tid, const RayPosAndW* in_rayPosAndNear, const RayDirAndT* in_rayDirAndFar,
                                            float4* rayPosAndNear, float4* rayDirAndFar, float4* accumColor, float4* accumuThoroughput, 
                                            RandomGen* gen, uint* rayFlags, MisData* misData, float4* wavelengths)
{
  if(tid >= m_maxThreadId)
    return;

  *accumColor        = make_float4(0,0,0,0);
  *accumuThoroughput = make_float4(1,1,1,1);
  *rayFlags          = 0;
  *misData           = makeInitialMisData();  

  //const int x = int(tid) % m_winWidth;
  //const int y = int(tid) / m_winHeight;

  const RayPosAndW rayPosData = in_rayPosAndNear[tid];
  const RayDirAndT rayDirData = in_rayDirAndFar[tid];

  float3 rayPos = float3(rayPosData.origin[0], rayPosData.origin[1], rayPosData.origin[2]);
  float3 rayDir = float3(rayDirData.direction[0], rayDirData.direction[1], rayDirData.direction[2]);
  transform_ray3f(m_worldViewInv, &rayPos, &rayDir);

  if(KSPEC_SPECTRAL_RENDERING !=0 && m_spectral_mode != 0)
  {
    *wavelengths = float4(rayPosData.wave);
    //const uint2 wavesXY = unpackXY1616(rayPosData.waves01);
    //const uint2 wavesZW = unpackXY1616(rayDirData.waves23);
    //const float scale = (1.0f/65535.0f)*(LAMBDA_MAX - LAMBDA_MIN);
    //*wavelengths = float4(float(wavesXY[0])*scale + LAMBDA_MIN,
    //                      float(wavesXY[1])*scale + LAMBDA_MIN,
    //                      float(wavesZW[0])*scale + LAMBDA_MIN,
    //                      float(wavesZW[1])*scale + LAMBDA_MIN);
  }
  else
    *wavelengths = float4(0,0,0,0);

  *rayPosAndNear = to_float4(rayPos, 0.0f);
  *rayDirAndFar  = to_float4(rayDir, FLT_MAX);
  *gen           = m_randomGens[tid];
}


void Integrator::kernel_RayTrace2(uint tid, uint bounce, const float4* rayPosAndNear, const float4* rayDirAndFar,
                                  float4* out_hit1, float4* out_hit2, float4* out_hit3, uint* out_instId, uint* rayFlags)
{
  if(tid >= m_maxThreadId)
    return;
  uint currRayFlags = *rayFlags;
  if(isDeadRay(currRayFlags))
    return;

  const float4 rayPos = *rayPosAndNear;
  const float4 rayDir = *rayDirAndFar ;

  const CRT_Hit hit   = m_pAccelStruct->RayQuery_NearestHit(rayPos, rayDir);
  RecordRayHitIfNeeded(bounce, hit);

  if(hit.geomId != uint32_t(-1))
  {
    const float2 uv     = float2(hit.coords[0], hit.coords[1]);
    
    // slightly undershoot the intersection to prevent self-intersection and other bugs
    const float3 hitPos = to_float3(rayPos) + hit.t * (1.f - 1e-6f) * to_float3(rayDir);
    // const float3 overHit  = to_float3(rayPos) + hit.t * (1.f + 1e-6f) * to_float3(rayDir);

    // alternative, you may consider Johannes Hanika solution from  Ray Tracing Gems2  
    /////////////////////////////////////////////////////////////////////////////////
    // // get distance vectors from triangle vertices
    // vec3 tmpu = P - A, tmpv = P - B, tmpw = P - C
    // // project these onto the tangent planes
    // // defined by the shading normals
    // float dotu = min (0.0, dot(tmpu , nA))
    // float dotv = min (0.0, dot(tmpv , nB))
    // float dotw = min (0.0, dot(tmpw , nC))
    // tmpu -= dotu*nA
    // tmpv -= dotv*nB
    // tmpw -= dotw*nC
    // // finally P' is the barycentric mean of these three
    // vec3 Pp = P + u*tmpu + v*tmpv + w*tmpw
    /////////////////////////////////////////////////////////////////////////////////

    const uint triOffset  = m_matIdOffsets[hit.geomId];
    const uint vertOffset = m_vertOffset  [hit.geomId];
  
    const uint A = m_triIndices[(triOffset + hit.primId)*3 + 0];
    const uint B = m_triIndices[(triOffset + hit.primId)*3 + 1];
    const uint C = m_triIndices[(triOffset + hit.primId)*3 + 2];

    const float4 data1 = (1.0f - uv.x - uv.y)*m_vNorm4f[A + vertOffset] + uv.y*m_vNorm4f[B + vertOffset] + uv.x*m_vNorm4f[C + vertOffset];
    const float4 data2 = (1.0f - uv.x - uv.y)*m_vTang4f[A + vertOffset] + uv.y*m_vTang4f[B + vertOffset] + uv.x*m_vTang4f[C + vertOffset];

    float3 hitNorm     = to_float3(data1);
    float3 hitTang     = to_float3(data2);
    float2 hitTexCoord = float2(data1.w, data2.w);

    // transform surface point with matrix and flip normal if needed
    //
    hitNorm                = normalize(mul3x3(m_normMatrices[hit.instId], hitNorm));
    hitTang                = normalize(mul3x3(m_normMatrices[hit.instId], hitTang));
    const float flipNorm   = dot(to_float3(rayDir), hitNorm) > 0.001f ? -1.0f : 1.0f; // beware of transparent materials which use normal sign to identity "inside/outside" glass for example
    hitNorm                = flipNorm * hitNorm;
    hitTang                = flipNorm * hitTang; // do we need this ??
    
    if (flipNorm < 0.0f) currRayFlags |=  RAY_FLAG_HAS_INV_NORMAL;
    else                 currRayFlags &= ~RAY_FLAG_HAS_INV_NORMAL;

    const uint midOriginal = m_matIdByPrimId[m_matIdOffsets[hit.geomId] + hit.primId];
    const uint midRemaped  = RemapMaterialId(midOriginal, hit.instId);

    *rayFlags              = packMatId(currRayFlags, midRemaped);
    *out_hit1              = to_float4(hitPos,  hitTexCoord.x); 
    *out_hit2              = to_float4(hitNorm, hitTexCoord.y);
    *out_hit3              = to_float4(hitTang, hit.t);
    *out_instId            = hit.instId;
  }
  else
    *rayFlags              = currRayFlags | (RAY_FLAG_IS_DEAD | RAY_FLAG_OUT_OF_SCENE);
}

float4 Integrator::GetLightSourceIntensity(uint a_lightId, const float4* a_wavelengths, float3 a_rayPos, float3 a_rayDir)
{
  float4 lightColor = m_lights[a_lightId].intensity;  
  const uint specId = m_lights[a_lightId].specId;
  if(KSPEC_SPECTRAL_RENDERING !=0 && m_spectral_mode != 0 && specId < 0xFFFFFFFF)
  {
    const uint2 data  = m_spec_offset_sz[specId];
    const uint offset = data.x;
    const uint size   = data.y;
    lightColor = SampleUniformSpectrum(m_spec_values.data() + offset, *a_wavelengths, size);
  }
  lightColor *= m_lights[a_lightId].mult;
  
  uint iesId = m_lights[a_lightId].iesId;
  if(iesId != uint(-1))
  {
    if((m_lights[a_lightId].flags & LIGHT_FLAG_POINT_AREA) != 0)
      a_rayDir = normalize(to_float3(m_lights[a_lightId].pos) - a_rayPos);
    const float3 dirTrans = to_float3(m_lights[a_lightId].iesMatrix*to_float4(a_rayDir, 0.0f));
    float sintheta        = 0.0f;
    const float2 texCoord = sphereMapTo2DTexCoord((-1.0f)*dirTrans, &sintheta);
    const float4 texColor = m_textures[iesId]->sample(texCoord);
    lightColor *= texColor;
  }

  return lightColor;
}


void Integrator::kernel_SampleLightSource(uint tid, const float4* rayPosAndNear, const float4* rayDirAndFar, const float4* wavelengths,
                                          const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_hitPart3,
                                          const uint* rayFlags, uint bounce, RandomGen* a_gen, float4* out_shadeColor)
{
  if(tid >= m_maxThreadId)
    return;
  const uint currRayFlags = *rayFlags;
  if(isDeadRay(currRayFlags))
    return;
    
  const uint32_t matId = extractMatId(currRayFlags);
  const float3 ray_dir = to_float3(*rayDirAndFar);
  
  const float4 data1  = *in_hitPart1;
  const float4 data2  = *in_hitPart2;
  const float4 lambda = *wavelengths;

  SurfaceHit hit;
  hit.pos  = to_float3(data1);
  hit.norm = to_float3(data2);
  hit.tang = to_float3(*in_hitPart3);
  hit.uv   = float2(data1.w, data2.w);

  const float2 rands = rndFloat2_Pseudo(a_gen); // don't use single rndFloat4 (!!!)
  const float rndId  = rndFloat1_Pseudo(a_gen); // don't use single rndFloat4 (!!!)
  const int lightId  = std::min(int(std::floor(rndId * float(m_lights.size()))), int(m_lights.size() - 1u));
  RecordLightRndIfNeeded(bounce, lightId, rands);

  if(lightId < 0) // no lights or invalid light id
  {
    *out_shadeColor = float4(0.0f, 0.0f, 0.0f, 0.0f);
    return;
  }
  
  const LightSample lSam = LightSampleRev(lightId, rands, hit.pos);
  const float  hitDist   = std::sqrt(dot(hit.pos - lSam.pos, hit.pos - lSam.pos));

  const float3 shadowRayDir = normalize(lSam.pos - hit.pos); // explicitSam.direction;
  const float3 shadowRayPos = hit.pos + hit.norm * std::max(maxcomp(hit.pos), 1.0f)*5e-6f; // TODO: see Ray Tracing Gems, also use flatNormal for offset
  const bool   inShadow     = m_pAccelStruct->RayQuery_AnyHit(to_float4(shadowRayPos, 0.0f), to_float4(shadowRayDir, hitDist*0.9995f));
  const bool   inIllumArea  = (dot(shadowRayDir, lSam.norm) < 0.0f) || lSam.isOmni || lSam.hasIES;

  RecordShadowHitIfNeeded(bounce, inShadow);

  if(!inShadow && inIllumArea) 
  {
    const BsdfEval bsdfV    = MaterialEval(matId, lambda, shadowRayDir, (-1.0f)*ray_dir, hit.norm, hit.tang, hit.uv);
    float cosThetaOut       = std::max(dot(shadowRayDir, hit.norm), 0.0f);
    
    float      lgtPdfW      = LightPdfSelectRev(lightId) * LightEvalPDF(lightId, shadowRayPos, shadowRayDir, lSam.pos, lSam.norm);
    float      misWeight    = (m_intergatorType == INTEGRATOR_MIS_PT) ? misWeightHeuristic(lgtPdfW, bsdfV.pdf) : 1.0f;
    const bool isDirect     = (m_lights[lightId].geomType == LIGHT_GEOM_DIRECT); 
    const bool isPoint      = (m_lights[lightId].geomType == LIGHT_GEOM_POINT); 
    
    if(isDirect)
    {
      misWeight = 1.0f;
      lgtPdfW   = 1.0f;
    }
    else if(isPoint)
      misWeight = 1.0f;

    if(m_skipBounce >= 1 && int(bounce) < int(m_skipBounce)-1) // skip some number of bounces if this is set
      misWeight = 0.0f;
    
    const float4 lightColor = GetLightSourceIntensity(lightId, wavelengths, shadowRayPos, shadowRayDir);
    *out_shadeColor = (lightColor * bsdfV.val / lgtPdfW) * cosThetaOut * misWeight;
  }
  else
    *out_shadeColor = float4(0.0f, 0.0f, 0.0f, 0.0f);
}

void Integrator::kernel_NextBounce(uint tid, uint bounce, const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_hitPart3,
                                   const uint* in_instId, const float4* in_shadeColor, float4* rayPosAndNear, float4* rayDirAndFar,
                                   const float4* wavelengths, float4* accumColor, float4* accumThoroughput,
                                   RandomGen* a_gen, MisData* misPrev, uint* rayFlags)
{
  if(tid >= m_maxThreadId)
    return;
  const uint currRayFlags = *rayFlags;
  if(isDeadRay(currRayFlags))
    return;
    
  const uint32_t matId = extractMatId(currRayFlags);

  // process surface hit case
  //
  const float3 ray_dir = to_float3(*rayDirAndFar);
  const float3 ray_pos = to_float3(*rayPosAndNear);
  const float4 lambda  = *wavelengths;
  
  const float4 data1 = *in_hitPart1;
  const float4 data2 = *in_hitPart2;
  
  SurfaceHit hit;
  hit.pos  = to_float3(data1);
  hit.norm = to_float3(data2);
  hit.tang = to_float3(*in_hitPart3);
  hit.uv   = float2(data1.w, data2.w);

  const float hitDist = in_hitPart3->w;
  
  const MisData prevBounce = *misPrev;
  const float   prevPdfW   = prevBounce.matSamplePdf;

  // process light hit case
  //
  if(m_materials[matId].mtype == MAT_TYPE_LIGHT_SOURCE)
  {
    const uint   texId     = m_materials[matId].texid[0];
    const float2 texCoordT = mulRows2x4(m_materials[matId].row0[0], m_materials[matId].row1[0], hit.uv);
    const float4 texColor  = m_textures[texId]->sample(texCoordT);
    const uint   lightId   = m_instIdToLightInstId[*in_instId]; 
    
    const float4 emissColor = m_materials[matId].colors[EMISSION_COLOR];
    float4 lightIntensity   = emissColor * texColor;

    if(lightId != 0xFFFFFFFF)
    {
      const float lightCos = dot(to_float3(*rayDirAndFar), to_float3(m_lights[lightId].norm));
      const float lightDirectionAtten = (lightCos < 0.0f || m_lights[lightId].geomType == LIGHT_GEOM_SPHERE) ? 1.0f : 0.0f;
      lightIntensity = GetLightSourceIntensity(lightId, wavelengths, ray_pos, to_float3(*rayDirAndFar))*lightDirectionAtten;
    }

    float misWeight = 1.0f;
    if(m_intergatorType == INTEGRATOR_MIS_PT) 
    {
      if(bounce > 0)
      {
        if(lightId != 0xFFFFFFFF)
        {
          const float lgtPdf  = LightPdfSelectRev(lightId) * LightEvalPDF(lightId, ray_pos, ray_dir, hit.pos, hit.norm);
          misWeight           = misWeightHeuristic(prevPdfW, lgtPdf);
          if (prevPdfW <= 0.0f) // specular bounce
            misWeight = 1.0f;
        }
      }
    }
    else if(m_intergatorType == INTEGRATOR_SHADOW_PT && hasNonSpecular(currRayFlags))
      misWeight = 0.0f;
    
    if(m_skipBounce >= 1 && bounce < m_skipBounce) // skip some number of bounces if this is set
      misWeight = 0.0f;

    float4 currAccumColor      = *accumColor;
    float4 currAccumThroughput = *accumThoroughput;
    
    currAccumColor += currAccumThroughput * lightIntensity * misWeight;
   
    *accumColor = currAccumColor;
    *rayFlags   = currRayFlags | (RAY_FLAG_IS_DEAD | RAY_FLAG_HIT_LIGHT);
    return;
  }
  
  const uint bounceTmp    = bounce;
  const BsdfSample matSam = MaterialSampleAndEval(matId, bounceTmp, lambda, a_gen, (-1.0f)*ray_dir, hit.norm, hit.tang, hit.uv, misPrev, currRayFlags);
  const float4 bxdfVal    = matSam.val * (1.0f / std::max(matSam.pdf, 1e-20f));
  const float  cosTheta   = std::abs(dot(matSam.dir, hit.norm)); 

  MisData nextBounceData      = *misPrev;        // remember current pdfW for next bounce
  nextBounceData.matSamplePdf = (matSam.flags & RAY_EVENT_S) != 0 ? -1.0f : matSam.pdf; 
  nextBounceData.cosTheta     = cosTheta;   
  *misPrev                    = nextBounceData;

  if(m_intergatorType == INTEGRATOR_STUPID_PT)
  {
    *accumThoroughput *= cosTheta * bxdfVal; 
  }
  else if(m_intergatorType == INTEGRATOR_SHADOW_PT || m_intergatorType == INTEGRATOR_MIS_PT)
  {
    const float4 currThoroughput = *accumThoroughput;
    const float4 shadeColor      = *in_shadeColor;
    float4 currAccumColor        = *accumColor;

    currAccumColor += currThoroughput * shadeColor;
    *accumColor       = currAccumColor;
    *accumThoroughput = currThoroughput*cosTheta*bxdfVal; 
  }

  // compute point on the other side of the surface in case of transmission
  if((matSam.flags & RAY_EVENT_T) != 0)
  {
    hit.pos = hit.pos + hitDist * ray_dir * 2 * 1e-6f;
  }

  *rayPosAndNear = to_float4(OffsRayPos(hit.pos, hit.norm, matSam.dir), 0.0f); // todo: use flatNormal for offset
  *rayDirAndFar  = to_float4(matSam.dir, FLT_MAX);
  *rayFlags      = currRayFlags | matSam.flags;
}

void Integrator::kernel_HitEnvironment(uint tid, const uint* rayFlags, const float4* rayDirAndFar, const MisData* a_prevMisData, const float4* accumThoroughput,
                                       float4* accumColor)
{
  if(tid >= m_maxThreadId)
    return;
  const uint currRayFlags = *rayFlags;
  if(!isOutOfScene(currRayFlags))
    return;
  
  // TODO: HDRI maps
  const float4 envData  = GetEnvironmentColorAndPdf(to_float3(*rayDirAndFar));
  // const float3 envColor = to_float3(envData)/envData.w;    // explicitly account for pdf; when MIS will be enabled, need to deal with MIS weight also!

  const float4 envColor = envData;
  if(m_intergatorType == INTEGRATOR_STUPID_PT)     // todo: when explicit sampling will be added, disable contribution here for 'INTEGRATOR_SHADOW_PT'
    *accumColor = (*accumThoroughput) * envColor;
  else
    *accumColor += (*accumThoroughput) * envColor;
}


void Integrator::kernel_ContributeToImage(uint tid, const uint* rayFlags, uint channels, const float4* a_accumColor, const RandomGen* gen,
                                          const uint* in_pakedXY, const float4* wavelengths, float* out_color)
{
  
  if(tid >= m_maxThreadId) // don't contrubute to image in any "record" mode
    return;
  
  m_randomGens[tid] = *gen;
  if(m_disableImageContrib !=0)
    return;

  const uint XY = in_pakedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;
  
  float4 specSamples = *a_accumColor; 
  
  //////////////////////////////////////////// Max Trofimov RGB camera responce model
  float3 tmpColor = to_float3(specSamples);
  float3 rgbR     = to_float3(m_camRespoceRGB[0]);
  float3 rgbG     = to_float3(m_camRespoceRGB[1]);
  float3 rgbB     = to_float3(m_camRespoceRGB[2]);
  float3 rgb      = float3(dot(tmpColor, rgbR), dot(tmpColor, rgbG), dot(tmpColor, rgbB));
  //////////////////////////////////////////// Max Trofimov RGB camera responce model
  
  if(KSPEC_SPECTRAL_RENDERING!=0 && m_spectral_mode != 0) 
  {
    float4 waves = *wavelengths;
    
    if(m_camResponseSpectrumId[0] < 0)
    {
      const float3 xyz = SpectrumToXYZ(specSamples, waves, LAMBDA_MIN, LAMBDA_MAX, m_cie_x.data(), m_cie_y.data(), m_cie_z.data(),
                                       terminateWavelngths(*rayFlags));
      rgb = XYZToRGB(xyz);
    }
    else
    {
      float4 responceX, responceY, responceZ;
      {
        int specId = m_camResponseSpectrumId[0];
        if(specId >= 0)
        {
          const uint2 data  = m_spec_offset_sz[specId];
          const uint offset = data.x;
          const uint size   = data.y;
          responceX = SampleUniformSpectrum(m_spec_values.data() + offset, waves, size);
        }
        else
          responceX = float4(1,1,1,1);

        specId = m_camResponseSpectrumId[1];
        if(specId >= 0)
        {
          const uint2 data  = m_spec_offset_sz[specId];
          const uint offset = data.x;
          const uint size   = data.y;
          responceY = SampleUniformSpectrum(m_spec_values.data() + offset, waves, size);
        }
        else
          responceY = responceX;

        specId = m_camResponseSpectrumId[2];
        if(specId >= 0)
        {
          const uint2 data  = m_spec_offset_sz[specId];
          const uint offset = data.x;
          const uint size   = data.y;
          responceZ = SampleUniformSpectrum(m_spec_values.data() + offset, waves, size);
        }
        else
          responceZ = responceY;
      }

      float3 xyz = float3(0,0,0);
      for (uint32_t i = 0; i < SPECTRUM_SAMPLE_SZ; ++i) {
        xyz.x += specSamples[i]*responceX[i];
        xyz.y += specSamples[i]*responceY[i];
        xyz.z += specSamples[i]*responceZ[i]; 
      } 

      if(m_camResponseType == CAM_RESPONCE_XYZ)
        rgb = XYZToRGB(xyz);
      else
        rgb = xyz;
    }
  }

  float4 colorRes = m_exposureMult * to_float4(rgb, 1.0f);
  //if(x == 415 && (y == 256-130-1))
  //{
  //  int a = 2;
  //  //colorRes = float4(1,0,0,0);
  //}
  
  if(channels == 1)
  {
    const float mono = 0.2126f*colorRes.x + 0.7152f*colorRes.y + 0.0722f*colorRes.z;
    out_color[y*m_winWidth+x] += mono;
  }
  else if(channels <= 4)
  {
    out_color[(y*m_winWidth+x)*channels + 0] += colorRes.x;
    out_color[(y*m_winWidth+x)*channels + 1] += colorRes.y;
    out_color[(y*m_winWidth+x)*channels + 2] += colorRes.z;
  }
  else
  {
    auto waves = (*wavelengths);
    auto color = (*a_accumColor)*m_exposureMult;
    for(int i=0;i<4;i++) {
      const float t         = (waves[i] - LAMBDA_MIN)/(LAMBDA_MAX-LAMBDA_MIN);
      const int channelId   = std::min(int(float(channels)*t), int(channels)-1);
      const int offsetPixel = int(y)*m_winWidth + int(x);
      const int offsetLayer = channelId*m_winWidth*m_winHeight;
      out_color[offsetLayer + offsetPixel] += color[i];
    }
  }

}

void Integrator::kernel_CopyColorToOutput(uint tid, uint channels, const float4* a_accumColor, const RandomGen* gen, float* out_color)
{
  if(tid >= m_maxThreadId)
    return;
  
  const float4 color = *a_accumColor;

  if(channels == 4)
  {
    out_color[tid*4+0] += color[0];
    out_color[tid*4+1] += color[1];
    out_color[tid*4+2] += color[2];
    out_color[tid*4+3] += color[3];
  }
  else if(channels == 1)
    out_color[tid] += color[0];

  m_randomGens[tid] = *gen;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Integrator::NaivePathTrace(uint tid, uint channels, float* out_color)
{
  float4 accumColor, accumThroughput;
  float4 rayPosAndNear, rayDirAndFar;
  float4 wavelengths;
  RandomGen gen; 
  MisData   mis;
  uint      rayFlags;
  kernel_InitEyeRay2(tid, m_packedXY.data(), &rayPosAndNear, &rayDirAndFar, &wavelengths, &accumColor, &accumThroughput, &gen, &rayFlags, &mis);

  for(uint depth = 0; depth < m_traceDepth + 1; ++depth) // + 1 due to NaivePT uses additional bounce to hit light 
  {
    float4 shadeColor, hitPart1, hitPart2, hitPart3;
    uint instId = 0;
    kernel_RayTrace2(tid, depth, &rayPosAndNear, &rayDirAndFar, &hitPart1, &hitPart2, &hitPart3, &instId, &rayFlags);
    if(isDeadRay(rayFlags))
      break;
    
    kernel_NextBounce(tid, depth, &hitPart1, &hitPart2, &hitPart3, &instId, &shadeColor,
                      &rayPosAndNear, &rayDirAndFar, &wavelengths, &accumColor, &accumThroughput, &gen, &mis, &rayFlags);
    if(isDeadRay(rayFlags))
      break;
  }

  kernel_HitEnvironment(tid, &rayFlags, &rayDirAndFar, &mis, &accumThroughput,
                       &accumColor);

  kernel_ContributeToImage(tid, &rayFlags, channels, &accumColor, &gen, m_packedXY.data(), &wavelengths, 
                           out_color);
}

void Integrator::PathTrace(uint tid, uint channels, float* out_color)
{
  float4 accumColor, accumThroughput;
  float4 rayPosAndNear, rayDirAndFar;
  float4 wavelengths;
  RandomGen gen; 
  MisData   mis;
  uint      rayFlags;
  kernel_InitEyeRay2(tid, m_packedXY.data(), &rayPosAndNear, &rayDirAndFar, &wavelengths, &accumColor, &accumThroughput, &gen, &rayFlags, &mis);

  for(uint depth = 0; depth < m_traceDepth; depth++) 
  {
    float4   shadeColor, hitPart1, hitPart2, hitPart3;
    uint instId;
    kernel_RayTrace2(tid, depth, &rayPosAndNear, &rayDirAndFar, &hitPart1, &hitPart2, &hitPart3, &instId, &rayFlags);
    if(isDeadRay(rayFlags))
      break;
    
    kernel_SampleLightSource(tid, &rayPosAndNear, &rayDirAndFar, &wavelengths, &hitPart1, &hitPart2, &hitPart3, &rayFlags,
                             depth, &gen, &shadeColor);

    kernel_NextBounce(tid, depth, &hitPart1, &hitPart2, &hitPart3, &instId, &shadeColor,
                      &rayPosAndNear, &rayDirAndFar, &wavelengths, &accumColor, &accumThroughput, &gen, &mis, &rayFlags);

    if(isDeadRay(rayFlags))
      break;
  }

  kernel_HitEnvironment(tid, &rayFlags, &rayDirAndFar, &mis, &accumThroughput,
                        &accumColor);

  kernel_ContributeToImage(tid, &rayFlags, channels, &accumColor, &gen, m_packedXY.data(), &wavelengths, out_color);
}

void Integrator::PathTraceFromInputRays(uint tid, uint channels, const RayPosAndW* in_rayPosAndNear, const RayDirAndT* in_rayDirAndFar, float* out_color)
{
  float4 accumColor, accumThroughput;
  float4 rayPosAndNear, rayDirAndFar;
  float4 wavelengths;
  RandomGen gen; 
  MisData   mis;
  uint      rayFlags;
  kernel_InitEyeRayFromInput(tid, in_rayPosAndNear, in_rayDirAndFar, 
                             &rayPosAndNear, &rayDirAndFar, &accumColor, &accumThroughput, &gen, &rayFlags, &mis, &wavelengths);
  
  for(uint depth = 0; depth < m_traceDepth; depth++) 
  {
    float4 shadeColor, hitPart1, hitPart2, hitPart3;
    uint instId;
    kernel_RayTrace2(tid, depth, &rayPosAndNear, &rayDirAndFar, &hitPart1, &hitPart2, &hitPart3, &instId, &rayFlags);
    if(isDeadRay(rayFlags))
      break;
    
    kernel_SampleLightSource(tid, &rayPosAndNear, &rayDirAndFar, &wavelengths, &hitPart1, &hitPart2, &hitPart3, &rayFlags, 
                             depth, &gen, &shadeColor);

    kernel_NextBounce(tid, depth, &hitPart1, &hitPart2, &hitPart3, &instId, &shadeColor, &rayPosAndNear, &rayDirAndFar,
                      &wavelengths, &accumColor, &accumThroughput, &gen, &mis, &rayFlags);

    if(isDeadRay(rayFlags))
      break;
  }

  kernel_HitEnvironment(tid, &rayFlags, &rayDirAndFar, &mis, &accumThroughput,
                        &accumColor);
  
  //////////////////////////////////////////////////// same as for PathTrace

  kernel_CopyColorToOutput(tid, channels, &accumColor, &gen, out_color);
}
