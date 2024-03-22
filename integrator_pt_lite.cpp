#include "integrator_pt.h"
#include "include/crandom.h"

#include <chrono>
#include <string>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

void Integrator::kernel_InitRGen(uint tid, RandomGen* gen) 
{
  *gen = m_randomGens[tid];
  const float2 rands2 = rndFloat2_Pseudo(gen);
  const float4 rands4 = rndFloat4_Pseudo(gen);
}

void Integrator::PathTraceLite(uint tid, uint channels, float* out_color)
{
  float4 rayPosAndNear, rayDirAndFar;
  RandomGen gen;
  kernel_InitRGen(tid, &gen);
  {
    const uint XY = m_packedXY[tid];
    const uint x  = (XY & 0x0000FFFF);
    const uint y  = (XY & 0xFFFF0000) >> 16;
  
    const float4 pixelOffsets = GetRandomNumbersLens(tid, &gen);
  
    const float xCoordNormalized = (float(x) + pixelOffsets.x)/float(m_winWidth);
    const float yCoordNormalized = (float(y) + pixelOffsets.y)/float(m_winHeight);
  
    float3 rayDir = EyeRayDirNormalized(xCoordNormalized, yCoordNormalized, m_projInv);
    float3 rayPos = float3(0,0,0);
  
    transform_ray3f(m_worldViewInv, &rayPos, &rayDir);
  
    rayPosAndNear = to_float4(rayPos, 0.0f);
    rayDirAndFar  = to_float4(rayDir, FLT_MAX);
  }

  float4 accumColor      = make_float4(0,0,0,0);
  float4 accumThroughput = make_float4(1,1,1,1);
  float4 wavelengths     = float4(0.0f);
  MisData   mis          = makeInitialMisData();
  float     time         = 0.0f;
  uint      rayFlags     = 0;

  for(uint depth = 0; depth < m_traceDepth; depth++) 
  {
    float4 shadeColor, hitPart1, hitPart2, hitPart3;
    uint   instId;

    //kernel_RayTrace2(tid, depth, &rayPosAndNear, &rayDirAndFar, &time, 
    //                 &hitPart1, &hitPart2, &hitPart3, &instId, &rayFlags);
    {
      const CRT_Hit hit = m_pAccelStruct->RayQuery_NearestHit(rayPosAndNear, rayDirAndFar, 0.0f);

      if(hit.geomId != uint32_t(-1))
      {
        const float2 uv     = float2(hit.coords[0], hit.coords[1]);        
        const float3 hitPos = to_float3(rayPosAndNear) + hit.t * (1.f - 1e-6f) * to_float3(rayDirAndFar);

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
        hitNorm = mul3x3(m_normMatrices[hit.instId], hitNorm);
        hitTang = mul3x3(m_normMatrices[hit.instId], hitTang);
    
        hitNorm = normalize(hitNorm);
        hitTang = normalize(hitTang);
        
        const float flipNorm = dot(to_float3(rayDirAndFar), hitNorm) > 0.001f ? -1.0f : 1.0f; // beware of transparent materials which use normal sign to identity "inside/outside" glass for example
        hitNorm              = flipNorm * hitNorm;
        hitTang              = flipNorm * hitTang; // do we need this ??
    
        if (flipNorm < 0.0f) rayFlags |=  RAY_FLAG_HAS_INV_NORMAL;
        else                 rayFlags &= ~RAY_FLAG_HAS_INV_NORMAL;
        
        const uint midOriginal = m_matIdByPrimId[m_matIdOffsets[hit.geomId] + hit.primId];
        const uint midRemaped  = RemapMaterialId(midOriginal, hit.instId);
    
        rayFlags = packMatId(rayFlags, midRemaped);
        hitPart1 = to_float4(hitPos,  hitTexCoord.x); 
        hitPart2 = to_float4(hitNorm, hitTexCoord.y);
        hitPart3 = to_float4(hitTang, hit.t);
        instId   = hit.instId;
      }
      else
      {
        const uint flagsToAdd = (depth == 0) ? (RAY_FLAG_PRIME_RAY_MISS | RAY_FLAG_IS_DEAD | RAY_FLAG_OUT_OF_SCENE) : (RAY_FLAG_IS_DEAD | RAY_FLAG_OUT_OF_SCENE);
        rayFlags              = rayFlags | flagsToAdd;
      }
    }

    if(isDeadRay(rayFlags))
      break;
    
    //kernel_SampleLightSource(tid, &rayPosAndNear, &rayDirAndFar, &wavelengths, &hitPart1, &hitPart2, &hitPart3, &rayFlags, &time,
    //                         depth, &gen, &shadeColor);
    {
      const uint32_t matId = extractMatId(rayFlags);
      const float3 ray_dir = to_float3(rayDirAndFar);
      
      const float4 data1  = hitPart1;
      const float4 data2  = hitPart2;
      const float4 lambda = wavelengths;
    
      SurfaceHit hit;
      hit.pos  = to_float3(data1);
      hit.norm = to_float3(data2);
      hit.tang = to_float3(hitPart3);
      hit.uv   = float2(data1.w, data2.w);
      
      const int lightId   = 0;
      const float2 rands2 = rndFloat2_Pseudo(&gen);
    
      const LightSample lSam = areaLightSampleRev  (m_lights.data() + lightId, rands2);
      const float  hitDist   = std::sqrt(dot(hit.pos - lSam.pos, hit.pos - lSam.pos));
    
      const float3 shadowRayDir = normalize(lSam.pos - hit.pos); // explicitSam.direction;
      const float3 shadowRayPos = hit.pos + hit.norm * std::max(maxcomp(hit.pos), 1.0f)*5e-6f; // TODO: see Ray Tracing Gems, also use flatNormal for offset
    
      const bool   inIllumArea  = (dot(shadowRayDir, lSam.norm) < 0.0f) || lSam.isOmni || lSam.hasIES;
      const bool   needShade    = inIllumArea && !m_pAccelStruct->RayQuery_AnyHit(to_float4(shadowRayPos, 0.0f), to_float4(shadowRayDir, hitDist*0.9995f), 0.0f); /// (!!!) expression-way, RT pipeline bug work around, if change check test_213
    
      if(needShade) /// (!!!) expression-way to compute 'needShade', RT pipeline bug work around, if change check test_213
      {
        //const BsdfEval bsdfV    = MaterialEval(matId, lambda, shadowRayDir, (-1.0f)*ray_dir, hit.norm, hit.tang, hit.uv);
        const float4 color      = m_materials[matId].colors[GLTF_COLOR_BASE];
        const float4 bsdfVval   = lambertEvalBSDF(shadowRayDir, (-1.0f)*ray_dir, hit.norm)*color;
        const float bsdfVpdf    = lambertEvalPDF (shadowRayDir, (-1.0f)*ray_dir, hit.norm);
        float cosThetaOut       = std::max(dot(shadowRayDir, hit.norm), 0.0f);
        
        //float      lgtPdfW      = LightEvalPDF(lightId, shadowRayPos, shadowRayDir, lSam.pos, lSam.norm, lSam.pdf);
        const float hitDist   = length(shadowRayPos - lSam.pos);
        const float cosVal    = std::abs(dot(shadowRayDir, -1.0f*lSam.norm));
        const float lgtPdfW   = PdfAtoW(m_lights[0].pdfA, hitDist, cosVal);
        const float misWeight = misWeightHeuristic(lgtPdfW, bsdfVpdf);        
        
        const float4 lightColor = m_lights[0].intensity*m_lights[0].mult; //LightIntensity(lightId, lambda, shadowRayPos, shadowRayDir);
        shadeColor = (lightColor * bsdfVval / lgtPdfW) * cosThetaOut * misWeight;
      }
      else
        shadeColor = float4(0.0f, 0.0f, 0.0f, 0.0f);
    }
    
    kernel_NextBounce(tid, depth, &hitPart1, &hitPart2, &hitPart3, &instId, &shadeColor,
                      &rayPosAndNear, &rayDirAndFar, &wavelengths, &accumColor, &accumThroughput, &gen, &mis, &rayFlags);

    if(isDeadRay(rayFlags))
      break;
  }

  //kernel_HitEnvironment(tid, &rayFlags, &rayDirAndFar, &mis, &accumThroughput,
  //                      &accumColor);
  
  //kernel_ContributeToImage(tid, &rayFlags, channels, &accumColor, &gen, m_packedXY.data(), &wavelengths, out_color);
  {
    m_randomGens[tid] = gen;
    const uint XY = m_packedXY[tid];
    const uint x  = (XY & 0x0000FFFF);
    const uint y  = (XY & 0xFFFF0000) >> 16;
  
    float4 colorRes = m_exposureMult * accumColor;
    out_color[(y*m_winWidth+x)*channels + 0] += colorRes.x;
    out_color[(y*m_winWidth+x)*channels + 1] += colorRes.y;
    out_color[(y*m_winWidth+x)*channels + 2] += colorRes.z;
  }
}