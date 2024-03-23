#include "integrator_dr.h"

#include "include/cmaterial.h"
#include "include/cmat_gltf.h"
#include "include/cmat_conductor.h"
#include "include/cmat_glass.h"
#include "include/cmat_diffuse.h"
#include "include/cmat_plastic.h"

#include <chrono>
#include <string>
#include <omp.h>

#include "utils.h"
#include "imageutils.h"

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

void IntegratorDR::LoadSceneEnd()
{
  m_texAddressTable.resize(m_textures.size());
  for(auto& texInfo : m_texAddressTable)
    texInfo = {size_t(-1),0,0,0,0};

  m_gradSize = 0;
}

std::pair<size_t, size_t> IntegratorDR::PutDiffTex2D(uint32_t texId, uint32_t width, uint32_t height, uint32_t channels)
{
  if(texId >= m_texAddressTable.size())
  {
    std::cout << "[IntegratorDR::PutDiffTex2D]: bad tex id = " << texId << std::endl;
    return std::make_pair(size_t(-1), 0);
  }

  m_texAddressTable[texId].offset   = m_gradSize;
  m_texAddressTable[texId].width    = width;
  m_texAddressTable[texId].height   = height;
  m_texAddressTable[texId].channels = channels;
  m_texAddressTable[texId].fwidth   = float(width);
  m_texAddressTable[texId].fheight  = float(height);

  size_t oldOffset = m_gradSize;
  size_t currSize  = size_t(width)*size_t(height)*size_t(channels);
  
  m_gradSize += currSize;
  return std::make_pair(oldOffset, currSize);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline int4 bilinearOffsets(const float ffx, const float ffy, const int w, const int h)
{
	const int sx = (ffx > 0.0f) ? 1 : -1;
	const int sy = (ffy > 0.0f) ? 1 : -1;

	const int px = (int)(ffx);
	const int py = (int)(ffy);

	int px_w0, px_w1, py_w0, py_w1;
  // wrap
	{
		px_w0 = px        % w;
		px_w1 = (px + sx) % w;

		px_w0 = (px_w0 < 0) ? px_w0 + w : px_w0;
		px_w1 = (px_w1 < 0) ? px_w1 + w : px_w1;
	}

  // wrap
	{
		py_w0 = py        % h;
		py_w1 = (py + sy) % h;

		py_w0 = (py_w0 < 0) ? py_w0 + h : py_w0;
		py_w1 = (py_w1 < 0) ? py_w1 + h : py_w1;
	}

	const int offset0 = py_w0*w + px_w0;
	const int offset1 = py_w0*w + px_w1;
	const int offset2 = py_w1*w + px_w0;
	const int offset3 = py_w1*w + px_w1;

	return int4(offset0, offset1, offset2, offset3);
}

float4 IntegratorDR::Tex2DFetchAD(uint texId, float2 a_uv, const float* tex_data)
{
  const auto info = m_texAddressTable[texId];

  if(info.offset != size_t(-1) && tex_data != nullptr && m_gradMode != 0) 
  {
    const float m_fw     = info.fwidth;
    const float m_fh     = info.fheight;
    const int tex_width  = info.width;
    const int tex_height = info.height;
    
    float ffx = a_uv.x * m_fw - 0.5f; // a_texCoord should not be very large, so that the float does not overflow later. 
    float ffy = a_uv.y * m_fh - 0.5f; // This is left to the responsibility of the top level.
    
    auto sampler = m_textures[texId]->sampler();

    if ((sampler.addressU == Sampler::AddressMode::CLAMP) != 0 && ffx < 0) ffx = 0.0f;
    if ((sampler.addressV == Sampler::AddressMode::CLAMP) != 0 && ffy < 0) ffy = 0.0f;
    
    // Calculate the weights for each pixel
    //
    const int   px = (int)(ffx);
    const int   py = (int)(ffy);
    
    const float fx  = std::abs(ffx - (float)px);
    const float fy  = std::abs(ffy - (float)py);
    const float fx1 = 1.0f - fx;
    const float fy1 = 1.0f - fy;
    
    const float w1 = fx1 * fy1;
    const float w2 = fx  * fy1;
    const float w3 = fx1 * fy;
    const float w4 = fx  * fy;
    
    const int4 offsets = bilinearOffsets(ffx, ffy, tex_width, tex_height);

    // Calculate the weighted sum of pixels (for each color channel)
    //
    if(info.channels == 4)
    {
      const float4 f1    = float4(tex_data[info.offset+offsets.x*4+0], tex_data[info.offset+offsets.x*4+1], tex_data[info.offset+offsets.x*4+2], tex_data[info.offset+offsets.x*4+3]);
      const float4 f2    = float4(tex_data[info.offset+offsets.y*4+0], tex_data[info.offset+offsets.y*4+1], tex_data[info.offset+offsets.y*4+2], tex_data[info.offset+offsets.y*4+3]);
      const float4 f3    = float4(tex_data[info.offset+offsets.z*4+0], tex_data[info.offset+offsets.z*4+1], tex_data[info.offset+offsets.z*4+2], tex_data[info.offset+offsets.z*4+3]);
      const float4 f4    = float4(tex_data[info.offset+offsets.w*4+0], tex_data[info.offset+offsets.w*4+1], tex_data[info.offset+offsets.w*4+2], tex_data[info.offset+offsets.w*4+3]);
  
      const float outr = f1.x * w1 + f2.x * w2 + f3.x * w3 + f4.x * w4;
      const float outg = f1.y * w1 + f2.y * w2 + f3.y * w3 + f4.y * w4;
      const float outb = f1.z * w1 + f2.z * w2 + f3.z * w3 + f4.z * w4;
      const float outa = f1.w * w1 + f2.w * w2 + f3.w * w3 + f4.w * w4;
      
      return float4(outr, outg, outb, outa);
    }
    else
    {
      const float f1 = tex_data[info.offset+offsets.x];
      const float f2 = tex_data[info.offset+offsets.y];
      const float f3 = tex_data[info.offset+offsets.z];
      const float f4 = tex_data[info.offset+offsets.w];
  
      const float outVal = f1 * w1 + f2 * w2 + f3 * w3 + f4 * w4;
      
      return float4(outVal, outVal, outVal, outVal);
    }
  }
  else
    return m_textures[texId]->sample(a_uv);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void IntegratorDR::kernel_InitEyeRay(uint tid, const uint* packedXY, float4* rayPosAndNear, float4* rayDirAndFar, const float* a_data) // (tid,tidX,tidY,tidZ) are SPECIAL PREDEFINED NAMES!!!
{
  if(tid >= m_maxThreadId)
    return;
  const uint XY = packedXY[tid];

  const uint x = (XY & 0x0000FFFF);
  const uint y = (XY & 0xFFFF0000) >> 16;

  float3 rayDir = EyeRayDirNormalized((float(x)+0.5f)/float(m_winWidth), (float(y)+0.5f)/float(m_winHeight), m_projInv);
  float3 rayPos = float3(0,0,0);

  transform_ray3f(m_worldViewInv, 
                  &rayPos, &rayDir);
  
  *rayPosAndNear = to_float4(rayPos, 0.0f);
  *rayDirAndFar  = to_float4(rayDir, FLT_MAX);
}

bool IntegratorDR::kernel_RayTrace(uint tid, const float4* rayPosAndNear, float4* rayDirAndFar,
                                   Lite_Hit* out_hit, float2* out_bars, const float* a_data)
{
  if(tid >= m_maxThreadId)
    return false;
  const float4 rayPos = *rayPosAndNear;
  const float4 rayDir = *rayDirAndFar ;

  CRT_Hit hit = m_pAccelStruct->RayQuery_NearestHit(rayPos, rayDir);
  
  Lite_Hit res;
  res.primId = hit.primId;
  res.instId = hit.instId;
  res.geomId = hit.geomId;
  res.t      = hit.t;

  float2 baricentrics = float2(hit.coords[0], hit.coords[1]);
 
  *out_hit  = res;
  *out_bars = baricentrics;
  return (res.primId != -1);
}


void IntegratorDR::kernel_CalcRayColor(uint tid, const Lite_Hit* in_hit, const float2* bars, float4* finalColor, const uint* in_pakedXY, float* out_color, const float* a_data)
{ 
  if(tid >= m_maxThreadId)
    return;

  const Lite_Hit hit = *in_hit;
  if(hit.geomId == -1)
  {
    out_color[tid] = 0;
    return;
  }

  const uint32_t matId  = m_matIdByPrimId[m_matIdOffsets[hit.geomId] + hit.primId];
  const float4 mdata    = m_materials[matId].colors[GLTF_COLOR_BASE];
  const float2 uv       = *bars;

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

  const uint   texId     = m_materials[matId].texid[0];
  const float2 texCoordT = mulRows2x4(m_materials[matId].row0[0], m_materials[matId].row1[0], hitTexCoord);
  const float4 texColor  = Tex2DFetchAD(texId, texCoordT, a_data); 
  const float3 color     = mdata.w > 0.0f ? clamp(float3(mdata.w,mdata.w,mdata.w), 0.0f, 1.0f) : to_float3(mdata*texColor);

  const uint XY = in_pakedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;

  (*finalColor) = to_float4(color, 0);

  out_color[(y*m_winWidth+x)*4 + 0] = color.x;
  out_color[(y*m_winWidth+x)*4 + 1] = color.y;
  out_color[(y*m_winWidth+x)*4 + 2] = color.z;
  out_color[(y*m_winWidth+x)*4 + 3] = 0.0f;
}


float4 IntegratorDR::CastRayDR(uint tid, uint channels, float* out_color, const float* a_data)
{
  float4 rayPosAndNear, rayDirAndFar;
  kernel_InitEyeRay(tid, m_packedXY.data(), &rayPosAndNear, &rayDirAndFar, a_data);

  Lite_Hit hit; 
  float2   baricentrics; 
  if(!kernel_RayTrace(tid, &rayPosAndNear, &rayDirAndFar, &hit, &baricentrics, a_data))
    return float4(0,0,0,0);
  
  float4 finalColor;
  kernel_CalcRayColor(tid, &hit, &baricentrics, &finalColor, m_packedXY.data(), out_color, a_data);
  return finalColor;
}

extern float __enzyme_autodiff(void*, ...);
int enzyme_const, enzyme_dup, enzyme_out;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double RegLossImage1D(size_t a_size, const float* data)
{
  double summ = 0.0f;

  for(size_t i=1;i<a_size-1;i++) {
    float diffLeft  = data[i] - data[i-1];
    float diffRight = data[i] - data[i+1];
    summ += double(diffLeft*diffLeft + diffRight*diffRight);
  }

  return summ/double(a_size);
}

double RegLossImage2D(int w, int h, const float* data)
{
  double summ = 0.0f;
  std::vector<double> lines(h);

  for(int y=1;y<h-1;y++) {
    lines[y] = 0.0;
    for(int x=1;x<w-1;x++) {
      float diffTop    = data[y*w+x] - data[(y+1)*w+x];
      float diffBottom = data[y*w+x] - data[(y-1)*w+x];
      float diffLeft   = data[y*w+x] - data[y*w+x-1];
      float diffRight  = data[y*w+x] - data[y*w+x+1];
      lines[y] += std::sqrt(double(diffLeft*diffLeft + diffRight*diffRight + diffTop*diffTop + diffBottom*diffBottom));
    }
    summ += lines[y];
  }
 
  return summ;
}

using LiteMath::dot3;

double RegLossImage2D4f(int w, int h, const float* data)
{
  double summ = 0.0f;
  std::vector<double> lines(h);

  for(int y=1;y<h-1;y++) {
    lines[y] = 0.0;
    for(int x=1;x<w-1;x++) {
      float4 p0 = float4(data[(y*w+x)*4+0],     data[(y*w+x)*4+1],     data[(y*w+x)*4+2],     data[(y*w+x)*4+3]);
      float4 p1 = float4(data[((y+1)*w+x)*4+0], data[((y+1)*w+x)*4+1], data[((y+1)*w+x)*4+2], data[((y+1)*w+x)*4+3]);
      float4 p2 = float4(data[((y-1)*w+x)*4+0], data[((y-1)*w+x)*4+1], data[((y-1)*w+x)*4+2], data[((y-1)*w+x)*4+3]); 
      float4 p3 = float4(data[(y*w+x-1)*4+0], data[(y*w+x-1)*4+1], data[(y*w+x-1)*4+2], data[(y*w+x-1)*4+3]);
      float4 p4 = float4(data[(y*w+x+1)*4+0], data[(y*w+x+1)*4+1], data[(y*w+x+1)*4+2], data[(y*w+x+1)*4+3]);

      float4 diffTop    = p0 - p1;
      float4 diffBottom = p0 - p2;
      float4 diffLeft   = p0 - p3;
      float4 diffRight  = p0 - p4;

      lines[y] += std::sqrt(double(dot3(diffLeft,diffLeft) + dot3(diffRight,diffRight) + dot3(diffTop,diffTop) + dot3(diffBottom,diffBottom)));
    }
  }
  
  summ = 0.0;
  for(int sy=0; sy < h/2; sy++) {
    int index1 = h/2 + sy;
    int index2 = h/2 - sy;
    if(index1 < h-1)
      summ += lines[index1];
    if(sy!=0 && index2 >= 1)
      summ += lines[index2];
  }
 
  return summ; //double(w*h);
}

void Image1DRegularizer(size_t a_size, const float* data, float* grad)
{
  __enzyme_autodiff((void*)RegLossImage1D, 
                           enzyme_const, a_size,
                           enzyme_dup,   data, grad);

}

void Image2D4fRegularizer(int w, int h, const float* data, float* grad)
{
  __enzyme_autodiff((void*)RegLossImage2D4f, 
                           enzyme_const, w,
                           enzyme_const, h,
                           enzyme_dup,   data, grad);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float PixelLossRT(IntegratorDR* __restrict__ pIntegrator,
                  const float*  __restrict__ a_refImg,
                        float*  __restrict__ out_color,
                  const float*  __restrict__ a_data, 
                  const uint*   __restrict__ in_pakedXY, 
                  uint tid, uint channels, uint pitch,
                  float*  __restrict__       outLoss)
{
  const uint XY = in_pakedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;

  const uint yRef  = pIntegrator->m_winHeight - y - 1; // in input images and when load data from HDD y has different direction
  float4 colorRend = pIntegrator->CastRayDR(tid, channels, out_color, a_data);
  float4 colorRef  = float4(a_refImg[(yRef*pitch+x)*channels + 0], 
                            a_refImg[(yRef*pitch+x)*channels + 1], 
                            a_refImg[(yRef*pitch+x)*channels + 2], 0.0f);

  float4 diff = colorRend - colorRef;
  float loss = LiteMath::dot3(diff, diff);
  (*outLoss) = loss;
  return loss;
}                      

float IntegratorDR::RayTraceDR(uint tid, uint channels, float* out_color, uint a_passNum,
                                const float* a_refImg, const float* a_data, float* a_dataGrad, size_t a_gradSize)
{
  memset(a_dataGrad, 0, sizeof(float)*a_gradSize);

  // init separate gradient for each thread
  //
  //std::vector<float> grads[MAXTHREADS_CPU];
  //for(int i=0;i<MAXTHREADS_CPU;i++)
  //   std::fill(grads[i].begin(), grads[i].end(), 0.0f);

  //double avgLoss = 0.0;
  auto start = std::chrono::high_resolution_clock::now();
  //#ifndef _DEBUG
  //#pragma omp parallel for default(shared) // num_threads(MAXTHREADS_CPU)
  //#endif

  float avgLoss = 0.0f;
  
  if(m_gradMode != 0)
  {
    for (int i = 0; i < int(tid); ++i) {
      float lossVal = 0.0f;
      __enzyme_autodiff((void*)PixelLossRT, 
                         enzyme_const, this,
                         enzyme_const, a_refImg,
                         enzyme_const, out_color,
                         enzyme_dup,   a_data, a_dataGrad,
                         enzyme_const, m_packedXY.data(),
                         enzyme_const, uint(i),
                         enzyme_const, channels,
                         enzyme_const, m_winWidth,
                         enzyme_const, &lossVal);
      avgLoss += float(lossVal)/float(a_passNum);
    }
  }
  else
  {
    for (int i = 0; i < int(tid); ++i) {
      float lossVal = PixelLossRT(this, a_refImg, out_color, a_data, m_packedXY.data(),
                                uint(i), channels, m_winWidth, &lossVal);
      avgLoss += float(lossVal)/float(a_passNum);
    }
  }

  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;

  // accumulate gradient from different threads (parallel reduction/hist)
  //

  //for(int i=0;i<MAXTHREADS_CPU;i++) 
  //  for(size_t j=0;j<a_gradSize; j++)
  //    a_dataGrad[j] += grads[i][j];

  //avgLoss /= float(m_winWidth*m_winHeight);
  //std::cout << "avgLoss = " << avgLoss << std::endl;
  
  //std::ofstream fout("z_grad.txt");
  //for(size_t i=0; i<a_gradSize; i++)
  //  fout << a_dataGrad[i]/float(a_passNum) << std::endl;
  //fout.close();

  return avgLoss;
}

float4 IntegratorDR::PathTraceReplay(uint tid, uint channels, uint cpuThreadId, float* out_color, 
                                     const float* drands, const float* dparams)
{
  float4  accumColor      = float4(0.0f);
  float4  accumThroughput = float4(1.0f);
  MisData mis             = makeInitialMisData();
  uint    rayFlags = 0;

  RandomGen gen = m_randomGens[RandomGenId(tid)];

  float4 rayPosAndNear, rayDirAndFar;
  float4 wavelengths;
  float  time;
  //kernel_InitEyeRay2(tid, &rayPosAndNear, &rayDirAndFar, &wavelengths, &accumColor, &accumThroughput, &gen, &rayFlags, &mis, &time);
  {
    const uint XY = m_packedXY[tid];
    const uint x  = (XY & 0x0000FFFF);
    const uint y  = (XY & 0xFFFF0000) >> 16;
  
    const float4 pixelOffsets = float4(drands[0], drands[1], drands[2], drands[3]);
    const float2 wt           = float2(drands[4], drands[5]);
  
    const float xCoordNormalized = (float(x) + pixelOffsets.x)/float(m_winWidth);
    const float yCoordNormalized = (float(y) + pixelOffsets.y)/float(m_winHeight);
  
    float3 rayDir = EyeRayDirNormalized(xCoordNormalized, yCoordNormalized, m_projInv);
    float3 rayPos = float3(0,0,0);
    
    if(KSPEC_SPECTRAL_RENDERING !=0 && m_spectral_mode != 0)
      wavelengths = SampleWavelengths(wt.x, LAMBDA_MIN, LAMBDA_MAX);
    else
      wavelengths = float4(0.0f);
  
    time = wt.y;
    transform_ray3f(m_worldViewInv, &rayPos, &rayDir);
  
    rayPosAndNear = to_float4(rayPos, 0.0f);
    rayDirAndFar  = to_float4(rayDir, FLT_MAX);
  }

  for(uint bounce = 0; bounce < m_traceDepth; bounce++) 
  {
    float4 shadeColor = float4(0.0f, 0.0f, 0.0f, 0.0f);
    float4 hitPart1, hitPart2, hitPart3;
    uint   instId;

    //kernel_RayTrace2(tid, bounce, &rayPosAndNear, &rayDirAndFar, &time, 
    //                 &hitPart1, &hitPart2, &hitPart3, &instId, &rayFlags);
    {
      //const CRT_Hit hit = m_recorded[cpuThreadId].perBounce[bounce].hit;
      const CRT_Hit hit = m_pAccelStruct->RayQuery_NearestHit(rayPosAndNear, rayDirAndFar, time);

      if(hit.geomId != uint32_t(-1))
      {
        const float2 uv     = float2(hit.coords[0], hit.coords[1]);
        
        // slightly undershoot the intersection to prevent self-intersection and other bugs
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
        const uint flagsToAdd = (bounce == 0) ? (RAY_FLAG_PRIME_RAY_MISS | RAY_FLAG_IS_DEAD | RAY_FLAG_OUT_OF_SCENE) : (RAY_FLAG_IS_DEAD | RAY_FLAG_OUT_OF_SCENE);
        rayFlags              = rayFlags | flagsToAdd;
      }

    }

    if(isDeadRay(rayFlags))
      break;
    
    //kernel_SampleLightSource(tid, &rayPosAndNear, &rayDirAndFar, &wavelengths, &hitPart1, &hitPart2, &hitPart3, &rayFlags, &time,
    //                         bounce, &gen, &shadeColor);
    {
      const uint32_t matId = extractMatId(rayFlags);
      const float3 ray_dir = to_float3(rayDirAndFar);
      const float4 lambda  = wavelengths;
    
      SurfaceHit hit;
      hit.pos  = to_float3(hitPart1);
      hit.norm = to_float3(hitPart2);
      hit.tang = to_float3(hitPart3);
      hit.uv   = float2(hitPart1.w, hitPart2.w);
      
      const int bounceTmp = int(bounce); 
      const float4 rands = GetRandomNumbersLgts(tid, &gen, bounceTmp); 
      const int lightId  = std::min(int(std::floor(rands.w * float(m_lights.size()))), int(m_lights.size() - 1u));
      RecordLightRndIfNeeded(bounce, rands); 

      if(lightId >= 0) // no lights or invalid light id
      {
        const LightSample lSam = LightSampleRev(lightId, to_float3(rands), hit.pos);
        const float  hitDist   = std::sqrt(dot(hit.pos - lSam.pos, hit.pos - lSam.pos));
      
        const float3 shadowRayDir = normalize(lSam.pos - hit.pos); // explicitSam.direction;
        const float3 shadowRayPos = hit.pos + hit.norm * std::max(maxcomp(hit.pos), 1.0f)*5e-6f; // TODO: see Ray Tracing Gems, also use flatNormal for offset
      
        const bool   inIllumArea  = (dot(shadowRayDir, lSam.norm) < 0.0f) || lSam.isOmni || lSam.hasIES;

        const bool   needShade    = inIllumArea && !m_pAccelStruct->RayQuery_AnyHit(to_float4(shadowRayPos, 0.0f), to_float4(shadowRayDir, hitDist*0.9995f), time); /// (!!!) expression-way, RT pipeline bug work around, if change check test_213
        RecordShadowHitIfNeeded(bounce, needShade);

        if(needShade) /// (!!!) expression-way to compute 'needShade', RT pipeline bug work around, if change check test_213
        {
          const BsdfEval bsdfV    = MaterialEval(matId, lambda, shadowRayDir, (-1.0f)*ray_dir, hit.norm, hit.tang, hit.uv);
          float cosThetaOut       = std::max(dot(shadowRayDir, hit.norm), 0.0f);
          
          float      lgtPdfW      = LightPdfSelectRev(lightId) * LightEvalPDF(lightId, shadowRayPos, shadowRayDir, lSam.pos, lSam.norm, lSam.pdf);
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
      
          const bool isDirectLight = !hasNonSpecular(rayFlags);
          if((m_renderLayer == FB_DIRECT   && !isDirectLight) || 
             (m_renderLayer == FB_INDIRECT && isDirectLight)) // skip some number of bounces if this is set
            misWeight = 0.0f;
            
          
          const float4 lightColor = LightIntensity(lightId, lambda, shadowRayPos, shadowRayDir);
          shadeColor = (lightColor * bsdfV.val / lgtPdfW) * cosThetaOut * misWeight;
        }
      }
    }

    kernel_NextBounce(tid, bounce, &hitPart1, &hitPart2, &hitPart3, &instId, &shadeColor,
                      &rayPosAndNear, &rayDirAndFar, &wavelengths, &accumColor, &accumThroughput, &gen, &mis, &rayFlags);

    if(isDeadRay(rayFlags))
      break;
  }

  kernel_HitEnvironment(tid, &rayFlags, &rayDirAndFar, &mis, &accumThroughput, &accumColor);

  return accumColor;
}

float PixelLossPT(IntegratorDR* __restrict__ pIntegrator,
                  uint tid, uint channels, uint pitch, uint cpuThreadId,
                  const float*  __restrict__ a_refImg,
                        float*  __restrict__ out_color,
                  const uint*   __restrict__ in_pakedXY, 
                  float*        __restrict__ outLoss,
                  const float*  __restrict__ a_drands,
                  const float*  __restrict__ a_dparams)
{
  const uint XY = in_pakedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;

  float4 colorRend = pIntegrator->PathTraceReplay(tid, channels, cpuThreadId, out_color, 
                                                  a_drands, a_dparams);

  const uint yRef  = pIntegrator->m_winHeight - y - 1; // in input images and when load data from HDD y has different direction
  float4 colorRef  = float4(a_refImg[(yRef*pitch+x)*channels + 0], 
                            a_refImg[(yRef*pitch+x)*channels + 1], 
                            a_refImg[(yRef*pitch+x)*channels + 2], 0.0f);

  float4 diff = colorRend - colorRef;
  float loss = LiteMath::dot3(diff, diff);
  (*outLoss) = loss;
  return loss;
}   


float IntegratorDR::PathTraceDR(uint size, uint channels, float* out_color, uint a_passNum,
                                const float* a_refImg, const float* a_data, float* a_dataGrad, size_t a_gradSize)
{
  m_disableImageContrib = 1;
  memset(a_dataGrad, 0, sizeof(float)*a_gradSize);

  // init separate gradient for each thread
  //
  std::vector<float> grads[MAXTHREADS_CPU];
  for(int i=0;i<MAXTHREADS_CPU;i++) {
    grads[i].resize(a_gradSize);
    std::fill(grads[i].begin(), grads[i].end(), 0.0f);
  }

  //double avgLoss = 0.0;
  auto start = std::chrono::high_resolution_clock::now();
  float avgLoss = 0.0f;
  
  if(m_gradMode != 0)
  {
    #ifndef _DEBUG
    #pragma omp parallel for default(shared) num_threads(MAXTHREADS_CPU)
    #endif
    for (int i = 0; i < int(size); ++i) 
    {
      float lossVal = 0.0f;
      for(int passId = 0; passId < int(a_passNum); passId++) 
      {
        // (1) record non differentiable data during common PT 
        //
        auto cpuThreadId = omp_get_thread_num();
        this->m_recorded[cpuThreadId].recordEnabled = true;
        this->PathTrace(i, channels, out_color);
        this->m_recorded[cpuThreadId].recordEnabled = false;
        
        // (2) perform path trace replay and differentiate actual function
        //
        float4 color = this->PathTraceReplay(i, channels, cpuThreadId, out_color, 
                                             this->m_recorded[cpuThreadId].perBounceRands.data(), grads[cpuThreadId].data());
       
        
        const uint XY = m_packedXY[i];
        const uint x  = (XY & 0x0000FFFF);
        const uint y  = (XY & 0xFFFF0000) >> 16;

        out_color[(y*m_winWidth+x)*channels + 0] += color.x;
        out_color[(y*m_winWidth+x)*channels + 1] += color.y;
        out_color[(y*m_winWidth+x)*channels + 2] += color.z;

        //__enzyme_autodiff((void*)PixelLossPT, 
        //                   enzyme_const, this,
        //                   enzyme_const, uint(i),
        //                   enzyme_const, channels,
        //                   enzyme_const, m_winWidth,
        //                   enzyme_const, cpuThreadId,
        //                   enzyme_const, a_refImg,
        //                   enzyme_const, out_color,
        //                   enzyme_const, m_packedXY.data(),
        //                   enzyme_const, &lossVal,
        //                   enzyme_const, this->m_recorded[cpuThreadId].perBounceRands.data(),
        //                   enzyme_dup,   a_data, grads[cpuThreadId].data());

        avgLoss += float(lossVal)/float(a_passNum);
      }
    }
  }
  else
  {
    //for (int i = 0; i < int(tid); ++i) {
    //  float lossVal = PixelLossPT(this, a_refImg, out_color, a_data, m_packedXY.data(),
    //                              uint(i), channels, m_winWidth, &lossVal);
    //  avgLoss += float(lossVal)/float(a_passNum);
    //}
  }
  
  const float normConst = 1.0f/float(a_passNum);
  SaveImage4fToBMP(out_color, m_winWidth, m_winHeight, 4, "z_render1.bmp", normConst, 2.4f);

  diffPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;

  // accumulate gradient from different threads (parallel reduction/hist)
  //
  for(int i=0;i<MAXTHREADS_CPU;i++) 
    for(size_t j=0;j<a_gradSize; j++)
      a_dataGrad[j] += grads[i][j];

  avgLoss /= float(m_winWidth*m_winHeight);

  //std::cout << "avgLoss = " << avgLoss << std::endl;
  //std::cout.flush();
  
  //std::ofstream fout("z_grad.txt");
  //for(size_t i=0; i<a_gradSize; i++)
  //  fout << a_dataGrad[i]/float(a_passNum) << std::endl;
  //fout.close();

  m_disableImageContrib = 0;
  return avgLoss;
}

