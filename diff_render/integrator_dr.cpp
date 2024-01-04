#include "integrator_dr.h"
#include "utils.h"

#include "include/cmaterial.h"
#include "include/cmat_gltf.h"
#include "include/cmat_conductor.h"
#include "include/cmat_glass.h"
#include "include/cmat_diffuse.h"
#include "include/cmat_plastic.h"

#include <chrono>
#include <string>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void IntegratorDR::LoadSceneEnd()
{
  m_matNonDiff.resize(m_materials.size());
  //m_matDiff.resize(m_materials.size());

  for(size_t matId=0;matId<m_materials.size();matId++) 
  {
    m_matNonDiff[matId].lambertTexId = as_uint(m_materials[matId].data[GLTF_UINT_TEXID0]);
    //m_matDiff[matId].color           = m_materials[matId].colors[GLTF_COLOR_BASE];
  }
}

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

  //const uint   texId     = as_uint(m_materials[matId].data[GLTF_UINT_TEXID0]);
  const uint   texId     = m_matNonDiff[matId].lambertTexId;

  const float2 texCoordT = mulRows2x4(m_materials[matId].row0[0], m_materials[matId].row1[0], hitTexCoord);
  const float4 texColor  = HydraTex2DFetch(texId, texCoordT, a_data); 
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


extern double __enzyme_autodiff(void*, ...);
int enzyme_const, enzyme_dup, enzyme_out;

float PixelLoss(IntegratorDR* __restrict__ pIntegrator,
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

  const uint yRef = pIntegrator->m_winHeight - y - 1; // in input images and when load data from HDD y has different direction

  //float4 colorRend = pIntegrator->PathTrace(tid, channels, out_color, a_data);
  float4 colorRend = pIntegrator->CastRayDR(tid, channels, out_color, a_data);

  float4 colorRef  = float4(a_refImg[(yRef*pitch+x)*channels + 0], 
                            a_refImg[(yRef*pitch+x)*channels + 1], 
                            a_refImg[(yRef*pitch+x)*channels + 2], 0.0f);

  float4 diff = colorRend - colorRef;
  float loss = LiteMath::dot3(diff, diff);
  (*outLoss) = loss;
  return loss;
}                     


float IntegratorDR::PathTraceDR(uint tid, uint channels, float* out_color, uint a_passNum,
                                const float* a_refImg, const float* a_data, float* a_dataGrad, size_t a_gradSize)
{
  memset(a_dataGrad, 0, sizeof(float)*a_gradSize);

  // init separate gradient for each thread
  //
  //std::vector<float> grads[MAXTHREADS];
  //for(int i=0;i<MAXTHREADS;i++)
  //   std::fill(grads[i].begin(), grads[i].end(), 0.0f);

  //double avgLoss = 0.0;
  auto start = std::chrono::high_resolution_clock::now();
  //#ifndef _DEBUG
  //#pragma omp parallel for default(shared) // num_threads(MAXTHREADS)
  //#endif

  float avgLoss = 0.0f;
  
  if(m_gradMode != 0)
  {
    for (int i = 0; i < int(tid); ++i) {
      float lossVal = 0.0f;
      __enzyme_autodiff((void*)PixelLoss, 
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
      float lossVal = PixelLoss(this, a_refImg, out_color, a_data, m_packedXY.data(),
                                uint(i), channels, m_winWidth, &lossVal);
      avgLoss += float(lossVal)/float(a_passNum);
    }
  }

  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;

  // accumulate gradient from different threads (parallel reduction/hist)
  //

  //for(int i=0;i<MAXTHREADS;i++) 
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

float4 IntegratorDR::HydraTex2DFetch(uint texId, float2 a_uv, const float* tex_data)
{
  if(texId == 1 && tex_data != nullptr && m_gradMode != 0) 
  {
    const float m_fw = 256.0f;
    const float m_fh = 256.0f;
    const int tex_width  = 256;
    const int tex_height = 256;
    
    float ffx = a_uv.x * m_fw - 0.5f; // a_uv should not be very large, so that the float does not overflow later. 
    float ffy = a_uv.y * m_fh - 0.5f; // This is left to the responsibility of the top level.
    
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
    
    const float4 f1    = float4(tex_data[offsets.x*4+0], tex_data[offsets.x*4+1], tex_data[offsets.x*4+2], tex_data[offsets.x*4+3]);
    const float4 f2    = float4(tex_data[offsets.y*4+0], tex_data[offsets.y*4+1], tex_data[offsets.y*4+2], tex_data[offsets.y*4+3]);
    const float4 f3    = float4(tex_data[offsets.z*4+0], tex_data[offsets.z*4+1], tex_data[offsets.z*4+2], tex_data[offsets.z*4+3]);
    const float4 f4    = float4(tex_data[offsets.w*4+0], tex_data[offsets.w*4+1], tex_data[offsets.w*4+2], tex_data[offsets.w*4+3]);
    // Calculate the weighted sum of pixels (for each color channel)
    //
    const float outr = f1.x * w1 + f2.x * w2 + f3.x * w3 + f4.x * w4;
    const float outg = f1.y * w1 + f2.y * w2 + f3.y * w3 + f4.y * w4;
    const float outb = f1.z * w1 + f2.z * w2 + f3.z * w3 + f4.z * w4;
    const float outa = f1.w * w1 + f2.w * w2 + f3.w * w3 + f4.w * w4;
    
    return float4(outr, outg, outb, outa);
  }
  else
    return m_textures[texId]->sample(a_uv);
    
}