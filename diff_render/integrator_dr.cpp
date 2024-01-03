#include "integrator_dr.h"
#include "utils.h"

#include <chrono>
#include <string>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

extern double __enzyme_autodiff(void*, ...);
int enzyme_const, enzyme_dup, enzyme_out;

float Loss(IntegratorDR* __restrict__ pIntegrator,
           const float*  __restrict__ a_refImg,
                 float*  __restrict__ out_color,
           const float*  __restrict__ a_data, 
           const uint* in_pakedXY, 
           uint tid, uint channels, uint pitch)
{
  const uint XY = in_pakedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;

  pIntegrator->tex_data = a_data; // hack (should it work?)

  float4 colorRend = pIntegrator->PathTrace(tid, channels, out_color);

  float4 colorRef  = float4(a_refImg[(y*pitch+x)*channels + 0], 
                            a_refImg[(y*pitch+x)*channels + 1], 
                            a_refImg[(y*pitch+x)*channels + 2], 0.0f);

  float4 diff = colorRend - colorRef;

  return LiteMath::dot3(diff, diff);
}                     


void IntegratorDR::PathTraceDR(uint tid, uint channels, float* out_color, uint a_passNum,
                               const float* a_refImg, const float* a_data, float* a_dataGrad, size_t a_gradSize)
{
  memset(a_dataGrad, 0, sizeof(float)*a_gradSize);

  // init separate gradient for each thread
  //
  //std::vector<float> grads[MAXTHREADS];
  //for(int i=0;i<MAXTHREADS;i++)
  //   std::fill(grads[i].begin(), grads[i].end(), 0.0f);

  double avgLoss = 0.0;

  ConsoleProgressBar progress(tid);
  progress.Start();
  auto start = std::chrono::high_resolution_clock::now();
  //#ifndef _DEBUG
  //#pragma omp parallel for default(shared) // num_threads(MAXTHREADS)
  //#endif
  for (int i = 0; i < int(tid); ++i) {
    float lossVal = 0.0f;
    for (int j = 0; j < int(a_passNum); ++j) {
      //PathTrace(uint(i), channels, out_color);
      lossVal += Loss(this, 
                      a_refImg, 
                      out_color, 
                      a_data,
                      m_packedXY.data(), 
                      uint(i), channels, m_winWidth);

      __enzyme_autodiff((void*)Loss, 
                         enzyme_const, this,
                         enzyme_const, a_refImg,
                         enzyme_const, out_color,
                         enzyme_dup,   a_data, a_dataGrad,
                         enzyme_const, m_packedXY.data(),
                         enzyme_const, uint(i),
                         enzyme_const, channels,
                         enzyme_const, m_winWidth);
    }
    avgLoss += double(lossVal)/double(a_passNum);
    progress.Update();
  }
  progress.Done();
  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;

  // accumulate gradient from different threads (parallel reduction/hist)
  //

  //for(int i=0;i<MAXTHREADS;i++) 
  //  for(size_t j=0;j<a_gradSize; j++)
  //    a_dataGrad[j] += grads[i][j];

  avgLoss /= double(m_winWidth*m_winHeight);
  std::cout << "avgLoss = " << avgLoss << std::endl;
  
  std::ofstream fout("grad.txt");
  for(size_t i=0; i<a_gradSize; i++)
    fout << a_dataGrad[i]/float(a_passNum) << std::endl;
  fout.close();
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

float4 IntegratorDR::HydraTex2DFetch(uint texId, float2 a_uv)
{
  if(texId == 1 && tex_data != nullptr) 
  {
    const float m_fw = 256.0f;
    const float m_fh = 256.0f;
    const int tex_width  = 256;
    const int tex_height = 256;
    
    float ffx = a_uv.x * m_fw - 0.5f; // a_uv should not be very large, so that the float does not overflow later. 
    float ffy = a_uv.y * m_fh - 0.5f; // This is left to the responsibility of the top level.
    
    if (ffx < 0) ffx = 0.0f;
    if (ffy < 0) ffy = 0.0f;
    
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

    //return float4(1.0, 0.1, 0.1,1);
  }
  else
    return m_textures[texId]->sample(a_uv);
}