#include "integrator_pt.h"
#include <memory>
#include <limits>
#include <omp.h>

#include "utils.h" // for progress bar

class IntegratorKMLT : public Integrator
{
public:

  IntegratorKMLT(int a_maxThreads, int a_spectral_mode, std::vector<uint32_t> a_features) : Integrator(a_maxThreads, a_spectral_mode, a_features)
  {

  }

  float  GetRandomNumbersSpec(uint tid, RandomGen* a_gen) override;
  float4 GetRandomNumbersLens(uint tid, RandomGen* a_gen) override;
  float4 GetRandomNumbersMats(uint tid, RandomGen* a_gen, int a_bounce) override;
  float4 GetRandomNumbersLgts(uint tid, RandomGen* a_gen, int a_bounce) override;
  float  GetRandomNumbersMatB(uint tid, RandomGen* a_gen, int a_bounce, int a_layer) override;

  void   kernel_InitEyeRay2(uint tid, const uint* packedXY, 
                            float4* rayPosAndNear, float4* rayDirAndFar, float4* wavelengths, 
                            float4* accumColor,    float4* accumuThoroughput,
                            RandomGen* gen, uint* rayFlags, MisData* misData) override;

  void   kernel_ContributeToImage(uint tid, const uint* rayFlags, uint channels, const float4* a_accumColor, const RandomGen* gen,
                                  const uint* in_pakedXY, const float4* wavelengths, float* out_color) override;

  void   PathTraceBlock(uint pixelsNum, uint channels, float* out_color, uint a_passNum) override;
  
  float4 F(const std::vector<float>& xVal, uint tid, int*pX, int* pY);

  static constexpr uint LGHT_ID      = 0;
  static constexpr uint MATS_ID      = 4;
  static constexpr uint BLND_ID      = 8;
  static constexpr uint PER_BOUNCE   = 10;
  static constexpr uint BOUNCE_START = 5;
  
  inline uint RandsPerThread() const { return PER_BOUNCE*m_traceDepth + BOUNCE_START; }

  std::vector<float>   m_allRandomNumbers;
  std::vector<float4>  m_allColors;
  std::vector<int2>    m_allXY;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float4 IntegratorKMLT::GetRandomNumbersLens(uint tid, RandomGen* a_gen) 
{ 
  return float4(m_allRandomNumbers[tid*RandsPerThread() + 0],
                m_allRandomNumbers[tid*RandsPerThread() + 1],
                m_allRandomNumbers[tid*RandsPerThread() + 2],
                m_allRandomNumbers[tid*RandsPerThread() + 3]);
}

float  IntegratorKMLT::GetRandomNumbersSpec(uint tid, RandomGen* a_gen) 
{ 
  return m_allRandomNumbers[tid*RandsPerThread() + 4]; 
}

float4 IntegratorKMLT::GetRandomNumbersMats(uint tid, RandomGen* a_gen, int a_bounce) 
{ 
  return float4(m_allRandomNumbers[tid*RandsPerThread() + BOUNCE_START + a_bounce*PER_BOUNCE + MATS_ID + 0],
                m_allRandomNumbers[tid*RandsPerThread() + BOUNCE_START + a_bounce*PER_BOUNCE + MATS_ID + 1],
                m_allRandomNumbers[tid*RandsPerThread() + BOUNCE_START + a_bounce*PER_BOUNCE + MATS_ID + 2],
                m_allRandomNumbers[tid*RandsPerThread() + BOUNCE_START + a_bounce*PER_BOUNCE + MATS_ID + 3]);
}

float4 IntegratorKMLT::GetRandomNumbersLgts(uint tid, RandomGen* a_gen, int a_bounce)
{
  return float4(m_allRandomNumbers[tid*RandsPerThread() + BOUNCE_START + a_bounce*PER_BOUNCE + LGHT_ID + 0],
                m_allRandomNumbers[tid*RandsPerThread() + BOUNCE_START + a_bounce*PER_BOUNCE + LGHT_ID + 1],
                m_allRandomNumbers[tid*RandsPerThread() + BOUNCE_START + a_bounce*PER_BOUNCE + LGHT_ID + 2],
                m_allRandomNumbers[tid*RandsPerThread() + BOUNCE_START + a_bounce*PER_BOUNCE + LGHT_ID + 3]);
}

float IntegratorKMLT::GetRandomNumbersMatB(uint tid, RandomGen* a_gen, int a_bounce, int a_layer) 
{ 
  return m_allRandomNumbers[tid*RandsPerThread() + BOUNCE_START + a_bounce*PER_BOUNCE + BLND_ID + a_layer];
}

void IntegratorKMLT::kernel_InitEyeRay2(uint tid, const uint* packedXY, 
                                       float4* rayPosAndNear, float4* rayDirAndFar, float4* wavelengths, 
                                       float4* accumColor,    float4* accumuThoroughput,
                                       RandomGen* gen, uint* rayFlags, MisData* misData) // 
{
  if(tid >= m_maxThreadId)
    return;

  *accumColor        = make_float4(0,0,0,0);
  *accumuThoroughput = make_float4(1,1,1,1);
  *rayFlags          = 0;
  *misData           = makeInitialMisData();
  
  ///////////////////////////////////////////////////////////////////////////////////// begin change
  const uint rgIndex = tid % uint(m_randomGens.size());
  RandomGen genLocal = m_randomGens[rgIndex];

  const float4 pixelOffsets = GetRandomNumbersLens(tid, &genLocal);

  uint x = uint(pixelOffsets.x*float(m_winWidth));
  uint y = uint(pixelOffsets.y*float(m_winHeight));
  
  if(x >= uint(m_winWidth-1))
    x = uint(m_winWidth-1);
  if(y >= uint(m_winHeight-1))
    y = uint(m_winHeight-1);
  
  m_allXY[tid] = int2(x,y);

  float3 rayDir = EyeRayDirNormalized(pixelOffsets.x, pixelOffsets.y, m_projInv);
  float3 rayPos = float3(0,0,0);
  ///////////////////////////////////////////////////////////////////////////////////// end change

  transform_ray3f(m_worldViewInv, &rayPos, &rayDir);
  
  float tmp = 0.0f;
  if(KSPEC_SPECTRAL_RENDERING !=0 && m_spectral_mode != 0)
  {
    float u = GetRandomNumbersSpec(tid, &genLocal);
    *wavelengths = SampleWavelengths(u, LAMBDA_MIN, LAMBDA_MAX);
    tmp = u;
  }
  else
  {
    const uint32_t sample_sz = sizeof((*wavelengths).M) / sizeof((*wavelengths).M[0]);
    for (uint32_t i = 0; i < sample_sz; ++i) 
      (*wavelengths)[i] = 0.0f;
  }

  RecordPixelRndIfNeeded(float2(pixelOffsets.x, pixelOffsets.y), tmp);
 
  *rayPosAndNear = to_float4(rayPos, 0.0f);
  *rayDirAndFar  = to_float4(rayDir, FLT_MAX);
  *gen           = genLocal;
}

void IntegratorKMLT::kernel_ContributeToImage(uint tid, const uint* rayFlags, uint channels, const float4* a_accumColor, const RandomGen* gen,
                                              const uint* in_pakedXY, const float4* wavelengths, float* out_color)
{
  
  if(tid >= m_maxThreadId) // don't contrubute to image in any "record" mode
    return;
  
  const float4 pixelOffsets = GetRandomNumbersLens(tid, const_cast<RandomGen*>(gen));
  const uint   rgenIndex    = tid % uint(m_randomGens.size()); ////////////////// change
  m_randomGens[rgenIndex]   = *gen;
  if(m_disableImageContrib !=0)
    return;
  
  /////////////////////////////////////////////////////////////////////////////// change
  uint x = uint(pixelOffsets.x*float(m_winWidth));
  uint y = uint(pixelOffsets.y*float(m_winHeight));
  
  if(x >= uint(m_winWidth-1))
    x = uint(m_winWidth-1);
  if(y >= uint(m_winHeight-1))
    y = uint(m_winHeight-1);
  /////////////////////////////////////////////////////////////////////////////// change   

  float4 specSamples = *a_accumColor; 
  float4 tmpVal      = specSamples*m_camRespoceRGB;
  float3 rgb         = to_float3(tmpVal);
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

  float4 colorRes  = m_exposureMult * to_float4(rgb, 1.0f);
  m_allColors[tid] = colorRes;

  // /////////////////////////////////////////////////////////////////////////////// change, add atomics
  // if(channels == 1)
  // {
  //   const float mono = 0.2126f*colorRes.x + 0.7152f*colorRes.y + 0.0722f*colorRes.z;
  //   #pragma omp atomic
  //   out_color[y*m_winWidth+x] += mono;
  // }
  // else if(channels <= 4)
  // { 
  //   #pragma omp atomic
  //   out_color[(y*m_winWidth+x)*channels + 0] += colorRes.x;
  //   #pragma omp atomic
  //   out_color[(y*m_winWidth+x)*channels + 1] += colorRes.y;
  //   #pragma omp atomic
  //   out_color[(y*m_winWidth+x)*channels + 2] += colorRes.z;
  // }
  // else
  // {
  //   auto waves = (*wavelengths);
  //   auto color = (*a_accumColor)*m_exposureMult;
  //   for(int i=0;i<4;i++) {
  //     const float t         = (waves[i] - LAMBDA_MIN)/(LAMBDA_MAX-LAMBDA_MIN);
  //     const int channelId   = std::min(int(float(channels)*t), int(channels)-1);
  //     const int offsetPixel = int(y)*m_winWidth + int(x);
  //     const int offsetLayer = channelId*m_winWidth*m_winHeight;
  //     #pragma omp atomic
  //     out_color[offsetLayer + offsetPixel] += color[i];
  //   }
  // }
  // /////////////////////////////////////////////////////////////////////////////// end change, atomics
}

float4 IntegratorKMLT::F(const std::vector<float>& xVal, uint tid, int*pX, int* pY)
{
  float* pRands = m_allRandomNumbers.data() + m_traceDepth*PER_BOUNCE*tid;
  assert(pRands + xVal.size() <= m_allRandomNumbers.data() + m_allRandomNumbers.size()); // check if we are out of bounds
  memcpy(pRands, xVal.data(), xVal.size()*sizeof(float));
  PathTrace(tid, 4, nullptr);
  (*pX) = m_allXY[tid].x;
  (*pY) = m_allXY[tid].y;
  return m_allColors[tid];
}

static inline float contribFunc(float4 color) // WORKS ONLY FOR RGB(!!!)
{
  return std::max(0.333334f*(color.x + color.y + color.z), 0.0f);
}

constexpr float MUTATE_COEFF_SCREEN = 128.0f;
constexpr float MUTATE_COEFF_BSDF   = 64.0f;

/**
\brief mutate random number in primary sample space -- interval [0,1]
\param valueX    - input original value
\param rands     - input pseudo random numbers
\param p2        - parameter of step size. The greater parameter is, the smaller step we gain. Default = 64.0f;
\param p1        - parameter of step size. The greater parameter is, the smaller step we gain. Default = 1024.0f;
\return mutated random float in range [0,1]

*/
static inline float MutateKelemen(float valueX, float2 rands, float p2, float p1)   // mutate in primary space
{
  const float s1    = 1.0f / p1;
  const float s2    = 1.0f / p2;
  const float power = -std::log(s2 / s1);
  const float dv    = std::max(s2*( std::exp(power*std::sqrt(rands.x)) - std::exp(power) ), 0.0f);

  if (rands.y < 0.5f)
  {
    valueX += dv;
    if (valueX > 1.0f)
      valueX -= 1.0f;
  }
  else
  {
    valueX -= dv;
    if (valueX < 0.0f)
      valueX += 1.0f;
  }

  return valueX;
}

std::vector<float> MutatePrimarySpace(const std::vector<float>& v2, RandomGen* pGen, int bounceNum)
{
  std::vector<float> res(v2.size());

  res[0] = MutateKelemen(v2[0], rndFloat2_Pseudo(pGen), MUTATE_COEFF_SCREEN*1.0f, 1024.0f); 
  res[1] = MutateKelemen(v2[1], rndFloat2_Pseudo(pGen), MUTATE_COEFF_SCREEN*1.0f, 1024.0f);
  res[2] = MutateKelemen(v2[2], rndFloat2_Pseudo(pGen), MUTATE_COEFF_BSDF, 1024.0f);
  res[3] = MutateKelemen(v2[3], rndFloat2_Pseudo(pGen), MUTATE_COEFF_BSDF, 1024.0f);
  res[4] = MutateKelemen(v2[4], rndFloat2_Pseudo(pGen), MUTATE_COEFF_BSDF, 1024.0f);

  for(int bounce=0; bounce < bounceNum; bounce++) 
  {
    const uint matIndex = IntegratorKMLT::BOUNCE_START + bounce*IntegratorKMLT::PER_BOUNCE + IntegratorKMLT::MATS_ID;
    const uint ltgIndex = IntegratorKMLT::BOUNCE_START + bounce*IntegratorKMLT::PER_BOUNCE + IntegratorKMLT::LGHT_ID;
    
    res[matIndex+0] = MutateKelemen(v2[matIndex+0], rndFloat2_Pseudo(pGen), MUTATE_COEFF_BSDF, 1024.0f);
    res[matIndex+1] = MutateKelemen(v2[matIndex+1], rndFloat2_Pseudo(pGen), MUTATE_COEFF_BSDF, 1024.0f);
    res[matIndex+2] = MutateKelemen(v2[matIndex+2], rndFloat2_Pseudo(pGen), MUTATE_COEFF_BSDF, 1024.0f);
    res[matIndex+3] = MutateKelemen(v2[matIndex+3], rndFloat2_Pseudo(pGen), MUTATE_COEFF_BSDF, 1024.0f);

    res[ltgIndex+0] = MutateKelemen(v2[ltgIndex+0], rndFloat2_Pseudo(pGen), MUTATE_COEFF_BSDF, 1024.0f);
    res[ltgIndex+1] = MutateKelemen(v2[ltgIndex+1], rndFloat2_Pseudo(pGen), MUTATE_COEFF_BSDF, 1024.0f);
    res[ltgIndex+2] = MutateKelemen(v2[ltgIndex+2], rndFloat2_Pseudo(pGen), MUTATE_COEFF_BSDF, 1024.0f);
    res[ltgIndex+3] = MutateKelemen(v2[ltgIndex+3], rndFloat2_Pseudo(pGen), MUTATE_COEFF_BSDF, 1024.0f);
  }

  return res;
}

void IntegratorKMLT::PathTraceBlock(uint pixelsNum, uint channels, float* out_color, uint a_passNum)
{
  this->SetFrameBufferLayer(FB_INDIRECT);
  
  uint maxThreads = 1;
  #ifndef _DEBUG
  #pragma omp parallel
  {
    #pragma omp single
    maxThreads = omp_get_max_threads();
  }
  #endif

  m_maxThreadId = maxThreads;
  m_allRandomNumbers.resize(m_traceDepth*RandsPerThread()*maxThreads);
  m_allColors.resize(maxThreads);
  m_allXY.resize(maxThreads);
  m_randomGens.resize(maxThreads);

  const size_t samplesPerPass = size_t(pixelsNum)*size_t(a_passNum) / maxThreads;

  ConsoleProgressBar progress(pixelsNum*a_passNum);
  progress.Start();
  auto start = std::chrono::high_resolution_clock::now();
  
  double avgBrightnessOut = 0.0f;

  #ifndef _DEBUG
  #pragma omp parallel
  #endif
  {
    int tid = omp_get_thread_num();
    
    const int NRandomisation = (clock() % 9) + (clock() % 4) + 1;
    RandomGen gen1 = RandomGenInit(tid*7 + 1);
    RandomGen gen2 = RandomGenInit(tid);
    for (int i = 0; i < NRandomisation; i++) 
    {
      NextState(&gen1);
      NextState(&gen2);
    }

    // (1) Initial State
    //
    int xScr = 0, yScr = 0;
    std::vector<float> xVec(m_traceDepth*RandsPerThread());
    for(size_t i=0;i<xVec.size();i++)
      xVec[i] = rndFloat1_Pseudo(&gen2);

    float4 yColor = F(xVec, tid, &xScr, &yScr);
    float  y      = contribFunc(yColor);
   
    // (2) Markov Chain
    //
    size_t accept     = 0;
    size_t largeSteps = 0;
    double accumBrightness = 0.0;

    for(size_t i=0;i<samplesPerPass;i++) 
    {
      auto xOld = xVec;

      const float plarge     = 0.33f;                           // 33% for large step;
      const bool isLargeStep = true; // (rndFloat1_Pseudo(&gen1) < plarge);
      
      std::vector<float> xNew(xOld.size());
      {
        if (isLargeStep)                                        // large step
        {
          for(size_t i=0;i<xNew.size();i++)
            xNew[i] = rndFloat1_Pseudo(&gen2);
        }
        else
          xNew = MutatePrimarySpace(xOld, &gen2, m_traceDepth); // small step
      }

      float  yOld      = y;
      float4 yOldColor = yColor;
  
      int xScrOld = xScr, yScrOld = yScr;
      int xScrNew = 0,    yScrNew = 0;
  
      float4 yNewColor = F(xNew, tid, &xScrNew, &yScrNew);
      float  yNew      = contribFunc(yNewColor);
      
      //{
      //  const int offset = yScrNew*m_winWidth + xScrNew;
      //  #pragma omp atomic
      //  out_color[offset*4+0] += yNewColor.x;
      //  #pragma omp atomic
      //  out_color[offset*4+1] += yNewColor.y;
      //  #pragma omp atomic
      //  out_color[offset*4+2] += yNewColor.z;
      //}
      
      float a = (yOld == 0.0f) ? 1.0f : std::min(1.0f, yNew / yOld);
      float p = rndFloat1_Pseudo(&gen1);

      if (p <= a) // accept 
      {
        xVec   = xNew;
        y      = yNew;
        yColor = yNewColor;
        xScr   = xScrNew;
        yScr   = yScrNew;
        accept++;
      }
      else        // reject
      {
        //x      = x;
        //y      = y;
        //yColor = yColor;
      }

      if(isLargeStep)
      {
        accumBrightness += double(yNew);
        largeSteps++;
      }

      // (5) contrib to image
      //
      float3 contribAtX = to_float3(yOldColor)*(1.0f / std::max(yOld, 1e-6f))*(1.0f - a);
      float3 contribAtY = to_float3(yNewColor)*(1.0f / std::max(yNew, 1e-6f))*a;

      if (dot(contribAtX, contribAtX) > 1e-12f)
      { 
        const int offset = yScrOld*m_winWidth + xScrOld;
        #pragma omp atomic
        out_color[offset*4+0] += contribAtX.x;
        #pragma omp atomic
        out_color[offset*4+1] += contribAtX.y;
        #pragma omp atomic
        out_color[offset*4+2] += contribAtX.z;
      }

      if (dot(contribAtY, contribAtY) > 1e-12f)
      { 
        const int offset = yScrNew*m_winWidth + xScrNew;
        #pragma omp atomic
        out_color[offset*4+0] += contribAtY.x;
        #pragma omp atomic
        out_color[offset*4+1] += contribAtY.y;
        #pragma omp atomic
        out_color[offset*4+2] += contribAtY.z;
      }
      

      progress.Update();
    }

    double avgBrightness = accumBrightness / double(largeSteps);
    #pragma omp atomic
    avgBrightnessOut += avgBrightness;
  }

  avgBrightnessOut = avgBrightnessOut / double(m_maxThreadId);

  double actualBrightness = 0.0;
  {
    for(uint i=0;i<pixelsNum;i++)
    {
      float4 color = float4(out_color[i*4+0], out_color[i*4+1], out_color[i*4+2], out_color[i*4+3]);
      actualBrightness += contribFunc(color);
    }
    actualBrightness /= double(pixelsNum);
  }
  
  const float normConst = float(a_passNum)*float(avgBrightnessOut/actualBrightness);
  for(uint i=0;i<pixelsNum*4;i++)
    out_color[i] *= normConst;

  progress.Done();
  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<Integrator> CreateIntegratorKMLT(int a_maxThreads = 1, int a_spectral_mode = 0, std::vector<uint32_t> a_features = {}) 
{ 
  return std::make_shared<IntegratorKMLT>(a_maxThreads, a_spectral_mode, a_features); 
}