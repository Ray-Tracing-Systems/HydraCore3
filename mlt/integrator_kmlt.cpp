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
                            RandomGen* gen, uint* rayFlags, MisData* misData, int* pX, int* pY);


  void   PathTraceBlock(uint pixelsNum, uint channels, float* out_color, uint a_passNum) override;
  
  float4 PathTraceF(uint tid, int*pX, int* pY);

  static constexpr uint LGHT_ID      = 0;
  static constexpr uint MATS_ID      = 4;
  static constexpr uint BLND_ID      = 8;
  static constexpr uint PER_BOUNCE   = 10;
  static constexpr uint BOUNCE_START = 5;
  
  inline uint RandsPerThread() const { return PER_BOUNCE*m_traceDepth + BOUNCE_START; }

  std::vector<float>    m_allRands;
  std::vector<int>      m_allLargeSteps;
  uint                  m_randsPerThread = 0;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

#define LAZY_EVALUATION 0

float4 IntegratorKMLT::GetRandomNumbersLens(uint tid, RandomGen* a_gen) 
{ 
  float* data  = m_allRands.data() + tid*m_randsPerThread;
#if LAZY_EVALUATION  
  float4 rands = rndFloat4_Pseudo(a_gen);
  
  if(m_allLargeSteps[tid] == 0)
  {
    const float2 xy = float2(rands.x, rands.y);
    const float2 zw = float2(rands.z, rands.w);
    rands.x = MutateKelemen(data[0], xy, MUTATE_COEFF_SCREEN*1.0f, 1024.0f); // screen
    rands.y = MutateKelemen(data[1], zw, MUTATE_COEFF_SCREEN*1.0f, 1024.0f); // screen
    //rands.z = MutateKelemen(data[2], rndFloat2_Pseudo(a_gen), MUTATE_COEFF_BSDF, 1024.0f);        // lens
    //rands.w = MutateKelemen(data[3], rndFloat2_Pseudo(a_gen), MUTATE_COEFF_BSDF, 1024.0f);        // lens
  }

  data[0] = rands.x;
  data[1] = rands.y;
  //data[2] = rands.z;
  //data[3] = rands.w;

  return rands;
#else
  return float4(data[0], data[1], data[2], data[3]);
#endif
}

float  IntegratorKMLT::GetRandomNumbersSpec(uint tid, RandomGen* a_gen) 
{ 
  float* data   = m_allRands.data() + tid*m_randsPerThread;
#if LAZY_EVALUATION  
  float2 rndVal = rndFloat2_Pseudo(a_gen);
  
  if(m_allLargeSteps[tid] == 0)
  {
    float2 copy = rndVal;
    rndVal.x = MutateKelemen(data[4], copy, MUTATE_COEFF_BSDF, 1024.0f);
  }

  data[4] = rndVal.x;
  return rndVal.x; 
#else
  return data[4]; 
#endif
}

float4 IntegratorKMLT::GetRandomNumbersMats(uint tid, RandomGen* a_gen, int a_bounce) 
{ 
  float* data  = m_allRands.data() + tid*m_randsPerThread + BOUNCE_START + a_bounce*PER_BOUNCE + MATS_ID;
#if LAZY_EVALUATION  
  float4 rands = rndFloat4_Pseudo(a_gen);
  
  if(m_allLargeSteps[tid] == 0)
  {
    const float2 xy = float2(rands.x, rands.y);
    const float2 zw = float2(rands.z, rands.w);
    const float4 r2 = rndFloat4_Pseudo(a_gen);
    rands.x = MutateKelemen(data[0], xy,                 MUTATE_COEFF_BSDF, 1024.0f); // 
    rands.y = MutateKelemen(data[1], zw,                 MUTATE_COEFF_BSDF, 1024.0f); // 
    rands.z = MutateKelemen(data[2], float2(r2.x, r2.y), MUTATE_COEFF_BSDF, 1024.0f); // 
    rands.w = MutateKelemen(data[3], float2(r2.z, r2.w), MUTATE_COEFF_BSDF, 1024.0f); // 
  }
  
  data[0] = rands.x;
  data[1] = rands.y;
  data[2] = rands.z;
  data[3] = rands.w;

  return rands;
#else
  return float4(data[0], data[1], data[2], data[3]);
#endif
}

float4 IntegratorKMLT::GetRandomNumbersLgts(uint tid, RandomGen* a_gen, int a_bounce)
{
  float* data  = m_allRands.data() + tid*m_randsPerThread + BOUNCE_START + a_bounce*PER_BOUNCE + LGHT_ID;
#if LAZY_EVALUATION  
  float4 rands = rndFloat4_Pseudo(a_gen);
  
  if(m_allLargeSteps[tid] == 0)
  {
    const float2 xy = float2(rands.x, rands.y);
    const float2 zw = float2(rands.z, rands.w);
    const float4 r2 = rndFloat4_Pseudo(a_gen);
    rands.x = MutateKelemen(data[0], xy,                 MUTATE_COEFF_BSDF, 1024.0f); // 
    rands.y = MutateKelemen(data[1], zw,                 MUTATE_COEFF_BSDF, 1024.0f); // 
    rands.z = MutateKelemen(data[2], float2(r2.x, r2.y), MUTATE_COEFF_BSDF, 1024.0f); // 
    rands.w = MutateKelemen(data[3], float2(r2.z, r2.w), MUTATE_COEFF_BSDF, 1024.0f); // 
  }

  data[0] = rands.x;
  data[1] = rands.y;
  data[2] = rands.z;
  data[3] = rands.w;

  return rands;
#else
  return float4(data[0], data[1], data[2], data[3]);
#endif  
}

float IntegratorKMLT::GetRandomNumbersMatB(uint tid, RandomGen* a_gen, int a_bounce, int a_layer) 
{ 
  float* data   = m_allRands.data() + tid*m_randsPerThread + BOUNCE_START + a_bounce*PER_BOUNCE + BLND_ID;
#if LAZY_EVALUATION  
  float2 rndVal = rndFloat2_Pseudo(a_gen);
  
  if(m_allLargeSteps[tid] == 0)
  {
    float2 copy = rndVal;
    rndVal.x = MutateKelemen(data[a_layer], copy, MUTATE_COEFF_BSDF, 1024.0f);
  }
 
  data[a_layer] = rndVal.x;
  return rndVal.x;
#else
  return data[a_layer];
#endif   
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void IntegratorKMLT::kernel_InitEyeRay2(uint tid, const uint* packedXY, 
                                       float4* rayPosAndNear, float4* rayDirAndFar, float4* wavelengths, 
                                       float4* accumColor,    float4* accumuThoroughput,
                                       RandomGen* gen, uint* rayFlags, MisData* misData, int* pX, int* pY) // 
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

  (*pX) = int(x);
  (*pY) = int(y);

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

float4 IntegratorKMLT::PathTraceF(uint tid, int*pX, int* pY)
{
  float4 accumColor, accumThroughput;
  float4 rayPosAndNear, rayDirAndFar;
  float4 wavelengths;
  RandomGen gen; 
  MisData   mis;
  uint      rayFlags;
  kernel_InitEyeRay2(tid, m_packedXY.data(), &rayPosAndNear, &rayDirAndFar, &wavelengths, &accumColor, &accumThroughput, &gen, &rayFlags, &mis, pX, pY);

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

  kernel_HitEnvironment(tid, &rayFlags, &rayDirAndFar, &mis, &accumThroughput, &accumColor);

  return accumColor*m_exposureMult;
}

static inline float contribFunc(float4 color) // WORKS ONLY FOR RGB(!!!)
{
  return std::max(0.333334f*(color.x + color.y + color.z), 0.0f);
}

uint32_t AlignedSize(uint32_t a_size, uint32_t a_alignment)
{
  if (a_size % a_alignment == 0)
    return a_size;
  else
  {
    uint32_t sizeCut = a_size - (a_size % a_alignment);
    return sizeCut + a_alignment;
  }
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

  m_randsPerThread = AlignedSize(RandsPerThread()*m_traceDepth*maxThreads, uint32_t(16)); 

  m_maxThreadId = maxThreads;
  m_randomGens.resize(maxThreads);
  m_allRands.resize(maxThreads*m_randsPerThread);
  m_allLargeSteps.resize(maxThreads);

  const size_t samplesPerPass = (size_t(pixelsNum)*size_t(a_passNum)) / size_t(maxThreads);

  ConsoleProgressBar progress(pixelsNum*a_passNum);
  progress.Start();
  auto start = std::chrono::high_resolution_clock::now();
  
  double avgBrightnessOut = 0.0f;
  float  avgAcceptanceRate = 0.0f;

  std::cout << "[IntegratorKMLT]: state size = " << m_traceDepth*RandsPerThread() << std::endl;

  #ifndef _DEBUG
  #pragma omp parallel 
  //#pragma omp parallel for schedule(static)
  //for(uint tid=0;tid<maxThreads;tid++)
  #endif
  {
    int tid = omp_get_thread_num();
    
    RandomGen gen1 = RandomGenInit(tid*7 + 1);
    RandomGen gen2 = RandomGenInit(tid);
    for (int i = 0; i < 10 + tid % 17; i++) 
    {
      NextState(&gen1);
      NextState(&gen2);
    }

    // (1) Initial State
    //
    int xScr = 0, yScr = 0;
    std::vector<float> xVec(m_traceDepth*RandsPerThread());
    float* xNew = m_allRands.data() + m_randsPerThread*tid; //(m_traceDepth*RandsPerThread());
    
  #if LAZY_EVALUATION==0  
    for(size_t i=0;i<xVec.size();i++)
      xVec[i] = rndFloat1_Pseudo(&gen2);
    for(size_t i=0;i<xVec.size();i++) // eval F(xVec)
      xNew[i] = xVec[i];
  #else
    m_allLargeSteps[tid] = 1; // [lazy evaluation]: evaluate large step inside PathTraceF 
  #endif

    float4 yColor = PathTraceF(tid, &xScr, &yScr);
    float  y      = contribFunc(yColor);
    
  #if LAZY_EVALUATION
    for(size_t i=0;i<xVec.size();i++) // [lazy evaluation]: save initial state
      xVec[i] = xNew[i];
  #endif

    // (2) Markov Chain
    //
    size_t accept     = 0;
    size_t largeSteps = 0;
    double accumBrightness = 0.0;

    for(size_t i=0;i<samplesPerPass;i++) 
    {
      auto xOld = xVec;

      const float plarge     = 0.25f;                         // 25% of large step;
      const bool isLargeStep = (rndFloat1_Pseudo(&gen1) < plarge);
      
    #if LAZY_EVALUATION
      m_allLargeSteps[tid] = isLargeStep ? 1 : 0;             // pass step type to intergator for lazy evaluation
    #else
      if (isLargeStep)                                      // large step
      {
        for(size_t i=0;i<xVec.size();i+=4) {
          const float4 r1 = rndFloat4_Pseudo(&gen2);
          xNew[i+0] = r1.x;
          xNew[i+1] = r1.y;
          xNew[i+2] = r1.z;
          xNew[i+3] = r1.w;
        }
      }
      else
      {
        const float4 r1 = rndFloat4_Pseudo(&gen2);
        const float4 r2 = rndFloat4_Pseudo(&gen2);

        xNew[0] = MutateKelemen(xOld[0], float2(r1.x, r1.y),      MUTATE_COEFF_SCREEN*1.0f, 1024.0f); // screen 
        xNew[1] = MutateKelemen(xOld[1], float2(r1.z, r1.w),      MUTATE_COEFF_SCREEN*1.0f, 1024.0f); // screen
        xNew[2] = MutateKelemen(xOld[2], rndFloat2_Pseudo(&gen2), MUTATE_COEFF_BSDF, 1024.0f);        // lens
        xNew[3] = MutateKelemen(xOld[3], rndFloat2_Pseudo(&gen2), MUTATE_COEFF_BSDF, 1024.0f);        // lens
       
        for(size_t i = 4; i < xVec.size(); i+=2) { 
          const float4 r1 = rndFloat4_Pseudo(&gen2);
          xNew[i+0] = MutateKelemen(xOld[i+0], float2(r1.x, r1.y), MUTATE_COEFF_BSDF, 1024.0f);
          xNew[i+1] = MutateKelemen(xOld[i+1], float2(r1.z, r1.w), MUTATE_COEFF_BSDF, 1024.0f);
        }
      }
    #endif

      float  yOld      = y;
      float4 yOldColor = yColor;
  
      int xScrOld = xScr, yScrOld = yScr;
      int xScrNew = 0,    yScrNew = 0;
  
      float4 yNewColor = PathTraceF(tid, &xScrNew, &yScrNew); // eval F(xNew)
      float  yNew      = contribFunc(yNewColor);
      
      float a = (yOld == 0.0f) ? 1.0f : std::min(1.0f, yNew / yOld);
      float p = rndFloat1_Pseudo(&gen1);

      if (p <= a) // accept 
      { 
      #if LAZY_EVALUATION==0
        for(size_t i=0;i<xVec.size();i++)
          xVec[i] = xNew[i];
      #endif
        y      = yNew;
        yColor = yNewColor;
        xScr   = xScrNew;
        yScr   = yScrNew;
        accept++;

      #if LAZY_EVALUATION
        // [lazy evaluation]: save new state
        for(size_t i=0;i<xVec.size();i++) 
          xVec[i] = xNew[i];
      #endif
      }
      else        // reject
      {
        //x      = x;
        //y      = y;
        //yColor = yColor;

      #if LAZY_EVALUATION
        // [lazy evaluation]: restore previous state
        for(size_t i=0;i<xVec.size();i++) // save initial state
          xNew[i] = xVec[i];
      #endif
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

    #pragma omp atomic
    avgAcceptanceRate += float(accept);

    double avgBrightness = accumBrightness / double(largeSteps);
    #pragma omp atomic
    avgBrightnessOut += avgBrightness;
  }

  progress.Done();

  avgBrightnessOut  = avgBrightnessOut  / double(m_maxThreadId);
  avgAcceptanceRate = avgAcceptanceRate / float(pixelsNum*a_passNum);

  std::cout << "[IntegratorKMLT]: acceptance rate = " << avgAcceptanceRate << std::endl;
  std::cout.flush();

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
  for(uint i=0;i<pixelsNum*channels;i++)
    out_color[i] *= normConst;

  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<Integrator> CreateIntegratorKMLT(int a_maxThreads = 1, int a_spectral_mode = 0, std::vector<uint32_t> a_features = {}) 
{ 
  return std::make_shared<IntegratorKMLT>(a_maxThreads, a_spectral_mode, a_features); 
}