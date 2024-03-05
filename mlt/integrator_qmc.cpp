#include "integrator_qmc.h"
#include "utils.h" // for progress bar

float  IntegratorQMC::GetRandomNumbersSpec(uint tid, RandomGen* a_gen) 
{ 
  return rndFloat1_Pseudo(a_gen); 
}

float4 IntegratorQMC::GetRandomNumbersLens(uint tid, RandomGen* a_gen) 
{ 
  float4 rands = rndFloat4_Pseudo(a_gen);
  rands.x = qmc::rndFloat(tid, 0, m_qmcTable[0]);
  rands.y = qmc::rndFloat(tid, 1, m_qmcTable[0]);
  return rands; 
}

float4 IntegratorQMC::GetRandomNumbersMats(uint tid, RandomGen* a_gen, int a_bounce) 
{ 
  float4 rands = rndFloat4_Pseudo(a_gen);
  rands.x = qmc::rndFloat(tid, 2, m_qmcTable[0]);
  rands.y = qmc::rndFloat(tid, 3, m_qmcTable[0]);
  return rands; 
}

float4 IntegratorQMC::GetRandomNumbersLgts(uint tid, RandomGen* a_gen, int a_bounce)
{
  float4 rands = rndFloat4_Pseudo(a_gen); // don't use single rndFloat4 (!!!)
  float  rndId = rndFloat1_Pseudo(a_gen); // don't use single rndFloat4 (!!!)
  rndId   = qmc::rndFloat(tid, 4, m_qmcTable[0]);
  rands.x = qmc::rndFloat(tid, 5, m_qmcTable[0]);
  rands.y = qmc::rndFloat(tid, 6, m_qmcTable[0]);
  return float4(rands.x, rands.y, rands.z, rndId);
}

float IntegratorQMC::GetRandomNumbersMatB(uint tid, RandomGen* a_gen, int a_bounce, int a_layer) 
{ 
  return rndFloat1_Pseudo(a_gen); 
}

void IntegratorQMC::kernel_InitEyeRay2(uint tid, const uint* packedXY, 
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

void IntegratorQMC::kernel_ContributeToImage(uint tid, const uint* rayFlags, uint channels, const float4* a_accumColor, const RandomGen* gen,
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

  float4 colorRes = m_exposureMult * to_float4(rgb, 1.0f);
  
  /////////////////////////////////////////////////////////////////////////////// change, add atomics
  if(channels == 1)
  {
    const float mono = 0.2126f*colorRes.x + 0.7152f*colorRes.y + 0.0722f*colorRes.z;
    #pragma omp atomic
    out_color[y*m_winWidth+x] += mono;
  }
  else if(channels <= 4)
  { 
    #pragma omp atomic
    out_color[(y*m_winWidth+x)*channels + 0] += colorRes.x;
    #pragma omp atomic
    out_color[(y*m_winWidth+x)*channels + 1] += colorRes.y;
    #pragma omp atomic
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
      #pragma omp atomic
      out_color[offsetLayer + offsetPixel] += color[i];
    }
  }
   /////////////////////////////////////////////////////////////////////////////// end change, atomics
}

void IntegratorQMC::PathTraceBlock(uint pixelsNum, uint channels, float* out_color, uint a_passNum)
{
  ConsoleProgressBar progress(pixelsNum*a_passNum);
  progress.Start();
  auto start = std::chrono::high_resolution_clock::now();
  
  using IndexType = unsigned; // or int
  const size_t samplesNum = std::min<size_t>(size_t(std::numeric_limits<IndexType>::max()), size_t(pixelsNum)*size_t(a_passNum));
  m_maxThreadId = uint(samplesNum);

  #ifndef _DEBUG
  #pragma omp parallel for default(shared) 
  #endif
  for (IndexType i = 0; i < IndexType(samplesNum); ++i) 
  {
    PathTrace(i, channels, out_color);
    progress.Update();
  }
  progress.Done();
  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<Integrator> CreateIntegratorQMC(int a_maxThreads = 1, int a_spectral_mode = 0, std::vector<uint32_t> a_features = {}) 
{ 
  return std::make_shared<IntegratorQMC>(a_maxThreads, a_spectral_mode, a_features); 
}