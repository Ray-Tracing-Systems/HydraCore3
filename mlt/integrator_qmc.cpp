#include "integrator_qmc.h"
#include "utils.h" // for progress bar

#include <cassert>

IntegratorQMC:: IntegratorQMC(int a_maxThreads, int a_spectral_mode, std::vector<uint32_t> a_features) : Integrator(a_maxThreads, a_spectral_mode, a_features)
{
  qmc::init(m_qmcTable);
  
  // TODO: move this to LoadScene
  //

  // init optic parameters
  //
  lines = {
    { .curvatureRadius=0.0302249007f,  .thickness=0.00083350006f, .eta=1.62f,        .apertureRadius=0.0151700005f},   // 0
    { .curvatureRadius=0.0113931f,     .thickness=0.00741360011f, .eta=1.0f,         .apertureRadius=0.0103400005f},   // 1
    { .curvatureRadius=0.0752018988f,  .thickness=0.00106540008f, .eta=1.63900006f,  .apertureRadius=0.00889999978f}, // 2
    { .curvatureRadius=0.00833490025f,  .thickness=0.0111549003f,   .eta=1.0f,         .apertureRadius=0.00671000034f},
    { .curvatureRadius=0.00958819967f,  .thickness=0.00200540014f,  .eta=1.65400004f,  .apertureRadius=0.00451000035f},
    { .curvatureRadius=0.0438676998f,   .thickness=0.00538950041f,  .eta=1.0f,         .apertureRadius=0.00407000026f},
    { .curvatureRadius=0.0f,            .thickness=0.00141630007f,  .eta=0.0f,         .apertureRadius=0.00275000022f},
    { .curvatureRadius=0.0294541009f,   .thickness=0.00219339994f,  .eta=1.51699996f,  .apertureRadius=0.00298000011f},
    { .curvatureRadius=-0.00522650033f, .thickness=0.000971400063f, .eta=1.80499995f,  .apertureRadius=0.00292000012f},
    { .curvatureRadius=-0.0142884003f,  .thickness=6.27000045e-05f, .eta=1.0f,         .apertureRadius=0.00298000011f},
    { .curvatureRadius=-0.0223726016f,  .thickness=0.000940000056f, .eta=1.67299998f,  .apertureRadius=0.00298000011f},
    { .curvatureRadius=-0.0150404004f,  .thickness=0.0233591795f,   .eta=1.0f,         .apertureRadius=0.00326000014f},
  };

  std::reverse(lines.begin(), lines.end()); // if you need this ...  

  m_diagonal = 0.035f;
  m_aspect   = 1.0f;

  // CalcPhysSize();
  //
  m_physSize.x = 2.0f*std::sqrt(m_diagonal * m_diagonal / (1.0f + m_aspect * m_aspect));
  m_physSize.y = m_aspect * m_physSize.x;
}

void IntegratorQMC::EnableQMC()
{
  const bool dof    = (m_camLensRadius > 0.0f) || (m_enableOpticSim != 0);
  const bool spd    = (m_spectral_mode != 0);
  const bool motion = (m_normMatrices2.size() != 0);

  assert(qmc::QRNG_DIMENSIONS >= 11); // code of this function assume >= 11 dimentions

  if(dof && spd && motion) // (0,1): pixel id; (2,3): lens; (4,5): spd and motion; (6,7): mat; (8,9,10): light.
  {
    m_qmcMotionDim = 4;
    m_qmcSpdDim    = 5;
    m_qmcMatDim    = 0; // 6; // may set to zero to enable OMC here
    m_qmcLgtDim    = 0; // 8; // may set to zero to enable OMC here
  }
  else if(dof && spd)    // (0,1): pixel id; (2,3): lens; (4): spd; (5,6): mat; (7,8,9): light.
  {
    m_qmcMotionDim = 0;
    m_qmcSpdDim    = 4;
    m_qmcMatDim    = 5;
    m_qmcLgtDim    = 7;
  }
  else if(spd && motion) // (0,1): pixel id; (2,3): motion and spd; (4,5): mat; (6,7,8): light.
  {
    m_qmcMotionDim = 2;
    m_qmcSpdDim    = 3;
    m_qmcMatDim    = 4;
    m_qmcLgtDim    = 6;

  }
  else if(dof && motion) // (0,1): pixel id; (2,3): lens; (4): motion; (5,6): mat; (7,8,9): light.
  {
    m_qmcMotionDim = 4;
    m_qmcSpdDim    = 0;
    m_qmcMatDim    = 5; // may set to zero to enable OMC here
    m_qmcLgtDim    = 7; // may set to zero to enable OMC here
  }
  else if(dof)  // (0,1): pixel id; (2,3): lens; (4,5): mat; (6,7,8): light.
  {
    m_qmcMotionDim = 0;
    m_qmcSpdDim    = 0;
    m_qmcMatDim    = 4;
    m_qmcLgtDim    = 6;
  }
  else if(spd)  // (0,1): pixel id; (2,3): mat; (4) : spd; (5,6,7): light.
  {
    m_qmcMotionDim = 0;
    m_qmcSpdDim    = 4;
    m_qmcMatDim    = 2;
    m_qmcLgtDim    = 5;
  }
  else if(motion) // (0,1): pixel id; (2,3): mat; (4) : motion; (5,6,7): light.
  {
    m_qmcMotionDim = 4;
    m_qmcSpdDim    = 0;
    m_qmcMatDim    = 2;
    m_qmcLgtDim    = 5;
  }
  else // (0,1): pixel id; (2,3): mat; (4,5,6): light.
  {
    m_qmcMotionDim = 0;
    m_qmcSpdDim    = 0;
    m_qmcMatDim    = 2;
    m_qmcLgtDim    = 4;
  }

  m_qmcDofDim = 2;
}

float  IntegratorQMC::GetRandomNumbersSpec(uint tid, RandomGen* a_gen) 
{ 
  if(m_qmcSpdDim != 0)
    return qmc::rndFloat(tid, m_qmcSpdDim, m_qmcTable[0]);
  else 
    return rndFloat1_Pseudo(a_gen);
}

float  IntegratorQMC::GetRandomNumbersTime(uint tid, RandomGen* a_gen) 
{ 
  if(m_qmcMotionDim != 0)
    return qmc::rndFloat(tid, m_qmcMotionDim, m_qmcTable[0]);
  else
    return rndFloat1_Pseudo(a_gen);
}

float4 IntegratorQMC::GetRandomNumbersLens(uint tid, RandomGen* a_gen) 
{ 
  float4 rands = rndFloat4_Pseudo(a_gen);
  rands.x = qmc::rndFloat(tid, 0, m_qmcTable[0]);
  rands.y = qmc::rndFloat(tid, 1, m_qmcTable[0]);
  if((m_camLensRadius > 0.0f || (m_enableOpticSim != 0)) && m_qmcDofDim != 0)
  {
    rands.z = qmc::rndFloat(tid, 2, m_qmcTable[0]);
    rands.w = qmc::rndFloat(tid, 3, m_qmcTable[0]);
  }
  return rands; 
}

float4 IntegratorQMC::GetRandomNumbersMats(uint tid, RandomGen* a_gen, int a_bounce) 
{ 
  float4 rands = rndFloat4_Pseudo(a_gen);
  if(a_bounce == 0 && m_qmcMatDim != 0)
  {
    rands.x = qmc::rndFloat(tid, m_qmcMatDim + 0, m_qmcTable[0]);
    rands.y = qmc::rndFloat(tid, m_qmcMatDim + 1, m_qmcTable[0]);
  }
  return rands; 
}

float4 IntegratorQMC::GetRandomNumbersLgts(uint tid, RandomGen* a_gen, int a_bounce)
{
  float4 rands = rndFloat4_Pseudo(a_gen); // don't use single rndFloat4 (!!!)
  float  rndId = rndFloat1_Pseudo(a_gen); // don't use single rndFloat4 (!!!)
  if(a_bounce == 0 && m_qmcLgtDim != 0)
  {
    rands.x = qmc::rndFloat(tid, m_qmcLgtDim + 0, m_qmcTable[0]);
    rands.y = qmc::rndFloat(tid, m_qmcLgtDim + 1, m_qmcTable[0]);
    rndId   = qmc::rndFloat(tid, m_qmcLgtDim + 2, m_qmcTable[0]);
  }
  return float4(rands.x, rands.y, rands.z, rndId);
}

float IntegratorQMC::GetRandomNumbersMatB(uint tid, RandomGen* a_gen, int a_bounce, int a_layer) 
{ 
  return rndFloat1_Pseudo(a_gen); 
}

uint IntegratorQMC::RandomGenId(uint tid) { return tid % uint(m_randomGens.size()); }

Integrator::EyeRayData IntegratorQMC::SampleCameraRay(RandomGen* pGen, uint tid)
{
  const float4 pixelOffsets = GetRandomNumbersLens(tid, pGen);
  float3 rayDir = EyeRayDirNormalized(pixelOffsets.x, pixelOffsets.y, m_projInv);
  float3 rayPos = float3(0,0,0);

  if (m_camLensRadius > 0.0f)
  {
    const float tFocus         = m_camTargetDist / (-rayDir.z);
    const float3 focusPosition = rayPos + rayDir*tFocus;
    const float2 xy            = m_camLensRadius*2.0f*MapSamplesToDisc(float2(pixelOffsets.z - 0.5f, pixelOffsets.w - 0.5f));
    rayPos.x += xy.x;
    rayPos.y += xy.y;
    rayDir = normalize(focusPosition - rayPos);
  }
  else if(m_enableOpticSim != 0) // not nessesary part of QMC. Just implemented here for test cases, could be moved in main class further  
  {
    const float2 xy = 0.25f*m_physSize*float2(2.0f*pixelOffsets.x - 1.0f, 2.0f*pixelOffsets.y - 1.0f);
    
    rayPos = float3(xy.x, xy.y, 0);
    
    const float2 rareSam  = LensRearRadius()*2.0f*MapSamplesToDisc(float2(pixelOffsets.z - 0.5f, pixelOffsets.w - 0.5f));
    const float3 shootTo  = float3(rareSam.x, rareSam.y, LensRearZ());
    const float3 ray_dirF = normalize(shootTo - rayPos);
    
    float cosTheta = std::abs(ray_dirF.z);
    rayDir         = ray_dirF;
    bool success   = TraceLensesFromFilm(rayPos, rayDir);
    
    if (!success) 
    {
      rayPos = float3(0,-10000000.0,0.0); // shoot ray under the floor
      rayDir = float3(0,-1,0);
    }
    else
    {
      rayDir = float3(-1,-1,-1)*normalize(rayDir);
      rayPos = float3(-1,-1,-1)*rayPos;
    }
  }
  
  /////////////////////////////////////////////////////
  uint x = uint(pixelOffsets.x*float(m_winWidth));
  uint y = uint(pixelOffsets.y*float(m_winHeight));
  if(x >= uint(m_winWidth-1))
    x = uint(m_winWidth-1);
  if(y >= uint(m_winHeight-1))
    y = uint(m_winHeight-1);
  /////////////////////////////////////////////////////

  EyeRayData res;
  {
    res.rayPos = rayPos;
    res.rayDir = rayDir;
    res.x      = x;
    res.y      = y;
    res.timeSam = 0.0f;
    res.waveSam = 1.0f;
    if(m_normMatrices2.size() != 0)
      res.timeSam = GetRandomNumbersTime(tid, pGen);
    if(KSPEC_SPECTRAL_RENDERING !=0 && m_spectral_mode != 0)
      res.waveSam = GetRandomNumbersSpec(tid, pGen);
    res.cosTheta = 1.0f;
  }
  
  RecordPixelRndIfNeeded(float2(pixelOffsets.x, pixelOffsets.y), res.waveSam);

  return res;
}

void IntegratorQMC::kernel_ContributeToImage(uint tid, const uint* rayFlags, uint channels, const float4* a_accumColor, const RandomGen* gen,
                                             const uint* in_pakedXY, const float4* wavelengths, float* out_color)
{
  
  if(tid >= m_maxThreadId) // don't contrubute to image in any "record" mode
    return;
  
  const float4 pixelOffsets = GetRandomNumbersLens(tid, const_cast<RandomGen*>(gen));
  m_randomGens[RandomGenId(tid)] = *gen;
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
  using IndexType = unsigned; // or int
  const size_t samplesNum = std::min<size_t>(size_t(std::numeric_limits<IndexType>::max()), size_t(pixelsNum)*size_t(a_passNum));
  m_maxThreadId = uint(samplesNum);

  EnableQMC();

  std::cout << "[IntegratorQMC]: max spp = " << samplesNum/pixelsNum << std::endl;
  std::cout.flush();

  ConsoleProgressBar progress(pixelsNum*a_passNum);
  progress.Start();
  auto start = std::chrono::high_resolution_clock::now();

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline bool Quadratic(float A, float B, float C, float *t0, float *t1) {
  // Find quadratic discriminant
  double discrim = (double)B * (double)B - 4. * (double)A * (double)C;
  if (discrim < 0.) 
    return false;
  double rootDiscrim = std::sqrt(discrim);
  float floatRootDiscrim   = float(rootDiscrim);
  // Compute quadratic _t_ values
  float q;
  if ((float)B < 0)
      q = -.5 * (B - floatRootDiscrim);
  else
      q = -.5 * (B + floatRootDiscrim);
  *t0 = q / A;
  *t1 = C / q;
  if ((float)*t0 > (float)*t1) 
  {
    // std::swap(*t0, *t1);
    float temp = *t0;
    *t0 = *t1;
    *t1 = temp;
  }
  return true;
}

static inline bool Refract(const float3 wi, const float3 n, float eta, float3 *wt) {
  // Compute $\cos \theta_\roman{t}$ using Snell's law
  float cosThetaI  = dot(n, wi);
  float sin2ThetaI = std::max(float(0), float(1.0f - cosThetaI * cosThetaI));
  float sin2ThetaT = eta * eta * sin2ThetaI;
  // Handle total internal reflection for transmission
  if (sin2ThetaT >= 1) return false;
  float cosThetaT = std::sqrt(1 - sin2ThetaT);
  *wt = eta * (-1.0f)*wi + (eta * cosThetaI - cosThetaT) * n;
  return true;
}

static inline float3 faceforward(const float3 n, const float3 v) { return (dot(n, v) < 0.f) ? (-1.0f)*n : n; }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool IntegratorQMC::IntersectSphericalElement(float radius, float zCenter, float3 rayPos, float3 rayDir, 
                                              float *t, float3 *n) const
{
  // Compute _t0_ and _t1_ for ray--element intersection
  const float3 o = rayPos - float3(0, 0, zCenter);
  const float  A = rayDir.x * rayDir.x + rayDir.y * rayDir.y + rayDir.z * rayDir.z;
  const float  B = 2 * (rayDir.x * o.x + rayDir.y * o.y + rayDir.z * o.z);
  const float  C = o.x * o.x + o.y * o.y + o.z * o.z - radius * radius;
  float  t0, t1;
  if (!Quadratic(A, B, C, &t0, &t1)) 
    return false;
  
  // Select intersection $t$ based on ray direction and element curvature
  bool useCloserT = (rayDir.z > 0.0f) != (radius < 0.0);
  *t = useCloserT ? std::min(t0, t1) : std::max(t0, t1);
  if (*t < 0.0f) 
    return false;
  
  // Compute surface normal of element at ray intersection point
  *n = normalize(o + (*t)*rayDir);
  *n = faceforward(*n, -1.0f*rayDir);
  return true;
}

bool IntegratorQMC::TraceLensesFromFilm(float3& inoutRayPos, float3& inoutRayDir) const
{
  float elementZ = 0;
  // Transform _rCamera_ from camera to lens system space
  // 
  float3 rayPosLens = float3(inoutRayPos.x, inoutRayPos.y, -inoutRayPos.z);
  float3 rayDirLens = float3(inoutRayDir.x, inoutRayDir.y, -inoutRayDir.z);

  for(int i=0; i<lines.size(); i++)
  {
    const auto element = lines[i];                                  
    // Update ray from film accounting for interaction with _element_
    elementZ -= element.thickness;
    
    // Compute intersection of ray with lens element
    float t;
    float3 n;
    bool isStop = (element.curvatureRadius == 0.0f);
    if (isStop) 
    {
      // The refracted ray computed in the previous lens element
      // interface may be pointed towards film plane(+z) in some
      // extreme situations; in such cases, 't' becomes negative.
      if (rayDirLens.z >= 0.0f) 
        return false;
      t = (elementZ - rayPosLens.z) / rayDirLens.z;
    } 
    else 
    {
      const float radius  = element.curvatureRadius;
      const float zCenter = elementZ + element.curvatureRadius;
      if (!IntersectSphericalElement(radius, zCenter, rayPosLens, rayDirLens, &t, &n))
        return false;
    }

    // Test intersection point against element aperture
    const float3 pHit = rayPosLens + t*rayDirLens;
    const float r2    = pHit.x * pHit.x + pHit.y * pHit.y;
    if (r2 > element.apertureRadius * element.apertureRadius) 
      return false;
    
    rayPosLens = pHit;
    // Update ray path for from-scene element interface interaction
    if (!isStop) 
    {
      float3 wt;
      float etaI = lines[i+0].eta;                                                      
      float etaT = (i == lines.size()-1) ? 1.0f : lines[i+1].eta;
      if(etaT == 0.0f)
        etaT = 1.0f;                                                          
      if (!Refract(normalize((-1.0f)*rayDirLens), n, etaI / etaT, &wt))
        return false;
      rayDirLens = wt;
    }

  }

  // Transform _rLens_ from lens system space back to camera space
  //
  inoutRayPos = float3(rayPosLens.x, rayPosLens.y, -rayPosLens.z);
  inoutRayDir = float3(rayDirLens.x, rayDirLens.y, -rayDirLens.z);
  return true;  
}
