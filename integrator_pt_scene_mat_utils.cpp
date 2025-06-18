#include "integrator_pt_scene.h"
#include "spectral/spec/conversions.h"
#include "spectral/spec/spectral_util.h"
#include "include/transfer_matrix.h"
#include "include/airy_reflectance.h"
#include "fourier/fourier.h"


static inline void save_to_file(const char* name, float *arr, int x_samples, int y_samples)
{
  std::ofstream precomp_file;
  precomp_file.open(name);
  for (int i = 0; i < x_samples; ++i)
  {
    for (int j = 0; j < y_samples; ++j)
    {
      precomp_file << arr[i * y_samples + j] << " ";
    }
  }
  precomp_file.close();
}

#include <chrono>
using namespace std::chrono;

ThinFilmPrecomputed precomputeThinFilmSpectral(
        const float extIOR, const uint* eta_id_vec, const uint* k_id_vec, const std::vector<float> &spec_values, 
        const std::vector<uint2> &spec_offsets, const float* eta_vec, const float* k_vec,
        const float* a_thickness, int layers)
{
  ThinFilmPrecomputed res;
  res.ext_reflectivity.resize(FILM_ANGLE_RES * FILM_LENGTH_RES);
  res.ext_transmittivity.resize(FILM_ANGLE_RES * FILM_LENGTH_RES);
  res.int_reflectivity.resize(FILM_ANGLE_RES * FILM_LENGTH_RES);
  res.int_transmittivity.resize(FILM_ANGLE_RES * FILM_LENGTH_RES);

  auto start = high_resolution_clock::now();
  for (size_t i = 0; i < FILM_LENGTH_RES; ++i)
  {
    float wavelength = (LAMBDA_MAX - LAMBDA_MIN - 1) / (FILM_LENGTH_RES - 1) * i + LAMBDA_MIN;
    std::vector<complex> ior;
    ior.reserve(layers + 1);
    ior[0] = complex(extIOR, 0.f);
    float eta, k;
    uint2 data;
    uint offset;
    uint size;
    for (size_t layer = 0; layer < layers; ++layer)
    {
      eta = eta_vec[layer];
      uint eta_id = eta_id_vec[layer];
      if (eta_id < 0xFFFFFFFF)
      {
        data  = spec_offsets[eta_id];
        offset = data.x;
        size   = data.y;
        eta = SampleUniformSpectrum(spec_values.data() + offset, {wavelength, 0, 0, 0}, size)[0];
      }

      k = k_vec[layer];
      uint k_id = k_id_vec[layer];
      if (k_id < 0xFFFFFFFF)
      {
        data  = spec_offsets[k_id];
        offset = data.x;
        size   = data.y;
        k = SampleUniformSpectrum(spec_values.data() + offset, {wavelength, 0, 0, 0}, size)[0];
      }

      ior[layer + 1] = complex(eta, k);
    }
    for (int j = 0; j < FILM_ANGLE_RES; ++j)
    {
      float theta = M_PI_2 / float(FILM_ANGLE_RES - 1) * j;
      
      float cosTheta = clamp(cosf(theta), 1e-3f, 1.f);
      FrReflRefr forward;
      FrReflRefr backward;
      if (layers == 2)
      {
        forward = FrFilm(cosTheta, ior[0], ior[1], ior[2], a_thickness[0], wavelength);
        backward = FrFilm(cosTheta, ior[2], ior[1], ior[0], a_thickness[0], wavelength);
        //forward = TransferMatrixForward(cosTheta, ior.data(), a_thickness, layers, wavelength);
        //backward = TransferMatrixBackward(cosTheta, ior.data(), a_thickness, layers, wavelength);
      }
      else
      {
        forward = multFrFilm(cosTheta, ior.data(), a_thickness, layers, wavelength);
        backward = multFrFilm_r(cosTheta, ior.data(), a_thickness, layers, wavelength);
        //forward = TransferMatrixForward(cosTheta, ior.data(), a_thickness, layers, wavelength);
        //backward = TransferMatrixBackward(cosTheta, ior.data(), a_thickness, layers, wavelength);
      }
      res.ext_reflectivity[i * FILM_ANGLE_RES + j] = forward.refl;
      res.ext_transmittivity[i * FILM_ANGLE_RES + j] = forward.refr;
      res.int_reflectivity[i * FILM_ANGLE_RES + j] = backward.refl;
      res.int_transmittivity[i * FILM_ANGLE_RES + j] = backward.refr;  
      #ifdef _DEBUG
      if(forward.refl < 0 || std::isnan(forward.refl) || std::isinf(forward.refl))
      {
        std::cout << "WARNING! Precomputed film external reflectance is " << forward.refl << std::endl;
      }
      if(forward.refr < 0 || std::isnan(forward.refr) || std::isinf(forward.refr))
      {
        std::cout << "WARNING! Precomputed film external transmittance is " << forward.refr << std::endl;
      }
      if(backward.refl < 0 || std::isnan(backward.refl) || std::isinf(backward.refl))
      {      
        std::cout << "WARNING! Precomputed film internal reflectance is " << backward.refl << std::endl;
      }
      if(backward.refr < 0 || std::isnan(backward.refr) || std::isinf(backward.refr))
      {
        std::cout << "WARNING! Precomputed film internal transmittance is " << backward.refr << std::endl;
      }
      #endif
    }
  }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  std::cout << duration.count() << std::endl;
  save_to_file("../precomputed_film_refl_ext.txt", res.ext_reflectivity.data(), FILM_LENGTH_RES, FILM_ANGLE_RES);
  save_to_file("../precomputed_film_refl_int.txt", res.int_reflectivity.data(), FILM_LENGTH_RES, FILM_ANGLE_RES);
  save_to_file("../precomputed_film_refr_ext.txt", res.ext_transmittivity.data(), FILM_LENGTH_RES, FILM_ANGLE_RES);
  save_to_file("../precomputed_film_refr_int.txt", res.int_transmittivity.data(), FILM_LENGTH_RES, FILM_ANGLE_RES);
  return res;
}

ThinFilmPrecomputed precomputeThinFilmRGB(
        const float extIOR, const uint* eta_id_vec, const uint* k_id_vec, const std::vector<float> &spec_values, 
        const std::vector<uint2> &spec_offsets, const float* eta_vec, const float* k_vec, const float* a_thickness, int layers, 
        const std::vector<float> &m_cie_x, const std::vector<float> &m_cie_y, const std::vector<float> &m_cie_z, 
        const uint thickness_res, const float thickness_min, const float thickness_max)
{
  ThinFilmPrecomputed res;
  res.ext_reflectivity.resize(FILM_ANGLE_RES * 3 * thickness_res);
  res.ext_transmittivity.resize(FILM_ANGLE_RES * 3 * thickness_res);
  res.int_reflectivity.resize(FILM_ANGLE_RES * 3 * thickness_res);
  res.int_transmittivity.resize(FILM_ANGLE_RES * 3 * thickness_res);

  for (uint t = 0; t < thickness_res; ++t)
  {
    float thickness = thickness_res == 1 ? a_thickness[0] : (thickness_max - thickness_min) / (thickness_res - 1) * t + thickness_min;

    spec::BasicSpectrum spec_refl_ext[FILM_ANGLE_RES];
    spec::BasicSpectrum spec_refr_ext[FILM_ANGLE_RES];
    spec::BasicSpectrum spec_refl_int[FILM_ANGLE_RES];
    spec::BasicSpectrum spec_refr_int[FILM_ANGLE_RES];

    for (size_t i = 0; i < FILM_LENGTH_RES; ++i)
    {
      float wavelength = (LAMBDA_MAX - LAMBDA_MIN) / (FILM_LENGTH_RES - 1) * i + LAMBDA_MIN;
      std::vector<complex> ior;
      ior.reserve(layers + 1);
      ior[0] = complex(extIOR, 0.f);
      float eta, k;
      uint2 data;
      uint offset;
      uint size;
      for (size_t layer = 0; layer < layers; ++layer)
      {
        eta = eta_vec[layer];
        uint eta_id = eta_id_vec[layer];
        if (eta_id < 0xFFFFFFFF)
        {
          data  = spec_offsets[eta_id];
          offset = data.x;
          size   = data.y;
          eta = SampleUniformSpectrum(spec_values.data() + offset, {wavelength, 0, 0, 0}, size)[0];
        }

        k = k_vec[layer];
        uint k_id = k_id_vec[layer];
        if (k_id < 0xFFFFFFFF)
        {
          data  = spec_offsets[k_id];
          offset = data.x;
          size   = data.y;
          k = SampleUniformSpectrum(spec_values.data() + offset, {wavelength, 0, 0, 0}, size)[0];
        }

        ior[layer + 1] = complex(eta, k);
      }
      for (int j = 0; j < FILM_ANGLE_RES; ++j)
      {
        float theta = M_PI_2 / float(FILM_ANGLE_RES - 1) * j;
        
        float cosTheta = clamp(cosf(theta), 1e-3f, 1.f);
        FrReflRefr forward;
        FrReflRefr backward;
        if (layers == 2)
        {
          forward = FrFilm(cosTheta, ior[0], ior[1], ior[2], thickness, wavelength);
          backward = FrFilm(cosTheta, ior[2], ior[1], ior[0], thickness, wavelength);
          //forward = TransferMatrixForward(cosTheta, ior.data(), a_thickness, layers, wavelength);
          //backward = TransferMatrixBackward(cosTheta, ior.data(), a_thickness, layers, wavelength);
        }
        else
        {
          forward = multFrFilm(cosTheta, ior.data(), a_thickness, layers, wavelength);
          backward = multFrFilm_r(cosTheta, ior.data(), a_thickness, layers, wavelength);
          //forward = TransferMatrixForward(cosTheta, ior.data(), a_thickness, layers, wavelength);
          //backward = TransferMatrixBackward(cosTheta, ior.data(), a_thickness, layers, wavelength);
        }

        if(forward.refl < 0 || std::isnan(forward.refl) || std::isinf(forward.refl))
        {
          std::cout << "WARNING! Precomputed film external reflectance is " << forward.refl << std::endl;
        }
        if(forward.refr < 0 || std::isnan(forward.refr) || std::isinf(forward.refr))
        {
          std::cout << "WARNING! Precomputed film external transmittance is " << forward.refr << std::endl;
        }
        if(backward.refl < 0 || std::isnan(backward.refl) || std::isinf(backward.refl))
        {      
          std::cout << "WARNING! Precomputed film internal reflectance is " << backward.refl << std::endl;
        }
        if(backward.refr < 0 || std::isnan(backward.refr) || std::isinf(backward.refr))
        {
          std::cout << "WARNING! Precomputed film internal transmittance is " << backward.refr << std::endl;
        }

        spec_refl_ext[j].set(wavelength, forward.refl);
        spec_refr_ext[j].set(wavelength, forward.refr);
        spec_refl_int[j].set(wavelength, backward.refl);
        spec_refr_int[j].set(wavelength, backward.refr);
      }
    }

    for (uint i = 0; i < FILM_ANGLE_RES; ++i)
    {
      auto rgb = spec::xyz2rgb(spec::spectre2xyz(spec_refl_ext[i]));
      res.ext_reflectivity[(FILM_ANGLE_RES * t + i) * 3]     = rgb.x;
      res.ext_reflectivity[(FILM_ANGLE_RES * t + i) * 3 + 1] = rgb.y;
      res.ext_reflectivity[(FILM_ANGLE_RES * t + i) * 3 + 2] = rgb.z;

      rgb = spec::xyz2rgb(spec::spectre2xyz(spec_refr_ext[i]));
      res.ext_transmittivity[(FILM_ANGLE_RES * t + i) * 3]     = rgb.x;
      res.ext_transmittivity[(FILM_ANGLE_RES * t + i) * 3 + 1] = rgb.y;
      res.ext_transmittivity[(FILM_ANGLE_RES * t + i) * 3 + 2] = rgb.z;

      rgb = spec::xyz2rgb(spec::spectre2xyz(spec_refl_int[i]));
      res.int_reflectivity[(FILM_ANGLE_RES * t + i) * 3]     = rgb.x;
      res.int_reflectivity[(FILM_ANGLE_RES * t + i) * 3 + 1] = rgb.y;
      res.int_reflectivity[(FILM_ANGLE_RES * t + i) * 3 + 2] = rgb.z;

      rgb = spec::xyz2rgb(spec::spectre2xyz(spec_refr_int[i]));
      res.int_transmittivity[(FILM_ANGLE_RES * t + i) * 3]     = rgb.x;
      res.int_transmittivity[(FILM_ANGLE_RES * t + i) * 3 + 1] = rgb.y;
      res.int_transmittivity[(FILM_ANGLE_RES * t + i) * 3 + 2] = rgb.z;
    }
  }

  return res;
}

std::vector<float> precomputeConductorFourier(
        uint eta_id, uint k_id, float eta, float k,
        const std::vector<float> &spec_values, const std::vector<uint2> &spec_offsets)
{
  std::vector<float> res;
  res.resize(FOURIER_ANGLE_RES * FourierSpec::SIZE);

  if (eta_id < 0xFFFFFFFF && k_id < 0xFFFFFFFF)
  {
    uint2 data  = spec_offsets[eta_id];
    uint offset = data.x;
    uint size   = data.y;
    FourierSpec etaSpecF = FourierSpec(spec_values.data() + offset, size);

    data   = spec_offsets[k_id];
    offset = data.x;
    size   = data.y;
    FourierSpec kSpecF = FourierSpec(spec_values.data() + offset, size);

    auto etaSpec = fourier::to_std_spectrum(etaSpecF);
    auto kSpec   = fourier::to_std_spectrum(kSpecF);

    auto std_samples = etaSpec.size();

    std::vector<float> spec;
    spec.resize(std_samples);

    for (uint i = 0; i < FILM_ANGLE_RES; ++i)
    {
      float theta = M_PI_2 / float(FILM_ANGLE_RES - 1) * i;
      float cosTheta = clamp(cosf(theta), 1e-6f, 1.f);

      for (uint j = 0; j < std_samples; ++j)
      {
        spec[j] = FrComplexConductor(cosTheta, {etaSpec[j], kSpec[j]});
      }

      FourierSpec specF = fourier::std_spectrum_to_fourier(spec);
      
      std::copy(specF.v, specF.v + FourierSpec::SIZE, res.begin() + i * FourierSpec::SIZE);
    }
  }

  // TO-DO constant ior

  return res;
}

ThinFilmPrecomputed precomputeThinFilmFourier(
        const float extIOR, const uint* eta_id_vec, const uint* k_id_vec, const std::vector<float> &spec_values, 
        const std::vector<uint2> &spec_offsets, const float* eta_vec, const float* k_vec,
        const float* a_thickness, int layers)
{
  ThinFilmPrecomputed res;
  res.ext_reflectivity.resize(FILM_ANGLE_RES * FourierSpec::SIZE);
  res.ext_transmittivity.resize(FILM_ANGLE_RES * FourierSpec::SIZE);
  res.int_reflectivity.resize(FILM_ANGLE_RES * FourierSpec::SIZE);
  res.int_transmittivity.resize(FILM_ANGLE_RES * FourierSpec::SIZE);

  std::vector<std::vector<float>> etaSpecs;
  etaSpecs.resize(layers);
  std::vector<std::vector<float>> kSpecs;
  kSpecs.resize(layers);

  ////
  for (size_t layer = 0; layer < layers; ++layer)
  {
    uint eta_id = eta_id_vec[layer];
    if (eta_id < 0xFFFFFFFF)
    {
      uint2 data  = spec_offsets[eta_id];
      uint offset = data.x;
      uint size   = data.y;
      FourierSpec etaSpecF = FourierSpec(spec_values.data() + offset, size);
      etaSpecs[layer] = fourier::to_std_spectrum(etaSpecF);
    }

    uint k_id = k_id_vec[layer]; 
    if (k_id < 0xFFFFFFFF)
    {   
      uint2 data   = spec_offsets[k_id];
      uint offset = data.x;
      uint size   = data.y;
      FourierSpec kSpecF = FourierSpec(spec_values.data() + offset, size);
      kSpecs[layer] = fourier::to_std_spectrum(kSpecF);
    }
  }
  ///

  ThinFilmPrecomputed realSpec;
  realSpec.ext_reflectivity.resize(nCIESamples);
  realSpec.ext_transmittivity.resize(nCIESamples);
  realSpec.int_reflectivity.resize(nCIESamples);
  realSpec.int_transmittivity.resize(nCIESamples);
    
  for (int j = 0; j < FILM_ANGLE_RES; ++j)
  {
    float theta = M_PI_2 / float(FILM_ANGLE_RES - 1) * j;
    float cosTheta = clamp(cosf(theta), 1e-3f, 1.f);

    for (size_t i = 0; i < nCIESamples; ++i)
    {
      float wavelength = (LAMBDA_MAX - LAMBDA_MIN - 1) / (nCIESamples - 1) * i + LAMBDA_MIN;
      std::vector<complex> ior;
      ior.reserve(layers + 1);
      ior[0] = complex(extIOR, 0.f);
      float eta, k;
      uint2 data;
      uint offset;
      uint size;
      for (size_t layer = 0; layer < layers; ++layer)
      {
        eta = eta_vec[layer];
        uint eta_id = eta_id_vec[layer];
        if (eta_id < 0xFFFFFFFF)
          eta = etaSpecs[layer][i];

        k = k_vec[layer];
        uint k_id = k_id_vec[layer];
        if (k_id < 0xFFFFFFFF)
          k = kSpecs[layer][i];

        ior[layer + 1] = complex(eta, k);
      }

      FrReflRefr forward;
      FrReflRefr backward;
      if (layers == 2)
      {
        forward = FrFilm(cosTheta, ior[0], ior[1], ior[2], a_thickness[0], wavelength);
        backward = FrFilm(cosTheta, ior[2], ior[1], ior[0], a_thickness[0], wavelength);
        //forward = TransferMatrixForward(cosTheta, ior.data(), a_thickness, layers, wavelength);
        //backward = TransferMatrixBackward(cosTheta, ior.data(), a_thickness, layers, wavelength);
      }
      else
      {
        forward = multFrFilm(cosTheta, ior.data(), a_thickness, layers, wavelength);
        backward = multFrFilm_r(cosTheta, ior.data(), a_thickness, layers, wavelength);
        //forward = TransferMatrixForward(cosTheta, ior.data(), a_thickness, layers, wavelength);
        //backward = TransferMatrixBackward(cosTheta, ior.data(), a_thickness, layers, wavelength);
      }
      realSpec.ext_reflectivity[i] = forward.refl;
      realSpec.ext_transmittivity[i] = forward.refr;
      realSpec.int_reflectivity[i] = backward.refl;
      realSpec.int_transmittivity[i] = backward.refr;  
    }

    FourierSpec extReflF = fourier::std_spectrum_to_fourier(realSpec.ext_reflectivity);
    FourierSpec extRefrF = fourier::std_spectrum_to_fourier(realSpec.ext_transmittivity);
    FourierSpec intReflF = fourier::std_spectrum_to_fourier(realSpec.int_reflectivity);
    FourierSpec intRefrF = fourier::std_spectrum_to_fourier(realSpec.int_transmittivity);

    /*
    std::cout << "______________________________________________________________________________________________________________\n";
    std::cout << "Film: ";
    for (int i = 0; i < FourierSpec::SIZE; ++i) 
    {
      std::cout << extReflF[i] << " " << extRefrF[i] << " ";
    }
    std::cout << std::endl;
    */

    std::copy(extReflF.v, extReflF.v + FourierSpec::SIZE, res.ext_reflectivity.begin() + j * FourierSpec::SIZE);
    std::copy(extRefrF.v, extRefrF.v + FourierSpec::SIZE, res.ext_transmittivity.begin() + j * FourierSpec::SIZE);
    std::copy(intReflF.v, intReflF.v + FourierSpec::SIZE, res.int_reflectivity.begin() + j * FourierSpec::SIZE);
    std::copy(intRefrF.v, intRefrF.v + FourierSpec::SIZE, res.int_transmittivity.begin() + j * FourierSpec::SIZE);
  }

  return res;
}