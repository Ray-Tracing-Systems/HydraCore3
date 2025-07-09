#ifndef FOURIER_FOURIER_H_
#define FOURIER_FOURIER_H_
#include "LiteMath.h"
#include "include/cglobals.h"
#include <complex>
#include <vector>
#include "fspec.h"


//CPU-only

namespace fourier
{

    using ffunc_t = std::vector<float> (*)(const FourierSpec &);

    
    FourierSpec real_fourier_moments_of(const std::vector<float> &phases, const std::vector<float> &values);

    inline FourierSpec spectrum_to_fourier(const std::vector<float> &wavelenghts, const std::vector<float> &values)
    {
        return real_fourier_moments_of(wl_to_phases(wavelenghts), values);
    }

    FourierSpec std_spectrum_to_fourier(const std::vector<float> &values);

    std::vector<float> precompute_mese_coeffs(const FourierSpec &spec);

    float mese_precomp(float phase, const std::vector<float> &q);

    std::vector<float> mese(const std::vector<float> &phases, const FourierSpec &spec);

    std::vector<float> fourier_series(const std::vector<float> &phases, const FourierSpec &spec);

    std::vector<float> mese(const FourierSpec &spec);

    std::vector<float> fourier_series(const std::vector<float> &phases, const FourierSpec &spec);
    std::vector<float> fourier_series(const FourierSpec &spec);
    std::vector<float> fourier_series_lut(const FourierSpec &spec);



    std::vector<float> to_std_spectrum(const FourierSpec &spec);

    void set_calc_func(ffunc_t function);


    extern std::vector<float> fourier_series(const std::vector<float> &phases, const FourierSpec &spec);



}


#endif