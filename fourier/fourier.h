#ifndef FOURIER_FOURIER_H_
#define FOURIER_FOURIER_H_
#include "LiteMath.h"
#include <complex>
#include <vector>

namespace fourier
{

    using Complex = std::complex<float>;

    constexpr Complex I{0.0f, 1.0f};

    float to_phase(float wl, float start, float end);

    inline std::vector<float> wl_to_phases(std::vector<float> wavelenghts, float start, float end)
    {
        for(unsigned i = 0; i < wavelenghts.size(); ++i) {
            wavelenghts[i] = to_phase(wavelenghts[i], start, end);
        }
        return wavelenghts;
    }

    std::vector<float> real_fourier_moments_of(const std::vector<float> &phases, const std::vector<float> &values, int n);

    std::vector<Complex> fourier_moments_of(const std::vector<float> &phases, const std::vector<float> &values, int n);

    std::vector<Complex> precompute_mese_coeffs(const std::vector<Complex> &gamma);
    std::vector<float> precompute_mese_coeffs(const std::vector<float> &gamma);

    std::vector<Complex> lagrange_multipliers(const std::vector<float> &moments);

    float bounded_mese_l(float phase, const std::vector<Complex> &lagrange_m);

    inline float bounded_mese_m(float phase, const std::vector<float> &moments)
    {
        return bounded_mese_l(phase, lagrange_multipliers(moments));
    }

    float mese_precomp(float phase, const std::vector<float> &moments);

    inline std::vector<float> mese(std::vector<float> phases, const std::vector<float> &moments)
    {
        std::vector<float> q = precompute_mese_coeffs(moments);

        for(unsigned k = 0; k < phases.size(); ++k) {
            phases[k] = mese_precomp(phases[k], q);
        }
        return phases;
    }

    inline float mese(float phase, const std::vector<float> &moments)
    {
        return mese_precomp(phase, precompute_mese_coeffs(moments));
    }


}


#endif