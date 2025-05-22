#ifndef FOURIER_FOURIER_H_
#define FOURIER_FOURIER_H_
#include "LiteMath.h"
#include <complex>
#include <vector>
#include "vecn.h"

namespace fourier
{   
    constexpr uint FOURIER_M = 8;
    constexpr uint VEC_SIZE  = FOURIER_M + 1; 

    using Complex = std::complex<float>;
    using fvec = vec<float, VEC_SIZE>;
    using fvec_c = vec<Complex, VEC_SIZE>;

    constexpr Complex I{0.0f, 1.0f};

    float to_phase(float wl, float start, float end);

    inline std::vector<float> wl_to_phases(std::vector<float> wavelenghts, float start, float end)
    {
        for(unsigned i = 0; i < wavelenghts.size(); ++i) {
            wavelenghts[i] = to_phase(wavelenghts[i], start, end);
        }
        return wavelenghts;
    }

    fvec real_fourier_moments_of(const std::vector<float> &phases, const std::vector<float> &values, int n);

    fvec_c fourier_moments_of(const std::vector<float> &phases, const std::vector<float> &values, int n);

    fvec_c precompute_mese_coeffs(const fvec_c &gamma);
    fvec precompute_mese_coeffs(const fvec &gamma);

    fvec_c lagrange_multipliers(const fvec &moments);

    float bounded_mese_l(float phase, const fvec_c &lagrange_m);

    inline float bounded_mese_m(float phase, const fvec &moments)
    {
        return bounded_mese_l(phase, lagrange_multipliers(moments));
    }

    float mese_precomp(float phase, const fvec &moments);

    inline std::vector<float> mese(std::vector<float> phases, const fvec &moments)
    {
        fvec q = precompute_mese_coeffs(moments);

        for(unsigned k = 0; k < phases.size(); ++k) {
            phases[k] = 0.5f * mese_precomp(phases[k], q);
        }
        return phases;
    }

    inline float mese(float phase, const fvec &moments)
    {
        return 0.5f * mese_precomp(phase, precompute_mese_coeffs(moments));
    }


}


#endif