#ifndef FOURIER_FOURIER_H_
#define FOURIER_FOURIER_H_
#include "LiteMath.h"
#include <complex>
#include <vector>
#include "fspec.h"

//CPU-only

namespace fourier
{

    float to_phase(float wl, float start, float end);

    inline std::vector<float> wl_to_phases(std::vector<float> wavelenghts, float start, float end)
    {
        for(unsigned i = 0; i < wavelenghts.size(); ++i) {
            wavelenghts[i] = to_phase(wavelenghts[i], start, end);
        }
        return wavelenghts;
    }

    FourierSpec spectrum_to_fourier(const std::vector<float> &phases, const std::vector<float> &values, int n);


    std::vector<float> precompute_mese_coeffs(const FourierSpec &spec);

    float mese_precomp(float phase, const std::vector<float> &q);

    inline std::vector<float> mese(std::vector<float> phases, const FourierSpec &spec)
    {
        std::vector<float> q = precompute_mese_coeffs(spec);

        for(unsigned k = 0; k < phases.size(); ++k) {
            phases[k] = mese_precomp(phases[k], q);
        }
        return phases;
    }

    inline float mese(float phase, const FourierSpec &spec)
    {
        return mese_precomp(phase, precompute_mese_coeffs(spec));
    }


}


#endif