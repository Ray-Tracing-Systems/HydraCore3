#include "fourier.h"
#include "levinson.h"
#include <cassert>
#include "LiteMath.h"
#include <numeric>


#include <iostream>

namespace fourier
{
    ffunc_t fourier_function = fourier_series;


    using Complex = std::complex<float>;
    static constexpr Complex I{0.0f, 1.0f};

    float to_phase(float wl, float start, float end)
    {
        return std::fma(LiteMath::M_PI, (wl - start) / (end - start), -LiteMath::M_PI);
    }

    FourierSpec real_fourier_moments_of(const std::vector<float> &phases, const std::vector<float> &values)
    {
        const size_t N = phases.size();
        assert(N == values.size());
        FourierSpec spec;
        for(int i = 0; i <= FourierSpec::M; ++i) {
            float val = 0.0f;
            for(size_t j = 0; j < N; ++j) {
                val += values[j] * std::cos(float(i) * phases[j]);
            }
            spec[i] = 2.0f * val / float(N);
        }
        spec[0] *= 0.5f;
        return spec;
    }

    std::vector<float> precompute_mese_coeffs(const FourierSpec &spec)
    {
        std::vector<float> e0(FourierSpec::M + 1, 0.0f);
        e0[0] = 1.0f;
        std::vector<float> data(2 * FourierSpec::M + 1); //Matrix values
        data[FourierSpec::M] = LiteMath::INV_TWOPI * spec[0];
        for(int i = 1; i <= FourierSpec::M; ++i) {
            data[FourierSpec::M + i] = LiteMath::INV_TWOPI * spec[i];
            data[FourierSpec::M - i] = LiteMath::INV_TWOPI * spec[i]; //[std::conj(spec[i])] in fact;
        }
        return levinson<float>(data, e0);
    }


    float mese_precomp(float phase, const std::vector<float> &q)
    {
        Complex t = 0.0f;
        for(int i = 0; i <= FourierSpec::M; ++i) t += LiteMath::INV_TWOPI * q[i] * std::exp(-I * float(i) * phase); 

        float div = std::fabs(t);
        div *= div;

        return (LiteMath::INV_TWOPI * std::real(q[0])) / div;
    }

    std::vector<float> mese(std::vector<float> phases, const FourierSpec &spec)
    {
        std::vector<float> q = precompute_mese_coeffs(spec);

        for(unsigned k = 0; k < phases.size(); ++k) {
            phases[k] = 0.5f * mese_precomp(phases[k], q);
        }
        return phases;
    }

    std::vector<float> fourier_series(std::vector<float> phases, const FourierSpec &spec)
    {
        return {};
    }


    void to_std_spectrum(const FourierSpec &spec, float *out_spectrum)
    {
        std::vector<float> phases(size_t(LAMBDA_MAX - LAMBDA_MIN + 1));
        std::iota(phases.begin(), phases.end(), LAMBDA_MIN);
        phases = wl_to_phases(phases, LAMBDA_MIN, LAMBDA_MAX);

        std::vector<float> spectrum = fourier_function(phases, spec);
        std::copy(spectrum.begin(), spectrum.end(), out_spectrum);

    }

    void set_calc_func(ffunc_t function)
    {
        fourier_function = function;
    }

}