#include "fourier.h"
#include "levinson.h"
#include <cassert>
#include "LiteMath.h"


#include <iostream>

namespace fourier
{

    using Complex = std::complex<float>;
    static constexpr Complex I{0.0f, 1.0f};

    float to_phase(float wl, float start, float end)
    {
        return std::fma(LiteMath::M_PI, (wl - start) / (end - start), -LiteMath::M_PI);
    }

    FourierSpec spectrum_to_fourier(const std::vector<float> &phases, const std::vector<float> &values, int n)
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

}