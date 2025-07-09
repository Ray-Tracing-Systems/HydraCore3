#include "fourier.h"
#include "levinson.h"
#include <cassert>
#include "LiteMath.h"

#include "../spectrum.h"



namespace fourier
{
    ffunc_t fourier_function = fourier_series;
    using Complex = std::complex<float>;
    static constexpr Complex I{0.0f, 1.0f};


    static std::vector<float> compute_lut()
    {   
        std::vector<float> res;
        const std::vector<float> &wavelenghts = Get_CIE_lambda();
        res.reserve(wavelenghts.size() * FourierSpec::SIZE);
        
        for(int i = 0; i < wavelenghts.size(); ++i) {
            for(int j = 0; j < FourierSpec::SIZE; ++j) {
                float phase = to_phase(wavelenghts[i], LAMBDA_MIN, LAMBDA_MAX);
                res.push_back(std::cos(j * phase));
            }
        }
        return res;
    }

    static const std::vector<float> FOURIER_PRECOMP = compute_lut();

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
        data[FourierSpec::M] = LiteMath::INV_PI * spec[0];
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

    std::vector<float> mese(const std::vector<float> &phases, const FourierSpec &spec)
    {
        std::vector<float> q = precompute_mese_coeffs(spec);

        std::vector<float> res(phases.size());

        for(unsigned k = 0; k < phases.size(); ++k) {
            res[k] = 0.5f * mese_precomp(phases[k], q);
        }
        return res;
    }

    static float fourier_series(float phase, const FourierSpec &spec)
    {
        float res = 0.0f;
        for(int i = 0; i < FourierSpec::SIZE; ++i) {
            res += spec[i] * std::cos(phase * i);
        }
        return res;
    }

    std::vector<float> fourier_series(const std::vector<float> &phases, const FourierSpec &spec)
    {
        std::vector<float> res(phases.size());
        for(size_t i = 0; i < phases.size(); ++i) {
            res[i] = fourier_series(phases[i], spec);
        }
        return res;
    }

    std::vector<float> mese(const FourierSpec &spec) { return mese(wl_to_phases(Get_CIE_lambda()), spec); }
    std::vector<float> fourier_series(const FourierSpec &spec) { return fourier_series(wl_to_phases(Get_CIE_lambda()), spec); }

    std::vector<float> fourier_series_lut(const FourierSpec &spec)
    {   
        //compute_lut_if_needed();
        static std::vector<float> phases = wl_to_phases(Get_CIE_lambda());

        std::vector<float> result(phases.size());

        for(int i = 0; i < phases.size(); ++i) {
            const size_t offset = i * FourierSpec::SIZE;
            float val = 0.0f;

            for(int j = 0; j < FourierSpec::SIZE; ++j) {
                val += FOURIER_PRECOMP[offset + j] * spec[j];
            }
            result[i] = val;
        }
        return result;
    }

    FourierSpec std_spectrum_to_fourier(const std::vector<float> &values)
    {
        return real_fourier_moments_of(wl_to_phases(Get_CIE_lambda()), values);
    }

    std::vector<float> to_std_spectrum(const FourierSpec &spec)
    {
        std::vector<float> phases = wl_to_phases(Get_CIE_lambda());

        return fourier_function(spec);
    }

    void set_calc_func(ffunc_t function)
    {
        fourier_function = function;
    }

}