#include "fourier.h"
#include "levinson.h"
#include <cassert>
#include "LiteMath.h"
#include <numeric>

#include "../spectrum.h"


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
        float t = 0.0f;
        for(int i = 0; i <= FourierSpec::M; ++i) t += LiteMath::INV_TWOPI * q[i] * std::cos(float(i) * phase); 

        float div = std::fabs(t);
        div *= div;

        return (LiteMath::INV_TWOPI * q[0]) / div;
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
        for(int i = 0; i < phases.size(); ++i) {
            res[i] = fourier_series(phases[i], spec);
        }
        return res;
    }

    std::vector<float> mese(const FourierSpec &spec) { return mese(wl_to_phases(Get_CIE_lambda()), spec); }


    std::vector<Complex> exponential_moments(const FourierSpec &spec, float &gamma0)
    {
        gamma0 = 0.5f * INV_TWO_PI * std::exp(PI * I * (moments[0] - 0.5f));
        std::vector<Complex> res(M + 1);
        res[0] = 2.0f * gamma0.real();
        for(int i = 1; i <= FourierSpec::M; ++i) {
            Complex sm{0.0f, 0.0f};
            for(int j = 1; j <= i - 1; ++j) {
                sm += float(i - j) * res[j] * spec[i - j]; 
            }

            res[i] = TWO_PI * I * (float(i) * gamma0 * spec[i] + sm) / float(i);
        }
        return res;
    }

    std::vector<float> lagrange_multipliers(const FourierSpec &spec)
    {      
        Complex gamma0;
        std::vector<Complex> gamma = exponential_moments(spec, gamma0);

        //Calculating q
        std::vector<float> q = precompute_mese_coeffs(std::vector<float>(gamma.begin(), gamma.end()));


        std::vector<float> lambda(FourierSpec::SIZE);
        for(int i = 0; i <= FourierSpec::M; ++i) {
            Complex t{0.0f, 0.0f};

            for(int k = 0; k <= FourierSpec::M - i; ++k) {
                Complex ts{0.0f, 0.0f};
                for(int j = 0; j <= FourierSpec::M - k - i; ++j) {
                    ts += std::conj(q[j + k + i]) * q[j];
                }
                t += ts * (k == 0 ? gamma0 : gamma[k]);
            }
            lambda[i] = t / (PI * I * q[0]);
        }
        return std::real(lambda);
    }
    
    float bounded_mese_l(float phase, const std::vector<float> &lagrange_m)
    {
        const int M = lagrange_m.size() - 1;
        float sum = std::real(lagrange_m[0]);
        for(int i = 1; i <= M; ++i) {
            sum += 2.0f * std::real(lagrange_m[i] * std::cos(float(i) * phase));
        } 

        return INV_PI * std::atan(sum) + 0.5f;
    }



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