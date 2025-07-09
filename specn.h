#ifndef SPECN_H_
#define SPECN_H_

#include <include/cglobals.h>

constexpr int SPECN_SIZE = 10;

template<int N>
struct SpecN_
{
    static constexpr int SIZE = N;
 
    float v[SIZE];

    SpecN_() : v{} 
    {
        for(int i = 0; i < SIZE; ++i) {
            v[i] = 0.0f;
        }
    }

    explicit SpecN_(const float *ptr) : v{} 
    { 
        for(int i = 0; i < SIZE; ++i) {
            v[i] = ptr[i];
        }
    }

    SpecN_(const SpecN_ &other) = default;

    explicit SpecN_(float x) : v{} {
        for(int i = 0; i < SIZE; ++i) {
            v[i] = x;
        }
    }

    explicit SpecN_(int x) : SpecN_(float(x)) {}

    SpecN_ &operator=(const SpecN_ &other) = default;

    SpecN_ operator+(const SpecN_ &other) const 
    {
        SpecN_ c(*this);
        c += other;
        return c;
    }

    SpecN_ operator-(const SpecN_ &other) const 
    {
        SpecN_ c(*this);
        c -= other;
        return c;
    }

    SpecN_ operator*(float f) const 
    {
        SpecN_ c(*this);
        c *= f;
        return c;
    }

    SpecN_ operator/(float f) const 
    {
        SpecN_ c(*this);
        c /= f;
        return c;
    }
    //Convolution in time domain
    SpecN_ operator*(const SpecN_ &other) const
    {
        SpecN_ c(*this);
        c *= other;
        return c;
    }

    float &operator[](int i) 
    {
        return v[i];
    }

    float operator[](int i) const 
    {
        return v[i];
    }

    SpecN_ operator-() const 
    {
        SpecN_ res;
        for(uint i = 0; i < SIZE; ++i) {
            res.v[i] = -v[i];
        }
        return res;
    }

    SpecN_ &operator+=(const SpecN_ &other) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] += other.v[i];
        }
        return *this;
    }

    SpecN_ &operator-=(const SpecN_ &other) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] -= other.v[i];
        }
        return *this;
    }

    SpecN_ &operator*=(float f) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] *= f;
        }
        return *this;
    }

    SpecN_ &operator/=(float f) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] /= f;
        }
        return *this;
    }

    SpecN_ &operator*=(const SpecN_ &other) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] *= other.v[i];
        }
        return *this;
    }

    bool operator==(const SpecN_ &other) const 
    {
        for(int i = 0; i < SIZE; ++i) {
            if(v[i] != other.v[i]) return false;
        }
        return true;
    }

    bool operator!=(const SpecN_ &other) const 
    {
        for(int i = 0; i < SIZE; ++i) {
            if(v[i] != other.v[i]) return true;
        }
        return false;
    }


};

template<int N>
inline SpecN_<N> operator*(float f, const SpecN_<N> &spec)
{
    return spec * f;
}

template<typename Spec>
struct BsdfEvalN_
{
  Spec val;
  float  pdf; 
};

template<typename Spec>
struct BsdfSampleN_
{
  Spec val;
  float3 dir;
  float  pdf; 
  uint   flags;
  float  ior;
};



using SpecN = SpecN_<SPECN_SIZE>;
using BsdfSampleN = BsdfSampleN_<SpecN>;
using BsdfEvalN = BsdfEvalN_<SpecN>;

static inline SpecN SampleWavelengthsN(float u, float a, float b) 
{
  // pdf is 1.0f / (b - a)
  SpecN res;

  res[0] = lerp(a, b, u);

  float delta = (b - a) / float(SpecN::SIZE);
  for (uint32_t i = 1; i < SpecN::SIZE; ++i) 
  {
      res[i] = res[i - 1] + delta;
      if (res[i] > b)
        res[i] = a + (res[i] - b);
  }

  return res;
}

static inline SpecN SampleUniformSpectrumN(const float* a_spec_values, SpecN a_wavelengths, uint32_t a_sz)
{
  SpecN res; 
  for(int i = 0; i < SpecN::SIZE; ++i) {
    size_t idx = size_t(a_wavelengths[i] - LAMBDA_MIN);
    res[i] = a_spec_values[idx];
  }

  return res;
}

#endif