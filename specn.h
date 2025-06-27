#ifndef SPECN_H_
#define SPECN_H_
#include <algorithm>
#include <numeric>




struct SpecN
{
    static constexpr int SIZE = 10;
 
    float v[SIZE];

    SpecN() : v{} 
    {
        for(int i = 0; i < SIZE; ++i) {
            v[i] = 0.0f;
        }
    }

    explicit SpecN(const float *ptr) : v{} 
    { 
        for(int i = 0; i < SIZE; ++i) {
            v[i] = ptr[i];
        }
    }

    SpecN(const SpecN &other) = default;

    explicit SpecN(float x) : v{} {
        for(int i = 0; i < SIZE; ++i) {
            v[i] = x;
        }
    }

    explicit SpecN(int x) : SpecN(float(x)) {}

    SpecN &operator=(const SpecN &other) = default;

    SpecN operator+(const SpecN &other) const 
    {
        SpecN c(*this);
        c += other;
        return c;
    }

    SpecN operator-(const SpecN &other) const 
    {
        SpecN c(*this);
        c -= other;
        return c;
    }

    SpecN operator*(float f) const 
    {
        SpecN c(*this);
        c *= f;
        return c;
    }

    SpecN operator/(float f) const 
    {
        SpecN c(*this);
        c /= f;
        return c;
    }
    //Convolution in time domain
    SpecN operator*(const SpecN &other) const
    {
        SpecN c(*this);
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

    SpecN operator-() const 
    {
        SpecN res;
        for(uint i = 0; i < SIZE; ++i) {
            res.v[i] = -v[i];
        }
        return res;
    }

    SpecN &operator+=(const SpecN &other) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] += other.v[i];
        }
        return *this;
    }

    SpecN &operator-=(const SpecN &other) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] -= other.v[i];
        }
        return *this;
    }

    SpecN &operator*=(float f) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] *= f;
        }
        return *this;
    }

    SpecN &operator/=(float f) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] /= f;
        }
        return *this;
    }

    SpecN &operator*=(const SpecN &other) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] *= other.v[i];
        }
        return *this;
    }

    bool operator==(const SpecN &other) const 
    {
        for(int i = 0; i < SIZE; ++i) {
            if(v[i] != other.v[i]) return false;
        }
        return true;
    }

    bool operator!=(const SpecN &other) const 
    {
        for(int i = 0; i < SIZE; ++i) {
            if(v[i] != other.v[i]) return true;
        }
        return false;
    }


};

inline SpecN operator*(float f, const SpecN &spec)
{
    return spec * f;
}


struct BsdfEvalN
{
  SpecN val;
  float  pdf; 
};

struct BsdfSampleN
{
  SpecN val;
  float3 dir;
  float  pdf; 
  uint   flags;
  float  ior;
};

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
  const int  WAVESN = int(LAMBDA_MAX-LAMBDA_MIN);

  SpecN res;
  for(int i = 0; i < SpecN::SIZE; ++i) {
    size_t idx = size_t(a_wavelengths[i] - LAMBDA_MIN);
    res[i] = a_spec_values[idx];
  }

  return res;
}


#endif