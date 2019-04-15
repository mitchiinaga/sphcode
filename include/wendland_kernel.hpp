#pragma once

#include <cmath>

#include "defines.hpp"
#include "kernel_function.hpp"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971693993751
#endif

// Wendland (1995), Dehnen & Aly (2012)
namespace sph
{
namespace wendland
{

// Wendland C4 Kernel
#if DIM == 2
    constexpr real sigma = 9.0 / M_PI;
#else
    constexpr real sigma = 495. / (32 * M_PI);
#endif

class C4Kernel : public KernelFunction {
public:
    real w(const real r, const real h)
    {
        const real q = r / h;
        return sigma / powh(h) * pow6(0.5 * (1.0 - q + std::abs(1.0 - q))) * (1.0 + 6.0 * q + 35.0 / 3.0 * q * q);
    }

    vec_t dw(const vec_t &rij, const real r, const real h)
    {
        const real q = r / h;
        const real cc = -56.0 / 3.0 * sigma / (powh(h) * sqr(h)) * pow5(0.5 * (1.0 - q + std::abs(1.0 - q))) * (1.0 + 5.0 * q);
        return rij * cc;
    }
};

}
}
