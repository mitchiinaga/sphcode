#pragma once

#include <cmath>

#include "defines.hpp"
#include "kernel_function.hpp"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971693993751
#endif

// cubic spline kernel
namespace sph
{

#if DIM == 1
    constexpr real sigma = 2.0 / 3.0;
#elif DIM == 2
    constexpr real sigma = 10.0 / (7.0 * M_PI);
#else
    constexpr real sigma = 1.0 / M_PI;
#endif

class CubicSpline : public KernelFunction {
public:
    CubicSpline()
    {
    }

    real w(const real r, const real h)
    {
        const real h_ = h * 0.5;
        const real q = r / h_;
        return sigma / powh(h_) * (0.25 * pow3(0.5 * (2.0 - q + std::abs(2.0 - q))) - pow3(0.5 * (1.0 - q + std::abs(1.0 - q))));
    }

    vec_t dw(const vec_t &rij, const real r, const real h)
    {
        const real h_ = h * 0.5;
        const real q = r / h_;
        const real c = -sigma / (powh(h_) * h_ * r) * (0.75 * sqr(0.5 * (2.0 - q + std::abs(2.0 - q))) - 3.0 * sqr(0.5 * (1.0 - q + std::abs(1.0 - q))));
        return rij * c;
    }
};

}
