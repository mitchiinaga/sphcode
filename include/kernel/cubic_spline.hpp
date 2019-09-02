#pragma once

#include <cmath>

#include "defines.hpp"
#include "kernel_function.hpp"

// cubic spline kernel
namespace sph
{
namespace Spline
{
#if DIM == 1
    constexpr real sigma_cubic = 2.0 / 3.0;
#elif DIM == 2
    constexpr real sigma_cubic = 10.0 / (7.0 * M_PI);
#else
    constexpr real sigma_cubic = 1.0 / M_PI;
#endif

class Cubic : public KernelFunction {
public:
    Cubic()
    {
    }

    real w(const real r, const real h) const
    {
        const real h_ = h * 0.5;
        const real q = r / h_;
        return sigma_cubic / powh(h_) * (0.25 * pow3(0.5 * (2.0 - q + std::abs(2.0 - q))) - pow3(0.5 * (1.0 - q + std::abs(1.0 - q))));
    }

    vec_t dw(const vec_t &rij, const real r, const real h) const
    {
        if(r == 0.0) {
            return vec_t(0);
        }
        const real h_ = h * 0.5;
        const real q = r / h_;
        const real c = -sigma_cubic / (powh(h_) * h_ * r) * (0.75 * sqr(0.5 * (2.0 - q + std::abs(2.0 - q))) - 3.0 * sqr(0.5 * (1.0 - q + std::abs(1.0 - q))));
        return rij * c;
    }

    real dhw(const real r, const real h) const
    {
        const real h_ = h * 0.5;
        const real q = r / h_;
        return 0.5 * sigma_cubic / (powh(h_) * h_) * (sqr((std::abs(2.0 - q) + 2.0 - q) * 0.5) * ((3. + DIM) * 0.25 * q - 0.5 * DIM)
            + sqr((std::abs(1.0 - q) + 1.0 - q) * 0.5) * ((-3.0 - DIM) * q + DIM));
    }
};
}

}
