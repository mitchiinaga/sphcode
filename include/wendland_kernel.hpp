#pragma once

#include <cmath>
#include <cassert>

#include "defines.hpp"
#include "kernel_function.hpp"

// Wendland (1995), Dehnen & Aly (2012)
namespace sph
{
namespace Wendland {

// Wendland C4 Kernel
#if DIM == 1
    constexpr real sigma_c4 = 0.0;
#elif DIM == 2
    constexpr real sigma_c4 = 9.0 / M_PI;
#else
    constexpr real sigma_c4 = 495. / (32 * M_PI);
#endif

class C4Kernel : public KernelFunction {
public:
    C4Kernel()
    {
        assert(DIM != 1);
    }

    real w(const real r, const real h)
    {
        const real q = r / h;
        return sigma_c4 / powh(h) * pow6(0.5 * (1.0 - q + std::abs(1.0 - q))) * (1.0 + 6.0 * q + 35.0 / 3.0 * q * q);
    }

    vec_t dw(const vec_t &rij, const real r, const real h)
    {
        const real q = r / h;
        const real c = -56.0 / 3.0 * sigma_c4 / (powh(h) * sqr(h)) * pow5(0.5 * (1.0 - q + std::abs(1.0 - q))) * (1.0 + 5.0 * q);
        return rij * c;
    }

    real dhw(const real r, const real h)
    {
        const double q = r / h;
        return -sigma_c4 / (powh(h) * h * 3.0) * pow5(0.5 * (1.0 - q + std::abs(1.0 - q)))
            * (3.0 * DIM + 15.0 * DIM * q + (-56.0 + 17.0 * DIM) * q * q - 35.0 * (8.0 + DIM) * pow3(q));
    }
};

}
}
