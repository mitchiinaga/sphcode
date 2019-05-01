#pragma once

#include <memory>

#include "vector_type.hpp"
#include "parameters.hpp"

namespace sph
{

class Distance {
public:
    void initialize(std::shared_ptr<SPHParameters> param)
    {
    }

    vec_t calc_r_ij(const vec_t & r_i, const vec_t & r_j) const
    {
        return r_i - r_j;
    }
};

}