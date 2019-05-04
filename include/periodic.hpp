#pragma once

#include <memory>
#include <cassert>

#include "vector_type.hpp"
#include "parameters.hpp"

namespace sph
{

class Periodic {
    bool m_is_valid;
    vec_t m_max;
    vec_t m_min;
    vec_t m_range;

public:
    Periodic() : m_is_valid(false) {}
    void initialize(std::shared_ptr<SPHParameters> param)
    {
        if(param->periodic.is_valid) {
            m_is_valid = true;
            for(int i = 0; i < DIM; ++i) {
                m_max[i] = param->periodic.range_max[i];
                m_min[i] = param->periodic.range_min[i];
            }
            m_range = m_max - m_min;
        } else {
            m_is_valid = false;
        }
    }

    vec_t calc_r_ij(const vec_t & r_i, const vec_t & r_j) const
    {
        if(m_is_valid) {
            const vec_t r_ij1 = r_i - r_j;
            const vec_t r_ij2 = r_ij1 + m_range;
            const vec_t r_ij3 = r_ij1 - m_range;
            vec_t r_ij;

            #define R_IJ(a, b, c)\
                if(std::abs(r_ij##a[i]) <= std::abs(r_ij##b[i]) && std::abs(r_ij##a[i]) <= std::abs(r_ij##c[i])) {\
                    r_ij[i] = r_ij##a[i];\
                    continue;\
                }

            for(int i = 0; i < DIM; ++i) {
                R_IJ(1, 2, 3);
                R_IJ(2, 3, 1);
                R_IJ(3, 1, 2);
                assert(false);
            }

            return r_ij;
        } else {
            return r_i - r_j;
        }
    }

    void apply(vec_t & r) const
    {
        if(m_is_valid) {
            for(int i = 0; i < DIM; ++i) {
                if(r[i] < m_min[i]) {
                    r[i] += m_range[i];
                } else if(r[i] > m_max[i]) {
                    r[i] -= m_range[i];
                }
            }
        }
    }
};

}