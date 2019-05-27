#include <algorithm>

#include "vector_type.hpp"
#include "particle.hpp"
#include "periodic.hpp"
#include "exhaustive_search.hpp"

namespace sph {

// 全探索 (デバッグ用)
int exhaustive_search(
    SPHParticle & p_i,
    const real kernel_size,
    const std::vector<SPHParticle> & particles,
    const int num,
    std::vector<int> & neighbor_list,
    const int list_size,
    Periodic const * periodic,
    const bool is_ij)
{
    const real kernel_size_i2 = kernel_size * kernel_size;
    const vec_t & pos_i = p_i.pos;
    int count = 0;
    for(int j = 0; j < num; ++j) {
        const auto & p_j = particles[j];
        const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
        const real r2 = abs2(r_ij);
        const real kernel_size2 = is_ij ? std::max(kernel_size_i2, p_j.sml * p_j.sml) : kernel_size_i2;
        if(r2 < kernel_size2) {
            neighbor_list[count] = j;
            ++count;
        }
    }

    std::sort(neighbor_list.begin(), neighbor_list.begin() + count, [&](const int a, const int b) {
        const vec_t r_ia = periodic->calc_r_ij(pos_i, particles[a].pos);
        const vec_t r_ib = periodic->calc_r_ij(pos_i, particles[b].pos);
        return abs2(r_ia) < abs2(r_ib);
    });

    return count;
}

}
