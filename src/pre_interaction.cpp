#include <algorithm>

#include "parameters.hpp"
#include "pre_interaction.hpp"
#include "particle.hpp"
#include "openmp.hpp"

namespace sph
{

constexpr int neighbor_list_size = 10;
constexpr real kernel_ratio = 1.2;

void PreInteraction::initialize(std::shared_ptr<SPHParameters> param)
{
    use_time_dependent_av = param->av.use_time_dependent_av;
    use_balsara_switch = param->av.use_balsara_switch;
    gamma = param->physics.gamma;
    neighbor_number = param->physics.neighbor_number;
}

void PreInteraction::calculation(SPHParticle * particles, int num)
{
#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        std::shared_ptr<int[]> neighbor_list(new int(neighbor_number * neighbor_list_size));
        
        // neighbor search
        int const n_neighbor = exhaustive_search(p_i, p_i.sml * kernel_ratio, particles, num, neighbor_list, neighbor_number * neighbor_list_size);

        // density etc.
        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
        }
    }
}

int PreInteraction::exhaustive_search(
    SPHParticle & p_i,
    const real kernel_size,
    SPHParticle const * particles,
    const int num,
    std::shared_ptr<int[]> neighbor_list,
    const int list_size)
{
    const real kernel_size2 = kernel_size * kernel_size;
    const vec_t & pos_i = p_i.pos;
    int count = 0;
    for(int j = 0; j < num; ++j) {
        const vec_t r_ij = pos_i - particles[j].pos;
        const real r2 = abs2(r_ij);
        if(r2 < kernel_size2) {
            neighbor_list[count] = j;
            ++count;
        }
    }

    std::sort(neighbor_list.get(), neighbor_list.get() + count, [&](const int a, const int b) {
        const vec_t r_ia = pos_i - particles[a].pos;
        const vec_t r_ib = pos_i - particles[b].pos;
        return abs2(r_ia) < abs2(r_ib);
    });

    return count;
}

}
