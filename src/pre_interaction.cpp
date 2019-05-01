#include <algorithm>

#include "parameters.hpp"
#include "pre_interaction.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "distance.hpp"
#include "openmp.hpp"
#include "kernel/kernel_function.hpp"
#include "exception.hpp"

namespace sph
{

constexpr real kernel_ratio = 1.2;

void PreInteraction::initialize(std::shared_ptr<SPHParameters> param)
{
    m_use_time_dependent_av = param->av.use_time_dependent_av;
    m_use_balsara_switch = param->av.use_balsara_switch;
    m_gamma = param->physics.gamma;
    m_neighbor_number = param->physics.neighbor_number;
}

void PreInteraction::calculation(std::shared_ptr<Simulation> sim)
{
    auto * particles = sim->get_particles().get();
    auto * distance = sim->get_distance().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        std::shared_ptr<int[]> neighbor_list(new int(m_neighbor_number * neighbor_list_size));
        
        // neighbor search
        int const n_neighbor = exhaustive_search(p_i, p_i.sml * kernel_ratio, particles, num, neighbor_list, m_neighbor_number * neighbor_list_size, distance);

        // smoothing length
        constexpr real A = DIM == 1 ? 2.0 :
                           DIM == 2 ? M_PI :
                                      4.0 * M_PI / 3.0;
        p_i.sml = std::pow(m_neighbor_number * p_i.mass / (p_i.dens * A), 1.0 / DIM);

        // density etc.
        real dens_i = 0.0;
        const vec_t & pos_i = p_i.pos;
        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = distance->calc_r_ij(pos_i, p_j.pos);
            const real r = abs(r_ij);

            if(r >= p_i.sml) {
                break;
            }

            dens_i += p_j.mass * kernel->w(r, p_i.sml);
        }

        p_i.dens = dens_i;
        p_i.pres = (m_gamma - 1.0) * dens_i * p_i.ene;
        p_i.sound = std::sqrt(m_gamma * p_i.pres / dens_i);
    }

    // signal velocity
    THROW_ERROR("signal velocity");
}

int PreInteraction::exhaustive_search(
    SPHParticle & p_i,
    const real kernel_size,
    SPHParticle const * particles,
    const int num,
    std::shared_ptr<int[]> neighbor_list,
    const int list_size,
    Distance const * distance)
{
    const real kernel_size2 = kernel_size * kernel_size;
    const vec_t & pos_i = p_i.pos;
    int count = 0;
    for(int j = 0; j < num; ++j) {
        const vec_t r_ij = distance->calc_r_ij(pos_i, particles[j].pos);
        const real r2 = abs2(r_ij);
        if(r2 < kernel_size2) {
            neighbor_list[count] = j;
            ++count;
        }
    }

    std::sort(neighbor_list.get(), neighbor_list.get() + count, [&](const int a, const int b) {
        const vec_t r_ia = distance->calc_r_ij(pos_i, particles[a].pos);
        const vec_t r_ib = distance->calc_r_ij(pos_i, particles[b].pos);
        return abs2(r_ia) < abs2(r_ib);
    });

    return count;
}

}
