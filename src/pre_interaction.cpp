#include <algorithm>

#include "parameters.hpp"
#include "pre_interaction.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "periodic.hpp"
#include "openmp.hpp"
#include "kernel/kernel_function.hpp"
#include "exception.hpp"

namespace sph
{

void PreInteraction::initialize(std::shared_ptr<SPHParameters> param)
{
    m_use_time_dependent_av = param->av.use_time_dependent_av;
    m_use_balsara_switch = param->av.use_balsara_switch;
    m_gamma = param->physics.gamma;
    m_neighbor_number = param->physics.neighbor_number;
    m_kernel_ratio = 1.0;
}

void PreInteraction::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();

    omp_real h_per_v_sig(std::numeric_limits<real>::max());

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

        // guess smoothing length
        constexpr real A = DIM == 1 ? 2.0 :
                           DIM == 2 ? M_PI :
                                      4.0 * M_PI / 3.0;
        p_i.sml = std::pow(m_neighbor_number * p_i.mass / (p_i.dens * A), 1.0 / DIM);
        
        // neighbor search
        int const n_neighbor = exhaustive_search(p_i, p_i.sml * m_kernel_ratio, particles, num, neighbor_list, m_neighbor_number * neighbor_list_size, periodic);
        p_i.neighbor = n_neighbor;

        // smoothing length

        // density etc.
        real dens_i = 0.0;
        real dh_dens_i = 0.0;
        real v_sig_max = p_i.sound * 2.0;
        const vec_t & pos_i = p_i.pos;
        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
            const real r = abs(r_ij);

            if(r >= p_i.sml) {
                break;
            }

            dens_i += p_j.mass * kernel->w(r, p_i.sml);
            dh_dens_i += p_j.mass * kernel->dhw(r, p_i.sml);

            if(i != j) {
                const real v_sig = p_i.sound + p_j.sound - inner_product(pos_i - p_j.pos, p_i.vel - p_j.vel) / r;
                if(v_sig > v_sig_max) {
                    v_sig_max = v_sig;
                }
            }
        }

        p_i.dens = dens_i;
        p_i.pres = (m_gamma - 1.0) * dens_i * p_i.ene;
        p_i.gradh = 1.0 / (1.0 + p_i.sml / (DIM * dens_i) * dh_dens_i);

        const real h_per_v_sig_i = p_i.sml / v_sig_max;
        if(h_per_v_sig.get() > h_per_v_sig_i) {
            h_per_v_sig.get() = h_per_v_sig_i;
        }
    }

    sim->set_h_per_v_sig(h_per_v_sig.min());

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        p_i.sound = std::sqrt(m_gamma * p_i.pres / p_i.dens);
    }
}

int PreInteraction::exhaustive_search(
    SPHParticle & p_i,
    const real kernel_size,
    const std::vector<SPHParticle> & particles,
    const int num,
    std::vector<int> & neighbor_list,
    const int list_size,
    Periodic const * periodic)
{
    const real kernel_size2 = kernel_size * kernel_size;
    const vec_t & pos_i = p_i.pos;
    int count = 0;
    for(int j = 0; j < num; ++j) {
        const vec_t r_ij = periodic->calc_r_ij(pos_i, particles[j].pos);
        const real r2 = abs2(r_ij);
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

void PreInteraction::initial_smoothing(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

        // guess smoothing length
        constexpr real A = DIM == 1 ? 2.0 :
                           DIM == 2 ? M_PI :
                                      4.0 * M_PI / 3.0;
        p_i.sml = std::pow(m_neighbor_number * p_i.mass / (p_i.dens * A), 1.0 / DIM);
        
        // neighbor search
        int const n_neighbor = exhaustive_search(p_i, p_i.sml, particles, num, neighbor_list, m_neighbor_number * neighbor_list_size, periodic);
        p_i.neighbor = n_neighbor;

        // density
        real dens_i = 0.0;
        const vec_t & pos_i = p_i.pos;
        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
            const real r = abs(r_ij);

            if(r >= p_i.sml) {
                break;
            }

            dens_i += p_j.mass * kernel->w(r, p_i.sml);
        }

        p_i.dens = dens_i;
    }
}

}
