#include <algorithm>

#include "parameters.hpp"
#include "pre_interaction.hpp"
#include "simulation.hpp"
#include "periodic.hpp"
#include "openmp.hpp"
#include "kernel/kernel_function.hpp"
#include "exception.hpp"
#include "bhtree.hpp"

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph
{

void PreInteraction::initialize(std::shared_ptr<SPHParameters> param)
{
    m_use_time_dependent_av = param->av.use_time_dependent_av;
    if(m_use_time_dependent_av) {
        m_alpha_max = param->av.alpha_max;
        m_alpha_min = param->av.alpha_min;
        m_epsilon = param->av.epsilon;
    }
    m_use_balsara_switch = param->av.use_balsara_switch;
    m_gamma = param->physics.gamma;
    m_neighbor_number = param->physics.neighbor_number;
    m_iteration = param->iterative_sml;
    if(m_iteration) {
        m_kernel_ratio = 1.2;
    } else {
        m_kernel_ratio = 1.0;
    }
}

void PreInteraction::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();
    const real dt = sim->get_dt();
    auto * tree = sim->get_tree().get();

    omp_real h_per_v_sig(std::numeric_limits<real>::max());

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

        // guess smoothing length
        constexpr real A = DIM == 1 ? 2.0 :
                           DIM == 2 ? M_PI :
                                      4.0 * M_PI / 3.0;
        p_i.sml = std::pow(m_neighbor_number * p_i.mass / (p_i.dens * A), 1.0 / DIM) * m_kernel_ratio;
        
        // neighbor search
#ifdef EXHAUSTIVE_SEARCH
        const int n_neighbor_tmp = exhaustive_search(p_i, p_i.sml, particles, num, neighbor_list, m_neighbor_number * neighbor_list_size, periodic, false);
#else
        const int n_neighbor_tmp = tree->neighbor_search(p_i, neighbor_list, particles, false);
#endif
        // smoothing length
        if(m_iteration) {
            p_i.sml = newton_raphson(p_i, particles, neighbor_list, n_neighbor_tmp, periodic, kernel);
        }

        // density etc.
        real dens_i = 0.0;
        real dh_dens_i = 0.0;
        real v_sig_max = p_i.sound * 2.0;
        const vec_t & pos_i = p_i.pos;
        int n_neighbor = 0;
        for(int n = 0; n < n_neighbor_tmp; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
            const real r = std::abs(r_ij);

            if(r >= p_i.sml) {
                break;
            }

            ++n_neighbor;
            dens_i += p_j.mass * kernel->w(r, p_i.sml);
            dh_dens_i += p_j.mass * kernel->dhw(r, p_i.sml);

            if(i != j) {
                const real v_sig = p_i.sound + p_j.sound - 3.0 * inner_product(r_ij, p_i.vel - p_j.vel) / r;
                if(v_sig > v_sig_max) {
                    v_sig_max = v_sig;
                }
            }
        }

        p_i.dens = dens_i;
        p_i.pres = (m_gamma - 1.0) * dens_i * p_i.ene;
        p_i.gradh = 1.0 / (1.0 + p_i.sml / (DIM * dens_i) * dh_dens_i);
        p_i.neighbor = n_neighbor;

        const real h_per_v_sig_i = p_i.sml / v_sig_max;
        if(h_per_v_sig.get() > h_per_v_sig_i) {
            h_per_v_sig.get() = h_per_v_sig_i;
        }

        // Artificial viscosity
        if(m_use_balsara_switch && DIM != 1) {
#if DIM != 1
            // balsara switch
            real div_v = 0.0;
#if DIM == 2
            real rot_v = 0.0;
#else
            vec_t rot_v = 0.0;
#endif
            for(int n = 0; n < n_neighbor; ++n) {
                int const j = neighbor_list[n];
                auto & p_j = particles[j];
                const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
                const real r = std::abs(r_ij);
                const vec_t dw = kernel->dw(r_ij, r, p_i.sml);
                const vec_t v_ij = p_i.vel - p_j.vel;
                div_v -= p_j.mass * inner_product(v_ij, dw);
                rot_v += vector_product(v_ij, dw) * p_j.mass;
            }
            div_v /= p_i.dens;
            rot_v /= p_i.dens;
            p_i.balsara = std::abs(div_v) / (std::abs(div_v) + std::abs(rot_v) + 1e-4 * p_i.sound / p_i.sml);

            // time dependent alpha
            if(m_use_time_dependent_av) {
                const real tau_inv = m_epsilon * p_i.sound / p_i.sml;
                const real dalpha = (-(p_i.alpha - m_alpha_min) * tau_inv + std::max(-div_v, 0.0) * (m_alpha_max - p_i.alpha)) * dt;
                p_i.alpha += dalpha;
            }
#endif
        } else if(m_use_time_dependent_av) {
            real div_v = 0.0;
            for(int n = 0; n < n_neighbor; ++n) {
                int const j = neighbor_list[n];
                auto & p_j = particles[j];
                const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
                const real r = std::abs(r_ij);
                const vec_t dw = kernel->dw(r_ij, r, p_i.sml);
                const vec_t v_ij = p_i.vel - p_j.vel;
                div_v -= p_j.mass * inner_product(v_ij, dw);
            }
            div_v /= p_i.dens;
            const real tau_inv = m_epsilon * p_i.sound / p_i.sml;
            const real dalpha = (-(p_i.alpha - m_alpha_min) * tau_inv + std::max(-div_v, 0.0) * (m_alpha_max - p_i.alpha)) * dt;
            p_i.alpha += dalpha;
        }
    }

    sim->set_h_per_v_sig(h_per_v_sig.min());

#ifndef EXHAUSTIVE_SEARCH
    tree->set_kernel();
#endif
}

void PreInteraction::initial_smoothing(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();
    auto * tree = sim->get_tree().get();

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        const vec_t & pos_i = p_i.pos;
        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

        // guess smoothing length
        constexpr real A = DIM == 1 ? 2.0 :
                           DIM == 2 ? M_PI :
                                      4.0 * M_PI / 3.0;
        p_i.sml = std::pow(m_neighbor_number * p_i.mass / (p_i.dens * A), 1.0 / DIM);
        
        // neighbor search
#ifdef EXHAUSTIVE_SEARCH
        int const n_neighbor = exhaustive_search(p_i, p_i.sml, particles, num, neighbor_list, m_neighbor_number * neighbor_list_size, periodic, false);
#else
        int const n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, false);
#endif

        // density
        real dens_i = 0.0;
        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
            const real r = std::abs(r_ij);

            if(r >= p_i.sml) {
                break;
            }

            dens_i += p_j.mass * kernel->w(r, p_i.sml);
        }

        p_i.dens = dens_i;
    }
}

inline real powh_(const real h) {
#if DIM == 1
    return 1;
#elif DIM == 2
    return h;
#elif DIM == 3
    return h * h;
#endif
}

real PreInteraction::newton_raphson(
    const SPHParticle & p_i,
    const std::vector<SPHParticle> & particles,
    const std::vector<int> & neighbor_list,
    const int n_neighbor,
    const Periodic * periodic,
    const KernelFunction * kernel
)
{
    real h_i = p_i.sml / m_kernel_ratio;
    constexpr real A = DIM == 1 ? 2.0 :
                       DIM == 2 ? M_PI :
                                  4.0 * M_PI / 3.0;
    const real b = p_i.mass * m_neighbor_number / A;

    // f = rho h^d - b
    // f' = drho/dh h^d + d rho h^{d-1}

    constexpr real epsilon = 1e-5;
    constexpr int max_iter = 10;
    const auto & r_i = p_i.pos;
    for(int i = 0; i < max_iter; ++i) {
        const real h_b = h_i;

        real dens = 0.0;
        real ddens = 0.0;
        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
            const real r = std::abs(r_ij);

            if(r >= h_i) {
                break;
            }

            dens += p_j.mass * kernel->w(r, h_i);
            ddens += p_j.mass * kernel->dhw(r, h_i);
        }

        const real f = dens * powh(h_i) - b;
        const real df = ddens * powh(h_i) + DIM * dens * powh_(h_i);

        h_i -= f / df;

        if(std::abs(h_i - h_b) < (h_i + h_b) * epsilon) {
            return h_i;
        }
    }

    return p_i.sml / m_kernel_ratio;
}

}
