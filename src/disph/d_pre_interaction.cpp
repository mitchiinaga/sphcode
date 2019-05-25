#include <algorithm>

#include "disph/d_pre_interaction.hpp"
#include "parameters.hpp"
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
namespace disph
{

void PreInteraction::calculation(std::shared_ptr<Simulation> sim)
{
    if(m_first) {
        initial_smoothing(sim);
        m_first = false;
    }

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
        real pres_i = 0.0;
        real dh_pres_i = 0.0;
        real n_i = 0.0;
        real dh_n_i = 0.0;
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
            const real w_ij = kernel->w(r, p_i.sml);
            const real dhw_ij = kernel->dhw(r, p_i.sml);
            dens_i += p_j.mass * w_ij;
            n_i += w_ij;
            pres_i += p_j.mass * p_j.ene * w_ij;
            dh_pres_i += p_j.mass * p_j.ene * dhw_ij;
            dh_n_i += dhw_ij;

            if(i != j) {
                const real v_sig = p_i.sound + p_j.sound - 3.0 * inner_product(r_ij, p_i.vel - p_j.vel) / r;
                if(v_sig > v_sig_max) {
                    v_sig_max = v_sig;
                }
            }
        }

        p_i.dens = dens_i;
        p_i.pres = (m_gamma - 1.0) * pres_i;
        // f_ij = 1 - p_i.gradh / (p_j.mass * p_j.ene)
        p_i.gradh = p_i.sml / (DIM * n_i) * dh_pres_i / (1.0 + p_i.sml / (DIM * n_i) * dh_n_i);
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
                div_v -= p_j.mass * p_j.ene * inner_product(v_ij, dw);
                rot_v += vector_product(v_ij, dw) * (p_j.mass * p_j.ene);
            }
            const real p_inv = (m_gamma - 1.0) / p_i.pres;
            div_v *= p_inv;
            rot_v *= p_inv;
            p_i.balsara = std::abs(div_v) / (std::abs(div_v) + std::abs(rot_v) + 1e-4 * p_i.sound / p_i.sml);

            // time dependent alpha
            if(m_use_time_dependent_av) {
                const real tau_inv = m_epsilon * p_i.sound / p_i.sml;
                const real dalpha = (-(p_i.alpha - m_alpha_min) * tau_inv + std::max(-div_v, (real)0.0) * (m_alpha_max - p_i.alpha)) * dt;
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
                div_v -= p_j.mass * p_j.ene * inner_product(v_ij, dw);
            }
            const real p_inv = (m_gamma - 1.0) / p_i.pres;
            div_v *= p_inv;
            const real tau_inv = m_epsilon * p_i.sound / p_i.sml;
            const real dalpha = (-(p_i.alpha - m_alpha_min) * tau_inv + std::max(-div_v, (real)0.0) * (m_alpha_max - p_i.alpha)) * dt;
            p_i.alpha += dalpha;
        }
    }

    sim->set_h_per_v_sig(h_per_v_sig.min());

#ifndef EXHAUSTIVE_SEARCH
    tree->set_kernel();
#endif
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
    const real b = m_neighbor_number / A;

    // f = n h^d - b
    // f' = dn/dh h^d + d n h^{d-1}

    constexpr real epsilon = 1e-4;
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

            dens += kernel->w(r, h_i);
            ddens += kernel->dhw(r, h_i);
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
}
