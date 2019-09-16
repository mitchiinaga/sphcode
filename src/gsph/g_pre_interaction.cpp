#include <algorithm>

#include "parameters.hpp"
#include "gsph/g_pre_interaction.hpp"
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
namespace gsph
{

void PreInteraction::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::PreInteraction::initialize(param);
    m_is_2nd_order = param->gsph.is_2nd_order;
}

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
    auto * tree = sim->get_tree().get();

    // for MUSCL
    auto & grad_d = sim->get_vector_array("density");
    auto & grad_p = sim->get_vector_array("pressure");
    vec_t * grad_v[DIM] = {
        sim->get_vector_array("velocity_0").data(),
#if DIM == 2
        sim->get_vector_array("velocity_1").data(),
#elif DIM == 3
        sim->get_vector_array("velocity_1").data(),
        sim->get_vector_array("velocity_2").data(),
#endif
    };

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

            if(i != j) {
                const real v_sig = p_i.sound + p_j.sound - 3.0 * inner_product(r_ij, p_i.vel - p_j.vel) / r;
                if(v_sig > v_sig_max) {
                    v_sig_max = v_sig;
                }
            }
        }

        p_i.dens = dens_i;
        p_i.pres = (m_gamma - 1.0) * dens_i * p_i.ene;
        p_i.neighbor = n_neighbor;
        p_i.v_sig = v_sig_max;

        // MUSCL法のための勾配計算
        if(!m_is_2nd_order) {
            continue;
        }

        vec_t dd, du; // dP = (gamma - 1) * (rho * du + drho * u)
        vec_t dv[DIM];
        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
            const real r = std::abs(r_ij);
            const vec_t dw_ij = kernel->dw(r_ij, r, p_i.sml);
            dd += dw_ij * p_j.mass;
            du += dw_ij * (p_j.mass * (p_j.ene - p_i.ene));
            for(int k = 0; k < DIM; ++k) {
                dv[k] += dw_ij * (p_j.mass * (p_j.vel[k] - p_i.vel[k]));
            }
        }
        grad_d[i] = dd;
        grad_p[i] = (dd * p_i.ene + du) * (m_gamma - 1.0);
        const real rho_inv = 1.0 / p_i.dens;
        for(int k = 0; k < DIM; ++k) {
            grad_v[k][i] = dv[k] * rho_inv;
        }
    }

#ifndef EXHAUSTIVE_SEARCH
    tree->set_kernel();
#endif
}

}
}
