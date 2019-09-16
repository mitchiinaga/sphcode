#include "defines.hpp"
#include "particle.hpp"
#include "periodic.hpp"
#include "simulation.hpp"
#include "bhtree.hpp"
#include "kernel/kernel_function.hpp"
#include "gsph/g_fluid_force.hpp"

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph
{
namespace gsph
{

void FluidForce::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::FluidForce::initialize(param);
    m_is_2nd_order = param->gsph.is_2nd_order;
    m_gamma = param->physics.gamma;

    hll_solver();
}

// van Leer (1979) limiter
inline real limiter(const real dq1, const real dq2)
{
    const real dq1dq2 = dq1 * dq2;
    if(dq1dq2 <= 0) {
        return 0.0;
    } else {
        return 2.0 * dq1dq2 / (dq1 + dq2);
    }
}

// Cha & Whitworth (2003)
void FluidForce::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();
    auto * tree = sim->get_tree().get();
    const real dt = sim->get_dt();
    auto timeid = sim->get_timeid();

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

#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        if(!(p_i.timeid & timeid)) {
            continue;
        }
        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);
        
        // neighbor search
#ifdef EXHAUSTIVE_SEARCH
        int const n_neighbor = exhaustive_search(p_i, p_i.sml, particles, num, neighbor_list, m_neighbor_number * neighbor_list_size, periodic, true);
#else
        int const n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, true);
#endif

        // fluid force
        const vec_t & r_i = p_i.pos;
        const vec_t & v_i = p_i.vel;
        const real h_i = p_i.sml;
        const real rho2_inv_i = 1.0 / sqr(p_i.dens);

        vec_t acc(0.0);
        real dene = 0.0;

        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
            const real r = std::abs(r_ij);

            if(r >= std::max(h_i, p_j.sml) || r == 0.0) {
                continue;
            }

            const real r_inv = 1.0 / r;
            const vec_t e_ij = r_ij * r_inv;
            const real ve_i = inner_product(v_i, e_ij);
            const real ve_j = inner_product(p_j.vel, e_ij);
            real vstar, pstar;

            if(m_is_2nd_order) {
                // Murante et al. (2011)

                real right[4], left[4];
                const real delta_i = 0.5 * (1.0 - p_i.sound * dt * r_inv);
                const real delta_j = 0.5 * (1.0 - p_j.sound * dt * r_inv);

                // velocity
                const real dv_ij = ve_i - ve_j;
                vec_t dv_i, dv_j;
                for(int k = 0; k < DIM; ++k) {
                    dv_i[k] = inner_product(grad_v[k][i], e_ij);
                    dv_j[k] = inner_product(grad_v[k][j], e_ij);
                }
                const real dve_i = inner_product(dv_i, e_ij) * r;
                const real dve_j = inner_product(dv_j, e_ij) * r;
                right[0] = ve_i - limiter(dv_ij, dve_i) * delta_i;
                left[0] = ve_j + limiter(dv_ij, dve_j) * delta_j;

                // density
                const real dd_ij = p_i.dens - p_j.dens;
                const real dd_i = inner_product(grad_d[i], e_ij) * r;
                const real dd_j = inner_product(grad_d[j], e_ij) * r;
                right[1] = p_i.dens - limiter(dd_ij, dd_i) * delta_i;
                left[1] = p_j.dens + limiter(dd_ij, dd_j) * delta_j;

                // pressure
                const real dp_ij = p_i.pres - p_j.pres;
                const real dp_i = inner_product(grad_p[i], e_ij) * r;
                const real dp_j = inner_product(grad_p[j], e_ij) * r;
                right[2] = p_i.pres - limiter(dp_ij, dp_i) * delta_i;
                left[2] = p_j.pres + limiter(dp_ij, dp_j) * delta_j;

                // sound speed
                right[3] = std::sqrt(m_gamma * right[2] / right[1]);
                left[3] = std::sqrt(m_gamma * left[2] / left[1]);

                m_solver(left, right, pstar, vstar);
            } else {
                const real right[4] = {
                    ve_i,
                    p_i.dens,
                    p_i.pres,
                    p_i.sound,
                };
                const real left[4] = {
                    ve_j,
                    p_j.dens,
                    p_j.pres,
                    p_j.sound,
                };

                m_solver(left, right, pstar, vstar);
            }

            const vec_t dw_i = kernel->dw(r_ij, r, h_i);
            const vec_t dw_j = kernel->dw(r_ij, r, p_j.sml);
            const vec_t v_ij = e_ij * vstar;
            const real rho2_inv_j = 1.0 / sqr(p_j.dens);
            const vec_t f = dw_i * (p_j.mass * pstar * rho2_inv_i) + dw_j * (p_j.mass * pstar * rho2_inv_j);

            acc -= f;
            dene -= inner_product(f, v_ij - v_i);
        }

        p_i.acc = acc;
        p_i.dene = dene;
    }
}

void FluidForce::hll_solver()
{
    m_solver = [&](const real left[], const real right[], real & pstar, real & vstar) {
        const real u_l   = left[0];
        const real rho_l = left[1];
        const real p_l   = left[2];
        const real c_l   = left[3];

        const real u_r   = right[0];
        const real rho_r = right[1];
        const real p_r   = right[2];
        const real c_r   = right[3];

        const real roe_l = std::sqrt(rho_l);
        const real roe_r = std::sqrt(rho_r);
        const real roe_inv = 1.0 / (roe_l + roe_r);

        const real u_t = (roe_l * u_l + roe_r * u_r) * roe_inv;
        const real c_t = (roe_l * c_l + roe_r * c_r) * roe_inv;
        const real s_l = std::min(u_l - c_l, u_t - c_t);
        const real s_r = std::max(u_r + c_r, u_t + c_t);

        const real c1 = rho_l * (s_l - u_l);
        const real c2 = rho_r * (s_r - u_r);
        const real c3 = 1.0 / (c1 - c2);
        const real c4 = p_l - u_l * c1;
        const real c5 = p_r - u_r * c2;
        
        vstar = (c5 - c4) * c3;
        pstar = (c1 * c5 - c2 * c4) * c3;
    };
}

}
}