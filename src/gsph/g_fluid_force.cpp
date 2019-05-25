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

void FluidForce::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();
    auto * tree = sim->get_tree().get();

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
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

            const vec_t e_ij = r_ij / r;
            const real ve_i = inner_product(v_i, e_ij);
            const real ve_j = inner_product(p_j.vel, e_ij);

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
            real vstar, pstar;
            m_solver(left, right, pstar, vstar);

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