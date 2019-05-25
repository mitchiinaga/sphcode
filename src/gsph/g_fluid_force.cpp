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
        const real p_per_rho2_i = p_i.pres / sqr(p_i.dens);
        const real h_i = p_i.sml;
        const real gradh_i = p_i.gradh;

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

            const vec_t dw_i = kernel->dw(r_ij, r, h_i);
            const vec_t dw_j = kernel->dw(r_ij, r, p_j.sml);
            const vec_t dw_ij = (dw_i + dw_j) * 0.5;
            const vec_t v_ij = v_i - p_j.vel;

            const real pi_ij = artificial_viscosity(p_i, p_j, r_ij);
            const real dene_ac = m_use_ac ? artificial_conductivity(p_i, p_j, r_ij, dw_ij) : 0.0;

#if 0
            acc -= dw_ij * (p_j.mass * (p_per_rho2_i + p_j.pres / sqr(p_j.dens) + pi_ij));
            dene += p_j.mass * (p_per_rho2_i + 0.5 * pi_ij) * inner_product(v_ij, dw_ij);
#else
            acc -= dw_i * (p_j.mass * (p_per_rho2_i * gradh_i + 0.5 * pi_ij)) + dw_j * (p_j.mass * (p_j.pres / sqr(p_j.dens) * p_j.gradh + 0.5 * pi_ij));
            dene += p_j.mass * p_per_rho2_i * gradh_i * inner_product(v_ij, dw_i) + 0.5 * p_j.mass * pi_ij * inner_product(v_ij, dw_ij) + dene_ac;
#endif
        }

        p_i.acc = acc;
        p_i.dene = dene;
    }
}

}
}