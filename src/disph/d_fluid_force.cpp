#include "defines.hpp"
#include "disph/d_fluid_force.hpp"
#include "particle.hpp"
#include "periodic.hpp"
#include "simulation.hpp"
#include "bhtree.hpp"
#include "kernel/kernel_function.hpp"

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph
{
namespace disph
{

void FluidForce::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::FluidForce::initialize(param);
    m_gamma = param->physics.gamma;
}

// Hopkins (2013)
// pressure-energy formulation
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
        const real gamma2_u_i = sqr(m_gamma - 1.0) * p_i.ene;
        const real gamma2_u_per_pres_i = gamma2_u_i / p_i.pres;
        const real m_u_inv = 1.0 / (p_i.mass * p_i.ene);
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
            const real f_ij = 1.0 - gradh_i / (p_j.mass * p_j.ene);
            const real f_ji = 1.0 - p_j.gradh * m_u_inv;
            const real u_per_pres_j = p_j.ene / p_j.pres;

            const real pi_ij = artificial_viscosity(p_i, p_j, r_ij);
            const real dene_ac = m_use_ac ? artificial_conductivity(p_i, p_j, r_ij, dw_ij) : 0.0;

            acc -= dw_i * (p_j.mass * (gamma2_u_per_pres_i * p_j.ene * f_ij + 0.5 * pi_ij)) + dw_j * (p_j.mass * (gamma2_u_i * u_per_pres_j * f_ji + 0.5 * pi_ij));
            dene += p_j.mass * gamma2_u_per_pres_i * p_j.ene * f_ij * inner_product(v_ij, dw_i) + 0.5 * p_j.mass * pi_ij * inner_product(v_ij, dw_ij) + dene_ac;
        }

        p_i.acc = acc;
        p_i.dene = dene;
    }
}

}
}