#include "defines.hpp"
#include "fluid_force.hpp"
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

void FluidForce::initialize(std::shared_ptr<SPHParameters> param)
{
    m_neighbor_number = param->physics.neighbor_number;
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

#ifdef WITHOUT_GRADH
            acc -= dw_ij * (p_j.mass * (p_per_rho2_i + p_j.pres / sqr(p_j.dens) + pi_ij));
            dene += p_j.mass * (p_per_rho2_i + 0.5 * pi_ij) * inner_product(v_ij, dw_ij);
#else
            acc -= dw_i * (p_j.mass * (p_per_rho2_i * gradh_i + 0.5 * pi_ij)) + dw_j * (p_j.mass * (p_j.pres / sqr(p_j.dens) * p_j.gradh + 0.5 * pi_ij));
            dene += p_j.mass * p_per_rho2_i * gradh_i * inner_product(v_ij, dw_i) + 0.5 * p_j.mass * pi_ij * inner_product(v_ij, dw_ij);
#endif
        }

        p_i.acc = acc;
        p_i.dene = dene;
    }
}

real FluidForce::artificial_viscosity(const SPHParticle & p_i, const SPHParticle & p_j, const vec_t & r_ij)
{
    // Monaghan (1997)
    const auto v_ij = p_i.vel - p_j.vel;
    const real vr = inner_product(v_ij, r_ij);

    if(vr < 0) {
        const real alpha = 0.5 * (p_i.alpha + p_j.alpha);
        const real balsara = 0.5 * (p_i.balsara + p_j.balsara);
        const real w_ij = vr / std::abs(r_ij);
        const real v_sig = p_i.sound + p_j.sound - 3.0 * w_ij;
        const real rho_ij_inv = 2.0 / (p_i.dens + p_j.dens);
        
        const real pi_ij = -0.5 * balsara * alpha * v_sig * w_ij * rho_ij_inv;
        return pi_ij;
    } else {
        return 0;
    }
}

}