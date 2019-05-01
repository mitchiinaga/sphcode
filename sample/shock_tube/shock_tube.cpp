#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"

namespace sph
{

void Solver::make_shock_tube()
{
    if(DIM != 1) {
        THROW_ERROR("DIM != 1");
    }

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real dx_r = 0.5 / N;
    const real dx_l = dx_r * 0.25;

    const int num = N * 10;
    m_sim->set_particles(std::shared_ptr<SPHParticle>(new SPHParticle[num], std::default_delete<SPHParticle[]>()));
    m_sim->set_particle_num(num);

    auto * p = m_sim->get_particles().get();
    real x = -0.5 + dx_l * 0.5;
    real dx = dx_l;
    real dens = 1.0;
    real pres = 1.0;
    const real mass = real(N);
    const real gamma = m_param->physics.gamma;
    bool left = true;

    for(int i = 0; i < num; ++i) {
        auto & p_i = p[i];
        p_i.pos[0] = x;
        p_i.vel[0] = 0.0;
        p_i.dens = dens;
        p_i.pres = pres;
        p_i.mass = mass;
        p_i.ene = p_i.pres / ((gamma - 1.0) * p_i.dens);

        x += dx;
        if(x > 0.5 && left) {
            x = 0.5 + dx_r * 0.5;
            dx = dx_r;
            dens = 0.25;
            pres = 0.1795;
            left = false;
        }
    }
}

}