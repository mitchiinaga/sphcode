#include <cassert>

#include "sph.hpp"
#include "parameters.hpp"
#include "particle.hpp"

namespace sph
{

SPH::SPH(std::shared_ptr<SPHParameters> param)
{
    m_timestep.initialize(param);
}

void SPH::initialize()
{
    // calc_tree();
    // pre_interaction();
    // calc_force();
}

void SPH::integrate(real * time)
{
    real const dt = m_timestep.calculation(m_particles.get(), m_particle_num);
    predict(dt);
    // calc_tree();
    // pre_interaction();
    // calc_force();
    correct(dt);
    *time += dt;
}

int SPH::get_particle_num()
{
    return m_particle_num;
}

const SPHParticle * SPH::get_particles()
{
    const SPHParticle * p = m_particles.get();
    assert(p);
    return p;
}

void SPH::predict(const real dt)
{
    SPHParticle * p = m_particles.get();
    assert(p);

#pragma omp parallel for
    for(int i = 0; i < m_particle_num; ++i) {
        // k -> k+1/2
        p[i].vel_i = p[i].vel + p[i].acc * (0.5 * dt);
        p[i].ene_i = p[i].ene + p[i].dene * (0.5 * dt);

        // k -> k+1
        p[i].pos += p[i].vel_i * dt;
        p[i].vel += p[i].acc * dt;
        p[i].ene += p[i].dene * dt;
    }
}

void SPH::correct(const real dt)
{
    SPHParticle * p = m_particles.get();
    assert(p);

#pragma omp parallel for
    for(int i = 0; i < m_particle_num; ++i) {
        p[i].vel = p[i].vel_i + p[i].acc * (0.5 * dt);
        p[i].ene = p[i].ene_i + p[i].dene * (0.5 * dt);
    }
}

}