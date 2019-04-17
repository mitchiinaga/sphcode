#include "sph.hpp"
#include "parameters.hpp"

namespace sph
{

SPH::SPH(std::shared_ptr<SPHParameters> & param)
{
}

void SPH::output_particles(const real time)
{

}

void SPH::output_energy(const real time)
{

}

void SPH::integrate(real * time)
{
    *time += 0.01;
}

int SPH::get_particle_num()
{
    return 0;
}

}