#pragma once

#include <memory>

#include "particle.hpp"

namespace sph
{

struct Simulation {
    std::shared_ptr<SPHParticle[]> particles;
    int particle_num;
    real time;
    real dt;
};

}