#pragma once

#include <memory>

namespace sph
{

class SPHParticle;
class KernelFunction;
class Distance;

struct Simulation {
    std::shared_ptr<SPHParticle[]> particles;
    int particle_num;
    real time;
    real dt;
    real h_per_v_sig_max; // h / v_sig
    std::shared_ptr<KernelFunction> kernel;
    std::shared_ptr<Distance> distance;
};

}