#pragma once

#include <memory>

namespace sph
{
struct SPHParameters;
class SPHParticle;

class TimeStep {
    real c_sound; // dt_s = c_sound * h / c
    real c_force; // dt_f = c_force * sqrt(h / a)
public:
    void initialize(std::shared_ptr<SPHParameters> param);
    real calculation(SPHParticle * particles, int num);
};
}