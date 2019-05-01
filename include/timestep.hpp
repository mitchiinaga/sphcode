#pragma once

#include "module.hpp"

namespace sph
{

class TimeStep : public Module {
    real c_sound; // dt_s = c_sound * h / c
    real c_force; // dt_f = c_force * sqrt(h / a)
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};
}