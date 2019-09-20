#pragma once

#include <vector>
#include "module.hpp"

namespace sph
{

class TimeStep : public Module {
protected:
    real c_sound; // dt_s = c_sound * h / c
    real c_force; // dt_f = c_force * sqrt(h / a)
public:
    virtual void initialize(std::shared_ptr<SPHParameters> param) override;
};

namespace fixed {
class Timestep : public sph::TimeStep {
public:
    void calculation(std::shared_ptr<Simulation> sim) override;
};
}

namespace indivisual {
class Timestep : public sph::TimeStep {
    std::vector<int> m_timeids;
    int m_current;
    int m_neighbor_number;
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};
}

}
