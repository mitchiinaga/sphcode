#pragma once

#include <fstream>
#include <memory>

#include "defines.hpp"

namespace sph
{
class SPHParticle;
class Simulation;

class Output {
    int m_count;
    std::ofstream m_out_energy;
public:
    Output(int count = 0);
    ~Output();
    void output_particle(std::shared_ptr<Simulation> sim);
    void output_energy(std::shared_ptr<Simulation> sim);
};

}
