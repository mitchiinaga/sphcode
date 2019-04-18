#pragma once

#include <fstream>
#include "defines.hpp"

namespace sph
{
class SPHParticle;

class Output {
    int m_count;
    std::ofstream m_out_energy;
public:
    Output(int count = 0);
    ~Output();
    void output_particle(const SPHParticle * particles, const int num, const real time);
    void output_energy(const SPHParticle * particles, const int num, const real time);
};

}
