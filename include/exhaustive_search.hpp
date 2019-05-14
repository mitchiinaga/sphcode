#pragma once

#include <vector>

#include "defines.hpp"

namespace sph
{

struct SPHParticle;
class Periodic;

int exhaustive_search(
    SPHParticle & p_i,
    const real kernel_size,
    const std::vector<SPHParticle> & particles,
    const int num,
    std::vector<int> & neighbor_list,
    const int list_size,
    Periodic const * periodic,
    const bool is_ij);

}
