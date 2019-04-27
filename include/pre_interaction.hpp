#pragma once

#include "module.hpp"

namespace sph
{
struct SPHParameters;
class SPHParticle;

class PreInteraction : public Module {
    bool use_balsara_switch;
    bool use_time_dependent_av;
    real gamma;
    int neighbor_number;

    int exhaustive_search(
        SPHParticle & p_i,
        const real kernel_size,
        SPHParticle const * particles,
        const int num,
        std::shared_ptr<int[]> neighbor_list,
        const int list_size); // for debug

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(SPHParticle * particles, int num);
};
}
