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
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(SPHParticle * particles, int num);
};
}
