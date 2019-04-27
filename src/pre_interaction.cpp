#include "parameters.hpp"
#include "pre_interaction.hpp"
#include "particle.hpp"
#include "openmp.hpp"

namespace sph
{

void PreInteraction::initialize(std::shared_ptr<SPHParameters> param)
{
    use_time_dependent_av = param->av.use_time_dependent_av;
    use_balsara_switch = param->av.use_balsara_switch;
    gamma = param->physics.gamma;
    neighbor_number = param->physics.neighbor_number;
}

void PreInteraction::calculation(SPHParticle * particles, int num)
{
}

}
