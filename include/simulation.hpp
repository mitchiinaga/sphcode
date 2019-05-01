#pragma once

#include <memory>

namespace sph
{

struct SPHParameters;
class SPHParticle;
class KernelFunction;
class Distance;

#define ADD_MEMBER(type, name)\
public:\
    void set_##name(type v) { m_##name = v; }\
    type & get_##name() { return m_##name; }\
private:\
    type m_##name

class Simulation {
    ADD_MEMBER(std::shared_ptr<SPHParticle[]>, particles);
    ADD_MEMBER(int, particle_num);
    ADD_MEMBER(real, time);
    ADD_MEMBER(real, dt);
    ADD_MEMBER(real, h_per_v_sig_max);
    ADD_MEMBER(std::shared_ptr<KernelFunction>, kernel);
    ADD_MEMBER(std::shared_ptr<Distance>, distance);

public:
    Simulation(std::shared_ptr<SPHParameters> param);
    void update_time();
};

#undef ADD_MEMBER

}