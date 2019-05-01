#pragma once

#include <memory>

namespace sph
{

class SPHParticle;
class KernelFunction;
class Distance;

#define ADD_MEMBER(type, name)\
public:\
    void set_##name(type v) { name = v; }\
    type & get_##name() { return name; }\
private:\
    type name

class Simulation {
    ADD_MEMBER(std::shared_ptr<SPHParticle[]>, particles);
    ADD_MEMBER(int, particle_num);
    ADD_MEMBER(real, time);
    ADD_MEMBER(real, dt);
    ADD_MEMBER(real, h_per_v_sig_max);
    ADD_MEMBER(std::shared_ptr<KernelFunction>, kernel);
    ADD_MEMBER(std::shared_ptr<Distance>, distance);

public:
    void update_time()
    {
        time += dt;
    }
};

#undef ADD_MEMBER

}