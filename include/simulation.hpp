#pragma once

#include <memory>
#include <vector>
#include <unordered_map>

#include "particle.hpp"

namespace sph
{

struct SPHParameters;
class KernelFunction;
class Periodic;
class BHTree;

#define ADD_MEMBER(type, name)\
public:\
    void set_##name(type v) { m_##name = v; }\
    type & get_##name() { return m_##name; }\
private:\
    type m_##name

class Simulation {
    ADD_MEMBER(std::vector<SPHParticle>, particles);
    ADD_MEMBER(int, particle_num);
    ADD_MEMBER(real, time);
    ADD_MEMBER(real, dt);
    ADD_MEMBER(real, h_per_v_sig);
    ADD_MEMBER(std::shared_ptr<KernelFunction>, kernel);
    ADD_MEMBER(std::shared_ptr<Periodic>, periodic);
    ADD_MEMBER(std::shared_ptr<BHTree>, tree);
    ADD_MEMBER(int, timeid);

    std::unordered_map<std::string, std::vector<real>> additional_scalar_array;
    std::unordered_map<std::string, std::vector<vec_t>> additional_vector_array;

public:
    Simulation(std::shared_ptr<SPHParameters> param);
    void update_time();
    void make_tree();
    void add_scalar_array(const std::vector<std::string> & names);
    void add_vector_array(const std::vector<std::string> & names);
    std::vector<real> & get_scalar_array(const std::string & name);
    std::vector<vec_t> & get_vector_array(const std::string & name);
};

#undef ADD_MEMBER

}