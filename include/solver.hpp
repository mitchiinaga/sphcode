#pragma once

#include <memory>
#include <unordered_map>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/any.hpp>

#include "defines.hpp"
#include "timestep.hpp"
#include "pre_interaction.hpp"
#include "fluid_force.hpp"

namespace sph
{

struct SPHParameters;
class Simulation;
class Output;

enum struct Sample {
    ShockTube,
    GreshoChanVortex,
    HydroStatic,
    KHI,
    DoNotUse,
};

class Solver {
    std::shared_ptr<SPHParameters> m_param;
    std::shared_ptr<Output>        m_output;
    std::string                    m_output_dir;
    std::shared_ptr<Simulation>    m_sim;

    // modules
    TimeStep       m_timestep;
    PreInteraction m_pre;
    FluidForce     m_fforce;

    void read_parameterfile(const char * filename);
    void make_initial_condition();
    void initialize();
    void predict();
    void correct();
    void integrate();

    // for sample
    Sample                                      m_sample;
    std::unordered_map<std::string, boost::any> m_sample_parameters;

    void make_shock_tube();
    void make_gresho_chan_vortex();
    void make_hydrostatic();
    void make_khi();

public:
    Solver(int argc, char * argv[]);
    void run();
};

}