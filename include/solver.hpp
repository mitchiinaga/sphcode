#pragma once

#include <memory>
#include <unordered_map>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/any.hpp>

#include "defines.hpp"

namespace sph
{

struct SPHParameters;
class Simulation;
class Output;

class Module;

enum struct Sample {
    ShockTube,
    GreshoChanVortex,
    PairingInstability,
    HydroStatic,
    KHI,
    Evrard,
    DoNotUse,
};

class Solver {
    std::shared_ptr<SPHParameters>  m_param;
    std::shared_ptr<Output>         m_output;
    std::string                     m_output_dir;
    std::shared_ptr<Simulation>     m_sim;

    // modules
    std::shared_ptr<Module> m_timestep;
    std::shared_ptr<Module> m_pre;
    std::shared_ptr<Module> m_fforce;
    std::shared_ptr<Module> m_gforce;

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
    void make_pairing_instability();
    void make_hydrostatic();
    void make_khi();
    void make_evrard();

public:
    Solver(int argc, char * argv[]);
    void run();
};

}