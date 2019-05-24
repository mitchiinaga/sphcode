#include <cassert>

#include <iostream>
#include <chrono>

#include "solver.hpp"
#include "parameters.hpp"
#include "particle.hpp"
#include "logger.hpp"
#include "exception.hpp"
#include "output.hpp"
#include "simulation.hpp"
#include "periodic.hpp"
#include "bhtree.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sph
{

Solver::Solver(int argc, char * argv[])
{
    std::cout << "--------------SPH simulation-------------\n\n";
    if(argc == 1) {
        std::cerr << "how to use\n" << std::endl;
        std::cerr << "sph <paramter.json>" << std::endl;
        std::exit(EXIT_FAILURE);
    } else {
        read_parameterfile(argv[1]);
    }

    Logger::open(m_output_dir);

#ifdef _OPENMP
    WRITE_LOG << "Open MP is valid.";
    int num_threads;
    if(argc == 3) {
        num_threads = std::atoi(argv[2]);
        omp_set_num_threads(num_threads);
    } else {
        num_threads = omp_get_max_threads();
    }
    WRITE_LOG << "the number of threads = " << num_threads << "\n";
#else
    WRITE_LOG << "OpenMP is invalid.\n";
#endif
    WRITE_LOG << "parameters";

    WRITE_LOG << "output directory     = " << m_output_dir;

    WRITE_LOG << "time";
    WRITE_LOG << "* start time         = " << m_param->time.start;
    WRITE_LOG << "* end time           = " << m_param->time.end;
    WRITE_LOG << "* output time        = " << m_param->time.output;
    WRITE_LOG << "* enerty output time = " << m_param->time.energy;

    WRITE_LOG << "CFL condition";
    WRITE_LOG << "* sound speed = " << m_param->cfl.sound;
    WRITE_LOG << "* force       = " << m_param->cfl.force;

    WRITE_LOG << "Artificial Viscosity";
    WRITE_LOG << "* alpha = " << m_param->av.alpha;
    if(m_param->av.use_balsara_switch) {
        WRITE_LOG << "* use Balsara switch";
    }
    if(m_param->av.use_time_dependent_av) {
        WRITE_LOG << "* use time dependent AV";
        WRITE_LOG << "* alpha max = " << m_param->av.alpha_max;
        WRITE_LOG << "* alpha min = " << m_param->av.alpha_min;
        WRITE_LOG << "* epsilon = " << m_param->av.epsilon;
    }

    if(m_param->ac.is_valid) {
        WRITE_LOG << "Artificial Conductivity";
        WRITE_LOG << "* alpha = " << m_param->ac.alpha;
    }

    WRITE_LOG << "Tree";
    WRITE_LOG << "* max tree level       = " << m_param->tree.max_level;
    WRITE_LOG << "* leaf particle number = " << m_param->tree.leaf_particle_num;

    WRITE_LOG << "Physics";
    WRITE_LOG << "* Neighbor number = " << m_param->physics.neighbor_number;
    WRITE_LOG << "* gamma           =" << m_param->physics.gamma;

    WRITE_LOG << "Kernel";
    if(m_param->kernel == KernelType::CUBIC_SPLINE) {
        WRITE_LOG << "* Cubic Spline";
    } else if(m_param->kernel == KernelType::WENDLAND) {
        WRITE_LOG << "* Wendland";
    } else {
        THROW_ERROR("kernel is unknown.");
    }

    if(m_param->iterative_sml) {
        WRITE_LOG << "iterative calculation for smoothing length is valid.";
    }

    if(m_param->periodic.is_valid) {
        WRITE_LOG << "Periodic boundary condition is valid.";
    }
    
    if(m_param->gravity.is_valid) {
        WRITE_LOG << "Gravity is valid.";
        WRITE_LOG << "G     = " << m_param->gravity.constant;
        WRITE_LOG << "theta = " << m_param->gravity.theta;
    }

    switch(m_sample) {
#define WRITE_SAMPLE(a, b) case a: WRITE_LOG << "Sample: " b " test"; break
        WRITE_SAMPLE(Sample::ShockTube, "shock tube");
        WRITE_SAMPLE(Sample::GreshoChanVortex, "Gresho-Chan vortex");
        WRITE_SAMPLE(Sample::HydroStatic, "Hydro static");
        WRITE_SAMPLE(Sample::KHI, "Kelvin-Helmholtz Instability");
        WRITE_SAMPLE(Sample::Evrard, "Evrard collapse");
#undef WRITE_SAMPLE
    }

    WRITE_LOG;

    m_output = std::make_shared<Output>();
}

void Solver::read_parameterfile(const char * filename)
{
    namespace pt = boost::property_tree;

    m_param = std::make_shared<SPHParameters>();

    pt::ptree input;

    std::string name_str = filename;
    if(name_str == "shock_tube") {
        pt::read_json("sample/shock_tube/shock_tube.json", input);
        m_sample = Sample::ShockTube;
        m_sample_parameters["N"] = input.get<int>("N", 100);
    } else if(name_str == "gresho_chan_vortex") {
        pt::read_json("sample/gresho_chan_vortex/gresho_chan_vortex.json", input);
        m_sample = Sample::GreshoChanVortex;
        m_sample_parameters["N"] = input.get<int>("N", 64);
    } else if(name_str == "hydrostatic") {
        pt::read_json("sample/hydrostatic/hydrostatic.json", input);
        m_sample = Sample::HydroStatic;
        m_sample_parameters["N"] = input.get<int>("N", 32);
    } else if(name_str == "khi") {
        pt::read_json("sample/khi/khi.json", input);
        m_sample = Sample::KHI;
        m_sample_parameters["N"] = input.get<int>("N", 128);
    } else if(name_str == "evrard") {
        pt::read_json("sample/evrard/evrard.json", input);
        m_sample = Sample::Evrard;
        m_sample_parameters["N"] = input.get<int>("N", 20);
    } else {
        pt::read_json(filename, input);
        m_sample = Sample::DoNotUse;
    }

    m_output_dir = input.get<std::string>("outputDirectory");

    // time
    m_param->time.start = input.get<real>("startTime", real(0));
    m_param->time.end   = input.get<real>("endTime");
    if(m_param->time.end < m_param->time.start) {
        THROW_ERROR("endTime < startTime");
    }
    m_param->time.output = input.get<real>("outputTime", (m_param->time.end - m_param->time.start) / 100);
    m_param->time.energy = input.get<real>("energyTime", m_param->time.output);

    // CFL
    m_param->cfl.sound = input.get<real>("cflSound", 0.3);
    m_param->cfl.force = input.get<real>("cflForce", 0.125);

    // Artificial Viscosity
    m_param->av.alpha = input.get<real>("avAlpha", 1.0);
    m_param->av.use_balsara_switch = input.get<bool>("useBalsaraSwitch", true);
    m_param->av.use_time_dependent_av = input.get<bool>("useTimeDependentAV", false);
    if(m_param->av.use_time_dependent_av) {
        m_param->av.alpha_max = input.get<real>("alphaMax", 2.0);
        m_param->av.alpha_min = input.get<real>("alphaMin", 0.1);
        if(m_param->av.alpha_max < m_param->av.alpha_min) {
            THROW_ERROR("alphaMax < alphaMin");
        }
        m_param->av.epsilon = input.get<real>("epsilonAV", 0.2);
    }

    // Artificial Conductivity
    m_param->ac.is_valid = input.get<bool>("useArtificialConductivity", false);
    if(m_param->ac.is_valid) {
        m_param->ac.alpha = input.get<real>("alphaAC", 1.0);
    }

    // Tree
    m_param->tree.max_level = input.get<int>("maxTreeLevel", 20);
    m_param->tree.leaf_particle_num = input.get<int>("leafParticleNumber", 4);

    // Physics
    m_param->physics.neighbor_number = input.get<int>("neighborNumber", 32);
    m_param->physics.gamma = input.get<real>("gamma");

    // Kernel
    std::string kernel_name = input.get<std::string>("kernel", "cubic_spline");
    if(kernel_name == "cubic_spline") {
        m_param->kernel = KernelType::CUBIC_SPLINE;
    } else if(kernel_name == "wendland") {
        m_param->kernel = KernelType::WENDLAND;
    } else {
        THROW_ERROR("kernel is unknown.");
    }

    // smoothing length
    m_param->iterative_sml = input.get<bool>("iterativeSmoothingLength", true);

    // periodic
    m_param->periodic.is_valid = input.get<bool>("periodic", false);
    if(m_param->periodic.is_valid) {
        {
            auto & range_max = input.get_child("rangeMax");
            if(range_max.size() != DIM) {
                THROW_ERROR("rangeMax != DIM");
            }
            int i = 0;
            for(auto & v : range_max) {
                m_param->periodic.range_max[i] = std::stod(v.second.data());
                ++i;
            }
        }

        {
            auto & range_min = input.get_child("rangeMin");
            if(range_min.size() != DIM) {
                THROW_ERROR("rangeMax != DIM");
            }
            int i = 0;
            for(auto & v : range_min) {
                m_param->periodic.range_min[i] = std::stod(v.second.data());
                ++i;
            }
        }
    }

    // gravity
    m_param->gravity.is_valid = input.get<bool>("useGravity", false);
    if(m_param->gravity.is_valid) {
        m_param->gravity.constant = input.get<real>("G", 1.0);
        m_param->gravity.theta = input.get<real>("theta", 0.5);
    }
}

void Solver::run()
{
    initialize();
    assert(m_sim->get_particles().size() == m_sim->get_particle_num());

    const real t_end = m_param->time.end;
    real t_out = m_param->time.output;
    real t_ene = m_param->time.energy;

    m_output->output_particle(m_sim);
    m_output->output_energy(m_sim);

    const auto start = std::chrono::system_clock::now();
    auto t_cout_i = start;
    int loop = 0;

    real t = m_sim->get_time();
    while(t < t_end) {
        integrate();
        const real dt = m_sim->get_dt();
        const int num = m_sim->get_particle_num();
        ++loop;

        m_sim->update_time();
        t = m_sim->get_time();
        
        // 1�b���Ƃɉ�ʏo�͂���
        const auto t_cout_f = std::chrono::system_clock::now();
        const real t_cout_s = std::chrono::duration_cast<std::chrono::seconds>(t_cout_f - t_cout_i).count();
        if(t_cout_s >= 1.0) {
            WRITE_LOG << "loop: " << loop << ", time: " << t << ", dt: " << dt << ", num: " << num;
            t_cout_i = std::chrono::system_clock::now();
        } else {
            WRITE_LOG_ONLY << "loop: " << loop << ", time: " << t << ", dt: " << dt << ", num: " << num;
        }

        if(t > t_out) {
            m_output->output_particle(m_sim);
            t_out += m_param->time.output;
        }

        if(t > t_ene) {
            m_output->output_energy(m_sim);
            t_ene += m_param->time.energy;
        }
    }
    const auto end = std::chrono::system_clock::now();
    const real calctime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    WRITE_LOG << "\ncalculation is finished";
    WRITE_LOG << "calclation time: " << calctime << " ms";
}

void Solver::initialize()
{
    m_sim = std::make_shared<Simulation>(m_param);

    make_initial_condition();

    m_timestep.initialize(m_param);
    m_pre.initialize(m_param);
    m_fforce.initialize(m_param);
    m_gforce.initialize(m_param);

    auto & p = m_sim->get_particles();
    const int num = m_sim->get_particle_num();
    const real gamma = m_param->physics.gamma;
    const real c_sound = gamma * (gamma - 1.0);

    assert(p.size() == num);
    const real alpha = m_param->av.alpha;
#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        p[i].alpha = alpha;
        p[i].balsara = 1.0;
        p[i].sound = std::sqrt(c_sound * p[i].ene);
    }

#ifndef EXHAUSTIVE_SEARCH
    auto tree = m_sim->get_tree();
    tree->resize(num);
    tree->make(p, num);
#endif

    m_pre.initial_smoothing(m_sim);
    m_pre.calculation(m_sim);
    m_fforce.calculation(m_sim);
    m_gforce.calculation(m_sim);
}

void Solver::integrate()
{
    m_timestep.calculation(m_sim);

    predict();
#ifndef EXHAUSTIVE_SEARCH
    m_sim->make_tree();
#endif
    m_pre.calculation(m_sim);
    m_fforce.calculation(m_sim);
    m_gforce.calculation(m_sim);
    correct();
}

void Solver::predict()
{
    auto & p = m_sim->get_particles();
    const int num = m_sim->get_particle_num();
    auto * periodic = m_sim->get_periodic().get();
    const real dt = m_sim->get_dt();
    const real gamma = m_param->physics.gamma;
    const real c_sound = gamma * (gamma - 1.0);

    assert(p.size() == num);

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        // k -> k+1/2
        p[i].vel_p = p[i].vel + p[i].acc * (0.5 * dt);
        p[i].ene_p = p[i].ene + p[i].dene * (0.5 * dt);

        // k -> k+1
        p[i].pos += p[i].vel_p * dt;
        p[i].vel += p[i].acc * dt;
        p[i].ene += p[i].dene * dt;
        p[i].sound = std::sqrt(c_sound * p[i].ene);

        periodic->apply(p[i].pos);
    }
}

void Solver::correct()
{
    auto & p = m_sim->get_particles();
    const int num = m_sim->get_particle_num();
    const real dt = m_sim->get_dt();
    const real gamma = m_param->physics.gamma;
    const real c_sound = gamma * (gamma - 1.0);

    assert(p.size() == num);

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        p[i].vel = p[i].vel_p + p[i].acc * (0.5 * dt);
        p[i].ene = p[i].ene_p + p[i].dene * (0.5 * dt);
        p[i].sound = std::sqrt(c_sound * p[i].ene);
    }
}

void Solver::make_initial_condition()
{
    switch(m_sample) {
#define MAKE_SAMPLE(a, b) case a: make_##b(); break
        MAKE_SAMPLE(Sample::ShockTube, shock_tube);
        MAKE_SAMPLE(Sample::GreshoChanVortex, gresho_chan_vortex);
        MAKE_SAMPLE(Sample::HydroStatic, hydrostatic);
        MAKE_SAMPLE(Sample::KHI, khi);
        MAKE_SAMPLE(Sample::Evrard, evrard);
        case Sample::DoNotUse:
            // make distribution
            break;
        default:
            THROW_ERROR("unknown sample type.");
#undef MAKE_SAMPLE
    }
}

}