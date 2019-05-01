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
    }

    WRITE_LOG << "Tree";
    WRITE_LOG << "* max tree level       = " << m_param->tree.max_level;
    WRITE_LOG << "* leaf particle number = " << m_param->tree.leaf_particle_num;

    WRITE_LOG << "Physics";
    WRITE_LOG << "* Neighbor number = " << m_param->physics.neighbor_number;
    WRITE_LOG << "* gamma           =" << m_param->physics.gamma;

    WRITE_LOG;

    m_output = std::make_shared<Output>();
}

void Solver::read_parameterfile(const char * filename)
{
    namespace pt = boost::property_tree;

    m_param = std::make_shared<SPHParameters>();

    pt::ptree input;
    pt::read_json(filename, input);

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
    m_param->cfl.force = input.get<real>("cflForce", 0.25);

    // Artificial Viscosity
    m_param->av.alpha = input.get<real>("avAlpha", 1.0);
    m_param->av.use_balsara_switch = input.get<bool>("useBalsaraSwitch", true);
    m_param->av.use_time_dependent_av = input.get<bool>("useTimeDependentAV", false);

    // Tree
    m_param->tree.max_level = input.get<int>("maxTreeLevel", 20);
    m_param->tree.leaf_particle_num = input.get<int>("leafParticleNumber", 1);

    // Physics
    m_param->physics.neighbor_number = input.get<int>("neighborNumber", 32);
    m_param->physics.gamma = input.get<real>("gamma");
}

void Solver::run()
{
    initialize();

    real t = m_param->time.start;
    real t_b = t;
    const real t_end = m_param->time.end;
    real t_out = m_param->time.output;
    real t_ene = m_param->time.energy;

//    m_output->output_particle(p->get_particles(), p->get_particle_num(), t);
//    m_output->output_energy(p->get_particles(), p->get_particle_num(), t);

    const auto start = std::chrono::system_clock::now();
    auto t_cout_i = start;
    int loop = 0;

    while(t < t_end) {
        integrate(&t);
        ++loop;
        
        // 1•b‚²‚Æ‚É‰æ–Êo—Í‚·‚é
        const auto t_cout_f = std::chrono::system_clock::now();
        const real t_cout_s = std::chrono::duration_cast<std::chrono::seconds>(t_cout_f - t_cout_i).count();
        if(t_cout_s >= 1.0) {
//            WRITE_LOG << "loop: " << loop << ", time: " << t << ", dt: " << t - t_b << ", num: " << p->get_particle_num();
            t_cout_i = std::chrono::system_clock::now();
        } else {
//            WRITE_LOG_ONLY << "loop: " << loop << ", time: " << t << ", dt: " << t - t_b << ", num: " << p->get_particle_num();
        }

        if(t > t_out) {
//            m_output->output_particle(p->get_particles(), p->get_particle_num(), t);
            t_out += m_param->time.output;
        }

        if(t > t_ene) {
//            m_output->output_energy(p->get_particles(), p->get_particle_num(), t);
            t_ene += m_param->time.energy;
        }

        t_b = t;
    }
    const auto end = std::chrono::system_clock::now();
    const real calctime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    WRITE_LOG << "\ncalculation is finished";
    WRITE_LOG << "calclation time: " << calctime << " ms";
}

void Solver::initialize()
{
    m_sim = std::make_shared<Simulation>();

    m_timestep.initialize(m_param);
    m_pre.initialize(m_param);
    // calc_tree();
    m_pre.calculation(m_sim);
    // calc_force();
}

void Solver::integrate(real * time)
{
    m_timestep.calculation(m_sim);
    real const dt = m_sim->dt;
    predict(dt);
    // calc_tree();
//    m_pre.calculation(m_particles.get(), m_particle_num);
    // calc_force();
    correct(dt);
    *time += dt;
}

void Solver::predict(const real dt)
{
    SPHParticle * p = m_sim->particles.get();
    const int num = m_sim->particle_num;
    assert(p);

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        // k -> k+1/2
        p[i].vel_i = p[i].vel + p[i].acc * (0.5 * dt);
        p[i].ene_i = p[i].ene + p[i].dene * (0.5 * dt);

        // k -> k+1
        p[i].pos += p[i].vel_i * dt;
        p[i].vel += p[i].acc * dt;
        p[i].ene += p[i].dene * dt;
    }
}

void Solver::correct(const real dt)
{
    SPHParticle * p = m_sim->particles.get();
    const int num = m_sim->particle_num;
    assert(p);

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        p[i].vel = p[i].vel_i + p[i].acc * (0.5 * dt);
        p[i].ene = p[i].ene_i + p[i].dene * (0.5 * dt);
    }
}

}