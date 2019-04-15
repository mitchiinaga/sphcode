#include <iostream>
#include <string>
#include <cstdlib>
#include <chrono>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "core.hpp"
#include "vector_type.hpp"
#include "logger.hpp"
#include "parameters.hpp"
#include "exception.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sph
{

Core::Core(int argc, char * argv[])
{
    std::cout << "--------------SPH simulation-------------\n\n";
    if(argc == 1) {
        std::cerr << "how to use\n" << std::endl;
        std::cerr << "sph <paramter.json>" << std::endl;
        std::exit(EXIT_FAILURE);
    } else {
        read_parameterfile(argv[1]);
    }

    Logger::open(output_dir);
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

    WRITE_LOG << "output directory     = " << output_dir;

    WRITE_LOG << "time";
    WRITE_LOG << "* start time         = " << param->time.start;
    WRITE_LOG << "* end time           = " << param->time.end;
    WRITE_LOG << "* output time        = " << param->time.output;
    WRITE_LOG << "* enerty output time = " << param->time.energy;

    WRITE_LOG << "CFL condition";
    WRITE_LOG << "* sound speed = " << param->cfl.sound;
    WRITE_LOG << "* force       = " << param->cfl.force;

    WRITE_LOG << "Artificial Viscosity";
    WRITE_LOG << "* alpha = " << param->av.alpha;
    if(param->av.use_balsala_switch) {
        WRITE_LOG << "* use Balsala switch";
    }
    if(param->av.use_time_dependent_av) {
        WRITE_LOG << "* use time dependent AV";
    }

    WRITE_LOG << "Tree";
    WRITE_LOG << "* Neighbor number      = " << param->tree.neighbor_number;
    WRITE_LOG << "* max tree level       = " << param->tree.max_level;
    WRITE_LOG << "* leaf particle number = " << param->tree.leaf_particle_num;

    WRITE_LOG << "Physics";
    WRITE_LOG << "* gamma =" << param->physics.gamma;

    WRITE_LOG;
}

void Core::read_parameterfile(const char * filename)
{
    namespace pt = boost::property_tree;

    param = std::make_shared<SPHParameters>();

    pt::ptree input;
    pt::read_json(filename, input);

    output_dir = input.get<std::string>("outputDirectory");

    // time
    param->time.start = input.get<real>("startTime", real(0));
    param->time.end   = input.get<real>("endTime");
    if(param->time.end < param->time.start) {
        THROW_ERROR("endTime < startTime");
    }
    param->time.output = input.get<real>("outputTime", (param->time.end - param->time.start) / 100);
    param->time.energy = input.get<real>("energyTime", param->time.output);

    // CFL
    param->cfl.sound = input.get<real>("cflSound", 0.3);
    param->cfl.force = input.get<real>("cflForce", 0.25);

    // Artificial Viscosity
    param->av.alpha = input.get<real>("avAlpha", 1.0);
    param->av.use_balsala_switch = input.get<bool>("useBalsalaSwitch", true);
    param->av.use_time_dependent_av = input.get<bool>("useTimeDependentAV", false);

    // Tree
    param->tree.neighbor_number = input.get<int>("neighborNumber", 32);
    param->tree.max_level = input.get<int>("maxTreeLevel", 20);
    param->tree.leaf_particle_num = input.get<int>("leafParticleNumber", 1);

    // Physics
    param->physics.gamma = input.get<real>("gamma");
}

void Core::run()
{
//    auto p = std::make_unique<SPH>(&param);
//    p->settings();
//    output_particle(p.get());
//    real t_out = outtime;
//    const auto start = std::chrono::system_clock::now();
//    real t = 0.0;
//    real t_b = 0.0;
//    int loop = 0;
//    auto t_cout_i = start;
//    while(t < endtime) {
//        p->integrate();
//        t = p->get_time();
//        ++loop;
//        
//        // 1ïbÇ≤Ç∆Ç…âÊñ èoóÕÇ∑ÇÈ
//        const auto t_cout_f = std::chrono::system_clock::now();
//        const real t_cout_s = std::chrono::duration_cast<std::chrono::seconds>(t_cout_f - t_cout_i).count();
//        if(t_cout_s >= 1.0) {
//            WRITE_LOG << "loop: " << loop << ", time: " << t << ", dt: " << t - t_b << ", num: " << p->size();
//            t_cout_i = std::chrono::system_clock::now();
//        } else {
//            WRITE_LOG_ONLY << "loop: " << loop << ", time: " << t << ", dt: " << t - t_b << ", num: " << p->size();
//        }
//
//        t_b = t;
//        if(t > t_out) {
//            const int next_count = output_particle(p.get());
//            t_out = outtime * next_count;
//        }
//    }
//    const auto end = std::chrono::system_clock::now();
//    const real calctime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    WRITE_LOG << "\ncalculation is finished";
//    WRITE_LOG << "calclation time: " << calctime << " ms";
}

}
