#pragma once

#include "defines.hpp"

namespace sph
{
enum struct KernelType {
    CUBIC_SPLINE,
    WENDLAND,
    UNKNOWN,
};

struct SPHParameters {

    struct Time {
        real start;
        real end;
        real output;
        real energy;
    } time;

    struct CFL {
        real sound;
        real force;
    } cfl;

    struct ArtificialViscosity {
        real alpha;
        bool use_balsara_switch;
        bool use_time_dependent_av;
        real alpha_max;
        real alpha_min;
        real epsilon; // tau = h / (epsilon * c)
    } av;

    struct Tree {
        int max_level;
        int leaf_particle_num;
    } tree;

    struct Physics {
        int neighbor_number;
        real gamma;
    } physics;

    KernelType kernel;

    struct Periodic {
        bool is_valid;
        real range_max[DIM];
        real range_min[DIM];
    } periodic;
};

}