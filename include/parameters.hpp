#pragma once

#include "defines.hpp"

namespace sph
{

enum struct SPHType {
    SSPH,
    DISPH,
    GSPH,
};

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

    SPHType type;

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

    struct ArtificialConductivity {
        real alpha;
        bool is_valid;
    } ac;

    struct Tree {
        int max_level;
        int leaf_particle_num;
    } tree;

    struct Physics {
        int neighbor_number;
        real gamma;
    } physics;

    KernelType kernel;

    bool iterative_sml;

    struct Periodic {
        bool is_valid;
        real range_max[DIM];
        real range_min[DIM];
    } periodic;

    struct Gravity {
        bool is_valid;
        real constant;
        real theta;
    } gravity;

    struct GSPH {
        bool is_2nd_order;
    } gsph;
};

}