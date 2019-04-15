#pragma once

#include "defines.hpp"

namespace sph
{

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
        bool use_balsala_switch;
        bool use_time_dependent_av;
    } av;

    struct Tree {
        int neighbor_number;
        int max_level;
        int leaf_particle_num;
    } tree;

    struct Physics {
        real gamma;
    } physics;
};

}