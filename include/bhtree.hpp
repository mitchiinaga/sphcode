#pragma once

#include <memory>
#include <vector>
#include <algorithm>

#include "vector_type.hpp"
#include "particle.hpp"

namespace sph
{
struct SPHParameters;

constexpr int NCHILD = DIM == 1 ? 2 : DIM == 2 ? 4 : 8;

class BHTree
{
    struct BHNode
    {
        SPHParticle * next;
        real mass;
        int num;
        BHNode * childs[NCHILD];
        vec_t pos;
        real edge;

        void clear() {
            next = nullptr;
            mass = 0.0;
            num = 0;
            std::fill(childs, childs + NCHILD, nullptr);
            pos = 0.0;
            edge = 0.0;
        }
    };

    int  m_max_level;
    int  m_leaf_particle_num;
    bool m_periodic;
    vec_t m_range_max;
    vec_t m_range_min;
    std::vector<BHNode> m_nodes;
public:
    void initialize(std::shared_ptr<SPHParameters> param);
    void assign(const int particle_num, const int tree_size = 5);
    void make(std::vector<SPHParticle> & particles, const int particle_num);
};

}