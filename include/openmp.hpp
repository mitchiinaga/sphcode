#pragma once

#include <memory>

#include "defines.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

class omp_real {
    int                     m_threads;
    std::unique_ptr<real[]> m_values;

public:
    omp_real(real const v = 0.0)
    {
#ifdef _OPENMP
        m_threads = omp_get_max_threads();
#else
        m_threads = 1;
#endif
        m_values = std::make_unique<real[]>(m_threads);

        for(int i = 0; i < m_threads; ++i) {
            m_values[i] = v;
        }
    }

    real & get()
    {
#ifdef _OPENMP
        return m_values[omp_get_thread_num()];
#else
        return m_values[0];
#endif
    }

    real min()
    {
        real v = m_values[0];
        for(int i = 1; i < m_threads; ++i) {
            if(v > m_values[i]) {
                v = m_values[i];
            }
        }
        return v;
    }

    real max()
    {
        real v = m_values[0];
        for(int i = 1; i < m_threads; ++i) {
            if(v < m_values[i]) {
                v = m_values[i];
            }
        }
        return v;
    }

    real sum()
    {
        real v = 0;
        for(int i = 0; i < m_threads; ++i) {
            v += m_values[i];
        }
        return v;
    }
};
