target_sources(sph
    PRIVATE
        main.cpp
        bhtree.cpp
        exhaustive_search.cpp
        fluid_force.cpp
        gravity_force.cpp
        logger.cpp
        output.cpp
        pre_interaction.cpp
        simulation.cpp
        solver.cpp
        timestep.cpp
    )

add_subdirectory(disph)
add_subdirectory(gsph)
add_subdirectory(sample)