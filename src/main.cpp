#include <iostream>

#include "solver.hpp"
#include "exception.hpp"

int main(int argc, char *argv[])
{
    std::ios_base::sync_with_stdio(false);
    sph::exception_handler([&]() {
        sph::Solver solver(argc, argv);
        solver.run();
    });
    return 0;
}