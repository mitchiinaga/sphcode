#include <iostream>

#include "core.hpp"
#include "exception.hpp"

int main(int argc, char *argv[])
{
    std::ios_base::sync_with_stdio(false);
    sph::exception_handler([&]() {
        sph::Core core(argc, argv);
        core.run();
    });
    return 0;
}