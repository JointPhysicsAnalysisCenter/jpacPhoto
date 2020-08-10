#include "toy_monte_carlo.hpp"
#include <iostream>

using namespace jpacPhoto;

int main()
{
    toy_monte_carlo mc;

    mc.generate(10., 220000);
    std::cout << "we did it! \n";
    return 1;
};