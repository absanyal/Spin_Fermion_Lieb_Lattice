#ifndef RNGESUS_HPP
#define RNGESUS_HPP

// #include <stdint-gcc.h>
#include <iostream>
#include <cmath>

class xorshift64
{
public:
    uint64_t seed;
    bool seed_set = false;
    void set_seed(uint64_t);
    uint64_t get_seed();
    double random();
};

void xorshift64::set_seed(uint64_t tempseed)
{
    seed = tempseed;
    seed_set = true;
}

uint64_t xorshift64::get_seed()
{
    return seed;
}

double xorshift64::random()
{
    if (seed_set == false)
    {
        std::cout << "\nSeed was never set." << std::endl;
	return 0;
    }
    else
    {
        seed ^= seed << 13;
        seed ^= seed >> 17;
        seed ^= seed << 5;
        set_seed(seed);
        return seed * pow(2, -64);
    }
}

#endif
