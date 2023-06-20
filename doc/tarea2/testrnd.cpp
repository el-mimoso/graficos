
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <random>

using namespace std;

int main(int argc, char const *argv[])
{
    unsigned int seed = 95;
    std::mt19937 gen(seed);

    FILE *f = fopen("RandomTests.txt", "w");

    // std::random_device rd;
    // std::mt19937 gen(rd());

    std::uniform_real_distribution<double> dis(0.0, 1.0);

    for (size_t i = 0; i < 10000; i++)
    {
        double a  = rand() / (RAND_MAX + 1.0);
        // double a = drand48();
        // double a = rand() / (32767 + 1.0);
        // double a  = rand();

        // double a = dis(gen);
        fprintf(f, "%g \n",a);
    }

    fclose(f);

    return 0;
}


