#include "../src/vec.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <math.h>

double pi = 3.141592;

inline double random_double()
{
    return rand() / (RAND_MAX + 1.0);
}

inline Vector random_cosine_direction()
{
    auto r1 = random_double();
    auto r2 = random_double();
    auto z = sqrt(1 - r2);

    auto phi = 2 * pi * r1;
    auto x = cos(phi) * sqrt(r2);
    auto y = sin(phi) * sqrt(r2);

    return Vector(x, y, z);
}

int main()
{
    int N = 100;

    auto sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        auto v = random_cosine_direction();
        std::cout << v.x << " " << v.y << " " << v.z << '\n';
        sum += v.z * v.z * v.z / (v.z / pi);
    }

    std::cout << std::fixed << std::setprecision(12);
    std::cout << "Pi/2     = " << pi / 2 << '\n';
    std::cout << "Estimate = " << sum / N << '\n';
}