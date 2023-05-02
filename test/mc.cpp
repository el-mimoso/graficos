#include "../src/vec.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>

double pi = 3.141592;

inline double pdf(const Vector &p)
{
    return 1 / (4 * pi);
}

inline double random_double()
{
    return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max)
{
    return min + (max - min) * random_double();
}

Vector random_in_unit_sphere()
{
    while (true)
    {
        Vector p = Vector(random_double(-1,1), random_double(-1,1), random_double(-1,1));
        if (p.x*p.x + p.y*p.y + p.z*p.z >= 1)
            continue;
        return p;
    }
}

Vector random_unit_vector()
{
    return (random_in_unit_sphere()).normalize();
}


int main()
{
    int N = 100;
    auto sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        Vector d = random_unit_vector();
        std::cout << d.x << " " << d.y << " " << d.z << '\n';
        auto cosine_squared = d.z * d.z;
        sum += cosine_squared / pdf(d);
    }
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "I = " << sum / N << '\n';
}