#include "vec.h"
#include "ray.h"
#include "sphere.h"

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
using namespace std;

const double pi = 3.141592;


// limita el valor de x a [0,1]
inline double
clamp(const double x)
{
    if (x < 0.0)
        return 0.0;
    else if (x > 1.0)
        return 1.0;
    return x;
}

// convierte un valor de color en [0,1] a un entero en [0,255]
inline int toDisplayValue(const double x)
{
    return int(pow(clamp(x), 1.0 / 2.2) * 255 + .5);
}

Vector cross(Vector a, Vector b)
{
    return a % b;
}

void coordinateSystem(const Vector &n, Vector &s, Vector &t)
{
    if (std::abs(n.x) > std::abs(n.y))
    {
        float invLen = 1.0f / std::sqrt(n.x * n.x + n.z * n.z);
        t = Vector(n.z * invLen, 0.0f, -n.x * invLen);
    }
    else
    {
        float invLen = 1.0f / std::sqrt(n.y * n.y + n.z * n.z);
        t = Vector(0.0f, n.z * invLen, -n.y * invLen);
    }
    s = cross(t, n);
}

// Retorna real aleatorio en el rango  [0,1).
inline double random_double()
{
    return drand48();
}

// Retorna real aleatorio en el rango [min,max).
inline double random_double(double min, double max)
{
    return min + (max - min) * random_double();
}

// Retorna un vector aleatorio en una esfera unitaria.
inline Vector random_in_sphere()
{
    auto r1 = random_double();
    auto r2 = random_double();

    auto x = cos(2 * pi * r1) * 2 * sqrt(r2 * (1 - r2));
    auto y = sin(2 * pi * r1) * 2 * sqrt(r2 * (1 - r2));
    auto z = 1 - 2 * r2;

    return Vector(x, y, z);
}

// Retorna un vector aleatorio en un hemisferio unitario.
inline Vector random_in_hemisphere()
{
    auto r1 = random_double();
    auto r2 = random_double();

    auto theta = acos(1 - r1);
    auto phi = 2 * pi * r2;

    auto x = sin(theta) * cos(phi);
    auto y = sin(theta) * sin(phi);
    auto z = 1 - r1;

    return Vector(x, y, z);
}

// Retorna un vector aleatorio en coseno hemisferico unitario.
inline Vector random_cosine_hemisphere(double &theta)
{
    auto r1 = random_double();
    auto r2 = random_double();

    auto phi = 2 * pi * r1;

    auto z = sqrt(1 - r2);
    auto x = cos(phi) * sqrt(r2);
    auto y = sin(phi) * sqrt(r2);
    
    theta = acos(z);
    
    return Vector(x, y, z);
}