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

// int arraySize(Sphere array[])
// {
//     return sizeof(array) / sizeof(array[0]);
// }

inline bool intersect(const Ray &r, double &t, int &id, Sphere spheres[], int spheresLength)
{
    // arreglo de pares donde almacenamos la id y la distancia del rayo
    pair<int, double> spheresData[spheresLength];

    for (int i = 0; i < spheresLength; i++)
    {
        // asignación de valores al par
        spheresData[i].first = i;
        spheresData[i].second = spheres[i].intersect(r);
        if (spheresData[i].second > 0)
        {
            id = spheresData[i].first;
            t = spheresData[i].second;
        }
    }
    // ordenamiento de pares por la menor distancia
    for (int i = 0; i < spheresLength; i++)
    {
        if (t > spheresData[i].second && spheresData[i].second > 0.01)
        {
            // actualizamos valores de t e id solo si el valor almacenado en el par es menor al valor de t y que sea positivo
            id = spheresData[i].first;
            t = spheresData[i].second;
        }
    }
    if (t > 0)
    {
        return true;
    }
    return false;
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
// transforma un vector local a global
Vector makeGlobal(Vector &target, const Vector &n, const Vector &s, const Vector &t)
{
    Vector globalized(
        Vector(s.x, t.x, n.x).dot(target),
        Vector(s.y, t.y, n.y).dot(target),
        Vector(s.z, t.z, n.z).dot(target));
    return globalized;
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

int randint(int Min, int Max)
{
    return rand() % (Max + 1 - Min) + Min;
}

// Retorna un vector aleatorio en una esfera unitaria.
inline Vector random_in_sphere()
{
    double r1 = random_double();
    double r2 = random_double();

    double theta = acos(1.0 - (2.0 * r1));
    double phi = 2.0 * pi * r2;

    double x = cos(phi) * sin(theta);
    double y = sin(phi) * sin(theta);
    // auto z = cos(theta);
    double z = 1.0 - 2.0 * r1;

    return Vector(x, y, z);
}
// retorna un vector dado un angulo theta y phi
Vector createVec(double &theta, double &phi)
{
    double x = sin(theta) * cos(phi);
    double y = sin(theta) * sin(phi);
    double z = cos(theta);
    return Vector(x, y, z);
}
// retorna un vector dado un angulo theta, phi y radio de la esfera
Vector createVec(double &theta, double &phi, Sphere s)
{
    double radio = s.r;
    double x = radio * sin(theta) * cos(phi);
    double y = radio * sin(theta) * sin(phi);
    double z = radio * cos(theta);
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

// funciones para la importancia de luz según: área de la superficie
// obtener phi y theta para area
void areaPhiTheta(double &phi, double &theta)
{
    double r1 = random_double();
    double r2 = random_double();
    theta = acos(1.0 - 2.0 * r1);
    phi = 2.0 * pi * r2;
}
// prob Area
double areaProb(Point &x, Point &xprim, Sphere fuenteLuz)
{
    double r = fuenteLuz.r;
    Vector wi = x - xprim;
    double areaInv = 1.0 / (4.0 * pi * r * r);
    double distanceSquared = wi.dot(wi);
    Vector normal = (xprim - fuenteLuz.p).normalize();
    wi.normalize();
    double cosThetao = normal.dot(wi);
    return areaInv * (distanceSquared / cosThetao);
}
// evaluacionArea -> radiancia LE en la ecuacion de rendering
Color areaEval(Point &x, Vector &wi, int idLuz, Sphere fuenteLuz, Sphere spheres[], int spheresLength)
{
    Color Le;
    double t;
    double auxdist = sqrt(wi.dot(wi));
    wi.normalize();
    Ray ray(x, wi);
    int id = 0;
    if (intersect(ray, t, id, spheres, spheresLength))
    {
        double resta = auxdist - t;
        if (resta < 0.01)
        {
            Le = fuenteLuz.l;
        }
        else
        {
            Le = Color(0.0, 0.0, 0.0);
        }
    }
    else
    {
        Le = Color(0.0, 0.0, 0.0);
    }
    return Le;
}

// funciones para la importancia de luz según: angulo solido
// muestreo Solid Angle
void solidAnglePhiTheta(double &phi, double &theta, Point p, Sphere fuenteLuz, double &cosTethaMax)
{
    double r = fuenteLuz.r;
    double r1 = random_double();
    double r2 = random_double();

    Vector wc = (fuenteLuz.p - p);
    double pcp = sqrt(wc.dot(wc));
    double sinTethaMax = r / pcp;

    cosTethaMax = sqrt(1.0 - (sinTethaMax * sinTethaMax)); 
    theta = acos(1.0 - r1 + (r1 * cosTethaMax));
    phi = 2.0 * pi * r2;
}
// probabilidad de muestreo Solid Angle
double solidAngleProb(double cosTethaMax)
{
    return 1.0 / (2.0 * pi * (1.0 - cosTethaMax));
}
// evaluacion Solid Angle
Color solidAngleEval(Point &x, Vector &wi, Sphere fuenteLuz, Sphere spheres[], int spheresLength)
{
    Color Le;
    double t;
    // wi.normalize();
    Ray ray(x, wi);
    int id = 0;

    if (intersect(ray, t, id, spheres, spheresLength))
    {
       Sphere obj = spheres[id];
        if (obj.l.dot(obj.l) > 0.01)
        {
            Le = fuenteLuz.l;
        }
        else
        {
            Le = Color(0.0, 0.0, 0.0);
        }
    }
    else
    {
        Le = Color(0.0, 0.0, 0.0);
    }
    return Le;
}