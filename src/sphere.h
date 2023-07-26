#ifndef SPHEREH
#define SPHEREH

#include "vec.h"
#include "ray.h"
#include "material.h"
#pragma once

class Sphere
{
public:
    double r; // radio de la esfera
    Point p;  // posicion
    Material *m; // material
    // Color c;  // color
    // Color l;  // luz

    // Sphere(double r_, Point p_, Color c_, Color l_) : r(r_), p(p_), c(c_), l(l_) {}
    Sphere(double r_, Point p_, Material* m_) : r(r_), p(p_), m(m_) {}

    double intersect(const Ray &ray) const
    {
        // regresar distancia si hay intersección
        // regresar 0.0 si no hay interseccion
        Vector op = ray.o - p;
        Vector d = ray.d;
        double t;
        double tol = 0.00001;
        double b = op.dot(d);
        double ce = op.dot(op) - r * r;
        double disc = b * b - ce;
        if (disc < 0)
            return 0.0;
        else
            disc = sqrt(disc);
        t = -b - disc;
        if (t > tol)
            return t;
        else
            t = -b + disc;
        if (t > tol)
            return t;
        else
            return 0;
    }

    void get_UV( const Point &hitpoint, double &u, double &v ) const
    {
        // Vector d = x - p;
        // u = 0.5 + atan2(d.z, d.x) / (2 * M_PI);
        // v = 0.5 - asin(d.y / r) / M_PI;
        auto theta = acos(-hitpoint.y);
        auto phi = atan2(-hitpoint.z, hitpoint.x) + pi;

        u = phi / (2 * pi);
        v = theta / pi;
    }


};
#endif