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

    // determina si el rayo intersecta a esta esfera
    // double intersect(const Ray &ray) const
    // {
    //     // se asignan los valores de la formula general cuadratica "la del chicharronero"

    //     Vector oc = ray.o - p;
    //     Vector a = ray.d;
    //     double t = 0.0;
    //     double tol = 0.00001;
    //     double b = oc.dot(ray.d);
    //     double c = oc.dot(oc) - r * r;
    //     double discriminant = b * b - c;
    //     // // solo si el valor de discriminant es positivo hacemos la raiz cuadrada para ahorrar tiempo de computo
    //     // // de otra forma retorna -1
    //     if (discriminant < tol)
    //     {
    //         return 0.0;
    //     }
    //     else
    //     {
    //         discriminant = sqrt(discriminant);
    //     }
    //     // tminus
    //     t = -b - discriminant;

    //     // asegurarnos que t sea positivo
    //     if (t > tol)
    //         return t;
    //     else
    //         // tplus
    //         t = -b + discriminant;
    //     // asegurarnos que t sea positivo
    //     if (t > tol)
    //         return t;
    //     else
    //         return 0.0;
    // }

};
#endif