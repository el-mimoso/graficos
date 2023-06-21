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

    // determina si el rayo intersecta a esta esfera
    double intersect(const Ray &ray) const
    {
        // se asignan los valores de la formula general cuadratica "la del chicharronero"

        Vector oc = ray.o - p;
        double a = ray.d.dot(ray.d);
        double b = oc.dot(ray.d);
        double c = oc.dot(oc) - r * r;
        double discriminant = b * b - a * c;
        // solo si el valor de discriminant es positivo hacemos la raiz cuadrada para ahorrar tiempo de computo
        // de otra forma retorna -1
        if (discriminant < 0)
        {
            return 0.0;
        }
        else
        {
            double tplus = (-b + sqrt(discriminant)) / a;
            double tminus = (-b - sqrt(discriminant)) / a;
            double t;

            // ambos positivos
            if (tminus > 0 && tplus > 0)
            {
                t = std::min(tminus, tplus);
            }
            // tminus positivo, tplus negativo
            else if (tminus > 0 && tplus < 0)
            {
                t = tminus;
            }
            // tminus negativo, tplus positivo
            else if (tminus < 0 && tplus > 0)
            {
                t = tplus;
            }
            else
            {
                t = 0;
            }
            return t;
        }
    }
};
#endif