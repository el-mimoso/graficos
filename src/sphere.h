#include "vec.h"
#include "ray.h"
#pragma once

class Sphere
{
public:
    double r; // radio de la esfera
    Point p;  // posicion
    Color c;  // color
    Color l;  // luz

    Sphere(double r_, Point p_, Color c_, Color l_) : r(r_), p(p_), c(c_), l(l_) {}

    // determina si el rayo intersecta a esta esfera
    // -1 si no toca, t de lo contrario.
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
            return (-b - sqrt(discriminant)) / a;
        }
    }
};