#ifndef MATERIALH
#define MATERIALH

#pragma once
#include "vec.h"
#include "ray.h"
#include "utils.h"
#include "sphere.h"

class Material
{
public:
    // probabilidad
    virtual double probability() { return 0.0; };
    // muestreo de direcciones
    virtual Vector sampling() { return Vector(0, 0, 0); };
    // evaluacion de la funcion de distribucion
    virtual Color eval_f() { return Color(0, 0, 0); };
    // retorna la emmitancia
    virtual Color emmitance() { return Color(0, 0, 0); };

    // destructor
    // ~Material();
};

class DifusseOG : public Material
{
private:
    double tetha;
    Color albedo;
    Color emmitanceColor;

public:
    DifusseOG(Color _albedo, Color _emmitanceColor) : albedo(_albedo), emmitanceColor(_emmitanceColor){};
    virtual double probability() override
    {
        return (1.0 / pi) * cos(tetha);
    }
    virtual Vector sampling() override
    {
        return random_cosine_hemisphere(tetha);
    }
    virtual Color eval_f() override
    {
        return albedo * (1.0 / pi);
    }
    virtual Color emmitance() override
    {
        return emmitanceColor;
    }
};
class DifusseON : public Material
{
private:
    double tetha;
    Color albedo;
    Color emmitanceColor;
    double sigma;

public:
    DifusseON(Color _albedo, Color _emmitanceColor, double _sigma) : albedo(_albedo), emmitanceColor(_emmitanceColor), sigma(_sigma){};

    virtual double probability() override
    {
        return (1.0 / pi) * cos(tetha);
    }
    virtual Vector sampling() override
    {
        return random_cosine_hemisphere(tetha);
    }
    virtual Color eval_f() override
    {
        double A, B;
        double sigma2 = sigma * sigma;
        double alpha,betha;

        A = 1.0 - (sigma2 / (2.0 * (sigma2 + 0.33)));
        B = (0.45 * sigma2) / (sigma2 + 0.09);

        Color fr = albedo * (1.0 / pi);
        fr = fr * (A + B * sin(alpha) * tan(betha));
        return fr;
    }
    virtual Color emmitance() override
    {
        return emmitanceColor;
    }
};

// class Metal : Material
// {
// };

// class Dielectric : Material
// {
// };
#endif

// virtual Color BDRF(const Ray &wi, Ray &wo, const registro &rec) override
// {

//     wo.o = rec.x;
//     wo.d = random_coseno(theta);
//     double pdf = (1.0 / pi) * cos(theta);

//     double ro, phio, thetao;
//     CarteEsfericas(wo.d, thetao, phio, ro);

//     double ri, thetai, phii;
//     CarteEsfericas(wi.d, thetai, phii, ri);

//     double A = 1.0 - sigma / (2.0 * (sigma + 0.33));
//     double cosio = cos((phii - phio));
//     double alpha = (thetai > thetao) ? thetai : thetao;
//     double beta = (thetai < thetao) ? thetai : thetao;

//     double B = 0;
//     if (cosio = !0)
//     {
//         B = (0.45 * sigma) / (sigma + 0.09) * cosio * sin(alpha) * tan(beta);
//     }

//     wo.d = LocalGlobal(rec.n, wo.d);

//     return albedo * (1.0 / pi)(A + B)(1.0 / pdf); //(1.0/pdf)(1.0/pdf)
// }