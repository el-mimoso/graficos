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
    // wi direccion de entrada
    // wo direccion de salida
    // N normal
    virtual Color eval_f(const Ray wi, const Ray wo, const Vector N) { return Color(0, 0, 0); };
    // retorna la emmitancia
    virtual Color emmitance() { return Color(0, 0, 0); };

    // destructor
    // ~Material();
};

class DifusseOG : public Material
{
private:
    double tetha;
    double pdf;
    Color albedo;
    Color emmitanceColor;


public:
    DifusseOG(Color _albedo, Color _emmitanceColor) : albedo(_albedo), emmitanceColor(_emmitanceColor){};
    virtual double probability() override
    {
        pdf = (1.0 / pi) * cos(tetha);
        return pdf;
    }
    virtual Vector sampling() override
    {
        return random_cosine_hemisphere(tetha);
    }
    virtual Color eval_f(const Ray wi, const Ray wo, const Vector N) override
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
    double pdf;
    Color albedo;
    Color emmitanceColor;

    double A, B;

    double sigma;

public:
    DifusseON(Color _albedo, Color _emmitanceColor, double _sigma) : albedo(_albedo), emmitanceColor(_emmitanceColor), sigma(_sigma){
        //precalculamos los valores de A y B para el modelo de Oren-Nayar
        double sigma2 = sigma * sigma;
        A = 1.0 - (sigma2 / (2.0 * (sigma2 + 0.33)));
        B = (0.45 * sigma2) / (sigma2 + 0.09);
    };

    virtual double probability() override
    {
        pdf = (1.0 / pi) * cos(tetha);
        return pdf;
    }
    virtual Vector sampling() override
    {
        return random_cosine_hemisphere(tetha);
    }
    virtual Color eval_f(const Ray wi, const Ray wo, const Vector N) override
    {
        double alpha,betha;
        
        double NdotL = N.dot(wi.d);
        double NdotV = N.dot(wo.d);

        double thetaI = acos(NdotL);
        double thetaO = acos(NdotV);

        alpha = max(thetaI, thetaO);
        betha = min(thetaI, thetaO);


        Point l = (wi.d - N*NdotL).normalize();
        Point v = (wo.d - N*NdotV).normalize();

        double cosPhi = l.dot(v);
        double cosMax = max(0.0, cosPhi);

        Color fr = albedo * (1.0 / pi);
        fr = fr * (A + B * cosMax * sin(alpha) * tan(betha));
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