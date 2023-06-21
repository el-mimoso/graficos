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
    virtual double probability(){return 0.0;};
    // muestreo de direcciones
    virtual Vector sampling(){ return Vector(0,0,0);};
    // evaluacion de la funcion de distribucion
    virtual Color eval_f(){return Color(0,0,0);};
    // retorna la emmitancia
    virtual Color emmitance(){ return Color(0,0,0);};
    
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
    DifusseOG(Color _albedo, Color _emmitanceColor) : albedo(_albedo), emmitanceColor(_emmitanceColor) {};
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

// class DiffuseON : Material
// {
// };

// class Metal : Material
// {
// };

// class Dielectric : Material
// {
// };
#endif

