#ifndef MATERIALH
#define MATERIALH

#pragma once
#include "vec.h"
#include "ray.h"
#include "utils.h"
#include "sphere.h"
// #include "rt.cpp"

class Material
{
public:
    Color emmitanceColor;
    // probabilidad
    virtual double probability() { return 0.0; };
    // muestreo de direcciones
    virtual Vector sampling(const Ray wi, const Vector N, Point x) { return Vector(0, 0, 0); };
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
    

public:
    Color emmitanceColor;
    DifusseOG(Color _albedo, Color _emmitanceColor) : albedo(_albedo), emmitanceColor(_emmitanceColor){};
    virtual double probability() override
    {
        pdf = (1.0 / pi) * cos(tetha);
        return pdf;
    }
    virtual Vector sampling(const Ray wi, const Vector N, Point x) override
    {
        Vector s, t;
        Vector wo = random_cosine_hemisphere(tetha);
        coordinateSystem(N, s, t);
        wo = makeGlobal(wo, N, s, t).normalize();
        return wo;
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
    DifusseON(Color _albedo, Color _emmitanceColor, double _sigma) : albedo(_albedo), emmitanceColor(_emmitanceColor), sigma(_sigma)
    {
        // precalculamos los valores de A y B para el modelo de Oren-Nayar
        double sigma2 = sigma * sigma;
        A = 1.0 - (sigma2 / (2.0 * (sigma2 + 0.33)));
        B = (0.45 * sigma2) / (sigma2 + 0.09);
    };

    virtual double probability() override
    {
        pdf = (1.0 / pi) * cos(tetha);
        return pdf;
    }
    virtual Vector sampling(const Ray wi, const Vector N, Point X) override
    {
        Vector s, t;
        Vector wo = random_cosine_hemisphere(tetha);
        coordinateSystem(N, s, t);
        wo = makeGlobal(wo, N, s, t).normalize();
        return wo;
    }
    virtual Color eval_f(const Ray wi, const Ray wo, const Vector N) override
    {
        double alpha, betha;

        double NdotL = N.dot(wi.d);
        double NdotV = N.dot(wo.d);

        double thetaI = acos(NdotL);
        double thetaO = acos(NdotV);

        alpha = max(thetaI, thetaO);
        betha = min(thetaI, thetaO);

        Point l = (wi.d - N * NdotL).normalize();
        Point v = (wo.d - N * NdotV).normalize();

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

Vector reflect(const Ray r, Vector N)
{
    // Vector wr = (r.d *-1) - N * 2.0 * (N.dot(r.d));
    Vector wr = r.d - N * 2 * (N.dot(r.d));
    return wr.normalize();
}

Vector refract(const Vector wi, Vector N, double etta, double cosTi, double cosTt)
{
    Vector wt = ((wi * -1) * etta) + N * (etta * cosTi - cosTt);
    return wt.normalize();
}

class Specular : public Material
{
private:
    Color emmitanceColor;
    Color albedo;

public:
    Specular(Color _albedo, Color _emmitanceColor) : albedo(_albedo), emmitanceColor(_emmitanceColor){};
    virtual double probability() override
    {
        return 1.0;
    }
    virtual Vector sampling(const Ray wi, const Vector N, Point x) override
    {
        Vector wo = reflect(wi, N);
        return wo.normalize();
    }
    virtual Color eval_f(const Ray wi, const Ray wo, const Vector N) override
    {
        return albedo;
    }
    virtual Color emmitance() override
    {
        return emmitanceColor;
    }
};


class Glass : public Material
{
private:
    Color emmitanceColor;
    Color albedo;
    double NI = 1.0;
    double NT;
    double F, n2;
    double rll, rPerp;
    double fr, ft;
    double pr, pt;
    double nvwr, nvwt;
    Vector wr, wt;
    bool bounce;

public:
    Glass(Color _albedo, Color _emmitanceColor, double _nt) : albedo(_albedo), emmitanceColor(_emmitanceColor), NT(_nt) {}

    double probability() override
    {
        if (bounce)
        {
            return pr;
        }
    
        return pt;
    }

    Vector sampling(const Ray r, const Vector n, Vector nv) override
    {
        double ni = NI;
        double nt = NT;
        Vector wi = r.d * -1;

        double cosTin = (wi).dot(n);
        cosTin = clamp(cosTin, -1, 1);
        bool into = cosTin > 0.f;
        if (!into)
        {
            swap(ni, nt);
            cosTin = fabs(cosTin);
        }

        n2 = ni / nt;

        double sinTi = sqrt(max(0.0, 1.0 - cosTin * cosTin));
        double sinTt = n2 * sinTi;
        double cosTrans = sqrt(max(0.0, 1 - sinTt * sinTt));

        if (sinTt >= 1)
        {
            F = 0.0;
        }
        else
        {
            rll = ((nt * cosTin) - (ni * cosTrans)) / ((nt * cosTin) + (ni * cosTrans));
            rPerp = ((ni * cosTin) - (nt * cosTrans)) / ((ni * cosTin) + (nt * cosTrans));
            F = (rll * rll + rPerp * rPerp) * 0.5;
        }

        wr = reflect(r, nv);

        wr.normalize();
        double cosTr = wr.dot(nv);
        wt = refract(wi, nv, n2, cosTin, cosTrans);
        wt.normalize();

        bounce = random_double() < F;

        pr = F;
        pt = 1.0 - F;

        fr = F / fabs(cosTin);
        ft = (((nt * nt) / (ni * ni)) * (1 - F)) / fabs(cosTrans);

        pr = ((fr * fabs(nv.dot(wr)) / pr));
        pt = ((ft * fabs(nv.dot(wt)) / pt));

        // cout << "pr: " << pr << " pt: " << pt << endl;

        if (bounce)
        {
            return wr;
        }

        return wt;
    }

    Color eval_f(const Ray wi, const Ray wo, const Vector N) override
    {
        return albedo; 
    }

    Color emmitance() override
    {
        return emmitanceColor;
    }
};

#endif