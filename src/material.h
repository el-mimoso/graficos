#ifndef MATERIALH
#define MATERIALH

#pragma once
#include "vec.h"
#include "ray.h"
#include "utils.h"
#include "sphere.h"
#include "texture.h"
// #include "rt.cpp"

struct hitInfo
{
    // int id;
    double t;
    double u;
    double v;
    Point x;
    Vector N;
    Vector nv;
    // Ray wi;
};

class Material
{
public:
    // probabilidad
    virtual double probability(hitInfo &hitInfo) { return 0.0; };
    // muestreo de direcciones
    virtual Vector sampling(const Ray wi, hitInfo &hitInfo) { return Vector(0, 0, 0); };
    // evaluacion de la funcion de distribucion
    // wi direccion de entrada
    // wo direccion de salida
    // N normal
    virtual Color eval_f(const Ray wi, const Ray wo, hitInfo &hitInfo) { return Color(0, 0, 0); };
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
    Color emmitanceColor = Color(0, 0, 0);

public:
    DifusseOG(Color _albedo, Color _emmitanceColor) : albedo(_albedo), emmitanceColor(_emmitanceColor){};
    virtual double probability(hitInfo &hitInfo) override
    {
        pdf = (1.0 / pi) * cos(tetha);
        return pdf;
    }
    virtual Vector sampling(const Ray wi, hitInfo &hitInfo) override
    {
        Vector s, t;
        Vector wo = random_cosine_hemisphere(tetha);
        coordinateSystem(hitInfo.N, s, t);
        wo = makeGlobal(wo, hitInfo.N, s, t).normalize();
        return wo;
    }
    virtual Color eval_f(const Ray wi, const Ray wo, hitInfo &hitInfo) override
    {
        return albedo * (1.0 / pi);
    }
    virtual Color emmitance() override
    {
        return emmitanceColor;
    }
};

class DifusseTx : public Material
{
private:
    double tetha;
    double pdf;
    // Color albedo;
    Texture *albedo;
    Texture *normalMap;
    Color emmitanceColor = Color(0, 0, 0);

public:
    // Cuando solo tenemos un mapa de textura
    DifusseTx(Texture *_albedo) : albedo(_albedo){};
    // Cuando tenemos un mapa de textura y un mapa de normales
    DifusseTx(Texture *_albedo, Texture *_normalMap) : albedo(_albedo), normalMap(_normalMap){};
    virtual double probability(hitInfo &hitInfo) override
    {
        pdf = (1.0 / pi) * cos(tetha);
        return pdf;
    }
    virtual Vector sampling(const Ray wi, hitInfo &hitInfo) override
    {
        Vector s, t;
        Vector wo = random_cosine_hemisphere(tetha);
        // obedecer al mapa de normales
        if (normalMap != nullptr)
        {
            // cout << "normalMap" << endl;
            //
 
            Vector normal = normalMap->value(hitInfo.u, hitInfo.v, hitInfo.x);
            normal = normal * 2.0 - Vector(1, 1, 1);
            normal = normal.normalize();    
            hitInfo.N = normal;
            cout<< hitInfo.N.x<<" "<<hitInfo.N.y<<" "<<hitInfo.N.z<<endl;
        }
        // obedecer a la normal de la superficie.

        coordinateSystem(hitInfo.N, s, t);
        wo = makeGlobal(wo, hitInfo.N, s, t).normalize();
        return wo;
    }
    virtual Color eval_f(const Ray wi, const Ray wo, hitInfo &hitInfo) override
    {
        // return hitInfo.N * (1.0 / pi);
        // return Color(hitInfo.N.x + 1, hitInfo.N.y + 1, hitInfo.N.z + 1) * .5;
        return albedo->value(hitInfo.u, hitInfo.v, hitInfo.x) * (1.0 / pi);
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

    virtual double probability(hitInfo &hitInfo) override
    {
        pdf = (1.0 / pi) * cos(tetha);
        return pdf;
    }
    virtual Vector sampling(const Ray wi, hitInfo &hitInfo) override
    {
        Vector s, t;
        Vector wo = random_cosine_hemisphere(tetha);
        coordinateSystem(hitInfo.N, s, t);
        wo = makeGlobal(wo, hitInfo.N, s, t).normalize();
        return wo;
    }
    virtual Color eval_f(const Ray wi, const Ray wo, hitInfo &hitInfo) override
    {
        double alpha, betha;

        double NdotL = hitInfo.N.dot(wi.d);
        double NdotV = hitInfo.N.dot(wo.d);

        double thetaI = acos(NdotL);
        double thetaO = acos(NdotV);

        alpha = max(thetaI, thetaO);
        betha = min(thetaI, thetaO);

        Point l = (wi.d - hitInfo.N * NdotL).normalize();
        Point v = (wo.d - hitInfo.N * NdotV).normalize();

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
    virtual double probability(hitInfo &hitInfo) override
    {
        return 1.0;
    }
    virtual Vector sampling(const Ray wi, hitInfo &hitInfo) override
    {
        Vector wo = reflect(wi, hitInfo.N);
        return wo.normalize();
    }
    virtual Color eval_f(const Ray wi, const Ray wo, hitInfo &hitInfo) override
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

    double probability(hitInfo &hitInfo) override
    {
        if (bounce)
        {
            return pr;
        }

        return pt;
    }

    Vector sampling(const Ray r, hitInfo &hitInfo) override
    {
        double ni = NI;
        double nt = NT;
        Vector wi = r.d * -1;

        double cosTin = (wi).dot(hitInfo.N);
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

        wr = reflect(r, hitInfo.nv);

        wr.normalize();
        double cosTr = wr.dot(hitInfo.nv);
        wt = refract(wi, hitInfo.nv, n2, cosTin, cosTrans);
        wt.normalize();

        bounce = random_double() < F;

        pr = F;
        pt = 1.0 - F;

        fr = F / fabs(cosTin);
        ft = (((nt * nt) / (ni * ni)) * (1 - F)) / fabs(cosTrans);

        pr = ((fr * fabs(hitInfo.nv.dot(wr)) / pr));
        pt = ((ft * fabs(hitInfo.nv.dot(wt)) / pt));

        // cout << "pr: " << pr << " pt: " << pt << endl;

        if (bounce)
        {
            return wr;
        }

        return wt;
    }

    Color eval_f(const Ray wi, const Ray wo, hitInfo &hitInfo) override
    {
        return albedo;
    }

    Color emmitance() override
    {
        return emmitanceColor;
    }
};

#endif