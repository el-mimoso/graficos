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
    Color emmitanceColor;

public:
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
double fresnel(const double etap, const double kappap, const double cos)
{
    // double sente = sin(acos(cos));
    // double nksin = etap - kappap - sente * sente;

    // double ab = sqrt(nksin * nksin + etap * kappap * 4.0);
    // double a = sqrt((ab + nksin) / 2.0);

    // double rpernum = ab + cos * cos - 2.0 * a * cos;
    // double rperdem = ab + cos * cos + 2.0 * a * cos;
    // double rper = rpernum / rperdem;

    // double rparnum = ab * cos * cos + sente * sente * sente * sente - 2.0 * a * cos * sente * sente;
    // double rpardem = ab * cos * cos + sente * sente * sente * sente + 2.0 * a * cos * sente * sente;
    // double rpar = rper * rparnum / rpardem;

    // return (1.0 / 2.0) * (rper + rpar);

    double coste = Clamp(cos, -1.0, 1.0);

    double sente = sqrt(1.0 - coste * coste);
    double nksin = etap - kappap - sente * sente;

    double ab = sqrt(nksin * nksin + etap * kappap * 4.0);
    double a = sqrt((ab + nksin) * 0.5);

    double rpernum = ab + coste * coste - 2.0 * a * coste;
    double rperdem = ab + coste * coste + 2.0 * a * coste;
    double rper = rpernum / rperdem;

    double rparnum = ab * coste * coste + sente * sente * sente * sente - 2.0 * a * coste * sente * sente;
    double rpardem = ab * coste * coste + sente * sente * sente * sente + 2.0 * a * coste * sente * sente;
    double rpar = rper * (rparnum / rpardem);

    return (0.5) * (rper + rpar);
};

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
        return Color(0, 0, 0);
    }
};

class Metal : public Material
{
private:
    Color eta;
    Color kappa;
    double alpha;
    Vector wh;
    double air_eta = 1.00029;
    Vector n;
    Vector wo;

public:
    Metal(Color _eta, Color _kappa, double _alpha) : eta(_eta), kappa(_kappa), alpha(_alpha){};
    virtual double probability() override
    {
        double cos_H = n.dot(wh);
        double p = cos_H / fabs(wo.dot(wh));
        return p;
    }
    virtual Vector sampling(const Ray wi, const Vector N, Point x) override
    {
        n = N;
        Vector dir = random_dir();
        Vector s, t;
        coordinateSystem(N, s, t);
        wh = s * dir.x + t * dir.y + N * dir.z;
        Vector v = (wi.d * (-1.0)).normalize();
        wo = reflect(wi, wh);
        return wo.normalize();
    }
    virtual Color eval_f(const Ray wi, const Ray wo, const Vector N) override
    {

        Vector v = wi.d * (-1.0);
        double cosv1 = v.dot(N);
        double cosv2 = wo.d.dot(N);

        double cosih = v.dot(wh);
        double cosh = wh.dot(N);

        double smith = 0.0;

        if ((v.dot(wh) / cosv1) > 0.f && (wo.d.dot(wh) / cosv2) > 0.f)
        {
            smith = Smith_G(cosv1) * Smith_G(cosv2);
        }

        double cos_tetha = wh.dot(v);

        Color eta_aux = eta.mult(eta) * (1.0 / (air_eta * air_eta));      //*(1.0/(1.00029*1.00029))
        Color kapa_aux = kappa.mult(kappa) * (1.0 / (air_eta * air_eta)); //

        double R = fresnel(eta_aux.x, kapa_aux.x, cos_tetha);
        double G = fresnel(eta_aux.y, kapa_aux.y, cos_tetha);
        double B = fresnel(eta_aux.z, kapa_aux.z, cos_tetha);

        // double R = fresnel(eta.x, kappa.x, cos_tetha);
        // double G = fresnel(eta.y, kappa.y, cos_tetha);
        // double B = fresnel(eta.z, kappa.z, cos_tetha);

        // double delta = (wo.d.dot(N) > 0.0001) ? 1.0 : 0.0;

        Color fr = Color(R, G, B) * smith * (1.0 / (fabs(cosv1) * fabs(cosv2)));

        return fr;
        // return Color(1,1,1);
    }
    virtual Color emmitance() override
    {
        return Color(0, 0, 0);
    }
    Vector random_dir()
    {
        double r1 = random_double();
        double r2 = random_double();

        double theta = atan(sqrt(-alpha * alpha * log10(1.0 - r1)));
        double phi = 2.0 * pi * r1;

        double x = sin(theta) * cos(phi);
        double y = sin(theta) * sin(phi);
        double z = cos(theta);

        return Vector(x, y, z);
    }
    double Smith_G(const double cosV)
    {
        // Smith G shading term
        double cosv2 = cosV * cosV;
        double a = (cosV / (sqrt(1.0 - cosv2) * alpha));
        double a2 = a * a;

        if (a < 1.6)
        {
            return ((3.535 * a) + (2.181 * a2)) / (1 + (2.276 * a) + (2.577 * a2));
        }
        else
            return 1.0;
    }
};

// class Dielectric : Material
// {
// };
#endif