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
    virtual double probability(const Ray wi, const Vector wo, hitInfo &hitInfo) { return 0.0; };
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
    virtual double probability(const Ray wi, const Vector wo, hitInfo &hitInfo) override
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
    virtual double probability(const Ray wi, const Vector wo, hitInfo &hitInfo) override
    {
        pdf = (1.0 / pi) * fabs(cos(tetha));
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
            // cout<< hitInfo.N.x<<" "<<hitInfo.N.y<<" "<<hitInfo.N.z<<endl;
        }
        // obedecer a la normal de la superficie.

        coordinateSystem(hitInfo.N, s, t);
        wo = makeGlobal(wo, hitInfo.N, s, t).normalize();
        return wo;
    }
    virtual Color eval_f(const Ray wi, const Ray wo, hitInfo &hitInfo) override
    {
        // return hitInfo.N * (1.0 / pi);

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

    virtual double probability(const Ray wi, const Vector wo, hitInfo &hitInfo) override
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
    virtual double probability(const Ray wi, const Vector wo, hitInfo &hitInfo) override
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

    double probability(const Ray wi, const Vector wo, hitInfo &hitInfo) override
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

class Metal : public Material
{
private:
    Color eta;
    Color kappa;
    double alpha;
    Vector wh;
    double air_eta = 1.00029;
    // Vector n;
    Vector wo;

public:
    Metal(Color _eta, Color _kappa, double _alpha) : eta(_eta.mult(_eta)), kappa(_kappa.mult(kappa)), alpha(_alpha){};
    virtual double probability(const Ray wi, const Vector wo, hitInfo &hitInfo) override
    {
        // wh = ((wi.d * -1.0) + wo);

        // wh.normalize();
        double cos_H = hitInfo.N.dot(wh);
        // double d = D(cos_H);
        double p = D(wh, hitInfo) * cos_H / fabs(wo.dot(wh));
        return p;
    }
    virtual Vector sampling(const Ray wi, hitInfo &hitInfo) override
    {
        // n = hitInfo.N;
        Vector dir = random_dir();
        Vector s, t;
        coordinateSystem(hitInfo.N, s, t);
        wh = s * dir.x + t * dir.y + hitInfo.N * dir.z;
        // Vector v = (wi.d * (-1.0)).normalize();
        wh.normalize();
        wo = reflect(wi, wh);
        return wo.normalize();
    }
    virtual Color eval_f(const Ray wi, const Ray wo, hitInfo &hitInfo) override
    {

        Vector v = wi.d * (-1.0);
        double cosv1 = v.dot(hitInfo.N);
        double cosv2 = wo.d.dot(hitInfo.N);

        double cosih = v.dot(wh);
        double cosh = wh.dot(hitInfo.N);

        double G_ = 0.0;

        if ((v.dot(wh) / cosv1) > 0.f && (wo.d.dot(wh) / cosv2) > 0.f)
        {
            G_ = G(cosv1) * G(cosv2);
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

        Color fr = Color(R, G, B) * D(wh, hitInfo) * G_ * (1.0 / (fabs(cosv1) * fabs(cosv2)));

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
        double phi = 2.0 * pi * r2;

        double x = sin(theta) * cos(phi);
        double y = sin(theta) * sin(phi);
        double z = cos(theta);

        return Vector(x, y, z);
    }
    double G(const double cosV)
    {
        // G shading term
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
    double D(Vector wh, hitInfo &hitInfo)
    {
        // float tan2 = tan2Theta(wh);
        // if (std::isinf(tan2))
        //     return 0;

        // float cos4 = cos2Theta(wh) * cos2Theta(wh);
        // return std::exp(-tan2 / (alpha * alpha)) / (M_PI * alpha * alpha * cos4);

        double costh = hitInfo.N.dot(wh);
        if (costh > 0)
        {
            double alpha2 = alpha * alpha;
            double costh2 = costh * costh;
            double costh4 = costh2 * costh2;
            double tang = sqrt(1.0 - costh2) / costh;
            return exp(-((tang * tang) / alpha2)) / (pi * alpha2 * costh4);
        }
        else
            return 0.0;
    }
    double fresnel(const double etap, const double kappap, const double cos)
    {
        double sente = sin(acos(cos));
        double nksin = etap - kappap - sente * sente;

        double ab = sqrt(nksin * nksin + etap * kappap * 4.0);
        double a = sqrt((ab + nksin) / 2.0);

        double rpernum = ab + cos * cos - 2.0 * a * cos;
        double rperdem = ab + cos * cos + 2.0 * a * cos;
        double rper = rpernum / rperdem;

        double rparnum = ab * cos * cos + sente * sente * sente * sente - 2.0 * a * cos * sente * sente;
        double rpardem = ab * cos * cos + sente * sente * sente * sente + 2.0 * a * cos * sente * sente;
        double rpar = rper * rparnum / rpardem;

        return (1.0 / 2.0) * (rper + rpar);

        // double coste = clamp(cos, -1.0, 1.0);

        // double sente = sqrt(1.0 - coste * coste);
        // double nksin = etap - kappap - sente * sente;

        // double ab = sqrt(nksin * nksin + etap * kappap * 4.0);
        // double a = sqrt((ab + nksin) * 0.5);

        // double rpernum = ab + coste * coste - 2.0 * a * coste;
        // double rperdem = ab + coste * coste + 2.0 * a * coste;
        // double rper = rpernum / rperdem;

        // double rparnum = ab * coste * coste + sente * sente * sente * sente - 2.0 * a * coste * sente * sente;
        // double rpardem = ab * coste * coste + sente * sente * sente * sente + 2.0 * a * coste * sente * sente;
        // double rpar = rper * (rparnum / rpardem);

        // return (0.5) * (rper + rpar);
    };
};

// inline double fresnel(const Vector &wo, const Vector &n, float n1, float n2)
// {
//     const float f0 = std::pow((n1 - n2) / (n1 + n2), 2.0f);
//     return f0 + (1.0f - f0) * std::pow(1.0f - wo.dot(n), 5.0f);
// }

// class Phong : public Material
// {
// public:
//     double diff;
//     double spec;
//     Texture *albedo;
//     Texture *normalMap;

//     Phong(Texture *_albedo, double _diff, double _spec) : albedo(_albedo), diff(_diff), spec(_spec){};
//     Phong(Texture *_albedo, Texture *_normalMap) : albedo(_albedo), normalMap(_normalMap) {}

//     virtual double probability(const Ray wi, const Vector wo, hitInfo &hitInfo) override
//     {
//         return 1.0/(2.0 * pi);
//     }
//     virtual Vector sampling(const Ray wi, hitInfo &hitInfo) override
//     {
//         double r1 = random_double();
//         double r2 = random_double();

//         double theta = acos(pow(1 - r1, 1.0 / (spec + 1)));
//         double phi = 2 * pi * r2;

//         Vector wr = Vector(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
//         wr = wr.normalize();

//         Vector u = Vector(0.00424, 1, 0.00764)%wr;
//         u = u.normalize();
//         Vector v = u % wr;
//         v = v.normalize();

//         Vector WI = wr * cos(theta) + u * sin(theta) * cos(phi) + v * sin(theta) * sin(phi);
//         WI.normalize();

//         return WI;
//     }
//     virtual Color eval_f(const Ray wi, const Ray wo, hitInfo &hitInfo) override
//     {
//         Vector wr = reflect(wi, hitInfo.N);
//         wr.normalize();

//         Color albedoColor = albedo->value(hitInfo.u, hitInfo.v, hitInfo.x);
//         Color diffuseColor = albedoColor * diff * max(wi.d.dot(wo.d), (double)0.0);
//         Color specularColor = albedoColor * spec * pow(cosPhi, 100);

//         return diffuseColor + specularColor;
//     }
//     virtual Color emmitance() override
//     {
//         return Color(0, 0, 0);
//     }

// };

#endif