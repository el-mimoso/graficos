// rt: un lanzador de rayos minimalista
// g++ -O3 -fopenmp rt.cpp -o rt
#include "vec.h"
#include "ray.h"
#include "sphere.h"
#include "utils.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <utility>
#include <typeinfo>

using namespace std;

// materiales para la escena

DifusseOG m1(Color(0.75, 0.25, 0.25), Color());  // pared izq
DifusseOG m2(Color(.25, .25, .75), Color());     // pared der
DifusseOG m3(Color(.25, .75, .25), Color());     // pared detras
DifusseOG m4(Color(.25, .25, .30), Color());     // suelo
DifusseOG m5(Color(.25, .25, .25), Color());     // techo
// Specular m6(Color(1, 1, 1), Color());
Specular blueMirror(Color(.70, .70, 1), Color());
Specular goldMirror(Color(.7, 1, 0.8), Color());

Specular m6(Color(1, 1, 1), Color());                           // esfera abajo-izq Espejo
Glass redGlass(Color(1, .76, .76), Color(), 1.5);               // esfera abajo-der Glass
Glass glass(Color(1, 1, 1), Color(), 1.5);               // esfera abajo-der Glass
DifusseON m8(Color(1, 1, 1), Color(10, 10, 10), 0.5); // esfera arriba LUZ

Sphere spheres[] = {
    // Escena: radio, posicion, color, luz

    Sphere(1e5, Point(-1e5 - 150, 0, 0), &m2),  // pared izq
    Sphere(1e5, Point(1e5 + 150, 0, 0), &m6),   // pared der
    Sphere(1e5, Point(0, 0, -1e5 - 350), &m6),  // pared detras
    Sphere(1e5, Point(0, -1e5 - 40.8, 0), &m4), // suelo
    Sphere(1e5, Point(0, 1e5 + 150.8, 0), &m5), // techo
    Sphere(1e5, Point(0, 0, 1e5 + 350), &m4),  // pared detras

    Sphere(16.5, Point(-50, -24.3, -150), &redGlass),
    Sphere(16.5, Point(0.0, -24.3, -150), &m2),
    Sphere(16.5, Point(50, -24.3, -150), &m3),

    Sphere(16.5, Point(-50, -24.3, 0), &goldMirror),
    Sphere(16.5, Point(0, -24.3, 0), &m6),
    Sphere(16.5, Point(50, -24.3, 0), &blueMirror),

    Sphere(16.5, Point(15, 24.3, 0), &glass),  // esfera arriba
    Sphere(30.5, Point(-150, 60, -80), &m8), // LUZ
    Sphere(40.5, Point(135, 100, -80), &m8)    // LUZs

};

const int spheresLength = sizeof(spheres) / sizeof(spheres[0]);

// calcula la intersección del rayo r con todas las esferas
// regresar true si hubo una intersección, falso de otro modo
// almacenar en t la distancia sobre el rayo en que sucede la interseccion
// almacenar en id el indice de spheres[] de la esfera cuya interseccion es mas cercana
inline bool intersect(const Ray &r, double &t, int &id)
{
    // arreglo de pares donde almacenamos la id y la distancia del rayo
    pair<int, double> spheresData[spheresLength];
    double t_maxDistance = t = 100000000000;

    for (int i = 0; i < spheresLength; i++)
    {
        // asignación de valores al par
        spheresData[i].first = i;
        spheresData[i].second = spheres[i].intersect(r);

        if (spheresData[i].second && spheresData[i].second > 0.01 && spheresData[i].second < t)
        {
            // actualizamos valores de t e id solo si el valor almacenado en el par es menor al valor de t y que sea positivo
            id = spheresData[i].first;
            t = spheresData[i].second;
        }
    }
    if (t < t_maxDistance)
    {
        return true;
    }
    return false;
}

// Calcula el valor de color para el rayo dado, seleccionar muestreo aleatorio
// regresar el color resultante
Color shade(const Ray &r, int depth)
{
    double t;
    int id = 0;

    if (depth <= 0)
        return Color(.3, .3, .3);

    if (!intersect(r, t, id))
        // el rayo no intersecto objeto, return Vector() == negro
        return Color(.3, .3, .3);

    const Sphere &obj = spheres[id];

    // determinar coordenadas del punto de interseccion
    Point x = r.o + r.d * t;

    // determinar la dirección normal en el punto de interseccion
    Vector n = (x - obj.p).normalize();
    Vector nv = checkNormalDir(n, r);

    Material *material = obj.m;
    Color emittance = material->emmitance();

    // terminamos el camino si pegamos con una fuente de luz.
    if (emittance.dot(emittance) > 0.01)
    {
        return emittance;
    }

    Point wo;
    double p;
    double tetha;

    wo = material->sampling(r, n, x);
    p = material->probability();
    // direccion saliente
    Ray newRay(x, wo);
    double cos_theta = newRay.d.dot(n);
    Color BRDF = material->eval_f(r, newRay, n);

    // si el material es especular perfecto
    if (typeid(*material) == typeid(Specular))
    {
        Vector wr;
        wr = material->sampling(r, n, x);
        p = material->probability();
        newRay = Ray(x, wr);
        BRDF = material->eval_f(r, newRay, n);
        cos_theta = newRay.d.dot(n);
        return BRDF.mult(shade(Ray(x, wr), depth - 1));
    }
    // si el material es dielectrico
    if (typeid(*material) == typeid(Glass))
    {
        Vector wrt;
        // Se mandan los valores de n y nv para que se haga la reflexion interna
        wrt = material->sampling(r, nv, n);
        p = material->probability();
        // direccion del nuevo rayo
        newRay = Ray(x, wrt);
        // cout<<wrt.x<<" "<<wrt.y<<" "<<wrt.z<<endl;
        // evaluacion de la BRDF
        BRDF = material->eval_f(r, newRay, n);
        return BRDF.mult(shade(newRay, depth - 1));
    }

    Color incomingColor = shade(newRay, depth - 1);
    Color colorvalue = (incomingColor.mult(BRDF) * (fabs(cos_theta) / p));
    return colorvalue;
}

// funcion principal
int main(int argc, char *argv[])
{
    int w = 1024, h = 768; // image resolution
    // Numero de muestras por pixel.
    const int pixel_samples = 1000;
    // Numero de rebotes.
    const int depth = 10;

    // fija la posicion de la camara y la dirección en que mira
    Ray camera(Point(20,50.2, 200), Vector(-0.1, -0.3, -1).normalize());

    // parametros de la camara
    Vector cx = Vector(w * 0.5095 / h, 0., 0.);
    Vector cy = (cx % camera.d).normalize() * 0.5095;

    // auxiliar para valor de pixel y matriz para almacenar la imagen
    Color *pixelColors = new Color[w * h];

// se define bloque paralelo
#pragma omp parallel
    {
        // siguiente bucle for se ejecuta de manera paralela,
        // asignando un chunk a cada nucleo del CPU disponible en ese momento
#pragma omp for schedule(dynamic, 1)
        for (int y = 0; y < h; y++)
        {
            // recorre todos los pixeles de la imagen
            fprintf(stderr, "\r%5.2f%%", 100. * y / (h - 1));

            for (int x = 0; x < w; x++)
            {
                int idx = (h - y - 1) * w + x; // index en 1D para una imagen 2D x,y son invertidos
                Color pixelValue = Color();    // pixelValue en negro por ahora
                // para el pixel actual, computar la dirección que un rayo debe tener
                Vector cameraRayDir = cx * (double(x) / w - .5) + cy * (double(y) / h - .5) + camera.d;
                for (int i = 0; i < pixel_samples; i++)
                {
                    // computar el color del pixel para el punto que intersectó el rayo desde la camara
                    pixelValue = pixelValue + shade(Ray(camera.o, cameraRayDir.normalize()), depth);
                }
                pixelValue = pixelValue * (1.0 / pixel_samples);
                // limitar los tres valores de color del pixel a [0,1]
                pixelColors[idx] = Color(clamp(pixelValue.x), clamp(pixelValue.y), clamp(pixelValue.z));
            }
        }
    }

    fprintf(stderr, "\n");

    FILE *f = fopen("image.ppm", "w");
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int p = 0; p < w * h; p++)
    {
        // escribe todos los valores de los pixeles
        fprintf(f, "%d %d %d ", toDisplayValue(pixelColors[p].x), toDisplayValue(pixelColors[p].y),
                toDisplayValue(pixelColors[p].z));
    }
    fclose(f);

    delete[] pixelColors;

    return 0;
}