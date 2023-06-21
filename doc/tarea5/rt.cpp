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
using namespace std;

// materiales para la escena

DifusseOG m1(Color(0.75, 0.25, 0.25), Color());
DifusseOG m2(Color(.25, .25, .75), Color());     // pared der
DifusseOG m3(Color(.25, .75, .25), Color());     // pared detras
DifusseOG m4(Color(.25, .75, .75), Color());     // suelo
DifusseOG m5(Color(.75, .75, .25), Color());     // techo
DifusseOG m6(Color(.2, .3, .4), Color());        // esfera abajo-izq
DifusseOG m7(Color(.4, .3, .2), Color());        // esfera abajo-der
DifusseOG m8(Color(1, 1, 1), Color(10, 10, 10)); // esfera arriba

Sphere spheres[] = {
    // Escena: radio, posicion, color, luz
    Sphere(1e5, Point(-1e5 - 49, 0, 0), &m1),    // pared izq
    Sphere(1e5, Point(1e5 + 49, 0, 0), &m2),     // pared der
    Sphere(1e5, Point(0, 0, -1e5 - 81.6), &m3),  // pared detras
    Sphere(1e5, Point(0, -1e5 - 40.8, 0), &m4),  // suelo
    Sphere(1e5, Point(0, 1e5 + 40.8, 0), &m5),   // techo
    Sphere(16.5, Point(-23, -24.3, -34.6), &m6), // esfera abajo-izq
    Sphere(16.5, Point(23, -24.3, -3.6), &m7),   // esfera abajo-der
    Sphere(10.5, Point(0, 24.3, 0), &m8)         // esfera arriba
    // Sphere(0, Point(0, 24.3, 0), Color(1, 1, 1), Color(10, 10, 10)) // esfera arriba
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

    for (int i = 0; i < spheresLength; i++)
    {
        // asignación de valores al par
        spheresData[i].first = i;
        spheresData[i].second = spheres[i].intersect(r);
        if (spheresData[i].second > 0)
        {
            id = spheresData[i].first;
            t = spheresData[i].second;
        }
    }
    // ordenamiento de pares por la menor distancia
    for (int i = 0; i < spheresLength; i++)
    {
        if (t > spheresData[i].second && spheresData[i].second > 0.01)
        {
            // actualizamos valores de t e id solo si el valor almacenado en el par es menor al valor de t y que sea positivo
            id = spheresData[i].first;
            t = spheresData[i].second;
        }
    }
    if (t > 0)
    {
        return true;
    }
    return false;
}

// Calcula el valor de color para el rayo dado, seleccionar muestreo aleatorio
// 0 para esfera unitaria
// 1 para esfera hemisferica
// 2 para coseno hemisferica
// regresar el color resultante
Color shade(const Ray &r, int depth, int mode)
{
    double t;
    int id = 0;

    if (depth <= 0)
        return Color(0, 0, 0);

    if (!intersect(r, t, id))
        // el rayo no intersecto objeto, return Vector() == negro
        return Color(0, 0, 0);

    const Sphere &obj = spheres[id];

    // determinar coordenadas del punto de interseccion
    Point x = r.o + r.d * t;

    // determinar la dirección normal en el punto de interseccion
    Vector n = (x - obj.p).normalize();

    Point target;
    double p;
    double tetha;

    target = obj.m->sampling();
    p = obj.m->probability();

    // se arma el sistema de coordenadas
    Vector s;
    Vector ta;
    coordinateSystem(n, s, ta);

    Point dir = makeGlobal(target, n, s, ta);

    Ray newRay(x, dir);

    double cos_theta = newRay.d.dot(n);
    Color BRDF = obj.m->eval_f();
    Color emittance = obj.m->emmitance();

    // Color BRDF = obj.c * (1 / pi)
    // Color emittance = obj.l;

    Color incomingColor = shade(newRay, depth - 1, mode);
    Color colorValue = emittance + (incomingColor.mult(BRDF) * cos_theta) * (1.0 / p);
    return colorValue;
}

// funcion principal
int main(int argc, char *argv[])
{
    int w = 1024, h = 768; // image resolution
    // Numero de muestras por pixel.
    const int pixel_samples = 32;
    // Numero de rebotes.
    const int depth = 3;
    // Modo de muestreo. 0 para esfera unitaria, 1 para esfera hemisferica, 2 para coseno hemisferica
    const int mode = 2;

    // fija la posicion de la camara y la dirección en que mira
    Ray camera(Point(0, 11.2, 214), Vector(0, -0.042612, -1).normalize());

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
                    pixelValue = pixelValue + shade(Ray(camera.o, cameraRayDir.normalize()), depth, mode);
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