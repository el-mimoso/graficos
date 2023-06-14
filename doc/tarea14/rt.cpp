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

// Escena Funky con mas fuentes de luz
// Sphere spheres[] = {
//     // Escena: radio, posición, color, luz
//     Sphere(1e5, Point(-1e5 - 49, 0, 0), Color(.75, .25, .25), Color()),            // pared izq
//     Sphere(1e5, Point(1e5 + 49, 0, 0), Color(.25, .25, .75), Color()),             // pared der
//     Sphere(1e5, Point(0, 0, -1e5 - 81.6), Color(.25, .75, .25), Color()),          // pared detras
//     Sphere(1e5, Point(0, -1e5 - 40.8, 0), Color(.25, .75, .75), Color()),          // suelo
//     Sphere(1e5, Point(0, 1e5 + 40.8, 0), Color(.75, .75, .25), Color(.1, .1, .1)), // techo
//     Sphere(16.5, Point(-23, -24.3, -34.6), Color(.2, .3, .4), Color(0, 0, 6)),     // esfera abajo-izq
//     Sphere(16.5, Point(23, -24.3, -3.6), Color(.4, .3, .2), Color(7, 0, 0)),       // esfera abajo-der
//     Sphere(10.5, Point(0, 24.3, 0), Color(1, 1, 1), Color(.2, .2, .2))             // esfera arriba
// };

// Escena OG
Sphere spheres[] = {
    Sphere(1e5, Point(-1e5 - 49, 0, 0), Color(.75, .25, .25), Color()),   // pared izq
    Sphere(1e5, Point(1e5 + 49, 0, 0), Color(.25, .25, .75), Color()),    // pared der
    Sphere(1e5, Point(0, 0, -1e5 - 81.6), Color(.25, .75, .25), Color()), // pared detras
    Sphere(1e5, Point(0, -1e5 - 40.8, 0), Color(.25, .75, .75), Color()), // suelo
    Sphere(1e5, Point(0, 1e5 + 40.8, 0), Color(.75, .75, .25), Color()),  // techo
    Sphere(16.5, Point(-23, -24.3, -34.6), Color(.2, .3, .4), Color()),   // esfera abajo-izq
    Sphere(16.5, Point(23, -24.3, -3.6), Color(.4, .3, .2), Color()),     // esfera abajo-der
    Sphere(10.5, Point(0, 24.3, 0), Color(1, 1, 1), Color(10, 10, 10))    // esfera arriba
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
// regresar el color resultante
Color shade(const Ray &r, int depth)
{
    double t;
    int id = 0;
    //condición de salida de la recursion introduciendo sesgo
    if (depth <= 0)
        return Color(0, 0, 0);

    if (!intersect(r, t, id))
        // el rayo no intersecto objeto, return Vector() == negro
        return Color(0, 0, 0);

    const Sphere &obj = spheres[id];

    // determinar coordenadas del punto de interseccion
    Point x = r.o + r.d * t;
    Vector n = (x - obj.p).normalize();

    Point target;
    double p;
    double tetha;

    // muestreo de coseno de hemisferio.
    target = random_cosine_hemisphere(tetha);
    p = (1.0 / pi) * cos(tetha);

    // se arma el sistema de coordenadas
    Vector s;
    Vector ta;
    coordinateSystem(n, s, ta);

    Point dir(
        Point(s.x, ta.x, n.x).dot(target),
        Point(s.y, ta.y, n.y).dot(target),
        Point(s.z, ta.z, n.z).dot(target));

    Ray newRay(x, dir);
    double cos_theta = newRay.d.dot(n);

    Color BRDF = obj.c * (1 / pi);
    Color emittance = obj.l;
    // llamada Recursiva y decremento de numero de rebotes
    Color incomingColor = shade(newRay, depth - 1);
    Color colorValue = emittance + (incomingColor.mult(BRDF) * cos_theta) * (1.0 / p);
    return colorValue;
}

// funcion principal
int main(int argc, char *argv[])
{
    int w = 1024, h = 768; // image resolution
    // Numero de muestras por pixel.
    const int pixel_samples = 32;
    // Numero de rebotes para introducir sesgo
    const int depth = 10;

    // fija la posicion de la camara y la dirección
    Ray camera(Point(0, 11.2, 214), Vector(0, -0.042612, -1).normalize());

    // parametros de la camara
    Vector cx = Vector(w * 0.5095 / h, 0., 0.);
    Vector cy = (cx % camera.d).normalize() * 0.5095;
    // auxiliar para valor de pixel y matriz para almacenar la imagen
    Color *pixelColors = new Color[w * h];

// se define bloque paralelo
#pragma omp parallel
    {
#pragma omp for schedule(dynamic, 1)
        for (int y = 0; y < h; y++)
        {
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