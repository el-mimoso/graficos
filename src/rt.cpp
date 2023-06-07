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

Sphere spheres[] = {
    // Escena: radio, posicion, color, luz
    Sphere(1e5, Point(-1e5 - 49, 0, 0), Color(.75, .25, .25), Color()),   // pared izq
    Sphere(1e5, Point(1e5 + 49, 0, 0), Color(.25, .25, .75), Color()),    // pared der
    Sphere(1e5, Point(0, 0, -1e5 - 81.6), Color(.25, .75, .25), Color()), // pared detras
    Sphere(1e5, Point(0, -1e5 - 40.8, 0), Color(.25, .75, .75), Color()), // suelo
    Sphere(1e5, Point(0, 1e5 + 40.8, 0), Color(.75, .75, .25), Color()),  // techo
    Sphere(16.5, Point(-23, -24.3, -34.6), Color(.2, .3, .4), Color()),   // esfera abajo-izq
    Sphere(16.5, Point(23, -24.3, -3.6), Color(.4, .3, .2), Color()),     // esfera abajo-der
    // Sphere(10.5, Point(0, 24.3, 0), Color(1, 1, 1), Color(10, 10, 10))    // esfera arriba
    Sphere(0, Point(0, 24.3, 0), Color(1, 1, 1), Color(5555, 5555, 5555)) // esfera arriba
    // sleek daft punk reference ヾ(⌐■_■)//ノ♪ 
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
        if (t > spheresData[i].second && spheresData[i].second > 0.0001)
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

// calcula el indice de la esfera que es luz
int getLightIndx()
{
    for (int i = 0; i < spheresLength; i++)
    {
        if (spheres[i].l.x > 0 || spheres[i].l.y > 0 || spheres[i].l.z > 0 && spheres[i].r == 0)
        {
            return i;
        }
    }
    return -1;
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

    int lightIndx = getLightIndx(); // indice de la esfera que es luz puntual
    double tlight;

    // distancia de la luz

    if (depth <= 0)
        return Color(0, 0, 0);

    if (!intersect(r, t, id))
        // el rayo no intersecto objeto, return Vector() == negro
        return Color(0, 0, 0);

    const Sphere &obj = spheres[id];

    const Sphere &light = spheres[lightIndx];

    // determinar coordenadas del punto de interseccion
    Point x = r.o + r.d * t;
    // printf("T  %g\n", t);

    // determinar la dirección normal en el punto de interseccion
    Vector n = (x - obj.p).normalize();

    // se lanza un rayo desde la fuente de luz hacia el punto de interseccion
    Ray newRay(light.p, x - light.p);

    Color emittance = light.l;

    double auxDistance = sqrt(newRay.d.dot(newRay.d));

    lightIndx = 0;
    // determinar si hay linea de vista entre la luz y la esfera
    if (intersect(newRay, tlight, lightIndx))
    {

        Color xprim = newRay.o + newRay.d * tlight;
        //revisamos que sea el mismo objeto con el que intersectamos. 
        if (id == lightIndx)
        {
        
            double xpNormSqr = xprim.dot(xprim);
            emittance = emittance * (1 / xpNormSqr);
        }
        else
            emittance = Color();
    }
    else
        emittance = Color();

    Vector lightD = newRay.d.normalize();
    double cos_theta = n.dot(lightD);

    Color BRDF = obj.c * (1 / pi);

    // Color incomingColor = shade(newRay, depth - 1, mode);

    Color colorValue = emittance.mult(BRDF * (-cos_theta));
    return colorValue;
}

// funcion principal
int main(int argc, char *argv[])
{
    int w = 1024, h = 768; // image resolution
    // Numero de muestras por pixel.
    const int pixel_samples = 32;
    // Numero de rebotes.
    const int depth = 2;
    // Modo de muestreo. 0 para esfera unitaria, 1 para esfera hemisferica, 2 para coseno hemisferica
    const int mode = 0;

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
                // ciclo de muestreo para antiAlias.
                for (int i = 0; i < pixel_samples; i++)
                {
                    // auto u = cx * (double(x) / w - .5);
                    // auto v = cy * (double(y) / h - .5);
                    // se le agrega un valor aleatorio a cada pixel en X Y para muestrear mas puntos, no solo el centro del pixel
                    auto u = cx * (double(x + random_double()) / w - 0.5);
                    auto v = cy * (double(y + random_double()) / h - 0.5);
                    // se lanza rayo con las variaciones en direccion
                    Vector cameraRayDir = u + v + camera.d;
                    // incrementamos valores de ese pixel

                    pixelValue = pixelValue + shade(Ray(camera.o, cameraRayDir.normalize()), depth, mode);
                }
                // promedio y escalado del valor del pixel
                auto scale = 1.0 / pixel_samples;
                pixelValue = pixelValue * scale;

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