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
#include <vector>
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
    Sphere(16.5, Point(23, -24.3, -3.6), Color(.4, .3, .2), Color()), // esfera abajo-der
    Sphere(10.5, Point(0, 24.3, 0), Color(1, 1, 1), Color(10,10,10)),
};
// almacenamos las luces.
vector<int> lights;
const int spheresLength = sizeof(spheres) / sizeof(spheres[0]);
// int lightsLength [] = arraySize(spheres);

int getLightsIndx()
{
    int lightsCount = 0;
    for (int i = 0; i < spheresLength; i++)
    {
        if (spheres[i].l.x > 0.001 || spheres[i].l.y > 0.001 || spheres[i].l.z > 0.001)
        {
            // lightsPos[i].first = i;
            lights.push_back(i);
        }
    }
    return lights.size();
    // cout << areaProb(0);
}
const int totalLights = getLightsIndx();

// Calcula el valor de color para el rayo dado, seleccionar muestreo aleatorio
// 0 para muestreo de area
// 1 para angulo solido
// regresar el color resultante
Color shade(const Ray &r, int mode)
{
    double t;
    int id = 0;
    // seleccionar una luz aleatoria entre las luces de la escena
    int thisLight = randint(0, totalLights - 1);
    thisLight = lights[thisLight];

    if (!intersect(r, t, id, spheres, spheresLength))
        // el rayo no intersecto objeto, return Vector() == negro
        return Color(0, 0, 0);

    const Sphere &obj = spheres[id];

    // determinar coordenadas del punto de interseccion
    Point x = r.o + r.d * t;

    // determinar la dirección normal en el punto de interseccion
    Vector n = (x - obj.p).normalize();

    // Point target;
    double p;
    double theta;
    double phi;
    Vector s;
    Vector ta;
    Color emittance;
    double cosTheta;
    Color colorValue, Le, BRDF;
    Vector wc, wi;

    Sphere lightSrc = spheres[thisLight];

    if (id == thisLight)
    {
        colorValue = Color(0, 0, 0);
    }
    else
    {
        // para muestreo de area
        if (mode == 0)
        {
            // double theta, phi;
            wc = (lightSrc.p - x).normalize();
            // actualizamos valores de teta y phi por referencia
            areaPhiTheta(phi, theta);
            // vector unitario d
            Vector d = createVec(theta, phi, lightSrc).normalize();
            Point x1 = lightSrc.p + d * lightSrc.r;
            wi = x1 - x;
            Le = areaEval(x, wi, 7, lightSrc, spheres, spheresLength);
            BRDF = obj.c * (1.0 / pi);
            cosTheta = n.dot(wi.normalize());
            p = areaProb(x, x1, lightSrc);
            colorValue = Le.mult(BRDF) * (cosTheta / p);
        }
        // para muestreo de angulo solido
        else if (mode == 1)
        {
            // // armamos sistema de coordenadas
            wc = (lightSrc.p - x).normalize();
            coordinateSystem(wc, s, ta);
            // // actualizamos valores de teta, phi y Cos theta Max por referencia
            double cosThetaMax;
            solidAnglePhiTheta(phi, theta, x, lightSrc, cosThetaMax);
            wi = createVec(theta, phi).normalize();
            Vector wiGlobal = makeGlobal(wi, wc, s, ta);
            wiGlobal.normalize();
            Le = solidAngleEval(x, wiGlobal, lightSrc, spheres, spheresLength);
            BRDF = obj.c * (1 / pi);
            cosTheta = n.dot(wiGlobal);
            p = solidAngleProb(cosThetaMax);
            colorValue = colorValue + Le.mult(BRDF) * (cosTheta / p);
        }
    }
    return (obj.l + colorValue);
}

// funcion principal
int main(int argc, char *argv[])
{
    int w = 1024, h = 768; // image resolution
    // Numero de muestras por pixel.
    const int pixel_samples = 32;
    // Modo de muestreo. 0 para area, 1 para angulo solido.
    const int mode = 1;
    cout << "Total lights : " << totalLights << endl;

    for (int i = 0; i < totalLights; i++)
    {
        cout << "indices de luces: " <<lights[i] << endl;
    }
    
    // getLightsIndx();

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
            // ciclo de muestreo
            for (int x = 0; x < w; x++)
            {
                int idx = (h - y - 1) * w + x; // index en 1D para una imagen 2D x,y son invertidos
                Color pixelValue = Color();    // pixelValue en negro por ahora
                // para el pixel actual, computar la dirección que un rayo debe tener
                Vector cameraRayDir = cx * (double(x) / w - .5) + cy * (double(y) / h - .5) + camera.d;
                for (int i = 0; i < pixel_samples; i++)
                {
                    // computar el color del pixel para el punto que intersectó el rayo desde la camara
                    pixelValue = pixelValue + shade(Ray(camera.o, cameraRayDir.normalize()), mode);
                }
                // se divide por el numero de muestras.
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