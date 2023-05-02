// rt: un lanzador de rayos minimalista
// g++ -O3 -fopenmp rt.cpp -o rt
#include "vec.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <utility>
using namespace std;

// TODO: pasar las clases a su propio archivo.

class Ray
{
public:
    Point o;
    Vector d;                                  // origen y direcccion del rayo
    Ray(Point o_, Vector d_) : o(o_), d(d_) {} // constructor

    Point at(double t) const
    {
        return o + d * t;
    }
};

class Sphere
{
public:
    double r; // radio de la esfera
    Point p;  // posicion
    Color c;  // color
    Color l;  // luz

    Sphere(double r_, Point p_, Color c_, Color l_) : r(r_), p(p_), c(c_), l(l_) {}

    // determina si el rayo intersecta a esta esfera
    // -1 si no toca, t de lo contrario.
    double intersect(const Ray &ray) const
    {
        // se asignan los valores de la formula general cuadratica "la del chicharronero"
        auto oc = ray.o - p;
        auto a = ray.d.dot(ray.d);
        auto b = oc.dot(ray.d);
        auto c = oc.dot(oc) - r * r;
        auto discriminant = b * b - a * c;
        // solo si el valor de discriminant es positivo hacemos la raiz cuadrada para ahorrar tiempo de computo
        // de otra forma retorna -1
        if (discriminant < 0)
        {
            return -1.0;
        }
        else
        {
            return (-b - sqrt(discriminant)) / a;
        }
    }
};

Sphere spheres[] = {
    // Escena: radio, posicion, color, luz
    Sphere(1e5, Point(-1e5 - 49, 0, 0), Color(.75, .25, .25), Color()),   // pared izq
    Sphere(1e5, Point(1e5 + 49, 0, 0), Color(.25, .25, .75), Color()),    // pared der
    Sphere(1e5, Point(0, 0, -1e5 - 81.6), Color(.25, .75, .25), Color()), // pared detras
    Sphere(1e5, Point(0, -1e5 - 40.8, 0), Color(.25, .75, .75), Color()), // suelo
    Sphere(1e5, Point(0, 1e5 + 40.8, 0), Color(.75, .75, .25), Color()),  // techo
    Sphere(16.5, Point(-23, -24.3, -34.6), Color(.2, .3, .4), Color()),   // esfera abajo-izq
    Sphere(16.5, Point(23, -24.3, -3.6), Color(.4, .3, .2), Color()),     // esfera abajo-der
    Sphere(10.5, Point(0, 24.3, 0), Color(1, 1, 1), Color(10, 10, 10))    // esfera arriba
    // Sphere(0, Point(0, 24.3, 0), Color(1, 1, 1), Color(10, 10, 10)) // esfera arriba
};
const int spheresLength = sizeof(spheres) / sizeof(spheres[0]);

// TODO: pasar las siguientes funciones utiles a su propio archivo.
// limita el valor de x a [0,1]
inline double
clamp(const double x)
{
    if (x < 0.0)
        return 0.0;
    else if (x > 1.0)
        return 1.0;
    return x;
}

// convierte un valor de color en [0,1] a un entero en [0,255]
inline int toDisplayValue(const double x)
{
    return int(pow(clamp(x), 1.0 / 2.2) * 255 + .5);
}

Vector cross(Vector a, Vector b)
{
    return a % b;
}

void coordinateSystem(const Vector &n, Vector &s, Vector &t)
{
    if (std::abs(n.x) > std::abs(n.y))
    {
        float invLen = 1.0f / std::sqrt(n.x * n.x + n.z * n.z);
        t = Vector(n.z * invLen, 0.0f, -n.x * invLen);
    }
    else
    {
        float invLen = 1.0f / std::sqrt(n.y * n.y + n.z * n.z);
        t = Vector(0.0f, n.z * invLen, -n.y * invLen);
    }
    s = cross(t, n);
}

// Returna real aleatorio en el rango  [0,1).
inline double random_double()
{
    return rand() / (RAND_MAX + 1.0);
}

// Returna real aleatorio en el rango [min,max).
inline double random_double(double min, double max)
{
    return min + (max - min) * random_double();
}

Vector random_in_unit_sphere()
{
    while (true)
    {
        Vector p = Vector(random_double(-1, 1), random_double(-1, 1), random_double(-1, 1));
        if (p.dot(p) >= 1)
            continue;
        return p;
    }
}
Vector random_unit_vector()
{
    return random_in_unit_sphere().normalize();
}

Vector random_in_hemisphere(const Vector &normal)
{
    Vector in_unit_sphere = random_in_unit_sphere();
    if (in_unit_sphere.dot(normal) > 0.0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return in_unit_sphere * -1;
}

const double pi = 3.1415926535897932385;

// calcula la intersección del rayo r con todas las esferas
// regresar true si hubo una intersección, falso de otro modo
// almacenar en t la distancia sobre el rayo en que sucede la interseccion
// almacenar en id el indice de spheres[] de la esfera cuya interseccion es mas cercana
inline bool intersect(const Ray &r, double &t, int &id)
{
    // arreglo de pares donde almacenamos la id y la distancia del rayo
    pair<int, int> spheresData[spheresLength];

    for (int i = 0; i < spheresLength; i++)
    {
        // asignación de valores al par
        spheresData[i].first = i;
        spheresData[i].second = spheres[i].intersect(r);
        if (spheresData[i].second > 0)
        {
            t = spheres[i].intersect(r);
            id = i;
        }
    }
    // ordenamiento de pares por la menor distancia
    for (int i = 0; i < spheresLength; i++)
    {
        if (spheresData[i].second < t && spheresData[i].second > 0)
        {
            // actualizamos valores de t e id solo si el valor almacenado en el par es menor al valor de t y que sea positivo
            t = spheresData[i].second;
            id = spheresData[i].first;
        }
    }
    if (t > 0)
    {
        return true;
    }
    return false;
}

// Calcula el valor de color para el rayo dado
Color shade(const Ray &r, int depth)
{
    double t;
    int id = 0;

    if (depth <= 0)
        return Color(0, 0, 0);

    if (!intersect(r, t, id))
        // el rayo no intersecto objeto, return Vector() == negro
        return Color();

    const Sphere &obj = spheres[id];

    if (obj.l.x * obj.l.x + obj.l.y * obj.l.y + obj.l.z * obj.l.z > 0)
    {
        return obj.l;
    }

    // determinar coordenadas del punto de interseccion
    Point x = r.at(t);

    // determinar la dirección normal en el punto de interseccion
    Vector n = (x - obj.p).normalize();

    // const double p = .25 / pi;
    const double p = .5 / pi;

    // Point target = x + n + random_unit_vector();
    Point target = x + random_in_hemisphere(n);
    Ray newRay(x, target - x);

    double cos_theta = newRay.d.dot(n);

    Color BRDF = obj.c * (1 / pi);
    Color emittance = obj.l;

    Color incomingColor = shade(newRay, depth - 1);

    Color colorValue = emittance + (incomingColor.mult(BRDF) * cos_theta) * (1 / p);
    return colorValue;
    // return Color(n.x + 1, n.y + 1, n.z + 1) * .5;
}

int main(int argc, char *argv[])
{
    int w = 1024, h = 768; // image resolution
    // Numero de muestras por pixel.
    const int pixel_samples = 32;
    const int depth = 32;

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
                    auto u = cx * (double(x) / w - .5);
                    auto v = cy * (double(y) / h - .5);
                    // se le agrega un valor aleatorio a cada pixel en X Y para muestrear mas puntos, no solo el centro del pixel
                    // auto u = cx * (double(x + random_double()) / w - 0.5);
                    // auto v = cy * (double(y + random_double()) / h - 0.5);
                    // se lanza rayo con las variaciones en direccion
                    Vector cameraRayDir = u + v + camera.d;
                    // incrementamos valores de ese pixel
                    pixelValue = pixelValue + shade(Ray(camera.o, cameraRayDir.normalize()), depth);
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