// rt: un lanzador de rayos minimalista
// g++ -O3 -fopenmp rt.cpp -o rt
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <utility>
using namespace std;

// TODO: pasar las clases a su propio archivo.

class Vector
{
public:
    double x, y, z; // coordenadas x,y,z

    // Constructor del vector, parametros por default en cero
    Vector(double x_ = 0, double y_ = 0, double z_ = 0)
    {
        x = x_;
        y = y_;
        z = z_;
    }

    // operador para suma y resta de vectores
    Vector operator+(const Vector &b) const { return Vector(x + b.x, y + b.y, z + b.z); }
    Vector operator-(const Vector &b) const { return Vector(x - b.x, y - b.y, z - b.z); }
    // operator multiplicacion vector y escalar
    Vector operator*(double b) const { return Vector(x * b, y * b, z * b); }

    // operator % para producto cruz
    Vector operator%(Vector &b) { return Vector(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }

    // producto punto con vector b
    double dot(const Vector &b) const { return x * b.x + y * b.y + z * b.z; }

    // producto elemento a elemento (Hadamard product)
    Vector mult(const Vector &b) const { return Vector(x * b.x, y * b.y, z * b.z); }

    // normalizar vector
    Vector &normalize() { return *this = *this * (1.0 / sqrt(x * x + y * y + z * z)); }
};
typedef Vector Point;
typedef Vector Color;

class Ray
{
public:
    Point o;
    Vector d;                                  // origen y direcccion del rayo
    Ray(Point o_, Vector d_) : o(o_), d(d_) {} // constructor
};

class Sphere
{
public:
    double r; // radio de la esfera
    Point p;  // posicion
    Color c;  // color

    Sphere(double r_, Point p_, Color c_) : r(r_), p(p_), c(c_) {}

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
            double tplus = (-b + sqrt(discriminant)) / a;
            double tminus = (-b - sqrt(discriminant)) / a;
            double t;

            // ambos positivos
            if (tminus > 0 && tplus > 0)
            {
                t = min(tminus, tplus);
            }
            // tminus positivo, tplus negativo
            else if (tminus > 0 && tplus < 0)
            {
                t = tminus;
            }
            // tminus negativo, tplus positivo
            else if (tminus < 0 && tplus > 0)
            {
                t = tplus;
            }
            else
            {
                t = 0;
            }
            return t;
        }
    }
};

Sphere spheres[] = {
    // Escena: radio, posicion, color
    Sphere(1e5, Point(-1e5 - 49, 0, 0), Color(.75, .25, .25)),       // pared izq
    Sphere(1e5, Point(1e5 + 49, 0, 0), Color(.25, .25, .75)),        // pared der
    Sphere(1e5, Point(0, 0, -1e5 - 81.6), Color(.75, .75, .75)),     // pared detras
    Sphere(1e5, Point(0, -1e5 - 40.8, 0), Color(.75, .75, .75)),     // suelo
    Sphere(1e5, Point(0, 1e5 + 40.8, 0), Color(.75, .75, .75)),      // techo
    Sphere(16.5, Point(-23, -24.3, -34.6), Color(.999, .999, .999)), // esfera abajo-izq
    Sphere(16.5, Point(23, -24.3, -3.6), Color(.999, .999, .999)),   // esfera abajo-der
    Sphere(10.5, Point(0, 24.3, 0), Color(1, 1, 1))                  // esfera arriba
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

// Retorna real aleatorio en el rango  [0,1).
inline double random_double()
{
    return drand48();
}

// Retorna real aleatorio en el rango [min,max).
inline double random_double(double min, double max)
{
    return min + (max - min) * random_double();
}

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
Color shade(const Ray &r)
{
    double t;
    int id = 0;
    if (!intersect(r, t, id))
        // el rayo no intersecto objeto, return Vector() == negro
        return Color();

    const Sphere &obj = spheres[id];

    // determinar coordenadas del punto de interseccion
    Point x = r.o + r.d * t;

    // determinar la dirección normal en el punto de interseccion
    Vector n = (x - obj.p).normalize();

    // determinar el color que se regresara
    Color colorValue;

    // color de valores de normales en punto de interseccion
    colorValue = Color(n.x + 1, n.y + 1, n.z + 1) * .5;

    return colorValue;
}

int main(int argc, char *argv[])
{
    int w = 1024, h = 768; // image resolution
    // Numero de muestras por pixel.
    const int pixel_samples = 32;

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

                // Vector cameraRayDir = cx * (double(x) / w - .5) + cy * (double(y) / h - .5) + camera.d;
                // pixelValue = shade(Ray(camera.o, cameraRayDir.normalize()));

                // ciclo de muestreo para antiAlias.
                for (int i = 0; i < pixel_samples - 1; i++)
                {
                    // se le agrega un valor aleatorio a cada pixel en X Y para muestrear mas puntos, no solo el centro del pixel
                    auto u = cx * (double(x + random_double(-.5, .5)) / w - 0.5);
                    auto v = cy * (double(y + random_double(-.5, .5)) / h - 0.5);

                    // se lanza rayo con las variaciones en direccion
                    Vector cameraRayDir = u + v + camera.d;
                    // incrementamos valores de ese pixel
                    pixelValue = pixelValue + shade(Ray(camera.o, cameraRayDir.normalize()));
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