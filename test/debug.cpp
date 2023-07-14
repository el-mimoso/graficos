// rt: un lanzador de rayos minimalista
// g++ -O3 -fopenmp rt.cpp -o rt
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <stdio.h>
#include <time.h> // for clock_t, clock(), CLOCKS_PER_SEC
#define NUM_THREADS 8
#include <cstdlib>

auto pi = 3.1416;

inline double random_double()
{
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max)
{
    // Returns a random real in [min,max).
    return min + (max - min) * random_double();
}

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

inline Vector cross(const Vector &u, const Vector &v)
{
    return Vector(u.y * v.z - u.z * v.y,
                  u.z * v.x - u.x * v.z,
                  u.x * v.y - u.y * v.x);
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
inline Vector random_esfera()
{ // Esta funcion crea direcciones aleatorias esfericas.
    auto r1 = random_double();
    auto r2 = random_double();

    auto x = cos(2 * pi * r1) * 2 * sqrt(r2 * (1 - r2));
    auto y = sin(2 * pi * r1) * 2 * sqrt(r2 * (1 - r2));
    auto z = 1 - 2 * r2;

    return Vector(x, y, z);
}

inline Vector random_hemisferio()
{ // Esta funcion crea direcciones aleatorias en un hemisferio.
    auto r1 = random_double();
    auto r2 = random_double();

    auto theta = acos(1 - r1);
    auto phi = 2 * pi * r2;

    auto x = sin(theta) * cos(phi);
    auto y = sin(theta) * sin(phi);
    auto z = 1 - r1;

    return Vector(x, y, z);
}

inline Vector random_coseno(double &theta)
{ // Esta funcion crea direcciones aleatorias con distribucion coseno hemisferico.
    auto r1 = random_double();
    auto r2 = random_double();

    auto phi = 2 * pi * r1;

    auto z = sqrt(1 - r2);
    auto x = cos(phi) * sqrt(r2);
    auto y = sin(phi) * sqrt(r2);

    return Vector(x, y, z);
}

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
    Color l;

    Sphere(double r_, Point p_, Color c_, Color m_) : r(r_), p(p_), c(c_), l(m_) {}

    // PROYECTO 1
    // determina si el rayo intersecta a esta esfera
    double intersect(const Ray &ray) const
    {
        Vector oc = ray.o - p;

        double a = ray.d.dot(ray.d);
        double b = oc.dot(ray.d);
        double c = oc.dot(oc) - r * r;
        double discriminant = b * b - a * c;
        // regresar distancia si hay intersección
        // regresar 0.0 si no hay interseccion
        if (discriminant < 0)
        {
            return 0.0;
        }
        else
        {
            return (-b - sqrt(discriminant) / a);
        }
    }
};

Sphere spheres[] = {
    // Escena: radio, posicion, color
    Sphere(1e5, Point(-1e5 - 49, 0, 0), Color(.75, .25, .25), Color()),   // pared izq
    Sphere(1e5, Point(1e5 + 49, 0, 0), Color(.25, .25, .75), Color()),    // pared der
    Sphere(1e5, Point(0, 0, -1e5 - 81.6), Color(.25, .75, .25), Color()), // pared detras
    Sphere(1e5, Point(0, -1e5 - 40.8, 0), Color(.25, .75, .75), Color()), // suelo
    Sphere(1e5, Point(0, 1e5 + 40.8, 0), Color(.75, .75, .25), Color()),  // techo
    Sphere(16.5, Point(-23, -24.3, -34.6), Color(.2, .3, .4), Color()),   // esfera abajo-izq
    Sphere(16.5, Point(23, -24.3, -3.6), Color(.4, .3, .2), Color()),     // esfera abajo-der
    Sphere(10.5, Point(0, 24.3, 0), Color(1, 1, 1), Color(10, 10, 10))    // esfera arriba
};

// limita el valor de x a [0,1]
inline double clamp(const double x)
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

// PROYECTO 1
// calcular la intersección del rayo r con todas las esferas
// regresar true si hubo una intersección, falso de otro modo
// almacenar en t la distancia sobre el rayo en que sucede la interseccion
// almacenar en id el indice de spheres[] de la esfera cuya interseccion es mas cercana
inline bool intersect(const Ray &r, double &t, int &id)
{
    int NS = 8;
    double aux[NS];

    for (int i = 0; i < NS; i++)
    {
        aux[i] = spheres[i].intersect(r);
        if (aux[i] > 0)
        {
            t = aux[i];
            id = i;
        };
    };

    for (int i = 0; i < NS; i++)
    {
        if (t > aux[i] && aux[i] > 0.01)
        {
            t = aux[i];
            id = i;
        };
    };
    if (t > 0)
    {
        // printf("%f  %d\t ",t,id);
        return true;
    };

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

    // if (obj.l.x * obj.l.x + obj.l.y * obj.l.y + obj.l.z * obj.l.z > 0)
    // {
    //     return obj.l;
    // }

    // determinar coordenadas del punto de interseccion
    Point x = r.o + r.d * t;

    // determinar la dirección normal en el punto de interseccion
    Vector n = (x - obj.p).normalize();

    const double p = 1 / (4 * pi);
    // const double p = .5 /  pi;
    // double tetha;
    // Point target = random_coseno(tetha);
    // const double p = (1 / pi) * cos(tetha);
    Vector s;
    Vector ta;
    coordinateSystem(n, s, ta);

    Point target =  random_esfera();
    // Point target = random_hemisferio();

    Point dir(Point(s.x, ta.x, n.x).dot(target), Point(s.y, ta.y, n.y).dot(target), Point(s.z, ta.z, n.z).dot(target));

    Ray newRay(x, dir);

    double cos_theta = newRay.d.dot(n);

    Color BRDF = obj.c * (1 / pi);
    Color emittance = obj.l;

    Color incomingColor = shade(newRay, depth - 1);

    Color colorValue = emittance + (incomingColor.mult(BRDF) * cos_theta) * (p);
    return colorValue;
    // return Color(n.x + 1, n.y + 1, n.z + 1) * .5;
}

int main(int argc, char *argv[])
{
    double time_spent = 0.0;
    int muestras = 64;

    clock_t begin = clock();
    // sleep(3);

    int w = 1024, h = 768; // image resolution

    // fija la posicion de la camara y la dirección en que mira
    Ray camera(Point(0, 11.2, 214), Vector(0, -0.042612, -1).normalize());

    // parametros de la camara
    Vector cx = Vector(w * 0.5095 / h, 0., 0.);
    Vector cy = (cx % camera.d).normalize() * 0.5095;

    // auxiliar para valor de pixel y matriz para almacenar la imagen
    Color *pixelColors = new Color[w * h];

    // PROYECTO 1
    // usar openmp para paralelizar el ciclo: cada hilo computara un renglon (ciclo interior),

    omp_set_num_threads(NUM_THREADS); // Lineas de Codigo para paralelizar
#pragma omp parallel for              // schedule(dynamic,1) //Con la paralelizacion se reduce en un 60 % aproxiamadamente el tiempo de ejecucion.

    for (int y = 0; y < h; y++)
    {
        // recorre todos los pixeles de la imagen
        fprintf(stderr, "\r%5.2f%%", 100. * y / (h - 1));
        for (int x = 0; x < w; x++)
        {

            int idx = (h - y - 1) * w + x; // index en 1D para una imagen 2D x,y son invertidos

            Color pixelValue = Color(); // pixelValue en negro por ahora
            // para el pixel actual, computar la dirección que un rayo debe tener

            Vector cameraRayDir = cx * (double(x) / w - .5) + cy * (double(y) / h - .5) + camera.d;
            // computar el color del pixel para el punto que intersectó el rayo desde la camara
            for (int k = 0; k < muestras; k++)
            {
                pixelValue = pixelValue + shade(Ray(camera.o, cameraRayDir.normalize()), 6);
            }
            // limitar los tres valores de color del pixel a [0,1]
            pixelColors[idx] = Color(clamp(pixelValue.x), clamp(pixelValue.y), clamp(pixelValue.z));
        }
    }

    fprintf(stderr, "\n");
    clock_t end = clock();

    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;

    printf("The elapsed time is %f seconds", time_spent);

    // PROYECTO 1
    // Investigar formato ppm
    FILE *f = fopen("image.ppm", "w");
    // escribe cabecera del archivo ppm, ancho, alto y valor maximo de color
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int p = 0; p < w * h; p++)
    { // escribe todos los valores de los pixeles
        fprintf(f, "%d %d %d \n", toDisplayValue(pixelColors[p].x), toDisplayValue(pixelColors[p].y),
                toDisplayValue(pixelColors[p].z));
    }
    fclose(f);

    delete[] pixelColors;

    return 0;
}
