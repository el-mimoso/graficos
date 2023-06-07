// rt: un lanzador de rayos minimalista
// g++ -O3 -fopenmp rt.cpp -o rt
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <utility>
#include <bits/stdc++.h>

using namespace std;

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

            //ambos positivos
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

// Calcular el valor mas grande de un arreglo
void largest(double arr[], int n, double &min, double &max)
{
    int i;
    max = arr[0];
    min = arr[0];

    for (i = 1; i < n; i++)
    {
        if (arr[i] > max)
        {
            max = arr[i];
        }
        if (arr[i] < min)
        {
            min = arr[i];
        }
    }

    printf("\n");
    printf("El valor maximo es: %f\n", max);
    printf("El valor minimo es: %f\n", min);

    // return max, min;
}

double fineTuneT(double Max, double Min, double t)
{
    double tF = (t - Min) / (Max - Min);
    return tF;
}

// Calcula el valor de color para el rayo dado
Color shade(const Ray &r, int iterador, double ts[])
{
    double t;
    int id = 0;
    if (!intersect(r, t, id))
        // el rayo no intersecto objeto, return Vector() == negro
        return Color();

    const Sphere &obj = spheres[id];

    // PROYECTO 1
    // determinar coordenadas del punto de interseccion
    Point x = r.o + r.d * t;

    // determinar la dirección normal en el punto de interseccion
    Vector n = (x - obj.p).normalize();

    // determinar el color que se regresara
    Color colorValue;

    // obtener solo el color de los objetos
    //  colorValue = obj.c;

    // color de valores de normales en punto de interseccion
    // colorValue = Color(n.x + 1, n.y + 1, n.z + 1) * .5;

    // color(grises) de acuerdo a la profundidad

    colorValue = Color(t, t, t);
    ts[iterador] = t;
    return colorValue;
}

int main(int argc, char *argv[])
{
    int w = 1024, h = 768; // image resolution

    // vector <double> tArray; // para almacenar las distancias t
    double tArray[w * h];
    // tVector.reserve(w*h);

    // fija la posicion de la camara y la dirección en que mira
    Ray camera(Point(0, 11.2, 214), Vector(0, -0.042612, -1).normalize());

    // parametros de la camara
    Vector cx = Vector(w * 0.5095 / h, 0., 0.);
    Vector cy = (cx % camera.d).normalize() * 0.5095;

    // auxiliar para valor de pixel y matriz para almacenar la imagen
    Color *pixelColors = new Color[w * h];

// PROYECTO 1
// usar openmp para paralelizar el ciclo: cada hilo computara un renglon (ciclo interior),
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

                // computar el color del pixel para el punto que intersectó el rayo desde la camara
                pixelValue = shade(Ray(camera.o, cameraRayDir.normalize()), idx, tArray);

                pixelColors[idx] = pixelValue;
            }
        }
    }

    // tamanio de tArray
    int n = w * h;
    // mayor t
    double maxT, minT;
    largest(tArray, n, minT, maxT);

    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    std::cout << maxT << "\n ";
    std::cout << minT << "\n ";

    // Recorremos la imagen y normalizamos la nueva t.
    for (size_t i = 0; i < w * h; i++)
    {
        double newT = fineTuneT(maxT, minT, pixelColors[i].x);
        pixelColors[i] = Color(newT, newT, newT);
        pixelColors[i] = Color(clamp(pixelColors[i].x), clamp(pixelColors[i].y), clamp(pixelColors[i].z));
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