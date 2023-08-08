// rt: un lanzador de rayos minimalista
// g++ -O3 -fopenmp rt.cpp -o rt
#include "vec.h"
#include "ray.h"
#include "sphere.h"
#include "utils.h"
#include "texture.h"
#include "sky.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <utility>
#include <typeinfo>

using namespace std;

// colores y Texturas
SolidColor red(Color(0.75, 0.25, 0.25));
SolidColor white(Color(0.75, 0.75, 0.75));
SolidColor green(Color(0.25, 0.75, 0.25));
SolidColor blue(Color(0.25, 0.25, 0.75));

CheckerTexture checker(&white, &blue);

ImageTexture earth("assets/earthmap.jpg");

ImageTexture golfNormals("assets/golfNormals.jpg");

ImageTexture rockAlbedo("assets/rock_albedo.jpeg");
ImageTexture rockNormals("assets/rockNormals.jpeg");

ImageTexture brickAlbedo("assets/brick_albedo.jpeg");
ImageTexture brickNormals("assets/brickNormals.jpeg");

AmbientMap skybox("assets/grass.hdr");
// AmbientMap skybox("assets/cafe.hdr");
// materiales para la escena

DifusseOG m1(Color(0.75, 0.25, 0.25), Color()); // pared izq
DifusseOG m2(Color(0.25, 0.25, 0.75), Color()); // pared der
DifusseOG m3(Color(0.25, 0.75, 0.25), Color()); // pared detras
// DifusseOG m4(Color(0.25, 0.75, 0.75), Color()); // suelo
DifusseTx m4(&white);                           // suelo
DifusseOG m5(Color(0.75, 0.75, 0.25), Color()); // techo
DifusseTx m6(&rockAlbedo, &rockNormals);        // esfera abajo-izq
// DifusseTx m7(&white, &golfNormals);           // esfera abajo-der
// Specular m7(Color(0.999, 0.999, 0.999), Color()); // esfera abajo-der
Metal m7(Color(0.143245, 0.377423, 1.43919), Color(3.98479, 2.3847, 1.60434), 0.3); // esfera abajo-der
// DifusseTx m7(&brickAlbedo, &brickNormals);                // esfera abajo-der
DifusseOG m8(Color(1.00, 1.00, 1.00), Color(10, 10, 10)); // esfera arriba LUZ

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
    hitInfo info;

    if (depth <= 0)
        return Color(0, 0, 0);

    if (!intersect(r, t, id))
        // el rayo no intersecto objeto, return Vector() == negro
        // return skybox.value(r);
        return Color(0, 0, 0);

    const Sphere &obj = spheres[id];

    // determinar coordenadas del punto de interseccion
    Point x = r.o + r.d * t;
    info.x = x;

    // determinar la dirección normal en el punto de interseccion
    Vector n = (x - obj.p).normalize();
    info.N = n;

    // actualizar las coordenadas u v
    info.nv = checkNormalDir(n, r);
    obj.get_UV(info.N, info.u, info.v);

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

    wo = material->sampling(r, info);
    p = material->probability(r, wo, info);
    // direccion saliente
    Ray newRay(x, wo);
    double cos_theta = newRay.d.dot(n);
    Color BRDF = material->eval_f(r, newRay, info);

    // si el material es especular perfecto
    if (typeid(*material) == typeid(Specular))
    {
        Vector wr;
        wr = material->sampling(r, info);
        p = material->probability(r, wo, info);
        newRay = Ray(x, wr);
        BRDF = material->eval_f(r, newRay, info);
        cos_theta = newRay.d.dot(n);
        return BRDF.mult(shade(Ray(x, wr), depth - 1));
    }
    // si el material es dielectrico
    if (typeid(*material) == typeid(Glass))
    {
        Vector wrt;
        // Se mandan los valores de n y nv para que se haga la reflexion interna
        wrt = material->sampling(r, info);
        p = material->probability(r, wo, info);
        // direccion del nuevo rayo
        newRay = Ray(x, wrt);
        // cout<<wrt.x<<" "<<wrt.y<<" "<<wrt.z<<endl;
        // evaluacion de la BRDF
        BRDF = material->eval_f(r, newRay, info);
        return BRDF.mult(shade(newRay, depth - 1));
    }
    // if (typeid(*material) == typeid(Metal))
    // {
    //     wo = material->sampling(r, info);
    //     // cout<<wo.x<<" "<<wo.y<<" "<<wo.z<<endl;
    //     p = material->probability(r, wo, info);
    //     Ray newRay(x, wo);
    //     double cos_theta = newRay.d.dot(n);
    //     Color BRDF = material->eval_f(r, newRay, info);
    //     // cout<<p<<endl;
    //     return BRDF.mult(shade(newRay, depth - 1)) * ( p);
    // }

    Color incomingColor = shade(newRay, depth - 1);
    Color colorvalue = (incomingColor.mult(BRDF) * (fabs(cos_theta) / p));
    return colorvalue;
}

// funcion principal
int main(int argc, char *argv[])
{
    int w = 1024, h = 768; // image resolution
    // Numero de muestras por pixel.
    const int pixel_samples = 32;
    // Numero de rebotes.
    const int depth = 3;
    // super sampling
    const int ssa = 4;

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
                    // super sampling
                    for (int j = 0; j < ssa; j++)
                    {
                        auto u = cx * (double(x + random_double(-.5, .5)) / w - 0.5);
                        auto v = cy * (double(y + random_double(-.5, .5)) / h - 0.5);
                        // se lanza rayo con las variaciones en direccion
                        Vector cameraRayDir = u + v + camera.d;
                        // incrementamos valores de ese pixel
                        pixelValue = pixelValue + shade(Ray(camera.o, cameraRayDir.normalize()), depth);
                    }
                    // computar el color del pixel para el punto que intersectó el rayo desde la camara
                    // pixelValue = pixelValue + shade(Ray(camera.o, cameraRayDir.normalize()), depth);
                }
                pixelValue = pixelValue * (1.0 / (pixel_samples * ssa));
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