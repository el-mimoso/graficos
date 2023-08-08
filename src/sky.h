#ifndef SKY_H
#define SKY_H

#include "vec.h"
#include "ray.h"
#include "utils.h"
#include "image.h"
#include <iostream>

class Sky
{
public:
    virtual Color value(const Ray ray) const = 0;
};

class AmbientMap : public Sky
{
private:
    float *data;
    int width, height;
    int bytes_per_scanline;

public:
    // constructor con datos vacios
    AmbientMap() : data(nullptr), width(0), height(0), bytes_per_scanline(0) {}
    // constructor con valores del skybox
    AmbientMap(const char *filename)
    {
        int n;
        // carga la imagen con valores flotantes
        data = stbi_loadf(
            filename, &width, &height, &n, 0);
        if (!data)
        {
            std::cerr << "ERROR: Could not load skybox HDRI '" << filename << "'.\n";
            width = height = 0;
        }
    }
    // funcion que retorna el color del skybox
    virtual Color value(const Ray ray) const override
    {
        // se obtiene la direccion del rayo
        Vector unit_direction = ray.d;
        unit_direction.normalize();
        // se obtiene el angulo theta
        // double theta = acos(-unit_direction.y);
        double theta = acos(unit_direction.y);
        // se obtiene el angulo phi
        // double phi = atan2(-unit_direction.z, unit_direction.x) + pi;
        double phi = atan2(unit_direction.z, unit_direction.x) + pi;
        // se obtiene la posicion en x de la imagen
        double u = phi / (2 * pi);
        // se obtiene la posicion en y de la imagen
        double v = (theta / pi) ;
        // se obtiene la posicion en la imagen
        int i = u * width;
        int j = v * height;
        // si la posicion es mayor que el tamaño de la imagen
        // if (i >= width)
        //     i = width - 1;
        // if (j >= height)
        //     j = height - 1;
        int aux = 3 * i + 3 * width * j;

        // se obtiene el color de la imagen
        float r = data[aux];
        float g = data[aux + 1];
        float b = data[aux + 2];
        // se retorna el color
        return Color(r, g, b);
    }

    ~AmbientMap()
    {
        stbi_image_free(data);
    }
};

#endif