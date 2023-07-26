#ifndef TEXTURE_H
#define TEXTURE_H

#include "vec.h"
#include "utils.h"
#include "image.h"

#include <iostream>

class Texture
{
public:
    virtual Color value(double u, double v, const Vector &p) const = 0;
};

class SolidColor : public Texture
{
public:
    SolidColor() {}
    SolidColor(Color c) : color_value(c) {}

    SolidColor(double red, double green, double blue)
        : SolidColor(Color(red, green, blue)) {}

    virtual Color value(double u, double v, const Vector &p) const override
    {
        return color_value;
    }

private:
    Color color_value;
};

class CheckerTexture : public Texture
{
public: // variables
    Texture *even;
    Texture *odd;

public:
    CheckerTexture() {}
    CheckerTexture(Texture *_even, Texture *_odd) : even(_even), odd(_odd) {}

    virtual Color value(double u, double v, const Vector &p) const override
    {
        double sines = sin(10 * p.x) * sin(10 * p.y) * sin(10 * p.z);
        if (sines < 0)
            return odd->value(u, v, p);
        else
            return even->value(u, v, p);
    }
};

class ImageTexture : public Texture
{
private:
    unsigned char *data;
    int width, height;
    int bytes_per_scanline;

public:
    const static int bytes_per_pixel = 3;
    ImageTexture() : data(nullptr), width(0), height(0), bytes_per_scanline(0) {}
    ImageTexture(const char *filename)
    {
        auto components_per_pixel = bytes_per_pixel;

        data = stbi_load(
            filename, &width, &height, &components_per_pixel, components_per_pixel);

        if (!data)
        {
            std::cerr << "ERROR: Could not load texture image file '" << filename << "'.\n";
            width = height = 0;
        }

        bytes_per_scanline = bytes_per_pixel * width;
    }
    // ~imageTexture(){
    //     delete data;
    // }

    virtual Color value(double u, double v, const Vector &p) const override
    {
        // If we have no texture data, then return solid cyan as a debugging aid.
        if (data == nullptr)
            return Color(0, 1, 1);

        // Clamp input texture coordinates to [0,1] x [1,0]
        u = clamp(u, 0.0, 1.0);
        v = 1.0 - clamp(v, 0.0, 1.0); // Flip V to image coordinates

        auto i = static_cast<int>(u * width);
        auto j = static_cast<int>(v * height);

        // Clamp integer mapping, since actual coordinates should be less than 1.0
        if (i >= width)
            i = width - 1;
        if (j >= height)
            j = height - 1;

        const auto color_scale = 1.0 / 255.0;
        auto pixel = data + j * bytes_per_scanline + i * bytes_per_pixel;

        return Color(pixel[0] * color_scale, pixel[1] * color_scale, pixel[2] * color_scale);
    }
};
#endif