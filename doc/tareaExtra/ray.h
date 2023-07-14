#pragma once
#include "vec.h"

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