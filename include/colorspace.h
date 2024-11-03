#pragma once
#include "whitepoint.h"
namespace colorama {
    class RGBColorspace {
        public:
            const Whitepoint wp;
            double red_x;
            double red_y;
            double green_x;
            double green_y;
            double blue_x;
            double blue_y;

            RGBColorspace(const Whitepoint& wp, double rx, double ry, double gx, double gy, double bx, double by);
            static const RGBColorspace& sRGB();
            static const RGBColorspace& Adobe();
            static const RGBColorspace& Rec709();
            static const RGBColorspace& ACES();
            static const RGBColorspace& Rec2020();
            static const RGBColorspace& AppleDisplayP3();
    };
}