#include "colorspace.h"

namespace colorama {
    RGBColorspace::RGBColorspace(const Whitepoint& wp, double rx, double ry, double gx, double gy, double bx, double by)
            : wp(wp), red_x(rx), red_y(ry), green_x(gx), green_y(gy), blue_x(bx), blue_y(by) {}

    const RGBColorspace& RGBColorspace::sRGB() {
        static const RGBColorspace srgb(Whitepoint::D65(), 0.640, 0.330, 0.300, 0.600, 0.150, 0.060);
        return srgb;
    }

    const RGBColorspace& RGBColorspace::Adobe() {
        static const RGBColorspace adobe(Whitepoint::D65(), 0.640, 0.330, 0.210, 0.710, 0.150, 0.060);
        return adobe;
    }

    const RGBColorspace& RGBColorspace::Rec709() {
        static const RGBColorspace rec709(Whitepoint::D65(), 0.640, 0.330, 0.300, 0.600, 0.150, 0.060);
        return rec709;
    }
    const RGBColorspace& RGBColorspace::ACES() {
        static const RGBColorspace aces(Whitepoint::D65(), 0.7347, 0.2653, 0.0000, 1.0000, 0.0001, -0.0774);
        return aces;
    }

    const RGBColorspace& RGBColorspace::Rec2020() {
        static const RGBColorspace rec2020(Whitepoint::D65(), 0.708, 0.292, 0.170, 0.797, 0.131, 0.046);
        return rec2020;
    }
    const RGBColorspace& RGBColorspace::AppleDisplayP3() {
        static const RGBColorspace displayP3(Whitepoint::D65(), 0.680, 0.320, 0.265, 0.690, 0.150, 0.060);
        return displayP3;
    }
}