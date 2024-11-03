#pragma once
#include <vec3d.h>
#include "colorspace.h"

#define CIELAB_E 0.008856
#define CIELAB_K 903.3
#define CIELAB_KE pow((CIELAB_E * CIELAB_K + 16.0)/116.0, 3.0)

namespace colorama {
    class Color : public spacely::Vec3d {
        public:
            enum ColorSpace {
                RGB_LIN,
                SRGB,
                XYZ,
                LAB,
                HSL,
                NO_COLORSPACE
            };

            using spacely::Vec3d::Vec3d;

            double r() const;
            double g() const;
            double b() const;

            // Add conversion constructor from vec3
            Color(const Vec3d& v);
            Color(double r, double g, double b, ColorSpace color_space);
            const Color to_rgbDisplayP3();
            const Color to_rgb();
            const Color to_xyz(const Whitepoint& wp = Whitepoint::D65()) ;
            const Color to_srgb();
            // From any color space to LAB
            const Color to_lab(Whitepoint wp);
            // From any to HSL
            const Color to_hsl();
            // Add assignment operator from vec3
            Color& operator=(const Vec3d& v);

            operator Vec3d() const;

            void set_color_space(ColorSpace color_space) ;
            ColorSpace get_color_space() const ;
            void clamp() ;

        private:
            ColorSpace color_space_ = ColorSpace::RGB_LIN;
            const spacely::mat3x3 xyz_rgb_matrix(const RGBColorspace& cs);
            const spacely::mat3x3 rgb_xyz_matrix(const RGBColorspace& cs);
            
            double linear_to_srgb(double c);
            double srgb_to_linear(double c);
            double cielab_cube(double c, double threshold);
            double cielab_cuberoot(double c);
            
    };
}