#include "color.h"
#include <iostream>
namespace colorama {
    Color::Color(const Vec3d& v) : Vec3d(v) {}
    
    Color::Color(double r, double g, double b, ColorSpace color_space) : Vec3d(r, g, b), color_space_(color_space) {}
    double Color::r() const {
        return x();
    }
    double Color::g() const {
        return y();
    }
    double Color::b() const {
        return z();
    }
    const Color Color::to_rgbDisplayP3() {
        if (color_space_ == ColorSpace::XYZ) {
            const spacely::mat3x3& m = xyz_rgb_matrix(RGBColorspace::AppleDisplayP3());
            Color result = spacely::Vec3d::mat3x3_times(m, Color(x(), y(), z()));
            result.set_color_space(ColorSpace::RGB_LIN);
            return result;
        } else {
            throw std::range_error("Color space not supported");
        }
    }
    const Color Color::to_rgb() {
        if(color_space_ == ColorSpace::RGB_LIN) {
            return *this;
        } else if (color_space_ == ColorSpace::SRGB) {
            Color result(srgb_to_linear(x()), srgb_to_linear(y()), srgb_to_linear(z()));
            result.set_color_space(ColorSpace::RGB_LIN);
            return result;

        } else if (color_space_ == ColorSpace::XYZ) {
            const spacely::mat3x3& m = xyz_rgb_matrix(RGBColorspace::sRGB());
            Color result = spacely::Vec3d::mat3x3_times(m, Color(x(), y(), z()));
            result.set_color_space(ColorSpace::RGB_LIN);
            return result;

        } else if (color_space_ == ColorSpace::LAB) {
            Color xyz_color = to_xyz();
            Color result = xyz_color.to_rgb();
            result.set_color_space(ColorSpace::RGB_LIN);
            return result;

        } else if (color_space_ == ColorSpace::HSL) {
            Color srgb_color = to_srgb();
            Color result = srgb_color.to_rgb();
            result.set_color_space(ColorSpace::RGB_LIN);
            return result;

        } else {
            std::cerr << "Color space not supported" << std::endl;
            return Color(0, 0, 0);
        }
    }
    const Color Color::to_xyz(const Whitepoint& wp) {
        if(color_space_ == ColorSpace::RGB_LIN) {
            const spacely::mat3x3& m = rgb_xyz_matrix(RGBColorspace::sRGB());
            Color result = spacely::Vec3d::mat3x3_times(m, Color(x(), y(), z()));
            result.set_color_space(ColorSpace::XYZ);
            return result;
        } else if (color_space_ == ColorSpace::SRGB) {
            Color rgb_lin = to_rgb();
            Color result = rgb_lin.to_xyz(wp);
            result.set_color_space(ColorSpace::XYZ);
            return result;
        } else if (color_space_ == ColorSpace::XYZ) {
            return *this;
        } else if (color_space_ == ColorSpace::LAB) {
            double L = x();
            double a = y();
            double b = z();
            double Y_ = (L + 16.0)/116.0;
            double X_ = a/500.0 + Y_; 
            double Z_ = Y_ - b/200.0;

            X_ = cielab_cube(X_, CIELAB_E)  * wp.x;
            Y_ = cielab_cube(Y_, CIELAB_KE) * wp.y;
            Z_ = cielab_cube(Z_, CIELAB_E)  * wp.z;
            Color result(X_, Y_, Z_);
            result.set_color_space(ColorSpace::XYZ);
            return result;
        } else {
            std::cerr << "Color space not supported" << std::endl;
            return Color(0, 0, 0);
        }
    }
    // From any color space to sRGB
    const Color Color::to_srgb(){
        if(color_space_ == ColorSpace::RGB_LIN) {
            Color result(linear_to_srgb(x()), linear_to_srgb(y()), linear_to_srgb(z()));
            result.set_color_space(ColorSpace::SRGB);
            return result;
        } else if (color_space_ == ColorSpace::SRGB) {
            return *this;
        } else if (color_space_ == ColorSpace::XYZ) {
            Color rgb_color = to_rgb();
            Color result = rgb_color.to_srgb();
            result.set_color_space(ColorSpace::SRGB);
            return result;
        } else if (color_space_ == ColorSpace::LAB) {
            Color rgb_color = to_rgb();
            Color result = rgb_color.to_srgb();
            result.set_color_space(ColorSpace::SRGB);
            return result;
        } else if (color_space_ == ColorSpace::HSL) {
            const double H = x();
            const double S = y() / 100.0;
            const double L = z() / 100.0;
            
            // Use epsilon for floating point comparisons
            const double epsilon = 1e-10;
            
            // Early return for grayscale
            if (S < epsilon) {
                Color result(L, L, L);
                result.set_color_space(ColorSpace::SRGB);
                return result;
            }
            
            const double C = (1.0 - std::abs(2.0 * L - 1.0)) * S;
            const double m = L - C/2.0;
            
            // Normalize hue to [0, 6)
            double h_prime = std::fmod(H, 360.0) / 60.0;
            if (h_prime < 0.0) h_prime += 6.0;
            
            const double X = C * (1.0 - std::abs(std::fmod(h_prime, 2.0) - 1.0));
            
            double R = 0.0;
            double G = 0.0;
            double B = 0.0;

            // Using h_prime instead of H makes the comparisons more stable
            const int sector = static_cast<int>(h_prime);
            switch (sector) {
                case 0:  // 0 <= h_prime < 1
                    R = C; G = X; B = 0;
                    break;
                case 1:  // 1 <= h_prime < 2
                    R = X; G = C; B = 0;
                    break;
                case 2:  // 2 <= h_prime < 3
                    R = 0; G = C; B = X;
                    break;
                case 3:  // 3 <= h_prime < 4
                    R = 0; G = X; B = C;
                    break;
                case 4:  // 4 <= h_prime < 5
                    R = X; G = 0; B = C;
                    break;
                default: // 5 <= h_prime < 6
                    R = C; G = 0; B = X;
                    break;
            }

            // Add the lightness adjustment and ensure values are in [0, 1]
            R = std::max(0.0, std::min(1.0, R + m));
            G = std::max(0.0, std::min(1.0, G + m));
            B = std::max(0.0, std::min(1.0, B + m));

            Color result(R, G, B);
            result.set_color_space(ColorSpace::SRGB);
            return result;
        } else {
            std::cerr << "Color space not supported" << std::endl;
            return Color(0, 0, 0);
        }
    }
    // From any color space to LAB
    const Color Color::to_lab(Whitepoint wp){
        if (color_space_ == ColorSpace::RGB_LIN) {
            Color color_xyz = to_xyz();
            Color result = color_xyz.to_lab(wp);
            result.set_color_space(ColorSpace::LAB);
            return result;
        } else if(color_space_ == ColorSpace::SRGB) {
            Color color_xyz = to_xyz();
            Color result = color_xyz.to_lab(wp);
            result.set_color_space(ColorSpace::LAB);
            return result;
        } else if(color_space_ == ColorSpace::XYZ){
            double X_ = cielab_cuberoot(x()/wp.x);
            double Y_ = cielab_cuberoot(y()/wp.y);
            double Z_ = cielab_cuberoot(z()/wp.z);
            double L = 116.0 * Y_ - 16.0;
            double a = 500.0 * (X_ - Y_);
            double b = 200.0 * (Y_ - Z_);
            Color result(L, a, b);  
            result.set_color_space(ColorSpace::LAB);
            return result;
        } else if (color_space_ == ColorSpace::LAB) {
            return *this;
        } else if (color_space_ == ColorSpace::HSL) {
            Color color_xyz = to_xyz();
            Color result = color_xyz.to_lab(wp);
            result.set_color_space(ColorSpace::LAB);
            return result;
        } else {
            std::cerr << "Color space not supported" << std::endl;
            return Color(0, 0, 0);
        }
    }
    // From any to HSL
    const Color Color::to_hsl() {
        if (color_space_ == ColorSpace::RGB_LIN) {
            Color srgb_color = to_srgb();
            Color result = srgb_color.to_hsl();
            result.set_color_space(ColorSpace::HSL);
            return result;

        } else if (color_space_ == ColorSpace::SRGB) {
            const double R = x();
            const double G = y();
            const double B = z();

            const double cmin = std::min(R, std::min(G, B));
            const double cmax = std::max(R, std::max(G, B));
            const double delta = cmax - cmin;

            double H = 0.0;
            double S = 0.0;
            double L = (cmax + cmin) / 2.0;

            // Use small epsilon for floating point comparisons
            const double epsilon = 1e-10;

            if (delta > epsilon) {  // If not grayscale
                if (std::abs(cmax - R) < epsilon) {
                    H = ((G - B) / delta);
                    if (H < 0.0) H += 6.0;
                } else if (std::abs(cmax - G) < epsilon) {
                    H = ((B - R) / delta) + 2.0;
                } else { // cmax == B
                    H = ((R - G) / delta) + 4.0;
                }

                H *= 60.0;
                
                // Calculate Saturation
                if (L < 0.5) {
                    S = delta / (cmax + cmin);
                } else {
                    S = delta / (2.0 - cmax - cmin);
                }
                S *= 100.0;
            }

            L *= 100.0;

            Color result(H, S, L);
            result.set_color_space(ColorSpace::HSL);
            return result;
        } else if (color_space_ == ColorSpace::XYZ) {
            Color rgb_color = to_rgb();
            Color result = rgb_color.to_hsl();
            result.set_color_space(ColorSpace::HSL);
            return result;
        } else if (color_space_ == ColorSpace::LAB) {
            Color rgb_color = to_rgb();
            Color result = rgb_color.to_hsl();
            result.set_color_space(ColorSpace::HSL);
            return result;
        } else {
            std::cerr << "Color space not supported" << std::endl;
            return Color(0, 0, 0);
        }
    }
    // Add assignment operator from vec3
    Color& Color::operator=(const Vec3d& v) {
        Vec3d::operator=(v);
        return *this;
    }


    // Add conversion operator to vec3
    Color::operator Vec3d() const { return Vec3d(x(), y(), z()); }

    void Color::set_color_space(ColorSpace color_space) {
        color_space_ = color_space;
    }
    Color::ColorSpace Color::get_color_space() const {
        return color_space_;
    }
    void Color::clamp() {
        e[0] = std::max(0.0, std::min(e[0], 1.0));
        e[1] = std::max(0.0, std::min(e[1], 1.0));
        e[2] = std::max(0.0, std::min(e[2], 1.0));
    }


    const spacely::mat3x3 Color::xyz_rgb_matrix(const RGBColorspace& cs) {
        const spacely::mat3x3 m = rgb_xyz_matrix(cs);
        spacely::mat3x3 result = spacely::Vec3d::mat3x3_inverse(m);
        return result;
    }
    const spacely::mat3x3 Color::rgb_xyz_matrix(const RGBColorspace& cs) {
        spacely::mat3x3 m ={{     
                {cs.red_x/cs.red_y, cs.green_x/cs.green_y, cs.blue_x/cs.blue_y},
                {1.0, 1.0, 1.0},
                {(1 - cs.red_x - cs.red_y)/cs.red_y, (1 - cs.green_x - cs.green_y)/cs.green_y, (1 - cs.blue_x - cs.blue_y)/cs.blue_y}
        }};
        spacely::mat3x3 wpXYZ = {{
                {cs.wp.x, cs.wp.x, cs.wp.x},
                {cs.wp.y, cs.wp.y, cs.wp.y},
                {cs.wp.z, cs.wp.z, cs.wp.z}
        }};
        spacely::mat3x3 inv_m = spacely::Vec3d::mat3x3_inverse(m);
        spacely::mat3x3 minv_x_wp  = spacely::Vec3d::mat3x3_mult(inv_m, wpXYZ);
        spacely::mat3x3 minv_x_wp_t = spacely::Vec3d::mat3x3_transpose(minv_x_wp);
        spacely::mat3x3 result = spacely::Vec3d::mat3x3_times(m, minv_x_wp_t);
        return result;
    }
            
    double Color::linear_to_srgb(double c){
        if(c <= 0.0031308) {
            return 12.92*c;
        } else {
            return 1.055*pow(c, 1/2.4) - 0.055;
        }
    }
    double Color::srgb_to_linear(double c){
        if(c <= 0.04045) {
            return c/12.92;
        } else {
            return pow((c + 0.055)/1.055, 2.4);
        }
    }
    double Color::cielab_cube(double c, double threshold){
        if(c * c * c <= threshold)
            return (116.0 * c - 16.0) / CIELAB_K;
        else
            return c * c * c;
    }
    double Color::cielab_cuberoot(double c){
        if (c <= CIELAB_E)
            return (CIELAB_K * c + 16.0)/116.0;
        else
            return pow(c, 1.0/3.0);
    }
}