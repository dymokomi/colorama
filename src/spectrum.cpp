#include "spectrum.h"

namespace colorama {

    void Spectrum::setLookupTable(const std::vector<std::vector<std::vector<spacely::Vec3d>>>& table, double step) { 
        Spectrum::lookup_table = table; Spectrum::step = step; 
    }   
    const std::vector<std::vector<std::vector<spacely::Vec3d>>>& Spectrum::getLookupTable() { return lookup_table; }

    Spectrum::Spectrum() : data_(RESPONSE_SAMPLES) {}
    Spectrum::Spectrum(const std::vector<double>& data) : data_(data) { }
    Spectrum::Spectrum(const double* data) : data_(data, data + RESPONSE_SAMPLES) {}
    Spectrum::Spectrum(double r, double g, double b, double coeff_a, double coeff_b, double coeff_c) : data_(RESPONSE_SAMPLES) {
        // Convert RGB to XYZ
        data_ = std::vector<double>(RESPONSE_SAMPLES);
        for(int i = 0; i < RESPONSE_SAMPLES; i++){
            data_[i] = spectrum_function(START_WAVELENGTH + i * (END_WAVELENGTH - START_WAVELENGTH) / RESPONSE_SAMPLES, coeff_a, coeff_b, coeff_c);
        }
    }
    Spectrum::Spectrum(double r, double g, double b) : Spectrum(Color(r, g, b, Color::ColorSpace::RGB_LIN)) {}
    Spectrum::Spectrum(Color c)
        : data_(RESPONSE_SAMPLES) {
        // Convert RGB to XYZ
        c.set_color_space(Color::ColorSpace::RGB_LIN);
        data_ = std::vector<double>(RESPONSE_SAMPLES);

        spacely::Vec3d rgb_coeffs;

        int size = static_cast<int>(1.0f / Spectrum::step + 0.5f) + 1;
        double r_scaled = c.x() * (size - 1);
        double g_scaled = c.y() * (size - 1);
        double b_scaled = c.z() * (size - 1);

        int r_index = std::min(static_cast<int>(r_scaled), size - 1);
        int g_index = std::min(static_cast<int>(g_scaled), size - 1);
        int b_index = std::min(static_cast<int>(b_scaled), size - 1);

        double r_frac = r_scaled - r_index;
        double g_frac = g_scaled - g_index;
        double b_frac = b_scaled - b_index;

        // Check if we need to interpolate
        if (r_frac < 1e-6 && g_frac < 1e-6 && b_frac < 1e-6) {
            // Exact match, no interpolation needed
            rgb_coeffs = Spectrum::lookup_table[r_index][g_index][b_index];
        } else {
            // Trilinear interpolation
            spacely::Vec3d c000 = Spectrum::lookup_table[r_index][g_index][b_index];
            spacely::Vec3d c001 = (b_index < size - 1) ? Spectrum::lookup_table[r_index][g_index][b_index + 1] : c000;
            spacely::Vec3d c010 = (g_index < size - 1) ? Spectrum::lookup_table[r_index][g_index + 1][b_index] : c000;
            spacely::Vec3d c011 = (g_index < size - 1 && b_index < size - 1) ? Spectrum::lookup_table[r_index][g_index + 1][b_index + 1] : c010;
            spacely::Vec3d c100 = (r_index < size - 1) ? Spectrum::lookup_table[r_index + 1][g_index][b_index] : c000;
            spacely::Vec3d c101 = (r_index < size - 1 && b_index < size - 1) ? Spectrum::lookup_table[r_index + 1][g_index][b_index + 1] : c100;
            spacely::Vec3d c110 = (r_index < size - 1 && g_index < size - 1) ? Spectrum::lookup_table[r_index + 1][g_index + 1][b_index] : c100;
            spacely::Vec3d c111 = (r_index < size - 1 && g_index < size - 1 && b_index < size - 1) ? Spectrum::lookup_table[r_index + 1][g_index + 1][b_index + 1] : c110;

            spacely::Vec3d c00 = c000 * (1 - r_frac) + c100 * r_frac;
            spacely::Vec3d c01 = c001 * (1 - r_frac) + c101 * r_frac;
            spacely::Vec3d c10 = c010 * (1 - r_frac) + c110 * r_frac;
            spacely::Vec3d c11 = c011 * (1 - r_frac) + c111 * r_frac;

            spacely::Vec3d c0 = c00 * (1 - g_frac) + c10 * g_frac;
            spacely::Vec3d c1 = c01 * (1 - g_frac) + c11 * g_frac;

            rgb_coeffs = c0 * (1 - b_frac) + c1 * b_frac;
        }

        for (int i = 0; i < RESPONSE_SAMPLES; i++) {
            double wavelength = START_WAVELENGTH + i * (END_WAVELENGTH - START_WAVELENGTH) / RESPONSE_SAMPLES;
            data_[i] = spectrum_function(wavelength, rgb_coeffs.x(), rgb_coeffs.y(), rgb_coeffs.z());
        }
    }
    double& Spectrum::operator[](int index) { return data_[index]; }
    const double& Spectrum::operator[](int index) const { return data_[index]; }
    
    int Spectrum::num_wavelengths() const { return static_cast<int>(data_.size()); }

    Color Spectrum::to_rgb(Observer * observer) const {
        Color xyz = to_XYZ(observer);
        Color rgb = xyz.to_rgb();

        return rgb;
    }
    Color Spectrum::to_XYZ(Observer * observer) const {
        Color xyz;
        
        std::vector<double> x_bar(observer->get_length());
        std::vector<double> y_bar(observer->get_length());
        std::vector<double> z_bar(observer->get_length());

        // Multiply spectrum by each cone fundamental
        int i;
        // Calculate the normalization factor for illuminant
        double illuminant_norm = 0.0;
        for (int i = 0; i < observer->get_length(); i++) {
            illuminant_norm +=  observer->y_bar[i];// * d65_spd[i];
        }
        //illuminant_norm /= 16.0;

        for(i = 0; i < observer->get_length()-1; i++){
        
            x_bar[i] = data_[i] * observer->x_bar[i];// * d65_spd[i];
            y_bar[i] = data_[i] * observer->y_bar[i];// * d65_spd[i];
            z_bar[i] = data_[i] * observer->z_bar[i];// * d65_spd[i];
        }

        // Integrate
        // Normalize
        double X = integrate(x_bar.data(), observer->get_length()) / illuminant_norm;
        double Y = integrate(y_bar.data(), observer->get_length()) / illuminant_norm;
        double Z = integrate(z_bar.data(), observer->get_length()) / illuminant_norm;

        Color result = Color(X, Y, Z);
        result.set_color_space(Color::ColorSpace::XYZ);
        return result;
    }
        
    std::vector<std::vector<std::vector<spacely::Vec3d>>> Spectrum::load_lookup_tables(double step) {
        const int size = static_cast<int>(1.0 / step + 0.5) + 1;
        const int total_size = size * size * size;
        std::cout << "Loading lookup table with dimensions: " << size << "x" << size << "x" << size << std::endl;
        
        std::ifstream infile("lookup_table.bin", std::ios::binary);
        std::vector<std::vector<std::vector<spacely::Vec3d>>> lookup_table(size, std::vector<std::vector<spacely::Vec3d>>(size, std::vector<spacely::Vec3d>(size)));

        if (infile.is_open()) {
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    for (int k = 0; k < size; ++k) {
                        infile.read(reinterpret_cast<char*>(&lookup_table[i][j][k]), sizeof(spacely::Vec3d));
                    }
                }
            }
            infile.close();

            if (infile.gcount() == sizeof(spacely::Vec3d)) {
                std::cout << "Lookup table loaded from lookup_table.bin" << std::endl;
            } else {
                std::cerr << "Error: Incomplete data in lookup_table.bin" << std::endl;
                lookup_table.clear();
            }
        } else {
            std::cerr << "Unable to open file for reading" << std::endl;
        }

        return lookup_table;
    }

    Spectrum& Spectrum::operator+=(const Spectrum& v) {
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] += v.data_[i];
        }
        return *this;
    }
    Spectrum& Spectrum::operator-=(const Spectrum& v) {
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] -= v.data_[i];
        }
        return *this;
    }

    Spectrum& Spectrum::operator*=(double t) {
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] *= t;
        }
        return *this;
    }
    Spectrum& Spectrum::operator*=(const Spectrum& v) {
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] *= v.data_[i];
        }
        return *this;
    }
    Spectrum Spectrum::operator*(double t) const {
        Spectrum result = *this;
        for (size_t i = 0; i < result.data_.size(); ++i) {
            result.data_[i] *= t;
        }
        return result;
    }
    Spectrum Spectrum::operator/(double t) const {
        Spectrum result = *this;
        for (size_t i = 0; i < result.data_.size(); ++i) {
            result.data_[i] /= t;
        }
        return result;
    }
    Spectrum Spectrum::operator*(const Spectrum& v) const {
        Spectrum result = *this;
        for (size_t i = 0; i < result.data_.size(); ++i) {
            result.data_[i] *= v.data_[i];
        }
        return result;
    }
    Spectrum Spectrum::operator+(const Spectrum& v) const {
        Spectrum result = *this;
        for (size_t i = 0; i < result.data_.size(); ++i) {
            result.data_[i] += v.data_[i];
        }
        return result;
    }
    Spectrum Spectrum::operator-(const Spectrum& v) const {
        Spectrum result = *this;
        for (size_t i = 0; i < result.data_.size(); ++i) {
            result.data_[i] -= v.data_[i];
        }
        return result;
    }
    Spectrum& Spectrum::operator/=(double t) {
        return *this *= 1/t;
    }
    Spectrum& Spectrum::operator/=(const Spectrum& v) {
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] /= v.data_[i];
        }
        return *this;
    }


    double Spectrum::get_wavelength(int index) const {
        if (index < 0 || index >= static_cast<int>(data_.size())) {
            throw std::out_of_range("Index out of range");
        }
        return START_WAVELENGTH + index * (END_WAVELENGTH - START_WAVELENGTH) / (data_.size() - 1);
    }



    Spectrum Spectrum::d65(){
        Spectrum s;
        for(int i = 0; i < RESPONSE_SAMPLES; i++){
            s[i] = d65_spd[i];
        }
        return s;
    }
    Spectrum Spectrum::d65_cb(){
        Spectrum s;
        for(int i = 0; i < RESPONSE_SAMPLES; i++){
            s[i] = d65_cb_spd[i];
        }
        return s;
    }
    Spectrum Spectrum::d50_cb(){
        Spectrum s;
        for(int i = 0; i < RESPONSE_SAMPLES; i++){
            s[i] = d50_cb_spd[i];
        }
        return s;
    }
    Spectrum Spectrum::d50(){
        Spectrum s;
        for(int i = 0; i < RESPONSE_SAMPLES; i++){
            s[i] = d50_spd[i];
        }
        return s;
    }
    Spectrum Spectrum::studio_led(){
        Spectrum s;
        for(int i = 0; i < RESPONSE_SAMPLES; i++){
            s[i] = studio_led_spd[i];
        }
        return s;
    }
    const std::vector<double>& Spectrum::get_data() const {
        return data_;
    }

    double Spectrum::integrate(const double* y, int length) const {
        double sum = 0;
        int i;
        for(i = 0; i < length-1; i++){
            sum += y[i];
        }
        return sum;
    }

    // modesl a spectrum function from a polynomial
    double Spectrum::spectrum_function(double lambda, double coeff_a, double coeff_b, double coeff_c) {
        // Normalize lambda to 0-1 range
        double x = (lambda - START_WAVELENGTH) / (END_WAVELENGTH- START_WAVELENGTH);
        double y = coeff_a * x * x + coeff_b * x + coeff_c;
        return 0.5 * y / std::sqrt(1.0 + y * y) + 0.5;
    }
}