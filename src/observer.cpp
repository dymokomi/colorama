#include "observer.h"

namespace colorama {
    Observer::Observer(enum Observer::StandardObserver stdobs, int samples, double start, double end) {
        switch (stdobs) {
            case CIE1931_2Deg:
                length = 471;
                wavelengths = CIE1931_2Deg_wave;
                x_bar = CIE1931_2Deg_x;
                y_bar = CIE1931_2Deg_y;
                z_bar = CIE1931_2Deg_z;
                break;

            case CIE1931_2Deg_Judd:
                length = 41;
                wavelengths = CIE1931_2Deg_Judd_wave;
                x_bar = CIE1931_2Deg_Judd_x;
                y_bar = CIE1931_2Deg_Judd_y;
                z_bar = CIE1931_2Deg_Judd_z;
                break;
            
            case CIE1931_2Deg_JuddVos:
                length = 90;
                wavelengths = CIE1931_2Deg_JuddVos_wave;
                x_bar = CIE1931_2Deg_JuddVos_x;
                y_bar = CIE1931_2Deg_JuddVos_y;
                z_bar = CIE1931_2Deg_JuddVos_z;
                break;

            case CIE1964_10Deg:
                length = 471;
                wavelengths = CIE1964_10Deg_wave;
                x_bar = CIE1964_10Deg_x;
                y_bar = CIE1964_10Deg_y;
                z_bar = CIE1964_10Deg_z;
                break;

            case StilesBurch1955_2Deg:
                length = 69;
                wavelengths = StilesBurch1955_2Deg_wave;
                x_bar = StilesBurch1955_2Deg_x;
                y_bar = StilesBurch1955_2Deg_y;
                z_bar = StilesBurch1955_2Deg_z;
                break;

            case StilesBurch1955_10Deg:
                length = 89;
                wavelengths = StilesBurch1955_10Deg_wave;
                x_bar = StilesBurch1955_10Deg_x;
                y_bar = StilesBurch1955_10Deg_y;
                z_bar = StilesBurch1955_10Deg_z;
                break;

            case CIE2006_2Deg:
                length = 441;
                wavelengths = CIE2006_2Deg_wave;
                x_bar = CIE2006_2Deg_x;
                y_bar = CIE2006_2Deg_y;
                z_bar = CIE2006_2Deg_z;
                break;

            case CIE2006_10Deg:
                length = 441;
                wavelengths = CIE2006_10Deg_wave;
                x_bar = CIE2006_10Deg_x;
                y_bar = CIE2006_10Deg_y;
                z_bar = CIE2006_10Deg_z;
                break;
            
            default:
                break;
        }
        interpolateXYZSimple(samples, start, end);
    }
    const int Observer::get_length() const {
        return length;
    }

    void Observer::interpolateXYZSimple(int numSamples, double startWavelength, double endWavelength) {
        std::vector<double> newWavelengths(numSamples);
        std::vector<double> newX(numSamples);
        std::vector<double> newY(numSamples);
        std::vector<double> newZ(numSamples);
        // length is already a member variable, so we don't need to redefine it here
        // wavelengths is a double*, so we can't use .size() on it
        // We'll use the existing length member variable instead
        double step = (endWavelength - startWavelength) / (numSamples - 1);
        for (int i = 0; i < numSamples; ++i) {
            newWavelengths[i] = startWavelength + i * step;
            for (int j = 0; j < length; j++) {
                if( int(newWavelengths[i]) == int(wavelengths[j])) {
                    newX[i] = x_bar[j];
                    newY[i] = y_bar[j];
                    newZ[i] = z_bar[j];
                }
            }
        }
        std::copy(newX.begin(), newX.end(), x_bar);
        std::copy(newY.begin(), newY.end(), y_bar);
        std::copy(newZ.begin(), newZ.end(), z_bar);
        length = numSamples;
    }
    void Observer::interpolateXYZ(int numSamples, double startWavelength, double endWavelength) {
        if (numSamples <= 0 || startWavelength >= endWavelength) {
            return;  // Invalid input, do nothing
        }

        std::vector<double> newWavelengths(numSamples);
        std::vector<double> newX(numSamples);
        std::vector<double> newY(numSamples);
        std::vector<double> newZ(numSamples);

        double step = (endWavelength - startWavelength) / (numSamples - 1);

        

        for (int i = 0; i < numSamples; ++i) {
            double wavelength = newWavelengths[i];
            
            // Find the two nearest points in the original data
            auto it = std::lower_bound(wavelengths, wavelengths + length, wavelength);
            int index = static_cast<int>(std::distance(wavelengths, it));

            if (index == 0) {
                newX[i] = x_bar[0];
                newY[i] = y_bar[0];
                newZ[i] = z_bar[0];
            } else if (index == length) {
                newX[i] = x_bar[length - 1];
                newY[i] = y_bar[length - 1];
                newZ[i] = z_bar[length - 1];
            } else {
                double t = (wavelength - wavelengths[index - 1]) / (wavelengths[index] - wavelengths[index - 1]);
                newX[i] = x_bar[index - 1] + t * (x_bar[index] - x_bar[index - 1]);
                newY[i] = y_bar[index - 1] + t * (y_bar[index] - y_bar[index - 1]);
                newZ[i] = z_bar[index - 1] + t * (z_bar[index] - z_bar[index - 1]);
            }
        }

        // Overwrite the original arrays with the new interpolated values
        length = numSamples;
        std::copy(newWavelengths.begin(), newWavelengths.end(), wavelengths);
        std::copy(newX.begin(), newX.end(), x_bar);
        std::copy(newY.begin(), newY.end(), y_bar);
        std::copy(newZ.begin(), newZ.end(), z_bar);
    }
}