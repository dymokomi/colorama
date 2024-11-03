#pragma once

#include <vector>
#include <thread>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vec3d.h>
#include "color.h"
#include "observer.h"

namespace colorama {
    class Spectrum {
        public:
            static constexpr int RESPONSE_SAMPLES = 31; 
            static constexpr double START_WAVELENGTH = 400.0;
            static constexpr double END_WAVELENGTH = 700.0;
            static void setLookupTable(const std::vector<std::vector<std::vector<spacely::Vec3d>>>& table, double step = 0.01f);
            static const std::vector<std::vector<std::vector<spacely::Vec3d>>>& getLookupTable();

            Spectrum();
            Spectrum(const std::vector<double>& data);
            Spectrum(const double* data);
            Spectrum(double r, double g, double b, double coeff_a, double coeff_b, double coeff_c);
            Spectrum(double r, double g, double b);
            Spectrum(Color c);
            double& operator[](int index);
            const double& operator[](int index) const ;
            
            int num_wavelengths() const ;

            Color to_rgb(Observer * observer) const;
            Color to_XYZ(Observer * observer) const;
               
            static std::vector<std::vector<std::vector<spacely::Vec3d>>> load_lookup_tables(double step = 0.01);

            Spectrum& operator+=(const Spectrum& v);
            Spectrum& operator-=(const Spectrum& v);
            Spectrum& operator*=(double t);
            Spectrum& operator*=(const Spectrum& v);
            Spectrum operator*(double t) const;
            Spectrum operator/(double t) const;
            Spectrum operator*(const Spectrum& v) const;
            Spectrum operator+(const Spectrum& v) const;
            Spectrum operator-(const Spectrum& v) const;
            Spectrum& operator/=(double t);
            Spectrum& operator/=(const Spectrum& v);


            double get_wavelength(int index) const;
            friend std::ostream& operator<<(std::ostream& os, const Spectrum& s);
            static Spectrum d65();
            static Spectrum d65_cb();
            static Spectrum d50_cb();
            static Spectrum d50();
            static Spectrum studio_led();
            const std::vector<double>& get_data() const;

        private:
            std::vector<double> data_;
            static inline std::vector<std::vector<std::vector<spacely::Vec3d>>> lookup_table;
            static inline double step;
            // D50 from 400 to 700nm in 10nm step
            static inline double d50_spd[31] = {
                18.0306416435159065, 20.6652060223712901, 21.9527407119419564, 21.1422996064334434, 27.3614468598853620, 
                31.9039467673700941, 33.1344364664343658, 33.4108486992791427, 34.7786120485748356, 33.6282778746323672, 
                35.0035740880590254, 35.3288768031519211, 35.5174543775374971, 37.3348492673887336, 36.8433847337951619, 
                37.4145659848912615, 36.5673016066648415, 35.7390522252738805, 36.1716434032807328, 34.1899881946123472, 
                35.7217924589155302, 36.3000311992217277, 36.2168040207649611, 35.0028793093285131, 36.1494104839038783, 
                34.9829135626512695, 35.9053603129809957, 37.6654176739129838, 36.2502631017350510, 31.9528372496181987, 
                33.4969281272612278
            };
            static inline double d65_spd[31] = {
                28.8166638832520299, 31.8569814237367765, 32.5345420827917948, 30.1842513703383979, 36.5157768073820819, 
                40.7441759660340708, 41.0241424424860455, 39.9965540444639700, 40.3663605096281231, 37.8898411308639851, 
                38.0789229675722112, 37.5384901672569811, 36.4896605315936569, 37.4991416450690807, 36.3555969825463876, 
                36.2305870757724406, 34.8217010512392875, 33.5452071341029594, 33.3550110029610920, 30.8818345074978744, 
                31.3416898915805469, 31.1999307466009554, 30.5381791398231996, 29.0025073017624955, 29.1454852062788845, 
                27.8666930568731601, 27.9320882114473967, 28.6505295475365607, 27.2598900943542723, 24.2781426550377049, 
                24.9355067274829914
            };
            static inline double d65_cb_spd[31] = {
                27.4641780627958099, 32.6944719012318643, 35.4812882740959665, 37.4376353241936073, 38.2055015413569308, 
                38.9342408175183081, 39.0643378963497980, 38.8344671179633281, 38.2866899439359827, 37.9316129543432652, 
                37.4200282007427276, 36.3440373231890277, 34.9931796850965995, 33.4750543742208322, 33.8291531902884941, 
                36.5563009781246180, 39.4966905944213806, 36.7646519389600215, 32.9116264237926970, 28.7289564306839367, 
                29.1534837405551208, 29.8215762581634714, 29.8127726964380386, 28.3200798972135175, 26.2218976859837980, 
                24.8690837008412835, 25.0881945704522131, 26.7833692893618220, 29.3598783543404203, 32.3609147291902133, 
                34.1734702711056855
            };
            static inline double d50_cb_spd[31] = {
                19.8664266423976557, 23.6741337318269345, 26.1948328629854927, 28.3615677800741857, 29.7344750036840431, 
                31.2421552125673578, 32.2299967639673284, 32.9408872357394245, 33.4133331951046344, 33.9331718525879751, 
                34.3774783974455289, 34.3345287647759605, 34.1449579723034091, 33.7480441255639931, 34.5818594080800068, 
                37.2151161972691469, 39.9994372117098465, 38.5717321809008880, 36.1221220969195471, 33.1541543772710483, 
                33.7598923000935329, 34.7136703497210846, 34.9817352984518095, 34.0590587069642865, 32.4625171891094553, 
                31.4435741795694454, 31.8908427680593967, 33.4873842859142243, 35.9814250243813163, 38.9271774167869395, 
                40.7710495779458100
            };
            static inline double studio_led_spd[31] = {
                0.0076918600433055, 0.1269077935276378, 0.5329653496946619, 2.3022668977667857, 9.3757140290895560, 
                37.9950088041393386, 63.2875346245249659, 45.2798004513542409, 31.0892347451363449, 25.8727005533851084, 
                27.3306791993759290, 30.7703147492545632, 33.7864399599110428, 36.2222798769514114, 38.2913744342267748, 
                39.8424451194275449, 41.1835138105096377, 42.2873825957799028, 43.4254462003392092, 44.9514354199367361, 
                46.8206837654504255, 48.2153901499761943, 48.3188117096549590, 47.5153045310449471, 45.4908290797702719, 
                42.0532783872339735, 37.3181882978236601, 32.1100725073108038, 27.0223572283176772, 22.3252684022835552, 
                17.7050507935610675
            };
            double integrate(const double* y, int length) const;
            double spectrum_function(double lambda, double coeff_a, double coeff_b, double coeff_c);
    
    };
}