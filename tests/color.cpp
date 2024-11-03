#include "doctest/doctest.h"
#include "../include/color.h"

using namespace colorama;

TEST_CASE("Color Constructor") {
    Color c(1, 2, 3);
    CHECK(c.r() == 1);
    CHECK(c.g() == 2);
    CHECK(c.b() == 3);
}

TEST_CASE("Color Conversion from vec3d") {
    spacely::Vec3d v(1, 2, 3);
    Color c(v);
    CHECK(c.r() == 1);
    CHECK(c.g() == 2);
    CHECK(c.b() == 3);
}

TEST_CASE("Color Conversion to sRGB") {
    Color c(0.25, 0.5, 0.75); // Color in linear RGB
    Color srgb = c.to_srgb();
    CHECK(srgb.r() == doctest::Approx(0.537099));
    CHECK(srgb.g() == doctest::Approx(0.735357));
    CHECK(srgb.b() == doctest::Approx(0.880825));
}
TEST_CASE("Color Conversion to XYZ") {
    Color c(0.25, 0.5, 0.75); // Color in linear RGB
    Color xyz = c.to_xyz();
    CHECK(xyz.x() == doctest::Approx(0.41723));
    CHECK(xyz.y() == doctest::Approx(0.464876));
    CHECK(xyz.z() == doctest::Approx(0.777158));
}
TEST_CASE("Color Conversion to HSL") {
    Color c(0.25, 0.5, 0.75); // Color in linear RGB
    Color hsl = c.to_hsl();
    CHECK(hsl.r() == doctest::Approx(205.393));
    CHECK(hsl.g() == doctest::Approx(59.0518));
    CHECK(hsl.b() == doctest::Approx(70.8962));
}
TEST_CASE("Color Conversion to LAB") {
    Color c(0.25, 0.5, 0.75); // Color in linear RGB
    Color lab = c.to_lab(colorama::Whitepoint::D65());
    CHECK(lab.r() == doctest::Approx(73.8608));
    CHECK(lab.g() == doctest::Approx(-7.33197));
    CHECK(lab.b() == doctest::Approx(-23.804));
}

TEST_CASE("Color Conversion to sRGB and back") {
    Color c(0.25, 0.5, 0.75); // Color in linear RGB
    Color srgb = c.to_srgb();
    Color original = srgb.to_rgb();
    CHECK(original.r() == doctest::Approx(0.25));
    CHECK(original.g() == doctest::Approx(0.5));
    CHECK(original.b() == doctest::Approx(0.75));
}

TEST_CASE("Color Conversion to XYZ and back") {
    Color c(0.25, 0.5, 0.75); // Color in linear RGB
    Color xyz = c.to_xyz();
    Color original = xyz.to_rgb();
    CHECK(original.r() == doctest::Approx(0.25));
    CHECK(original.g() == doctest::Approx(0.5));
    CHECK(original.b() == doctest::Approx(0.75));
}
TEST_CASE("sRGB to HSL and back") {
    Color c(0.25, 0.5, 0.75, Color::ColorSpace::SRGB); // Color in sRGB
    Color hsl = c.to_hsl();
    Color original = hsl.to_srgb();
    CHECK(original.r() == doctest::Approx(0.25));
    CHECK(original.g() == doctest::Approx(0.5));
    CHECK(original.b() == doctest::Approx(0.75));
}
TEST_CASE("Manual Color Conversion to HSL and back") {
    Color c(0.25, 0.5, 0.75); // Color in linear RGB
    Color srgb = c.to_srgb();
    Color hsl = srgb.to_hsl();
    Color originalsRGB = hsl.to_srgb();
    Color original = originalsRGB.to_rgb();
    CHECK(original.r() == doctest::Approx(0.25).epsilon(0.01));
    CHECK(original.g() == doctest::Approx(0.5).epsilon(0.01));
    CHECK(original.b() == doctest::Approx(0.75).epsilon(0.01));
}
TEST_CASE("Color Conversion to HSL and back") {
    Color c(0.25, 0.5, 0.75); // Color in linear RGB
    Color hsl = c.to_hsl();
    Color original = hsl.to_rgb();
    CHECK(original.r() == doctest::Approx(0.25).epsilon(0.01));
    CHECK(original.g() == doctest::Approx(0.5).epsilon(0.01));
    CHECK(original.b() == doctest::Approx(0.75).epsilon(0.01));
}

TEST_CASE("Color Conversion to LAB and back") {
    Color c(0.25, 0.5, 0.75); // Color in linear RGB
    Color lab = c.to_lab(colorama::Whitepoint::D65());
    Color original = lab.to_rgb();
    CHECK(original.r() == doctest::Approx(0.25));
    CHECK(original.g() == doctest::Approx(0.5));
    CHECK(original.b() == doctest::Approx(0.75));
}
