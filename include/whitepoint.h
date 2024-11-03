#pragma once

namespace colorama {
    class Whitepoint {
    public:
        double x;
        double y;
        double z;

        Whitepoint(double x, double y, double z);

        static const Whitepoint& A();
        static const Whitepoint& B();
        static const Whitepoint& C();
        static const Whitepoint& D50();
        static const Whitepoint& D55();
        static const Whitepoint& D65();
        static const Whitepoint& D75();
        static const Whitepoint& E();
        static const Whitepoint& F2();
        static const Whitepoint& F7();
        static const Whitepoint& F11();
    };
}