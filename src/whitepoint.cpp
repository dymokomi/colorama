
#include "whitepoint.h"

namespace colorama {
    Whitepoint::Whitepoint(double x, double y, double z) : x(x), y(y), z(z) {}

    const Whitepoint& Whitepoint::A() {
        static const Whitepoint a(1.09850, 1.0, 0.35585);
        return a;
    }

    const Whitepoint& Whitepoint::B() {
        static const Whitepoint b(0.99072, 1.0, 0.85223);
        return b;
    }

    const Whitepoint& Whitepoint::C() {
        static const Whitepoint c(0.98074, 1.0, 1.18232);
        return c;
    }

    const Whitepoint& Whitepoint::D50() {
        static const Whitepoint d50(0.96422, 1.0, 0.82521);
        return d50;
    }

    const Whitepoint& Whitepoint::D55() {
        static const Whitepoint d55(0.95682, 1.0, 0.92149);
        return d55;
    }

    const Whitepoint& Whitepoint::D65() {
        static const Whitepoint d65(0.95047, 1.0, 1.08883);
        return d65;
    }

    const Whitepoint& Whitepoint::D75() {
        static const Whitepoint d75(0.94972, 1.0, 1.22638);
        return d75;
    }

    const Whitepoint& Whitepoint::E() {
        static const Whitepoint e(1.00000, 1.0, 1.00000);
        return e;
    }

    const Whitepoint& Whitepoint::F2() {
        static const Whitepoint f2(0.99186, 1.0, 0.67393);
        return f2;
    }

    const Whitepoint& Whitepoint::F7() {
        static const Whitepoint f7(0.95041, 1.0, 1.08747);
        return f7;
    }

    const Whitepoint& Whitepoint::F11() {
        static const Whitepoint f11(1.00962, 1.0, 0.64350);
        return f11;
    }
}
