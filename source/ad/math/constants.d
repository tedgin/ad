/// This module extends the `std.math.constants` to support `GDN` objects.
module ad.math.constants;

static import std.math.constants;

import ad.core;

/// Euler's constant $(MATH e)
enum GDN!Deg E(ulong Deg) = GDN!Deg.mkConst(std.math.constants.E);

/// $(MATH π)
enum GDN!Deg PI(ulong Deg) = GDN!Deg.mkConst(std.math.constants.PI);

/// $(MATH π/2)
enum GDN!Deg PI_2(ulong Deg) = GDN!Deg.mkConst(std.math.constants.PI_2);

/// $(MATH π/4)
enum GDN!Deg PI_4(ulong Deg) = GDN!Deg.mkConst(std.math.constants.PI_4);

/// $(MATH 1/π)
enum GDN!Deg M_1_PI(ulong Deg) = GDN!Deg.mkConst(std.math.constants.M_1_PI);

/// $(MATH 2/π)
enum GDN!Deg M_2_PI(ulong Deg) = GDN!Deg.mkConst(std.math.constants.M_2_PI);

/// $(MATH 2/√π)
enum GDN!Deg M_2_SQRTPI(ulong Deg) = GDN!Deg.mkConst(std.math.constants.M_2_SQRTPI);

/// $(MATH ln(10))
enum GDN!Deg LN10(ulong Deg) = GDN!Deg.mkConst(std.math.constants.LN10);

/// $(MATH ln(2))
enum GDN!Deg LN2(ulong Deg) = GDN!Deg.mkConst(std.math.constants.LN2);

/// $(MATH log(2))
enum GDN!Deg LOG2(ulong Deg) = GDN!Deg.mkConst(std.math.constants.LOG2);

/// $(MATH lg(e))
enum GDN!Deg LOG2E(ulong Deg) = GDN!Deg.mkConst(std.math.constants.LOG2E);

/// $(MATH lg(10))
enum GDN!Deg LOG2T(ulong Deg) = GDN!Deg.mkConst(std.math.constants.LOG2T);

/// $(MATH log(e))
enum GDN!Deg LOG10E(ulong Deg) = GDN!Deg.mkConst(std.math.constants.LOG10E);

/// $(MATH √2)
enum GDN!Deg SQRT2(ulong Deg) = GDN!Deg.mkConst(std.math.constants.SQRT2);

/// $(MATH √½)
enum GDN!Deg SQRT1_2(ulong Deg) = GDN!Deg.mkConst(std.math.constants.SQRT1_2);

unittest
{
    assert(E!1 == std.math.constants.E, "E broken");
    assert(PI!1 == std.math.constants.PI, "PI broken");
    assert(PI_2!1 == std.math.constants.PI_2, "PI_2 broken");
    assert(PI_4!1 == std.math.constants.PI_4, "PI_4 broken");
    assert(M_1_PI!1 == std.math.constants.M_1_PI, "M_1_PI broken");
    assert(M_2_PI!1 == std.math.constants.M_2_PI, "M_2_PI broken");
    assert(M_2_SQRTPI!1 == std.math.constants.M_2_SQRTPI, "SQRT1_2 broken");
    assert(LN10!2 == std.math.constants.LN10, "LN10 broken");
    assert(LN2!1 == std.math.constants.LN2, "LN2 broken");
    assert(LOG2!1 == std.math.constants.LOG2, "LOG2 broken");
    assert(LOG2E!1 == std.math.constants.LOG2E, "LOG2E broken");
    assert(LOG2T!1 == std.math.constants.LOG2T, "LOG2T broken");
    assert(LOG10E!1 == std.math.constants.LOG10E, "LOG10E broken");
    assert(SQRT2!1 == std.math.constants.SQRT2, "SQRT2 broken");
    assert(SQRT1_2!1 == std.math.constants.SQRT1_2, "SQRT1_2 broken");
}
