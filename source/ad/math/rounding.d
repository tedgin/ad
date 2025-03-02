/// It extends `std.math.rounding` to support `GDN` objects.
module ad.math.rounding;

static import std.math.rounding;

import std.math: isFinite;

import ad.core;
import ad.math.internal: dirac;
import ad.math.traits: isFinite;

// TODO: finish implement this module.


/**
 * Returns the value of `g` rounded upward to the nearest integer.
 *
 * If $(MATH f(x) = ⌈g(x)⌉), then $(MATH f' = g'∑$(SUB i∊ℤ)𝛿(g-i))
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` to round
 *
 * Returns:
 *   the rounded `GDN`
 */
pure nothrow @nogc @safe GDN!Deg ceil(ulong Deg)(in GDN!Deg g)
{
    static if (Deg == 1) {
        const f_red = std.math.rounding.ceil(g.reduce());
        const f_val = f_red;
    } else {
        const f_red = ceil(g.reduce());
        const f_val = f_red.val;
    }

    GDN!Deg.DerivType!1 df;
    if (isFinite(g)) {
        df = g.d * dirac(g.reduce() - f_red);
    }

    return GDN!Deg(f_val, df);
}

///
unittest
{
    assert(ceil(GDN!1(1)) is GDN!1(1, real.infinity));
    assert(ceil(GDN!1(-1.4)) is GDN!1(-1, 0));
}

unittest
{
    import std.math: LN2;

    assert(ceil(GDN!1.nan) is GDN!1.nan);
    assert(ceil(GDN!1(real.infinity)) is GDN!1(real.infinity, real.nan));
    assert(ceil(GDN!1(-real.infinity)) is GDN!1(-real.infinity, real.nan));
    assert(ceil(GDN!2(2.3)) is GDN!2(3, 0, 0));
    assert(ceil(GDN!1(0, -1/LN2)) is GDN!1(0, -real.infinity));
}

/**
 * Returns the value of `g` rounded downward to the nearest integer.
 *
 * If $(MATH f(x) = ⌊g(x)⌋), then $(MATH f' = g'∑$(SUB i∊ℤ)𝛿(g-i))
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` to round
 *
 * Returns:
 *   the rounded `GDN`
 */
pure nothrow @nogc @safe GDN!Deg floor(ulong Deg)(in GDN!Deg g)
{
    static if (Deg == 1) {
        const f_red = std.math.rounding.floor(g.reduce());
        const f_val = f_red;
     } else {
        const f_red = floor(g.reduce());
        const f_val = f_red.val;
     }

    GDN!Deg.DerivType!1 df;
    if (isFinite(g)) {
        df = g.d * dirac(g.reduce() - f_red);
    }

    return GDN!Deg(f_val, df);
}

///
unittest
{
    assert(floor(GDN!1(1)) is GDN!1(1, real.infinity));
    assert(floor(GDN!1(-1.4)) is GDN!1(-2, 0));
}

unittest
{
    import std.math: LN2;

    assert(floor(GDN!1.nan) is GDN!1.nan);
    assert(floor(GDN!1(real.infinity)) is GDN!1(real.infinity, real.nan));
    assert(floor(GDN!1(-real.infinity)) is GDN!1(-real.infinity, real.nan));
    assert(floor(GDN!2(2.3)) is GDN!2(2, 0, 0));
    assert(floor(GDN!1(0, -1/LN2)) is GDN!1(0, -real.infinity));
}
