/// It extends `std.math.rounding` to support `GDN` objects.
module ad.math.rounding;

static import core.math;
static import std.math.rounding;

import std.math: FloatingPointControl, isFinite, isNaN, nearbyint, signbit;

import ad.core;
import ad.math.internal: dirac;
import ad.math.operations: nextDown, nextUp;
import ad.math.traits: asReal, isFinite, isInfinity;


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
    static if (Deg == 1)
        const f_red = std.math.rounding.ceil(g.reduce());
    else
        const f_red = ceil(g.reduce());

    GDN!Deg.DerivType!1 df;
    if (isFinite(g)) {
        df = g.d * dirac(g.reduce() - f_red);
    }

    return GDN!Deg(asReal(f_red), df);
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
    static if (Deg == 1)
        const f_red = std.math.rounding.floor(g.reduce());
    else
        const f_red = floor(g.reduce());

    GDN!Deg.DerivType!1 df;
    if (isFinite(g)) {
        df = g.d * dirac(g.reduce() - f_red);
    }

    return GDN!Deg(asReal(f_red), df);
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


/**
 * Rounds `g` to the nearest integer value, using the current rounding mode.
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` object to be rounded.
 *
 * Returns:
 *   An integer representing the rounded value of `g`.
 */
pure nothrow @nogc @safe long lrint(ulong Deg)(in GDN!Deg g)
{
    return std.math.rounding.lrint(g.val);
}

///
unittest
{
    assert(lrint(GDN!1(1.9)) == 2L);
}


/**
 * Returns the value of a `GDN` rounded to the nearest integer.
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` object to be rounded.
 * Returns:
 *   An integer representing the rounded value of `g`.
 */
nothrow @nogc @safe long lround(ulong Deg)(in GDN!Deg g)
{
    return std.math.rounding.lround(g.val);
}

///
unittest
{
    assert(lround(GDN!1(1.5)) == 2L);
    assert(lround(GDN!1(-0.5)) == -1L);
}


private pragma(inline, true) pure nothrow @nogc @safe
GDN!Deg nearbyint_impl(string impl, ulong Deg)(in GDN!Deg g)
{
    mixin("const f = " ~ impl ~ "(g.val);");

    auto dfdg = GDN!Deg.mkZeroDeriv();
    if (isInfinity(g)) {
        dfdg = GDN!Deg.mkNaNDeriv();
    } else {
        auto fn = std.math.rounding.nearbyint(nextDown(g).val);
        auto fp = std.math.rounding.nearbyint(nextUp(g).val);

        if (f == 0) {
            if (signbit(f) == 1) {
                fp = f;
            } else {
                fn = f;
            }
        }

        if (fn != fp) {
            dfdg = dirac(g.reduce() - g.val);
        }
    }

    return GDN!Deg(f, dfdg*g.d);
}

unittest
{
    enum impl = "std.math.rounding.nearbyint";

    assert(nearbyint_impl!impl(GDN!1.infinity) is GDN!1(real.infinity, real.nan));

    const f = nearbyint_impl!impl(GDN!2(1.5));
    assert(f == 2 && f.d == real.infinity && isNaN(f.d!2));

    assert(nearbyint_impl!impl(GDN!1(-0.)) is GDN!1(-0., 0));
    assert(nearbyint_impl!impl(GDN!1(+0.)) is GDN!1(+0., 0));
}


/**
 * Rounds `g` to the nearest integer value, using the current rounding mode.
 *
 * If $(MATH f(g(x)) = nearbyint(g(x))), then $(MATH f' = (df/dg)g'), where
 * $(MATH df/dg = 𝛿(g - m)), $(MATH 𝛿) is the Dirac delta function, and $(MATH m) is a rounding mode
 * split point.
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` object to be rounded.
 *
 * Returns:
 *   A `GDN` object representing the rounded value of `g`.
 */
pure nothrow @nogc @safe GDN!Deg nearbyint(ulong Deg)(in GDN!Deg g)
{
    return nearbyint_impl!"std.math.rounding.nearbyint"(g);
}

///
unittest
{
    const e = nearbyint(GDN!2(1.5));
    assert(e == 2 && e.d == real.infinity && isNaN(e.d!2));
}


/**
 * Rounds `g` to the nearest integer value, using the current rounding mode. If the return value is
 * not identical to `g`, the `FE_INEXACT` exception is raised.
 *
 * If $(MATH f(g(x)) = rint(g(x))), then $(MATH f' = (df/dg)g'), where $(MATH df/dg = 𝛿(g - m)),
 * $(MATH 𝛿) is the Dirac delta function, and $(MATH m) is a rounding mode split point.
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` object to be rounded.
 *
 * Returns:
 *   A `GDN` object representing the rounded value of `g`.
 */
pragma(inline, true) pure nothrow @nogc @safe GDN!Deg rint(ulong Deg)(in GDN!Deg g)
{
    return nearbyint_impl!"std.math.rounding.rint"(g);
}

///
unittest
{
    import std.math: ieeeFlags, resetIeeeFlags;

    resetIeeeFlags();
    const e = rint(GDN!2(1.5));
    assert(ieeeFlags.inexact);
    assert(e == 2 && e.d == real.infinity && isNaN(e.d!2));
}

unittest
{
    import std.math: ieeeFlags, resetIeeeFlags;

    resetIeeeFlags();
    const w = rint(GDN!1.one);
    assert(!ieeeFlags.inexact);
    assert(w is GDN!1.one);
}


// TODO: Implement quantize
// TODO: Implement rndtol
// TODO: Implement round
// TODO: Implement trunc