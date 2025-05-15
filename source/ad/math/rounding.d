/// It extends `std.math.rounding` to support `GDN` objects.
module ad.math.rounding;

static import core.math;
static import std.math.rounding;

import std.math: abs, FloatingPointControl, isFinite, isInfinity, isNaN, nearbyint, signbit;
import std.traits: arity, isIntegral, Parameters, ReturnType;

static import ad.math.internal;

import ad.core;
import ad.math.internal:
    areAll, asGDN, asReal, CommonGDN, dirac, isGDN, isGDNOrReal, isOne, nextDown, nextUp, pow;


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
    return ad.math.internal.ceil(g);
}

///
unittest
{
    assert(ceil(GDN!1(1)) is GDN!1(1, real.infinity));
    assert(ceil(GDN!1(-1.4)) is GDN!1(-1, 0));
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
    return ad.math.internal.floor(g);
}

///
unittest
{
    assert(floor(GDN!1(1)) is GDN!1(1, real.infinity));
    assert(floor(GDN!1(-1.4)) is GDN!1(-2, 0));
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
    if (isInfinity(g.val)) {
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
 * If $(MATH f(x) = nearbyint(g(x))), then $(MATH f' = (df/dg)g'), where $(MATH df/dg = 𝛿(g - m)),
 * $(MATH 𝛿) is the Dirac delta function, and $(MATH m) is a rounding mode split point.
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
 * If $(MATH f(x) = rint(g(x))), then $(MATH f' = (df/dg)g'), where $(MATH df/dg = 𝛿(g - m)),
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


/**
 * Round val to a multiple of unit. rfunc specifies the rounding function to use.
 *
 * Params:
 *   V = the value type
 *   U = the unit type
 *   round = the rounding function to use
 *   val = the value to round
 *   unit = the unit to round to
 *
 * Returns:
 *   The rounded value of `val` to the nearest multiple of `unit`.
 */
CommonGDN!(U, V) quantize(alias round=rint, U, V)(in V val, in U unit)
if (is(typeof(round(CommonGDN!(U, V).init)) : CommonGDN!(U, V)))
{
    return quantize_impl!round(val, unit);
}

///
unittest
{
    const q = quantize!nearbyint(GDN!1(5), GDN!1(3));
    assert(q is GDN!1(6, 2));
}

/**
 * Round `g` to a multiple of `base ^^ exp`. `rfunc` specifies the rounding function to use.
 *
 * Params:
 *   G = the type of `g`
 *   I = the exponent type
 *   round = the rounding function to use
 *   base = the base of the number to round to
 *   g = the `GDN` object to round
 *   exp = the exponent of the number to round to
 *
 * Returns:
 *   The rounded `GDN` object.
 */
CommonGDN!(G, typeof(base)) quantize(alias base, alias round=rint, G, I)(in G g, in I exp)
if (isOne!(isGDN, G, typeof(base))
    && areAll!(isGDNOrReal, G, typeof(base))
    && (is(typeof(round(G.init)) : G) || is(typeof(round(typeof(base).init)) : typeof(base)))
    && isIntegral!I)
{
    alias Deg = typeof(return).DEGREE;
    enum b = asGDN!Deg(base);
    return quantize_impl!round(asGDN!Deg(g), pow(b, exp));
}

/// ditto
CommonGDN!(G, typeof(base)) quantize(alias base, long exp=1, alias round=rint, G)(in G g)
if (isOne!(isGDN, G, typeof(base))
    && areAll!(isGDNOrReal, G, typeof(base))
    && (is(typeof(round(G.init)) : G) || is(typeof(round(typeof(base).init)) : typeof(base))))
{
    alias Deg = typeof(return).DEGREE;
    enum unit = pow(asGDN!Deg(base), exp);
    return quantize_impl!round(asGDN!Deg(g), unit);
}

///
unittest
{
    import ad.math.operations: isClose;

    const f = quantize!10(GDN!1(345.678_9), -2);
    assert(isClose(f, 345.68) && f.d == 0);

    const g = quantize!(GDN!1(2))(GDN!1(1.6), -1);
    assert(g is GDN!1(1.5, -0.75));

    assert(quantize!22(GDN!1(12_345.678_9)) is GDN!1(12_342, 0));
}

unittest
{
}

private pragma(inline, true)
GDN!Deg quantize_impl(alias round, ulong Deg)(in GDN!Deg val, in GDN!Deg unit)
if (is(typeof(round(GDN!Deg.init)) : GDN!Deg))
{
    return round(val / unit) * unit;
}

unittest
{
    const f = quantize_impl!rint(GDN!1(1.5), GDN!1(0.5, 0));
    assert(f is GDN!1(1.5, 0));
}


/**
 * This function rounds `g` to a `long` using the current rounding mode. All of the derivative terms
 * are lost.
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` object to be rounded.
 *
 * Returns:
 *   the rounded value of `g`.
 */
pragma(inline, true) pure nothrow @nogc @safe long rndtol(ulong Deg)(in GDN!Deg g)
{
    return core.math.rndtol(g.val);
}

///
unittest
{
    assert(rndtol(GDN!1(1.1)) == 1L);
}


/**
 * Returns a value g rounded to the nearest integer.
 *
 * If $(MATH f(x) = round(g(x))), then $(MATH f' = g'∑$(SUB i∊ℤ)𝛿(g-i-½),
 *
 * Params:
 *   Deg = the degree of g
 *   g = the `GDN` to round.
 *
 * Returns:
 *   The rounded `GDN`,
 */
nothrow @nogc @trusted GDN!Deg round(ulong Deg)(in GDN!Deg g)
{
    static if (Deg == 1)
        const f_red = std.math.rounding.round(g.reduce());
    else
        const f_red = round(g.reduce());

    GDN!Deg.DerivType!1 df;
    if (isFinite(g.val)) {
        df = g.d * dirac(abs(g.reduce() - f_red) - 0.5);
    }

    return GDN!Deg(asReal(f_red), df);
}

///
unittest
{
    assert(round(GDN!1(4.5)) is GDN!1(5, real.infinity));
    assert(round(GDN!1(-4.5)) is GDN!1(-5, real.infinity));
}

unittest
{
    assert(round(GDN!2(5.4)) is GDN!2(5, 0, 0));
}


/**
 * Returns a value g truncated to the integer part.
 *
 * Where $(MATH g(x) < 0), $(MATH f(x) = ⌈g(x)⌉), otherwise $(MATH f(x) = ⌊g(x)⌋).
 *
 * Params:
 *   Deg = the degree of g
 *   g = the `GDN` to round.
 *
 * Returns:
 *   The truncated `GDN`,
 */

pure nothrow @nogc @trusted GDN!Deg trunc(ulong Deg)(in GDN!Deg g)
{
    static if (Deg == 1)
        const f_red = std.math.rounding.trunc(g.reduce());
    else
        const f_red = trunc(g.reduce());

    GDN!Deg.DerivType!1 df;
    if (isFinite(g.val)) {
        df = g.d * dirac(abs(g.reduce() - f_red));
    }

    return GDN!Deg(asReal(f_red), df);
}

///
unittest
{
    assert(trunc(GDN!1(0.01)) is GDN!1(+0., 0));
    assert(trunc(GDN!1(-0.49)) is GDN!1(-0., 0));
}

unittest
{
    assert(trunc(GDN!2(0.5)) is GDN!2(0, 0, 0));
}
