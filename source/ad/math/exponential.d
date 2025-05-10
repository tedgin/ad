/// It extends `std.math.exponential` to support `GDN` objects.
module ad.math.exponential;

static import std.math.exponential;

import std.math: E, isInfinity, isNaN, LN2, sgn;
import std.traits: isIntegral;

static import ad.math.internal;

import ad.core;
import ad.math.internal: areAll, asReal, CommonGDN, dirac, floor, isGDN, isGDNOrReal, isOne, sgn;


/**
 * This function raises $(MATH e) to a given power.
 *
 * If $(MATH f(x) = e$(SUP g(x))), then $(MATH f' = e$(SUP g)g').
 *
 * Params:
 *   Deg = the degree of the GDN
 *   g = the power $(MATH e) is raised to.
 *
 * Returns:
 *   A GDN object representing $(MATH e) raised to the power of `g`.
 */
nothrow pure @nogc @safe GDN!Deg exp(ulong Deg)(in GDN!Deg g)
{
    static if (Deg == 1)
        alias exp_fn = std.math.exponential.exp;
    else
        alias exp_fn = exp;

    const f_red = exp_fn(g.reduce());
    return GDN!Deg(asReal(f_red), f_red * g.d);
}

///
unittest
{
    assert(exp(GDN!1(0, 3)) is GDN!1(1, 3));
}

unittest
{
    import std.format;

    assert(exp(GDN!1.nan) is GDN!1.nan);
    assert(exp(GDN!1.zero) is GDN!1.one);

    const e = exp(-GDN!2.infinity);
    assert(e is GDN!2(0, -0., 0), format("exp(-inf) != %s", e));
}


/**
 * This function computes the base-2 exponential of a given `GDN` object `g`.
 *
 * $(MATH f(x) = 2$(SUP g(x))), then $(MATH f' = 2$(SUP g)g'ln(2))
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` object for which the base-2 exponential is to be calculated.
 *
 * Returns:
 *   A `GDN` object representing the base-2 exponential of `g`.
 */
pure nothrow @nogc @safe GDN!Deg exp2(ulong Deg)(in GDN!Deg g)
{
    static if (Deg == 1)
        alias exp2_fn = std.math.exponential.exp2;
    else
        alias exp2_fn = exp2;

    const f_red = exp2_fn(g.reduce());
    return GDN!Deg(asReal(f_red), f_red * g.d * LN2);
}

///
unittest
{
    assert(exp2(GDN!1(0)) is GDN!1(1, LN2));
}

unittest
{
    assert(exp2(GDN!2(2)) is GDN!2(4, 4*LN2, 4*LN2^^2));
    // f = 2^2 = 4
    // <f',f"> = 2^<2,1> * <1,0> * ln(2)
    //         = <4,4*ln(2)> * <ln(2),0>
    //         = <4*ln(2), 4*ln(2)^2>
}


/**
 * Calculates the value of $(MATH e$(SUP g) - 1).
 *
 * If $(MATH f(x) = e$(SUP g(x)) - 1), then $(MATH f' = e$(SUP g)g').
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the exponent
 *
 * Returns:
 *   A `GDN` object representing $(MATH e$(SUP g) - 1).
 */
pure nothrow @nogc @safe GDN!Deg expm1(ulong Deg)(in GDN!Deg g)
{
    static if (Deg == 1)
        alias exp_fn = std.math.exponential.exp;
    else
        alias exp_fn = exp;

    return GDN!Deg(std.math.exponential.expm1(g.val), exp_fn(g.reduce()) * g.d);
}

///
unittest
{
    assert(expm1(GDN!1(0)) is GDN!1(0, 1));
}

unittest
{
    assert(expm1(GDN!2(1)) is GDN!2(E-1, E, E));
}


/**
 * Separates a `GDN` into its significand and exponent.
 *
 * If $(MATH g(x) = f(x)2$(SUP e(x))) where $(MATH ½ ≤ f < 1) and $MATH(e : ℝ ↦ ℤ), then
 * $(MATH f = g⋅2$(SUP -e)) and $(MATH f' = g'⋅2$(SUP -e) - ln(2)g⋅2$(SUP -e)e'). $(MATH e) is
 * constant when $(MATH f ≠ ½). When $(MATH f = ½), $(MATH e) jumps by $(MATH sgn(g')). This means
 * that $(MATH e' = sgn(g')𝛿(g⋅2$(SUP -e) - ½)). Thus
 * $(MATH f' = g'⋅2$(SUP -e) - ln(2)g⋅2$(SUP -e)sgn(g')𝛿(g⋅2$(SUP -e) - ½)), which simplifies to
 * $(MATH f' = (g' - ln(2)sgn(g')g⋅𝛿(g⋅2$(SUP -e) - ½))2$(SUP -e)).
 *
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` object to be separated.
 *   e = the exponent of `g`
 *
 * Returns:
 *   A `GDN` object representing the significand of `g`.
 *
 * Special cases:
 *   - If `g` is ±0, then `f` is ±0, `f.d` is NaN, and `e` is 0.
 *   - If `g` is +∞, then `f` is +∞, `f.d` is NaN, and `e` is `int.max`.
 *   - If `g` is -∞, then `f` is -∞, `f.d` is NaN, and `e` is `int.min`.
 *   - If `g` is ±NaN, then `f` is ±NaN, `f.d` is NaN, and `e` is `int.min`.
 */
pure nothrow @nogc @safe GDN!Deg frexp(ulong Deg)(in GDN!Deg g, out int e)
{
    if (isNaN(g.val)) {
        e = int.min;
        return GDN!Deg(g.val, GDN!Deg.mkNaNDeriv);
    }

    if (isInfinity(g.val)) {
        e = sgn(g) == 1 ? int.max : int.min;
        return GDN!Deg(g.val, GDN!Deg.mkNaNDeriv);
    }

    if (g == 0) {
        e = 0;
        return GDN!Deg(g.val, GDN!Deg.mkNaNDeriv);
    }

    const f_val = std.math.exponential.frexp(g.val, e);
    const _2r_e = 2.0L^^(-e);

    auto df = g.d * _2r_e;
    if (f_val == 0.5 && g.d != 0) {
        const g_red_scale = g.reduce() * _2r_e;
        df = df - LN2 * sgn(g.d) * g_red_scale * dirac(g_red_scale - 0.5);
    }

    return GDN!Deg(f_val, df);
}

///
unittest
{
    int e;

    const f = frexp(GDN!1(1.1), e);
    assert(f is GDN!1(0.55, 0.5) && e == 1);

    const g = frexp(GDN!1(+0.), e);
    assert(g is GDN!1(0, real.nan) && e == 0);

    const h = frexp(GDN!1(-0.), e);
    assert(h is GDN!1(-0., real.nan) && e == 0);

    const i = frexp(GDN!1(real.infinity), e);
    assert(i is GDN!1(real.infinity, real.nan) && e == int.max);

    const j = frexp(GDN!1(-real.infinity), e);
    assert(j is GDN!1(-real.infinity, real.nan) && e == int.min);

    const k = frexp(GDN!1(real.nan), e);
    assert(k is GDN!1(real.nan, real.nan) && e == int.min);
}

//  TODO: test frexp
unittest
{
    import std.format: format;

    int e;
    const q = frexp(GDN!1(-real.nan), e);
    assert(q is GDN!1(-real.nan, real.nan) && e == int.min);

    const w = frexp(GDN!1(0.5, 0),e);
    assert(w is GDN!1(0.5, 0) && e == 0);

    const r = frexp(GDN!1(0.5, 1), e);
    // f' = 1 * 2^0 - ln(2) * sgn(1) * 0.5 * 2^0 * dirac(0.5 * 2^0 - 0.5)
    //    = 1       - ln(2)          * 0.5       * dirac(0.5       - 0.5)
    //    = 1       - ln(2)/2                    * inf
    //    = -inf
    assert(r is GDN!1(0.5, -real.infinity) && e == 0);

    const t = frexp(GDN!2(3), e);
    // <f',f"> = <1,0> * 2^(-2) = <0.25, 0>
    assert(t is GDN!2(0.75, 0.25, 0) && e == 2, format("frexp(GDN!2(1), %s) != %s", e, t));
}


/**
 * Extracts the exponent of g as a signed integral value.
 *
 * Params:
 *   Deg = the degree of g
 *   g = the GDN to find the exponent of
 *
 * Returns:
 *   the integral exponent of g
 */
pure nothrow @nogc @safe int ilogb(ulong Deg)(in GDN!Deg g)
{
    return std.math.exponential.ilogb(g.val);
}

///
unittest
{
    assert(ilogb(GDN!1.one) == 0);
}


// TODO: implement ldexp


/+ TODO: implement log
/**
This function computes the natural logarithm of its argument. It is analogous
to std.math.log().

Params:
    x = the argument
*/
nothrow pure @safe real log(in real x)
{
    return std.math.log(x);
}

/// ditto
nothrow pure @safe GDN!Deg log(ulong Deg)(in GDN!Deg x)
{
    return x.log();
}
+/


// TODO: implement log10
// TODO: implement log1p


/**
 * This function computes the base-2 logarithm of the given `GDN` object `g`. It is analogous to
 * `std.math.exponential.log2()`.
 *
 * If $(MATH f(x) = lg(g(x))), then $(MATH f' = g'/[g⋅ln(2)])
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` object for which the base-2 logarithm is to be calculated.
 *
 * Returns:
 *   A `GDN` object representing the base-2 logarithm of `g`.
 */
pure nothrow @nogc @safe GDN!Deg log2(ulong Deg)(in GDN!Deg g)
{
    return ad.math.internal.log2(g);
}

///
unittest
{
    import std.math: LN2;

    assert(log2(GDN!1(2)) is GDN!1(1, 1/(2 * LN2)));
    assert(log2(GDN!1(-0.)) is GDN!1(-real.infinity, real.nan));
    assert(log2(GDN!1(+0.)) is GDN!1(-real.infinity, real.infinity));
}



// TODO: implement logb


// TODO: implement pow
/**
 * This function determines the value of a `GDN` raised to an integer power.
 *
 * If $(MATH f(x) = g(x)$(SUP n)), where $(MATH n ∈ ℤ), then if $(MATH n = 0), $(MATH f' = 0g'),
 * otherwise $(MATH f' = ng$(SUP n-1)g').
 *
 * Params:
 *   I = the integer type of the exponent
 *   Deg = the degree of the `GDN` `g`
 *   g = the GDN base
 *   n = the exponent
 *
 * Returns:
 *   It returns a `GDN` representing `g` raised of `n`.
 */
pure nothrow @nogc @safe GDN!Deg pow(I, ulong Deg)(in GDN!Deg g, in I n) if (isIntegral!I)
{
    return ad.math.internal.pow(g, n);
}

///
unittest
{
    assert(pow(GDN!1(2), 3) is GDN!1(8, 12));
}

// TODO: implement
// pure nothrow @nogc @safe GDN!Deg pow(I, ulong Deg)(in I x, in GDN!Deg g) if (isIntegral!I);

// TODO: implement
// pure nothrow @nogc @safe
// CommonGDN!(G, H) pow(G, H)(G g, H h) if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H));


// TODO: implement powmod
// TODO: implement scalbn