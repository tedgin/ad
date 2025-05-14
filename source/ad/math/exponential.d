/// It extends `std.math.exponential` to support `GDN` objects.
module ad.math.exponential;

static import core.math;
static import std.math.exponential;

import std.math: E, isInfinity, isNaN, LN10, LN2, sgn;
import std.traits: isIntegral, Select;

static import ad.math.internal;

import ad.core;
import ad.math.internal:
    areAll, areNone, asReal, CommonGDN, dirac, floor, isGDN, isGDNOrReal, isInfinity, isNaN, isOne,
    sgn, signbit;


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
    std.math.exponential.frexp(g.val, e);

    if (g == 0 || isInfinity(g) || isNaN(g)) {
        return GDN!Deg(g.val, GDN!Deg.mkNaNDeriv());
    }

    return g * 2.0L ^^ -e;
}

///
unittest
{
    int e;
    const f = frexp(GDN!1(1.1), e);
    assert(f is GDN!1(0.55, 0.5) && e == 1);
}

unittest
{
    import std.format: format;

    int e;

    const q = frexp(GDN!1(+0.), e);
    assert(q == +0. && isNaN(q.d) && e == 0);

    const w = frexp(GDN!1(-0.), e);
    assert(w == -0. && isNaN(w.d) && e == 0);

    const r = frexp(GDN!1(+real.infinity), e);
    assert(r == +real.infinity && isNaN(r.d) && e == int.max);

    const t = frexp(GDN!1(-real.infinity), e);
    assert(t == -real.infinity && isNaN(r.d) && e == int.min);

    const y = frexp(GDN!1(real.nan), e);
    assert(signbit(y) == 0 && isNaN(y) && isNaN(y.d) && e == int.min);

    const u = frexp(GDN!1(-real.nan), e);
    assert(u is GDN!1(-real.nan, real.nan) && e == int.min);

    const i = frexp(GDN!2(3), e);
    assert(i is GDN!2(0.75, 0.25, 0) && e == 2, format("frexp(GDN!2(1), %s) != %s", e, i));
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


/**
 * This function computes $(MATH 2$(SUP c)g).
 *
 * If $(MATH f(x) = 2$(SUP c)g(x)), then $(MATH f' = 2$(SUP c)g').
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the generalized dual number being scaled.
 *   c = the power of $(MATH 2) used to scale `g`,
 *
 * Returns:
 *   A `GDN` object resulting from the computation.
 */
pragma(inline, true) pure nothrow @nogc @safe GDN!Deg ldexp(ulong Deg)(in GDN!Deg g, in int c)
{
    alias ldexp_red = Select!(Deg == 1, core.math.ldexp, ldexp);
    return GDN!Deg(core.math.ldexp(g.val, c), ldexp_red(g.d, c));
}

///
unittest
{
    assert(ldexp(GDN!2(1), 2) is GDN!2(4, 4, 0));
}


/**
 * This function computes the natural logarithm of its argument.
 *
 * If $(MATH f(x) = ln(g(x))), then $(MATH f' = g'/g).
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the argument
 *
 * Returns:
 *   the natural logarithm of `g`.
 */
pure nothrow @nogc @safe GDN!Deg log(ulong Deg)(in GDN!Deg g)
{
    return g.log();
}

///
unittest
{
    assert(log(GDN!2(1)) is GDN!2(0, 1, -1));
}


/**
 * This function computes the logarithm of the given `GDN` object `g`.
 *
 * If $(MATH f(x) = log(g(x))), then $(MATH f' = g'/[g⋅ln(10)])
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` object for which the logarithm is to be calculated.
 *
 * Returns:
 *   A `GDN` object representing the logarithm of `g`.
 */
pure nothrow @nogc @safe GDN!Deg log10(ulong Deg)(in GDN!Deg g)
{
    const df = signbit(g) == 1 ? GDN!Deg.mkNaNDeriv : g.d / (LN10*g.reduce());
    return GDN!Deg(std.math.exponential.log10(g.val), df);
}

///
unittest
{
    assert(log10(GDN!1(1)) is GDN!1(0, 1/LN10));
}

unittest
{
    import std.format: format;
    import std.math: LOG2;

    const q = log10(GDN!1(-1));
    assert(isNaN(q) && isNaN(q.d), format("log10(-1) = %s", q));

    assert(log10(GDN!2(2)) is GDN!2(LOG2, 0.5/LN10, -0.25/LN10));
    // f = log(2)
    // <f',f"> = <1,0>/(ln(10)<2,1>)
    //         = <1,0>/<2,1>/ln(10)
    //         = <0.5,-0.25>/ln(10)
    //         = <1/[2ln(10)],-1/[4ln(10)]>
}


// TODO: implement log1p
/**
 * Calculates the natural logarithm of 1 + g.
 *
 * If $(MATH f(x) = ln(1 + g(x))), then $(MATH f' = g'/(1 + g)).
 *
 * Params:
 *   Deg = the degree of g
 *   g = the value to compute the shifted logarithm of
 *
 * Returns:
 *  $(MATH ln(1 + g)) as a `GDN`
 */
pure nothrow @nogc @safe GDN!Deg log1p(ulong Deg)(in GDN!Deg g)
{
    const df = g <= -1 ? GDN!Deg.mkNaNDeriv() : g.d/(1+g.reduce());
    return GDN!Deg(std.math.exponential.log1p(g.val), df);
}

///
unittest
{
    import std.math: log;

    assert(log1p(GDN!1(3)) is GDN!1(log(4), 0.25));
    assert(log1p(GDN!1(-1)) is GDN!1(-real.infinity, real.nan));
}

unittest
{
    assert(isNaN(log1p(GDN!1(-2))));
}


/**
 * This function computes the base-2 logarithm of the given `GDN` object `g`.
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


/**
 * Extracts the exponent of g as a signed integral valued real
 *
 * Params:
 *   Deg = the degree of g
 *   g = the `GDN` to extract the exponent from its value
 *
 * Returns:
 *   The exponent of g
 */
nothrow @nogc @safe real logb(ulong Deg)(in GDN!Deg g)
{
    return std.math.exponential.logb(g.val);
}

///
unittest
{
    assert(logb(GDN!1(1)) == 0);
}


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
 *   It returns a `GDN` representing `g` raised to `n`.
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

/**
 * This function determine the value of a integer raised to a `GDN` power;
 *
 * If $(MATH f(x) = n$(SUP g(x))), then $(MATH f' = n$(SUP g)g'ln(n)).
 *
 * Params:
 *   I = the integer type of the base
 *   Deg = the degree of g
 *   n = the base
 *   g = the exponent
 *
 * Returns:
 *   A `GDN` representing `n` raised to `g`.
 */
pure nothrow @nogc @safe GDN!Deg pow(I, ulong Deg)(in I n, in GDN!Deg g) if (isIntegral!I)
{
    static if (Deg == 1) {
        const f_val = std.math.exponential.pow(n, g.val);
        return GDN!Deg(f_val, f_val * g.d * std.math.exponential.log(n));
    } else {
        const f_red = pow(n, g.reduce());
        return GDN!Deg(f_red.val(), f_red * g.d * std.math.exponential.log(n));
    }
}

///
unittest
{
    import std.math: log;

    assert(pow(2, GDN!1(5)) is GDN!1(32, 32*log(2)));
}

unittest
{
    const ln3 = std.math.exponential.log(3);
    assert(pow(3, GDN!2(1, 2, 4)) is GDN!2(3, 6*ln3, 12*ln3*(ln3 + 1)));
    // f = 3
    // <f',f"> = 3^<1,2> * <2,4> * ln(3)
    //         = <3,3*2ln(3)> * <2,4> * ln(3)
    //         = <1,2ln(3)> * <2,4> * 3ln(3)
    //         = <2,4ln(3) + 4> * 3ln(3)
    //         = <6ln(3), 12ln(3)(ln(3)+1)>
}

/**
 * This function calculates $(MATH g$(SUP h)). It is the same as `g ^^ h`.
 *
 * Params:
 *   G = the type of g, either real or a GDN
 *   H = the type of h, either real or a GDN
 *   g = the base
 *   h = the exponent
 *
 * Returns:
 *   It returns `g ^^ h`.
 */
pragma(inline, true) pure nothrow @nogc @safe
CommonGDN!(G, H) pow(G, H)(in G g, in H h)
if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H) && areNone!(isIntegral, G, H))
{
    return g^^h;
}

unittest
{
    assert(pow(GDN!1(2), 3.) is GDN!1(8, 12));
}

// TODO: implement powmod
// TODO: implement scalbn