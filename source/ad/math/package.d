/** TODO: finish documentation after all of the submodules are completed.
 * This module extends the `core.math` and `std.math` libraries to support `GDN` objects. It is
 * decomposed into submodules in the same way that std.math is. It also exports all of the symbols
 * from `ad.math`, `core.math` and `std.math` to make it easier to work with real and generalized
 * dual numbers together.
 */
module ad.math;

public import core.math: toPrec, yl2x, yl2xp1;
public import std.math;
public import ad.core;
public import ad.math.algebraic;
public import ad.math.constants;
public import ad.math.exponential;
public import ad.math.operations;
public import ad.math.remainder;
public import ad.math.rounding;
public import ad.math.traits;
public import ad.math.trigonometry;

static import core.math;

import std.algorithm: min;
import std.math.constants: LN2;
import std.traits: isFloatingPoint, Select;

import ad.math.internal: areAll, asGDN, dirac, CommonGDN, isGDN, isGDNOrReal, isOne, signbit;


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
 * This function rounds the value of a `GDN` to a given floating point type removing all derivative
 * information.
 *
 * Params:
 *   F = the float point type to be converted to
 *   Deg = the degree of `g`
 *   g = the generalized dual number to be converted
 *
 * Returns:
 *  the rounded valued with precision determined by `F`.
 *
 */
pragma(inline, true) pure nothrow @nogc @safe
F toPrec(F, ulong Deg)(in GDN!Deg g) if (isFloatingPoint!F)
{
    return core.math.toPrec!F(g.val);
}

///
unittest
{
    assert(typeid(toPrec!float(GDN!1.zero)) == typeid(float));
}


/**
 * This function computes $(MATH h⋅lg(g)). It either `g` or `h` has type `real`, it is converted to
 * a constant generalized dual number with the same degree as the other parameter.
 *
 * If $(MATH f(x) = h(x)lg(g(x))), then $(MATH f' = h'lg(g) + hg'/(ln(2)g))
 *
 * Params:
 *   G = the type of `g`
 *   H = the type of `h`
 *   g = the argument of logarithm
 *   h = the multiplier of the logarithm
 *
 * Returns:
 *   The resulting generalized dual number will have a degree equal to the lesser of the degrees of
 *   `g` and `h`.
 */
pragma(inline, true) pure nothrow @nogc @safe
CommonGDN!(G, H) yl2x(G, H)(in G g, in H h) if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H))
{
    alias Deg = typeof(return).DEGREE;
    return yl2x_impl(asGDN!Deg(g), asGDN!Deg(h));
}

///
unittest
{
    import std.math: LN2;

    assert(yl2x(GDN!1(2), GDN!1(3)) is GDN!1(3, 1 + 1.5/LN2));
    assert(yl2x(GDN!1(+0., -1), GDN!1(1)) is GDN!1(-real.infinity,  -real.infinity));
    assert(typeof(yl2x(GDN!2(1), GDN!1(2))).DEGREE == 1);
    assert(yl2x(GDN!1(1), 2.) is GDN!1(0, 2/LN2));
}

unittest
{
    assert(yl2x(1, GDN!1(2)) is GDN!1(0, 0));
}

private pure nothrow @nogc @safe GDN!Deg yl2x_impl(ulong Deg)(in GDN!Deg g, in GDN!Deg h)
{
    alias yl2x_red = Select!(Deg == 1, core.math.yl2x, yl2x_impl);

    GDN!Deg.DerivType!1 df;
    if (signbit(g) == 0) {
        const g_red = g.reduce();
        df = yl2x_red(g_red, h.d) + h.reduce()*g.d/(LN2 * g_red);
    }

    return GDN!Deg(core.math.yl2x(g.val, h.val), df);
}

unittest
{
    assert(isNaN(yl2x(GDN!1(-1), GDN!1(1))));

    const f = yl2x(GDN!1(0), GDN!1(1));
    assert(f == -real.infinity && isNaN(f.d));

    assert(yl2x(GDN!2(1), GDN!2(2)) is GDN!2(0, 2/LN2, 0));
    // f = 0
    // <f',f"> = h'lg(g) + hg'/(ln(2)g)
    //    = <1,0>lg(<1,1>) + <2,1><1,0>/ln(2)<1,1>
    //    = <0,0+1/ln(2)> + <2,1>/<1,1>/ln(2)
    //    = <0,1/ln(2)> + <2,-1>/ln(2)
    //    = <2/ln(2),0>
}

/**
 * Computes $(MATH h⋅lg(g + 1)), for $(MATH -(1 - √½) ≤ x ≤ +(1 - √½)). When $(MATH g) is outside of
 * this interval, the results are undefined. It either `g` or `h` has type `real`, it is converted
 * to a constant generalized dual number with the same degree as the other parameter.
 *
 * If $(MATH f(x) = h(x)lg(g(x) + 1)), then $(MATH f' = h'lg(g + 1) + hg'/[ln(2)(g + 1)])
 *
 * Params:
 *   G = the type of `g`
 *   H = the type of `h`
 *   g = the argument of logarithm
 *   h = the multiplier of the logarithm
 *
 * Returns:
 *   The resulting generalized dual number will have a degree equal to the lesser of the degree of
 *   `g` and `h`.
 */
pragma(inline, true) pure nothrow @nogc @safe
CommonGDN!(G, H) yl2xp1(G, H)(in G g, in H h) if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H))
{
    alias Deg = typeof(return).DEGREE;
    return yl2xp1_impl(asGDN!Deg(g), asGDN!Deg(h));
}

///
unittest
{
    import std.math: LN2;

    assert(yl2xp1(GDN!1(0), GDN!1(3)) is GDN!1(0, 3/LN2));
    assert(typeof(yl2xp1(GDN!2(0), GDN!1(1))).DEGREE == 1);
    assert(yl2xp1(GDN!1(0), 1) is GDN!1(0, 1/LN2));
}

unittest
{
    assert(yl2xp1(0, GDN!1(1)) is GDN!1(0, 0));
}

private pure nothrow @nogc @safe GDN!Deg yl2xp1_impl(ulong Deg)(in GDN!Deg g, in GDN!Deg h)
{
    alias yl2xp1_red = Select!(Deg == 1, core.math.yl2xp1, yl2xp1_impl);

    GDN!Deg.DerivType!1 df;
    if (g > -1) {
        const g_red = g.reduce();
        df = yl2xp1_red(g_red, h.d) + h.reduce()*g.d/(LN2 * (g_red + 1));
    }

    return GDN!Deg(core.math.yl2xp1(g.val, h.val), df);
}

unittest
{
    import std.math: isNaN, LN2;

    assert(yl2xp1_impl(GDN!2(0), GDN!2(1)) is GDN!2(0, 1/LN2, 1/LN2));
    // f = 0
    // <f',f"> = h'lg(g+1) + hg'/[ln(2)(g+1)]
    //    = <1,0>lg(<0,1>+1) + <1,1><1,0>/[ln(2)(<0,1>+1)]
    //    = <1,0>lg<1,1> + <1,1>/[ln(2)<1,1>]
    //    = <0,1/ln(2)> + <1,0>/ln(2)
    //    = <1/ln(2),1/ln(2)>

    assert(isNaN(yl2xp1_impl(GDN!1(-1), GDN!1(0)).d));

    const q = yl2xp1_impl(GDN!1(-2), GDN!1(0));
    assert(isNaN(q.d), "yl2xp1(2,0) should not have a derivative");
}
