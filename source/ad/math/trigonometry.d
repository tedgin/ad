/// It extends the `std.math.trigonometry` to support `GDN` objects.
module ad.math.trigonometry;

static import core.math;
static import std.math.trigonometry;

import std.math: isFinite, pow, sqrt;

import ad.core;
import ad.math.internal:
    areAll, asGDN, CommonGDN, isFinite, isGDN, isGDNOrReal, isInfinity, isOne, pow, sqrt;


/**
 * This function computes the sine of its argument.
 *
 * If $(MATH f(x) = sin(g(x))), then $(MATH f' = cos(g)g').
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` to compute the sine of`
 *
 * Returns:
 *   the sine expressed as a `GDN`.
 */
pragma(inline, true) pure nothrow @nogc @safe GDN!Deg sin(ulong Deg)(in GDN!Deg g)
{
    static if (Deg == 1)
        alias cosine = core.math.cos;
    else
        alias cosine = cos;

    return GDN!Deg(core.math.sin(g.val), cosine(g.reduce())*g.d);
}

///
unittest
{
    assert(sin(GDN!2(0)) is GDN!2(0, 1, 0));
}

unittest
{
    import std.math: isClose, PI_2;
    import ad.math.traits: isNaN;

    assert(sin(GDN!1.zero) is GDN!1.zero);

    const g = sin(GDN!1(PI_2));
    assert(g == 1 && isClose(g.d, 0., 0., real.epsilon));

    assert(isNaN(sin(GDN!1.infinity)));
    assert(isNaN(sin(-GDN!1.infinity)));

    assert(sin(GDN!2(0)) is GDN!2(0, 1, 0));
    // f = 0
    // <f',f"> = cos(<0,1>)<1,0> = <1,0><1,0> = <1, 0>
}


/**
 * This function computes the cosine of the argument.
 *
 * If $(MATH f(x) = cos(g(x))), then $(MATH f' = -sin(g)g').
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` to compute the cosine of
 *
 * Returns:
 *   the cosine of `g`
 */
pragma(inline, true) pure nothrow @nogc @safe GDN!Deg cos(ulong Deg)(in GDN!Deg g)
{
    static if (Deg == 1)
        alias sine = core.math.sin;
    else
        alias sine = sin;

    return GDN!Deg(core.math.cos(g.val), -sine(g.reduce())*g.d);
}

///
unittest
{
    const g = GDN!2(0);
    const f = cos(g);
    assert(f == 1 && f.d == 0 && f.d!2 == -1);
}

unittest
{
    import std.math: isClose, PI_2;
    import ad.math.traits: isNaN;

    const g = cos(GDN!1(PI_2));
    assert(isClose(g.val, 0., 0., real.epsilon) && g.d == -1);

    assert(isNaN(cos(GDN!1.infinity)));
    assert(isNaN(cos(-GDN!1.infinity)));

    assert(cos(GDN!2(0)) is GDN!2(1, -0., -1));
    // f = 1
    // <f',f"> = -sin(<0,1>)<1,0> = -<0,1><1,0> = <-0,-1>
}


/**
 * This function computes the tangent of its argument
 *
 * If $(MATH f(x) = tan(g(x))), then $(MATH f' = g'sec$(SUP 2)(g)).
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` to take the tangent of
 *
 * Returns:
 *   the tangent of `g`
 */
pure nothrow @nogc @safe GDN!Deg tan(ulong Deg)(in GDN!Deg g)
{
    static if (Deg == 1)
        alias cosine = core.math.cos;
    else
        alias cosine = cos;

    return GDN!Deg(std.math.trigonometry.tan(g.val), g.d/pow(cosine(g.reduce()), 2));
}

///
unittest
{
    import std.math: isClose, PI_4;

    const f = tan(GDN!1(-PI_4));
    assert(f == -1 && isClose(f.d, 2));
}

unittest
{
    import std.math: isClose, PI_4;
    import ad.math.traits: isNaN;

    const w = tan(GDN!1(PI_4));
    assert(w.val == 1 && isClose(w.d, 2.));

    assert(isNaN(tan(GDN!1(real.nan))));
    assert(isNaN(tan(GDN!1(real.infinity))));

    assert(tan(GDN!2(0)) is GDN!2(0, 1, 0));
    // f = 0
    // <f',f"> = <1,0>/cos(<0,1>)^2
    //    = <1,0>/<1,-0>^2
    //    = <1,0>/<1,2*1*-0>
    //    = <1,0>/<1,-0>
    //    = <1,0>
}


/**
 * This function computes the arcsine (inverse sine) of its argument `g`.
 *
 * If $(MATH f(x) = sin$(SUP -1)g(x)), then $(MATH f' = g'/√(1 - g$(SUP 2)), |g| ≤ 1).
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` to compute the arcsine of
 *
 * Returns:
 *   the arcsine of `g`
 */
pure nothrow @nogc @safe GDN!Deg asin(ulong Deg)(in GDN!Deg g)
{
    return GDN!Deg(std.math.trigonometry.asin(g.val), g.d/sqrt(1 - pow(g.reduce(), 2)));
}

///
unittest
{
    import std.math: PI_2;

    assert(asin(GDN!1(0)) is GDN!1(0, 1));
    assert(asin(GDN!1(1)) is GDN!1(PI_2, real.infinity));
}

unittest
{
    import std.format: format;
    import std.math: PI_2;
    import ad.math.traits: isNaN;

    assert(asin(GDN!1(-1)) is GDN!1(-PI_2, real.infinity));
    assert(isNaN(asin(GDN!1(2))));

    const g = GDN!2(0);
    const f = asin(g);
    assert(f is GDN!2(0, 1, 0), format("asin(%s) != %s", g, f));
    // f = 0
    // <f',f"> = <1,0>/sqrt(1 - <0,1>^^2)
    //    = <1,0>/sqrt(1 - <0,0>)
    //    = <1,0>/sqrt(<1,0>)
    //    = <1,0>/<1,0>
    //    = <1,0>
}


/**
 * This function computes the arccosine (inverse cosine) of its argument `g`.
 *
 * If $(MATH f(x) = cos$(SUP -1)g(x)), then $(MATH f' = -g'/√(1 - g$(SUP 2)), |g| ≤ 1).
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` to compute the arccosine of
 *
 * Returns:
 *   the arccosine of `g`
 */
pure nothrow @nogc @safe GDN!Deg acos(ulong Deg)(in GDN!Deg g)
{
    return GDN!Deg(std.math.trigonometry.acos(g.val), -g.d/sqrt(1 - pow(g.reduce(), 2)));
}

///
unittest
{
    import std.math: PI_2;

    assert(acos(GDN!1(0)) is GDN!1(PI_2, -1));
    assert(acos(GDN!1(1)) is GDN!1(0, -real.infinity));
}

unittest
{
    import std.math: PI, PI_2;
    import ad.math.traits: isNaN;

    assert(acos(GDN!1(-1)) is GDN!1(PI, -real.infinity));
    assert(isNaN(acos(GDN!1(2))));
    assert(acos(GDN!2(0)) is GDN!2(PI_2, -1, 0));
}


/**
 * This function computes the arctangent (inverse tangent) of its argument `g`.
 *
 * If $(MATH f(x) = tan$(SUP -1)g(x)), then $(MATH f' = g'/(1 + g$(SUP 2))).
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` to compute the arctangent of
 *
 * Returns:
 *   the arctangent of `g`
 */
pure nothrow @nogc @safe GDN!Deg atan(ulong Deg)(in GDN!Deg g)
{
    return GDN!Deg(std.math.trigonometry.atan(g.val), g.d/(1 + pow(g.reduce(), 2)));
}

///
unittest
{
    assert(atan(GDN!1(0)) is GDN!1(0, 1));
}

unittest
{
    import std.math: PI_2;

    assert(atan(GDN!1(real.infinity)) is GDN!1(PI_2, +0.));
    assert(atan(GDN!1(-real.infinity)) is GDN!1(-PI_2, +0.));
    assert(atan(GDN!2(0)) is GDN!2(0, 1, 0));
}


/**
 * This function computes the arctangent (inverse tangent) of $(MATH g/h).
 *
 * If $(MATH f(x) = tan$(SUP -1)(g(x)/h(x))), then $(MATH f' = (g'h - gh')/(h$(SUP 2) + g$(SUP 2))).
 *
 * Params:
 *   G = the type of g
 *   H = the type of h
 *   g = the numerator of the arctangent argument
 *   h = the denominator of the arctangent argument
 *
 * Returns:
 *   It returns the angle resulting from the arctan(g/h).
 */
pure nothrow @nogc @safe
CommonGDN!(G, H) atan2(G, H)(in G g, in H h) if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H))
{
    alias RGDN = typeof(return);

    const gg = asGDN!(RGDN.DEGREE)(g);
    const hh = asGDN!(RGDN.DEGREE)(h);
    const gg_red = gg.reduce();
    const hh_red = hh.reduce();

    RGDN.DerivType!1 df;
    if (isFinite(gg) && isInfinity(hh)) {
        if (isFinite(gg.d) && isFinite(hh.d)) {
            df = gg.d / hh_red;
        }
    } else if (isInfinity(gg) && isFinite(hh)) {
        if (isFinite(gg.d) && isFinite(hh.d)) {
            df = -hh.d / gg_red;
        }
    } else {
        df =  (gg.d*hh_red - gg_red*hh.d) / (pow(hh_red, 2) + pow(gg_red, 2));
    }

    return RGDN(std.math.trigonometry.atan2(gg.val, hh.val), df);
}

///
unittest
{
    import std.math: PI_4;

    assert(atan2(GDN!1(1), GDN!1(1)) is GDN!1(PI_4, 0));
    assert(atan2(GDN!1(1), -1) is GDN!1(3*PI_4, -0.5));
}

unittest
{
    import std.format: format;
    import std.math: isNaN, PI, PI_2, PI_4;

    assert(atan2(GDN!2(-1), GDN!2(1)) is GDN!2(-PI_4, 1, 0));
    // <f',f"> = (<1,0><1,1> - <-1,1><1,0>) / (<1,1>^2 + <-1,1>^2)
    //         = (<1,1>      - <-1,1>)      / (<1,2>   + <1,-2>)
    //         = <2,0>                      / <2,0>
    //         = <1,0>

    assert(atan2(GDN!1.nan, GDN!1(2)) is GDN!1.nan);
    assert(atan2(GDN!1(-1), GDN!1.nan) is GDN!1.nan);

    const nz = GDN!1(-0.);
    const pz = GDN!1(+0.);

    assert(atan2(pz, GDN!1(1)) is GDN!1(+0., 1));
    assert(atan2(nz, GDN!1(3)) is GDN!1(-0., 1./3));

    const q = atan2(pz, pz);
    assert(q.val is +0. && isNaN(q.d), format("atan2(+0, +0) != %s", q));

    const w = atan2(nz, pz);
    assert(w.val is -0. && isNaN(w.d), format("atan2(-0, +0) != %s", w));

    assert(atan2(pz, GDN!1(-1)) is GDN!1(PI, -1));
    assert(atan2(nz, GDN!1(-1)) is GDN!1(-PI, -1));

    const e = atan2(pz, nz);
    assert(e == PI && isNaN(e.d));

    const r = atan2(nz, nz);
    assert(r == -PI && isNaN(r.d));

    assert(atan2(GDN!1(1), pz) is GDN!1(PI_2, -1));
    assert(atan2(GDN!1(1), nz) is GDN!1(PI_2, -1));
    assert(atan2(GDN!1(-1), pz) is GDN!1(-PI_2, 1));
    assert(atan2(GDN!1(-1), nz) is GDN!1(-PI_2, 1));

    const ni = GDN!1(-real.infinity);
    const pi = GDN!1(real.infinity);

    assert(atan2(GDN!1(1), pi) is GDN!1(+0., +0.));
    assert(atan2(GDN!1(-1), pi) is GDN!1(-0., +0.));

    const t = atan2(ni, GDN!1(0));
    assert(t is GDN!1(-PI_2, +0.), format("atan2(-inf, 0) != %s", t));

    assert(atan2(pi, GDN!1(0)) is GDN!1(PI_2,-0.));

    const y = atan2(GDN!1(1), ni);
    assert(y is GDN!1(PI, -0.), format("atan(1, -inf) != %s", y));

    const u = atan2(ni, pi);
    assert(u == -PI_4 && isNaN(u.d));

    const i = atan2(pi, pi);
    assert(i == PI_4 && isNaN(i.d));

    const o = atan2(ni, ni);
    assert(o == -3*PI_4 && isNaN(o.d));

    const p = atan2(pi, ni);
    assert(p == 3*PI_4 && isNaN(p.d));
}


/**
This function calculates the hyperbolic sine of its argument `g`.

If $(MATH f(x) = sinh(g(x))), then $(MATH f' = g'cosh(g)).
*/
// TODO: implement
GDN!Deg sinh(ulong Deg)(in GDN!Deg g) pure nothrow @nogc @safe;

/**
This function calculates the hyperbolic cosine of its argument `g`.

If $(MATH f(x) = cosh(g(x))), then $(MATH f' = g'sinh(g)).
*/
// TODO: implement
GDN!Deg cosh(ulong Deg)(in GDN!Deg g) pure nothrow @nogc @safe;

/**
This function calculates the hyperbolic tangent of its argument `g`.

If $(MATH f(x) = tanh(g(x))), then $(MATH f' = g'/cosh$(SUP 2)(g)).
*/
// TODO: implement
GDN!Deg tanh(ulong Deg)(in GDN!Deg g) pure nothrow @nogc @safe;

/**
This function calculates the inverse hyperbolic sine of its argument `g`.

If $(MATH f(x) = sinh$(SUP -1)g(x)), then $(MATH f' = g'/√(1 + g$(SUP 2))).
*/
// TODO: implement
GDN!Deg asinh(ulong Deg)(in GDN!Deg g) pure nothrow @nogc @safe;

/**
This function calculates the inverse hyperbolic cosine of its argument `g`.

If $(MATH f(x) = cosh$(SUP -1)g(x)), then $(MATH f' = g'/√(g$(SUP 2) - 1)).
*/
// TODO: implement
GDN!Deg acosh(ulong Deg)(in GDN!Deg g) pure nothrow @nogc @safe;

/**
This function calculates the inverse hyperbolic tangent of its argument `g`.

If $(MATH f(x) = tanh$(SUP -1)g(x)), then $(MATH f' = g'/(1 - g$(SUP 2)), |g| < 1).
*/
// TODO: implement
GDN!Deg atanh(ulong Deg)(in GDN!Deg g) pure nothrow @nogc @safe;
