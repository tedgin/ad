/// It extends the `std.math.trigonometry` to support `GDN` objects.
module ad.math.trigonometry;

static import core.math;
static import std.math.trigonometry;

import std.math: sqrt;

import ad.core;
import ad.math.internal: sqrt;


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

    return GDN!Deg(std.math.trigonometry.tan(g.val), g.d/cosine(g.reduce())^^2);
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
    return GDN!Deg(std.math.trigonometry.asin(g.val), g.d/sqrt(1 - g.reduce()^^2));
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
    return GDN!Deg(std.math.trigonometry.acos(g.val), -g.d/sqrt(1 - g.reduce()^^2));
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
    return GDN!Deg(std.math.trigonometry.atan(g.val), g.d/(1 + g.reduce()^^2));
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


/+
/**
This function computes the arctangent (inverse tangent) of $(MATH g/h).

If $(MATH f(x) = tan$(SUP -1)(g(x)/h(x))), the $(MATH f' = (g'h - gh')/(h$(SUP 2) + g$(SUP 2))).
*/
GDN!(GDeg < HDeg ? GDeg : HDeg) atan2(ulong GDeg, ulong HDeg)(in GDN!GDeg g, in GDN!HDeg h)
    pure nothrow @nogc @safe
{
    return atan2_impl(cast(typeof(return)) g, cast(typeof(return)) h);
}

/// ditto
GDN!Deg atan2(ulong Deg)(in GDN!Deg g, in real c) pure nothrow @nogc @safe
{
    return atan2_impl(g, GDN!Deg.mkConst(c));
}

/// ditto
GDN!Deg atan2(ulong Deg)(in real c, in GDN!Deg h) pure nothrow @nogc @safe
{
    return atan2_impl(GDN!Deg.mkConst(c), h);
}

///
unittest
{
    import std.math : PI_4;

    const f = atan2(GDN!1(1), GDN!1(1));
    assert(f == PI_4 && f.d == 0);

    const e = atan2(GDN!2(-1), GDN!1(1));
    assert(e == -PI_4 && e.d == 1);

    const q = atan2(GDN!1(1), -1);
    assert(q == 3 * PI_4 && q.d == -0.5);

    const w = atan2(-1, GDN!1(-1));
    assert(w == -3 * PI_4 && w.d == 0.5);
}

unittest
{
    import std.math : PI;

    assert(atan2(GDN!1(+0.), GDN!2(-1)).same(GDN!1(PI, -1)));
}

package pragma(inline, true)
GDN!Deg atan2_impl(ulong Deg)(in GDN!Deg g, in GDN!Deg h) pure nothrow @nogc @safe
{
    const g_red = g.reduce();
    const h_red = h.reduce();
    const gdh = g.d * h_red;
    const hdg = h.d * g_red;

    return GDN!Deg(std.math.atan2(g.val, h.val), (g.d*h_red - g_red*h.d)/(h_red^^2 + g_red^^2));
}

unittest
{
    import std.format;
    import std.math : PI;

    const nz = GDN!1(-0.);
    const pz = GDN!1(+0.);

    assert(atan2_impl(GDN!1(-0.), GDN!1(-1)).same(GDN!1(-PI, -1)));
    assert(atan2_impl(GDN!1.nan, GDN!1(2)).same(GDN!1.nan));
    assert(atan2_impl(GDN!1(-1), GDN!1.nan).same(GDN!1.nan));
    assert(atan2_impl(GDN!1(+0.), GDN!1(1)).same(GDN!1(+0., 1)));

    const h = GDN!1(3);
    const f = atan2_impl(nz, h);
    assert(f.same(GDN!1(-0., 1./3)), format("atan2(%s, %s) != %s", nz, h, f));

    const w = atan2_impl(pz, pz);
    assert(w.same(GDN!1(+0., real.nan)), format("atan2(%s, %s) != %s", pz, pz, w));

    const e = atan2_impl(nz, pz);
    assert(e.same(GDN!1(-0., real.infinity)), format("atan2(%s, %s) != %s", nz, pz, e));
    // lim g->0-,h->0+ (h-g)/(hh + gg) => h-g > 0

    // TODO: test + /-0, -0
    // TODO: test >0, + /-0
    // TODO: test <0, + /-0
    // TODO: test >0, inf
    // TODO: test + /-inf, finite
    // TODO: test >0, -inf
    // TODO: test + /-inf, inf
    // TODO: test + /-inf, -inf
}
+/

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
