/// It extends the `std.math.trigonometry` to support `GDN` objects.
module ad.math.trigonometry;

public import std.math.trigonometry;

static import core.math;
import std.math: sqrt;

import ad.core;
import ad.math.algebraic: sqrt;

/**
This function computes the sine of its argument.

If $(MATH f(x) = sin(g(x))), then $(MATH f' = cos(g)g').
*/
pragma(inline, true) GDN!Deg sin(ulong Deg)(in GDN!Deg g) nothrow pure @nogc @safe
{
    static if (Deg == 1)
        const cos_g = std.math.cos(g.reduce());
    else
        const cos_g = cos(g.reduce());

    return GDN!Deg(core.math.sin(g.val), cos_g*g.d);
}

///
unittest
{
    const f = sin(GDN!2(0));
    assert(f == 0 && f.d == 1 && f.d!2 == 0);
}

unittest
{
    import std.math : isClose, PI_2;

    assert(sin(GDN!1.zero).same(GDN!1.zero));

    const g = sin(GDN!1(PI_2));
    assert(g.val == 1 && isClose(g.d, 0., 0., real.epsilon));

    assert(sin(GDN!1.infinity).same(GDN!1.nan));
    assert(sin(-GDN!1.infinity).same(GDN!1.nan));

    assert(sin(GDN!2(0)).same(GDN!2(0, 1, 0)));
    // f = 0
    // <f',f"> = cos(<0,1>)<1,0> = <1,0><1,0> = <1, 0>
}

/**
This function computes the cosine of the argument.

If $(MATH f(x) = cos(g(x))), then $(MATH f' = -sin(g)g').
*/
pragma(inline, true) GDN!Deg cos(ulong Deg)(in GDN!Deg g) nothrow pure @nogc @safe
{
    static if (Deg == 1)
        const sin_g = std.math.sin(g.reduce());
    else
        const sin_g = sin(g.reduce());

    return GDN!Deg(core.math.cos(g.val), -sin_g*g.d);
}

///
unittest
{
    const f = cos(GDN!2(0));
    assert(f == 1 && f.d == 0 && f.d!2 == -1);
}

unittest
{
    import std.math : isClose, PI_2;

    const g = cos(GDN!1(PI_2));
    assert(isClose(g.val, 0., 0., real.epsilon) && g.d == -1);

    assert(cos(GDN!1.infinity).same(GDN!1.nan));
    assert(cos(-GDN!1.infinity).same(GDN!1.nan));

    const q = GDN!2(0);
    const w = cos(q);
    assert(w.same(GDN!2(1, -0., -1)));
    // f = 1
    // <f',f"> = -sin(<0,1>)<1,0> = -<0,1><1,0> = <-0,-1>
}


/**
This function computes the tangent of its argument

If $(MATH f(x) = tan(g(x))), then $(MATH f' = g'sec$(SUP 2)(g)).
*/
GDN!Deg tan(ulong Deg)(in GDN!Deg g) nothrow pure @nogc @safe
{
    static if (Deg == 1)
        const cos_g = std.math.cos(g.reduce());
    else
        const cos_g = cos(g.reduce());

    return GDN!Deg(std.math.tan(g.val), g.d/cos_g^^2);
}

///
unittest
{
    import std.math : isClose, PI_4;

    const g = GDN!1(-PI_4);
    const f = tan(g);
    assert(f == -1 && isClose(f.d, 2));
}

unittest
{
    import std.math : isClose, PI_4;

    const q = GDN!1(PI_4);
    const w = tan(q);
    assert(w.val == 1 && isClose(w.d, 2.));

    assert(tan(GDN!1(real.nan)).same(GDN!1.nan));
    assert(tan(GDN!1(real.infinity)).same(GDN!1.nan));

    assert(tan(GDN!2(0)).same(GDN!2(0, 1, 0)));
    // f = 0
    // <f',f"> = <1,0>/cos(<0,1>)^2
    //    = <1,0>/<1,-0>^2
    //    = <1,0>/<1,2*1*-0>
    //    = <1,0>/<1,-0>
    //    = <1,0>
}


/**
This function computes the arcsine (inverse sine) of its argument `g`.

If $(MATH f(x) = sin$(SUP -1)g(x)), then $(MATH f' = g'/√(1 - g$(SUP 2)), |g| ≤ 1).
*/
GDN!Deg asin(ulong Deg)(in GDN!Deg g) nothrow pure @nogc @safe
{
    return GDN!Deg(std.math.asin(g.val), g.d/sqrt(1 - g.reduce()^^2));
}

///
unittest
{
    import std.math : PI_2;

    const f = asin(GDN!1(0));
    assert(f == 0 && f.d == 1);

    const q = asin(GDN!1(1));
    assert(q == PI_2 && q.d == real.infinity);
}

unittest
{
    import std.format;
    import std.math : PI_2;

    assert(asin(GDN!1(-1)).same(GDN!1(-PI_2, real.infinity)));
    assert(asin(GDN!1(2)).same(GDN!1.nan));

    const g = GDN!2(0);
    const f = asin(g);
    const fp = GDN!1(1,0) / sqrt(1-GDN!1(0,1)^^2);
    assert(fp.same(GDN!1(1,0)), format("f' = %s", GDN!1(0,1) * GDN!1(0,1)));
    assert(f.same(GDN!2(0, 1, 0)), format("asin(%s) != %s", g, f));
    // f = 0
    // <f',f"> = <1,0>/sqrt(1 - <0,1>^^2)
    //    = <1,0>/sqrt(1 - <0,0>)
    //    = <1,0>/sqrt(<1,0>)
    //    = <1,0>/<1,0>
    //    = <1,0>
}

/**
This function computes the arccosine (inverse cosine) of its argument `g`.

If $(MATH f(x) = cos$(SUP -1)g(x)), then $(MATH f' = -g'/√(1 - g$(SUP 2)), |g| ≤ 1).
*/
GDN!Deg acos(ulong Deg)(in GDN!Deg g) nothrow pure @nogc @safe
{
    return GDN!Deg(std.math.acos(g.val), -g.d/sqrt(1 - g.reduce()^^2));
}

///
unittest
{
    import std.math : PI_2;

    const f = acos(GDN!1(0));
    assert(f == PI_2 && f.d == -1);

    const q = acos(GDN!1(1));
    assert(q == 0 && q.d == -real.infinity);
}

unittest
{
    import std.math : PI, PI_2;

    assert(acos(GDN!1(-1)).same(GDN!1(PI, -real.infinity)));
    assert(acos(GDN!1(2)).same(GDN!1.nan));
    assert(acos(GDN!2(0)).same(GDN!2(PI_2, -1, 0)));
}

/**
This function computes the arctangent (inverse tangent) of its argument `g`.

If $(MATH f(x) = tan$(SUP -1)g(x)), then $(MATH f' = g'/(1 + g$(SUP 2))).
*/
GDN!Deg atan(ulong Deg)(in GDN!Deg g) nothrow pure @nogc @safe
{
    return GDN!Deg(std.math.atan(g.val), g.d/(1 + g.reduce()^^2));
}

///
unittest
{
    const f = atan(GDN!1(0));
    assert(f == 0 && f.d == 1);
}

unittest
{
    import std.math : PI_2;

    assert(atan(GDN!1(real.infinity)).same(GDN!1(PI_2, +0.)));
    assert(atan(GDN!1(-real.infinity)).same(GDN!1(-PI_2, +0.)));
    assert(atan(GDN!2(0)).same(GDN!2(0, 1, 0)));
}

/**
This function computes the arctangent (inverse tangent) of $(MATH g/h).

If $(MATH f(x) = tan$(SUP -1)(g(x)/h(x))), the $(MATH f' = (g'h - gh')/(h$(SUP 2) + g$(SUP 2))).
*/
GDN!(GDeg < HDeg ? GDeg : HDeg) atan2(ulong GDeg, ulong HDeg)(in GDN!GDeg g, in GDN!HDeg h)
    nothrow pure @nogc @safe
{
    return atan2_impl(cast(typeof(return)) g, cast(typeof(return)) h);
}

/// ditto
GDN!Deg atan2(ulong Deg)(in GDN!Deg g, in real c) nothrow pure @nogc @safe
{
    return atan2_impl(g, GDN!Deg.mkConst(c));
}

/// ditto
GDN!Deg atan2(ulong Deg)(in real c, in GDN!Deg h) nothrow pure @nogc @safe
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
GDN!Deg atan2_impl(ulong Deg)(in GDN!Deg g, in GDN!Deg h) nothrow pure @nogc @safe
{
    const g_red = g.reduce();
    const h_red = h.reduce();
    return GDN!Deg(std.math.atan2(g.val, h.val), (g.d*h_red - g_red*h.d)/(h_red^^2 + g_red^^2));
}

unittest
{
    import std.math : PI;

    assert(atan2_impl(GDN!1(-0.), GDN!1(-1)).same(GDN!1(-PI, -1)));

    // TODO: test NaN, anything
    // TODO: test anything, NaN
    // TODO: test +/-0, >0
    // TODO: test +/-0, +0
    // TODO: test +/-0, <0
    // TODO: test +/-0, -0
    // TODO: test >0, +/-0
    // TODO: test <0, +/-0
    // TODO: test >0, inf
    // TODO: test +/-inf, finite
    // TODO: test >0, -inf
    // TODO: test +/-inf, inf
    // TODO: test +/-inf, -inf
}

/**
This function calculates the hyperbolic sine of its argument `g`.

If $(MATH f(x) = sinh(g(x))), then $(MATH f' = g'cosh(g)).
*/
// TODO: implement
GDN!Deg sinh(ulong Deg)(in GDN!Deg g) nothrow pure @nogc @safe;

/**
This function calculates the hyperbolic cosine of its argument `g`.

If $(MATH f(x) = cosh(g(x))), then $(MATH f' = g'sinh(g)).
*/
// TODO: implement
GDN!Deg cosh(ulong Deg)(in GDN!Deg g) nothrow pure @nogc @safe;

/**
This function calculates the hyperbolic tangent of its argument `g`.

If $(MATH f(x) = tanh(g(x))), then $(MATH f' = g'/cosh$(SUP 2)(g)).
*/
// TODO: implement
GDN!Deg tanh(ulong Deg)(in GDN!Deg g) nothrow pure @nogc @safe;

/**
This function calculates the inverse hyperbolic sine of its argument `g`.

If $(MATH f(x) = sinh$(SUP -1)g(x)), then $(MATH f' = g'/√(1 + g$(SUP 2))).
*/
// TODO: implement
GDN!Deg asinh(ulong Deg)(in GDN!Deg g) nothrow pure @nogc @safe;

/**
This function calculates the inverse hyperbolic cosine of its argument `g`.

If $(MATH f(x) = cosh$(SUP -1)g(x)), then $(MATH f' = g'/√(g$(SUP 2) - 1)).
*/
// TODO: implement
GDN!Deg acosh(ulong Deg)(in GDN!Deg g) nothrow pure @nogc @safe;

/**
This function calculates the inverse hyperbolic tangent of its argument `g`.

If $(MATH f(x) = tanh$(SUP -1)g(x)), then $(MATH f' = g'/(1 - g$(SUP 2)), |g| < 1).
*/
// TODO: implement
GDN!Deg atanh(ulong Deg)(in GDN!Deg g) nothrow pure @nogc @safe;
