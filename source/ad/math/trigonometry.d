/// It extends the `std.math.trigonometry` to support `GenDualNum` objects.
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
pragma(inline, true) GenDualNum!Deg sin(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe
{
    static if (Deg == 1)
        const cos_g = std.math.cos(g.reduce());
    else
        const cos_g = cos(g.reduce());

    return GenDualNum!Deg(core.math.sin(g.val), cos_g*g.d);
}

///
unittest
{
    const f = sin(GenDualNum!2(0));
    assert(f == 0 && f.d == 1 && f.d!2 == 0);
}

unittest
{
    import std.math : isClose, PI_2;

    alias GDN = GenDualNum;

    assert(sin(GDN!1.zero).same(GenDualNum!1.zero));

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
pragma(inline, true) GenDualNum!Deg cos(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe
{
    static if (Deg == 1)
        const sin_g = std.math.sin(g.reduce());
    else
        const sin_g = sin(g.reduce());

    return GenDualNum!Deg(core.math.cos(g.val), -sin_g*g.d);
}

///
unittest
{
    const f = cos(GenDualNum!2(0));
    assert(f == 1 && f.d == 0 && f.d!2 == -1);
}

unittest
{
    import std.math : isClose, PI_2;

    alias GDN = GenDualNum;

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
GenDualNum!Deg tan(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe
{
    static if (Deg == 1)
        const cos_g = std.math.cos(g.reduce());
    else
        const cos_g = cos(g.reduce());

    return GenDualNum!Deg(std.math.tan(g.val), g.d/cos_g^^2);
}

///
unittest
{
    import std.math : isClose, PI_4;

    const g = GenDualNum!1(-PI_4);
    const f = tan(g);
    assert(f == -1 && isClose(f.d, 2));
}

unittest
{
    import std.math : isClose, PI_4;

    alias GDN = GenDualNum;

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
GenDualNum!Deg asin(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe
{
    return GenDualNum!Deg(std.math.asin(g.val), g.d/sqrt(1 - g.reduce()^^2));
}

///
unittest
{
    import std.math : PI_2;

    const f = asin(GenDualNum!1(0));
    assert(f == 0 && f.d == 1);

    const q = asin(GenDualNum!1(1));
    assert(q == PI_2 && q.d == real.infinity);
}

unittest
{
    import std.format;
    import std.math : PI_2;

    alias GDN = GenDualNum;

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
GenDualNum!Deg acos(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe
{
    return GenDualNum!Deg(std.math.acos(g.val), -g.d/sqrt(1 - g.reduce()^^2));
}

///
unittest
{
    import std.math : PI_2;

    const f = acos(GenDualNum!1(0));
    assert(f == PI_2 && f.d == -1);

    const q = acos(GenDualNum!1(1));
    assert(q == 0 && q.d == -real.infinity);
}

unittest
{
    import std.math : PI, PI_2;

    alias GDN = GenDualNum;

    assert(acos(GDN!1(-1)).same(GDN!1(PI, -real.infinity)));
    assert(acos(GDN!1(2)).same(GDN!1.nan));
    assert(acos(GDN!2(0)).same(GDN!2(PI_2, -1, 0)));
}

/**
This function computes the arctangent (inverse tangent) of its argument `g`.

If $(MATH f(x) = tan$(SUP -1)g(x)), then $(MATH f' = g'/(1 + g$(SUP 2))).
*/
GenDualNum!Deg atan(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe
{
    return GenDualNum!Deg(std.math.atan(g.val), g.d/(1 + g.reduce()^^2));
}

///
unittest
{
    const f = atan(GenDualNum!1(0));
    assert(f == 0 && f.d == 1);
}

unittest
{
    import std.math : PI_2;

    alias GDN = GenDualNum;

    assert(atan(GDN!1(real.infinity)).same(GDN!1(PI_2, +0.)));
    assert(atan(GDN!1(-real.infinity)).same(GDN!1(-PI_2, +0.)));
    assert(atan(GDN!2(0)).same(GDN!2(0, 1, 0)));
}

/**
This function computes the arctangent (inverse tangent) of $(MATH g/h).

If $(MATH f(x) = tan$(SUP -1)(g(x)/h(x))), the $(MATH f' = (g'h - gh')/(h$(SUP 2) + g$(SUP 2))).
*/
GenDualNum!(GDeg < HDeg ? GDeg : HDeg)
atan2(ulong GDeg, ulong HDeg)(in GenDualNum!GDeg g, in GenDualNum!HDeg h)
    nothrow pure @nogc @safe
{
    return atan2_impl(cast(typeof(return)) g, cast(typeof(return)) h);
}

/// ditto
GenDualNum!Deg atan2(ulong Deg)(in GenDualNum!Deg g, in real c) nothrow pure @nogc @safe
{
    return atan2_impl(g, GenDualNum!Deg.mkConst(c));
}

/// ditto
GenDualNum!Deg atan2(ulong Deg)(in real c, in GenDualNum!Deg h) nothrow pure @nogc @safe
{
    return atan2_impl(GenDualNum!Deg.mkConst(c), h);
}

///
unittest
{
    import std.math : PI_4;

    const f = atan2(GenDualNum!1(1), GenDualNum!1(1));
    assert(f == PI_4 && f.d == 0);

    const e = atan2(GenDualNum!2(-1), GenDualNum!1(1));
    assert(e == -PI_4 && e.d == 1);

    const q = atan2(GenDualNum!1(1), -1);
    assert(q == 3 * PI_4 && q.d == -0.5);

    const w = atan2(-1, GenDualNum!1(-1));
    assert(w == -3 * PI_4 && w.d == 0.5);
}

unittest
{
    import std.math : PI;

    alias GDN = GenDualNum;

    assert(atan2(GDN!1(+0.), GDN!2(-1)).same(GDN!1(PI, -1)));
}

package pragma(inline, true)
GenDualNum!Deg atan2_impl(ulong Deg)(in GenDualNum!Deg g, in GenDualNum!Deg h) nothrow pure
    @nogc @safe
{
    const g_red = g.reduce();
    const h_red = h.reduce();
    return GenDualNum!Deg(
        std.math.atan2(g.val, h.val),
        (g.d*h_red - g_red*h.d) / (h_red^^2 + g_red^^2));
}

unittest
{
    import std.math : PI;

    alias GDN = GenDualNum;

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
GenDualNum!Deg sinh(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe;

/**
This function calculates the hyperbolic cosine of its argument `g`.

If $(MATH f(x) = cosh(g(x))), then $(MATH f' = g'sinh(g)).
*/
// TODO: implement
GenDualNum!Deg cosh(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe;

/**
This function calculates the hyperbolic tangent of its argument `g`.

If $(MATH f(x) = tanh(g(x))), then $(MATH f' = g'/cosh$(SUP 2)(g)).
*/
// TODO: implement
GenDualNum!Deg tanh(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe;

/**
This function calculates the inverse hyperbolic sine of its argument `g`.

If $(MATH f(x) = sinh$(SUP -1)g(x)), then $(MATH f' = g'/√(1 + g$(SUP 2))).
*/
// TODO: implement
GenDualNum!Deg asinh(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe;

/**
This function calculates the inverse hyperbolic cosine of its argument `g`.

If $(MATH f(x) = cosh$(SUP -1)g(x)), then $(MATH f' = g'/√(g$(SUP 2) - 1)).
*/
// TODO: implement
GenDualNum!Deg acosh(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe;

/**
This function calculates the inverse hyperbolic tangent of its argument `g`.

If $(MATH f(x) = tanh$(SUP -1)g(x)), then $(MATH f' = g'/(1 - g$(SUP 2)), |g| < 1).
*/
// TODO: implement
GenDualNum!Deg atanh(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe;
