/// This module extends the core.math and std.math libraries to supporting GenDualNum objects.
module ad.math;

public import core.math : toPrec, yl2x, yl2xp1;
public import std.math;
public import ad.core;

static import core.math;
import std.algorithm.iteration : map;
import std.array : array;
import std.traits : isFloatingPoint, isImplicitlyConvertible, Select;

static import ad.math.impl;

// core.math implementations
pragma(inline, true)
{
    /**
    This function computes the cosine of the argument.

    If $(MATH f(g) = cos(g)), then $(MATH f' = -sin(g)g').
    */
    GenDualNum!Degree cos(ulong Degree)(in GenDualNum!Degree x) nothrow pure @nogc @safe
    {
        return ad.math.impl.cos(x);
    }

    ///
    unittest
    {
        const f = cos(GenDualNum!2(0));
        assert(f == 1 && f.d == 0 && f.d!2 == -1);
    }

    /**
    This function computes $(MATH 2$(SUP exp)x).

    If $(MATH f(g;c) = 2$(SUP c)g), then $(MATH f' = 2$(SUP c)g').

    Params:
        x = the generalized dual number being scaled by $(MATH 2$(SUP exp)).
        exp = the power of $(MATH 2) used to scale `x`,
    */
    GenDualNum!Degree ldexp(ulong Degree)(in GenDualNum!Degree x, in int exp) nothrow pure
    @nogc @safe
    {
        return ad.math.impl.ldexp(x, exp);
    }

    ///
    unittest
    {
        const f = ldexp(GenDualNum!2(1), 2);
        assert(f == 4 && f.d == 4 && f.d!2 == 0);
    }

    /**
    Rounds `x` to the nearest integer value, using the current rounding mode. If the return value is
    not identical to `x`, the `FE_INEXACT` exception is raised.

    If $(MATH f(g) = rint(g)), then $(MATH f' = (df/dg)g'), where $(MATH df/dg = 𝛿(g - m)),
    $(MATH 𝛿) is the Dirac delta function, and $(MATH m) is a rounding mode split point.
    */
    GenDualNum!Degree rint(ulong Degree)(in GenDualNum!Degree x) nothrow pure @nogc @safe
    {
        const f = core.math.rint(x.val);

        auto df = GenDualNum!Degree.mkZeroDeriv();
        if (isInfinity(x.val))
        {
            df = GenDualNum!Degree.mkNaNDeriv();
        }
        else
        {
            auto fn = nearbyint(nextDown(x.val));
            auto fp = nearbyint(nextUp(x.val));

            if (f == 0)
            {
                if (signbit(f) == 1)
                    fp = f;
                else
                    fn = f;
            }

            if (fn != fp)
            {
                static if (Degree == 1)
                    df = real.infinity;
                else
                    df = (x.reduce() - x.val).dirac();
            }
        }

        return GenDualNum!Degree(f, df * x.d);
    }

    ///
    unittest
    {
        resetIeeeFlags();
        const e = rint(GenDualNum!2(1.5));
        assert(ieeeFlags.inexact);
        assert(e == 2 && e.d == real.infinity && isNaN(e.d!2));
    }

    unittest
    {
        import std.format;

        alias GDN = GenDualNum;

        const q = rint(GDN!1.infinity);
        assert(q.same(GDN!1(real.infinity, real.nan)), format("rint(inf) != %s", q));

        resetIeeeFlags();
        const w = rint(GDN!1.one);
        assert(!ieeeFlags.inexact);
        assert(w.same(GDN!1(1, 0)), format("rint(1) != %s", w));

        assert(rint(GDN!2(1.5)).same(GDN!2(2, real.infinity, real.nan)));

        const r = rint(GDN!1(-0.));
        assert(r.same(GDN!1(-0., 0)), format("rint(-0) != %s", r));

        assert(rint(GDN!1(+0.)).same(GDN!1(+0., 0)));
    }

    /**
    This function rounds `x` to a `long` using the current rounding mode. All of the derivative
    terms are lost.
    */
    long rndtol(ulong Degree)(in GenDualNum!Degree x) nothrow pure @nogc @safe
    {
        return core.math.rndtol(x.val);
    }

    ///
    unittest
    {
        assert(rndtol(GenDualNum!1(1.1)) == 1L);
    }

    /**
    This function computes the sine of its argument.

    If $(MATH f(g) = sin(g)), then $(MATH f' = cos(g)g').
    */
    GenDualNum!Degree sin(ulong Degree)(in GenDualNum!Degree x) nothrow pure @nogc @safe
    {
        return ad.math.impl.sin(x);
    }

    ///
    unittest
    {
        const f = sin(GenDualNum!2(0));
        assert(f == 0 && f.d == 1 && f.d!2 == 0);
    }

    /**
    This function rounds `x` to a given floating point type removing all derivative information.

    Params:
        T = the float point type to be converted to
        x = the generalized dual number to be converted
    */
    T toPrec(T, ulong Degree)(in GenDualNum!Degree x) nothrow pure @nogc @safe if (isFloatingPoint!T)
    {
        return core.math.toPrec!T(x.val);
    }

    ///
    unittest
    {
        assert(typeid(toPrec!float(GenDualNum!1.zero)) == typeid(float));
    }

    /**
    This function computes $(MATH y⋅lg(x)). It either `x` or `y` has type `real`, it is converted to
    a constant generalized dual number with the same degree as the other parameter.

    If $(MATH f(x) = h(x)lg(g(x))), then $(MATH f' = h'lg(g) + hg'/(ln(2)g))

    Params:
        x = the argument of logarithm
        y = the multiplier of the logarithm

    Returns:
        The resulting generalized dual number will have a degree equal to the lesser of the degrees
        of `x` and `y`.
    */
    GenDualNum!(XDegree < YDegree ? XDegree : YDegree)
    yl2x(ulong XDegree, ulong YDegree)(in GenDualNum!XDegree x, in GenDualNum!YDegree y)
    nothrow pure @nogc @safe
    {
        return ad.math.impl.yl2x(x, y);
    }

    /// ditto
    GenDualNum!Degree yl2x(T, ulong Degree)(in GenDualNum!Degree x, in T c) nothrow pure @nogc @safe
            if (isImplicitlyConvertible!(T, real))
    {
        return ad.math.impl.yl2x(x, c);
    }

    /// ditto
    GenDualNum!Degree yl2x(T, ulong Degree)(in T c, in GenDualNum!Degree y) nothrow pure @nogc @safe
            if (isImplicitlyConvertible!(T, real))
    {
        return ad.math.impl.yl2x(c, y);
    }

    ///
    unittest
    {
        const f = yl2x(GenDualNum!1(2), GenDualNum!1(3));
        assert(f == 3 && f.d == 1 + 1.5 / std.math.LN2);

        const x = yl2x(GenDualNum!1(+0., -1), GenDualNum!1(1));
        assert(x == -real.infinity && x.d == -real.infinity);

        assert(typeof(yl2x(GenDualNum!2(1), GenDualNum!1(2))).DEGREE == 1);

        const z = yl2x(GenDualNum!1(1), 2.);
        assert(z == 0 && z.d == 2 / std.math.LN2);
    }

    unittest
    {
        assert(yl2x(1, GenDualNum!1(2)).same(GenDualNum!1(0, 0)));
    }

    /**
    Computes $(MATH y⋅lg(x + 1)), for $(MATH -(1 - √½) ≤ x ≤ +(1 - √½)). When $(MATH x) is outside
    of this interval, the results are undefined. It either `x` or `y` has type `real`, it is
    converted to a constant generalized dual number with the same degree as the other parameter.

    If $(MATH f(x) = h(x)lg(g(x) + 1)), then $(MATH f' = h'lg(g + 1) + hg'/[ln(2)(g + 1)])

    Params:
        x = the argument of logarithm
        y = the multiplier of the logarithm

    Returns:
        The resulting generalized dual number will have a degree equal to the lesser of the degrees
        of `x` and `y`.
    */
    GenDualNum!(XDegree < YDegree ? XDegree : YDegree)
    yl2xp1(ulong XDegree, ulong YDegree)(in GenDualNum!XDegree x, in GenDualNum!YDegree y)
    nothrow pure @nogc @safe
    {
        return ad.math.impl.yl2xp1(x, y);
    }

    /// ditto
    GenDualNum!Degree yl2xp1(T, ulong Degree)(in GenDualNum!Degree x, in T c) nothrow pure
    @nogc @safe if (isImplicitlyConvertible!(T, real))
    {
        return ad.math.impl.yl2xp1(x, c);
    }

    /// ditto
    GenDualNum!Degree yl2xp1(T, ulong Degree)(in T c, in GenDualNum!Degree y) nothrow pure
    @nogc @safe if (isImplicitlyConvertible!(T, real))
    {
        return ad.math.impl.yl2xp1(c, y);
    }

    ///
    unittest
    {
        const f = yl2xp1(GenDualNum!1(0), GenDualNum!1(3));
        assert(f == 0 && f.d == 3 / std.math.LN2);

        assert(typeof(yl2xp1(GenDualNum!2(0), GenDualNum!1(1))).DEGREE == 1);

        const y = yl2xp1(GenDualNum!1(0), 1);
        assert(y == 0 && y.d == 1 / std.math.LN2);
    }

    unittest
    {
        assert(yl2xp1(0, GenDualNum!1(1)).same(GenDualNum!1(0, 0)));
    }
}

// CONSTANTS

/// Euler's constant $(MATH e)
enum GenDualNum!Degree E(ulong Degree) = GenDualNum!Degree.mkConst(std.math.E);

/// $(MATH π)
enum GenDualNum!Degree PI(ulong Degree) = GenDualNum!Degree.mkConst(std.math.PI);

/// $(MATH π/2)
enum GenDualNum!Degree PI_2(ulong Degree) = GenDualNum!Degree.mkConst(std.math.PI_2);

/// $(MATH π/4)
enum GenDualNum!Degree PI_4(ulong Degree) = GenDualNum!Degree.mkConst(std.math.PI_4);

/// $(MATH 1/π)
enum GenDualNum!Degree M_1_PI(ulong Degree) = GenDualNum!Degree.mkConst(std.math.M_1_PI);

/// $(MATH 2/π)
enum GenDualNum!Degree M_2_PI(ulong Degree) = GenDualNum!Degree.mkConst(std.math.M_2_PI);

/// $(MATH 2/√π)
enum GenDualNum!Degree M_2_SQRTPI(ulong Degree) = GenDualNum!Degree.mkConst(std.math.M_2_SQRTPI);

/// $(MATH ln(10))
enum GenDualNum!Degree LN10(ulong Degree) = GenDualNum!Degree.mkConst(std.math.LN10);

/// $(MATH ln(2))
enum GenDualNum!Degree LN2(ulong Degree) = GenDualNum!Degree.mkConst(std.math.LN2);

/// $(MATH log(2))
enum GenDualNum!Degree LOG2(ulong Degree) = GenDualNum!Degree.mkConst(std.math.LOG2);

/// $(MATH lg(e))
enum GenDualNum!Degree LOG2E(ulong Degree) = GenDualNum!Degree.mkConst(std.math.LOG2E);

/// $(MATH lg(10))
enum GenDualNum!Degree LOG2T(ulong Degree) = GenDualNum!Degree.mkConst(std.math.LOG2T);

/// $(MATH log(e))
enum GenDualNum!Degree LOG10E(ulong Degree) = GenDualNum!Degree.mkConst(std.math.LOG10E);

/// $(MATH √2)
enum GenDualNum!Degree SQRT2(ulong Degree) = GenDualNum!Degree.mkConst(std.math.SQRT2);

/// $(MATH √½)
enum GenDualNum!Degree SQRT1_2(ulong Degree) = GenDualNum!Degree.mkConst(std.math.SQRT1_2);

unittest
{
    assert(E!1 == std.math.E, "E broken");
    assert(PI!1 == std.math.PI, "PI broken");
    assert(PI_2!1 == std.math.PI_2, "PI_2 broken");
    assert(PI_4!1 == std.math.PI_4, "PI_4 broken");
    assert(M_1_PI!1 == std.math.M_1_PI, "M_1_PI broken");
    assert(M_2_PI!1 == std.math.M_2_PI, "M_2_PI broken");
    assert(M_2_SQRTPI!1 == std.math.M_2_SQRTPI, "SQRT1_2 broken");
    assert(LN10!2 == std.math.LN10, "LN10 broken");
    assert(LN2!1 == std.math.LN2, "LN2 broken");
    assert(LOG2!1 == std.math.LOG2, "LOG2 broken");
    assert(LOG2E!1 == std.math.LOG2E, "LOG2E broken");
    assert(LOG2T!1 == std.math.LOG2T, "LOG2T broken");
    assert(LOG10E!1 == std.math.LOG10E, "LOG10E broken");
    assert(SQRT2!1 == std.math.SQRT2, "SQRT2 broken");
    assert(SQRT1_2!1 == std.math.SQRT1_2, "SQRT1_2 broken");
}

// ALGEBRAIC

// abs support
unittest
{
    assert(abs(GenDualNum!1(-1)).same(GenDualNum!1(1, -1)), "std.math.abs(GDN) not working");
}

/**
This function computes the absolute value of the argument.

If $(MATH f(g) = |g|), then $(MATH f' = sgn(g)g'), when $(MATH g ≠ 0)
*/
pragma(inline, true)
GenDualNum!Degree fabs(ulong Degree)(in GenDualNum!Degree x) nothrow pure @nogc @safe
{
    const df_val = signbit(x.val) == 0 ? 1.0L : -1.0L;

    static if (Degree == 1)
        const df = df_val;
    else
        const df = GenDualNum!(Degree - 1)(df_val, GenDualNum!(Degree - 1).mkZeroDeriv());

    return GenDualNum!Degree(core.math.fabs(x.val), df * x.d);
}

///
unittest
{
    const f = fabs(GenDualNum!2(-3));
    assert(f == 3 && f.d == -1 && f.d!2 == 0);

    const x = fabs(GenDualNum!1(+0.));
    assert(x == +0. && x.d == 1);

    const y = fabs(GenDualNum!1(-0.));
    assert(y == +0. && y.d == -1);
}

/**
This function computes the square root of its argument.

If $(MATH f(g) = √g), then $(MATH f' = g$(SUP -½)g'/2).
*/
pragma(inline, true)
GenDualNum!Degree sqrt(ulong Degree)(in GenDualNum!Degree x) nothrow pure @nogc @safe
{
    const df = signbit(x.val) == 1 ? GenDualNum!Degree.DerivType!1.nan : x.reduce() ^^ -0.5 / 2;
    return GenDualNum!Degree(core.math.sqrt(x.val), df * x.d);
}

///
unittest
{
    const f = sqrt(GenDualNum!2(1));
    // f = 1
    // <f',f"> = <1,0>/[2*sqrt(<1,1>)] = .5<1,0><1,1>^-.5 = <.5,-.25>
    assert(f == 1 && f.d == 0.5 && f.d!2 == -0.25);

    const x = sqrt(GenDualNum!1(+0.));
    assert(x == 0 && x.d == real.infinity);
}

unittest
{
    import std.format;

    alias GDN1 = GenDualNum!1;

    assert(sqrt(GDN1(-0.)).same(GDN1(-0., real.nan)), "sqrt(-0) incorrect");

    const x = sqrt(-GDN1.one);
    assert(x.same(GDN1.nan), format("sqrt(-1) = %s, should be %s", x, GDN1.nan));

    assert(sqrt(GDN1.infinity).same(GDN1(real.infinity, 0)), "sqrt(inf) incorrect");
}

/**
This function computes the cube root of the argument.

If $(MATH f(g) = ∛g), then $(MATH f' = g$(SUP -⅔)g'/3)
*/
GenDualNum!Degree cbrt(ulong Degree)(in GenDualNum!Degree x) nothrow @nogc @safe
{
    alias GDN = GenDualNum!Degree;

    if (x == 0)
        return GDN(0, GDN.DerivType!1.infinity);

    const p = -2.0L / 3;
    const x_pow = isInfinity(x.val) && x.val < 0 ? -(-x.reduce()) ^^ p : x.reduce() ^^ p;
    return GDN(std.math.cbrt(x.val), x.d * x_pow / 3);
}

///
unittest
{
    const f = cbrt(GenDualNum!2(8));
    assert(f == 2 && f.d == 1.0L / 12 && isClose(f.d!2, -1.0L / 144));
}

unittest
{
    import std.format;

    alias GDN1 = GenDualNum!1;

    assert(cbrt(GDN1(1)).same(GDN1(1, 1.0L / 3.0L)), "cbrt(1) incorrect");
    assert(cbrt(GDN1(0)).same(GDN1(0, real.infinity)), "cbrt(0) incorrect");
    assert(cbrt(-GDN1(0)).same(GDN1(-0, real.infinity)), "cbrt(-0) incorrect");

    const x = cbrt(GDN1.infinity);
    assert(x.same(GDN1(real.infinity, 0)), format("cbrt(%s) != %s", GDN1.infinity, x));

    const y = -GDN1.infinity;
    const cy = cbrt(y);
    assert(cy.same(GDN1(-real.infinity, 0)), format("cbrt(%s) != %s", y, cy));

}

/**
Calculates the distance of the point $(MATH (x, y)) from the origin $(MATH (0, 0)). It either `x` or
`y` has type `real`, it is converted to a constant generalized dual number with the same degree as
the other parameter.

If $(MATH f(x) = [g(x)$(SUP 2) + h(x)$(SUP 2)]$(SUP ½)), then
$(MATH f' = (gg' + hh')(g$(SUP 2) + h$(SUP 2))$(SUP -½)).
*/
GenDualNum!(XDegree < YDegree ? XDegree : YDegree)
hypot(ulong XDegree, ulong YDegree)(in GenDualNum!XDegree x, in GenDualNum!YDegree y)
nothrow pure @nogc @safe
{
    return ad.math.impl.hypot(x, y);
}

/// ditto
pragma(inline, true)
GenDualNum!Degree hypot(T, ulong Degree)(in GenDualNum!Degree x, in T c) nothrow pure @nogc @safe
        if (isImplicitlyConvertible!(T, real))
{
    return ad.math.impl.hypot(x, c);
}

/// ditto
pragma(inline, true)
GenDualNum!Degree hypot(T, ulong Degree)(in T c, in GenDualNum!Degree y) nothrow pure @nogc @safe
{
    return ad.math.impl.hypot(c, y);
}

///
unittest
{
    const f = hypot(GenDualNum!1(1), GenDualNum!1(2));
    assert(f == std.math.sqrt(5.0L) && f.d == 3 / std.math.sqrt(5.0L));
}

/**
Calculates the distance of the point $(MATH (x, y, z)) from the origin $(MATH (0, 0, 0)). If any of
`x`, `y` or `z` has type `real`, it is converted to a constant generalized dual number with the same
degree as the least degree of the other parameters.

If $(MATH f(x) = [g(x)$(SUP 2) + h(x)$(SUP 2) + i(x)$(SUP 2)]$(SUP ½)), then
$(MATH f' = (gg' + hh' + ii')(g$(SUP 2) + h$(SUP 2) + i$(SUP 2))$(SUP -½)).

Returns:
    The resulting generalized dual number will have a degree equal to the least of the degrees of
    `x`, `y`, and `z`.
*/
auto
hypot(ulong XDegree, ulong YDegree, ulong ZDegree)(
    in GenDualNum!XDegree x, in GenDualNum!YDegree y, in GenDualNum!ZDegree z)
nothrow pure @nogc @safe
{
    alias XYDeg = Select!(XDegree < YDegree, XDegree, YDegree);
    alias Deg = Select!(XYDeg < ZDegree, XYDeg, ZDegree);

    return ad.math.impl.hypot(cast(GenDualNum!Deg) x, cast(GenDualNum!Deg) y, cast(GenDualNum!Deg) z);
}

/// ditto
GenDualNum!(XDegree < YDegree ? XDegree : YDegree)
hypot(T, ulong XDegree, ulong YDegree)(in GenDualNum!XDegree x, in GenDualNum!YDegree y, in T c)
        if (isImplicitlyConvertible!(T, real))
{
    return hypot(x, y, GenDualNum!(XDegree < YDegree ? XDegree : YDegree).mkConst(c));
}

/// ditto
GenDualNum!(XDegree < ZDegree ? XDegree : ZDegree)
hypot(T, ulong XDegree, ulong ZDegree)(in GenDualNum!XDegree x, in T c, in GenDualNum!ZDegree z)
        if (isImplicitlyConvertible!(T, real))
{
    return hypot(x, GenDualNum!(XDegree < ZDegree ? XDegree : ZDegree).mkConst(c), z);
}

/// ditto
GenDualNum!(YDegree < ZDegree ? YDegree : ZDegree)
hypot(T, ulong YDegree, ulong ZDegree)(in T c, in GenDualNum!YDegree y, in GenDualNum!ZDegree z)
        if (isImplicitlyConvertible!(T, real))
{
    return hypot(GenDualNum!(YDegree < ZDegree ? YDegree : ZDegree).mkConst(c), y, z);
}

/// ditto
GenDualNum!Degree hypot(T, U, ulong Degree)(in GenDualNum!Degree x, in T c1, in U c2)
        if (isImplicitlyConvertible!(T, real) && isImplicitlyConvertible!(U, real))
{
    return hypot(x, GenDualNum!Degree.mkConst(c1), GenDualNum!Degree.mkConst(c2));
}

/// ditto
GenDualNum!Degree hypot(T, U, ulong Degree)(in T c1, in GenDualNum!Degree y, in U c2)
        if (isImplicitlyConvertible!(T, real) && isImplicitlyConvertible!(U, real))
{
    return hypot(GenDualNum!Degree.mkConst(c1), y, GenDualNum!Degree.mkConst(c2));
}

/// ditto
GenDualNum!Degree hypot(T, U, ulong Degree)(in T c1, in U c2, in GenDualNum!Degree z)
        if (isImplicitlyConvertible!(T, real) && isImplicitlyConvertible!(U, real))
{
    return hypot(GenDualNum!Degree.mkConst(c1), GenDualNum!Degree.mkConst(c2), z);
}

///
unittest
{
    const f = hypot(GenDualNum!1(1), GenDualNum!1(2), GenDualNum!1(3));
    assert(isClose(f.val, std.math.sqrt(14.0L)) && isClose(f.d, 6 / std.math.sqrt(14.0L)));
}

unittest
{
    import std.format;

    alias GDN = GenDualNum;

    assert(typeof(hypot(GDN!2.one, GDN!2.one, GDN!1.one)).DEGREE == 1);
    assert(typeof(hypot(GDN!2.one, GDN!1.one, GDN!3.one)).DEGREE == 1);
    assert(typeof(hypot(GDN!2.one, GDN!4.one, GDN!3.one)).DEGREE == 2);
    assert(typeof(hypot(GDN!5(0), GDN!6(1), 2)).DEGREE == 5);
    assert(typeof(hypot(GDN!7(3), 4, GDN!1(5))).DEGREE == 1);
    assert(typeof(hypot(6, GDN!2(7), GDN!3(8))).DEGREE == 2);

    const e = hypot(GDN!1(0), 1, 2);
    assert(e.same(GDN!1(std.math.sqrt(5.0L), 0)), format("hypot(0, 1, 2) != %s", e));

    const r = hypot(3, GDN!1(4), -1);
    const t = std.math.sqrt(26.0L);
    assert(r.same(GDN!1(t, 4 / t)));

    const y = hypot(-2, -3, GDN!1(-4));
    const u = std.math.sqrt(29.0L);
    assert(y.same(GDN!1(u, -4 / u)));
}

/**
This evaluates the polynomial $(MATH A(x) = a$(SUB 0) + a$(SUB 1)x + a$(SUB 2)x$(SUP 2) + ...) using
Horner's rule $(MATH A(x) = a$(SUB 0) + x(a$(SUB 1) + x(a$(SUB 2) + ...))). It either $(MATH x) or
$(MATH a$(SUB i)) have type `real`, they are converted to a constant generalized dual number with
the same degree as the other parameters.

If $(MATH f(x) = h$(SUB 0)(x) + h$(SUB 1)(x)g(x) + h$(SUB 2)(x)g(x)$(SUP 2) + ... + h$(SUB n)(x)g(x)$(SUP n)),
then
$(MATH f' = h$(SUB 0)' + h$(SUB 1)g' + g⋅(h$(SUB 1)' + 2h$(SUB 2)g' + g⋅(h$(SUB 2)' + 3h$(SUB 3)g' + g⋅(...(h$(SUB n)')...)))).

Params:
    x = the argument of the polynomial
    A = the coefficients of the polynomial

Returns:
    The resulting generalized dual number will have a degree equal to the lesser of the degrees of
    `x` and `A`.
*/
GenDualNum!(XDegree < ADegree ? XDegree : ADegree)
poly(ulong XDegree, ulong ADegree)(in GenDualNum!XDegree x, in GenDualNum!ADegree[] A) nothrow pure
@nogc @trusted
in
{
    assert(A.length > 0);
}
do
{
    return typeof(return)(poly_impl_base(x.val, A), poly_impl_deriv(x, A));
}

/// ditto
GenDualNum!Degree poly(T, ulong Degree)(in GenDualNum!Degree x, in T[] A) nothrow pure
@nogc @trusted if (isFloatingPoint!T)
in
{
    assert(A.length > 0);
}
do
{
    static if (Degree == 1)
        auto df_acc = 0.0L;
    else
        auto df_acc = GenDualNum!Degree.DerivType!1.zero;

    const g = x.reduce();

    ptrdiff_t n = A.length;
    n--;
    while (n > 0)
    {
        df_acc = n * A[n] * x.d + g * df_acc;
        n--;
    }

    return GenDualNum!Degree(std.math.poly(x.val, A), df_acc);
}

/// ditto
GenDualNum!Degree poly(T, ulong Degree)(in T x, in GenDualNum!Degree[] A) nothrow pure
@nogc @trusted if (isFloatingPoint!T)
in
{
    assert(A.length > 0);
}
do
{
    return GenDualNum!Degree(
        poly_impl_base(x, A),
        poly_impl_deriv(GenDualNum!Degree.mkConst(x), A));
}

/// ditto
pragma(inline, true)
GenDualNum!(XDegree < ADegree ? XDegree : ADegree)
poly(ulong XDegree, ulong ADegree, int N)(
    in GenDualNum!XDegree x, ref const GenDualNum!ADegree[N] A)
nothrow pure @nogc @safe if (N > 0 && N <= 10)
{
    return typeof(return)(poly_impl_base(x.val, A), poly_impl_deriv(x, A));
}

/// ditto
pragma(inline, true)
GenDualNum!Degree poly(T, ulong Degree, int N)(in GenDualNum!Degree x, const ref T[N] A)
nothrow pure @nogc @safe if (isFloatingPoint!T && N > 0 && N <= 10)
{
    static if (Degree == 1)
        auto df_acc = 0.0L;
    else
        auto df_acc = GenDualNum!Degree.DerivType!1.zero;

    static foreach (i; 1 .. N)
    {
        df_acc = (N - i) * A[N - i] * x.d + x.reduce() * df_acc;
    }

    return GenDualNum!Degree(std.math.poly(x.val, A), df_acc);
}

/// ditto
pragma(inline, true)
GenDualNum!Degree poly(T, ulong Degree, uint N)(in T x, const ref GenDualNum!Degree[N] A)
nothrow pure @safe if (isFloatingPoint!T && N > 0 && N <= 10)
{
    return GenDualNum!Degree(
        poly_impl_base(x, A),
        poly_impl_deriv(GenDualNum!Degree.mkConst(x), A));
}

///
unittest
{
    import std.format;

    const q = poly(GenDualNum!1(3), [GenDualNum!1(0), GenDualNum!1(1), GenDualNum!1(2)]);
    assert(q == 21 && q.d == 26, format("%s", q));

    const w = poly(GenDualNum!2(-1), [-2., -3., 4.]);
    assert(w == 5 && w.d == -11 && w.d!2 == 8, format("%s", w));
}

unittest
{
    import std.format;

    alias GDN = GenDualNum;

    assert(typeof(poly(GDN!1(0), [GDN!2(-1)])).DEGREE == 1);
    assert(typeof(poly(GDN!3(-2), [GDN!1(-3), GDN!1(4)])).DEGREE == 1);

    const w = poly(GDN!1(-4), [5., -5.]);
    // f = 25
    // f' = 0 + 0(-4) + -5*1 = -5
    assert(w.same(GDN!1(25, -5)));

    const q = poly(0., [GDN!1(-1), GDN!1(2)]);
    assert(q.same(GDN!1(-1)), format("%s", q));

    static GDN!3[2] e = [GDN!3(2), GDN!3(3)];

    const r = poly(GDN!2(-2), e);
    assert(typeof(r).DEGREE == 2, format("%s", r));

    assert(typeof(poly(-2., e)).DEGREE == 3);

    static real[3] y = [0., 1., 2.];

    const u = poly(GDN!1(-1), y);
    // f = 1
    // f' = 0 + 1*1 + -1(0 + 2*2*1 + -1(0))
    //    = 1 + -1(4)
    //    = 1 - 4
    //    = -3
    assert(u.same(GDN!1(1, -3)), format("poly(<-1,1>, %s) = %s", y, u));

    const i = poly(GDN!2(-2), y);
    // f = 6
    // <f',f"> = 0 + 1<1,0> + 2*2<-2,1><1,0>
    //    = <1,0> + 4<-2,1>
    //    = <1,0> + <-8,4>
    //    = <-7,4>
    assert(i.same(GDN!2(6, -7, 4)));
}

private pragma(inline, true)
GenDualNum!(GDeg < HDeg ? GDeg : HDeg).DerivType!1 poly_impl_deriv(ulong GDeg, ulong HDeg)(
    in GenDualNum!GDeg g, in GenDualNum!HDeg[] h)
nothrow pure @nogc @safe
in
{
    assert(h.length > 0);
}
do
{
    alias Deg = Select!(GDeg < HDeg, GDeg, HDeg);

    const gd = cast(GenDualNum!Deg) g;
    const gd_red = gd.reduce();

    ptrdiff_t n = h.length;
    n--;
    auto acc = (cast(GenDualNum!Deg) h[n]).d;
    while (n > 0)
    {
        acc = (cast(GenDualNum!Deg) h[n - 1]).d
            + n * (cast(GenDualNum!Deg) h[n]).reduce() * gd.d
            + gd_red * acc;
        n--;
    }

    return acc;
}

unittest
{
    alias GDN = GenDualNum;

    const q = poly_impl_deriv(GDN!2(0), [GDN!2(1), GDN!2(2)]);
    // f' = <h0',h0"> + <h1,h1'><g',g"> + <g,g'><h1',h1">
    //    = <1,0> + <2,1><1,0> + <0,1><1,0>
    //    = <1,0> + <2,1> + <0,1>
    //    = <3,2>
    assert(q.same(GDN!1(3, 2)));
}

// Taken from std.math.algebraic.polyImplBase
private real poly_impl_base(ulong Deg)(in real x, in GenDualNum!Deg[] h) nothrow pure @nogc @safe
in
{
    assert(h.length > 0);
}
do
{
    ptrdiff_t n = h.length;
    --n;
    auto acc = h[n].val;
    while (--n >= 0)
    {
        acc *= x;
        acc += h[n].val;
    }
    return acc;
}

unittest
{
    assert(poly_impl_base(0, [GenDualNum!1(-1)]) == -1);
    assert(poly_impl_base(2, [GenDualNum!3(-2), GenDualNum!3(-3), GenDualNum!3(4)]) == 8);
}

/+ TODO: Implement the following
nextPow2
truncPow2
+/

/+ TODO: Implement std.math.trigonometry support
// TODO: learn why remainder is not pure
DerivSeqType tan(DerivSeqType)(const DerivSeqType u)
in
{
    assert(remainder(u.val() - PI_2, PI) != 0);
}
do
{
    static if (DerivSeqType.Degree == 0) {
        return DerivSeqType(std.math.tan(u.val()));
    } else {
        return sin(u) / cos(u);
    }
}

pure DerivSeqType acos(DerivSeqType)(const DerivSeqType u)
in
{
    static if (DerivSeqType.Degree == 0) {
        assert(-1 <= u.val() && u.val() <= 1);
    } else {
        assert(-1 < u.val() && u.val() < 1);
    }
}
do
{
    static if (DerivSeqType.Degree == 0) {
        return std.math.acos(u.val());
    } else {
        return DerivSeq(
            std.math.acos(u.val()),
            -u.d() / sqrt(1 - pow(u.reduce(), 2)));
    }
}
+/

// TODO: Implement std.math.rounding support

/+ TODO: Implement std.math.exponential support
/**
This function computes the raises a given number to a given power. It is
analogous to std.math.pow().

Params:
    x = the base
    y = the exponent
*/
nothrow pure @safe real pow(in real x, real y)
{
    return std.math.pow(x, y);
}

/// ditto
nothrow pure @safe GenDualNum!Degree pow(ulong Degree)(in GenDualNum!Degree x, in real y)
{
    return x ^^ y;
}

/// ditto
nothrow pure @safe GenDualNum!Degree pow(ulong Degree)(in real x, in GenDualNum!Degree y)
{
    return x ^^ y;
}

/// ditto
nothrow pure @safe GenDualNum!Degree pow(ulong Degree)(in GenDualNum!Degree x, in GenDualNum!Degree y)
{
    return x ^^ y;
}

/**
This function raises e to a given power. It is analogous to std.math.exp().

Params:
    x = the power e is raised to.
*/
nothrow pure @safe real exp(in real x)
{
    return std.math.exp(x);
}

/// ditto
nothrow pure @safe GenDualNum!Degree exp(ulong Degree)(in GenDualNum!Degree x)
{
    return GenDualNum!Degree(exp(x.val), x.d * exp(x.reduce()));
}

unittest
{
    import std.format;

    assert(exp(GenDualNum!1.nan).same(GenDualNum!1.nan));

    assert(exp(GenDualNum!1.zero).same(GenDualNum!1.one));

    // XXX - derivSeq isn't implemented due to the following bug
    // https: //issues.dlang.org/show_bug.cgi?id=22621 is fixed.
    // assert(exp(GenDualNum!1.infinity).same(derivSeq(real.infinity, real.infinity)));

    const e = exp(-GenDualNum!1.infinity);
    assert(e.same(GenDualNum!1(0, -0.)), format("exp(-inf) != %s", e));
}

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
nothrow pure @safe GenDualNum!Degree log(ulong Degree)(in GenDualNum!Degree x)
{
    return x.log();
}
+/

// TODO: Implement std.math.remainder support
// TODO: Implement std.math.operations support

/+ TODO: Implement std.math.traits support
/**
This function determines whether its argument is a NaN.  It is analogous to
std.math.isNaN().

Params:
    x = the argument
*/
nothrow pure @safe bool isNaN(in real x)
{
    return std.math.isNaN(x);
}

/// ditto
nothrow pure @safe bool isNaN(ulong Degree)(in GenDualNum!Degree x)
{
    return std.math.isNaN(x.val);
}

unittest
{
    assert(isNaN(GenDualNum!1.nan));
}

/**
This function computes the sign of the argument. It is analogous to
std.math.sgn().

The derivate of sgn(x) evaluated at 0 is undefined or NaN. Otherwise it is 0 *
x'.

Params:
    x = the argument
*/
nothrow pure @safe real sgn(in real x) @nogc
{
    return std.math.sgn(x);
}

/// ditto
nothrow pure @safe GenDualNum!Degree sgn(ulong Degree)(in GenDualNum!Degree x)
{
    alias GDN = GenDualNum!Degree;

    if (isNaN(x))
        return GDN.nan;
    else if (x == 0)
        return GDN(0, GDN.DerivType!1.nan);
    else
        return GDN(std.math.sgn(x.val), GDN.DerivType!1.zero);
}

/// ditto
nothrow pure @safe GenDualNum!Degree sgn(ulong Degree : 1)(in GenDualNum!Degree x)
{
    alias GDN = GenDualNum!Degree;

    if (isNaN(x))
        return GDN.nan;
    else if (x == 0)
        return GDN(0, real.nan);
    else
        return GDN(std.math.sgn(x.val), 0);
}

unittest
{
    import std.format;

    assert(
        sgn(GenDualNum!1()).same(GenDualNum!1.nan),
        format("sgn(%s) is not %s", sgn(GenDualNum!1()), GenDualNum!1.nan));

    assert(sgn(GenDualNum!1(0, real.nan)).same(GenDualNum!1(0, real.nan)));

    // XXX - derivSeq isn't implemented due to the following bug
    // https: //issues.dlang.org/show_bug.cgi?id=22621 is fixed.
    // assert(sgn(derivSeq(1.0L, 2.0L, real.nan)).same(derivSeq(1.0L, 0.0L, real.nan)));

    assert(sgn(GenDualNum!1(real.infinity, 0)).same(GenDualNum!1(1, 0)));
    assert(sgn(GenDualNum!1(-real.infinity, 0)).same(GenDualNum!1(-1, 0)));
}
+/

// TODO: Implement std.math.hardware support
