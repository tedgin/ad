/// It extends the `std.math.algebraic` to support `GenDualNum` objects.
module ad.math.algebraic;

public import std.math.algebraic;

static import core.math;

import std.math : isInfinity, signbit;
import std.traits : isFloatingPoint, isImplicitlyConvertible, Select;

import ad.core;
import ad.math.exponential : log2;
import ad.math.rounding : ceil, floor;
import ad.math.traits : sgn;


// abs support
unittest
{
    assert(abs(GenDualNum!1(-1)).same(GenDualNum!1(1, -1)), "std.math.abs(GDN) not working");
}

/**
This function computes the absolute value of the argument.

If $(MATH f(x) = |g(x)|), then $(MATH f' = sgn(g)g'), when $(MATH g ≠ 0)
*/
pragma(inline, true) GenDualNum!Deg fabs(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe
{
    const df_val = signbit(g.val) == 0 ? 1.0L : -1.0L;

    static if (Deg == 1)
        const df = df_val;
    else
        const df = GenDualNum!Deg.DerivType!1.mkConst(df_val);

    return GenDualNum!Deg(core.math.fabs(g.val), df * g.d);
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

unittest
{
    alias GDN = GenDualNum ;

    assert(fabs(GDN!1(-1)).same(GDN!1(1, -1)));
}

/**
This function computes the square root of its argument.

If $(MATH f(x) = √g(x)), then $(MATH f' = g$(SUP -½)g'/2).
*/
pragma(inline, true) GenDualNum!Deg sqrt(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe
{
    const df = signbit(g.val) == 1 ? GenDualNum!Deg.DerivType!1.nan : g.reduce()^^-0.5/2;
    return GenDualNum!Deg(core.math.sqrt(g.val), df * g.d);
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

If $(MATH f(x) = ∛g(x)), then $(MATH f' = g$(SUP -⅔)g'/3)
*/
GenDualNum!Deg cbrt(ulong Deg)(in GenDualNum!Deg g) nothrow @nogc @safe
{
    alias GDN = GenDualNum!Deg;

    if (g == 0)
        return GDN(0, GDN.DerivType!1.infinity);

    const p = -2.0L / 3;
    const g_pow = isInfinity(g.val) && g < 0 ? -(-g.reduce())^^p : g.reduce()^^p;
    return GDN(std.math.cbrt(g.val), g.d * g_pow / 3);
}

///
unittest
{
    import std.math : isClose;

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
Calculates the distance of the point $(MATH (g, h)) from the origin $(MATH (0, 0)). It either `g` or
`h` has type `real`, it is converted to a constant generalized dual number with the same degree as
the other parameter.

If $(MATH f(x) = [g(x)$(SUP 2) + h(x)$(SUP 2)]$(SUP ½)), then
$(MATH f' = (gg' + hh')(g$(SUP 2) + h$(SUP 2))$(SUP -½)).
*/
GenDualNum!(GDeg < HDeg ? GDeg : HDeg)
hypot(ulong GDeg, ulong HDeg)(in GenDualNum!GDeg g, in GenDualNum!HDeg h)
nothrow pure @nogc @safe
{
    return hypot_impl(g, h);
}

/// ditto
GenDualNum!Deg hypot(T, ulong Deg)(in GenDualNum!Deg g, in T c) nothrow pure @nogc @safe
if (isImplicitlyConvertible!(T, real))
{
    return hypot_impl(g, GenDualNum!Deg.mkConst(c));
}

/// ditto
pragma(inline, true)
GenDualNum!Deg hypot(T, ulong Deg)(in T c, in GenDualNum!Deg h) nothrow pure @nogc @safe
{
    return hypot(h, c);
}

///
unittest
{
    const f = hypot(GenDualNum!1(1), GenDualNum!1(2));
    assert(f == std.math.sqrt(5.0L) && f.d == 3 / std.math.sqrt(5.0L));
}

private pragma(inline, true)
GenDualNum!Deg hypot_impl(ulong Deg)(in GenDualNum!Deg g, in GenDualNum!Deg h) nothrow pure
@nogc @safe
{
    const g_red = g.reduce();
    const h_red = h.reduce();

    static if (Deg == 1) {
        const f_red = std.math.hypot(g_red, h_red);
        const f_val = f_red;
    } else {
        const f_red = hypot_impl(g_red, h_red);
        const f_val = f_red.val;
    }

    const df = (g_red * g.d + h_red * h.d) / f_red;
    return GenDualNum!Deg(f_val, df);
}

unittest
{
    import std.math : isClose, sqrt;

    alias GDN = GenDualNum;

    const f = hypot_impl(GDN!2(1), GDN!2(2));
    // f = sqrt(5)
    // <f',f"> = (<1,1><1,0> + <2,1><1,0>) / <sqrt(5),3/sqrt(5)>
    //    = (<1,1> + <2,1>) / <sqrt(5),3/sqrt(5)>
    //    = <3,2>/<sqrt(5),3/sqrt(5)>
    //    = <3/sqrt(5) , (2sqrt(5) - 9/sqrt(5))/5>
    //    = <3/sqrt(5) , .4sqrt(5) - 1.8/sqrt(5)>
    const w = sqrt(5.0L);
    const q = GDN!2(w, 3/w, .4*w - 1.8/w);
    assert(isClose(f.d!2, q.d!2));
}

/**
Calculates the distance of the point $(MATH (g, h, i)) from the origin $(MATH (0, 0, 0)). If any of
`g`, `h` or `i` has type `real`, it is converted to a constant generalized dual number with the same
degree as the least degree of the other parameters.

If $(MATH f(x) = [g(x)$(SUP 2) + h(x)$(SUP 2) + i(x)$(SUP 2)]$(SUP ½)), then
$(MATH f' = (gg' + hh' + ii')(g$(SUP 2) + h$(SUP 2) + i$(SUP 2))$(SUP -½)).

Returns:
    The resulting generalized dual number will have a degree equal to the least of the degrees of
    `g`, `h`, and `i`.
*/
auto
hypot(ulong GDeg, ulong HDeg, ulong IDeg)(
    in GenDualNum!GDeg g, in GenDualNum!HDeg h, in GenDualNum!IDeg i)
nothrow pure @nogc @safe
{
    alias GHDeg = Select!(GDeg < HDeg, GDeg, HDeg);
    alias Deg = Select!(GHDeg < IDeg, GHDeg, IDeg);

    return hypot_impl(cast(GenDualNum!Deg) g, cast(GenDualNum!Deg) h, cast(GenDualNum!Deg) i);
}

/// ditto
GenDualNum!(GDeg < HDeg ? GDeg : HDeg)
hypot(T, ulong GDeg, ulong HDeg)(in GenDualNum!GDeg g, in GenDualNum!HDeg h, in T c)
nothrow pure @nogc @safe
if (isImplicitlyConvertible!(T, real))
{
    return hypot_impl(cast(typeof(return)) g, cast(typeof(return)) h, typeof(return).mkConst(c));
}

/// ditto
pragma(inline, true)
GenDualNum!(GDeg < IDeg ? GDeg : IDeg)
hypot(T, ulong GDeg, ulong IDeg)(in GenDualNum!GDeg g, in T c, in GenDualNum!IDeg i)
nothrow pure @nogc @safe
if (isImplicitlyConvertible!(T, real))
{
    return hypot(g, i, c);
}

/// ditto
pragma(inline, true)
GenDualNum!(HDeg < IDeg ? HDeg : IDeg)
hypot(T, ulong HDeg, ulong IDeg)(in T c, in GenDualNum!HDeg h, in GenDualNum!IDeg i)
nothrow pure @nogc @safe
if (isImplicitlyConvertible!(T, real))
{
    return hypot(h, i, c);
}

/// ditto
GenDualNum!Deg hypot(T, U, ulong Deg)(in GenDualNum!Deg g, in T c1, in U c2) nothrow pure
@nogc @safe
if (isImplicitlyConvertible!(T, real) && isImplicitlyConvertible!(U, real))
{
    return hypot_impl(g, GenDualNum!Deg.mkConst(c1), GenDualNum!Deg.mkConst(c2));
}

/// ditto
pragma(inline, true)
GenDualNum!Deg hypot(T, U, ulong Deg)(in T c1, in GenDualNum!Deg h, in U c2) nothrow pure
@nogc @safe
if (isImplicitlyConvertible!(T, real) && isImplicitlyConvertible!(U, real))
{
    return hypot(h, c1, c2);
}

/// ditto
pragma(inline, true)
GenDualNum!Deg hypot(T, U, ulong Deg)(in T c1, in U c2, in GenDualNum!Deg i) nothrow pure
@nogc @safe
if (isImplicitlyConvertible!(T, real) && isImplicitlyConvertible!(U, real))
{
    return hypot(i, c1, c2);
}

///
unittest
{
    import std.math : isClose;

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

private pragma(inline, true)
GenDualNum!Deg hypot_impl(ulong Deg)(in GenDualNum!Deg g, in GenDualNum!Deg h, in GenDualNum!Deg i)
nothrow pure @nogc @safe
{
    const g_red = g.reduce();
    const h_red = h.reduce();
    const i_red = i.reduce();

    static if (Deg == 1) {
        const f_red = std.math.hypot(g_red, h_red, i_red);
        const f_val = f_red;
    } else {
        const f_red = hypot_impl(g_red, h_red, i_red);
        const f_val = f_red.val;
    }

    const df = (g_red*g.d + h_red*h.d + i_red*i.d) / f_red;
    return GenDualNum!Deg(f_val, df);
}

unittest
{
    import std.math : isClose, sqrt;

    alias GDN = GenDualNum;

    const f = hypot_impl(GDN!2(1), GDN!2(2), GDN!2(3));
    // f = sqrt(14)
    // <f',f"> = (<1,1><1,0> + <2,1><1,0> + <3,1><1,0>) / <sqrt(14),6/sqrt(14)>
    //    = (<1,1> + <2,1> + <3,1>) / <sqrt(14),6/sqrt(14)>
    //    = <6,3>/<sqrt(14),6/sqrt(14)>
    //    = <6/sqrt(14) , (3sqrt(14) - 36/sqrt(14))/14>
    //    = <6/sqrt(14) , 3/sqrt(14) - 18/(7sqrt(14))>
    const w = sqrt(14.0L);
    const q = GDN!2(w, 6/w, 3/w - 18/(7*w));
    assert(isClose(f.d!2, q.d!2));
}

/**
This evaluates the polynomial $(MATH H(g) = h$(SUB 0) + h$(SUB 1)g + h$(SUB 2)g$(SUP 2) + ...) using
Horner's rule $(MATH H(g) = h$(SUB 0) + g⋅(h$(SUB 1) + g⋅(h$(SUB 2) + ...))). It either $(MATH g) or
$(MATH h$(SUB i)) have type `real`, they are converted to a constant generalized dual number with
the same degree as the other parameters.

If $(MATH f(x) = h$(SUB 0)(x) + h$(SUB 1)(x)g(x) + h$(SUB 2)(x)g(x)$(SUP 2) + ... + h$(SUB n)(x)g(x)$(SUP n)),
then
$(MATH f' = h$(SUB 0)' + h$(SUB 1)g' + g⋅(h$(SUB 1)' + 2h$(SUB 2)g' + g⋅(h$(SUB 2)' + 3h$(SUB 3)g' + g⋅(...(h$(SUB n)')...)))).

Params:
    g = the argument of the polynomial
    H = the coefficients of the polynomial

Returns:
    The resulting generalized dual number will have a degree equal to the lesser of the degrees of
    `g` and `H`.
*/
GenDualNum!(GDeg < HDeg ? GDeg : HDeg)
poly(ulong GDeg, ulong HDeg)(in GenDualNum!GDeg g, in GenDualNum!HDeg[] H)
nothrow pure @nogc @trusted
in(H.length > 0, "coefficient array cannot be empty")
{
    return typeof(return)(poly_impl_base(g.val, H), poly_impl_deriv(g, H));
}

/// ditto
GenDualNum!Deg poly(T, ulong Deg)(in GenDualNum!Deg g, in T[] H) nothrow pure @nogc @trusted
if (isFloatingPoint!T)
in(H.length > 0, "coefficient array cannot be empty")
{
    static if (Deg == 1)
        auto df_acc = 0.0L;
    else
        auto df_acc = GenDualNum!Deg.DerivType!1.zero;

    const g_red = g.reduce();

    ptrdiff_t n = H.length;
    n--;
    while (n > 0)
    {
        df_acc = n * H[n] * g.d + g_red * df_acc;
        n--;
    }

    return GenDualNum!Deg(std.math.poly(g.val, H), df_acc);
}

/// ditto
GenDualNum!Deg poly(T, ulong Deg)(in T g, in GenDualNum!Deg[] H) nothrow pure @nogc @trusted
if (isFloatingPoint!T)
in(H.length > 0, "coefficient array cannot be empty")
{
    return GenDualNum!Deg(
        poly_impl_base(g, H),
        poly_impl_deriv(GenDualNum!Deg.mkConst(g), H));
}

/// ditto
pragma(inline, true)
GenDualNum!(GDeg < HDeg ? GDeg : HDeg)
poly(ulong GDeg, ulong HDeg, int N)(in GenDualNum!GDeg g, ref const GenDualNum!HDeg[N] H)
nothrow pure @nogc @safe
if (N > 0 && N <= 10)
{
    return typeof(return)(poly_impl_base(g.val, H), poly_impl_deriv(g, H));
}

/// ditto
pragma(inline, true)
GenDualNum!Deg poly(T, ulong Deg, int N)(in GenDualNum!Deg g, const ref T[N] H) nothrow pure
@nogc @safe
if (isFloatingPoint!T && N > 0 && N <= 10)
{
    static if (Deg == 1)
        auto df_acc = 0.0L;
    else
        auto df_acc = GenDualNum!Deg.DerivType!1.zero;

    static foreach (i; 1 .. N)
    {
        df_acc = (N-i) * H[N-i] * g.d + g.reduce() * df_acc;
    }

    return GenDualNum!Deg(std.math.poly(g.val, H), df_acc);
}

/// ditto
pragma(inline, true)
GenDualNum!Deg poly(T, ulong Deg, uint N)(in T g, const ref GenDualNum!Deg[N] H) nothrow pure @nogc @safe
if (isFloatingPoint!T && N > 0 && N <= 10)
{
    return GenDualNum!Deg(
        poly_impl_base(g, H),
        poly_impl_deriv(GenDualNum!Deg.mkConst(g), H));
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

/**
Gives the next power of two after `g`.

This function is equivalent to $(MATH lim$(SUB 𝜀⟶0$(SUP +)) sgn(g)2$(SUP ⌈lg|g| + 𝜀⌉)).
*/
GenDualNum!Deg nextPow2(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe
out (f; f.val == std.math.nextPow2(g.val), "result doesn't agree with std.math.nextPow2")
{
    const lg_abs_g = log2(fabs(g));

    auto power = ceil(lg_abs_g);
    if (power == lg_abs_g) {
        power = power + 1;
    }

    return  sgn(g) * 2^^power;
}

///
unittest
{
    const q = nextPow2(GenDualNum!1(3));
    assert(q == 4 && q.d == 0);

    const w = nextPow2(GenDualNum!1(1));
    assert(w == 2 &&  w.d == real.infinity);
}

unittest
{
    import std.format;

    alias GDN = GenDualNum;

    const e = GDN!1(-1);
    const r = nextPow2(e);
    assert(r.same(GDN!1(-2, real.infinity)), format("nextPow2(%s) != %s", e, r));

    assert(nextPow2(GDN!1(+0.)).same(GDN!1(+0., real.nan)));
    assert(nextPow2(GDN!1(-0.)).same(GDN!1(-0., real.nan)));

    const q = GDN!2(0.1);
    const w = nextPow2(q);
    // h = sgn(g)
    // h = 1
    // <h',h"> = 2𝛿(<0.1,1>)<1,0> = 2<0,0><1,0> = <0,0>
    // <h,h',h"> = <1,0,0>
    // i = |g|
    // i = .1
    // <i',i"> = sgn(<.1,1>)<1,0>
    //    <1,0>*<1,0>
    //    <1,0>
    // <i,i',i"> = <.1,1,0>
    // j = lg(i)
    // j = lg(.1)
    // <j',j"> = <1,0>/(<.1,1>ln(2))
    //    = <1,0>/<.1ln(2),ln(2)>
    //    = <10/ln(2),(0*.1ln(2) - 1*ln(2))100/ln(2)^2)>
    //    = <10/ln(2),-ln(2)100/ln(2)^2>
    //    = <10/ln(2),-100/ln(2)>
    // <j,j',j"> = <lg(.1),10/ln(2),-100/ln(2)>
    // k = ⌈j⌉
    // k = -3
    // <k',k"> = <10/ln(2),-100/ln(2)>∑[i∊ℤ]𝛿(<lg(.1),10/ln(2)> - i)
    //    = <10/ln(2),-100/ln(2)><0,0>
    //    = <0,0>
    // <k,k',k"> = <-3,0,0>
    // l = 2^k
    // l = .125
    // <l',l"> = <2,0>^<-3,0>(<0,0>*<-3,0>/<2,0> + <0,0>ln<2,0>)
    //    = <.125,0*-3/2 + 0*ln(2)>(<0,0>/<2,0> + <0,0><ln(2),0>)
    //    = <.125,0>(<0,(0*2 - 0)/4> + <0,0+0>)
    //    = <.125,0>(<0,0> + <0,0>)
    //    = <.125,0><0,0>
    //    = <0,0 + 0>
    //    = <0,0>
    // <l,l',l"> = <.125,0,0>
    // f = h * l
    // f = 0.125
    // <f',f"> = <0,0><.125,0> + <1,0><0,0>
    //    = <0,0+0> + <0,0+0>
    //    = <0,0> + <0,0>
    //    = <0,0>
    assert(w.same(GDN!2(0.125,0,0)), format("nextPow(%s) != %s", q, w));
}

/**
Gives the previous power of two no larger than `g`.

This function is equivalent to $(MATH sgn(g)2$(SUP ⌊lg|g|⌋)).
*/
GenDualNum!Deg truncPow2(ulong Deg)(in GenDualNum!Deg g) nothrow pure @nogc @safe
out (f; f.val == std.math.truncPow2(g.val), "result doesn't agree with std.math.truncPow2")
{
    return  sgn(g) * 2^^floor(log2(fabs(g)));
}

///
unittest
{
    const q = truncPow2(GenDualNum!1(3));
    assert(q == 2 && q.d == 0);

    const w = truncPow2(GenDualNum!1(1));
    assert(w == 1 &&  w.d == real.infinity);
}

unittest
{
    alias GDN = GenDualNum;

    assert(truncPow2(GDN!1(-1)).same(GDN!1(-1, real.infinity)));
    assert(truncPow2(GDN!1(+0.)).same(GDN!1(+0., real.nan)));
    assert(truncPow2(GDN!1(-0.)).same(GDN!1(-0., real.nan)));
    assert(truncPow2(GDN!2(0.1)).same(GDN!2(0.0625,0,0)));
}
