/// It extends the `std.math.algebraic` to support `GDN` objects.
module ad.math.algebraic;

static import core.math;
static import std.math.algebraic;

import std.algorithm: min;
import std.traits: isFloatingPoint, isImplicitlyConvertible, Select;

import ad.core;
import ad.math.exponential: log2;
import ad.math.rounding: ceil, floor;
import ad.math.traits: isInfinity, sgn, signbit;

/**
 * This function computes the absolute value of the argument.
 *
 * If $(MATH f(x) = |g(x)|), then $(MATH f' = sgn(g)g'), when $(MATH g ≠ 0)
 *
 * Params:
 *   Deg = the degree of the `GDN` object to compute the absolute value of
 *   g = the `GDN` object to compute the absolute value of
 *
 * Returns:
 *   the absolute value of the `GDN` object
 */
pragma(inline, true) pure nothrow @nogc @safe GDN!Deg fabs(ulong Deg)(in GDN!Deg g)
{
    const df_val = signbit(g) == 0 ? 1.0L : -1.0L;

    static if (Deg == 1)
        const df = df_val;
    else
        const df = GDN!Deg.DerivType!1.mkConst(df_val);

    return GDN!Deg(core.math.fabs(g.val), df * g.d);
}

///
unittest
{
    assert(fabs(GDN!2(-3)) is GDN!2(3, -1, 0));
}

unittest
{
    assert(fabs(GDN!2(-3)) is GDN!2(3, -1, 0));
    assert(fabs(GDN!1(+0.)) is GDN!1(+0., 1));
    assert(fabs(GDN!1(-0.)) is GDN!1(+0., -1));
    assert(fabs(GDN!1(-1)) is GDN!1(1, -1));
    assert(fabs(GDN!1.nan) is GDN!1.nan);
}


// abs support
unittest
{
    assert(
        std.math.algebraic.abs(GDN!1(-1)) is GDN!1(1, -1),
        "std.math.algebraic.abs(GDN) not working");
}


/**
 * This function computes the square root of its argument.
 *
 * If $(MATH f(x) = √g(x)), then $(MATH f' = g$(SUP -½)g'/2).
 *
 * Params:
 *   Deg = the degree of the `GDN` object to compute the square root of
 *   g = the `GDN` object to compute the square root of
 *
 * Returns:
 *   the square root of the `GDN` object
 */
pragma(inline, true) pure nothrow @nogc @safe GDN!Deg sqrt(ulong Deg)(in GDN!Deg g)
{
    const dfdg = signbit(g) == 1 ? GDN!Deg.DerivType!1.nan : g.reduce()^^-0.5/2;
    return GDN!Deg(core.math.sqrt(g.val), dfdg * g.d);
}

///
unittest
{
    assert(sqrt(GDN!2(1)) is GDN!2(1, 0.5, -0.25));
    // f = 1
    // <f',f"> = <1,0>/[2*sqrt(<1,1>)] = .5<1,0><1,1>^-.5 = <.5,-.25>

    assert(sqrt(GDN!1(+0.)) is GDN!1(0, real.infinity));
}

unittest
{
    import std.format: format;
    import ad.math.traits: isNaN;

    assert(sqrt(GDN!1(-0.)) is GDN!1(-0., real.nan), "sqrt(-0) incorrect");

    const x = sqrt(-GDN!1.one);
    assert(isNaN(x), format("sqrt(-1) = %s, should be %s", x, GDN!1.nan));

    assert(sqrt(GDN!1.infinity) is GDN!1(real.infinity, 0), "sqrt(inf) incorrect");
}


/**
 * This function computes the cube root of the argument.
 *
 * If $(MATH f(x) = ∛g(x)), then $(MATH f' = g$(SUP -⅔)g'/3)
 *
 * Params:
 *   Deg = the degree of the `GDN` object to compute the cube root of
 *   g = the `GDN` object to compute the cube root of
 *
 * Returns:
 *   the cube root of the `GDN` object
 */
nothrow @nogc @safe GDN!Deg cbrt(ulong Deg)(in GDN!Deg g)
{
    if (g == 0) {
        return GDN!Deg(g.val, GDN!Deg.DerivType!1.infinity);
    }

    const p = -2.0L / 3;
    const g_pow = isInfinity(g) && g < 0 ? -(-g.reduce())^^p : g.reduce()^^p;
    return GDN!Deg(std.math.algebraic.cbrt(g.val), g.d * g_pow / 3);
}

///
unittest
{
    import std.math: isClose;

    const f = cbrt(GDN!2(8));
    assert(f == 2 && f.d == 1/12.0L && isClose(f.d!2, -1/144.0L));
}

unittest
{
    import std.format: format;

    assert(cbrt(GDN!1(1)) is GDN!1(1, 1/3.0L), "cbrt(1) incorrect");
    assert(cbrt(GDN!1(+0.)) is GDN!1(+0., real.infinity), "cbrt(0) incorrect");

    const g = GDN!1(-0.);
    const f = cbrt(g);
    const q = GDN!1(-0., real.infinity);
    assert(f is q, format("cbrt(%s) = %s, expected %s", g, f, q));

    const x = cbrt(GDN!1.infinity);
    assert(x is GDN!1(real.infinity, 0), format("cbrt(%s) != %s", GDN!1.infinity, x));

    const y = -GDN!1.infinity;
    const cy = cbrt(y);
    assert(cy is GDN!1(-real.infinity, 0), format("cbrt(%s) != %s", y, cy));
}


/**
 * Calculates the distance of the point $(MATH (g, h)) from the origin $(MATH (0, 0)). It either `g`
 * or `h` has type `real`, it is converted to a constant generalized dual number with the same
 * degree as the other parameter.
 *
 * If $(MATH f(x) = [g(x)$(SUP 2) + h(x)$(SUP 2)]$(SUP ½)), then
 * $(MATH f' = (gg' + hh')(g$(SUP 2) + h$(SUP 2))$(SUP -½)).
 *
 * Params:
 *   R = a type that is implicitly convertible to `real`
 *   GDeg = the degree of the `GDN` object `g`
 *   HDeg = the degree of the `GDN` object `h`
 *   g = the `GDN` object representing the x-coordinate
 *   h = the `GDN` object representing the y-coordinate
 *
 * Returns:
 *   The distance to the origin. The resulting `GDN` will have a degree equal to the lesser of the
 *   degrees of `g` and `h`.
 */
pure nothrow @nogc @safe
GDN!(min(GDeg, HDeg)) hypot(ulong GDeg, ulong HDeg)(in GDN!GDeg g, in GDN!HDeg h)
{
    return hypot_impl(g, h);
}

/// ditto
pure nothrow @nogc @safe GDN!GDeg hypot(R, ulong GDeg)(in GDN!GDeg g, in R h)
if (isImplicitlyConvertible!(R, real)) {
    return hypot_impl(g, GDN!GDeg.mkConst(h));
}

/// ditto
pure nothrow@nogc @safe GDN!HDeg hypot(R, ulong HDeg)(in R g, in GDN!HDeg h)
if (isImplicitlyConvertible!(R, real)) {
    return hypot_impl(GDN!HDeg.mkConst(g), h);
}

///
unittest
{
    import std.math: sqrt;

    assert(hypot(GDN!1(1), GDN!1(2)) is GDN!1(sqrt(5.0L), 3/sqrt(5.0L)));
}

private pragma(inline, true) pure nothrow @nogc @safe
GDN!Deg hypot_impl(ulong Deg)(in GDN!Deg g, in GDN!Deg h)
{
    const g_red = g.reduce();
    const h_red = h.reduce();

    static if (Deg == 1) {
        const f_red = std.math.algebraic.hypot(g_red, h_red);
        const f_val = f_red;
    } else {
        const f_red = hypot_impl(g_red, h_red);
        const f_val = f_red.val;
    }

    const df = (g_red * g.d + h_red * h.d) / f_red;
    return GDN!Deg(f_val, df);
}

unittest
{
    import std.math: isClose, sqrt;

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
 * Calculates the distance of the point $(MATH (g, h, i)) from the origin $(MATH (0, 0, 0)). If any
 * of `g`, `h` or `i` has type `real`, it is converted to a constant generalized dual number with
 * the same degree as the lesser degree of the other parameters.
 *
 * If $(MATH f(x) = [g(x)$(SUP 2) + h(x)$(SUP 2) + i(x)$(SUP 2)]$(SUP ½)), then
 * $(MATH f' = (gg' + hh' + ii')(g$(SUP 2) + h$(SUP 2) + i$(SUP 2))$(SUP -½)).
 *
 * Params:
 *   R = a type that is implicitly convertible to `real`
 *   S = a type that is implicitly convertible to `real`
 *   GDeg = the degree of the `GDN` object `g`
 *   HDeg = the degree of the `GDN` object `h`
 *   IDeg = the degree of the `GDN` object `i`
 *   g = the `GDN` object representing the x-coordinate
 *   h = the `GDN` object representing the y-coordinate
 *   i = the `GDN` object representing the z-coordinate
 *
 * Returns:
 *   The resulting generalized dual number will have a degree equal to the least of the degrees of
 *   `g`, `h`, and `i`.
 */
pure nothrow @nogc @safe
GDN!(min(GDeg, HDeg, IDeg))
hypot(ulong GDeg, ulong HDeg, ulong IDeg)(in GDN!GDeg g, in GDN!HDeg h, in GDN!IDeg i)
{
    alias Res = typeof(return);
    return hypot_impl(cast(Res) g, cast(Res) h, cast(Res) i);
}

/// ditto
pure nothrow @nogc @safe
GDN!(min(GDeg, HDeg)) hypot(R, ulong GDeg, ulong HDeg)(in GDN!GDeg g, in GDN!HDeg h, in R i)
if (isImplicitlyConvertible!(R, real)) {
    alias Res = typeof(return);
    return hypot_impl(cast(Res) g, cast(Res) h, Res.mkConst(i));
}

/// ditto
pure nothrow @nogc @safe
GDN!(min(GDeg, IDeg)) hypot(R, ulong GDeg, ulong IDeg)(in GDN!GDeg g, in R h, in GDN!IDeg i)
if (isImplicitlyConvertible!(R, real)) {
    alias Res = typeof(return);
    return hypot_impl(cast(Res) g, Res.mkConst(h), cast(Res) i);
}

/// ditto
pure nothrow @nogc @safe
GDN!(min(HDeg, IDeg)) hypot(R, ulong HDeg, ulong IDeg)(in R g, in GDN!HDeg h, in GDN!IDeg i)
if (isImplicitlyConvertible!(R, real)) {
    alias Res = typeof(return);
    return hypot_impl(Res.mkConst(g), cast(Res) h, cast(Res) i);
}

/// ditto
pure nothrow @nogc @safe GDN!GDeg hypot(R, S, ulong GDeg)(in GDN!GDeg g, in R h, in S i)
if (isImplicitlyConvertible!(R, real) && isImplicitlyConvertible!(S, real)) {
    return hypot_impl(g, GDN!GDeg.mkConst(h), GDN!GDeg.mkConst(i));
}

/// ditto
pure nothrow @nogc @safe GDN!HDeg hypot(R, S, ulong HDeg)(in R g, in GDN!HDeg h, in S i)
if (isImplicitlyConvertible!(R, real) && isImplicitlyConvertible!(S, real)) {
    return hypot_impl(GDN!HDeg.mkConst(g), h, GDN!HDeg.mkConst(i));
}

/// ditto
pure nothrow @nogc @safe GDN!IDeg hypot(R, S, ulong IDeg)(in R g, in S h, in GDN!IDeg i)
if (isImplicitlyConvertible!(R, real) && isImplicitlyConvertible!(S, real)) {
    return hypot_impl(GDN!IDeg.mkConst(g), GDN!IDeg.mkConst(h), i);
}

///
unittest
{
    import std.math: isClose, sqrt;

    const f = hypot(GDN!1(1), GDN!1(2), GDN!1(3));
    assert(isClose(f.val, sqrt(14.0L)) && isClose(f.d, 6/sqrt(14.0L)));
}

unittest
{
    import std.format: format;
    import std.math: sqrt;

    assert(typeof(hypot(GDN!2.one, GDN!2.one, GDN!1.one)).DEGREE == 1);
    assert(typeof(hypot(GDN!2.one, GDN!1.one, GDN!3.one)).DEGREE == 1);
    assert(typeof(hypot(GDN!2.one, GDN!4.one, GDN!3.one)).DEGREE == 2);
    assert(typeof(hypot(GDN!5(0), GDN!6(1), 2)).DEGREE == 5);
    assert(typeof(hypot(GDN!7(3), 4, GDN!1(5))).DEGREE == 1);
    assert(typeof(hypot(6, GDN!2(7), GDN!3(8))).DEGREE == 2);

    const e = hypot(GDN!1(0), 1, 2);
    assert(e is GDN!1(sqrt(5.0L), 0), format("hypot(0, 1, 2) != %s", e));

    const r = hypot(3, GDN!1(4), -1);
    const t = sqrt(26.0L);
    assert(r is GDN!1(t, 4/t));

    const y = hypot(-2, -3, GDN!1(-4));
    const u = sqrt(29.0L);
    assert(y is GDN!1(u, -4/u));
}

private pragma(inline, true) pure nothrow @nogc @safe
GDN!Deg hypot_impl(ulong Deg)(in GDN!Deg g, in GDN!Deg h, in GDN!Deg i)
{
    const g_red = g.reduce();
    const h_red = h.reduce();
    const i_red = i.reduce();

    static if (Deg == 1) {
        const f_red = std.math.algebraic.hypot(g_red, h_red, i_red);
        const f_val = f_red;
    } else {
        const f_red = hypot_impl(g_red, h_red, i_red);
        const f_val = f_red.val;
    }

    const df = (g_red*g.d + h_red*h.d + i_red*i.d) / f_red;
    return GDN!Deg(f_val, df);
}

unittest
{
    import std.math: isClose, sqrt;

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
 * This evaluates the polynomial $(MATH H(g) = h$(SUB 0) + h$(SUB 1)g + h$(SUB 2)g$(SUP 2) + ...)
 * using Horner's rule $(MATH H(g) = h$(SUB 0) + g⋅(h$(SUB 1) + g⋅(h$(SUB 2) + ...))). If either
 * $(MATH g) or $(MATH h$(SUB i)) have type `real`, they are converted to a constant GDN with the
 * same degree as the other parameters.
 *
 * If $(MATH f(x) = h$(SUB 0)(x) + h$(SUB 1)(x)g(x) + h$(SUB 2)(x)g(x)$(SUP 2) + ... + h$(SUB n)(x)g(x)$(SUP n)),
 * then
 * $(MATH f' = h$(SUB 0)' + h$(SUB 1)g' + g⋅(h$(SUB 1)' + 2h$(SUB 2)g' + g⋅(h$(SUB 2)' + 3h$(SUB 3)g' + g⋅(...(h$(SUB n)')...)))).
 *
 * Params:
 *   F = a floating point type
 *   GDeg = the degree of the argument of the polynomial
 *   HDeg = the degree of the coefficients of the polynomial
 *   g = the argument of the polynomial
 *   H = the coefficients of the polynomial
 *
 * Returns:
 *   the polynomial evaluated at `g`. The resulting GDN will have a degree equal to the lesser of
 *   the degrees of `g` and `H`.
 */
pure nothrow @nogc @safe
GDN!(min(GDeg, HDeg)) poly(ulong GDeg, ulong HDeg)(in GDN!GDeg g, in GDN!HDeg[] H)
in(H.length > 0, "coefficient array cannot be empty")
{
    return typeof(return)(poly_impl_base(g.val, H), poly_impl_deriv(g, H));
}

/// ditto
pure nothrow @nogc @safe GDN!GDeg poly(F, ulong GDeg)(in GDN!GDeg g, in F[] H)
if (isFloatingPoint!F)
in(H.length > 0, "coefficient array cannot be empty")
{
    static if (GDeg == 1)
        auto df_acc = 0.0L;
    else
        auto df_acc = GDN!GDeg.DerivType!1.zero;

    const g_red = g.reduce();

    ptrdiff_t n = H.length;
    n--;
    while (n > 0) {
        df_acc = n * H[n] * g.d + g_red * df_acc;
        n--;
    }

    return GDN!GDeg(std.math.algebraic.poly(g.val, H), df_acc);
}

/// ditto
pure nothrow @nogc @safe GDN!HDeg poly(F, ulong HDeg)(in F g, in GDN!HDeg[] H)
if (isFloatingPoint!F)
in(H.length > 0, "coefficient array cannot be empty")
{
    return GDN!HDeg(poly_impl_base(g, H), poly_impl_deriv(GDN!HDeg.mkConst(g), H));
}

/// ditto
pure nothrow @nogc @safe
GDN!(min(GDeg, HDeg)) poly(ulong GDeg, ulong HDeg, int N)(in GDN!GDeg g, ref const GDN!HDeg[N] H)
if (N > 0 && N <= 10) {
    return typeof(return)(poly_impl_base(g.val, H), poly_impl_deriv(g, H));
}

/// ditto
pure nothrow @nogc @safe GDN!GDeg poly(F, ulong GDeg, int N)(in GDN!GDeg g, const ref F[N] H)
if (isFloatingPoint!F && N > 0 && N <= 10) {
    static if (GDeg == 1)
        auto df_acc = 0.0L;
    else
        auto df_acc = GDN!GDeg.DerivType!1.zero;

    static foreach (i; 1 .. N) df_acc = (N-i) * H[N-i] * g.d + g.reduce() * df_acc;

    return GDN!GDeg(std.math.algebraic.poly(g.val, H), df_acc);
}

/// ditto
pure nothrow @nogc @safe GDN!HDeg poly(F, ulong HDeg, uint N)(in F g, const ref GDN!HDeg[N] H)
if (isFloatingPoint!F && N > 0 && N <= 10)
{
    return GDN!HDeg(poly_impl_base(g, H), poly_impl_deriv(GDN!HDeg.mkConst(g), H));
}

///
unittest
{
    assert(poly(GDN!1(3), [GDN!1(0), GDN!1(1), GDN!1(2)]) is GDN!1(21, 26));
    assert(poly(GDN!2(-1), [-2., -3., 4.]) is GDN!2(5, -11, 8));
}

unittest
{
    import std.format;

    assert(typeof(poly(GDN!1(0), [GDN!2(-1)])).DEGREE == 1);
    assert(typeof(poly(GDN!3(-2), [GDN!1(-3), GDN!1(4)])).DEGREE == 1);

    assert(poly(GDN!1(-4), [5., -5.]) is GDN!1(25, -5));
    // f = 25
    // f' = 0 + 0(-4) + -5*1 = -5

    const q = poly(0., [GDN!1(-1), GDN!1(2)]);
    assert(q is GDN!1(-1), format("%s", q));

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
    assert(u is GDN!1(1, -3), format("poly(<-1,1>, %s) = %s", y, u));

    assert(poly(GDN!2(-2), y) is GDN!2(6, -7, 4));
    // f = 6
    // <f',f"> = 0 + 1<1,0> + 2*2<-2,1><1,0>
    //    = <1,0> + 4<-2,1>
    //    = <1,0> + <-8,4>
    //    = <-7,4>
}

private pragma(inline, true) pure nothrow @nogc @safe
auto poly_impl_deriv(ulong GDeg, ulong HDeg)(in GDN!GDeg g, in GDN!HDeg[] h)
{
    alias Deg = Select!(GDeg < HDeg, GDeg, HDeg);

    const gd = cast(GDN!Deg) g;
    const gd_red = gd.reduce();

    ptrdiff_t n = h.length;
    n--;
    auto acc = (cast(GDN!Deg) h[n]).d;
    while (n > 0) {
        acc = (cast(GDN!Deg) h[n - 1]).d
            + n * (cast(GDN!Deg) h[n]).reduce() * gd.d
            + gd_red * acc;
        n--;
    }

    return acc;
}

unittest
{
    assert(poly_impl_deriv(GDN!2(0), [GDN!2(1), GDN!2(2)]) is GDN!1(3, 2));
    // f' = <h0',h0"> + <h1,h1'><g',g"> + <g,g'><h1',h1">
    //    = <1,0> + <2,1><1,0> + <0,1><1,0>
    //    = <1,0> + <2,1> + <0,1>
    //    = <3,2>
}

// Taken from std.math.algebraic.polyImplBase
private pure nothrow @nogc @safe real poly_impl_base(ulong Deg)(in real x, in GDN!Deg[] h)
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
    assert(poly_impl_base(0, [GDN!1(-1)]) == -1);
    assert(poly_impl_base(2, [GDN!3(-2), GDN!3(-3), GDN!3(4)]) == 8);
}


/**
 * Gives the next power of two after `g`.
 *
 * This function is equivalent to $(MATH lim$(SUB 𝜀⟶0$(SUP +)) sgn(g)2$(SUP ⌈lg|g| + 𝜀⌉)).
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the GDN to find the next power of 2 from
 *
 * Returns:
 *   The GDN object whose value is the next power of 2 after `g`.
 */
pure nothrow @nogc @safe GDN!Deg nextPow2(ulong Deg)(in GDN!Deg g)
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
    assert(nextPow2(GDN!1(3)) is GDN!1(4, 0));
    assert(nextPow2(GDN!1(1)) is GDN!1(2, real.infinity));
}

unittest
{
    import std.format: format;

    const e = GDN!1(-1);
    const r = nextPow2(e);
    assert(r is GDN!1(-2, real.infinity), format("nextPow2(%s) != %s", e, r));

    assert(nextPow2(GDN!1(+0.)) is GDN!1(+0., real.nan));
    assert(nextPow2(GDN!1(-0.)) is GDN!1(-0., real.nan));

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
    assert(w is GDN!2(0.125,0,0), format("nextPow(%s) != %s", q, w));
}

/**
 * Gives the previous power of two no larger than `g`.
 *
 * This function is equivalent to $(MATH sgn(g)2$(SUP ⌊lg|g|⌋)).
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the GDN to truncate to a poser of two
 *
 * Returns:
 *   `g` truncated to a power of 2
 */
pure nothrow @nogc @safe GDN!Deg truncPow2(ulong Deg)(in GDN!Deg g)
out (f; f.val == std.math.truncPow2(g.val), "result doesn't agree with std.math.truncPow2")
{
    return  sgn(g) * 2^^floor(log2(fabs(g)));
}

///
unittest
{
    assert(truncPow2(GDN!1(3)) is GDN!1(2, 0));
    assert(truncPow2(GDN!1(1)) is GDN!1(1, real.infinity));
}

unittest
{
    assert(truncPow2(GDN!1(-1)) is GDN!1(-1, real.infinity));
    assert(truncPow2(GDN!1(+0.)) is GDN!1(+0., real.nan));
    assert(truncPow2(GDN!1(-0.)) is GDN!1(-0., real.nan));
    assert(truncPow2(GDN!2(0.1)) is GDN!2(0.0625,0,0));
}
