/// It extends the `std.math.algebraic` module to support `GDN` objects.
module ad.math.algebraic;

static import core.math;
static import std.math.algebraic;

import std.algorithm: min;
import std.math: isInfinity;
import std.traits: isFloatingPoint, isImplicitlyConvertible, Select;

static import ad.math.internal;

import ad.core;
import ad.math.internal:
    areAll, asGDN, asReal, ceil, CommonGDN, floor, isGDN, isGDNOrReal, isOne, log2, sgn, signbit;


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
    return ad.math.internal.sqrt(g);
}

///
unittest
{
    assert(sqrt(GDN!2(1)) is GDN!2(1, 0.5, -0.25));
    assert(sqrt(GDN!1(+0.)) is GDN!1(0, real.infinity));
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
    const g_pow = isInfinity(g.val) && g < 0 ? -(-g.reduce())^^p : g.reduce()^^p;
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
 *   G = the type of `g`, either a `GDN` or a `real`
 *   H = the type of `h`, either a `GDN` or a `real`
 *   g = the `GDN` object representing the x-coordinate
 *   h = the `GDN` object representing the y-coordinate
 *
 * Returns:
 *   The distance to the origin. The resulting `GDN` will have a degree equal to the lesser of the
 *   degrees of `g` and `h`.
 */
pure nothrow @nogc @safe
CommonGDN!(G, H) hypot(G, H)(in G g, in H h) if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H))
{
    alias Deg = typeof(return).DEGREE;

    const gg = asGDN!Deg(g);
    const hh = asGDN!Deg(h);

    const g_red = gg.reduce();
    const h_red = hh.reduce();

    static if (Deg == 1)
        const f_red = std.math.algebraic.hypot(g_red, h_red);
    else
        const f_red = hypot(g_red, h_red);

    const df = (g_red * gg.d + h_red * hh.d) / f_red;
    return GDN!Deg(asReal(f_red), df);
}

///
unittest
{
    import std.math: sqrt;

    assert(hypot(GDN!1(1), GDN!1(2)) is GDN!1(sqrt(5.0L), 3/sqrt(5.0L)));
}

unittest
{
    import std.math: isClose, sqrt;

    const f = hypot(GDN!2(1), GDN!2(2));
    // f = sqrt(5)
    // <f',f"> = (<1,1><1,0> + <2,1><1,0>) / <sqrt(5),3/sqrt(5)>
    //    = (<1,1> + <2,1>) / <sqrt(5),3/sqrt(5)>
    //    = <3,2>/<sqrt(5),3/sqrt(5)>
    //    = <3/sqrt(5) , (2sqrt(5) - 9/sqrt(5))/5>
    //    = <3/sqrt(5) , .4sqrt(5) - 1.8/sqrt(5)>
    const w = sqrt(5.0L);
    const q = GDN!2(w, 3/w, .4*w - 1.8/w);
    assert(isClose(f.d!2, q.d!2));

    assert(hypot(GDN!2(0), GDN!1(1)) is GDN!1(1));
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
 *   G = the type of `g`, either a `GDN` or a `real`
 *   H = the type of `h`, either a `GDN` or a `real`
 *   I = the type of `i`, either a `GDN` or a `real`
 *   g = the `GDN` object representing the x-coordinate
 *   h = the `GDN` object representing the y-coordinate
 *   i = the `GDN` object representing the z-coordinate
 *
 * Returns:
 *   The resulting generalized dual number will have a degree equal to the least of the degrees of
 *   `g`, `h`, and `i`.
 */
pure nothrow @nogc @safe
CommonGDN!(G, H, I)
hypot(G, H, I)(in G g, in H h, in I i) if (isOne!(isGDN, G, H, I) && areAll!(isGDNOrReal, G, H, I))
{
    alias Deg = typeof(return).DEGREE;

    const gg = asGDN!Deg(g);
    const hh = asGDN!Deg(h);
    const ii = asGDN!Deg(i);

    const g_red = gg.reduce();
    const h_red = hh.reduce();
    const i_red = ii.reduce();

    static if (Deg == 1)
        const f_red = std.math.algebraic.hypot(g_red, h_red, i_red);
    else
        const f_red = hypot(g_red, h_red, i_red);

    const df = (g_red * gg.d + h_red * hh.d + i_red * ii.d) / f_red;
    return GDN!Deg(asReal(f_red), df);
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
    import std.math: isClose, sqrt;

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

    const f = hypot(GDN!2(1), GDN!2(2), GDN!2(3));
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
 *   G = the type of `g`
 *   H = the type of `h`
 *   g = the argument of the polynomial
 *   h = the coefficients of the polynomial
 *
 * Returns:
 *   the polynomial evaluated at `g`. The resulting GDN will have a degree equal to the lesser of
 *   the degrees of `G` and `H`.
 */
pure nothrow @nogc @safe
CommonGDN!(G, H) poly(G, H)(in G g, in H[] h) if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H))
in(h.length > 0, "coefficient array cannot be empty")
{
    return poly_impl_base(g, h);
}

/// ditto
pure nothrow @nogc @safe
CommonGDN!(G, H) poly(G, H, int N)(in G g, ref const H[N] h)
if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H) && N > 0 && N <= 10)
{
    return poly_impl_base(g, h);
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

    static GDN!3[2] e = [GDN!3(2), GDN!3(3)];

    const r = poly(GDN!2(-2), e);
    assert(typeof(r).DEGREE == 2, format("%s", r));
}

// Taken from std.math.algebraic.polyImplBase
private pragma(inline, true) pure nothrow @nogc @safe
CommonGDN!(G, H) poly_impl_base(G, H)(in G g, in H[] h)
if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H))
{
    alias Deg = typeof(return).DEGREE;

    ptrdiff_t n = h.length;
    --n;
    auto acc = asGDN!Deg(h[n]);
    while (--n >= 0) {
        acc = acc * asGDN!Deg(g) + asGDN!Deg(h[n]);
    }

    return acc;
}

unittest
{
    import std.format: format;

    const f = poly_impl_base(0, [GDN!1(-1)]);
    assert(f is GDN!1(-1), format("f = %s", f));
    // f = <-1,1>

    const w = poly_impl_base(GDN!1(1), [2]);
    assert(w is GDN!1(2, 0), format("w = %s", w));

    const q = poly_impl_base(GDN!2(2), [GDN!2(-2), GDN!2(-3), GDN!2(4)]);
    assert(q is GDN!2(8, 20, 18), format("q = %s", q));
    // q = h0 + h1*g + h2*g^2
    // q' = h'0 + h'1*g + h1*g' + h'2*g^2 + 2*h2*g*g'
    // q = -2 + -3*2 + 4*2^2
    // q = 8
    // <q',q"> = <1,0> + <1,0><2,1> + <-3,1><1,0> + <1,0><2,1>^2 + 2<4,1><2,1><1,0>
    //         = <1,0> + <2,1>      + <-3,1>      + <1,0><4,4>   + <8,2><2,1>
    //         = <0,2>                            + <4,4>        + <16,12>
    //         = <20,18>
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
out (f; f == std.math.nextPow2(g.val), "result doesn't agree with std.math.nextPow2")
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
out (f; f == std.math.truncPow2(g.val), "result doesn't agree with std.math.truncPow2")
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
