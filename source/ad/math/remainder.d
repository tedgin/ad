/// It extends the `std.math.remainder` to support `GDN` objects.
module ad.math.remainder;

static import std.math.remainder;

import std.math: sgn;

import ad.core;
import ad.math.internal:
    areAll, asGDN, CommonGDN, isFinite, isGDN, isGDNOrReal, isInfinity, isNaN, isOne, round, trunc;


/**
 * Returns the remainder of `g` divided by `h`.
 *
 * If $(MATH f(x) = g(x) (mod h(x))), then $(MATH f' = g' - (g - f)h'/h - δ(f)(g'h - gh')/h), where
 * $(MATH δ) is the Dirac Delta function.
 *
 * Params:
 *   G = the type of g, either `GDN` or `real`.
 *   H = the type of h, either `GDN` or `real`.
 *   g = the dividend.
 *   h = the divisor.
 *
 * Returns:
 *   The remainder of `g` divided by `h`.
 */
nothrow @nogc @safe
CommonGDN!(G, H) fmod(G, H)(in G g, in H h) if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H))
{
    return g % h;
}

///
unittest
{
    const f = fmod(GDN!1(5), GDN!1(3));
    assert(f is GDN!1(2,0));
}


/**
 * Breaks g into an integer and a fraction, each with the same sign as g.
 *
 * `f = modf(g, i)` can be expressed mathematically as follows. When $MATH(g ≥ 0), $(MATH i = ⌊g⌋),
 * and when $(MATH g < 0), $(MATH i = ⌈g⌉). This is the mathematical definition of `trunc(g)`.
 * $(MATH f = g - i).
 *
 * Params:
 *   Deg = the degree of g
 *   g = the GDN object to break into an integer and a fraction.
 *   i = the integer part of g.
 *
 * Returns:
 *   The fractional part of g.
 */
nothrow @nogc @safe GDN!Deg modf(ulong Deg)(in GDN!Deg g, out GDN!Deg i)
{
    i = trunc(g);

    if (isInfinity(g)) {
        return GDN!Deg(sgn(g.val) * 0., GDN!Deg.mkNaNDeriv());
    }

    return g - i;
}

///
unittest
{
    import std.math: isClose;

    const g = GDN!1(3.14159);
    GDN!1 i;
    const f = modf(GDN!1(g), i);
    assert(i is GDN!1(3, 0));
    assert(isClose(f.val, 0.14159));
    assert(f.d == 1);
}

unittest
{
    import std.format: format;

    GDN!1 i;

    const q = modf(GDN!1(-real.infinity), i);
    assert(
        q is GDN!1(-0., real.nan) && i is GDN!1(-real.infinity, real.nan),
        format("q: %s, i: %s", q, i));

    const w = modf(GDN!1(real.infinity), i);
    assert(w is GDN!1(0., real.nan) && i is GDN!1(real.infinity, real.nan));
}


/** Calculates integer quotient of `g` and `h` and its remainder using the definition of remainder
 * provided by IEC 60559.
 *
 * The integer quotient n is $(MATH n(x) = round(g(x)/h(x))), and the remainder is defined as
 * $(MATH f(g(x), h(x)) = g(x) - h(x)n(x)).
 *
 * Params:
 *   G = the type of g, either `GDN` or `real`.
 *   H = the type of h, either `GDN` or `real`.
 *   g = the dividend.
 *   h = the divisor.
 *   n = the integer quotient of g and h.
 *
 * Returns:
 *   The remainder of `g` divided by `h`.
 */
nothrow @nogc @safe
CommonGDN!(G, H) remquo(G, H)(in G g, in H h, out int n)
if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H))
{
    alias Deg = typeof(return).DEGREE;

    const gg = asGDN!Deg(g);
    const hh = asGDN!Deg(h);

    if (isNaN(gg) || isNaN(hh)) return GDN!Deg.nanCombine(gg, hh);

    if (gg == 0 && hh != 0) {
        n = 0;
        return gg;
    }

    if (isFinite(gg) && isInfinity(hh)) return gg;

    const n_gdn = round(gg / hh);
    n = cast(int) n_gdn.val;
    return gg - hh * n_gdn;
}

///
unittest
{
    import std.math: isClose;

    int n;
    const f = remquo(GDN!1(5.1), GDN!1(3), n);
    assert(n ==2 && isClose(f.val, -0.9) && f.d == -1);
}

unittest
{
    import ad.math.traits: isNaN;

    int n;

    n = int.min;
    const q = remquo(GDN!1(-0.), GDN!1(1), n);
    assert(n == 0 && q is GDN!1(-0., 1));

    n = int.min;
    const w = remquo(GDN!1(+0.), GDN!1(1), n);
    assert(n == 0 && w is GDN!1(+0., 1));

    assert(isNaN(remquo(GDN!1(-real.infinity), GDN!1(1), n)));
    assert(isNaN(remquo(GDN!1(real.infinity), GDN!1(1), n)));
    assert(isNaN(remquo(GDN!1(1), GDN!1(0), n)));
    assert(remquo(GDN!1(2), GDN!1(-real.infinity), n) is GDN!1(2));
}


/** Calculates the remainder of `g` divided by `h` using the definition of remainder provided by
 * IEC 60559.
 *
 * The remainder is defined as $(MATH f(g(x), h(x)) = g(x) - h(x)round(g(x)/h(x))).
 *
 * Params:
 *   G = the type of g, either `GDN` or `real`.
 *   H = the type of h, either `GDN` or `real`.
 *   g = the dividend.
 *   h = the divisor.
 *
 * Returns:
 *   The remainder of `g` divided by `h`.
 */
nothrow @nogc @safe
CommonGDN!(G, H) remainder(G, H)(in G g, in H h)
if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H))
{
    int _;
    return remquo(g, h, _);
}

///
unittest
{
    import std.math: isClose;

    const f = remainder(GDN!1(5.1), GDN!1(3));
    assert(isClose(f.val, -0.9));
    assert(f.d == -1);
}
