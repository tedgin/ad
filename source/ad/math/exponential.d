/// It extends `std.math.exponential` to support `GDN` objects.
module ad.math.exponential;

static import std.math.exponential;

import std.math: LN2;

import ad.core;
import ad.math.traits: signbit;


/**
 * This function computes the base-2 logarithm of the given `GDN` object `g`. It is analogous to
 * `std.math.exponential.log2()`.
 *
 * If $(MATH f(x) = lg(g(x))), then $(MATH f' = g'/[g⋅ln(2)])
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` object for which the base-2 logarithm is to be calculated.
 *
 * Returns:
 *   A `GDN` object representing the base-2 logarithm of `g`.
 */
pure nothrow @nogc @safe GDN!Deg log2(ulong Deg)(in GDN!Deg g)
{
    const df = signbit(g) == 1 ? GDN!Deg.mkNaNDeriv : g.d/(LN2 * g.reduce());
    return GDN!Deg(std.math.exponential.log2(g.val), df);
}

///
unittest
{
    import std.math: LN2;

    assert(log2(GDN!1(2)) is GDN!1(1, 1/(2 * LN2)));
    assert(log2(GDN!1(-0.)) is GDN!1(-real.infinity, real.nan));
    assert(log2(GDN!1(+0.)) is GDN!1(-real.infinity, real.infinity));
}

unittest
{
    assert(log2(GDN!2(1)) is GDN!2(0, 1/LN2, -1/LN2));
    assert(log2(GDN!1(1, -1)) is GDN!1(0, -1/LN2));
}


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
nothrow pure @safe GDN!Deg pow(ulong Deg)(in GDN!Deg x, in real y)
{
    return x ^^ y;
}

/// ditto
nothrow pure @safe GDN!Deg pow(ulong Deg)(in real x, in GDN!Deg y)
{
    return x ^^ y;
}

/// ditto
nothrow pure @safe GDN!Deg pow(ulong Deg)(in GDN!Deg x, in GDN!Deg y)
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
nothrow pure @safe GDN!Deg exp(ulong Deg)(in GDN!Deg x)
{
    return GDN!Deg(exp(x.val), x.d * exp(x.reduce()));
}

unittest
{
    import std.format;

    assert(exp(GDN!1.nan).same(GDN!1.nan));

    assert(exp(GDN!1.zero).same(GDN!1.one));

    // XXX - derivSeq isn't implemented due to the following bug
    // https: //issues.dlang.org/show_bug.cgi?id=22621 is fixed.
    // assert(exp(GDN!1.infinity).same(derivSeq(real.infinity, real.infinity)));

    const e = exp(-GDN!1.infinity);
    assert(e.same(GDN!1(0, -0.)), format("exp(-inf) != %s", e));
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
nothrow pure @safe GDN!Deg log(ulong Deg)(in GDN!Deg x)
{
    return x.log();
}
+/
