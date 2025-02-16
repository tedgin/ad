/// It extends `std.math.exponential` to support `GenDualNum` objects.
module ad.math.exponential;

public import std.math.exponential;

import std.math : LN2, signbit;

import ad.core;


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

/**
Calculates the base-2 logarithm of `g`.

If $(MATH f(x) = lg(g(x))), then $(MATH f' = g'/[g⋅ln(2)])
*/
GenDualNum!Degree log2(ulong Degree)(in GenDualNum!Degree g) nothrow pure @nogc @safe
{
    const df = signbit(g.val) == 1 ? GenDualNum!Degree.mkNaNDeriv : g.d / (LN2 * g.reduce());
    return GenDualNum!Degree(std.math.log2(g.val), df);
}

///
unittest
{
    import std.math : isNaN, LN2;

    const q = log2(GenDualNum!1(2));
    assert(q == 1 && q.d == 1 / (2 * LN2));

    const w = log2(GenDualNum!1(-0.));
    assert(w == -real.infinity && isNaN(w.d));

    const e = log2(GenDualNum!1(+0.));
    assert(e == -real.infinity && e.d == real.infinity);
}

unittest
{
    alias GDN = GenDualNum;

    const q = log2(GDN!2(1));
    assert(q.same(GDN!2(0, 1/LN2, -1/LN2)));

    assert(log2(GDN!1(1, -1)).same(GDN!1(0, -1/LN2)));
}
