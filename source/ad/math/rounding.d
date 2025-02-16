/// It extends `std.math.rounding` to support `GenDualNum` objects.
module ad.math.rounding;

public import std.math.rounding;

import std.math : isFinite;
static import std.math;

import ad.core;
import ad.math.dirac;


// TODO: finish implement this module.

/**
Returns the value of `g` rounded upward to the nearest integer.

If $(MATH f(x) = ⌈g(x)⌉), then $(MATH f' = g'∑$(SUB i∊ℤ)𝛿(g-i))
*/
GenDualNum!Degree ceil(ulong Degree)(in GenDualNum!Degree g) nothrow pure @nogc @safe
{
    static if (Degree == 1) {
        const f_red = std.math.ceil(g.reduce());
        const f_val = f_red;
    } else {
        const f_red = ceil(g.reduce());
        const f_val = f_red.val;
    }

    GenDualNum!Degree.DerivType!1 df;
    if (isFinite(g.val)) {
        df = g.d * dirac(g.reduce() - f_red);
    }

    return GenDualNum!Degree(f_val, df);
}

///
unittest
{
    const q = ceil(GenDualNum!1(1));
    assert(q == 1 && q.d == real.infinity);

    const w = ceil(GenDualNum!1(-1.4));
    assert(w == -1 && w.d == 0);
}

unittest
{
    import std.math : LN2;

    alias GDN = GenDualNum;

    assert(ceil(GDN!1.nan).same(GDN!1.nan));
    assert(ceil(GDN!1(real.infinity)).same(GDN!1(real.infinity, real.nan)));
    assert(ceil(GDN!1(-real.infinity)).same(GDN!1(-real.infinity, real.nan)));
    assert(ceil(GDN!2(2.3)).same(GDN!2(3, 0, 0)));
    assert(ceil(GDN!1(0, -1/LN2)).same(GDN!1(0, -real.infinity)));
}

/**
Returns the value of `g` rounded downward to the nearest integer.

If $(MATH f(x) = ⌊g(x)⌋), then $(MATH f' = g'∑$(SUB i∊ℤ)𝛿(g-i))
*/
GenDualNum!Degree floor(ulong Degree)(in GenDualNum!Degree g) nothrow pure @nogc @safe
{
    static if (Degree == 1) {
        const f_red = std.math.floor(g.reduce());
        const f_val = f_red;
     } else {
        const f_red = floor(g.reduce());
        const f_val = f_red.val;
     }

    GenDualNum!Degree.DerivType!1 df;
    if (std.math.isFinite(g.val)) {
        df = g.d * dirac(g.reduce() - f_red);
    }

    return GenDualNum!Degree(f_val, df);
}

///
unittest
{
    const q = floor(GenDualNum!1(1));
    assert(q == 1 && q.d == real.infinity);

    const w = floor(GenDualNum!1(-1.4));
    assert(w == -2 && w.d == 0);
}

unittest
{
    import std.math : LN2;

    alias GDN = GenDualNum;

    assert(floor(GDN!1.nan).same(GDN!1.nan));
    assert(floor(GDN!1(real.infinity)).same(GDN!1(real.infinity, real.nan)));
    assert(floor(GDN!1(-real.infinity)).same(GDN!1(-real.infinity, real.nan)));
    assert(floor(GDN!2(2.3)).same(GDN!2(2, 0, 0)));
    assert(floor(GDN!1(0, -1/LN2)).same(GDN!1(0, -real.infinity)));
}
