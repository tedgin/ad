/// It extends `std.math.traits` to support `GDN` objects.
module ad.math.traits;

public import std.math.traits;

static import std.math;

import ad.core;
import ad.math.dirac;

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
nothrow pure @safe bool isNaN(ulong Deg)(in GDN!Deg x)
{
    return std.math.isNaN(x.val);
}

unittest
{
    assert(isNaN(GDN!1.nan));
}
+/

/**
This function computes the sign of the argument.

If $(MATH f(x) = sgn(g(x))), then $(MATH f' = 2𝛿(g)g'), where $(MATH 𝛿) is the Dirac delta function.
*/
GDN!Deg sgn(ulong Deg)(in GDN!Deg g) nothrow pure @nogc @safe
{
    return GDN!Deg(std.math.sgn(g.val), 2*dirac(g.reduce())*g.d);
}

///
unittest
{
    const q = sgn(GDN!1(-2));
    assert(q == -1 && q.d == 0);

    const w = sgn(GDN!1(0));
    assert(w == 0 && w.d == real.infinity);
}

unittest
{
    assert(sgn(GDN!1()).same(GDN!1.nan));
    assert(sgn(GDN!1(0, real.nan)).same(GDN!1(0, real.nan)));
    assert(sgn(GDN!2(1.0L, 2.0L, real.nan)).same(GDN!2(1.0L, 0.0L, real.nan)));
    assert(sgn(GDN!1(real.infinity, 0)).same(GDN!1(1, 0)));
    assert(sgn(GDN!1(-real.infinity, 0)).same(GDN!1(-1, 0)));
    assert(sgn(GDN!1(-1)).same(GDN!1(-1, 0)));
}
