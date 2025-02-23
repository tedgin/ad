/// This module extends the `core.math` and `std.math` libraries to support `GDN` objects.
module ad.math;

public import core.math : toPrec, yl2x, yl2xp1;
public import std.math;
public import ad.core;
public import ad.math.algebraic;
public import ad.math.constants;
public import ad.math.exponential;
public import ad.math.operations;
public import ad.math.remainder;
public import ad.math.rounding;
public import ad.math.traits;
public import ad.math.trigonometry;

static import core.math;

import std.algorithm.iteration : map;
import std.array : array;
import std.traits : isFloatingPoint, isImplicitlyConvertible, Select;

import ad.math.dirac;

// core.math implementations
pragma(inline, true)
{
    /**
     * This function computes $(MATH 2$(SUP c)g).
     *
     * If $(MATH f(x) = 2$(SUP c)g(x)), then $(MATH f' = 2$(SUP c)g').
     *
     * Params:
     *   g = the generalized dual number being scaled.
     *   c = the power of $(MATH 2) used to scale `g`,
     */
    pure nothrow @nogc @safe GDN!Deg ldexp(ulong Deg)(in GDN!Deg g, in int c)
    {
        alias ldexp_red = Select!(Deg == 1, core.math.ldexp, ldexp);
        return GDN!Deg(core.math.ldexp(g.val, c), ldexp_red(g.d, c));
    }

    ///
    unittest
    {
        const f = ldexp(GDN!2(1), 2);
        assert(f == 4 && f.d == 4 && f.d!2 == 0);
    }

    /**
     * Rounds `g` to the nearest integer value, using the current rounding mode. If the return value
     * is not identical to `g`, the `FE_INEXACT` exception is raised.
     *
     * If $(MATH f(g(x)) = rint(g(x))), then $(MATH f' = (df/dg)g'), where $(MATH df/dg = 𝛿(g - m)),
     * $(MATH 𝛿) is the Dirac delta function, and $(MATH m) is a rounding mode split point.
     */
    pure nothrow @nogc @safe GDN!Deg rint(ulong Deg)(in GDN!Deg g)
    {
        const f = core.math.rint(g.val);

        auto df = GDN!Deg.mkZeroDeriv();
        if (isInfinity(g.val)) {
            df = GDN!Deg.mkNaNDeriv();
        } else {
            auto fn = nearbyint(nextDown(g.val));
            auto fp = nearbyint(nextUp(g.val));

            if (f == 0) {
                if (signbit(f) == 1)
                    fp = f;
                else
                    fn = f;
            }

            if (fn != fp) {
                df = dirac(g.reduce() - g.val);
            }
        }

        return GDN!Deg(f, df*g.d);
    }

    ///
    unittest
    {
        resetIeeeFlags();
        const e = rint(GDN!2(1.5));
        assert(ieeeFlags.inexact);
        assert(e == 2 && e.d == real.infinity && isNaN(e.d!2));
    }

    unittest
    {
        const q = rint(GDN!1.infinity);
        assert(q == real.infinity && isNaN(q.d));

        resetIeeeFlags();
        const w = rint(GDN!1.one);
        assert(!ieeeFlags.inexact);
        assert(w == 1 && w.d == 0);

        const e = rint(GDN!2(1.5));
        assert(e == 2 && e.d == real.infinity && isNaN(e.d!2));

        const r = rint(GDN!1(-0.));
        assert(r == 0 && signbit(r.val) == 1 && r.d == 0);

        const t = rint(GDN!1(+0.));
        assert(t == 0 && signbit(t.val) == 0 && t.d == 0);
    }

    /**
     * This function rounds `g` to a `long` using the current rounding mode. All of the derivative
     * terms are lost.
     */
    pure nothrow @nogc @safe long rndtol(ulong Deg)(in GDN!Deg g)
    {
        return core.math.rndtol(g.val);
    }

    ///
    unittest
    {
        assert(rndtol(GDN!1(1.1)) == 1L);
    }

    /**
     * This function rounds `g` to a given floating point type removing all derivative information.
     *
     * Params:
     *   T = the float point type to be converted to
     *   g = the generalized dual number to be converted
     */
    pure nothrow @nogc @safe T toPrec(T, ulong Deg)(in GDN!Deg g) if (isFloatingPoint!T)
    {
        return core.math.toPrec!T(g.val);
    }

    ///
    unittest
    {
        assert(typeid(toPrec!float(GDN!1.zero)) == typeid(float));
    }

    /**
     * This function computes $(MATH h⋅lg(g)). It either `g` or `h` has type `real`, it is converted
     * to a constant generalized dual number with the same degree as the other parameter.
     *
     * If $(MATH f(x) = h(x)lg(g(x))), then $(MATH f' = h'lg(g) + hg'/(ln(2)g))
     *
     * Params:
     *   g = the argument of logarithm
     *   h = the multiplier of the logarithm
     *
     * Returns:
     *   The resulting generalized dual number will have a degree equal to the lesser of the degrees
     *   of `g` and `h`.
     */
    pure nothrow @nogc @safe
    GDN!(GDeg < HDeg ? GDeg : HDeg) yl2x(ulong GDeg, ulong HDeg)(in GDN!GDeg g, in GDN!HDeg h)
    {
        return yl2x_impl(cast(typeof(return)) g, cast(typeof(return)) h);
    }

    /// ditto
    pure nothrow @nogc @safe GDN!Deg yl2x(ulong Deg)(in GDN!Deg g, in real c)
    {
        return yl2x_impl(g, GDN!Deg.mkConst(c));
    }

    /// ditto
    pure nothrow @nogc @safe GDN!Deg yl2x(ulong Deg)(in real c, in GDN!Deg h)
    {
        return yl2x_impl(GDN!Deg.mkConst(c), h);
    }

    ///
    unittest
    {
        import std.math: LN2;

        const f = yl2x(GDN!1(2), GDN!1(3));
        assert(f == 3 && f.d == 1 + 1.5/LN2);

        const x = yl2x(GDN!1(+0., -1), GDN!1(1));
        assert(x == -real.infinity && x.d == -real.infinity);

        assert(typeof(yl2x(GDN!2(1), GDN!1(2))).DEGREE == 1);

        const z = yl2x(GDN!1(1), 2.);
        assert(z == 0 && z.d == 2/LN2);
    }

    unittest
    {
        const q = yl2x(1, GDN!1(2));
        assert(q == 0 && q.d == 0);
    }

    private GDN!Deg yl2x_impl(ulong Deg)(in GDN!Deg g, in GDN!Deg h) nothrow pure @nogc @safe
    {
        alias yl2x_red = Select!(Deg == 1, core.math.yl2x, yl2x_impl);

        GDN!Deg.DerivType!1 df;
        if (std.math.signbit(g.val) == 0) {
            const g_red = g.reduce();
            df = yl2x_red(g_red, h.d) + h.reduce() * g.d / (std.math.LN2 * g_red);
        }

        return GDN!Deg(yl2x(g.val, h.val), df);
    }

    unittest
    {
        import std.math : LN2;

        const q = yl2x_impl(GDN!1(-1), GDN!1(1));
        assert(isNaN(q.val) && isNaN(q.d));

        const w = yl2x_impl(GDN!1(0), GDN!1(1));
        assert(w == -real.infinity && isNaN(w.d));

        const e = yl2x_impl(GDN!2(1), GDN!2(2));
        // f = 0
        // <f',f"> = h'lg(g) + hg'/(ln(2)g)
        //    = <1,0>lg(<1,1>) + <2,1><1,0>/ln(2)<1,1>
        //    = <0,0+1/ln(2)> + <2,1>/<1,1>/ln(2)
        //    = <0,1/ln(2)> + <2,-1>/ln(2)
        //    = <2/ln(2),0>
        assert(e.val == 0 && e.d == 2/LN2 && e.d!2 == 0);
    }

    /**
     * Computes $(MATH h⋅lg(g + 1)), for $(MATH -(1 - √½) ≤ x ≤ +(1 - √½)). When $(MATH g) is
     * outside of this interval, the results are undefined. It either `g` or `h` has type `real`, it
     * is converted to a constant generalized dual number with the same degree as the other
     * parameter.
     *
     * If $(MATH f(x) = h(x)lg(g(x) + 1)), then $(MATH f' = h'lg(g + 1) + hg'/[ln(2)(g + 1)])
     *
     * Params:
     *   g = the argument of logarithm
     *   h = the multiplier of the logarithm
     *
     * Returns:
     *   The resulting generalized dual number will have a degree equal to the lesser of the degrees
     *   of `g` and `h`.
     */
    pure nothrow @nogc @safe
    GDN!(GDeg < HDeg ? GDeg : HDeg) yl2xp1(ulong GDeg, ulong HDeg)(in GDN!GDeg g, in GDN!HDeg h)
    {
        return yl2xp1_impl(cast(typeof(return)) g, cast(typeof(return)) h);
    }

    /// ditto
    pure nothrow @nogc @safe GDN!Deg yl2xp1(ulong Deg)(in GDN!Deg g, in real c)
    {
        return yl2xp1_impl(g, GDN!Deg.mkConst(c));
    }

    /// ditto
    pure nothrow @nogc @safe GDN!Deg yl2xp1(ulong Deg)(in real c, in GDN!Deg h)
    {
        return yl2xp1_impl(GDN!Deg.mkConst(c), h);
    }

    ///
    unittest
    {
        import std.math: LN2;

        const f = yl2xp1(GDN!1(0), GDN!1(3));
        assert(f == 0 && f.d == 3/LN2);

        assert(typeof(yl2xp1(GDN!2(0), GDN!1(1))).DEGREE == 1);

        const y = yl2xp1(GDN!1(0), 1);
        assert(y == 0 && y.d == 1/LN2);
    }

    unittest
    {
        const q = yl2xp1(0, GDN!1(1));
        assert(q == 0 && q.d == 0);
    }

    private GDN!Deg yl2xp1_impl(ulong Deg)(in GDN!Deg g, in GDN!Deg h) nothrow pure @nogc @safe
    {
        alias yl2xp1_red = Select!(Deg == 1, core.math.yl2xp1, yl2xp1_impl);

        GDN!Deg.DerivType!1 df;
        if (g > -1) {
            const g_red = g.reduce();
            df = yl2xp1_red(g_red, h.d) + h.reduce() * g.d / (std.math.LN2 * (g_red + 1));
        }

        return GDN!Deg(core.math.yl2xp1(g.val, h.val), df);
    }

    unittest
    {
        import std.math : isNaN, LN2;

        const e = yl2xp1_impl(GDN!2(0), GDN!2(1));
        // f = 0
        // <f',f"> = h'lg(g+1) + hg'/[ln(2)(g+1)]
        //    = <1,0>lg(<0,1>+1) + <1,1><1,0>/[ln(2)(<0,1>+1)]
        //    = <1,0>lg<1,1> + <1,1>/[ln(2)<1,1>]
        //    = <0,1/ln(2)> + <1,0>/ln(2)
        //    = <1/ln(2),1/ln(2)>
        assert(e == 0 && e.d == 1/LN2 && e.d!2 == 1/LN2);

        assert(isNaN(yl2xp1_impl(GDN!1(-1), GDN!1(0)).d));

        const q = yl2xp1_impl(GDN!1(-2), GDN!1(0));
        assert(isNaN(q.d), "yl2xp1(2,0) should not have a derivative");
    }
}
