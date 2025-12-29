// This module implements functions used internally by the `ad.math` module.
module ad.math.internal;

static import core.math;
static import std.math.exponential;
static import std.math.operations;
static import std.math.rounding;
static import std.math.traits;

import std.algorithm: min;
import std.math: abs, isFinite, isNaN, LN2;
import std.traits: fullyQualifiedName, isImplicitlyConvertible, isIntegral, Select, TemplateOf;

import ad.core;


/*
 * General shared internals
 */
package pure nothrow @nogc @safe
{
    // Exposes the GDN.dirac method as a function to help with overload resolution.
    pragma(inline, true) GDN!Deg dirac(ulong Deg)(in GDN!Deg g)
    {
        return g.dirac();
    }

    // The Dirac delta function for a `real`.
    real dirac(in real g)
    {
        if (isNaN(g)) {
            return real.nan;
        }

        return g == 0 ? real.infinity : 0;
    }

    unittest
    {
        assert(dirac(0) == real.infinity);
        assert(dirac(1) == 0);
        assert(isNaN(dirac(real.nan)));
    }

    // Determines the minimum GDN degree for a sequence of types. Only GDN types are considered. If
    // no GDN types are present, the result is 0.
    enum minDeg(T...) = {
        bool foundGDN = false;
        ulong commonDeg = ulong.max;
        static foreach (U; T) {
            static if (isGDN!U) {
                foundGDN = true;
                commonDeg = min(commonDeg, U.DEGREE);
            }
        }
        return foundGDN ? commonDeg : 0;
    }();

    unittest
    {
        static assert(minDeg!(GDN!2) == 2);
        static assert(minDeg!char == 0);
        static assert(minDeg!(GDN!32, GDN!2, real) == 2);
    }
}


/*
 * std.math.traits shared internals
 */
package pure nothrow @nogc @safe
{
    // The implementation of areAll
    enum bool areAll(alias Test, TS...) = {
        auto res = true;
        static foreach(T; TS)
            static if (!Test!T) res = false;
        return res;
    }();

    // The implementation of areNone
    enum bool areNone(alias Test, TS...) = {
        auto res = true;
        static foreach(T; TS)
            static if (Test!T) res = false;
        return res;
    }();

    // The implementation of isGDN
    enum bool isGDN(T) = fullyQualifiedName!(TemplateOf!T) == "ad.core.GDN";

    // The implementation of isGDNOrReal
    enum bool isGDNOrReal(T) = isGDN!T || isImplicitlyConvertible!(T, real);

    // Converts a GDN to a GDN of the specified degree.
    pragma(inline, true) GDN!Deg asGDN(ulong Deg, T)(in T t) if (isGDN!T)
    {
        return cast(GDN!Deg) t;
    }

    unittest
    {
        assert(asGDN!1(GDN!2(3)) is GDN!1(3));
        assert(asGDN!2(GDN!2(3)) is GDN!2(3));
        assert(asGDN!3(GDN!2(3)) is GDN!3(3));
    }

    // Converts a real to a constant GDN of the specified degree.
    pragma(inline, true) GDN!Deg asGDN(ulong Deg)(in real t)
    {
        return GDN!Deg.mkConst(t);
    }

    unittest
    {
        assert(asGDN!2(3) is GDN!2(3, 0, 0));
    }


    // Converts a GDN or something implicitly convertible to a real to a real.
    pragma(inline, true) real asReal(F)(in F f) if (isGDNOrReal!F)
    {
        static if (isGDN!F)
            return f.val;
        else
            return f;
    }

    unittest
    {
        assert(asReal(GDN!2(3)) is 3.0L);
        assert(asReal(1) is 1.0L);
    }


    // The implementation of isFinite
    pragma(inline, true) bool isFinite(ulong Deg)(in GDN!Deg f)
    {
        return std.math.traits.isFinite(f.val);
    }


    // The implementation of isInfinity
    pragma(inline, true) bool isInfinity(ulong Deg)(in GDN!Deg f)
    {
        return std.math.traits.isInfinity(f.val);
    }


    // The implementation of isNaN
    pragma(inline, true) bool isNaN(ulong Deg)(in GDN!Deg f)
    {
        return std.math.traits.isNaN(f.val);
    }


    // The implementation of isOne
    enum bool isOne(alias Test, TS...) = {
        auto res = false;
        static foreach(T; TS)
            static if (Test!T) res = true;
        return res;
    }();


    // The implementation of CommonGDN
    template CommonGDN(G...) if (isOne!(isGDN, G) && areAll!(isGDNOrReal, G)) {
        alias CommonGDN = GDN!(minDeg!G);
    }


    // The implementation of sgn
    pragma(inline, true) GDN!Deg sgn(ulong Deg)(in GDN!Deg g)
    {
        return GDN!Deg(std.math.traits.sgn(g.val), 2*dirac(g.reduce())*g.d);
    }

    unittest
    {
        assert(sgn(GDN!1()) is GDN!1.nan);
        assert(sgn(GDN!1(0, real.nan)) is GDN!1(0, real.nan));
        assert(sgn(GDN!2(1.0L, 2.0L, real.nan)) is GDN!2(1.0L, 0.0L, real.nan));
        assert(sgn(GDN!1(real.infinity, 0)) is GDN!1(1, 0));
        assert(sgn(GDN!1(-real.infinity, 0)) is GDN!1(-1, 0));
        assert(sgn(GDN!1(-1)) is GDN!1(-1, 0));
    }

    // The implementation of signbit
    pragma(inline, true) int signbit(ulong Deg)(in GDN!Deg f)
    {
        return std.math.traits.signbit(f.val);
    }
}


/*
 * std.math.algebraic shared internals
 */
package pure @safe
{
    pragma(inline, true) nothrow @nogc GDN!Deg sqrt(ulong Deg)(in GDN!Deg g)
    {
        if (isNaN(g.val)) return g;

        const dfdg = signbit(g) == 1 ? GDN!Deg.DerivType!1.nan : g.reduce()^^-0.5/2;
        return GDN!Deg(core.math.sqrt(g.val), dfdg * g.d);
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
}


/*
 * std.math.exponential shared internals
 */
package pure nothrow @nogc @safe
{
    // The implementation of exp for GDN.
    pragma(inline, true) GDN!Deg exp(ulong Deg)(in GDN!Deg g)
    {
        static if (Deg == 1)
            alias exp_fn = std.math.exponential.exp;
        else
            alias exp_fn = exp;

        const f_red = exp_fn(g.reduce());
        return GDN!Deg(asReal(f_red), f_red * g.d);
    }

    unittest
    {
        assert(exp(GDN!1.nan) is GDN!1.nan);
        assert(exp(GDN!1.zero) is GDN!1.one);
        assert(exp(-GDN!2.infinity) is GDN!2(0, -0., 0));
    }

    // The implementation of log2 for GDN.
    pragma(inline, true) GDN!Deg log2(ulong Deg)(in GDN!Deg g)
    {
        if (isNaN(g.val)) return g;
        const df = signbit(g) == 1 ? GDN!Deg.mkNaNDeriv : 1.0L / (LN2*g.reduce());
        return GDN!Deg(std.math.exponential.log2(g.val), df * g.d);
    }

    unittest
    {
        assert(log2(GDN!2(1)) is GDN!2(0, 1/LN2, -1/LN2));
        assert(log2(GDN!1(1, -1)) is GDN!1(0, -1/LN2));
    }

    // The implement of pow for GDN where the exponent is an integer
    pragma(inline, true) GDN!Deg pow(I, ulong Deg)(in GDN!Deg g, in I n) if (isIntegral!I)
    {
        alias pow_red = Select!(Deg == 1, std.math.exponential.pow, pow);

        if (isNaN(g)) return g;

        const df = n == 0 ? GDN!Deg.DerivType!1(0.0L) : n * pow_red(g.reduce(), n-1);
        return GDN!Deg(std.math.exponential.pow(g.val, n), df * g.d);
    }

    unittest
    {
        import std.format: format;

        assert(isNaN(pow(GDN!1.nan, 1)));
        assert(pow(GDN!1(2, 3), 0) is GDN!1(1, 0));

        const q = pow(GDN!1(3, real.infinity), 0);
        assert(q == 1 && isNaN(q.d));

        assert(pow(GDN!1(-2, 2), -1) is GDN!1(-0.5, -0.5));
        assert(pow(GDN!2(3), 2) is GDN!2(9, 6, 2));
    }
}


/*
 * std.math.operations shared internals
 */
package pure nothrow @nogc @safe
{
    // The implementation of getNaNPayload for GDN.
    pragma(inline, true) ulong getNaNPayload(ulong Deg)(in GDN!Deg f)
    {
        return std.math.operations.getNaNPayload(f.val);
    }

    // The implementation of nextDown for GDN.
    pragma(inline, true) GDN!Deg nextDown(ulong Deg)(in GDN!Deg g)
    {
        return GDN!Deg(std.math.operations.nextDown(g.val), g.d);
    }

    unittest
    {
        assert(isNaN(nextDown(GDN!1.nan)));
    }


    // The implementation of nextUp for GDN.
    pragma(inline, true) GDN!Deg nextUp(ulong Deg)(in GDN!Deg g)
    {
        return GDN!Deg(std.math.operations.nextUp(g.val), g.d);
    }

    unittest
    {
        assert(isNaN(nextUp(GDN!1.nan)));
    }
}


/*
 * std.math.rounding shared internals
 */
package pure nothrow @nogc @safe
{
    pragma(inline,true) private
    GDN!Deg
    round_impl(string Fn, ulong Deg)(
        in GDN!Deg g,
        in GDN!Deg.DerivType!1 delegate(in GDN!Deg.DerivType!1, in GDN!Deg.DerivType!1)
            pure nothrow @nogc @safe
            fn_deriv)
    {
        if (isNaN(g.val)) return g;

        static if (Deg == 1)
            mixin("const f_red = std.math.rounding." ~ Fn ~ "(g.reduce());");
        else
            mixin("const f_red = " ~ Fn ~ "(g.reduce());");

        GDN!Deg.DerivType!1 df;
        if (isFinite(g)) {
            df = fn_deriv(f_red, g.reduce()) * g.d;
        }

        return GDN!Deg(asReal(f_red), df);
    }

    unittest
    {
        real sub(in real f, in real g) { return f - g; }
        assert(round_impl!("ceil", 1)(GDN!1(3.5), &sub) is GDN!1(4, 0.5));

        auto add(in GDN!1 f, in GDN!1 g) { return f + g; }
        const q = round_impl!("ceil", 2)(GDN!2(3.5), &add);
        assert(q is GDN!2(4, 7.5, 1));
    }


    // The implementation of ceil for GDN.
    pragma(inline, true) GDN!Deg ceil(ulong Deg)(in GDN!Deg g)
    {
        auto dfn(in GDN!Deg.DerivType!1 f_red, in GDN!Deg.DerivType!1 g_red) {
            return dirac(g_red - f_red);
        }

        return round_impl!("ceil", Deg)(g, &dfn);
    }

    unittest
    {
        import std.math: LN2;

        assert(ceil(GDN!1.nan) is GDN!1.nan);
        assert(ceil(GDN!1(real.infinity)) is GDN!1(real.infinity, real.nan));
        assert(ceil(GDN!1(-real.infinity)) is GDN!1(-real.infinity, real.nan));
        assert(ceil(GDN!2(2.3)) is GDN!2(3, 0, 0));
        assert(ceil(GDN!1(0, -1/LN2)) is GDN!1(0, -real.infinity));
    }


    // The implementation of floor for GDN
    pragma(inline, true) GDN!Deg floor(ulong Deg)(in GDN!Deg g)
    {
        auto dfn(in GDN!Deg.DerivType!1 f_red, in GDN!Deg.DerivType!1 g_red) {
            return dirac(g_red - f_red);
        }

        return round_impl!("floor", Deg)(g, &dfn);
    }

    unittest
    {
        import std.math: LN2;

        assert(floor(GDN!1.nan) is GDN!1.nan);
        assert(floor(GDN!1(real.infinity)) is GDN!1(real.infinity, real.nan));
        assert(floor(GDN!1(-real.infinity)) is GDN!1(-real.infinity, real.nan));
        assert(floor(GDN!2(2.3)) is GDN!2(2, 0, 0));
        assert(floor(GDN!1(0, -1/LN2)) is GDN!1(0, -real.infinity));
    }


    // The implementation of round for GDN
    pragma(inline, true) GDN!Deg round(ulong Deg)(in GDN!Deg g)
    {
        auto dfn(in GDN!Deg.DerivType!1 f_red, in GDN!Deg.DerivType!1 g_red) {
            return dirac(abs(g_red - f_red) - 0.5);
        }

        return round_impl!("round", Deg)(g, &dfn);
    }

    unittest
    {
        assert(round(GDN!2(5.4)) is GDN!2(5, 0, 0));
    }


    // The implementation of trunc for GDN
    pragma(inline, true) GDN!Deg trunc(ulong Deg)(in GDN!Deg g)
    {
        auto dfn(in GDN!Deg.DerivType!1 f_red, in GDN!Deg.DerivType!1 g_red) {
            return dirac(abs(g_red - f_red));
        }

        return round_impl!("trunc", Deg)(g, &dfn);
    }

    unittest
    {
        assert(trunc(GDN!2(0.5)) is GDN!2(0, 0, 0));
    }
}
