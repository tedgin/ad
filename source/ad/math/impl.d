/**
This module contains implementation logic used by ad.math. It is not intended to be used directly.
*/
module ad.math.impl;

import std.traits : isImplicitlyConvertible, Select;

static import core.math;
static import std.math;

import ad.core;

private
{
    mixin template InjectImpl(string name, string impl_name)
    {
        mixin("alias " ~ name ~ " = " ~ impl_name ~ ";");
    }

    mixin template InjectConvertImpl(string name, string impl_name)
    {
        mixin("
            GenDualNum!(GD < HD ? GD : HD) " ~
            name ~ "(ulong GD, ulong HD)(in GenDualNum!GD g, in GenDualNum!HD h) {
                alias D = Select!(GD < HD, GD, HD);
                return " ~ impl_name ~ "!D(cast(GenDualNum!D) g, cast(GenDualNum!D) h);
            }");

        mixin("
            GenDualNum!D " ~ name ~ "(T, ulong D)(in GenDualNum!D g, in T c)
            if (isImplicitlyConvertible!(T, real)) {
                return " ~ impl_name ~ "(g, GenDualNum!D.mkConst(c));
            }");

        mixin("
            GenDualNum!D " ~ name ~ "(T, ulong D)(in T c, in GenDualNum!D h)
            if (isImplicitlyConvertible!(T, real)) {
                return " ~ impl_name ~ "(GenDualNum!D.mkConst(c), h);
            }");
    }

    mixin template InjectCoreAlias(string name)
    {
        mixin("alias " ~ name ~ " = core.math." ~ name ~ ";");
    }

    mixin template InjectStdAlias(string name)
    {
        mixin("alias " ~ name ~ " = std.math." ~ name ~ ";");
    }
}

package pragma(inline, true) nothrow pure @nogc @safe
{
    // COS

    private GenDualNum!Degree cos_impl(ulong Degree)(in GenDualNum!Degree g)
    {
        return GenDualNum!Degree(cos(g.val), -sin(g.reduce()) * g.d);
    }

    unittest
    {
        alias GDN = GenDualNum;

        const g = cos_impl(GDN!1(std.math.PI_2));
        assert(std.math.isClose(g.val, 0., 0., real.epsilon) && g.d == -1);

        assert(cos_impl(GDN!1.infinity).same(GDN!1.nan));
        assert(cos_impl(-GDN!1.infinity).same(GDN!1.nan));
    }

    mixin InjectImpl!("cos", "cos_impl");
    mixin InjectCoreAlias!"cos";


    // LDEXP

    private GenDualNum!Degree ldexp_impl(ulong Degree)(in GenDualNum!Degree g, in int exp)
    {
        return GenDualNum!Degree(ldexp(g.val, exp), ldexp(g.d, exp));
    }

    mixin InjectImpl!("ldexp", "ldexp_impl");
    mixin InjectCoreAlias!"ldexp";


    // SIN

    private GenDualNum!Degree sin_impl(ulong Degree)(in GenDualNum!Degree g)
    {
        return GenDualNum!Degree(sin(g.val), cos(g.reduce()) * g.d);
    }

    unittest
    {
        alias GDN = GenDualNum;

        assert(sin_impl(GDN!1.zero).same(GenDualNum!1.zero));

        const g = sin_impl(GDN!1(std.math.PI_2));
        assert(g.val == 1 && std.math.isClose(g.d, 0., 0., real.epsilon));

        assert(sin_impl(GDN!1.infinity).same(GDN!1.nan));
        assert(sin_impl(-GDN!1.infinity).same(GDN!1.nan));
    }

    mixin InjectImpl!("sin", "sin_impl");
    mixin InjectCoreAlias!"sin";


    // YL2X

    private GenDualNum!Deg yl2x_impl(ulong Deg)(in GenDualNum!Deg g, in GenDualNum!Deg h)
    {
        GenDualNum!Deg.DerivType!1 df;
        if (std.math.signbit(g.val) == 1)
        {
            df = GenDualNum!Deg.mkNaNDeriv();
        }
        else
        {
            const g_red = g.reduce();
            df = yl2x(g_red, h.d) + h.reduce() * g.d / (std.math.LN2 * g_red);
        }

        return GenDualNum!Deg(yl2x(g.val, h.val), df);
    }

    unittest
    {
        alias GDN = GenDualNum;

        assert(yl2x_impl(GDN!1(-1), GDN!1(1)).same(GDN!1.nan));
        assert(yl2x_impl(GDN!1(-0.), GDN!1(1)).same(GDN!1(-real.infinity, real.nan)));

        const e = yl2x_impl(GDN!2(1), GDN!2(2));
        // f = 0
        // <f',f"> = h'lg(g) + hg'/(ln(2)g)
        //    = <1,0>lg(<1,1>) + <2,1><1,0>/ln(2)<1,1>
        //    = <0,0+1/ln(2)> + <2,1>/<1,1>/ln(2)
        //    = <0,1/ln(2)> + <2,-1>/ln(2)
        //    = <2/ln(2),0>
        assert(e.same(GDN!2(0, 2 / std.math.LN2, 0)));
    }

    mixin InjectConvertImpl!("yl2x", "yl2x_impl");
    mixin InjectCoreAlias!"yl2x";


    // YL2XP1

    private GenDualNum!Deg yl2xp1_impl(ulong Deg)(in GenDualNum!Deg g, in GenDualNum!Deg h)
    {
        GenDualNum!Deg.DerivType!1 df;
        if (g <= -1)
        {
            df = GenDualNum!Deg.mkNaNDeriv();
        }
        else
        {
            const g_red = g.reduce();
            df = yl2xp1(g_red, h.d) + h.reduce() * g.d / (std.math.LN2 * (g_red + 1));
        }

        return GenDualNum!Deg(yl2xp1(g.val, h.val), df);
    }

    unittest
    {
        alias GDN = GenDualNum;

        const e = yl2xp1_impl(GDN!2(0), GDN!2(1));
        // f = 0
        // <f',f"> = h'lg(g+1) + hg'/[ln(2)(g+1)]
        //    = <1,0>lg(<0,1>+1) + <1,1><1,0>/[ln(2)(<0,1>+1)]
        //    = <1,0>lg<1,1> + <1,1>/[ln(2)<1,1>]
        //    = <0,1/ln(2)> + <1,0>/ln(2)
        //    = <1/ln(2),1/ln(2)>
        assert(e.same(GDN!2(0, 1 / std.math.LN2, 1 / std.math.LN2)));

        assert(std.math.isNaN(yl2xp1_impl(GDN!1(-1), GDN!1(0)).d));

        const q = yl2xp1_impl(GDN!1(-2), GDN!1(0));
        assert(std.math.isNaN(q.d), "yl2xp1(2,0) should not have a derivative");
    }

    mixin InjectConvertImpl!("yl2xp1", "yl2xp1_impl");
    mixin InjectCoreAlias!"yl2xp1";


    // HYPOT

    private GenDualNum!Deg hypot2_impl(ulong Deg)(in GenDualNum!Deg g, in GenDualNum!Deg h)
    {
        const g_red = g.reduce();
        const h_red = h.reduce();
        const f_red = hypot(g_red, h_red);
        const df = (g_red * g.d + h_red * h.d) / f_red;

        static if (Deg == 1) return GenDualNum!1(f_red, df);
        else return GenDualNum!Deg(f_red.val, df);
    }

    unittest
    {
        alias GDN = GenDualNum;

        const f = hypot2_impl(GDN!2(1), GDN!2(2));
        // f = sqrt(5)
        // <f',f"> = (<1,1><1,0> + <2,1><1,0>) / <sqrt(5),3/sqrt(5)>
        //    = (<1,1> + <2,1>) / <sqrt(5),3/sqrt(5)>
        //    = <3,2>/<sqrt(5),3/sqrt(5)>
        //    = <3/sqrt(5) , (2sqrt(5) - 9/sqrt(5))/5>
        //    = <3/sqrt(5) , .4sqrt(5) - 1.8/sqrt(5)>
        const w = std.math.sqrt(5.0L);
        const q = GDN!2(w, 3 / w, .4 * w - 1.8 / w);
        assert(std.math.isClose(f.d!2, q.d!2));
    }

    private
    GenDualNum!Deg
    hypot3_impl(ulong Deg)(in GenDualNum!Deg g, in GenDualNum!Deg h, in GenDualNum!Deg i)
    {
        const g_red = g.reduce();
        const h_red = h.reduce();
        const i_red = i.reduce();
        const f_red = hypot(g_red, h_red, i_red);
        const df = (g_red * g.d + h_red * h.d + i_red * i.d) / f_red;

        static if (Deg == 1) return GenDualNum!1(f_red, df);
        else return GenDualNum!Deg(f_red.val, df);
    }

    unittest
    {
        alias GDN = GenDualNum;

        const f = hypot3_impl(GDN!2(1), GDN!2(2), GDN!2(3));
        // f = sqrt(14)
        // <f',f"> = (<1,1><1,0> + <2,1><1,0> + <3,1><1,0>) / <sqrt(14),6/sqrt(14)>
        //    = (<1,1> + <2,1> + <3,1>) / <sqrt(14),6/sqrt(14)>
        //    = <6,3>/<sqrt(14),6/sqrt(14)>
        //    = <6/sqrt(14) , (3sqrt(14) - 36/sqrt(14))/14>
        //    = <6/sqrt(14) , 3/sqrt(14) - 18/(7sqrt(14))>
        const w = std.math.sqrt(14.0L);
        const q = GDN!2(w, 6 / w, 3 / w - 18 / (7 * w));
        assert(std.math.isClose(f.d!2, q.d!2));
    }

    mixin InjectConvertImpl!("hypot", "hypot2_impl");
    mixin InjectImpl!("hypot", "hypot3_impl");
    mixin InjectStdAlias!"hypot";
}
