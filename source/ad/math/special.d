/// It extends the `std.mathspecial` module to support `GDN` objects.

module ad.math.special;

public import std.mathspecial;

import std.math: getNaNPayload, isInfinity, isNaN, signbit, trunc;
import std.traits: select;

static import ad.math.polygamma;

import ad.core;
import ad.math.internal:
    areAll, asGDN, asReal, CommonGDN, dirac, getNaNPayload, isGDN, isGDNOrReal, isOne, isNaN,
    signbit, sgn;


private pure nothrow @nogc @safe GDN!Deg polygamma(ulong N, ulong Deg)(in GDN!Deg g) if (N > 0)
{
    static if (Deg == 1)
        const df = ad.math.polygamma.polygamma!(N+1)(g.reduce());
    else
        const df = polygamma!(N+1)(g.reduce());

    return GDN!Deg(ad.math.polygamma.polygamma!N(g.val), df*g.d);
}

unittest
{
    import std.format: format;

    const q = GDN!1(ad.math.polygamma.polygamma!1(1), ad.math.polygamma.polygamma!2(1));
    assert(polygamma!1(GDN!1(1)) is q);

    const a_exp = GDN!1(ad.math.polygamma.polygamma!1(2), 3*ad.math.polygamma.polygamma!2(2));
    const a_act = polygamma!1(GDN!1(2, 3));
    assert(a_act is a_exp, format("Ψ₁(<2,3>) = %s != %s", a_act, a_exp));

    assert(polygamma!1(GDN!1(+0.)) is GDN!1(real.infinity, -real.infinity));
    assert(polygamma!1(GDN!1(-0.)) is GDN!1(real.infinity, real.infinity));

    const u = polygamma!1(GDN!1(-1));
    assert(u == real.infinity && isNaN(u.d));

    const w = polygamma!1(GDN!1(real.nan));
    assert(isNaN(w.val) && isNaN(w.d));

    const e = polygamma!1(GDN!1(1, real.nan));
    assert(e == ad.math.polygamma.polygamma!1(1) && isNaN(e.d));

    assert(polygamma!1(GDN!1(real.infinity)) == GDN!1(+0., -0.));

    const t = polygamma!1(GDN!1(-real.infinity));
    assert(t.val == real.infinity && isNaN(t.d), format("Ψ₁(-∞) != %s", t));

    const y = GDN!1(ad.math.polygamma.polygamma!2(2), ad.math.polygamma.polygamma!3(2));
    assert(polygamma!2(GDN!1(2)) is y);

    assert(polygamma!2(GDN!1(+0.)) is GDN!1(-real.infinity, real.infinity));
    assert(polygamma!2(GDN!1(-0.)) is GDN!1(real.infinity, real.infinity));

    const i = polygamma!2(GDN!1(-1));
    assert(isNaN(i.val) && isNaN(i.d));

    const o_exp = GDN!2(
        ad.math.polygamma.polygamma!1(3),
        ad.math.polygamma.polygamma!2(3),
        ad.math.polygamma.polygamma!3(3));
    // <f',f"> = <1,0>Ψ₂(<3,1>) = <1,0><Ψ₂(3),Ψ₃(3)> = <Ψ₂(3),Ψ₃(3)>
    const o_act = polygamma!1(GDN!2(3));
    assert(o_act is o_exp, format("Ψ₁(3) = %s != %s", o_exp, o_act));
}


/**
 * the gamma function, $(MATH Γ), of a generalized dual number
 *
 * If $(MATH f(x) = Γ(g(x))), then $(MATH f' = Γ(g)Ψ(g)g')
 *
 * Params:
 *   Deg = the degree of g
 *   g = the `GDN`argument
 *
 * Returns:
 *   $(MATH Γ(g)) as a `GDN`
 */
pure nothrow @nogc @safe GDN!Deg gamma(ulong Deg)(in GDN!Deg g)
{
    const g_red = g.reduce();

    static if (Deg == 1) {
        const f_red = std.mathspecial.gamma(g_red);
        const psi = std.mathspecial.digamma(g_red);
    } else {
        const f_red = gamma(g_red);
        const psi = digamma(g_red);
    }

    return GDN!Deg(asReal(f_red), f_red*psi*g.d);
}

///
unittest
{
    import std.mathspecial: digamma;

    assert(gamma(GDN!1(1)) is GDN!1(1, digamma(1)));
}

unittest
{
    import std.format: format;

    assert(gamma(GDN!1(2, 3)) is GDN!1(1, 3*std.mathspecial.digamma(2)));

    const q = gamma(GDN!1(real.nan));
    assert(isNaN(q.val) && isNaN(q.d));

// NB: Fails because of https://github.com/dlang/phobos/issues/10802, fixed on stable
//     const t = gamma(GDN!1(-0.));
//     assert(t is GDN!1(-real.infinity, -real.infinity), format("Γ(-0) = %s", t));

// NB: Fails because of https://github.com/dlang/phobos/issues/10802, fixed on stable
//     const y = gamma(GDN!1(+0.));
//     assert(y is GDN!1(real.infinity, -real.infinity), format("Γ(+0) = %s", y));

    const w = gamma(GDN!1(-2));
    assert(isNaN(w.val) && isNaN(w.d));

    assert(gamma(GDN!1(+real.infinity)) is GDN!1(+real.infinity, +real.infinity));

    const e = gamma(GDN!1(-real.infinity));
    assert(isNaN(e.val) && isNaN(e.d));

    const r = GDN!2(
        2,
        2*std.mathspecial.digamma(3),
        2*std.mathspecial.digamma(3)^^2 + 2*ad.math.polygamma.polygamma!1(3));
    assert(gamma(GDN!2(3)) is r);
    // <f',f"> = Γ(<3,1>)Ψ(<3,1>)<1,0>
    //         = <2,2Ψ(3)><Ψ(3),Ψ₁(3)><1,0>
    //         = <2Ψ(3),2Ψ²(3)+2Ψ₁(3)>
}


/**
 * Computes the natural logarithm of the gamma function for generalized dual number.
 *
 * If $(MATH f(x) = ln|Γ(g(x))|), $(MATH f' = Ψ(g)g').
 *
 * Params:
 *   Deg = the degree of g
 *   g = the `GDN` argument
 *
 * Returns:
 *   $(MATH ln|Γ(g)|) as a `GDN`
 */
pure nothrow @nogc @safe GDN!Deg logGamma(ulong Deg)(in GDN!Deg g)
{
    static if (Deg == 1)
        const df = std.mathspecial.digamma(g.reduce());
    else
        const df = digamma(g.reduce());

    return GDN!Deg(std.mathspecial.logGamma(g.val), df*g.d);
}

///
unittest
{
    import std.mathspecial: digamma;

    assert(logGamma(GDN!1(2)) is GDN!1(0, digamma(2)));
}

unittest
{
    import std.format: format;
    import std.math: log;

    assert(logGamma(GDN!1(3, 4)) is GDN!1(log(2.0L), 4*std.mathspecial.digamma(3)));

    const q = logGamma(GDN!1.nan);
    assert(isNaN(q.val) && isNaN(q.d));

    const w = logGamma(GDN!1(-1));
    assert(w == real.infinity && isNaN(w.d));

    const e = logGamma(GDN!1(-real.infinity));
    assert(e == real.infinity && isNaN(e.d), format("logGamma(-inf) != %s", e));

    assert(logGamma(GDN!1(real.infinity)) is GDN!1(real.infinity, real.infinity));

    const r = GDN!2(log(24.0L), std.mathspecial.digamma(5), ad.math.polygamma.polygamma!1(5));
    assert(logGamma(GDN!2(5)) == r);
    // <f',f"> = Ψ(<5,1>)*<1,0> = <Ψ(5), Ψ₁(5)>
}


/**
 * Computes the sign of the gamma function of a generalized dual number.
 *
 * If $(MATH f(x) = sgn(Γ(g(x)))), then $(MATH f' = 2𝛿(Γ(g))Γ(g)Ψ(g)g'). Since $(MATH Γ(g) ≠ 0, ∀g),
 * $(MATH f' = 0), if it exists. It doesn't exists when $(MATH g) is a non-positive integer or
 * $(MATH -∞) or when $(MATH g') is infinite.
 *
 * Params:
 *   Deg = the degree of g
 *   g = the `GDN` argument
 *
 * Returns:
 *   It returns $(MATH  sgn(Γ(g))) as a `GDN`.
 */
pure nothrow @nogc @safe GDN!Deg sgnGamma(ulong Deg)(in GDN!Deg g)
{
    const f = std.mathspecial.sgnGamma(g.val);

    real df;

    if (!isNaN(f)) {
        if (g.val < 0) {
            ulong ngz = cast(ulong) trunc(-g.val);
            if (ngz != -g.val) df = (ngz & 1) == 0 ? -0. : +0.;
        } else if (g.val is -0.0L) {
            df = -0.;
        } else {
            df = +0.;
        }
    }

    return GDN!Deg(f, df * g.d);
}

///
unittest
{
    const f = sgnGamma(GDN!1(1));
    assert(f == 1 && f.d == 0);
}

unittest
{
// NB: This fails because of https://github.com/dlang/phobos/issues/10801, fixed in stable
//     const g = sgnGamma(GDN!1(-0.5));
//     assert(g == -1 && g.d == 0);

    const h = sgnGamma(GDN!1(+0.));
    assert(h is GDN!1(1, 0));

    const i = sgnGamma(GDN!1(-1));
    assert(isNaN(i.val) && isNaN(i.d));
}


/**
 * Computes the beta function where at least one of the arguments is a generalized dual number.
 *
 * If $(MATH f(x) = B(g(x),h(x))), then $(MATH f' = $(SUP ∂B)/$(SUB ∂g)g' + $(SUP ∂B)/$(SUB ∂h)h').
 * $(MATH $(SUP ∂B(x$(SUB 1),x$(SUB 2)))/$(SUB ∂x$(SUB i)) = B⋅[Ψ(x$(SUB i)) - Ψ(x$(SUB 1)+x$(SUB 2))]),
 * so $(MATH f' = B⋅[Ψ(g) - Ψ(g+h)]g' + B⋅[Ψ(h) - Ψ(g+h)]h'). This reduces to
 * $(MATH f' = B(g,h)[Ψ(g)g' + Ψ(h)h' - Ψ(g+h)(g' + h')]).
 *
 * Params:
 *   G = the first `GDN` argument
 *   H = the second `GDN` argument
 *   g = the first `GDN` argument
 *   h = the second `GDN` argument
 *
 * Returns:
 *   $(MATH B(g,h)) as a `GDN`.
 */
pure nothrow @nogc @safe
CommonGDN!(G, H) beta(G, H)(in G g, in H h) if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H))
{
    alias Deg = typeof(return).DEGREE;

    static if (Deg == 1) {
        alias B = std.mathspecial.beta;
        alias psi = std.mathspecial.digamma;
    } else {
        alias B = beta;
        alias psi = digamma;
    }

    const gg = asGDN!Deg(g);
    const hh = asGDN!Deg(h);

    const g_red = gg.reduce();
    const h_red = hh.reduce();
    const f_red = B(g_red, h_red);

    return GDN!Deg(
        asReal(f_red), f_red*(gg.d*psi(g_red) + hh.d*psi(h_red) - (gg.d + hh.d)*psi(g_red+h_red)));
}

///
unittest
{
    assert(beta(GDN!1(2), GDN!1(1)) == GDN!1(0.5, -1));
}

unittest
{
    // f' = B(g,h)[g'Ψ(g) + h'Ψ(h) - (g' + h')Ψ(g+h)]
    // In the following m and n are positive integers
    // B(m,n) = (m - 1)!(n - 1)!/(m + n - 1)!
    // Ψ(n) = -γ + ∑ᵢ₌₁ⁿ⁻¹(1/i)
    // Ψ₁(n) = 𝜋²/6 - ∑ᵢ₌₁ⁿ⁻¹(1/i²)

    import std.format: format;
    import std.math: isClose;

    assert(beta(GDN!1(1), 2) is GDN!1(0.5, -0.75));

    const q_act = beta(3, GDN!1(4));
    const q_exp = GDN!1(1.0L/60, -37.0L/3_600);
    assert(q_act == q_exp);
    assert(isClose(q_act.d, q_exp.d), format("B'(3, <4,1>) = %s != %s", q_act.d, q_exp.d));
    // f = 2!3!/6! = 2/(6*5*4) = 1/(5*4*3) = 1/60
    // f' = 1/60[Ψ(4) - Ψ(7)]
    //    = [(1 + 1/2 + 1/3) - (1 + 1/2 + 1/3 + 1/4 + 1/5 + 1/6)]/60
    //    = -(1/4 + 1/5 + 1/6)/60
    //    = -(9/20 + 1/6)/60
    //    = -(27/60 + 10/60)/60
    //    = -37/60/60
    //    = -37/3600

    assert(beta(GDN!2(2), GDN!2(3)) == GDN!2(1.0L/12, -5.0L/36, 9.0L/32));
    // f = 1!2!/4! = 2/(4*3*2) = 1/12;
    // <f',f"> = B(<2,1>, <3,1>)[<1,0>Ψ(<2,1>) + <1,0>Ψ(<3,1>) - (<1,0> + <1,0>)Ψ(<2,1>+<3,1>)]
    // B(<2,1>, <3,1>) = <B(2,3), B(2,3)[Ψ(2) + Ψ(3) - (1 + 1)Ψ(2+3)]>
    //    = <B(2,3), B(2,3)[Ψ(2) + Ψ(3) - 2Ψ(5)]>
    // B(2,3) = 1/12
    // Ψ(<n,1>) = <Ψ(n),Ψ₁(n)>
    // <1,0>Ψ(<2,1>) + <1,0>Ψ(<3,1>) - (<1,0> + <1,0>)Ψ(<2,1>+<3,1>)
    //    = <Ψ(2),Ψ₁(2)> + <Ψ(3),Ψ₁(3)> - <2,0><Ψ(5),Ψ₁(5)>
    //    = <Ψ(2) + Ψ(3), Ψ₁(2) + Ψ₁(3)> - 2<Ψ(5),Ψ₁(5)>
    //    = <Ψ(2) + Ψ(3), Ψ₁(2) + Ψ₁(3)> - <2Ψ(5),2Ψ₁(5)>
    //    = <Ψ(2) + Ψ(3) - 2Ψ(5), Ψ₁(2) + Ψ₁(3) - 2Ψ₁(5)>
    // Ψ(2) + Ψ(3) - 2Ψ(5) = -γ + 1 + -γ + 1 + 1/2 - 2(-γ + 1 + 1/2 + 1/3 + 1/4)
    //    = -2γ + 2 + 1/2 - 2γ - 2 - 1 - 2/3 - 1/2 = -1 - 2/3
    //    = -5/3
    // Ψ₁(2) + Ψ₁(3) - 2Ψ₁(5) = 𝜋²/6 - 1 + 𝜋²/6 - (1 + 1/4) - 2[𝜋²/6 - (1 + 1/4 + 1/9 + 1/16)]
    //    = 𝜋²/3 - 1 - 1 - 1/4 - 2(𝜋²/6 - 1 - 1/4 - 1/9 - 1/16)
    //    = 𝜋²/3 - 2 - 1/4 - 𝜋²/3 + 2 + 1/2 + 2/9 + 1/8 = 1/2 - 1/4 + 1/8 + 2/9 = 1/4 + 25/72
    //    = 43/72
    // <f',f"> = <1/12,(-5/3)/12><-5/3,43/72> = <1/12,-5/36><-5/3,43,72> = <-5/36,25/108 + 43/864>
    //    = <-5/36,9/32>

    const w_act = beta(GDN!1(-0.5), GDN!1(1));
    // f = B(-0.5, 1) = Γ(-0.5)Γ(1)/Γ(0.5) = 𝜋/[sin(-𝜋/2)Γ(1.5)]/Γ(0.5) = -2𝜋Γ(1)/[√𝜋Γ(2)√𝜋]
    //   = -2
    // f' = B(-0.5,1)[Ψ(-0.5) + Ψ(1) - 2Ψ(0.5)]
    // Ψ(1.5) - Ψ(-0.5) = 𝜋⋅cot(-𝜋/2) =>
    // Ψ(-0.5) = Ψ(1.5) = -γ - 2ln(2) + ∑ᵢ₌₁¹2/(2i - 1)
    //         = 2 - γ - 2ln(2)
    // Ψ(0.5) = -γ - 2ln(2)
    // Ψ(1) = -γ
    // f' = -2{2 - γ - 2ln(2) + -γ - 2[-γ - 2ln(2)]} = -2[2 - 2γ - 2ln(2) + 2γ + 4ln(2)]
    //    = -4 - 4ln(2)
    const w_exp = GDN!1(-2, -4-4*LN2);
    assert(
        isClose(w_act.val, w_exp.val) && isClose(w_act.d, w_exp.d),
        format("B(-.5,1) = %s != %s", w_act, w_exp));

    const e = beta(GDN!1(-1), GDN!1(2));
    assert(isNaN(e.val) && isNaN(e.d));
    // f = B(-1, 2) = Γ(-1)Γ(1)/Γ(1) = Γ(-1), DNE

    const r = beta(GDN!1(-0.5), GDN!1(0.5));
    assert(r.val is -0. && isNaN(r.d));
    // f = B(-.5, .5) = Γ(-.5)Γ(.5)/Γ(0) =
    //   = -0
    // f' = B(-.5, .5)[Ψ(-.5) + Ψ(.5) - 2Ψ(0)]

    const t = beta(GDN!1(-1), GDN!1(1));
    assert(isNaN(t.val) && isNaN(t.d), format("B(-1,1) != %s", t));
    // lim(x->1) B(-x,x) = lim(x->1) Γ(-x)Γ(x)/Γ(0) = 0
    // lim(x->1) B(-x,1) = lim(x->1) Γ(-x)Γ(1)/Γ(1 - x) = lim(x->1) Γ(-x)/Γ(1 - x)
    // lim(x->1) 𝜋/[sin(-𝜋x)Γ(1+x)]/{𝜋/[sin(𝜋x)Γ(x)]} = lim(x->1) 𝜋sin(𝜋x)Γ(x)/[𝜋sin(-𝜋x)Γ(1+x)]
    // lim(x->1) sin(𝜋x)Γ(x)/[-sin(𝜋x)Γ(1+x)] = lim(x->1) -Γ(x)/Γ(1+x) = -Γ(1)/Γ(2) = -1
    // B(-1,1) DNE

    assert(beta(GDN!1(2), GDN!1(-1.5)) is beta(GDN!1(-1.5), GDN!1(2)));

    const y = std.mathspecial.beta(-0.5, -0.2);
    const u = std.mathspecial.digamma(-0.5);
    const i = std.mathspecial.digamma(-0.2);
    const o = std.mathspecial.digamma(-0.7);
    assert(beta(GDN!1(-0.5), GDN!1(-0.2)) is GDN!1(y, y*(u + i - 2*o)));
    // f = B(-.5,-.2)
    // f' = B(-.5,-.2)[Ψ(-.5) + Ψ(-.2) - 2Ψ(-.7)]

// NB: Fails because std.mathspecial.beta(-1.5, -.5) is NaN. Fixed in stable branch.
    // const p = beta(GDN!1(-1.5), GDN!1(-0.5));
    // // f = B(-1.5,-.5) = Γ(-1.5)Γ(-.5)/Γ(-2) = -0
    // // f' = f[Ψ(-1.5) + Ψ(-.5) - 2Ψ(-2)] DNE
    // assert(p.val is -0.0L && isNaN(p.d), format("B(-1.5,-0.5) != %s", p));

    const a = beta(GDN!1(-1), GDN!1(-0.5));
    assert(isNaN(a.val) && isNaN(a.d));
    // f = B(-1, -0.5) = Γ(-1)Γ(-.5)/Γ(-1.5) DNE

    const s = beta(GDN!1(-1), GDN!1(-2));
    assert(isNaN(s.val) && isNaN(s.d));

// NB: Fails because std.mathspecial.digamma(-0.) is NaN. Fixed in stable branch.
//     const d = beta(GDN!1(-0.), GDN!1(-0.5));
//     // f' = B(-0, -0.5)[Ψ(-0) + Ψ(-0.5) - 2Ψ(-0.5)] = -∞[∞ - Ψ(-0.5)] = -∞,
//     assert(d == -real.infinity && d.d == -real.infinity, format("B(-0,-.5) != %s", d));

    const f = beta(GDN!1(-0.), GDN!1(-1));
    assert(isNaN(f.val) && isNaN(f.d));

// NB: Fails because std.mathspecial.digamma(+0.) is NaN. Fixed in stable branch.
//     const g = beta(GDN!1(+0.), GDN!1(-1.5));
//     // f' = B(+0, -1.5)[Ψ(+0) + Ψ(-1.5) - 2Ψ(-1.5)] = ∞[-∞ - Ψ(-0.5)] = -∞,
//     assert(g == real.infinity && g.d == -real.infinity, format("B(+0, -1.5) = %s", g));

    const h = beta(GDN!1(+0.), GDN!1(-2));
    assert(isNaN(h.val) && isNaN(h.d));

// NB: Fails because std.mathspecial.digamma(+0.) = real.nan. Fixed in stable.
//    const j = beta(GDN!1(+0.), GDN!1(+0.));
//    assert(j == real.infinity && isNaN(j.d), format("B(+0,+0) = %s", j));

// NB: Fails because std.mathspecial.digamma(-0.) is NaN. Fixed in stable branch.
//     const k = beta(GDN!1(-0.), GDN!1(1));
//     // f' = B(-0, 1)[Ψ(-0) + Ψ(1) - 2Ψ(1)] = -∞[∞ - Ψ(1)] = -∞,
//     assert(k == -real.infinity && k.d == -real.infinity, format("B(-0, 1) = %s", k));

// NB: Fails because std.mathspecial.digamma(-0.) is NaN. Fixed in stable branch.
//     const l = beta(GDN!1(+0.), GDN!1(1));
//     // f' = B(+0, 1)[Ψ(+0) + Ψ(1) - 2Ψ(1)] = ∞[-∞ - Ψ(1)] = -∞,
//     assert(l == real.infinity && l.d == -real.infinity, format("B(+0,1) = %s", l));

// NB: Fails because std.mathspecial.beta(real.infinity, 1) is NaN. Fixed in
// stable branch.
//     const z = beta(GDN!1(real.infinity), GDN!1(1));
//     assert(z == 0 && isNaN(z.d), format("B(∞,1) != %s", z));

// NB: Fails because std.mathspecial.beta(real.infinity, real.infinity) is NaN.
// Fixed in stable branch.
//     const x = beta(GDN!1(real.infinity), GDN!1(real.infinity));
//     assert(x == 0 && isNaN(x.d), format("B(∞,∞) = %s", x));

    const c = beta(GDN!1(-real.infinity), GDN!1(1));
    assert(isNaN(c.val) && isNaN(c.d));

    const v = beta(GDN!1.nan, GDN!1(1));
    assert(isNaN(v.val) && isNaN(v.d));
}


/**
 * the digamma function,$(MATH Ψ), of a generalized dual number
 *
 * If $(MATH f(x) = Ψ(g(x))), then $(MATH f' = Ψ₁(g)g'), where $(MATH Ψ₁) is the polygamma function
 * of order $(MATH 1) (trigamma function).
 *
 * Params:
 *   Deg = the degree of g
 *   g = the `GDN`argument
 *
 * Returns:
 *   $(MATH Ψ(g)) as a `GDN`
 */
pure nothrow @nogc @safe GDN!Deg digamma(ulong Deg)(in GDN!Deg g)
{
    static if (Deg == 1)
        const df = ad.math.polygamma.polygamma!1(g.reduce());
    else
        const df = polygamma!1(g.reduce());

    return GDN!Deg(std.mathspecial.digamma(g.val), df*g.d);
}

///
unittest
{
    import ad.math.internal: isNaN;

    const trigamma_1 = PI^^2 / 6;  // Ψ₁(1)

    const f_act = digamma(GDN!1(1));
    const f_exp = GDN!1(std.mathspecial.digamma(1), trigamma_1);
    assert(f_act == f_exp);

    assert(isNaN(digamma(GDN!1(-1))));
    assert(digamma(GDN!1(real.infinity)) is GDN!1(real.infinity, +0.));
}

unittest
{
    import std.format: format;

    const e = GDN!1(std.mathspecial.digamma(2), 3*ad.math.polygamma.polygamma!1(2));
    assert(digamma(GDN!1(2, 3)) is e);

// NB: Fails because of https://github.com/dlang/phobos/issues/10802, fixed on stable
//     const f_nz = digamma(GDN!1(-0.));
//     assert(f_nz is GDN!1(real.infinity, real.infinity), format("Ψ₁(-0) = %s", f_nz));

// NB: Fails because of https://github.com/dlang/phobos/issues/10802, fixed on stable
//     const f_pz = digamma(GDN!1(+0.));
//     assert(f_pz is GDN!1(-real.infinity, real.infinity), format("Ψ₁(+0) = %s", f_pz));

    const q = digamma(GDN!1(-real.infinity));
    assert(isNaN(q.val) && isNaN(q.d));

    const g = digamma(GDN!1(real.nan));
    assert(isNaN(g.val) && isNaN(g.d));

    const w = GDN!2(
        std.mathspecial.digamma(2),
        ad.math.polygamma.polygamma!1(2),
        ad.math.polygamma.polygamma!2(2));

    assert(digamma(GDN!2(2)) is w);
}


/**
 * The natural logarithm of a `GDN` minus digamma of the same `GDN`.
 *
 * If $(MATH f(x) = ln(g(x)) - Ψ(g(x))), then $(MATH f' = [1/g - Ψ₁(g)]g'), where $(MATH Ψ₁) is the
 * polygamma function of order one (trigamma function).
 *
 * Params:
 *   Deg = the degree of g
 *   g = the argument, must be positive
 *
 * Returns:
 *   a GDN representing the natural logarithm of g minus digamma of g.
 */
pure nothrow @nogc @safe GDN!Deg logmdigamma(ulong Deg)(in GDN!Deg g)
in(signbit(g) == 0, "the argument must be positive")
do {
    const f = std.mathspecial.logmdigamma(g.val);
    const g_red = g.reduce();

    static if (Deg == 1)
        const dfdg = 1/g_red - ad.math.polygamma.polygamma!1(g_red);
    else
        const dfdg = g_red.inv() - polygamma!1(g_red);

    return GDN!Deg(f, dfdg*g.d);
}

///
unittest
{
    import std.math: isClose;

    // Euler-Mascheroni constant
    const γ = 0.577_215_664_901_532_860_607L;

    const f_act = logmdigamma(GDN!1(1));
    const f_exp = GDN!1(γ, 1-PI^^2/6);
    assert(isClose(f_act.val, f_exp.val) && isClose(f_act.d, f_exp.d));
}

unittest
{
    import std.format: format;
    import std.math: isClose;
    import ad.math: log;

    const w = logmdigamma(GDN!1(+0.));
    assert(w == real.infinity && isNaN(w.d), format("logmdigamma(+0) != %s", w));
    // ln(x) - 1/x ≤ Ψ(x) ≤ ln(x) - 1/(2x), x>0
    // -1/x ≤ Ψ(x) - ln(x) ≤ -1/(2x)
    // 1/x ≥ ln(x) - Ψ(x) ≥ 1/(2x)
    // lim{x⟶0⁺} ln(x) - Ψ(x) ≥ lim{x⟶0⁺} 1/(2x) = +∞

    const e_act = logmdigamma(GDN!2(2));
    const e_exp = log(GDN!2(2)) - digamma(GDN!2(2));
    assert(isClose(e_act.val, e_exp.val));
    assert(isClose(e_act.d.val, e_exp.d.val, 10*real.epsilon));
    assert(e_act.d!2 == e_exp.d!2);
}


/**
 * The inverse of the function formed from the natural logarithm of a `GDN` minus digamma of the
 * same `GDN`.
 *
 * If $(MATH f(x) = ln(g(x)) - Ψ(g(x))), then $(MATH g = f⁻¹(f)) and $(MATH f' = [1/g - Ψ₁(g)]g'),
 * where $(MATH Ψ₁) is the polygamma function of order one (trigamma function). Thus
 * $(MATH g' = f' / [1/g - Ψ₁(g)] = f'g / [1 - gΨ₁(g)])
 *
 * Params:
 *   Deg = the degree of f
 *   f = the `GDN` argument
 *
 * Return:
 *   It returns `g` such that the natural logarithm of `g` minus digamma of `g` is equal to f.
 */
pure nothrow @nogc @safe GDN!Deg logmdigammaInverse(ulong Deg)(in GDN!Deg f)
{
    static if (Deg == 1) {
        alias ln_m_digamma_inv = std.mathspecial.logmdigammaInverse;
        alias trigamma = ad.math.polygamma.polygamma!1;
    } else {
        alias ln_m_digamma_inv =  logmdigammaInverse;
        alias trigamma = polygamma!(1, Deg-1);
    }

    const g_red = ln_m_digamma_inv(f.reduce());

    GDN!Deg.DerivType!1 naive_derivative() { return f.d * g_red / (1 - g_red*trigamma(g_red)); }

    // Assume x > 0. (x + 1/2)/x² ≤ Ψ₁(x) ≤ (x + 1)/x² ⇒ -1/x ≤ 1 - xΨ₁(x) ≤ -1/(2x).
    // lim{x⟶0⁺} -1/x = -∞ = lim{x⟶0⁺} -1/(2x) ⇒ lim{x⟶0⁺} 1 - xΨ₁(x) = -∞
    // lim{x⟶+∞} -1/x = 0⁻ = lim{x⟶+∞} -1/(2x) ⇒ lim{x⟶+∞} 1 - xΨ₁(x) = 0⁻
    GDN!Deg.DerivType!1 dg;

    if (g_red == 0 && signbit(g_red) == 0) {
        static if (Deg == 1)
            dg = -0. * f.d;
        else
            dg = GDN!Deg.DerivType!1(asReal(g_red.d), naive_derivative().d);
    } else if (g_red == real.infinity) {
        static if (Deg == 1)
            dg = -real.infinity * f.d;
        else
            dg = GDN!Deg.DerivType!1(asReal(g_red.d), naive_derivative().d);
    } else {
        dg = naive_derivative();
    }

    return GDN!Deg(asReal(g_red), dg);
}

///
unittest
{
    import std.math: isClose;

    const f = GDN!1(1);
    const g = logmdigammaInverse(logmdigamma(f));
    assert(isClose(g.val, f.val) && isClose(g.d, f.d));
}

unittest
{
    import std.format: format;
    import std.math: isClose;

    const γ = 0.577_215_664_901_532_860_607L;
    const ζ3 = 1.202_056_903_159_594_285_400L;

    const q = logmdigammaInverse(GDN!1(+0.0L));
    // g = +∞
    // g' = lim{g⟶+∞} 1/[1/g - Ψ₁(g)] = lim{g⟶0⁺} 1/g - Ψ₁(g) = -lim{g⟶0⁺} Ψ₁(g) - 1/g
    //    ≤ -lim{g⟶0⁺} (g + 1)/g² - 1/g = -lim{g⟶0⁺} 1/g² = -∞
    // g' = -∞
    assert(q is GDN!1(real.infinity, -real.infinity), format("logmdigammaInverse(+0) != %s", q));

    const w = logmdigammaInverse(GDN!1(real.infinity));
    // lim{g⟶0⁺} Ψ(g) ~ lim{g⟶0⁺} ln(g) - 1/(2g) ⇒ lim{g⟶0⁺} Ψ(g) - ln(q) ~ lim{g⟶0⁺} -1/(2g)
    // ⇒ lim{g⟶0⁺} ln(g) - Ψ(g) ~ lim{g⟶0⁺} 1/(2g) = +∞
    // ⇒ g = 0⁺
    // g' = lim{g⟶0⁺} 1/[1/g - Ψ₁(g)] = lim{g⟶+∞} 1/g - Ψ₁(g) = -lim{g⟶+∞} Ψ₁(g) - 1/g
    // -lim{g⟶+∞} (g + 1/2)/g² - 1/g ≤ -lim{g⟶+∞} Ψ₁(g) - 1/g ≤ -lim{g⟶+∞} (g + 1)/g² - 1/g
    // ⇒ lim{g⟶+∞} 1/g - (g + 1/2)/g² ≤ lim{g⟶+∞} 1/g - Ψ₁(g) ≤ lim{g⟶+∞} 1/g - (g + 1)/g²
    // ⇒ lim{g⟶+∞} -1/(2g²) ≤ lim{g⟶+∞} 1/g - Ψ₁(g) ≤ lim{g⟶+∞} -1/g²
    // ⇒ 0⁻ ≤ lim{g⟶+∞} 1/g - Ψ₁(g) ≤ 0⁻ ⇒ lim{g⟶+∞} 1/g - Ψ₁(g) = 0⁻
    // g' = 0⁻
    assert(w is GDN!1(+0., -0.), format("logmdigammaInverse(∞) != %s", w));

// NB: Failed because std.mathspecial.logmdigammaInverse(-.0) is -real.infinity, fixed in stable
//     const e = logmdigammaInverse(GDN!1(-0.));
//     // Assume x > 0. ln(x) - 1/(2x) - Ψ₁(x) > 0 ⇒ ln(x) - Ψ₁(x) > 1/(2x) > 0.
//     // f(x) = ln(x) - Ψ₁(x) is defined for x ∊ ℝ⁺. f: ℝ⁺ ↦ ℝ⁺ ⇒ f⁻¹: ℝ⁺ ↦ ℝ⁺
//     // f(-0) DNE,
//     assert(isNaN(e), format("logmdigammaInverse(-0) = %s", e));

    const r = logmdigammaInverse(GDN!1(-1));
    assert(isNaN(r.val) && isNaN(r.d));

    const t = GDN!2(γ);
    const u = logmdigammaInverse(t);
    // ln(1) - Ψ(1) == γ ⇒ g = 1
    // g' = f'g/[1 - gΨ₁(g)] = 1*1/[1 - 1*Ψ₁(1)] = 1/(1 - 𝜋²/6) = 6/(6 - 𝜋²)
    // <g',g"> = <f',f"><g,g'>/[1 - <g,g'>Ψ₁(<g,g'>)]
    //         = <1,0><1,g'>/[1 - <1,g'>Ψ₁(<1,g'>)]
    //         = <1,g'>/[1 - <1,g'><Ψ₁(1),Ψ₂(1)g'>]
    // Ψ₂(1) = (-1)³2!ζ(3,1) = -2ζ(3)
    // <g',g"> = <1,g'>/[1 - <1,g'><𝜋²/6,-2ζ(3)g'>]
    //         = <1,g'>/[1 - <𝜋²/6,𝜋²g'/6 - 2ζ(3)g'>]
    //         = <1,g'>/<1 - 𝜋²/6,[2ζ(3) - 𝜋²/6]g'>
    //         = <1/(1 - 𝜋²/6),{(1 - 𝜋²/6)g' - [2ζ(3) - 𝜋²/6]g'}/(1 - 𝜋²/6)²>
    //         = <g',[1 - 2ζ(3)](g')³>
    // g" = [1 - 2ζ(3)](g')³ = [1 - 2ζ(3)][6/(6 - 𝜋²)]³
    // g" = 216[1 - 2ζ(3)]/(6 - 𝜋²)³
    const dg = 6/(6 - PI^^2);
    const d2g = 216*(1 - 2*ζ3)/(6 - PI^^2)^^3;
    assert(
        isClose(u.val, 1) && isClose(u.d.val, dg) && isClose(u.d!2, d2g),
        format("logmdigammaInverse(γ) != %s", u));
}


/* This function computes the derivative of I(g(x); a,b) with respect to x,
 * where a and b are constants.
 *
 * Let f(x) = I(g(x); a,b). Then f' = dI/dg⋅g' = g'(d/dg)B(g; a,b)/B(a,b)) where
 * B(g; a,b) = ∫₀ᵍtᵃ⁻¹(1-t)ᵇ⁻¹dt is the incomplete beta function.
 *
 * The integrand tᵃ⁻¹(1-t)ᵇ⁻¹ is Lebesgue integrable over 0 ≤ t ≤ 1. Therefore,
 * B' = gᵃ⁻¹(1-g)ᵇ⁻¹ almost everywhere. g = 0 and 1 are the only values where
 * gᵃ⁻¹(1-g)ᵇ⁻¹ doesn't exist for every a and b, but the one-sided limits do.
 * This algorithm defines
 * B'(0; a,b) = lim{g→0⁺} gᵃ⁻¹(1-g)ᵇ⁻¹ = { ∞, 0<a<1; 1, a=1; 0, a>1 }), and
 * B'(1; a,b) = lim{g→1⁻} gᵃ⁻¹(1-g)ᵇ⁻¹ = { ∞, 0<b<1; 1, b=1; 0, b>1 }).
 *
 * Thus I'(g; a,b) has the following form when 0 < a,b < ∞.
 *
 *    - I'(g; a,b) = gᵃ⁻¹(1-g)ᵇ⁻¹/B(a,b), 0 < g < 1
 *    - I'(0; a,b) = { ∞, a<1; 1/B(1,b), a=1; 0, a>1 }
 *    - I'(1; a,b) = { ∞, b<1; 1/B(a,1), b=1; 0, b>1 }
 *
 * Here are the degenerate cases of I'. Let H(x) = { 0, x<0; 1, x≥0 } be the
 * Heaviside step function in the following.
 *
 *    - I'(g; 0,b) = (d/dg)lim{a→0⁺} I(g; a,b) = (d/dg)(1 - H(-g)) = 𝛿(g)
 *    - I'(g; ∞,b) = (d/dg)lim{a→∞} I(g; a,b) = (d/dg)H(g-1) = 𝛿(g-1)
 *    - I'(g; a,0) = (d/dg)lim{b→0⁺} I(g; a,b) = (d/dg)H(g-1) = 𝛿(g-1)
 *    - I'(g; a,∞) = (d/dg)lim{b→∞} I(g; a,b) = (d/dg)[1 - H(-g)] = 𝛿(g)
 *    - I'(g; 0,0) = (d/dg)lim{a,b→0⁺} I(g; a,b), does not exist
 *    - I'(g; 0,∞) = (d/dg)lim{a→0⁺,b→∞) I(g; a,b) = (d/dg)[1 - H(-g)] = 𝛿(g)
 *    - I'(g; ∞,0) = (d/dg)lim{a→∞,b→0⁺} I(g; a,b) = (d/dg)H(g-1) = 𝛿(g-1)
 *    - I'(g; ∞,∞) = (d/dg)lim{a,b→∞) I(g; a,b), does not exist
 */
private pure nothrow @nogc @safe
GDN!Deg.DerivType!1 betaIncompleteDeriv(ulong Deg)(in real a, in real b, in GDN!Deg g)
{
    alias Deriv = typeof(return);

    Deriv dfdg;
    if ((a == 0 && b == 0) || (a == real.infinity && b == real.infinity)) {
        dfdg = Deriv.nan;
    } else if (a == 0 || b == real.infinity) {
        dfdg = dirac(g.reduce());
    } else if (a == real.infinity || b == 0) {
        dfdg = dirac(g.reduce() - 1);
    } else {
        const g_red = g.reduce();
        const numerator = g_red^^(a - 1) * (1 - g_red)^^(b - 1);
        const denominator = std.mathspecial.beta(a, b);
        if (numerator == real.infinity && denominator == real.infinity) {
            // In this case, the denominator is not really infinite. It's just
            // larger than real.max. Instead of returning ∞/∞ = NaN, return the
            // numerator. The derivatives still need to be scaled.
            static if (Deg == 1)
                dfdg = numerator;
            else
                dfdg = Deriv(numerator.val, numerator.d/denominator);
        } else if (numerator == 0 && denominator == 0) {
            // In this case, the denominator is not really zero. It's just too
            // small to represent. Instead of return 0/0 = NaN, return the
            // numerator. The derivatives still need to be scaled.
            static if (Deg == 1)
                dfdg = numerator;
            else
                dfdg = Deriv(numerator.val, numerator.d/denominator);
        } else {
            dfdg = numerator / denominator;
        }
    }

    return dfdg * g.d;
}

unittest
{
    import std.format: format;

    //
    // Tests of degree 1
    //

    // a = +0, b = +0

    assert(isNaN(betaIncompleteDeriv(+0., +0., GDN!1(0))));
    assert(isNaN(betaIncompleteDeriv(+0., +0., GDN!1(.5))));
    assert(isNaN(betaIncompleteDeriv(+0., +0., GDN!1(1))));

    // a = +0, b = .5

    assert(betaIncompleteDeriv(+0., .5, GDN!1(0)) is real.infinity);
    assert(betaIncompleteDeriv(+0., .5, GDN!1(.5)) == 0);
    assert(betaIncompleteDeriv(+0., .5, GDN!1(1)) == 0);

    // a = +0, b = 1

    assert(betaIncompleteDeriv(+0., 1, GDN!1(0)) is real.infinity);
    assert(betaIncompleteDeriv(+0., 1, GDN!1(.5)) == 0);
    assert(betaIncompleteDeriv(+0., 1, GDN!1(1)) == 0);

    // a = +0, b = 2

    assert(betaIncompleteDeriv(+0., 2, GDN!1(0)) is real.infinity);
    assert(betaIncompleteDeriv(+0., 2, GDN!1(.5)) == 0);
    assert(betaIncompleteDeriv(+0., 2, GDN!1(1)) == 0);

    // a = +0, b = ∞

    assert(betaIncompleteDeriv(+0., real.infinity, GDN!1(0)) is real.infinity);
    assert(betaIncompleteDeriv(+0., real.infinity, GDN!1(.5)) == 0);
    assert(betaIncompleteDeriv(+0., real.infinity, GDN!1(1)) == 0);

    // a = .5, b = +0

    assert(betaIncompleteDeriv(.5, +0., GDN!1(0)) == 0);
    assert(betaIncompleteDeriv(.5, +0., GDN!1(.5)) == 0);
    assert(betaIncompleteDeriv(.5, +0., GDN!1(1)) is real.infinity);

    // a = .5, b = .5

    assert(betaIncompleteDeriv(.5, .5, GDN!1(0)) is real.infinity);

    assert(isClose(betaIncompleteDeriv(.5, .5, GDN!1(.5)), M_2_PI));
    // f' = .5^-.5*.5^-.5/B(.5,.5) = 2Γ(1)/[Γ(.5)Γ(.5)] = 2/(√𝜋√𝜋) = 2/𝜋

    assert(betaIncompleteDeriv(.5, .5, GDN!1(1)) is real.infinity);

    // a = .5, b = 1

    assert(betaIncompleteDeriv(.5, 1, GDN!1(0)) is real.infinity);

    assert(betaIncompleteDeriv(.5, 1, GDN!1(.5)) == SQRT1_2);
    // f' = .5^-.5/B(.5,1) = √2*Γ(1.5)/[Γ(.5)Γ(1)] = √2(√𝜋/2)/√𝜋 = √2/2

    assert(betaIncompleteDeriv(.5, 1, GDN!1(1)) == .5);
    // f' = 1/B(.5,1) = Γ(1.5)/[Γ(.5)Γ(1)] = (√𝜋/2)/√𝜋 = 1/2

    // a = .5, b = 2

    assert(betaIncompleteDeriv(.5, 2, GDN!1(0)) is real.infinity);

    assert(betaIncompleteDeriv(.5, 2, GDN!1(.5)) == 3*SQRT2/8);
    // f' = .5^-.5*.5/B(.5,2) = (√2/2)Γ(5/2)/[Γ(1/2)Γ(2)] = √2(3√𝜋/4)/(2√𝜋) = 3√2/8

    assert(betaIncompleteDeriv(.5, 2, GDN!1(1)) == 0);

    // a = .5, b = ∞

    assert(betaIncompleteDeriv(.5, real.infinity, GDN!1(0)) is real.infinity);
    assert(betaIncompleteDeriv(.5, real.infinity, GDN!1(.5)) == 0);
        assert(betaIncompleteDeriv(.5, real.infinity, GDN!1(1)) == 0);

    // a = 1, b = +0

    assert(betaIncompleteDeriv(1, +0., GDN!1(0)) == 0);
    assert(betaIncompleteDeriv(1, +0., GDN!1(.5)) == 0);
    assert(betaIncompleteDeriv(1, +0., GDN!1(1)) is real.infinity);

    // a = 1, b = .5

    assert(betaIncompleteDeriv(1, .5, GDN!1(0)) == 0.5L);

    const zhh = betaIncompleteDeriv(1, .5, GDN!1(.5));
    // f' = .5^-.5/B(1,.5) = √2/2
    assert(zhh == SQRT1_2, format("I'(.5; 1,.5) = %s", zhh));

    assert(betaIncompleteDeriv(1, .5, GDN!1(1)) is real.infinity);

    // a = 1, b = 1

    assert(betaIncompleteDeriv(1, 1, GDN!1(0)) == 1);

    const ooh = betaIncompleteDeriv(1, 1, GDN!1(0.5L));
    assert(ooh == 1, format("I'(.5; 1,1) = %s", ooh));

    assert(betaIncompleteDeriv(1, 1, GDN!1(1)) == 1);

    // a = 1, b = 2

    assert(betaIncompleteDeriv(1, 2, GDN!1(0)) == 2);
    // f' = 1/B(1,2) = Γ(3)/[Γ(1)Γ(2)] = 2

    assert(betaIncompleteDeriv(1, 2, GDN!1(.5)) == 1);
    // f' = .5/B(1,2) = Γ(3)/[2Γ(1)Γ(2)] = 2/2 = 1

    assert(betaIncompleteDeriv(1, 2, GDN!1(1)) == 0);

    // a = 1, b = ∞

    assert(betaIncompleteDeriv(1, real.infinity, GDN!1(0)) is real.infinity);
    assert(betaIncompleteDeriv(1, real.infinity, GDN!1(.5)) == 0);
    assert(betaIncompleteDeriv(1, real.infinity, GDN!1(1)) == 0);

    // a = 2, b = +0

    assert(betaIncompleteDeriv(2, +0., GDN!1(0)) == 0);
    assert(betaIncompleteDeriv(2, +0., GDN!1(.5)) == 0);
    assert(betaIncompleteDeriv(2, +0., GDN!1(1)) is real.infinity);

    // a = 2, b = .5

    assert(betaIncompleteDeriv(2, .5, GDN!1(0)) == 0);

    assert(betaIncompleteDeriv(2, .5, GDN!1(.5)) == 3*SQRT2/8);
    // f' = .5*.5^-.5/B(2,.5) = (√2/2)Γ(5/2)/[Γ(1/2)Γ(2)] = 3√2/8

    assert(betaIncompleteDeriv(2, .5, GDN!1(1)) is real.infinity);

    // a = 2, b = 1

    assert(betaIncompleteDeriv(2, 1, GDN!1(0)) == 0);

    assert(betaIncompleteDeriv(2, 1, GDN!1(0.5)) == 1);
    // f' = .5^2/B(2,1) = 1

    assert(betaIncompleteDeriv(2, 1, GDN!1(1)) == 2);

    // a = 2, b = 2

    assert(betaIncompleteDeriv(2, 2, GDN!1(0)) == 0);

    assert(betaIncompleteDeriv(2, 2, GDN!1(0.5)) == 1.5L);
    // f' = .5*.5/B(2,2) = Γ(4)/[4Γ(2)Γ(2)] = 6/4 = 3/2

    assert(betaIncompleteDeriv(2, 2, GDN!1(1)) == 0);

    // a = 2, b = ∞

    assert(betaIncompleteDeriv(2, real.infinity, GDN!1(0)) is real.infinity);
    assert(betaIncompleteDeriv(2, real.infinity, GDN!1(.5)) == 0);
    assert(betaIncompleteDeriv(2, real.infinity, GDN!1(1)) == 0);

    // a = ∞, b = +0

    assert(betaIncompleteDeriv(real.infinity, +0., GDN!1(0)) == 0);
    assert(betaIncompleteDeriv(real.infinity, +0., GDN!1(.5)) == 0);
    assert(betaIncompleteDeriv(real.infinity, +0., GDN!1(1)) is real.infinity);

    // a = ∞, b = .5

    assert(betaIncompleteDeriv(real.infinity, .5, GDN!1(0)) == 0);
    assert(betaIncompleteDeriv(real.infinity, .5, GDN!1(.5)) == 0);
    assert(betaIncompleteDeriv(real.infinity, .5, GDN!1(1)) is real.infinity);

    // a = ∞, b = 1

    assert(betaIncompleteDeriv(real.infinity, 1, GDN!1(0)) == 0);
    assert(betaIncompleteDeriv(real.infinity, 1, GDN!1(.5)) == 0);
    assert(betaIncompleteDeriv(real.infinity, 1, GDN!1(1)) is real.infinity);

    // a = ∞, b = 2

    assert(betaIncompleteDeriv(real.infinity, 2, GDN!1(0)) == 0);
    assert(betaIncompleteDeriv(real.infinity, 2, GDN!1(.5)) == 0);
    assert(betaIncompleteDeriv(real.infinity, 2, GDN!1(1)) is real.infinity);

    // a = ∞, b = ∞

    assert(isNaN(betaIncompleteDeriv(real.infinity, real.infinity, GDN!1(0))));
    assert(isNaN(betaIncompleteDeriv(real.infinity, real.infinity, GDN!1(.5))));
    assert(isNaN(betaIncompleteDeriv(real.infinity, real.infinity, GDN!1(1))));

    // g' ≠ 1

    assert(betaIncompleteDeriv(2, 2, GDN!1(0.5, 0.5)) == 0.75L);
    // from earlier test: I'(GDN!1(0.5); 2,2) = 1.5 when g' = 1 so when g' = 0.5
    // expect I' = 1.5 * 0.5 = 0.75

    //
    // Tests of degree 2
    //

    // a = 0⁺

    assert(betaIncompleteDeriv(+0., 1, GDN!2(.5)) is GDN!1(0, 0));

    // a = 1, b = 0⁺

    assert(betaIncompleteDeriv(1, +0., GDN!2(.5)) is GDN!1(0, 0));
    // limit I(.5;1,0+) = 0 and derivative is dirac(g-1)==0

    // a = 1, b = 1

    assert(betaIncompleteDeriv(1, 1, GDN!2(0.0L)) is GDN!1(1, 0));

    assert(betaIncompleteDeriv(2, 2, GDN!2(.5)) is GDN!1(1.5, 0));
    // f' = <.5,1>^1(1 - <.5,1>)^1/B(2,2) = <.5,1><.5,-1>Γ(4)/[Γ(2)Γ(2)] = 6<.25,-.5 + .5> =<1.5,0>

    assert(betaIncompleteDeriv(1, 1, GDN!2(1.0L)) is GDN!1(1.0L, 0));

    // a = 1, b = ∞

    assert(betaIncompleteDeriv(1, real.infinity, GDN!2(.5)) is GDN!1(0, 0));

    // a = ∞, b = 1

    assert(betaIncompleteDeriv(real.infinity, 1, GDN!2(.5)) is GDN!1(0, 0));

    // g' ≠ 1

    assert(betaIncompleteDeriv(1, 1, GDN!2(.5, 2.0L, 0.0L)) is GDN!1(2.0L, 0.0L));

    // g" ≠ 0

    assert(betaIncompleteDeriv(1, 1, GDN!2(.5, 1.0L, 1.0L)) is GDN!1(1.0L, 1.0L));
}


/**
 * The regularized incomplete beta function $(MATH I$(SUB g)(a,b)).
 *
 * For fixed $(MATH a,b > 0), let $(MATH f(x) = I$(SUB g(x))(a,b)). Then
 * $(MATH f' = g'$(SUP dI$(SUB g))/$(SUB dg)).
 *
 * Params:
 *   Deg = the degree of g
 *   a = the first shape parameter, must be positive
 *   b = the second shape parameter, must be positive
 *   g = the argument, must belong to the interval $(MATH [0,1])
 *
 * Returns:
 *   the regularized incomplete beta function evaluated at g expressed as _a `GDN`
 */
pure nothrow @nogc @safe GDN!Deg betaIncomplete(ulong Deg)(in real a, in real b, in GDN!Deg g)
in {
    assert(isNaN(a) || signbit(a) == 0, "the first shape parameter must be positive");
    assert(isNaN(b) || signbit(b) == 0, "the second shape parameter must be positive");
    assert(isNaN!Deg(g) || (g >= 0 && g <= 1), "the argument must be in [0,1]");
}
do {
    return GDN!Deg(std.mathspecial.betaIncomplete(a, b, g.val), betaIncompleteDeriv(a, b, g));
}

///
unittest
{
    assert(betaIncomplete(1, 1, GDN!1(0.5)) is GDN!1(0.5, 1));
}


// NB: In master, but not released.
// /**
//  * The regularized incomplete beta complement function $(MATH I$(SUB g)$(SUP C)(a,b)).
//  *
//  * For fixed $(MATH a,b > 0), if $(MATH f(x) = I$(SUB g(x))$(SUP C)(a,b)), then
//  * $(MATH f' = g'$(SUP dI$(SUB g)$(SUP C))/$(SUB dg)). Since
//  * $(MATH I$(SUB g)$(SUP C) = 1 - I$(SUB g)), $(MATH I$(SUB g)$(SUP C)' = -I$(SUB g)'). Thus
//  * $(MATH f' = -g'$(SUP dI$(SUB g))$(SUB dg)).
//  *
//  * Params:
//  *   Deg = the degree of g
//  *   a = the first shape parameter, must be positive
//  *   b = the second shape parameter, must be positive
//  *   g = the argument, must belong to the interval $(MATH [0,1])
//  *
//  * Returns:
//  *   the regularized incomplete beta complement function evaluated at g expressed as _a `GDN`
//  */
// pure nothrow @nogc @safe GDN!Deg betaIncompleteCompl(ulong Deg)(in real a, in real b, in GDN!Deg g)
// in {
//     assert(isNaN(a) || signbit(a) == 0, "the first shape parameter must be positive");
//     assert(isNaN(b) || signbit(b) == 0, "the second shape parameter must be positive");
//     assert(isNaN!Deg(g) || (g >= 0 && g <= 1), "the argument must be in [0,1]");
// }
// do {
//     return GDN!Deg(std.mathspecial.betaIncompleteCompl(a, b, g.val), -betaIncompleteDeriv(a, b, g));
// }
//
// ///
// unittest
// {
//     const a = 1.0L;
//     const b = 2.0L;
//     const x = GDN!1(0.1L);
//     const i = betaIncomplete(a, b, x);
//     const icc = betaIncompleteCompl(b, a, 1.0L - x);
//     assert(isClose(icc, i) && icc.d == i.d);
// }


// TODO: implement the following functions.
// betaIncompleteInverse
// gammaIncomplete
// gammaIncompleteCompl
// gammaIncompleteComplInverse
// erf
// erfc
// normalDistribution
// normalDistributionInverse
