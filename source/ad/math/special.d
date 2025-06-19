module ad.math.special;

static import std.mathspecial;

import std.math: isNaN;
import std.traits: select;

static import ad.math.polygamma;

import ad.core;
import ad.math.internal: asReal, sgn;

/**
 * the gamma function, $MATH(Γ), of a generalized dual number
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
	import std.math: isNaN;

	assert(gamma(GDN!1(2, 3)) is GDN!1(1, 3*std.mathspecial.digamma(2)));

	const q = gamma(GDN!1(real.nan));
	assert(isNaN(q.val) && isNaN(q.d));

	const t = gamma(GDN!1(-0.));
	assert(t == GDN!1(-real.infinity) && isNaN(t.d), format("Γ(-0) != %s", t));

	const y = gamma(GDN!1(+0.));
	assert(y == GDN!1(real.infinity) && isNaN(y.d));

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
	//         = <2,2Ψ(3)><Ψ(3),Ψ⁽¹⁾(3)><1,0>
	//         = <2Ψ(3),2Ψ²(3)+2Ψ⁽¹⁾(3)>
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
	import std.math: isNaN, log;

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
	// <f',f"> = Ψ(<5,1>)*<1,0> = <Ψ(5), Ψ⁽¹⁾(5)>
}


/**
 * Computes the sign of the gamma function of a generalized dual number.
 *
 * If $(MATH f(x) = sgn(Γ(g(x))).
 *
 * Params:
 *   Deg = the degree of g
 *   g = the `GDN` argument
 *
 * Returns:
 *   It returns $(MATH  sgn(Γ(g))) as a `GDN`.
 */
pure nothrow @nogc @safe GDN!Deg sgnGamma(ulong Deg)(in GDN!Deg g)
out (f; f == std.mathspecial.sgnGamma(g.val), "result doesn't agree with std.mathspecial.sgnGamma")
{
	return sgn(gamma(g));
}

///
unittest
{
	const f = sgnGamma(GDN!1(1));
	assert(f == 1 && f.d == 0);

	// TODO: figure out why this fails
	import std.format: format;
	assert(std.mathspecial.sgnGamma(-0.5) == -1, format("sgnGamma(-0.5) = %s", std.mathspecial.sgnGamma(-0.5)));
	const g = sgnGamma(GDN!1(-0.5));
	assert(g == -1 && g.d == 0);

	// const h = sgnGamma(GDN!1(+0.));
	// assert(h == 1 && isNaN(h.d));

	const i = sgnGamma(GDN!1(-1));
	assert(isNaN(i.val) && isNaN(i.d));
}


// TODO: implement beta


/**
 * the digamma function,$(MATH Ψ), of a generalized dual number
 *
 * If $(MATH f(x) = Ψ(g(x))), then $(MATH f' = Ψ⁽¹⁾(g)g').
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
	import ad.math.traits: isNaN;

	assert(digamma(GDN!1(1)) is GDN!1(std.mathspecial.digamma(1), ad.math.polygamma.polygamma!1(1)));
	assert(isNaN(digamma(GDN!1(-1))));
	assert(digamma(GDN!1(real.infinity)) is GDN!1(real.infinity, +0.));
}

unittest
{
	import std.format: format;
	import std.math: isNaN;

	const e = GDN!1(std.mathspecial.digamma(2), 3*ad.math.polygamma.polygamma!1(2));
	assert(digamma(GDN!1(2, 3)) is e);

	const f_nz = digamma(GDN!1(-0.));
	assert(isNaN(f_nz.val) && isNaN(f_nz.d), format("Ψ⁽¹⁾(-0) = %s", f_nz));

	const f_pz = digamma(GDN!1(-0.));
	assert(isNaN(f_pz.val) && isNaN(f_pz.d), format("Ψ⁽¹⁾(+0) = %s", f_pz));

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


private pure nothrow @nogc @safe GDN!Deg polygamma(ulong n, ulong Deg)(in GDN!Deg g) if (n > 0)
{
	static if (Deg == 1)
		const df = ad.math.polygamma.polygamma!(n+1)(g.reduce());
	else
		const df = polygamma!(n+1)(g.reduce());

	return GDN!Deg(ad.math.polygamma.polygamma!n(g.val), df*g.d);
}

unittest
{
	import std.format: format;
	import std.math: isNaN;

	const q = GDN!1(ad.math.polygamma.polygamma!1(1), ad.math.polygamma.polygamma!2(1));
	assert(polygamma!1(GDN!1(1)) is q);

	const a_exp = GDN!1(ad.math.polygamma.polygamma!1(2), 3*ad.math.polygamma.polygamma!2(2));
	const a_act = polygamma!1(GDN!1(2, 3));
	assert(a_act is a_exp, format("Ψ⁽¹⁾(<2,3>) = %s != %s", a_act, a_exp));

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
	assert(isNaN(t.val) && isNaN(t.d), format("Ψ⁽¹⁾(-∞) != %s", t));

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
	// <f',f"> = <1,0>Ψ⁽²⁾(<3,1>) = <1,0><Ψ⁽²⁾(3),Ψ⁽³⁾(3)> = <Ψ⁽²⁾(3),Ψ⁽³⁾(3)>
	const o_act = polygamma!1(GDN!2(3));
	assert(o_act is o_exp, format("Ψ⁽¹⁾(3) = %s != %s", o_exp, o_act));
}
