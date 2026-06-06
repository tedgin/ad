/**
* This module extends the `core.math` Phobos module to support `GDN` objects.
*
* To emulate being implemented as intrinsics, all functions are declared to be inlined.
*/
module ad.math.core;

static import core.math;

import std.math : isNaN, LN2;
import std.meta : allSatisfy, anySatisfy;
import std.traits : isFloatingPoint, Select;

import ad;
import ad.internal : asGDN, CommonGDN, isConvertibleToGDN, isGDN, isNaN, signbit;


/**
 * This function computes the cosine of the argument.
 *
 * If $(MATH f(x) = cos(g(x))), then $(MATH f' = -sin(g)g').
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` to compute the cosine of
 *
 * Returns:
 *   the cosine of `g`
 */
pragma(inline, true) pure nothrow @nogc @safe GDN!Deg cos(ulong Deg)(in GDN!Deg g)
do {
	alias sine = Select!(Deg == 1, core.math.sin, sin);

	if (isNaN(g)) return g;
	return GDN!Deg(core.math.cos(g.val), -sine(g.reduce())*g.d);
}
/***/ unittest {
	const g = GDN!2(0);
	const f = cos(g);
	assert(f == 1 && f.d == 0 && f.d!2 == -1);
}
unittest {
	import std.math : isClose, NaN, PI_2;

	assert(cos(GDN!1(NaN(1), NaN(2))) is GDN!1(NaN(1), NaN(2)));

	const g = cos(GDN!1(PI_2));
	assert(isClose(g.val, 0., 0., real.epsilon) && g.d == -1);

	assert(isNaN(cos(GDN!1.infinity)));
	assert(isNaN(cos(-GDN!1.infinity)));

	assert(cos(GDN!2(0)) is GDN!2(1, -0., -1));
	// f = 1
	// <f',f"> = -sin(<0,1>)<1,0> = -<0,1><1,0> = <-0,-1>
}


/**
 * This function computes the sine of its argument.
 *
 * If $(MATH f(x) = sin(g(x))), then $(MATH f' = cos(g)g').
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` to compute the sine of`
 *
 * Returns:
 *   the sine expressed as a `GDN`.
 */
pragma(inline, true) pure nothrow @nogc @safe GDN!Deg sin(ulong Deg)(in GDN!Deg g)
do {
	alias cosine = Select!(Deg == 1, core.math.cos, cos);

	if (isNaN(g)) return g;
	return GDN!Deg(core.math.sin(g.val), cosine(g.reduce())*g.d);
}
/***/ unittest {
	assert(sin(GDN!2(0)) is GDN!2(0, 1, 0));
}
unittest {
	import std.math : isClose, NaN, PI_2;

	assert(sin(GDN!1(NaN(1), NaN(2))) is GDN!1(NaN(1), NaN(2)));
	assert(sin(GDN!1.zero) is GDN!1.zero);

	const g = sin(GDN!1(PI_2));
	assert(g == 1 && isClose(g.d, 0., 0., real.epsilon));

	assert(isNaN(sin(GDN!1.infinity)));
	assert(isNaN(sin(-GDN!1.infinity)));

	assert(sin(GDN!2(0)) is GDN!2(0, 1, 0));
	// f = 0
	// <f',f"> = cos(<0,1>)<1,0> = <1,0><1,0> = <1, 0>
}


/**
 * This function computes the absolute value of the argument.
 *
 * If $(MATH f(x) = |g(x)|), then $(MATH f' = sgn(g)g'), when $(MATH g ≠ 0)
 *
 * Params:
 *   Deg = the degree of the `GDN` object to compute the absolute value of
 *   g = the `GDN` object to compute the absolute value of
 *
 * Returns:
 *   the absolute value of the `GDN` object
 */
pragma(inline, true) pure nothrow @nogc @safe GDN!Deg fabs(ulong Deg)(in GDN!Deg g)
out(f; isNaN(f) || f >= 0)
do {
	const df_val = signbit(g) == 0 ? 1.0L : -1.0L;

	static if (Deg == 1)
		const df = df_val;
	else
		const df = GDN!Deg.DerivType!1.mkConst(df_val);

	return GDN!Deg(core.math.fabs(g.val), df * g.d);
}
/***/ unittest {
	assert(fabs(GDN!2(-3)) is GDN!2(3, -1, 0));
}
unittest {
	import std.math: NaN;

	assert(fabs(GDN!2(-3)) is GDN!2(3, -1, 0));
	assert(fabs(GDN!1(+0.)) is GDN!1(+0., 1));
	assert(fabs(GDN!1(-0.)) is GDN!1(+0., -1));
	assert(fabs(GDN!1(-1)) is GDN!1(1, -1));
	assert(fabs(GDN!1.nan) is GDN!1.nan);
	assert(fabs(GDN!1(-NaN(1))) is GDN!1(NaN(1)));
}


/**
 * This function computes $(MATH 2$(SUP c)g).
 *
 * If $(MATH f(x) = 2$(SUP c)g(x)), then $(MATH f' = 2$(SUP c)g').
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the generalized dual number being scaled.
 *   c = the power of $(MATH 2) used to scale `g`,
 *
 * Returns:
 *   A `GDN` object resulting from the computation.
 */
pragma(inline, true) pure nothrow @nogc @safe GDN!Deg ldexp(ulong Deg)(in GDN!Deg g, in int c)
do {
	alias ldexp_red = Select!(Deg == 1, core.math.ldexp, ldexp);
	if (isNaN(g)) return g;
	return GDN!Deg(core.math.ldexp(g.val, c), ldexp_red(g.d, c));
}
/***/ unittest {
	assert(ldexp(GDN!2(1), 2) is GDN!2(4, 4, 0));
}
unittest {
	import std.math : NaN;

	assert(ldexp(GDN!1(NaN(2)), 1) is GDN!1(NaN(2)));
}


/**
 * This function rounds `g` to a `long` using the current rounding mode. All of the derivative terms
 * are lost.
 *
 * Params:
 *   Deg = the degree of `g`
 *   g = the `GDN` object to be rounded.
 *
 * Returns:
 *   the rounded value of `g`.
 */
pragma(inline, true) pure nothrow @nogc @safe long rndtol(ulong Deg)(in GDN!Deg g)
do {
	return core.math.rndtol(g.val);
}
/***/ unittest {
	assert(rndtol(GDN!1(1.1)) == 1L);
}


/**
 * This function computes the square root of its argument.
 *
 * If $(MATH f(x) = √g(x)), then $(MATH f' = g$(SUP -½)g'/2).
 *
 * Params:
 *   Deg = the degree of the `GDN` object to compute the square root of
 *   g = the `GDN` object to compute the square root of
 *
 * Returns:
 *   the square root of the `GDN` object
 */
pragma(inline, true) pure nothrow @nogc @safe  GDN!Deg sqrt(ulong Deg)(in GDN!Deg g)
out(f; isNaN(f) || f >= 0)
do {
	if (isNaN(g)) return g;

	const dfdg = signbit(g) == 1 ? GDN!Deg.DerivType!1.nan : g.reduce()^^-0.5/2;
	return GDN!Deg(core.math.sqrt(g.val), dfdg * g.d);
}
/***/ unittest {
	assert(sqrt(GDN!2(1)) is GDN!2(1, 0.5, -0.25));
}
unittest {
	import std.format : format;

	assert(sqrt(GDN!1(-0.)) is GDN!1(-0., real.nan), "sqrt(-0) incorrect");

	const x = sqrt(-GDN!1.one);
	assert(isNaN(x), format("sqrt(-1) = %s, should be %s", x, GDN!1.nan));

	assert(sqrt(GDN!1.infinity) is GDN!1(real.infinity, 0), "sqrt(inf) incorrect");
}


/**
 * This function rounds the value of a `GDN` to a given floating point type removing all derivative
 * information.
 *
 * Params:
 *   F = the float point type to be converted to
 *   Deg = the degree of `g`
 *   g = the generalized dual number to be converted
 *
 * Returns:
 *  the rounded valued with precision determined by `F`.
 *
 */
pragma(inline, true) pure nothrow @nogc @safe
F toPrec(F, ulong Deg)(in GDN!Deg g) if (isFloatingPoint!F)
do {
	return core.math.toPrec!F(g.val);
}
/***/ unittest {
	import std.math : NaN;

	assert(toPrec!float(GDN!1(-NaN(1))) is float(-NaN(1)));
	assert(typeid(toPrec!float(GDN!1.zero)) == typeid(float));
}


/**
 * This function computes $(MATH h⋅lg(g)). It either `g` or `h` has type `real`, it is converted to
 * a constant generalized dual number with the same degree as the other parameter.
 *
 * If $(MATH f(x) = h(x)lg(g(x))), then $(MATH f' = h'lg(g) + hg'/(ln(2)g))
 *
 * Params:
 *   G = the type of `g`
 *   H = the type of `h`
 *   g = the argument of logarithm
 *   h = the multiplier of the logarithm
 *
 * Returns:
 *   The resulting generalized dual number will have a degree equal to the lesser of the degrees of
 *   `g` and `h`.
 */
pragma(inline, true) pure nothrow @nogc @safe
CommonGDN!(G, H) yl2x(G, H)(in G g, in H h)
if (anySatisfy!(isGDN, G, H) && allSatisfy!(isConvertibleToGDN, G, H))
do {
	alias Deg = typeof(return).DEGREE;

	const gg = asGDN!Deg(g);
	const hh = asGDN!Deg(h);

	if (isNaN(gg) || isNaN(hh)) return nanCombine(gg, hh);
	return yl2x_impl(gg, hh);
}
/***/ unittest {
	import std.math : LN2;

	assert(yl2x(GDN!1(2), GDN!1(3)) is GDN!1(3, 1 + 1.5/LN2));
	assert(yl2x(GDN!1(+0., -1), GDN!1(1)) is GDN!1(-real.infinity,  -real.infinity));
	assert(typeof(yl2x(GDN!2(1), GDN!1(2))).DEGREE == 1);
	assert(yl2x(GDN!1(1), 2.) is GDN!1(0, 2/LN2));
}
unittest {
	import std.math : NaN;

	assert(yl2x(GDN!1(-NaN(1), NaN(2)), GDN!1(NaN(1), NaN(3))) is GDN!1(-NaN(1), NaN(3)));
	assert(yl2x(1, GDN!1(2)) is GDN!1(0, 0));
}


private pure nothrow @nogc @safe GDN!Deg yl2x_impl(ulong Deg)(in GDN!Deg g, in GDN!Deg h)
do {
	alias yl2x_red = Select!(Deg == 1, core.math.yl2x, yl2x_impl);

	GDN!Deg.DerivType!1 df;
	if (signbit(g) == 0) {
		const g_red = g.reduce();
		df = yl2x_red(g_red, h.d) + h.reduce()*g.d/(LN2 * g_red);
	}

	return GDN!Deg(core.math.yl2x(g.val, h.val), df);
}
unittest {
	assert(isNaN(yl2x(GDN!1(-1), GDN!1(1))));

	const f = yl2x(GDN!1(0), GDN!1(1));
	assert(f == -real.infinity && isNaN(f.d));

	assert(yl2x(GDN!2(1), GDN!2(2)) is GDN!2(0, 2/LN2, 0));
	// f = 0
	// <f',f"> = h'lg(g) + hg'/(ln(2)g)
	//    = <1,0>lg(<1,1>) + <2,1><1,0>/ln(2)<1,1>
	//    = <0,0+1/ln(2)> + <2,1>/<1,1>/ln(2)
	//    = <0,1/ln(2)> + <2,-1>/ln(2)
	//    = <2/ln(2),0>
}


/**
 * Computes $(MATH h⋅lg(g + 1)), for $(MATH -(1 - √½) ≤ x ≤ +(1 - √½)). When $(MATH g) is outside of
 * this interval, the results are undefined. It either `g` or `h` has type `real`, it is converted
 * to a constant generalized dual number with the same degree as the other parameter.
 *
 * If $(MATH f(x) = h(x)lg(g(x) + 1)), then $(MATH f' = h'lg(g + 1) + hg'/[ln(2)(g + 1)])
 *
 * Params:
 *   G = the type of `g`
 *   H = the type of `h`
 *   g = the argument of logarithm
 *   h = the multiplier of the logarithm
 *
 * Returns:
 *   The resulting generalized dual number will have a degree equal to the lesser of the degree of
 *   `g` and `h`.
 */
pragma(inline, true) pure nothrow @nogc @safe
CommonGDN!(G, H) yl2xp1(G, H)(in G g, in H h)
if (anySatisfy!(isGDN, G, H) && allSatisfy!(isConvertibleToGDN, G, H))
do {
	alias Deg = typeof(return).DEGREE;

	const gg = asGDN!Deg(g);
	const hh = asGDN!Deg(h);

	if (isNaN(gg) || isNaN(hh)) return nanCombine(gg, hh);
	return yl2xp1_impl(gg, hh);
}
/***/ unittest {
	import std.math : LN2;

	assert(yl2xp1(GDN!1(0), GDN!1(3)) is GDN!1(0, 3/LN2));
	assert(typeof(yl2xp1(GDN!2(0), GDN!1(1))).DEGREE == 1);
	assert(yl2xp1(GDN!1(0), 1) is GDN!1(0, 1/LN2));
}
unittest {
	import std.math : NaN;

	assert(yl2xp1(GDN!1(-NaN(1), NaN(2)), GDN!1(NaN(1), NaN(3))) is GDN!1(-NaN(1), NaN(3)));
	assert(yl2xp1(0, GDN!1(1)) is GDN!1(0, 0));
}


private pure nothrow @nogc @safe GDN!Deg yl2xp1_impl(ulong Deg)(in GDN!Deg g, in GDN!Deg h)
do {
	alias yl2xp1_red = Select!(Deg == 1, core.math.yl2xp1, yl2xp1_impl);

	GDN!Deg.DerivType!1 df;
	if (g > -1) {
		const g_red = g.reduce();
		df = yl2xp1_red(g_red, h.d) + h.reduce()*g.d/(LN2 * (g_red + 1));
	}

	return GDN!Deg(core.math.yl2xp1(g.val, h.val), df);
}
unittest {
	assert(yl2xp1_impl(GDN!2(0), GDN!2(1)) is GDN!2(0, 1/LN2, 1/LN2));
	// f = 0
	// <f',f"> = h'lg(g+1) + hg'/[ln(2)(g+1)]
	//    = <1,0>lg(<0,1>+1) + <1,1><1,0>/[ln(2)(<0,1>+1)]
	//    = <1,0>lg<1,1> + <1,1>/[ln(2)<1,1>]
	//    = <0,1/ln(2)> + <1,0>/ln(2)
	//    = <1/ln(2),1/ln(2)>

	assert(isNaN(yl2xp1_impl(GDN!1(-1), GDN!1(0)).d));

	const q = yl2xp1_impl(GDN!1(-2), GDN!1(0));
	assert(isNaN(q.d), "yl2xp1(2,0) should not have a derivative");
}
