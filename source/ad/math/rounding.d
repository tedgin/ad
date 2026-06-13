/**
* This module extends `std.math.rounding` module to support `GDN` objects.
*/
module ad.math.rounding;

static import std.math.rounding;

import std.math : isInfinity, signbit;
import std.meta : allSatisfy, anySatisfy;
import std.traits : isIntegral;

static import ad.core.math;
static import ad.internal;

import ad;
import ad.internal :
	asGDN, CommonGDN, dirac, isConvertibleToGDN, isGDN, isNaN, nextDown, nextUp, pow;


/**
* This function determines the value of g rounded upward to the nearest integer.
*
* If $(MATH f(x) = ⌈g(x)⌉), then $(MATH f' = g'∑$(SUB i∊ℤ)𝛿(g-i))
*
* Params:
*   Deg = the degree of `g`
*   g = the `GDN` to round
*
* Returns:
*   the rounded `GDN`
*/
pure nothrow @nogc @safe GDN!Deg ceil(ulong Deg)(in GDN!Deg g)
do {
	return ad.internal.ceil(g);
}
/***/ unittest {
	assert(ceil(GDN!1(1)) is GDN!1(1, real.infinity));
	assert(ceil(GDN!1(-1.4)) is GDN!1(-1, 0));
}


/**
* This function determins the value of g rounded downward to the nearest integer.
*
* If $(MATH f(x) = ⌊g(x)⌋), then $(MATH f' = g'∑$(SUB i∊ℤ)𝛿(g-i))
*
* Params:
*   Deg = the degree of g
*   g = the `GDN` to round
*
* Returns:
*   the rounded `GDN`
*/
pure nothrow @nogc @safe GDN!Deg floor(ulong Deg)(in GDN!Deg g)
do {
	return ad.internal.floor(g);
}
/***/ unittest {
	assert(floor(GDN!1(1)) is GDN!1(1, real.infinity));
	assert(floor(GDN!1(-1.4)) is GDN!1(-2, 0));
}


/**
* This function rounds g to the nearest integer value, using the current rounding mode.
*
* Params:
*   Deg = the degree of g
*   g = the `GDN` object to be rounded.
*
* Returns:
*   an integer representing the rounded value of g.
*/
pure nothrow @nogc @safe long lrint(ulong Deg)(in GDN!Deg g)
do {
	return std.math.rounding.lrint(g.val);
}
/***/ unittest {
	assert(lrint(GDN!1(1.9)) == 2L);
}


/**
* This function determines the value of a `GDN` rounded to the nearest integer.
*
* Params:
*   Deg = the degree of g
*   g = the `GDN` object to be rounded.
* Returns:
*   An integer representing the rounded value of g.
*/
nothrow @nogc @safe long lround(ulong Deg)(in GDN!Deg g)
do {
	return std.math.rounding.lround(g.val);
}
/***/ unittest {
	assert(lround(GDN!1(1.5)) == 2L);
	assert(lround(GDN!1(-0.5)) == -1L);
}


private pure nothrow @nogc @safe GDN!Deg nearbyint_impl(string impl, ulong Deg)(in GDN!Deg g)
do {
	if (isNaN(g)) return g;

	mixin("const f = " ~ impl ~ "(g.val);");

	auto dfdg = GDN!Deg.mkZeroDeriv();
	if (isInfinity(g.val)) {
		dfdg = GDN!Deg.mkNaNDeriv();
	} else {
		auto fn = std.math.rounding.nearbyint(nextDown(g).val);
		auto fp = std.math.rounding.nearbyint(nextUp(g).val);

		if (f == 0.0L) {
			if (signbit(f) == 1) {
				fp = f;
			} else {
				fn = f;
			}
		}

		if (fn != fp) {
			dfdg = dirac(g.reduce() - g.val);
		}
	}

	return GDN!Deg(f, dfdg*g.d);
}
unittest {
	import std.math : isNaN, NaN;

	enum impl = "std.math.rounding.nearbyint";

	assert(nearbyint_impl!impl(GDN!1.infinity) is GDN!1(real.infinity, real.nan));

	const f = nearbyint_impl!impl(GDN!2(1.5));
	assert(f == 2 && f.d == real.infinity && isNaN(f.d!2));

	assert(nearbyint_impl!impl(GDN!1(-0.)) is GDN!1(-0., 0));
	assert(nearbyint_impl!impl(GDN!1(+0.)) is GDN!1(+0., 0));

	assert(nearbyint_impl!impl(GDN!1(NaN(1))) is GDN!1(NaN(1)));
}


/**
* This function rounds g to the nearest integer value, using the current rounding mode.
*
* If $(MATH f(x) = nearbyint(g(x))), then $(MATH f' = (df/dg)g'), where $(MATH df/dg = 𝛿(g - m)),
* $(MATH 𝛿) is the Dirac delta function, and $(MATH m) is a rounding mode split point.
*
* Params:
*   Deg = the degree of g
*   g = the `GDN` object to be rounded
*
* Returns:
*   a `GDN` object representing the rounded value of g
*/
pure nothrow @nogc @safe GDN!Deg nearbyint(ulong Deg)(in GDN!Deg g)
do {
	return nearbyint_impl!"std.math.rounding.nearbyint"(g);
}
/***/ unittest {
	import std.math : isNaN;

	const e = nearbyint(GDN!2(1.5));
	assert(e == 2 && e.d == real.infinity && isNaN(e.d!2));
}


/**
* This function rounds g to the nearest integer value, using the current rounding mode.
*
* If the return value is not identical to g, the `FE_INEXACT` exception is raised.
*
* If $(MATH f(x) = rint(g(x))), then $(MATH f' = (df/dg)g'), where $(MATH df/dg = 𝛿(g - m)),
* $(MATH 𝛿) is the Dirac delta function, and $(MATH m) is a rounding mode split point.
*
* Params:
*   Deg = the degree of g
*   g = the `GDN` object to be rounded.
*
* Returns:
*   a `GDN` object representing the rounded value of g
*/
pure nothrow @nogc @safe GDN!Deg rint(ulong Deg)(in GDN!Deg g)
do {
   return nearbyint_impl!"std.math.rounding.rint"(g);
}
/***/ unittest {
	import std.math : ieeeFlags, isNaN, resetIeeeFlags;

	resetIeeeFlags();
	const e = rint(GDN!2(1.5));
	assert(ieeeFlags.inexact);
	assert(e == 2 && e.d == real.infinity && isNaN(e.d!2));
}
unittest {
	import std.math : ieeeFlags, resetIeeeFlags;

	resetIeeeFlags();
	const w = rint(GDN!1.one);
	assert(!ieeeFlags.inexact);
	assert(w is GDN!1.one);
}


/**
* This function round val to a multiple of unit using the function rfunc.
*
* Params:
*   V = the value type
*   U = the unit type
*   rfunc = the rounding function to use
*   val = the value to round
*   unit = the unit to round to
*
* Returns:
*   the rounded value of val to the nearest multiple of unit
*/
pure nothrow @nogc @safe
CommonGDN!(U, V) quantize(alias rfunc=rint, U, V)(in V val, in U unit)
if (is(typeof(rfunc(CommonGDN!(U, V).init)) : CommonGDN!(U, V)))
do {
	return quantize_impl!rfunc(val, unit);
}
/***/ unittest {
	const q = quantize!nearbyint(GDN!1(5), GDN!1(3));
	assert(q is GDN!1(6, 2));
}


/**
* This function rounds g to a multiple of `base ^^ exp` using the function rfunc.
*
* Params:
*   G = the type of g
*   I = the exponent type
*   rfunc = the rounding function to use
*   base = the base of the number to round to
*   g = the `GDN` object to round
*   exp = the exponent of the number to round to
*
* Returns:
*   the rounded `GDN` object
*/
pure nothrow @nogc @safe
CommonGDN!(G, typeof(base)) quantize(alias base, alias rfunc=rint, G, I)(in G g, in I exp)
if (anySatisfy!(isGDN, G, typeof(base))
	&& allSatisfy!(isConvertibleToGDN, G, typeof(base))
	&& (is(typeof(rfunc(G.init)) : G) || is(typeof(rfunc(typeof(base).init)) : typeof(base)))
	&& isIntegral!I)
do {
	alias Deg = typeof(return).DEGREE;
	enum b = asGDN!Deg(base);

	const gg = asGDN!Deg(g);

	if (isNaN(b) || isNaN(gg)) return nanCombine(b, gg);
	return quantize_impl!rfunc(gg, pow(b, exp));
}
/// ditto
pure nothrow @nogc @safe
CommonGDN!(G, typeof(base)) quantize(alias base, long exp=1, alias rfunc=rint, G)(in G g)
if (anySatisfy!(isGDN, G, typeof(base))
	&& allSatisfy!(isConvertibleToGDN, G, typeof(base))
	&& (is(typeof(rfunc(G.init)) : G) || is(typeof(rfunc(typeof(base).init)) : typeof(base))))
do {
	alias Deg = typeof(return).DEGREE;
	enum b = asGDN!Deg(base);
	enum unit = pow(b, exp);

	const gg = asGDN!Deg(g);

	if (isNaN(b) || isNaN(gg)) return nanCombine(b, gg);
	return quantize_impl!rfunc(asGDN!Deg(g), unit);
}
/***/ unittest {
	import ad.math.operations : isClose;

	const f = quantize!10(GDN!1(345.678_9), -2);
	assert(isClose(f, 345.68) && f.d == 0);

	const g = quantize!(GDN!1(2))(GDN!1(1.6), -1);
	assert(g is GDN!1(1.5, -0.75));

	assert(quantize!22(GDN!1(12_345.678_9)) is GDN!1(12_342, 0));
}
unittest {
	import std.math : NaN;

	assert(quantize!(GDN!1(NaN(1), NaN(3)))(GDN!1(NaN(4), NaN(2)), 0) is GDN!1(NaN(4),NaN(3)));
	assert(quantize!(GDN!1(NaN(1), NaN(3)))(GDN!1(NaN(4), NaN(2))) is GDN!1(NaN(4), NaN(3)));
}


private pure nothrow @nogc @safe
GDN!Deg quantize_impl(alias round, ulong Deg)(in GDN!Deg val, in GDN!Deg unit)
if (is(typeof(round(GDN!Deg.init)) : GDN!Deg))
do {
	return round(val/unit) * unit;
}
unittest {
	const f = quantize_impl!rint(GDN!1(1.5), GDN!1(0.5, 0));
	assert(f is GDN!1(1.5, 0));
}


/**
* This function rounds g to a `long` using the current rounding mode.
*
* All of the derivative terms are lost.
*
* Params:
*   Deg = the degree of g
*   g = the `GDN` object to be rounded.
*
* Returns:
*   the rounded value of g.
*/
pure nothrow @nogc @safe long rndtol(ulong Deg)(in GDN!Deg g)
do {
	return ad.core.math.rndtol(g);
}
/***/ unittest {
	assert(rndtol(GDN!1(0.2)) == 0L);
}


/**
* This function returns a value g rounded to the nearest integer.
*
* If $(MATH f(x) = round(g(x))), then $(MATH f' = g'∑$(SUB i∊ℤ)𝛿(g-i-½)),
*
* Params:
*   Deg = the degree of g
*   g = the `GDN` to round.
*
* Returns:
*   the rounded `GDN`
*/
nothrow @nogc @trusted GDN!Deg round(ulong Deg)(in GDN!Deg g)
do {
	return ad.internal.round(g);
}
/***/ unittest {
	assert(round(GDN!1(4.5)) is GDN!1(5, real.infinity));
	assert(round(GDN!1(-4.5)) is GDN!1(-5, real.infinity));
}


/**
* This function truncates g to an integer.
*
* Where $(MATH g(x) < 0), $(MATH f(x) = ⌈g(x)⌉), otherwise $(MATH f(x) = ⌊g(x)⌋).
*
* Params:
*   Deg = the degree of g
*   g = the `GDN` to round.
*
* Returns:
*   the truncated `GDN`
*/
pure nothrow @nogc @trusted GDN!Deg trunc(ulong Deg)(in GDN!Deg g)
do {
	return ad.internal.trunc(g);
}
/***/ unittest {
	assert(trunc(GDN!1(0.01)) is GDN!1(+0., 0));
	assert(trunc(GDN!1(-0.49)) is GDN!1(-0., 0));
}
