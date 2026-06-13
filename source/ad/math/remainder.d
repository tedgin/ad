/**
* This extends the `std.math.remainder` Phobos module to support `GDN` objects.
*/
module ad.math.remainder;

static import std.math.remainder;

import std.math : copysign;
import std.meta : allSatisfy, anySatisfy;

import ad;
import ad.internal :
	asGDN, CommonGDN, isConvertibleToGDN, isFinite, isGDN, isInfinity, isNaN, round, trunc;


/**
* This function determines the remainder of g divided by h.
*
* It is the same as `g % h`.
*
* Params:
*   G = the type of g, either `GDN` or `real`
*   H = the type of h, either `GDN` or `real`
*   g = the dividend
*   h = the divisor
*
* Returns:
*   the remainder of g divided by h
*/
pure nothrow @nogc @safe
CommonGDN!(G, H) fmod(G, H)(in G g, in H h)
if (anySatisfy!(isGDN, G, H) && allSatisfy!(isConvertibleToGDN, G, H))
do {
	return g % h;
}
/***/ unittest {
	const f = fmod(GDN!1(5), GDN!1(3));
	assert(f is GDN!1(2,0));
}


/**
* This function breaks g into an integer and a fraction, each with the same sign as g.
*
* `f = modf(g, i)` can be expressed mathematically as $(MATH f = g - i). $(MATH i) is defined as
* follows. When $(MATH g ≥ 0), $(MATH i = ⌊g⌋), and when $(MATH g < 0), $(MATH i = ⌈g⌉). <em>This is
* the mathematical definition of the function `trunc` defined in `ad.math.rounding`.</em>
*
* Params:
*   Deg = the degree of g
*   g = the GDN object to break into an integer and a fraction
*   i = the integer part of g
*
* Returns:
*   the fractional part of g
*/
pure nothrow @nogc @safe GDN!Deg modf(ulong Deg)(in GDN!Deg g, out GDN!Deg i)
do {
	i = trunc(g);

	if (isNaN(g)) return g;

	if (isInfinity(g)) {
		return GDN!Deg(copysign(0.0L, g.val), GDN!Deg.mkNaNDeriv());
	}

	return g - i;
}
/***/ unittest {
	import std.math : isClose;

	const g = GDN!1(3.14159);
	GDN!1 i;
	const f = modf(GDN!1(g), i);
	assert(i is GDN!1(3, 0));
	assert(isClose(f.val, 0.14159));
	assert(f.d == 1);
}
unittest {
	import std.format : format;
	import std.math : NaN;

	GDN!1 i;

	const q = modf(GDN!1(-real.infinity), i);
	assert(
		q is GDN!1(-0., real.nan) && i is GDN!1(-real.infinity, real.nan),
		format("q: %s, i: %s", q, i));

	const w = modf(GDN!1(real.infinity), i);
	assert(w is GDN!1(0., real.nan) && i is GDN!1(real.infinity, real.nan));

	assert(modf(GDN!1(NaN(1)), i) is GDN!1(NaN(1)) && i is GDN!1(NaN(1)));
}


/**
* This function calculates the integer quotient of g and h and its remainder.
*
* It uses the definition of remainder provided by IEC 60559. The integer quotient n is
* `round(g/h)`, and the remainder is defined as `g - h*n`, where function `round` is defined in
* `ad.math.rounding`.
*
* If either g or h has type `real`, it is converted to a constant generalized dual number with the
* same degree as the other parameter. If g and h are `GDN` objects with different degrees, the one
* with the greater degree is converted to have the same degree as the lesser.
*
* Params:
*   G = the type of g, either `GDN` or `real`
*   H = the type of h, either `GDN` or `real`
*   g = the dividend
*   h = the divisor
*   n = the integer quotient of g and h
*
* Returns:
*   the remainder of g divided by h
*/
pure nothrow @nogc @safe
CommonGDN!(G, H) remquo(G, H)(in G g, in H h, out int n)
if (anySatisfy!(isGDN, G, H) && allSatisfy!(isConvertibleToGDN, G, H))
do {
	alias Deg = typeof(return).DEGREE;

	const gg = asGDN!Deg(g);
	const hh = asGDN!Deg(h);

	if (isNaN(gg) || isNaN(hh)) return nanCombine(gg, hh);

	if (gg == 0.0L && hh != 0.0L) {
		n = 0;
		return gg;
	}

	if (isFinite(gg) && isInfinity(hh)) return gg;

	const n_gdn = round(gg / hh);
	n = cast(int) n_gdn.val;
	return gg - hh * n_gdn;
}
/***/ unittest {
	import std.math : isClose;

	int n;
	const f = remquo(GDN!1(5.1), GDN!1(3), n);
	assert(n == 2 && isClose(f.val, -0.9) && f.d == -1);
}
unittest {
	import std.math : NaN;
	import ad.math.traits : isNaN;

	int n;

	n = int.min;
	const q = remquo(GDN!1(-0.), GDN!1(1), n);
	assert(n == 0 && q is GDN!1(-0., 1));

	n = int.min;
	const w = remquo(GDN!1(+0.), GDN!1(1), n);
	assert(n == 0 && w is GDN!1(+0., 1));

	assert(isNaN(remquo(GDN!1(-real.infinity), GDN!1(1), n)));
	assert(isNaN(remquo(GDN!1(real.infinity), GDN!1(1), n)));
	assert(isNaN(remquo(GDN!1(1), GDN!1(0), n)));
	assert(remquo(GDN!1(2), GDN!1(-real.infinity), n) is GDN!1(2));

	assert(remquo(GDN!1(NaN(1), NaN(3)), GDN!1(NaN(4), NaN(2)), n) is GDN!1(NaN(4), NaN(3)));
}


/**
* This function calculates the _remainder of g divided by h.
*
* It using the definition of _remainder provided by IEC 60559. In other words, it computes
* `g - h*round(g/h)`, where the `round` function is defined in `ad.math.rounding`.
*
* If either g or h has type `real`, it is converted to a constant generalized dual number with the
* same degree as the other parameter. If g and h are `GDN` objects with different degrees, the one
* with the greater degree is converted to have the same degree as the lesser.

* Params:
*   G = the type of g, either `GDN` or `real`
*   H = the type of h, either `GDN` or `real`
*   g = the dividend
*   h = the divisor
*
* Returns:
*   the remainder of g divided by h
*/
pure nothrow @nogc @safe
CommonGDN!(G, H) remainder(G, H)(in G g, in H h)
if (anySatisfy!(isGDN, G, H) && allSatisfy!(isConvertibleToGDN, G, H))
do {
	int _;
	return remquo(g, h, _);
}
/***/ unittest {
	import std.math : isClose;

	const f = remainder(GDN!1(5.1), GDN!1(3));
	assert(isClose(f.val, -0.9));
	assert(f.d == -1);
}
