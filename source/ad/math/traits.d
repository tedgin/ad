/**
* This module extends `std.math.traits` module to support `GDN` objects.
*/
module ad.math.traits;

static import std.math.traits;

import std.traits : isFloatingPoint, isIntegral, Select;

static import ad.internal;

import ad;
import ad.internal : asReal, isConvertibleToGDN;


/**
* This function checks if the given `GDN` object has a finite value.
*
* Params:
*   Deg = The degree of the `GDN` object
*   f = The `GDN` object to check
*
* Returns:
*   `true` if the value of the `GDN` object is finite, `false` otherwise
*/
pure nothrow @nogc @safe bool isFinite(ulong Deg)(in GDN!Deg f)
do {
	return ad.internal.isFinite(f);
}
/***/ unittest {
	assert(isFinite(GDN!1(1)));
}


/**
* This function checks if two `GDN` objects have the same binary representation.
*
* Params:
*   FDeg = The degree of the first `GDN` object.
*   GDeg = The degree of the second `GDN` object.
*   f = The first `GDN` object to compare.
*   g = The second `GDN` object to compare.
*
* Returns:
*   `true` if the two `GDN` objects are identical, `false` otherwise.
*/
pure nothrow @nogc @safe bool isIdentical(ulong FDeg, ulong GDeg)(in GDN!FDeg f, in GDN!GDeg g)
do {
	alias isIdenticalDeriv = Select!(FDeg == 1, std.math.traits.isIdentical, isIdentical);

	static if (FDeg == GDeg)
		return std.math.traits.isIdentical(f.val, g.val) && isIdenticalDeriv(f.d, g.d);
	else
		return false;
}
/***/ unittest {
	assert(isIdentical(GDN!1(1), GDN!1(1)));
	assert(!isIdentical(GDN!1(1), GDN!1(2)));
	assert(!isIdentical(GDN!1(1, 2), GDN!1(1, 1)));
	assert(!isIdentical(GDN!1(1), GDN!2(1)));
}


/**
* This function checks if the given `GDN` object has an infinite value.
*
* Params:
*   Deg = The degree of the `GDN` object
*   f = The `GDN` object to check
*
* Returns:
*   `true` if the value of the `GDN` object is infinite, `false` otherwise
*/
pure nothrow @nogc @safe bool isInfinity(ulong Deg)(in GDN!Deg f)
do {
	return ad.internal.isInfinity(f);
}
/***/ unittest {
	assert(isInfinity(GDN!1(real.infinity)));
}


/**
* This function determines whether the value of the given `GDN` object is NaN.
*
* Params:
*   Deg = the degree of the `GDN` object
*   f = the `GDN` object to check
*
* Returns:
*   `true` if the value of the `GDN` object is `NaN`, `false` otherwise
*/
pure nothrow @nogc @safe bool isNaN(ulong Deg)(in GDN!Deg f)
do {
	return ad.internal.isNaN(f);
}
/***/ unittest {
	assert(isNaN(GDN!1.nan));
}


/**
* This function checks if the given `GDN` object has a normal value.
*
* A normal value is one that is finite, non-zero, and not subnormal.
*
* Params:
*   Deg = The degree of the `GDN` object
*   f = The `GDN` object to check.
*
* Returns:
*   `true` if the value of the `GDN` object is normal, `false` otherwise
*/
pure nothrow @nogc @safe bool isNormal(ulong Deg)(in GDN!Deg f)
do {
	return std.math.traits.isNormal(f.val);
}
/***/ unittest {
	assert(isNormal(GDN!1(1)));
}


/**
* This function checks if the given `GDN` object is a power of 2.
*
* Params:
*   Deg = The degree of the `GDN` object
*   f = The `GDN` object to check
*
* Returns:
*   `true` if the value of the `GDN` object is a power of 2, `false` otherwise
*/
pure nothrow @nogc @safe bool isPowerOf2(ulong Deg)(in GDN!Deg f)
do {
	return std.math.traits.isPowerOf2(f.val);
}
/***/ unittest {
	assert(isPowerOf2(GDN!1(1)));
}


/**
* This function checks if the given `GDN` object has a subnormal value.
*
* Params:
*   Deg = The degree of the `GDN` object
*   f = The `GDN` object to check
*
* Returns:
*   `true` if the value of the `GDN` object is subnormal, `false` otherwise
*/
pure nothrow @nogc @safe bool isSubnormal(ulong Deg)(in GDN!Deg f)
do {
	return std.math.traits.isSubnormal(f.val);
}
/***/ unittest {
	assert(isSubnormal(GDN!1(real.min_normal / 2)));
}


/**
* This function checks if the sign bit of the value of a given `GDN` object is set.
*
* Params:
*   Deg = The degree of the `GDN` object
*   f = The `GDN` object to check
*
* Returns:
*   `1` if the sign bit of the `GDN` object's value is set, `0` otherwise
*/
pure nothrow @nogc @safe int signbit(ulong Deg)(in GDN!Deg f)
do {
	return ad.internal.signbit(f);
}
/***/ unittest {
	assert(signbit(GDN!1(-1.0)) == 1);
}


/**
* This function makes to have the same sign as from.
*
* Params:
*   G = GDN or implicitly convertible _to real
*   F = a floating-point type
*   I = an integral type
*   TDeg = the degree of the `GDN` object _to change the sign of
*   FDeg = the degree of the `GDN` object _to copy the sign _from
*   to = the value _to change the sign of
*   from = the value _to copy the sign _from
*
* Returns:
*   to with the same sign as from
*/
pure nothrow @nogc @safe
GDN!TDeg copysign(G, ulong TDeg)(in GDN!TDeg to, in G from) if (isConvertibleToGDN!G)
do {
	return GDN!TDeg(std.math.traits.copysign(to.val, asReal(from)), to.d);
}
/// ditto
pure nothrow @nogc @safe F copysign(F, ulong FDeg)(in F to, in GDN!FDeg from) if (isFloatingPoint!F)
do {
	return std.math.traits.copysign(to, from.val);
}
/// ditto
pure nothrow @nogc @safe real copysign(I, ulong FDeg)(in I to, in GDN!FDeg from) if (isIntegral!I)
do {
	return std.math.traits.copysign(to, from.val);
}
/***/ unittest {
	assert(isIdentical(copysign(GDN!1(-1), GDN!1(-2)), GDN!1(-1)));
	assert(isIdentical(copysign(GDN!1(-3), 4.), GDN!1(3)));
	assert(copysign(5., GDN!1(-6.)) is -5.);
	assert(copysign(7, GDN!1(8)) is 7);
}
unittest {
	assert(isIdentical(copysign(GDN!3(1), GDN!1(-1)), GDN!3(-1)));
}


/**
* This function computes the sign of a `GDN` Object.
*
* If $(MATH f(x) = sgn(g(x))), then $(MATH f' = 2𝛿(g)g'), where $(MATH 𝛿) is the Dirac delta
* function.
*
* To be in agreement with `std.math.traits.sgn`, the sign of  $(MATH sgn(±0) = ±0).
*
* Params:
*   Deg = the degree of the `GDN` object
*   g = the `GDN` object to compute the sign of
*
* Returns:
*   the sign of the `GDN` object
*/
pure nothrow @nogc @safe GDN!Deg sgn(ulong Deg)(in GDN!Deg g)
do {
	return ad.internal.sgn(g);
}
/***/ unittest {
	assert(isIdentical(sgn(GDN!1(-2)), GDN!1(-1, 0)));
	assert(isIdentical(sgn(GDN!1(0)), GDN!1(0, real.infinity)));
}
