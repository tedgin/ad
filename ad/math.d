/**
 * This module extends the std.math library to supporting PluralNum objects.
 */
module ad.math;

import std.math;
import std.traits;

public import ad.core;


/**
 * This function determines whether its argument is a NaN.  It is analogous to std.math.isNaN().
 *  
 * Params:
 *   T = The type of the argument. It must be a scalar type for a PluralNum.
 *   x = the argument
 */
@safe
export pure nothrow bool isNaN(T)(in T x) if (isScalarType!(T))
body {
	return std.math.isNaN(x);
}
@safe
export pure nothrow bool isNaN(ulong O, F)(in PluralNum!(O, F) x) 
body {
	return std.math.isNaN(x.val);
}
unittest {
	assert(isNaN(PluralNum!().nan));
	assert(isNaN(PluralNum!(1, float).nan));
}


/**
 * This function computes the sign of the argument. It is analogous to std.math.sgn().
 * 
 * The derivate of sgn(x) evaluated at 0 is undefined or NaN.  Otherwise it is 0*x'.
 * 
 * Params:
 *   T = The type of the argument. It must be a scalar type for a PluralNum.
 *   x = the argument
 */
@safe
export pure nothrow T sgn(T)(in T x) if (isScalarType!(T))
body { 
	return std.math.sgn(x); 
}
@safe
export pure nothrow PluralNum!(O, F) sgn(ulong O, F)(in PluralNum!(O, F) x)
body {
	alias PN = PluralNum!(O, F);
	if (x == 0) return PN(0, PN.DerivType!().nan);
	return PN(sgn(x.val), 0 * x.d);
}
unittest {
	assert(sgn(0) == 0);
	assert(sgn(2) == 1);
	assert(sgn(-3) == -1);

	assert(sgn(PluralNum!()()).same(PluralNum!()(real.nan, real.nan)));

	assert(sgn(PluralNum!()(0, real.nan)).same(PluralNum!()(0, real.nan)));
	assert(sgn(derivSeq(1.0L, 2.0L, real.nan)).same(derivSeq(1.0L, 0.0L, real.nan)));
	assert(sgn(PluralNum!()(real.infinity, 0)).same(PluralNum!().one));
	assert(sgn(PluralNum!()(-real.infinity, 0)).same(-PluralNum!().one));
}


/**
 * This function computes the absolute value of the argument. It is analogous to std.math.abs().
 * 
 * The derivative of the abs(x) is x*x'/abs(x).
 * 
 * Params:
 *   T = The type of the argument. It must be a scalar type for a PluralNum.
 *   x = the argument
 */
@safe
export pure nothrow T abs(T)(in T x) if (isScalarType!(T))
body {
	return std.math.abs(x);
}
@safe
export pure nothrow PluralNum!(O, F) abs(ulong O, F)(in PluralNum!(O, F) x)
body {
	alias PN = PluralNum!(O, F);
	const dx = x == PN.zero 
		? PN.DerivType!().nan 
		: x.d * sgn(x.reduce()); 
	return PN(std.math.abs(x.val), dx);
}
unittest {
	assert(abs(1) == 1);
	assert(abs(0) == 0);
	assert(abs(-1) == 1);

	assert(isNaN(abs(PluralNum!()())));

	// TODO fix this test to use same
	assert(abs(PluralNum!().var(-1)) == PluralNum!().var(1));
	// TODO more testing
}


/** 
 * TODO document
 * 
 * The derivative of cos(x) is -sin(x).
 */
@safe 
export pure nothrow T cos(T)(in T x) if (isScalarType!(T))
body {
	return std.math.cos(x);
}
@safe
export pure nothrow PluralNum!(O, F) cos(ulong O, F)(in PluralNum!(O, F) x)
body {
	return DerivSeqType(
		std.math.cos(x.val), 
		-x.d * sin(x.reduce()));
}
// TODO test


/+
export pure DerivSeqType sqrt(DerivSeqType)(const DerivSeqType u)
in {
	static if (DerivSeqType.Order == 0) {
		assert (u.val() >= 0);
	} else {
		assert (u.val() > 0);
	}
}
body {
	static if (DerivSeqType.Order == 0) {
		return DerivSeqType(std.math.sqrt(u.val()));
	} else {
		return DerivSeqType( 
			std.math.sqrt(u.val()), 
			u.d() / (2 * sqrt(u.reduce())));
	}		
}


export pure DerivSeqType pow(ValType, DerivSeqType)(const DerivSeqType u, const ValType k)
	/* TODO Depending on the values of _x and k, a complex or nan may arise.  Either way, this 
	 * should fail an assertion.
	 */
body {
	static if (DerivSeqType.Order == 0) {
		return DerivSeqType(std.math.pow(u.val(), k));
	} else {
		return DerivSeqType(
			std.math.pow(u.val(), k), 
			k * pow(u.reduce(), k - 1) * u.d());
	}
}


export pure DerivSeqType exp(DerivSeqType)(const DerivSeqType u)
body {
	static if (DerivSeqType.Order == 0) {
		return DerivSeqType(std.math.exp(u.val()));
	} else {
		return DerivSeqType( 
			std.math.exp(u.val()), 
			u.d() * exp(u.reduce()));
	}
}


export pure DerivSeqType log(DerivSeqType)(const DerivSeqType u)
body {
	return u.log();
}


export pure DerivSeqType sin(DerivSeqType)(const DerivSeqType u)
body {
	static if (DerivSeqType.Order == 0) {
		return DerivSeqType(std.math.sin(u.val()));
	} else {
		return DerivSeqType(
			std.math.sin(u.val()), 
			u.d() * cos(u.reduce()));
	}
}
+/

/+
// TODO learn why remainder is not pure
export DerivSeqType tan(DerivSeqType)(const DerivSeqType u)
in {
	assert(remainder(u.val() - PI_2, PI) != 0);
}
body {
	static if (DerivSeqType.Order == 0) {
		return DerivSeqType(std.math.tan(u.val()));
	} else {
		return sin(u) / cos(u);
	}
}


export pure DerivSeqType acos(DerivSeqType)(const DerivSeqType u)
in {
	static if (DerivSeqType.Order == 0) {
		assert(-1 <= u.val() && u.val() <= 1);
	} else {
		assert(-1 < u.val() && u.val() < 1);
	}
}
body {
	static if (DerivSeqType.Order == 0) {
		return std.math.acos(u.val());
	} else {
		return DerivSeq( 
			std.math.acos(u.val()), 
			-u.d() / sqrt(1 - pow(u.reduce(), 2)));
	}
}

unittest {
	const real c = 5;
	const auto p = DerivSeq!(1).param(1);
	const auto u = DerivSeq!(2).var(2);
	const auto v = DerivSeq!(2).var(3);

	writeln("c = ", c);
	writeln("p = ", p);
	writeln("u = ", u);
	writeln("v = ", v);
	writeln("u.x() = ", u.val());
	writeln("u.dx() = ", u.d!()());
	writeln("u.dx(2) = ", u.d!(2)());
	writeln("+u = ", +u);
	writeln("-u = ", -u);
	writeln("u.inv() = ", u.inv());
	writeln("u + v = ", u + v);
	writeln("u + c = ", u + c);
	writeln("u - v = ", u - v);
	writeln("u - c = ", u - c);
	writeln("c - u = ", c - u);
	writeln("u * u = ", u * u);
	writeln("u * c = ", u * c);
	writeln("u / v = ", u / v);
	writeln("u / c = ", u / c);
	writeln("c / u = ", c / u);
	writeln("sin(u) = ", sin(u));
	writeln("cos(u) = ", cos(u));
	writeln("tan(u) = ", tan(u));
	writeln("exp(u) = ", exp(u));
	writeln("log(u) = ", log(u));
	writeln("pow(u, c) = ", pow(u, c));
	writeln("abs(u) = ", abs(u));
	writeln("cos(2u) = ", cos(2*u));
	writeln("pow(cos(u), 2) = ", pow(cos(u), 2));
	writeln("cos(u)*cos(u) = ", cos(u)*cos(u));
}
+/