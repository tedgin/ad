/**
 * This module extends the std.math library to supporting PluralNum objects.
 */
module ad.math;

import std.math;
import std.traits;

import std.stdio;

public import ad.core;


/**
 * This function computes the absolute value of the argument. It is analogous to
 * std.math.abs().
 * 
 * The derivative of the abs(x) is x*x'/abs(x).
 * 
 * Params:
 *   x = the argument
 */
@safe
export pure nothrow abs(in real x)
body {
	return std.math.abs(x);
}
@safe
export pure nothrow PluralNum!O abs(ulong O)(in PluralNum!O x)
body {
	alias PN = PluralNum!O;
	const dx = x == 0 
		? PN.DerivType!().nan 
		: x.d * sgn(x.reduce()); 
	return PN(abs(x.val), dx);
}
unittest {
	assert(abs(PluralNum!()()).same(PluralNum!()(real.nan, real.nan)));
	
	assert(abs(PluralNum!().zero).same(derivSeq(0.0L, real.nan)));
	assert(abs(derivSeq(2.0L, 2.0L)).same(derivSeq(2.0L, 2.0L)));
	assert(abs(PluralNum!().var(-1)).same(derivSeq(1.0L, -1.0L)));
	
	assert(abs(PluralNum!().infinity).same(PluralNum!().infinity));
	assert(abs(-PluralNum!().infinity).same(PluralNum!().infinity));
}


/** 
 * This function computes the cosine of the argument. It is analogous to 
 * std.math.cos().
 * 
 * Params:
 *   x = the argument
 */
@safe 
export pure nothrow real cos(in real x)
body {
	return std.math.cos(x);
}
@safe
export pure nothrow PluralNum!O cos(ulong O)(in PluralNum!O x)
body {
	return PluralNum!O(
		cos(x.val), 
		-sin(x.reduce()) * x.d);
}
unittest {
	assert(cos(PluralNum!()()).same(PluralNum!()(real.nan, real.nan)));
	
	assert(cos(PluralNum!().zero).same(derivSeq(1.0L, 0.0L)));
	assert(cos(PluralNum!().var(PI_2)).same(derivSeq(cos(PI_2), -1.0L)));
	
	assert(cos(PluralNum!().infinity).same(PluralNum!().nan));
	assert(cos(-PluralNum!().infinity).same(PluralNum!().nan));
}


/**
 * This function raises e to a given power. It is analogous to std.math.exp().
 * 
 * Params:
 *   x = the power e is raised to.
 */
@safe
export pure nothrow real exp(in real x)
body {
	return std.math.exp(x);
}
@safe
export pure nothrow PluralNum!O exp(ulong O)(in PluralNum!O x)
body {
	return PluralNum!O( 
	    exp(x.val), 
	    x.d * exp(x.reduce()));
}
unittest {
	assert(exp(PluralNum!().nan).same(PluralNum!().nan));
	
	assert(exp(PluralNum!().var(0)).same(PluralNum!().var(1)));
	
	assert(exp(PluralNum!().infinity).same(derivSeq(real.infinity, real.infinity)));
	assert(exp(-PluralNum!().infinity).same(PluralNum!().zero));
}


/**
 * This function determines whether its argument is a NaN.  It is analogous to 
 * std.math.isNaN().
 *  
 * Params:
 *   x = the argument
 */
@safe
export pure nothrow bool isNaN(in real x)
body {
	return std.math.isNaN(x);
}
@safe
export pure nothrow bool isNaN(ulong O)(in PluralNum!O x) 
body {
	return std.math.isNaN(x.val);
}
unittest {
	assert(isNaN(PluralNum!().nan));
}


/**
 * This function computes the natural logarithm of its argument. It is analogous
 * to std.math.log().
 * 
 * Params:
 *   x = the argument
 */
@safe
export pure nothrow real log(in real x)
body {
	return std.math.log(x);
}
@safe
export pure nothrow PluralNum!O log(ulong O)(in PluralNum!O x)
body {
	return x.log();
}


/**
 * This function computes the raises a given number to a given power. It is 
 * analogous to std.math.pow().
 * 
 * Params:
 *   x = the base
 *   y = the exponent
 */
@safe
export pure nothrow real pow(in real x, real y)
body {
	return std.math.pow(x, y);
}
@safe
export pure nothrow PluralNum!O pow(ulong O)(in PluralNum!O x, in real y)
body {
	return x ^^ y;
}
@safe
export pure nothrow PluralNum!O pow(ulong O)(in real x, in PluralNum!O y)
body {
	return x ^^ y;
}
@safe
export pure nothrow PluralNum!O pow(ulong O)(in PluralNum!O x, in PluralNum!O y)
body {
	return x ^^ y;
}


/**
 * This function computes the sign of the argument. It is analogous to 
 * std.math.sgn().
 * 
 * The derivate of sgn(x) evaluated at 0 is undefined or NaN.  Otherwise it is 0
 * * x'.
 * 
 * Params:
 *   x = the argument
 */
@safe
export pure nothrow real sgn(in real x)
body { 
	return std.math.sgn(x); 
}
@safe
export pure nothrow PluralNum!O sgn(ulong O)(in PluralNum!O x)
body {
	alias PN = PluralNum!O;
	if (x == 0) return PN(0, PN.DerivType!().nan);
	return PN(sgn(x.val), 0 * x.d);
}
unittest {
	assert(sgn(PluralNum!()()).same(PluralNum!()(real.nan, real.nan)));

	assert(sgn(PluralNum!()(0, real.nan)).same(PluralNum!()(0, real.nan)));
	assert(sgn(derivSeq(1.0L, 2.0L, real.nan)).same(derivSeq(1.0L, 0.0L, real.nan)));
	assert(sgn(PluralNum!()(real.infinity, 0)).same(PluralNum!().one));
	assert(sgn(PluralNum!()(-real.infinity, 0)).same(-PluralNum!().one));
}


/**
 * This function computes the square root of its argument. It is analogous to 
 * std.math.sqrt().
 * 
 * Params:
 *   x = the argument
 */
@safe
export pure nothrow real sqrt(in real x) 
body {
	return std.math.sqrt(x);
}
@safe
export pure nothrow PluralNum!O sqrt(ulong O)(in PluralNum!O x)
body {
	alias PN = PluralNum!O;
	if (x == 0) return PN(0, PN.DerivType!().nan);
	return PN( 
		sqrt(x.val), 
		x.d / (2 * sqrt(x.reduce())));		
}
unittest {
	assert(sqrt(PluralNum!().var(1)).same(derivSeq(1.0L, 0.5L)));

	assert(sqrt(PluralNum!().var(0)).same(derivSeq(0.0L, real.nan)));
	assert(sqrt(PluralNum!().var(-1)).same(PluralNum!().nan));

	assert(sqrt(PluralNum!().infinity).same(derivSeq(real.infinity, 0.0L)));

	assert(sqrt(PluralNum!().nan).same(PluralNum!().nan));
}


/+
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
+/