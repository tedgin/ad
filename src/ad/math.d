module ad.math;

import std.math;
import std.stdio;

public import ad.core;


export pure DerivSeqType abs(DerivSeqType)(const DerivSeqType u)
in {
	static if (DerivSeqType.Order > 0) {
		assert (u.val() != 0);
	}
}
body {
	static if (DerivSeqType.Order == 0) {
		return DerivSeqType(std.math.abs(u.val()));
	} else {
		return DerivSeqType(std.math.abs(u.val()), u.reduce().sgn());
	}
}


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


export pure DerivSeqType cos(DerivSeqType)(const DerivSeqType u)
body {
	static if (DerivSeqType.Order == 0) {
		return DerivSeqType(std.math.cos(u.val()));
	} else {
		return DerivSeqType(
			std.math.cos(u.val()), 
			-u.d() * sin(u.reduce()));
	}
}


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
