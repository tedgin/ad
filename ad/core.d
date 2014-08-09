/** 
 * This module implements automatic differentiation using forward accumulation and operator overloading.
 * 
 * It can only differentiate functions of the form f:R->R. This is a completely unoptimized version.
 */
module ad.core;

import std.math;
import std.string;
import std.traits;


/**
 * This class implements a generalization of the dual number concept. 
 * 
 * This class provides the basic algebraic operations of this generalization. It overloads all operators that make sense
 * for fields.
 *
 * Params:
 *  Order = the number of derivatives represented 
 *  Field = the underlying type of real number
 */
export struct PluralNum(ulong Order = 1, Field = real) if (Order >= 1 && isFloatingPoint!Field && isScalarType!Field) {

	/**
	 * This is the type constructor of the derivatives of the plural number.
	 * 
	 * Params:
	 *  DerivOrder = the order of the derivative. This must be at least 1 and no more than Order, the order of the 
	 *    plural number.
	 */
	export template DerivType(ulong DerivOrder = 1) if (DerivOrder >= 1) {
		static if (DerivOrder == Order) {
			export alias Field DerivType;
		} else static if (1 <= DerivOrder && DerivOrder < Order) {
			export alias PluralNum!(Order - DerivOrder, Field) DerivType;
		}
	}
	unittest {
		assert (is(PluralNum!(2).DerivType!() == PluralNum!(1)));
		assert (is(PluralNum!(2).DerivType!(2) == real));
	}

	private Field _x;
	private DerivType!() _dx;

	/**
	 * A constant zero represented as a plural number.
	 */
	export static immutable PluralNum zero = param(0);
	
	/**
	 * A constant one represented as a plural number.
	 */
	export static immutable PluralNum one = param(1);
	
	export static immutable PluralNum nan = PluralNum(Field.nan, DerivType!().nan);

	export static immutable PluralNum infinity = PluralNum(Field.infinity, one.reduce());
	
	export static immutable ulong dig = Field.dig;

	export static immutable PluralNum epsilon = param(Field.epsilon);

	export static immutable ulong mant_dig = Field.mant_dig;

	export static immutable int max_10_exp = Field.max_10_exp;
	
	export static immutable int max_exp = Field.max_exp;

	export static immutable int min_10_exp = Field.min_10_exp;

	export static immutable int min_exp = Field.min_exp;

	export static immutable PluralNum max = param(Field.max);

	export static immutable PluralNum min_normal = param(Field.min_normal);

	/**
	 * This function constructs a plural number representing a constant with respect to the derivative.
	 *
	 * Params:
	 *  val = The value of the parameter or constant.
	 *
	 * Returns:
	 *  The plural number representation of the parameter or constant.
	 */
	@safe 
	export static pure nothrow PluralNum param(in Field val)
	body {
		static if (isScalarType!(DerivType!())) return PluralNum(val, 0);
		else                                    return PluralNum(val, DerivType!().zero);
	}
	unittest {
		const q = PluralNum!(1).param(1);
		assert(q._x is 1.0);
		assert(q._dx is 0.0);

		const w = PluralNum!(2).param(1);
		assert(w._x is 1.0);
		assert(w._dx == PluralNum!(1).zero);
	}

	/**
	 * This functions constructs a plural number representing the variable of differentiation.
	 *
	 * Params:
	 *  val = The value of the variable.
	 *
	 * Returns:
	 *  The plural number representation of the variable.
	 */
	@safe 
	export static pure nothrow PluralNum var(in Field val)
	body {
		static if (isScalarType!(DerivType!())) return PluralNum(val, 1);
		else                                    return PluralNum(val, DerivType!().one);
	}
	unittest {
		const q = PluralNum!(1).var(2);
		assert(q._x is 2.0);
		assert(q._dx is 1.0);
		
		const w = PluralNum!(2).var(2);
		assert(w._x is 2.0);
		assert(w._dx == PluralNum!(1).one);
	}

	/**
	 * Constructs a plural number from if value and derivatives
	 * 
	 * Params:
	 *   val = the value
	 *   derivs = the derivatives
	 */
	package pure nothrow this(in Field val, in DerivType!() derivs) 
	body {
		_x = val;
		_dx = derivs;
	}

	/**
	 * Constructs a plural number from its value and the values of its derivatives
	 * 
	 * Params:
	 *   derivVals = the an array of the derivative values where index is the order of the derivative
	 */
	package pure nothrow this(in Field[Order + 1] derivVals ...)
	body {
		_x = derivVals[0];
		static if (Order > 1) _dx = DerivType!()(derivVals[1 .. Order + 1]);
		else                  _dx = derivVals[1];
	}
	unittest {
		assert (PluralNum!()([1.0, 2.0]) == PluralNum!()(1.0, 2.0));
		assert (PluralNum!(2)(1.0L, 2.0L, 3.0L) == PluralNum!(2)(1.0, PluralNum!()(2.0, 3.0)));
	}


	@safe @property 
	export pure nothrow const PluralNum!(Order, typeof(Field.nan.re)) re() 
	body {
		return PluralNum!(Order, typeof(Field.nan.re))(_x.re, _dx.re);
	}
	
	@safe @property 
	export pure nothrow const PluralNum!(Order, typeof(Field.nan.im)) im() 
	body {
		return PluralNum!(Order, typeof(Field.nan.im))(_x.im, _dx.im);
	}
	
	@safe 
	export pure nothrow const Field opCast(Field)() 
	body {
		return _x;
	}
	unittest {
		const q = PluralNum!()(2, 1);
		assert(cast(real)(q) is 2.0);
	}
	
	/**
	 * Returns the value represented by the plural number.
	 * 
	 * Returns:
	 *  The value represented by the plural number
	 */
	@safe @property 
	export pure nothrow const Field val()
	body {
		return _x;
	}

	/**
	 * Returns the derivative of order DerivOrder of the plural number. DerivOrder must be at least one but no more than 
	 * Order. If DerivOrder == Order, the derivative will be a value of the type of the underlying field. Otherwise the 
	 * derivative will be a plural number but with order Order - DeriveOrder. If no value is provided for DerivOrder, a 
	 * value of 1 is assume.  
	 * 
	 * Params:
	 *  DerivOrder = The order of the derivate to compute.
	 * 
	 * Returns:
	 *  The derivative of the plural number.
	 */
	@safe @property 
	export pure nothrow const DerivType!(DerivOrder) d(ulong DerivOrder = 1)() 
	if (1 <= DerivOrder && DerivOrder <= Order)
	body {
		static if (DerivOrder == 1) return _dx;
		else                        return _dx.d!(DerivOrder - 1)();
	}
	unittest {
		const q = PluralNum!(3)(3, PluralNum!(2)(2, PluralNum!()(1, 0)));
		const dq = q.d!()();
		assert (is(typeof(dq) == const(PluralNum!(2))));
		assert (dq.val() == 2);
		
		const d2q = q.d!(2)();
		assert (is(typeof(d2q) == const(PluralNum!(1))));
		assert (d2q.val() == 1);

		const d3q = q.d!(3)();
		assert (is(typeof(d3q) == const(real)));
		assert (d3q == 0);
	}
	
	@safe 
	export pure nothrow const bool opEquals(in PluralNum that) 
	body {
		return this._x == that._x;
	}
	
	@safe
	export pure nothrow const bool opEquals(in Field val) 
	body{
		return _x == val;
	}
	
	@safe
	export pure nothrow const int opCmp(in PluralNum that) 
	body {
		return opCmp(that._x);
	}

	@safe 
	export pure nothrow const int opCmp(in Field val) 
	body {
		if (_x < val) return -1;
		if (_x > val) return 1;
		else          return 0;
	}
	unittest {
		const q = PluralNum!()(2, 1);
		assert (q < 3);
		assert (q !<> 2);
		assert (q >= 1);
	}

	@trusted
	export pure nothrow const hash_t toHash()
	body {
		return cast(hash_t)(_x);
	}

	@safe
	export pure nothrow const PluralNum opUnary(string op)()
	body {
		     static if (op == "+") return PluralNum(+_x, +_dx);
		else static if (op == "-") return PluralNum(-_x, -_dx);
		else static assert(false, "Operator " ~ op ~ " not implemented");
	}
	unittest {
		const q = PluralNum!()(2, 1);
		const w = +q;
		assert (w._x == 2);
		assert (w._dx == 1);

		const e = -q;
		assert (e._x == -2);
		assert (e._dx == -1);
	}

	/**
	 * Computes the inverse of the plural number. That is, given a plural number x, it computes the plural number y 
	 * where x * y is identically 1.
	 * 
	 * Returns:
	 *  It returns the inverted plural number.
	 */
	@safe
	export pure nothrow const PluralNum inv()
	body {
		const reducedX = reduce();
		return PluralNum(1 / _x, -_dx / (reducedX * reducedX));
	}
	unittest {
		import std.stdio;
		const x = PluralNum!(3)(1, PluralNum!(2)(2, PluralNum!(1)(3, 4)));
		const y = x.inv();
		assert (x * y == PluralNum!(3).one);
	}

	@safe
	export pure nothrow const PluralNum opBinaryRight(string op)(in Field val)
	body {
		     static if (op == "+")  return param(val) + this;
		else static if (op == "-")  return param(val) - this;
		else static if (op == "*")  return param(val) * this;
		else static if (op == "/")  return param(val) / this;
		else static if (op == "%")  return param(val) % this;
		else static if (op == "^^") return param(val) ^^ this;
		else static assert(false, "Operator " ~ op ~ " not implemented");
	}

	@safe
	export pure nothrow const PluralNum opBinary(string op)(in Field val)
	body {
		     static if (op == "+")  return this + param(val);
		else static if (op == "-")  return this - param(val);
		else static if (op == "*")  return this * param(val);
		else static if (op == "/")  return this / param(val);
		else static if (op == "%")  return this % param(val);
		else static if (op == "^^") return this ^^ param(val);
		else static assert(false, "Operator " ~ op ~ " not implemented");
	}

	@safe
	export pure nothrow const PluralNum opBinary(string op)(in PluralNum that) if (op == "+")
	body {
		return PluralNum(this._x + that._x, this._dx + that._dx);
	}
	unittest {
		const q = PluralNum!()(2, 1);
		const w = PluralNum!()(4, 3);
		const e = q + w;
		assert (e._x == 6 && e._dx == 4);

		const r = q + 1.0;
		assert (r._x == 3 && r._dx == 1);

		const t = 2 + q;
		assert (t._x == 4 && t._dx == 1);
	}

	@safe
	export pure nothrow const PluralNum opBinary(string op)(in PluralNum that) if (op == "-")
	body {
		return this + -that;
	} 

	@safe
	export pure nothrow const PluralNum opBinary(string op)(in PluralNum that) if (op == "*") 
	body {
		return PluralNum(this._x * that._x, this._dx * that.reduce() + this.reduce() * that._dx);
	} 
	unittest {
		const q = PluralNum!()(1, 2);
		const w = PluralNum!()(3, 5);
		const e = q * w;
		assert (e._x == 3 && e._dx == 11);
	}

	@safe
	export pure nothrow const PluralNum opBinary(string op)(in PluralNum that) if (op == "/") 
	body {
		return this * that.inv();
	} 

	@safe
	export pure nothrow const PluralNum opBinary(string op)(in PluralNum that) if (op == "%") 
	body {
		const thisModThat = this._x % that._x;
		const reduceThis = this.reduce();
		const reduceThat = that.reduce();
		const dThat_that = that._dx / reduceThat;
		const dx_term1 = (reduceThis % reduceThat) * dThat_that;
		const dx_term2a = this._dx - reduceThis * dThat_that;
		const dx_term2b = (thisModThat == 0 && this._x / that._x != 0) ? -infinity : one;
		return PluralNum!()(thisModThat, dx_term1 + dx_term2a * dx_term2b.reduce());
	} 
	unittest {
		const e = PluralNum!()(3, 1) % PluralNum!()(2, 4);
		assert (e._x == 1 && e._dx == -3);

		const q = PluralNum!()(-3, 1) % PluralNum!()(2, 4);
		assert (q._x == -1 && q._dx == 5);

		const r = PluralNum!()(3, 1) % PluralNum!()(-2, 4);
		assert (r._x == 1 && r._dx ==  5);

		const t = PluralNum!()(-3, 1) % PluralNum!()(-2, 4);
		assert (t._x == -1 && t._dx == -3);

		const w = PluralNum!()(4, 1) % PluralNum!()(2, 4);
		assert (w._x == 0 && w._dx == real.infinity);

		const y = PluralNum!()(0, 1) % PluralNum!()(2, 4);
		assert (y._x == 0 && y._dx == 1);

		const u = PluralNum!()(0, 1) % PluralNum!()(0, 4);
		assert (isnan(u._x) && isnan(u._dx));

		const i = PluralNum!()(real.infinity, 1) % PluralNum!()(2, 4);
		assert (isnan(i._x) && isnan(i._dx));

		const o = PluralNum!()(3, 1) % PluralNum!()(real.infinity, 1);
		assert (o._x == 3 && o._dx == 1);

		const p = PluralNum!()(real.nan, real.nan) % PluralNum!()(2, 4);
		assert (isnan(p._x) && isnan(p._dx));

		const a = PluralNum!()(3, 1) % PluralNum!()(real.nan, real.nan);
		assert (isnan(a._x) && isnan(a._dx));
	}

	@safe
	export pure nothrow const PluralNum opBinary(string op)(in PluralNum that) if (op == "^^") 
	body {
		const f = this.reduce();
		const fp = this._dx;
		const g = that.reduce();
		const gp = that._dx;
		const fug = f ^^ g;
		return PluralNum(this._x ^^ that._x, gp * fug * f.log() + fp * fug * g / f);
	}
	unittest {
		const q = PluralNum!()(2, 1) ^^ PluralNum!()(3, 4);
		assert(q._x == 8 && q._dx == (32*std.math.log(2) + 12));

		const w = PluralNum!()(0, 1) ^^ PluralNum!()(3, 4);
		assert(w._x == 0 && isnan(w._dx));

		const e = PluralNum!()(-1,1) ^^ PluralNum!()(3, 4);
		assert(e._x == -1 && isnan(e._dx));
	}
	
	const string toString()
	body {
		return format(0);
	}
	unittest {
		const q = PluralNum!(2)(1, PluralNum!(1)(2, 3));
		assert (q.toString() == "1.000000 + 2.000000d + 3.000000d2");
	}
	/**
	 * Computes the natural logarithm of the plural number. The logarithm is a method attached to this struct because it
	 * is required to compute the derivative of the ^^ operator.
	 * 
	 * The derivative of the logarithm is undefined when operand is non-positive.
	 * 
	 * returns:
	 *  A new plural number that is the natural logarithm of this dual number.
	 */
	@safe
	package pure nothrow const PluralNum log()
	body {
		const dlog = _x > 0 ? _dx / reduce() : nan.reduce();
		return PluralNum(std.math.log(_x), dlog);
	}
	unittest {
		const p = PluralNum!(2).var(1).log();
		assert(p._x == 0);
		assert(p._dx._x == 1);
		assert(p._dx._dx == -1);

		const z = PluralNum!(1).var(0).log();
		assert(z._x == -real.infinity);
		assert(isNaN(z._dx));

		const n = PluralNum!(1).var(-1).log();
		assert(isNaN(n._x));
		assert(isNaN(n._dx));
	}

	/**
	 * Reduces the order of the plural number by one. If the plural number is first order, the value returned will be a
	 * field. This truncates the highest order derivative.
	 * 
	 * returns:
	 *  The plural number of order one less.
	 */
	@safe
	package pure nothrow const DerivType!() reduce()
	body {
		static if(Order == 1) return _x;
		else                  return DerivType!()(_x, _dx.reduce());
	}
	unittest {
		const p2 = PluralNum!(2)(2, PluralNum!(1)(3, 4));
		const p1 = p2.reduce();
		assert(p2._x == p1._x);
		assert(p2._dx._x == p1._dx);

		const p0 = p1.reduce();
		assert(is(typeof(p0) == const(real)));
		assert(p2._x == p0);
	}
	
	private const string format(in ulong derivOrder)
	body {
		static if (isScalarType!(DerivType!()))
			const tail = std.string.format("%f%s", _dx, formatSuffix(derivOrder + 1));
		else                               
			const tail = _dx.format(derivOrder + 1);
		return std.string.format("%f%s + %s", _x, formatSuffix(derivOrder), tail);
	}
	unittest {
		assert(PluralNum!(1)(2.0, 3.0).format(1) == "2.000000d + 3.000000d2");

		const q = PluralNum!(2)(2.0, PluralNum!(1)(3.0, 4.0));
		assert(q.format(1) == "2.000000d + 3.000000d2 + 4.000000d3");
	}
}
unittest {
	// force the unit tests
	const PluralNum!() w;

	assert (isnan(w._x));
	assert (isnan(PluralNum!().init._x));
}


/**
 * Constructs a plural number from its value and the values of its derivatives. If there is only one value is the 
 * sequence, that value is returned.
 * 
 * Params:
 *   Len = the number of elements in the sequence (optional)
 *   Field = the underlying type of the real number (optional)
 *   derivVals = the an array of the derivative values where index is the order of the derivative
 */
@safe
export pure nothrow Field derivSeq(ulong Len, Field)(in Field[Len] val ...) if (Len == 1)
body {
	return val[0];
}
@safe
export pure nothrow PluralNum!(Len - 1, Field) derivSeq(ulong Len, Field)(in Field[Len] derivVals ...) if (Len > 1)
body {
	return PluralNum!(Len - 1, Field)(derivVals);
}
unittest {
	assert (derivSeq(1.0) == 1.0);
	assert (derivSeq(1.0L, 2.0L) == PluralNum!()(1.0, 2.0));
}


private string formatSuffix(in ulong derivOrder)
body {
	switch (derivOrder) {
		case 0:  return "";
		case 1:  return "d";
		default: return format("d%d", derivOrder);
	}
}
unittest {
	assert ("" == formatSuffix(0));
	assert ("d" == formatSuffix(1));
	assert ("d2" == formatSuffix(2));
	assert ("d10" == formatSuffix(10));
}
