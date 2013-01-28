/** 
 * This module implements automatic differentiation using forward accumulation and operator 
 * overloading.
 * 
 * TODO figure out how to describe this.
 *
 * It can only differentiate functions of the form f:R->R. This is a completely unoptimized version.
 */
module ad.core;

/* TODO document
 * TODO start using NaNs and infs
 * TODO implement % and ^^
 * TODO start using refs for struct params
 */

import std.math;
import std.string;


/**
 * This class implements a generalization of the dual number concept. 
 * 
 * This class provides the basic algebraic operations of this generalization.
 *
 * Params:
 *  Order = the number of derivatives represented 
 *  ValType = the underlying type of real number
 */
export struct DerivSeq(uint Order = 1, ValType = real) 
if (Order > 0) {
	private ValType _x;
	private TailSeq _dx;

	/**
	 * This type represents the sequence of the derivatives following the head of this sequence.
	 */
	export alias DerivSeq!(Order - 1, ValType) TailSeq;
	
	/**
	 * The number of derivatives represented.
	 */
	// TODO learn why this is better than an immutable static
	export enum { Order = Order }

	/**
	 * A constant zero represented as a DerivSeq.
	 */
	export static immutable DerivSeq zero = param(0);
	
	/**
	 * A constant one represented as a DerivSeq.
	 */
	export static immutable DerivSeq one = param(1);
	
	/**
	 * positive infinity
	 */
	export static immutable DerivSeq infinity = DerivSeq(ValType.infinity, TailSeq.nan);
	
	/**
	 * not a number
	 */
	export static immutable DerivSeq nan = DerivSeq(ValType.nan, TailSeq.nan);
	
	/**
	 * This function constructs a dual number representing a constant with respect to the 
	 * derivative.
	 *
	 * Params:
	 *	val = The value of the parameter or constant.
	 *
	 * Returns:
	 *	The dual number representation of the parameter or constant.
	 */
	export static pure DerivSeq param(const ValType val)
	body {
		return DerivSeq(val, TailSeq.zero);
	}
	
	/**
	 * This functions constructs a dual number representing the variable of differentiation.
	 *
	 * Params:
	 *  val = The value of the variable.
	 *
	 * Returns:
	 *  The dual number representation of the variable.
	 */
	export static pure DerivSeq var(const ValType val)
	body {
		return DerivSeq(val, TailSeq.one);
	}

	/**
	 * TODO document
	 */
	export pure this(const ValType val, const TailSeq derivs) 
	body {
		_x = val;
		_dx = derivs;
	}
	
	/**
	 * TODO document
	 */
	export pure const ValType opCast(ValType)() 
	body {
		return _x;
	}
	
	/**
	 * Returns the value represented by the dual number.
	 * 
	 * Returns:
	 * 	The value represented by the dual number
	 */
	export pure const ValType val()
	body {
		return _x;
	}

	/**
	 * Returns the derivative of order DerivOrder of the dual number. The derivative will be a dual 
	 * number but with order Order - DeriveOrder. DerivOrder must be at least one but no more than 
	 * Order. If no value is provided for DerivOrder, a value of 1 is assume.  
	 * 
	 * Returns:
	 * 	The derivative of the dual number as a dual number of reduced order.
	 */
	export pure const DerivSeq!(Order - DerivOrder, ValType) d(uint DerivOrder = 1)() 
	if (1 <= DerivOrder && DerivOrder <= Order)
	body {
		static if (DerivOrder == 1) return _dx;
		else                        return _dx.d!(DerivOrder - 1)();
	}
	
	export pure const bool opEquals(const DerivSeq that) 
	body {
		return this._x == that._x;
	}
	
	export pure const bool opEquals(const ValType val) 
	body{
		return _x == val;
	}
	
	export pure const int opCmp(ref const DerivSeq that) 
	body {
		return opCmp(that._x);
	}
	
	export pure const int opCmp(const ValType val) 
	body {
		if (_x < val) return -1;
		if (_x > val) return 1;
		else          return 0;
	}
	
	export pure const DerivSeq opUnary(string op)()
	body {
		     static if (op == "+") return DerivSeq(+_x, +_dx);
		else static if (op == "-") return DerivSeq(-_x, -_dx);
		else static assert(false, "Operator " ~ op ~ " not implemented");
	}

	export pure const DerivSeq inv()
	in { assert(_x != 0); }
	body {
		const reducedX = reduce();
		return DerivSeq(1 / _x, -_dx / (reducedX * reducedX));
	}

	export pure const DerivSeq sgn()
	in { assert(_x != 0); }
	body {
		return DerivSeq(_x < 0 ? -1 : 1, TailSeq.zero);
	}

	export pure const DerivSeq opBinary(string op)(const DerivSeq that)
	body {
		     static if (op == "+") return DerivSeq(this._x + that._x, this._dx + that._dx);
		else static if (op == "-") return DerivSeq(this._x - that._x, this._dx - that._dx);
		else static if (op == "*") 
			return DerivSeq(this._x * that._x, this._dx * that.reduce() + this.reduce() * that._dx);
		else static if (op == "/") return this * that.inv();
		else static if (op == "%") {
			const thisModThat = this._x % that._x;
			const reduceThis = this.reduce();
			const reduceThat = that.reduce();
			const dThat_that = that._dx / reduceThat;
			const dx_term1 = (reduceThis % reduceThat) * dThat_that;
			const dx_term2a = this._dx - reduceThis * dThat_that;
			const dx_term2b = (thisModThat == 0 && this._x / that._x != 0) ? -Inf : One;  
			return DerivSeq(thisModThat, dx_term1 + dx_term2a * dx_term2b);
		} else static if (op == "^^") { 
			const f = this.reduce();
			const fp = this._dx;
			const g = that.reduce();
			const gp = that._dx;
			return DerivSeq(thix._x ^^ that._x, f ^^ g * (gp * f.log() + fp * g / f));
		} else static assert(false, "Operator " ~ op ~ " not implemented");
	}
	
	export pure const DerivSeq opBinary(string op)(const ValType val)
	body {
		     static if (op == "+")  return this + param(val);
		else static if (op == "-")  return this - param(val);
		else static if (op == "*")  return this * param(val);
		else static if (op == "/")  return this / param(val);
		else static if (op == "%")  return this % param(val);
		else static if (op == "^^") return this ^^ param(val);
		else static assert(false, "Operator " ~ op ~ " not implemented");
	}

	export pure const DerivSeq opBinaryRight(string op)(const ValType val)
	body {
		     static if (op == "+")  return param(val) + this;
		else static if (op == "-")  return param(val) - this;
		else static if (op == "*")  return param(val) * this;
		else static if (op == "/")  return param(val) / this;
		else static if (op == "%")  return param(val) % this;
		else static if (op == "^^") return param(val) ^^ this;
		else static assert(false, "Operator " ~ op ~ " not implemented");
	}

	export pure const TailSeq reduce()
	body {
		static if(Order == 1) return TailSeq(_x);
		else                  return TailSeq(_x, _dx.reduce());
	}

	export const string toString()
	body {
		return format(0);
	}
	
	package pure const DerivSeq log()
	in {
		assert(_x > 0);
	}
	body {
		return DerivSeq(std.math.log(_x), _dx / reduce());
	}
	
	// TODO learn why std.string.format is not pure
	private const string format(const uint power)
	body {
		return std.string.format("%f%s + %s", _x, formatSuffix(power), _dx.format(power + 1));
	}
}


export struct DerivSeq(uint Order : 0, ValType = real) {		
	private ValType _x;
	
	export enum { Order = Order }

	export static immutable DerivSeq zero = param(0);
	
	export static immutable DerivSeq one = param(1);
	
	export static immutable DerivSeq infinity = DerivSeq(ValType.infinity);
	
	export static immutable DerivSeq nan = DerivSeq(ValType.nan);
	
	export static pure DerivSeq param(const ValType val)
	body {
		return DerivSeq(val);
	}

	export static pure DerivSeq var(const ValType val)
	body {
		return DerivSeq(val);
	}
	
	/**
	 * TODO document
	 */
	export pure const ValType opCast(ValType)() 
	body {
		return _x;
	}
	
	export pure const ValType val()
	body {
		return _x;
	}

	export pure const bool opEquals(const DerivSeq that) 
	body {
		return this._x == that._x;
	}
	
	export pure const bool opEquals(const ValType val) 
	body {
		return _x == val;
	}
	
	export pure const int opCmp(const DerivSeq that) 
	body {
		return opCmp(that._x);
	}
	
	export pure const int opCmp(const ValType val) 
	body {
		if (_x < val) return -1;
		if (_x > val) return 1;
		else          return 0;
	}
	
	export pure const DerivSeq opUnary(string op)()
	body {
		     static if (op == "+") return DerivSeq(+_x);
		else static if (op == "-") return DerivSeq(-_x);
		else static assert(false, "Operator " ~ op ~ " not implemented");
	}

	export pure const DerivSeq inv()
	in { assert(_x != 0); }
	body {
		return DerivSeq(1 / _x);
	}
	
	export pure const DerivSeq sgn()
	in { assert(_x != 0); }
	body {
		return DerivSeq(_x < 0 ? -1 : 1);
	}

	export pure const DerivSeq opBinary(string op)(const DerivSeq that)
	body {
		     static if (op == "+")  return DerivSeq(this._x + that._x);
		else static if (op == "-")  return DerivSeq(this._x - that._x);
		else static if (op == "*")  return DerivSeq(this._x * that._x);
		else static if (op == "/")  return DerivSeq(this._x / that._x);
		else static if (op == "%")  return DerivSeq(this._x % that._x);
		else static if (op == "^^") return DerivSeq(this._x ^^ that._x);
		else static assert(false, "Operator " ~ op ~ " not implemented");
	}
	
	export pure const DerivSeq opBinary(string op)(const ValType val)
	body {
		     static if (op == "+")  return DerivSeq(_x + val);
		else static if (op == "-")  return DerivSeq(_x - val);
		else static if (op == "*")  return DerivSeq(_x * val);
		else static if (op == "/")  return DerivSeq(_x / val);
		else static if (op == "%")  return DerivSeq(_x % val);
		else static if (op == "^^") return DerivSeq(_x ^^ val);
		else static assert(false, "Operator " ~ op ~ " not implemented");
	}

	export pure const DerivSeq opBinaryRight(string op)(const ValType val)
	body {
		     static if (op == "+")  return DerivSeq(val + _x);
		else static if (op == "-")  return DerivSeq(val - _x);
		else static if (op == "*")  return DerivSeq(val * _x);
		else static if (op == "/")  return DerivSeq(val / _x);
		else static if (op == "%")  return DerivSeq(val / _x);
		else static if (op == "^^") return DerivSeq(val ^^ _x);
		else static assert(false, "Operator " ~ op ~ " not implemented");
	}

	export const string toString()
	body {
		return format(0);
	}

	package pure const DerivSeq log()
	in {
		assert(_x > 0);
	}
	body {
		return DerivSeq(std.math.log(_x));
	}
	
	private const string format(const uint power) 
	body {
		return std.string.format("%f%s", _x, formatSuffix(power));
	}
}


private string formatSuffix(const uint power)
body {
		switch (power) {
		case 0:  return "";
		case 1:  return "d";
		default: return format("d%d", power);
		}
}
