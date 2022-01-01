/** 
This module implements automatic differentiation using forward accumulation and
operator overloading. It can only differentiate functions of the form *f:R->R*. 
This is a completely unoptimized version. The traditional definition of 
differentiation is used, not the generalized notion from distribution theory.
*/
module ad.core;

import std.format : format;
import std.math : isNaN, log;
import std.traits : isImplicitlyConvertible;

/**
This data structure implements a generalization of the dual number concept. It 
provides the basic algebraic operations of this generalization using operator 
overloading.

Params:
    Degree = the number of derivatives represented
*/
struct GenDualNum(ulong Degree = 1)
{
    /**
    This is the type constructor of the derivatives.
    
    Params:
        Order = the order of the derivative. This must be at least one but no 
            more  than Degree, the degree of the generalized dual number.
    */
    template DerivType(ulong Order = 1)
    {
        static assert(
            Order <= Degree, 
            "The order of derivative cannot be larger than the degree of the generalized dual" ~
            " number.");

        static if (Order < Degree) 
            alias DerivType = GenDualNum!(Degree - Order);
        else
            alias DerivType = real;
    }

    private static DerivType!1 mkNaNDeriv()
    {
        static if (Degree == 1)
            return real.nan;
        else
            return DerivType!1(real.nan, DerivType!1.mkNaNDeriv());
    }

    private static DerivType!1 mkZeroDeriv()
    {
        static if (Degree == 1)
            return 0;
        else
            return DerivType!1(0, DerivType!1.mkZeroDeriv());
    }

    // FIELDS

    private real _x;
    private DerivType!1 _dx;

    // CONSTRUCTORS

    package this(in real val, in DerivType!1 derivs)
    {
        _x = val;
        _dx = derivs;
    }

    private this(in real[Degree + 1] derivVals...)
    {
        _x = derivVals[0];

        static if (Degree == 1)
            _dx = derivVals[1];
        else
            _dx = DerivType!1(derivVals[1 .. Degree + 1]);
    }

    /**
    This function constructs a generalized dual number representing a parameter 
    or constant with respect to differentiation.
        
    Params:
        val = The value of the parameter or constant.
    
    Returns:
        The representation of the parameter or constant
    */
    static nothrow pure @safe GenDualNum param(in real val)
    {
        return GenDualNum(val, mkZeroDeriv());
    }

    /**
    This functions constructs a generalized dual number representing the 
    variable of differentiation.
    
    Params:
        val = The value of the variable.
        
    Returns:
        The representation of the variable.
    */
    static nothrow pure @safe GenDualNum var(in real val)
    {
        static if (Degree == 1)
            return GenDualNum(val, 1);
        else
            return GenDualNum(val, DerivType!1.one);
    }

    // TODO consider adding the following constructor for converting a 
    //   parameter to a variable.
    // static nothrow pure @safe GenDualNum var(in GenDualNum param)
    // {
    //     if (param.d!1.same(DerivType!1.zero) {
    //         return var(param.val);
    //     } else {
    //         return param;
    // }

    // CONSTANTS

    /// This is a constant _zero represented as a generalized dual number.
    static immutable GenDualNum zero = GenDualNum(0, mkZeroDeriv());

    /// This is a constant _one represented as a generalized dual number.
    static immutable GenDualNum one = GenDualNum(1, mkZeroDeriv());

    // PROPERTIES OF REAL

    /// This is a generalized dual number representing infinity.
    static immutable GenDualNum infinity = GenDualNum(real.infinity, mkZeroDeriv());

    /// This is a generalized dual number representing `NaN`.
    static immutable GenDualNum nan = GenDualNum(real.nan, mkNaNDeriv());

    /// This is the number of decimal digits of precision of the value.
    static immutable int dig = real.dig;

    /** 
    This is the smallest generalized dual number such that 
    `one + epsilon > one`.
    */
    static immutable GenDualNum epsilon = GenDualNum(real.epsilon, mkZeroDeriv());

    /// This is the number of bits in the mantissa of the value.
    static immutable int mant_dig = real.mant_dig;

    /**
    This is the maximum `int` value such that `10 ^^ max_10_exp` is 
    representable as a generalized dual number.
    */    
    static immutable int max_10_exp = real.max_10_exp;

    /** 
    This is the maximum `int` value such that `2 ^^ (max_exp - 1)` is 
    representable as a generalized dual number.
    */
    static immutable int max_exp = real.max_exp;

    /** 
    This is the minimum `int` value such that `10 ^^ min_10_exp` is 
    representable as a generalized dual number. 
    */
    static immutable int min_10_exp = real.min_10_exp;

    /** 
    This is the minimum `int` value such that `2 ^^ (min_exp - 1)` is 
    representable as a generalized dual number.
    */
    static immutable int min_exp = real.min_exp;

    /// This is the largest finite generalized dual number. 
    static immutable GenDualNum max = GenDualNum(real.max, mkZeroDeriv());

    /// This is the smallest positive generalized dual number. 
    static immutable GenDualNum min_normal = GenDualNum(real.min_normal, mkZeroDeriv());

    /// This is the real part.
    @property nothrow pure @safe GenDualNum re() const
    {
        return this;
    }

    /**
    This is the imaginary part. It is always has a zero value since GenDualNum 
    only supports real numbers.
    */
    @property nothrow pure @safe GenDualNum im() const
    {
        return zero;
    }

    // MEMBERS

    /**
    This is the value of the generalized dual number.
        
    Returns:
        the value
    */
    nothrow pure @safe real val() const
    {
        return _x;
    }

    /**
    This is the derivative of order `Order` of the generalized dual number. The 
    order must be at least one but no more that the degree of the generalized 
    dual number. If `Order == Degree`, the derivative will be a `real`. 
    Otherwise it will be a `GenDualNum` but with degree `Degree - Order`. 
    
    Params:
        Order = the order of the derivate to compute, default `1`
        
    Returns:
        the derivative
    */
    nothrow pure @safe DerivType!Order d(ulong Order = 1)() const 
    if (0 < Order && Order <= Degree)
    {
        static if (Order == 1)
            return _dx;
        else
            return _dx.d!(Order - 1);
    }

    /// ditto
    nothrow pure @safe GenDualNum d(ulong Order : 0)() const
    {
        return this;
    }

    /**
    This computes the inverse of a generalized dual number. That is, given a 
    generalized dual number *x*, it computes the generalized dual number *y* 
    where *xy = 1*.
    
    Returns:
        the inverted generalized dual number
    */
    nothrow pure @safe GenDualNum inv() const
    {
        const reducedX = reduce();
        return GenDualNum(1 / _x, -_dx / (reducedX * reducedX));
    }

    /*
    The logarithm is defined in this module instead of core, because it is 
    required to compute the derivative of the ^^ operator.
    */
    package GenDualNum log() const
    {
        static import std.math;

        const dlog = _x > 0 ? _dx / reduce() : nan.reduce();
        return GenDualNum(std.math.log(_x), dlog);
    }

    package DerivType!1 reduce() const
    {
        static if (Degree == 1)
            return DerivType!1(_x);
        else
            return DerivType!1(_x, _dx.reduce());
    }

    // OPERATOR OVERLOADING

    /**
    This provides support for using the operators `==` and `!=` to compare a 
    generalized dual number to a real or another generalized dual number. Two
    generalized dual numbers are equal if their values are equal regardless of 
    the values of their derivative terms.
    */
    nothrow pure @safe bool opEquals(ulong D)(in GenDualNum!D that) const
    {
        return this._x == that._x;
    }

    /// ditto
    nothrow pure @safe bool opEquals(in real val) const
    {
        return _x == val;
    }

    /**
    This provides support for using the operators `<`, `<=`, `>=`, and `>` to 
    compare a generalized dual number to a real or another generalized dual 
    number. If *x* and *y* are two generalized dual numbers, *x < y*, if the 
    value of *x* is less than the value of *y* regardless of the values of their 
    derivative terms.
    */
    nothrow pure @safe int opCmp(ulong D)(in GenDualNum!D that) const
    {
        return opCmp(that._x);
    }

    /// ditto
    nothrow pure @safe int opCmp(in real val) const
    {
        if (_x < val)
            return -1;
        if (_x > val)
            return 1;
        else
            return 0;
    }

    /// This provides support for using the prefix operators `+` and `-`.
    nothrow pure @safe GenDualNum opUnary(string op)() const 
    {
        static if (op == "+")
            return GenDualNum(+_x, +_dx);
        else static if (op == "-")
            return GenDualNum(-_x, -_dx);
        else
            static assert(false, "Operator " ~ op ~ " not implemented");
    }

    /**
    This provides support for using the arithmetic infix operators `+`, `-`, 
    `*`, `/`, `%`, and `^^` for combining a generalized dual number with a real 
    number or another generalized dual number. 
    */
    nothrow pure @safe GenDualNum opBinary(string op : "+")(in GenDualNum that) const
    {
        return GenDualNum(this._x + that._x, this._dx + that._dx);
    }

    /// ditto
    nothrow pure @safe GenDualNum opBinary(string op : "-")(in GenDualNum that) const
    {
        return this + -that;
    }

    /// ditto
    nothrow pure @safe GenDualNum opBinary(string op : "*")(in GenDualNum that) const
    {
        return GenDualNum(this._x * that._x, this._dx * that.reduce() + this.reduce() * that._dx);
    }

    /// ditto
    nothrow pure @safe GenDualNum opBinary(string op : "/")(in GenDualNum that) const
    {
        return this * that.inv();
    }

    /// ditto
    nothrow pure @safe GenDualNum opBinary(string op : "%")(in GenDualNum that) const
    {
        const thisModThat = this._x % that._x;
        const reduceThis = this.reduce();
        const reduceThat = that.reduce();
        const dThat_that = that._dx / reduceThat;
        const dxTerm1 = (reduceThis % reduceThat) * dThat_that;
        const dxTerm2a = this._dx - reduceThis * dThat_that;
        const dxTerm2b = (thisModThat == 0 && this._x / that._x != 0) ? -infinity : one;
        return GenDualNum(thisModThat, dxTerm1 + dxTerm2a * dxTerm2b.reduce());
    }

    /// ditto
    nothrow pure @safe GenDualNum opBinary(string op : "^^")(in GenDualNum that) const
    {
        const f = this.reduce();
        const fp = this._dx;
        const g = that.reduce();
        const gp = that._dx;
        const fug = f ^^ g;
        return GenDualNum(this._x ^^ that._x, gp * fug * f.log() + fp * fug * g / f);
    }

    /// ditto
    nothrow pure @safe GenDualNum opBinary(string op)(in real val) const
    {
        return mixin("this " ~ op ~ " param(val)");
    }

    /// ditto
    nothrow pure @safe GenDualNum opBinaryRight(string op)(in real val) const 
    {
        return mixin("param(val) " ~ op ~ " this");
    }

    /// This generates a hash for a generalized dual number.
    nothrow pure @trusted hash_t toHash() const
    {
        auto buf = cast(const(ubyte)*)&_x;

        hash_t res = 0;
        for (auto i = 0; i < real.sizeof; i += hash_t.sizeof)
        {
            for (auto j = 0; j < hash_t.sizeof; j++)
            {
                if (i + j < real.sizeof)
                {
                    res += cast(hash_t) buf[i + j] << 8 * j;
                }
            }
        }
        return res;
    }

    /**
    This generates a string version of a generalized dual number. The form of 
    the string will be *x₀*` + `*x₁*`d + `*x₂*`d2 + `*…*` + `*xₙ*`d`*n* for a 
    generalized dual number of degree *n* with value *x₀*, first derivative *x₁*, 
    second derivative *x₂*, etc.
    */
    pure @safe string toString() const
    {
        return toString(0);
    }

    private pure @safe string toString(in ulong derivOrder) const
    {
        static if (Degree == 1)
            const tail = format("%g%s", _dx, formatDerivOrder(derivOrder + 1));
        else
            const tail = _dx.toString(derivOrder + 1);

        return format!(char, real, string, string)(
            "%g%s + %s", _x, formatDerivOrder(derivOrder), tail);
    }
}

package struct GenDualNum(ulong Degree : 0)
{
    static assert(false, "The degree of generalized dual number must be greater than zero.");
}

// .init
unittest
{
    // force the unit tests
    const GenDualNum!1 w;
    assert(isNaN(w._x));
    assert(isNaN(GenDualNum!1.init._x));
}

// DerivType
unittest
{
    assert(is(GenDualNum!1.DerivType!1 == real));
    assert(is(GenDualNum!2.DerivType!1 == GenDualNum!1));
}

// param
unittest
{
    const q = GenDualNum!1.param(1);
    assert(q._x is 1.0);
    assert(q._dx is 0.0);

    const w = GenDualNum!2.param(1);
    assert(w._x is 1.0);
    assert(w._dx.same(GenDualNum!1.zero));
}

// var
unittest
{
    const q = GenDualNum!1.var(2);
    assert(q._x is 2.0);
    assert(q._dx is 1.0);

    const w = GenDualNum!2.var(2);
    assert(w._x is 2.0);
    assert(w._dx.same(GenDualNum!1.one));
}

// d
unittest
{
    const q = GenDualNum!3(3, GenDualNum!2(2, GenDualNum!1(1, 0)));
    const dq = q.d;
    assert(is(typeof(dq) == const(GenDualNum!2)));
    assert(dq.val == 2);

    const d2q = q.d!2;
    assert(is(typeof(d2q) == const(GenDualNum!1)));
    assert(d2q.val == 1);

    const d3q = q.d!3;
    assert(is(typeof(d3q) == const(real)));
    assert(d3q == 0);

    assert(q.same(q.d!0));
}

// reduce
unittest
{
    const p2 = GenDualNum!2(2, GenDualNum!1(3, 4));
    const p1 = p2.reduce();
    assert(p2._x == p1._x);
    assert(p2._dx._x == p1._dx);

    const p0 = p1.reduce();
    assert(is(typeof(p0) == const(real)));
    assert(p2._x == p0);
}

// inv
unittest
{
    import std.stdio;

    const x = GenDualNum!3(1, GenDualNum!2(2, GenDualNum!1(3, 4)));
    const y = x.inv();
    assert((x * y).same(GenDualNum!3.one));
}

// log
unittest
{
    const p = GenDualNum!2.var(1).log();
    assert(p._x == 0);
    assert(p._dx._x == 1);
    assert(p._dx._dx == -1);

    const z = GenDualNum!1.var(0).log();
    assert(z._x == -real.infinity);
    assert(isNaN(z.d!1));

    const n = GenDualNum!1.var(-1).log();
    assert(isNaN(n._x));
    assert(isNaN(n.d!1));
}

// opBinary(+)
unittest
{
    const q = GenDualNum!1(2, 1);
    const w = GenDualNum!1(4, 3);
    const e = q + w;
    assert(e._x == 6 && e._dx == 4);

    const r = q + 1.0;
    assert(r._x == 3 && r._dx == 1);

    const t = 2 + q;
    assert(t._x == 4 && t._dx == 1);
}

// opBinary(*)
unittest
{
    const q = GenDualNum!1(1, 2);
    const w = GenDualNum!1(3, 5);
    const e = q * w;
    assert(e._x == 3 && e._dx == 11);
}

// opBinary(%)
unittest
{
    const e = GenDualNum!1(3, 1) % GenDualNum!1(2, 4);
    assert(e._x == 1 && e._dx == -3);

    const q = GenDualNum!1(-3, 1) % GenDualNum!1(2, 4);
    assert(q._x == -1 && q._dx == 5);

    const r = GenDualNum!1(3, 1) % GenDualNum!1(-2, 4);
    assert(r._x == 1 && r._dx == 5);

    const t = GenDualNum!1(-3, 1) % GenDualNum!1(-2, 4);
    assert(t._x == -1 && t._dx == -3);

    const w = GenDualNum!1(4, 1) % GenDualNum!1(2, 4);
    assert(w._x == 0 && w._dx == real.infinity);

    const y = GenDualNum!1(0, 1) % GenDualNum!1(2, 4);
    assert(y._x == 0 && y._dx == 1);

    const u = GenDualNum!1(0, 1) % GenDualNum!1(0, 4);
    assert(isNaN(u.val) && isNaN(u.d!1));

    const i = GenDualNum!1(real.infinity, 1) % GenDualNum!1(2, 4);
    assert(isNaN(i.val) && isNaN(i.d!1));

    const o = GenDualNum!1(3, 1) % GenDualNum!1(real.infinity, 1);
    assert(o._x == 3 && o._dx == 1);

    const p = GenDualNum!1(real.nan, real.nan) % GenDualNum!1(2, 4);
    assert(isNaN(p.val) && isNaN(p.d!1));

    const a = GenDualNum!1(3, 1) % GenDualNum!1(real.nan, real.nan);
    assert(isNaN(a.val) && isNaN(a.d!1));
}

// opBinary(^^)
unittest
{
    const q = GenDualNum!1(2, 1) ^^ GenDualNum!1(3, 4);
    assert(q._x == 8 && q._dx == (32 * log(2) + 12));

    const w = GenDualNum!1(0, 1) ^^ GenDualNum!1(3, 4);
    assert(w.val == 0 && isNaN(w.d!1));

    const e = GenDualNum!1(-1, 1) ^^ GenDualNum!1(3, 4);
    assert(e.val == -1 && isNaN(e.d!1));
}

// opCmp
unittest
{
    const q = GenDualNum!1(2, 1);
    assert(q < 3);
    assert(q > 1);
    assert(q <= 2);
    assert(q >= 2);
}

// opEquals
unittest
{
    const q = GenDualNum!1(4, 2);
    assert(q == 4);
    assert(q != 2);
}

// opUnary
unittest
{
    const q = GenDualNum!1(2, 1);
    const w = +q;
    assert(w._x == 2);
    assert(w._dx == 1);

    const e = -q;
    assert(e._x == -2);
    assert(e._dx == -1);
}

// toHash
unittest
{
    auto q = GenDualNum!1(0.1L, 0.2L);
    auto w = GenDualNum!1(0.1L, 0.0L);
    assert(q.toHash() == w.toHash());

    assert(GenDualNum!1.min_normal.toHash() != (2 * GenDualNum!1.min_normal).toHash());
}

package template same(L, R)
if (isImplicitlyConvertible!(L, real) && isImplicitlyConvertible!(R, real))
{
    bool same(in L lhs, in R rhs)
    {
        return lhs == rhs || (isNaN(real(lhs)) && isNaN(real(rhs)));
    }
}

unittest
{
    assert(same(1, 1));
    assert(!same(2, 3));
    assert(same(real.nan, real.nan));
    assert(!same(real.nan, 5.0L));
    assert(same(real.infinity, real.infinity));
    assert(same(-real.infinity, -real.infinity));
    assert(!same(real.infinity, 2.9L));
}

package template same(T, ulong GDNDegree)
if (isImplicitlyConvertible!(T, real))
{
    bool same(in T, in GenDualNum!GDNDegree)
    {
        return false;
    }

    bool same(in GenDualNum!GDNDegree, in T)
    {
        return false;
    }
}

package template same(ulong LDegree, ulong RDegree)
{
    bool same(in GenDualNum!LDegree lhs, in GenDualNum!RDegree rhs)
    {
        static if (LDegree != RDegree)
            return false;
        else
            return same(lhs.val, rhs.val) && same(lhs.d, rhs.d);
    }
}

unittest
{
    const x = GenDualNum!2(1, 2, 3);
    assert(x.same(x));
    assert(!x.same(GenDualNum!2(2, 2, 3)));
    assert(!x.same(GenDualNum!2(1, 1, 3)));
    assert(!x.same(GenDualNum!2(1, 2, 2)));

    const n = GenDualNum!1();
    assert(n.same(n));
}

// XXX - This cannot be implemented until 
// https: //issues.dlang.org/show_bug.cgi?id=22621 is fixed.
// /**
// Constructs a generalized dual number from its value and the values of its 
// derivatives. If there is only one value is the sequence, that value is 
// returned.
//   
// Params:
//     Len = the number of elements in the sequence (optional)
//     derivVals = the an array of the derivative values where index is the 
//         order of the derivative
//  */
// nothrow pure @safe GenDualNum!(Len - 1) derivSeq(ulong Len)(in real[Len] derivVals...)
// {
//     return GenDualNum!(Len - 1)(derivVals);
// }
// 
// /// ditto
// nothrow pure @safe real derivSeq(ulong Len : 1)(in real[Len] val...) 
// {
//     return val[0];
// }
// 
// void derivSeq(ulong Len : 0)(in real[Len] derivVals...) 
// {
//     static assert(false, "there must be a least one element in the sequence");
// }
// 
// unittest
// {
//     assert(derivSeq(1.0) == 1.0, format("%g != 1", derivSeq(1.0)));
//     assert(derivSeq(1.0L, 2.0L).same(GenDualNum!1(1.0, 2.0)));
// }

private pure @safe string formatDerivOrder(in ulong derivOrder)
{
    switch (derivOrder)
    {
    case 0:
        return "";
    case 1:
        return "d";
    default:
        return format("d%d", derivOrder);
    }
}
