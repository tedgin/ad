/**
This module implements automatic differentiation using forward accumulation and
operator overloading. It can only differentiate functions of the form 
$(MATH f:R->R). This is a completely unoptimized version. The traditional 
definition of differentiation is used, not the generalized notion from 
distribution theory.

Macros:
    MATH = <var style="white-space: nowrap">$0</var>
    SUP  = <sup style="font-size: 70%; vertical-align: super;">$0</sup>
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
            more than Degree, the degree of the generalized dual number.
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

    // CONSTANTS

    /// This is a constant _zero represented as a generalized dual number.
    static immutable GenDualNum zero = GenDualNum(0, mkZeroDeriv());

    /// This is a constant _one represented as a generalized dual number.
    static immutable GenDualNum one = GenDualNum(1, mkZeroDeriv());

    // FIELDS

    private real _x;
    private DerivType!1 _dx;

    // CONSTRUCTORS

    /**
    This creates a generalized dual number from its value and the values of its 
    derivatives. 
    
    Params:
        derivVals = the array of the derivative values where the index is the 
            order of the derivative
    */
    this(in real[Degree + 1] derivVals...) nothrow pure @nogc @safe
    {
        _x = derivVals[0];

        static if (Degree == 1)
            _dx = derivVals[1];
        else
            _dx = DerivType!1(derivVals[1 .. Degree + 1]);
    }

    /**
    This creates a generalized dual number by copying an existing one. If the 
    new one's degree is less than the existing one's degree, the higher order 
    derivatives truncated. If the new one's degree is greater than the existing 
    one's degree, the higher order derivative are set to zero.

    Params:
        ThatDegree = the degree of generalized dual number being copied
        that = the generalized dual number being copied
    */
    this(ulong ThatDegree)(in GenDualNum!ThatDegree that) nothrow pure @nogc @safe
    {
        this._x = that._x;

        static if (ThatDegree < Degree && ThatDegree == 1)
            this._dx = DerivType!1.param(that._dx);
        else static if (ThatDegree > Degree && Degree == 1)
            this._dx = that._dx.val;
        else
            this._dx = DerivType!1(that._dx);
    }

    package this(real val, DerivType!1 derivs)
    {
        _x = val;
        _dx = derivs;
    }

    /**
    This function constructs a generalized dual number representing a parameter 
    or constant with respect to differentiation.
        
    Params:
        val = The value of the parameter or constant.
    
    Returns:
        The representation of the parameter or constant
    */
    static GenDualNum param(in real val) nothrow pure @nogc @safe
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
    static GenDualNum var(in real val) nothrow pure @nogc @safe
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
    //     if (same(param.d!1, DerivType!1.zero) {
    //         return var(param.val);
    //     } else {
    //         return param;
    // }

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
    This is the maximum `int` value such that $(MATH 10$(SUP max_10_exp)) is 
    representable as a generalized dual number.
    */    
    static immutable int max_10_exp = real.max_10_exp;

    /** 
    This is the maximum `int` value such that $(MATH 2$(SUP max_exp - 1)) is 
    representable as a generalized dual number.
    */
    static immutable int max_exp = real.max_exp;

    /** 
    This is the minimum `int` value such that $(MATH 10$(SUP min_10_exp)) is 
    representable as a generalized dual number. 
    */
    static immutable int min_10_exp = real.min_10_exp;

    /** 
    This is the minimum `int` value such that $(MATH 2$(SUP min_exp - 1)) is 
    representable as a generalized dual number.
    */
    static immutable int min_exp = real.min_exp;

    /// This is the largest finite generalized dual number. 
    static immutable GenDualNum max = GenDualNum(real.max, mkZeroDeriv());

    /// This is the smallest positive generalized dual number. 
    static immutable GenDualNum min_normal = GenDualNum(real.min_normal, mkZeroDeriv());

    /// This is the real part.
    GenDualNum re() const nothrow pure @safe @nogc
    {
        return this;
    }

    /**
    This is the imaginary part. It is always has a zero value since GenDualNum 
    only supports real numbers.
    */
    GenDualNum im() const nothrow pure @nogc @safe
    {
        return zero;
    }

    // MEMBERS

    /**
    This is the value of the generalized dual number.
        
    Returns:
        the value
    */
    real val() const nothrow pure @nogc @safe
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
    DerivType!Order d(ulong Order = 1)() const nothrow pure @nogc @safe
    if (0 < Order && Order <= Degree)
    {
        static if (Order == 1)
            return _dx;
        else
            return _dx.d!(Order - 1);
    }

    /// ditto
    GenDualNum d(ulong Order : 0)() const nothrow pure @nogc @safe
    {
        return this;
    }

    /**
    This computes the inverse of a generalized dual number. That is, given a 
    generalized dual number $(MATH x), it computes the generalized dual number 
    $(MATH y) where $(MATH xy = 1).
    
    Returns:
        the inverted generalized dual number
    */
    GenDualNum inv() const nothrow pure @nogc @safe
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
    bool opEquals(ulong D)(in GenDualNum!D that) const nothrow pure @nogc @safe
    {
        return this._x == that._x;
    }

    /// ditto
    bool opEquals(in real val) const nothrow pure @nogc @safe
    {
        return _x == val;
    }

    /**
    This provides support for using the operators `<`, `<=`, `>=`, and `>` to 
    compare a generalized dual number to a real or another generalized dual 
    number. If $(MATH x) and $(MATH y) are two generalized dual numbers, 
    $(MATH x < y), if the value of $(MATH x) is less than the value of $(MATH y) 
    regardless of the values of their derivative terms.
    */
    int opCmp(ulong D)(in GenDualNum!D that) const nothrow pure @nogc @safe
    {
        return opCmp(that._x);
    }

    /// ditto
    int opCmp(in real val) const nothrow pure @nogc @safe
    {
        if (_x < val)
            return -1;
        if (_x > val)
            return 1;
        else
            return 0;
    }

    /// This provides support for using the prefix operators `+` and `-`.
    GenDualNum opUnary(string op)() const nothrow pure @nogc @safe
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
    GenDualNum!(ThatDegree < Degree ? ThatDegree : Degree) 
    opBinary(string Op, ulong ThatDegree)(in GenDualNum!ThatDegree that) 
    const nothrow pure @nogc @safe
    if (ThatDegree != Degree)
    {
        static if (ThatDegree < Degree)
        {
            return GenDualNum!ThatDegree(this).opBinary!Op(that);
        }
        else
        {
            return this.opBinary!Op(GenDualNum!Degree(that));
        }
    }

    /// ditto
    GenDualNum opBinary(string Op : "+", ulong ThatDegree : Degree)(in GenDualNum!ThatDegree that) 
    const nothrow pure @nogc @safe
    {
        return GenDualNum(this.val + that.val, this.d + that.d);
    }

    /// ditto
    GenDualNum opBinary(string Op : "-", ulong ThatDegree : Degree)(in GenDualNum!ThatDegree that) 
    const nothrow pure @nogc @safe
    {
        return this + -that;
    }

    /// ditto
    GenDualNum opBinary(string Op : "*", ulong ThatDegree : Degree)(in GenDualNum!ThatDegree that) 
    const nothrow pure @nogc @safe
    {
        return GenDualNum(this._x * that._x, this._dx * that.reduce() + this.reduce() * that._dx);
    }

    /// ditto
    GenDualNum opBinary(string Op : "/", ulong ThatDegree : Degree)(in GenDualNum!ThatDegree that) 
    const nothrow pure @nogc @safe
    {
        return this * that.inv();
    }

    /// ditto
    GenDualNum opBinary(string Op : "%", ulong ThatDegree : Degree)(in GenDualNum!ThatDegree that) 
    const nothrow pure @nogc @safe
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
    GenDualNum opBinary(string Op : "^^", ulong ThatDegree : Degree)(in GenDualNum!ThatDegree that) 
    const nothrow pure @nogc @safe
    {
        const f = this.reduce();
        const fp = this._dx;
        const g = that.reduce();
        const gp = that._dx;
        const fug = f ^^ g;
        return GenDualNum(this._x ^^ that._x, gp * fug * f.log() + fp * fug * g / f);
    }

    /// ditto
    GenDualNum opBinary(string Op)(in real val) const nothrow pure @nogc @safe
    {
        return mixin("this " ~ Op ~ " param(val)");
    }

    /// ditto
    GenDualNum opBinaryRight(string Op)(in real val) const nothrow pure @nogc @safe
    {
        return mixin("param(val) " ~ Op ~ " this");
    }

    /// This generates a hash for a generalized dual number.
    hash_t toHash() const nothrow pure @nogc @trusted
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
    the string will be 
    $(MATH f(x₀))` + `$(MATH f⁽ⁱ⁾(x₀))`dx + `$(MATH f⁽²⁾(x₀))`(dx)² + `$(MATH …)` + `$(MATH f⁽ⁿ⁾(x₀))`(dx)`$(MATH ⁿ) 
    for a generalized dual number of degree $(MATH n) with value $(MATH f(x₀)), 
    first derivative $(MATH f⁽ⁱ⁾(x₀)), second derivative $(MATH f⁽²⁾(x₀)), etc.
    */
    string toString() const pure @safe
    {
        return toString(0);
    }

    private string toString(ulong derivOrder) const pure @safe
    {
        static if (Degree == 1)
            const tail = format("%g%s", _dx, formatDerivOrder(derivOrder + 1));
        else
            const tail = _dx.toString(derivOrder + 1);

        return format!(char, real, string, string)(
            "%g%s + %s", _x, formatDerivOrder(derivOrder), tail);
    }
}

/**
This fails with an explanation when an attempt is made to create a zero degree 
generalized dual number.
*/
struct GenDualNum(ulong Degree : 0)
{
    static assert(false, "The degree of generalized dual number must be greater than zero.");
}

// .init
unittest
{
    const GenDualNum!1 w;
    assert(isNaN(w._x), "_x doesn't default init to NaN");
}

// DerivType
unittest
{
    assert(is(GenDualNum!1.DerivType!1 == real), "The degree 1 derivative type should be real");

    assert(
        is(GenDualNum!2.DerivType!1 == GenDualNum!1), 
        "The degree 2 derivative type should be a degree 1 GenDualNum");

    assert(GenDualNum!1.mkNaNDeriv is real.nan, "The degree 1 derivative NaN should be a real NaN");

    assert(
        GenDualNum!2.mkNaNDeriv is GenDualNum!1(real.nan, real.nan), 
        "The degree 2 derivative NaN isn't correct");

    assert(GenDualNum!2.mkZeroDeriv is GenDualNum!1(0, 0),
            "The degree 2 derivative zero isn't correct");
}

// constructors
unittest
{
    const q = GenDualNum!2(0, 1, 2);
    assert(q._x is 0.0L, "value is set incorrectly");
    assert(same(q._dx, GenDualNum!1(1, 2)), "derivatives are set incorrectly");

    assert(same(GenDualNum!3(0, 1, 2, 0), GenDualNum!3(q)), "GenDualNum lifting not working");
    assert(same(GenDualNum!1(0, 1), GenDualNum!1(q)), "GenDualNum truncating not working");
}

// var
unittest
{
    const q = GenDualNum!1.var(2);
    assert(q._x is 2.0L, "GenDualNum!1.var(2)._x should be 2");
    assert(q._dx is 1.0L, "GenDualNum!1.var(2)._dx should be 1");

    const w = GenDualNum!2.var(2);
    assert(w._x is 2.0, "GenDualNum!2.var(2)._x should be 2");
    assert(same(w._dx, GenDualNum!1.one), "GenDualNum!2.var(2)._dx should be GenDualNum!1.one");
}

// real properties
unittest
{
    const q = GenDualNum!1();
    assert(q.re is q, "The real part of a GenDualNum should be the GenDualNum");
    assert(q.im is GenDualNum!1.zero, "The imaginary part of a GenDualNum should be zero");
}

// d
unittest
{
    const q = GenDualNum!3(3, GenDualNum!2(2, GenDualNum!1(1, 0)));
    assert(same(q, q.d!0));

    const d2q = q.d!2;

    assert(
        is(typeof(d2q) == const(GenDualNum!1)), 
        "The type of GenDualNum!3.d!2 should be GenDaulNum!1");

    assert(d2q.val == 1, "GenDualNum!3.d!2.val has incorrect value");
}

// inv
unittest
{
    const x = GenDualNum!3(1, GenDualNum!2(2, GenDualNum!1(3, 4)));
    assert(same(x * x.inv, GenDualNum!3.one), "x * x.inv should be one");
}

// log
unittest
{
    const p = GenDualNum!2.var(1).log();
    assert(p._x == 0);
    assert(p._dx._x == 1);
    assert(p._dx._dx == -1);
}

// comparison operations
unittest
{
    const q = GenDualNum!1(4, 2);
    assert(q == 4, "GenDualNum!1(4, 2) should equal 4");
    assert(q > 1, "GenDualNum!1(4, 2) should be greater than 1");
    assert(q <= 4, "GenDualNum!(4, 2) should be less than or equal to 4");

    const w = GenDualNum!2(4, 3, 1);
    assert(w == q, "Two GenDualNum objects with the same value should be equal");

    const e = GenDualNum!2(3, 3, 1);
    assert(
        e < w, 
        "If the value of 1 GenDualNum object is less than another, the GenDualNum object should" ~
        " be less than the other.");
}
 
// opUnary
unittest
{
    const q = GenDualNum!1(2, 1);
    assert(same(q, +q), "+q should be the identical to q");
}

// opBinary(+)
unittest 
{
    const gdn1 = GenDualNum!1(4, 5);
    const gdn2 = GenDualNum!2(1, 2, 3);
    const sum = GenDualNum!1(5, 7); 
    assert(same(sum, gdn2 + gdn1), "GenDualNum!2 + GenDualNum!1 not working");
    assert(same(sum, gdn1 + gdn2), "GenDualNum!1 + GenDualNum!2 not working");
}

// opBinary(-)
unittest
{
    assert(
        same(GenDualNum!1(1, -1), GenDualNum!1(2, 1) - GenDualNum!1(1, 2)), 
        "GenDualNum - GenDualNum not working");
}

// opBinary(%)
unittest
{
    assert(
        same(GenDualNum!1(1, -3), GenDualNum!1(3, 1) % GenDualNum!1(2, 4)),
        "GenDualNum % GenDualNum not working");
}

// opBinary(^^)
unittest
{
    const fx = 0.0L;
    const dfx = 1.0L;
    const f = GenDualNum!1(fx, dfx);

    const gx = 2.0L;
    const dgx = 3.0L;
    const g = GenDualNum!1(gx, dgx);

    const hx = fx ^^ gx;
    const dhx = hx * (dgx * log(fx) + gx * dfx / fx); 
    const h = GenDualNum!1(hx, dhx);

    assert(same(h, f ^^ g), "GenDualNum ^^ GenDualNum not working");
}

// toHash
unittest
{
    auto q = GenDualNum!1(0.1L, 0.2L);
    auto w = GenDualNum!1(0.1L, 0.0L);
    assert(q.toHash() == w.toHash(), "toHash isn't working correctly");
}

// toString
unittest
{
    assert(GenDualNum!1(0, 1).toString == "0 + 1dx", "GenDualNum.toString not working");

    assert(
        GenDualNum!2(0, 1, 2).toString(0) == "0 + 1dx + 2(dx)²", 
        "GendDualNum.toString(derivOrder) not working");
}

package template same(L, R)
if (isImplicitlyConvertible!(L, real) && isImplicitlyConvertible!(R, real))
{
    bool same(in L lhs, in R rhs)
    {
        return lhs == rhs || (isNaN(real(lhs)) && isNaN(real(rhs)));
    }
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
    const x = GenDualNum!1(0, 1);
    assert(!same(0, x), "0 should not be the same as GenDualNum!1(0, 1)");
    assert(!same(x, 0), "GenDualNum!1(0, 1) should not be the same as 0");

    const y = GenDualNum!2(0, 1, 0);
    assert(!same(x, y), "GenDualNum objects of different degree cannot be the same");
}

private string formatDerivOrder(ulong derivOrder) pure @safe
{
    switch (derivOrder)
    {
    case 0:
        return "";
    case 1:
        return "dx";
    default:
        return format("(dx)%s", formatPower(derivOrder));
    }
}

private string formatPower(ulong power) pure @safe
{
    switch (power)
    {
    case 0:
        return "\u2070";
    case 1:
        return "\u00B9";
    case 2:
        return "\u00B2";
    case 3:
        return "\u00B3";
    case 4:
        return "\u2074";
    case 5:
        return "\u2075";
    case 6:
        return "\u2076";
    case 7:
        return "\u2077";
    case 8:
        return "\u2078";
    case 9:
        return "\u2079";
    default:
        return formatPower(power / 10) ~ formatPower(power % 10);    
    }
}

unittest
{
    assert(formatPower(0) == "⁰", "formatPower(0) should be ⁰");
    assert(formatPower(1) == "¹", "formatPower(0) should be ¹");
    assert(formatPower(3) == "³", "formatPower(0) should be ³");
    assert(formatPower(4) == "⁴", "formatPower(0) should be ⁴");
    assert(formatPower(5) == "⁵", "formatPower(0) should be ⁵");
    assert(formatPower(6) == "⁶", "formatPower(0) should be ⁶");
    assert(formatPower(7) == "⁷", "formatPower(0) should be ⁷");
    assert(formatPower(8) == "⁸", "formatPower(0) should be ⁸");
    assert(formatPower(9) == "⁹", "formatPower(0) should be ⁹");
    assert(formatPower(10) == "¹⁰", "formatPower(0) should be ¹⁰");
}
