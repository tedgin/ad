/**
 * This module implements univariate automatic differentiation of arbitrary order using forward
 * accumulation. It supports differentiating functions of the form $(MATH f:ℝ→ℝ).
 */
module ad.core;

import std.algorithm: min;
import std.format: format;
import std.math: abs, isFinite, isInfinity, isNaN, LN2, log, signbit;
import std.traits: fullyQualifiedName, isImplicitlyConvertible, TemplateOf;

/**
 * This data structure implements a <em>generalized dual number</em>, a generalization of the dual
 * number that supports derivatives of arbitrary order. The <em>degree</em> of a generalized dual
 * number is the maximum order of derivative that it supports.
 *
 * `GDN` is intended to be a drop-in replacement for a `real`. When no operation is explicitly
 * defined for a `GDN` is will implicitly be converted to a `real`.
 *
 * In addition to differentiation, it supports basic algebraic operations mostly through operator
 * overloading. The `+` and `-` prefix operators are overloaded as well as the `+`, `-`, `*`, `/`,
 * `%`, and `^^` binary, infix operators. When two `GDN`s are combined by a binary operation, the
 * degree of the resulting `GDN` will be the lesser of the degrees of the two input `GDN`s.
 *
 * ```
 * auto x = GDN!2(3);
 * auto y = GDN!1(4);
 * auto z = x + y;
 * assert(typeof(z).DEGREE == 1);
 * ```
 *
 * The infix operators are also overloaded to support combining a `GDN` with a `real`. The `real` is
 * converted to a constant `GDN` of the same degree as the other input `GDN`. A <em>constant</em>
 * generalized dual number is one where the derivative of any order is 0.
 *
 * ```
 * auto x = GDN!1(2);
 * real y = 3;
 *
 * auto u = x + y;
 * assert(u == 5.0L && u.d == 1);
 *
 * auto v = y + x;
 * assert(v is u);
 * ```
 *
 * Params:
 *   Degree = the highest order of derivative that can be taken
*/
struct GDN(ulong Degree = 1) if (Degree > 0)
{
    /// The degree of the Generalized dual number
    enum ulong DEGREE = Degree;

    /**
     * This is the type constructor of the derivatives.
     *
     * Params:
     *   Order = the order of the derivative
     */
    template DerivType(ulong Order = 1) if (Order <= Degree)
    {
        static if (Order < Degree)
            alias DerivType = GDN!(Degree - Order);
        else
            alias DerivType = real;
    }

    // Constructs an object that has the same type as the derivative where all elements are NaN.
    package static pure nothrow @nogc @safe DerivType!1 mkNaNDeriv()
    {
        static if (Degree == 1)
            return real.nan;
        else
            return DerivType!1(real.nan, DerivType!1.mkNaNDeriv());
    }

    // Constructs an object that has the same type as the derivative where all elements are 0.
    package static pure nothrow @nogc @safe DerivType!1 mkZeroDeriv()
    {
        static if (Degree == 1)
            return 0;
        else
            return DerivType!1(0, DerivType!1.mkZeroDeriv());
    }

    // FIELDS

    // The value of the generalized dual number.
    private real _x;

    // The derivative of the generalized dual number.
    private DerivType!1 _dx;

    // CONSTRUCTORS

    /**
     * This constructs a generalized dual number representing the variable of differentiation. The
     * first derivative, $(MATH dx/dx = 1), and all higher order derivatives are $(MATH 0).
     *
     * Params:
     *   val = The value of the variable.
     */
    pure nothrow @nogc @safe this(in real val)
    {
        _x = val;

        static if (Degree == 1)
            _dx = 1;
        else
            _dx = GDN!(Degree - 1)(1, GDN!(Degree - 1).mkZeroDeriv());
    }

    /**
     * This creates a generalized dual number from its value and the values of its derivatives.
     *
     * Params:
     *   derivVals = the array of the derivative values where the index is the order of the
     *       derivative, i.e.,
     *       `[$(MATH f(x₀)), $(MATH f$(SUP (1))(x₀)), $(MATH f$(SUP (2))(x₀)), $(MATH …)]`.
     */
    pure nothrow @nogc @safe this(in real[Degree + 1] derivVals...)
    {
        _x = derivVals[0];

        static if (Degree == 1)
            _dx = derivVals[1];
        else
            _dx = DerivType!1(derivVals[1 .. Degree + 1]);
    }

    /**
     * This creates a generalized dual number by copying an existing one. If the new one's degree is
     * less than the existing one's degree, the higher order derivatives are truncated. If the new
     * one's degree is greater than the existing one's degree, the higher order derivative are set
     * to zero.
     *
     * Params:
     *   ThatDegree = the degree of generalized dual number being copied
     *   that = the generalized dual number being copied
     */
    pure nothrow @nogc @safe this(ulong ThatDegree)(in GDN!ThatDegree that)
    {
        this._x = that._x;

        static if (ThatDegree < Degree && ThatDegree == 1)
            this._dx = DerivType!1(that._dx, DerivType!1.mkZeroDeriv());
        else static if (ThatDegree > Degree && Degree == 1)
            this._dx = that._dx.val;
        else
            this._dx = DerivType!1(that._dx);
    }

    // This constructs a generalized dual number from its value and its first derivative.
    package pure nothrow @nogc @safe this(in real val, in DerivType!1 derivs)
    {
        _x = val;
        _dx = derivs;
    }

    /*
    This constructs a generalized dual number that has all of its derivative of each order set to
    zero.
    */
    package static pure nothrow @nogc @safe GDN mkConst(in real val)
    {
        return GDN(val, mkZeroDeriv());
    }

    // PROPERTIES

    /// This is a constant zero represented as a generalized dual number.
    enum auto zero = GDN.mkConst(0);

    /// This is a constant one represented as a generalized dual number.
    enum auto one = GDN.mkConst(1);

    /// This is a constant generalized dual number representing infinity.
    enum auto infinity = GDN.mkConst(real.infinity);

    /// This is a generalized dual number representing `NaN`. All derivatives are NaN as well.
    enum auto nan = GDN();

    /// This is the smallest generalized dual number such that `one + epsilon > one`.
    enum auto epsilon = GDN.mkConst(real.epsilon);

    /// This is the largest finite generalized dual number.
    enum auto max = GDN.mkConst(real.max);

    /// This is the smallest positive generalized dual number.
    enum auto min_normal = GDN.mkConst(real.min_normal);

    /// This is the number of decimal digits of precision of the value.
    enum int dig = real.dig;

    /// This is the number of bits in the mantissa of the value.
    enum int mant_dig = real.mant_dig;

    /**
     * This is the maximum `int` value such that $(MATH 10$(SUP max_10_exp)) is representable as a
     * generalized dual number.
     */
    enum int max_10_exp = real.max_10_exp;

    /**
     * This is the maximum `int` value such that $(MATH 2$(SUP max_exp - 1)) is representable as a
     * generalized dual number.
     */
    enum int max_exp = real.max_exp;

    /**
     * This is the minimum `int` value such that $(MATH 10$(SUP min_10_exp)) is representable as a
     * generalized dual number.
     */
    enum int min_10_exp = real.min_10_exp;

    /**
     * This is the minimum `int` value such that $(MATH 2$(SUP min_exp - 1)) is representable as a
     * generalized dual number.
     */
    enum int min_exp = real.min_exp;

    /**
     * This is the real part. It is identical to the generalized dual number, since `GDN` only
     * supports real numbers.
     */
    pure nothrow @safe @nogc GDN re() const
    {
        return this;
    }

    /**
     * This is the imaginary part. It is always has a zero value, since `GDN` only supports real
     * numbers.
     */
    pure nothrow @nogc @safe GDN im() const
    {
        return GDN.zero;
    }

    // MEMBERS

    /// This is the value of the generalized dual number.
    pure nothrow @nogc @safe real val() const
    {
        return _x;
    }

    /**
     * This is the derivative of order `Ord` of the generalized dual number. The order must be at
     * least one but no more that the degree of the generalized dual number. If `Ord == Degree`, the
     * derivative will be of type `real`. Otherwise it will be a `GDN` but with degree
     * `Degree - Ord`.
     *
     * Params:
     *   Ord = the order of the derivate to compute, default `1`
     *
     * Returns:
     *   The derivative of order `Ord` of the generalized dual number.
     *
     * Examples:
     *   ```
     *   auto x = GDN!3(2, 3, -1, -2);
     *   assert(x.d.val == 3);
     *   assert(x.d!3 == -2);
     *   assert(x.d!0 is x);
     *   ```
     */
    pure nothrow @nogc @safe DerivType!Ord d(ulong Ord = 1)() const if (0 < Ord && Ord <= Degree)
    {
        static if (Ord == 1)
            return _dx;
        else
            return _dx.d!(Ord - 1);
    }

    /// ditto
    pure nothrow @nogc @safe GDN d(ulong Ord : 0)() const
    {
        return this;
    }

    // This function evaluated the Dirac delta function of the generalized dual number.
    pure nothrow @nogc @safe package GDN dirac() const
    {
        static if (Degree == 1)
            const df = _x != 0
                ? 0.0L : (signbit(_x) == 1 ? real.infinity : -real.infinity);
        else
            const df = reduce().dirac() * (signbit(_x) == 1 ? 1 : -1);

        return GDN(_x == 0 ? real.infinity : 0, df * _dx);
    }

    /**
     * This computes the inverse of a generalized dual number. That is, given a generalized dual
     * number $(MATH g), it computes the generalized dual number $(MATH h) (or $(MATH g$(SUP -1)))
     * where $(MATH gh = 1).
     *
     * If $(MATH f(x) = g$(SUP -1)(x)), then $(MATH f' = -g$(SUP -2)g').
     *
     * Returns:
     *   the inverse of the generalized dual number
     *
     * Examples:
     *   ```
     *   auto x = GDN!1(2, 3);
     *   auto y = x.inv;
     *   assert(y == 0.5L && y.d == -0.75);
     *   ```
     */
    pure nothrow @nogc @safe GDN inv() const
    {
        const reduced = reduce();
        return GDN(1 / _x, -_dx / (reduced * reduced));
    }

    /*
    Computes the natural logarithm of a generalized dual number.

    It is defined in this module instead of core, because it is required to compute the derivative
    of the ^^ operator.
    */
    package pure nothrow @nogc @safe GDN log() const
    {
        static import std.math;

        if (signbit(_x) == 1)
            return nan;
        return GDN(std.math.log(_x), _dx / reduce());
    }

    /*
    This copies this generalized dual number with the copy have one degree less than the original.
    As a result the highest order derivative is removed.
    */
    package pure nothrow @nogc @safe auto reduce() const
    {
        static if (Degree > 1)
            return GDN!(Degree - 1)(_x, _dx.reduce);
         else
            return _x;
    }

    // OPERATOR OVERLOADING

    /**
     * This provides support for the `cast()` operator. It allows casting from a `GDN` to one with a
     * different degree.
     *
     * Params:
     *   T = the type of GDN to cast to
     *
     * Returns:
     *   a generalized dual number with the same value as the original but with a different degree
     *
     * Examples:
     *   ```
     *   auto x = GDN!1(2);
     *   auto y = cast(GDN!2) x;
     *   assert(y == x && typeof(y).DEGREE == 2);
     *   ```
     */
    pragma(inline, true) pure nothrow @nogc @safe T opCast(T)() const
    if (fullyQualifiedName!(TemplateOf!T) == "ad.core.GDN") {
        static if (T.DEGREE == Degree)
            return this;
        else
            return T(this);
    }

    /**
     * This provides support for using the operators `==` and `!=` to compare a generalized dual
     * number to a real or another generalized dual number. Two generalized dual numbers are equal
     * if their values are equal regardless of the values of their derivative terms.
     *
     * Params:
     *   ThatDeg = the degree of the generalized dual number being compared
     *   that = the generalized dual number or `real` being compared
     *
     * Returns:
     *   `true` if the values are equal, `false` otherwise
     *
     * Examples:
     *   ```
     *   auto x = GDN!1(2, 3);
     *   auto y = GDN!2(2, 1, 0);
     *   assert(x == y);
     *   assert(x == 2.0L);
     *   ```
     */
    pure nothrow @nogc @safe bool opEquals(ulong ThatDeg)(in GDN!ThatDeg that) const
    {
        return this._x == that._x;
    }

    /// ditto
    pure nothrow @nogc @safe bool opEquals(in real val) const
    {
        return _x == val;
    }

    /**
     * This provides support for using the operators `<`, `<=`, `>=`, and `>` to compare a
     * generalized dual number to a `real` or another generalized dual number. If $(MATH x) and
     * $(MATH y) are two generalized dual numbers, $(MATH x < y), if the value of $(MATH x) is less
     * than the value of $(MATH y) regardless of the values of their derivative terms.
     *
     * Params:
     *   ThatDeg = the degree of the generalized dual number being compared
     *   that = the generalized dual number or `real` being compared
     *
     * Examples:
     *   ```
     *   auto x = GDN!1(2, 3);
     *   auto y = GDN!2(3, 1, 0);
     *   auto z = GDN!1(2, 2);
     *   assert(x < y);
     *   assert(x <= z);
     *   assert(x > 1.0L);
     *   ```
     */
    pure nothrow @nogc @safe int opCmp(ulong ThatDeg)(in GDN!ThatDeg that) const
    {
        return opCmp(that._x);
    }

    /// ditto
    pure nothrow @nogc @safe int opCmp(in real val) const
    {
        if (_x < val)
            return -1;
        if (_x > val)
            return 1;
        else
            return 0;
    }

    /** <b>+g</b>
     *
     * This defines the identity operator for a generalized dual number.
     *
     * If $(MATH f(x) = +g(x)), then $(MATH f' = g').
     *
     * Returns:
     *   the generalized dual number
     */
    pure nothrow @nogc @safe GDN opUnary(string Op : "+")() const
    {
        return this;
    }

    /** <b>-g</b>
     *
     * This negates a generalized dual number.
     *
     * If $(MATH f(x) = -g(x)), then $(MATH f' = -g').
     *
     * Returns:
     *  the negated generalized dual number
     */
    pure nothrow @nogc @safe GDN opUnary(string Op : "-")() const
    {
        return GDN(-_x, -_dx);
    }

    /**
     * This ensures that when two generalized dual numbers are combined, the degree of the resulting
     * generalized dual number is the lesser of the degrees of the two being combined.
     *
     * Params:
     *   Op = the operator being applied
     *   ThatDegree = the degree of the generalized dual number being combined
     *   that = the generalized dual number being combined
     *
     * Returns:
     *   the GDN resulting from the of the operation
     */
    pure nothrow @nogc @safe
    GDN!(min(ThatDegree, Degree)) opBinary(string Op, ulong ThatDegree)(in GDN!ThatDegree that)
    const
    if (ThatDegree != Degree) {
        alias Res = typeof(return);
        return (cast(Res) this).opBinary!Op(cast(Res) that);
    }

    /** <b>g + h</b>
     *
     * This adds one generalized dual number to another. If either term has type `real`, that term
     * is converted to a constant generalized dual number with the same degree as the other.
     *
     * If $(MATH f(x) = g(x) + h(x)), then $(MATH f' = g' + h').
     *
     * Params:
     *   that = the addend
     *
     * Returns:
     *   the sum of the two generalized dual numbers with degree being the lesser of the degree of
     *   the augend and the degree of the addend
     *
     * Examples:
     *   ```
     *   auto x = GDN!2(2);
     *   auto y = GDN!3(3);
     *   auto z = x + y;
     *   auto w = 5 + x;
     *   assert(z == 5 && z.d == 2 && z.d!2 == 0);
     *   assert(w == 7 && w.d == 1 && w.d!2 == 0);
     *   ```
     */
    pure nothrow @nogc @safe
    GDN opBinary(string Op : "+", ulong ThatDegree : Degree)(in GDN!ThatDegree that) const
    {
        return GDN(this.val + that.val, this.d + that.d);
    }

    /** <b>g - h</b>
     *
     * This subtracts one generalized dual number from another. If either term has type `real`, it
     * is converted to a constant generalized dual number with the same degree as the other.
     *
     * If $(MATH f(x) = g(x) - h(x)), then $(MATH f' = g' - h').
     *
     * Params:
     *   that = the subtrahend
     *
     * Returns:
     *   the difference between the two generalized dual numbers with degree being the lesser of the
     *   degree of the minuend and the degree of the subtrahend
     *
     * Examples:
     *   ```
     *   auto x = GDN!1(2);
     *   auto y = GDN!1(3);
     *   auto z = x - y;
     *   auto w = 5 - x;
     *   assert(z == -1 && z.d == 0);
     *   assert(w == 3 && w.d == -1);
     *   ```
     */
    pure nothrow @nogc @safe
    GDN opBinary(string Op : "-", ulong ThatDegree : Degree)(in GDN!ThatDegree that) const
    {
        return this + -that;
    }

    /** <b>g * h</b>
     *
     * This multiplies one generalized dual number by another. If either factor has type `real`, it
     * is converted to a constant generalized dual number with the same degree as the other.
     *
     * If $(MATH f(x) = g(x)h(x)), then $(MATH f' = g'h + gh').
     *
     * Params:
     *   that = the multiplicand
     *
     * Returns:
     *   the product of the two generalized dual numbers with degree being the lesser of the degree
     *   of the multiplier and the degree of the multiplicand.
     *
     * Examples:
     *   ```
     *   auto x = GDN!1(2);
     *   auto y = GDN!1(3);
     *   auto z = x * y;
     *   auto w = 5 * x;
     *   assert(z == 6 && z.d == 5);
     *   assert(w == 10 && w.d == 5);
     *   ```
     */
    pure nothrow @nogc @safe
    GDN opBinary(string Op : "*", ulong ThatDegree : Degree)(in GDN!ThatDegree that) const
    {
        const prod = this._x * that._x;

        if (isNaN(prod))
            return nan;
        return GDN(prod, this._dx * that.reduce() + this.reduce() * that._dx);
    }

    /** <b>g / h</b>
     *
     * This divides one generalized dual number by another. If either the dividend or the divisor
     * has type `real`, it is converted to a constant generalized dual number with the same degree
     * as the divisor or dividend, respectively.
     *
     * If $(MATH f(x) = g(x)/h(x)), then $(MATH f' = (g'h - gh')/h$(SUP 2)).
     *
     * Params:
     *   that = the divisor
     *
     * Returns:
     *   the quotient of the two generalized dual numbers with degree being the lesser of the degree
     *   of the dividend and the degree of the divisor.
     *
     * Examples:
     *   ```
     *   auto x = GDN!1(2);
     *   auto y = GDN!1(-1);
     *   auto z = x / y;
     *   auto w = 5 / x;
     *   assert(z == -2 && z.d == -3);
     *   assert(w == 2.5 && w.d == -1.25);
     *  ```
     */
    pure nothrow @nogc @safe
    GDN opBinary(string Op : "/", ulong ThatDegree : Degree)(in GDN!ThatDegree that) const
    {
        return this * that.inv();
    }

    /** <b>g % h</b>
     *
    * This computes the modulus (remainder) of one generalized dual number divided by another. If
    * either the dividend or the divisor has type `real`, it is converted to a constant generalized
    * dual number with the same degree as the divisor or dividend, respectively.
    *
    * If $(MATH f(x) = g(x) (mod h(x))), then $(MATH f' = g' - (g - f)h'/h - δ(f)(g'h - gh')/h),
    * where $(MATH δ) is the Dirac Delta function.
    *
    * Params:
    *   that = the divisor
    *
    * Returns:
    *   the remainder of the two generalized dual numbers with degree being the lesser of the degree
    *   of the dividend and the degree of the divisor.
    *
    * Examples:
    *   ```
    *   auto x = GDN!1(2);
    *   auto y = GDN!1(-1);
    *   auto z = x % y;
    *   auto w = 5 % x;
    *   assert(z == 0 && z.d == -real.infinity);
    *   assert(w == 1 && w.d == -2);
    *   ```
    */
    /*
    Derivation:
    f (mod g) = f - ⌊f/g⌋g ⇒ ⌊f/g⌋ = {f - [f (mod g)]}/g
    d[f (mod g)] = f' - d(⌊f/g⌋)g - ⌊f/g⌋g' = f' - d⌊f/g⌋d(f/g)g - ⌊f/g⌋g'
           d⌊x⌋ = ∑ᵢ 𝛿(x-i), i∈ℤ, i.e., 0 when x ∉ ℤ, and a Dirac Delta function, 𝛿(x-x), otherwise.
           f/g ∈ ℤ when f (mod g) = 0, d⌊f/g⌋ = 𝛿(f (mod g))
       = f' - {f - [f (mod g)]}g'/g - 𝛿(f (mod g))gd(f/g)
    */
    pure nothrow @nogc @safe
    GDN opBinary(string Op : "%", ulong ThatDegree : Degree)(in GDN!ThatDegree that) const
    {
        const x = this.reduce(), dx = this.d;
        const y = that.reduce(), dy = that.d;
        const v = this / that, dv = v.d;
        const w = x % y;

        static if (Degree == 1)
            const diracW = w == 0 ? real.infinity : 0.0L;
        else
            const diracW = w.dirac();

        const term2 = isInfinity(that.val) ? zero.reduce() : diracW * y * dv;
        const term3 = (x - w) * dy / y;
        const dw = dx - term2 - term3;

        static if (Degree == 1)
            return GDN(w, dw);
        else
            return GDN(w.val, dw);
    }

    /** <b>g ^^ h</b>
     *
     * This raises one generalized dual number to the pow of another. If the base or exponent has
     * type `real`, it is converted to a constant generalized dual number with the same degree as
     * the exponent or base, respectively.
     *
     * If $(MATH f(x) = g(x)$(SUP h(x))), then $(MATH f' = g$(SUP h)[g'h/g + h'ln(g)]).
     *
     * Params:
     *   that = the exponent
     *
     * Returns:
     *   the power of the two generalized dual numbers with degree being the lesser of the degree of
     *   the base and the degree of the exponent.
     *
     * Examples:
     *   ```
     *   static import std.math;
     *
     *   auto x = GDN!1(2, -1);
     *   auto y = GDN!1(-2, 3);
     *   auto u = x ^^ y;
     *   auto w = 3 ^^ x;
     *   assert(u == 0.25 && u.d == 0.25 + 0.75*std.math.LN2);
     *   assert(w == 9 && std.math.isClose(w.d, -std.math.log(19_683)));
     *   ```
     */
    pure nothrow @nogc @safe
    GDN opBinary(string Op : "^^", ulong ThatDegree : Degree)(in GDN!ThatDegree that) const
    {
        const g = this.reduce();
        const gp = this._dx;
        const h = that.reduce();
        const hp = that._dx;

        static if (DEGREE == 1)
            const hpNaN = isNaN(hp);
        else
            const hpNaN = isNaN(hp.val);

        if (this._x == 0 && signbit(this._x) == 0 && that._x > 0 && !hpNaN)
        {
            return GDN(this._x ^^ that._x, g ^^ (h - 1) * gp * h);
        }

        return GDN(this._x ^^ that._x, g ^^ h * (gp * h / g + hp * g.log()));
    }

    /// ditto
    pure nothrow @nogc @safe GDN opBinary(string Op : "^^")(in real c) const
    {
        if (c == 0)
        {
            return GDN.one;
        }

        return GDN(_x ^^ c, c * reduce() ^^ (c - 1) * _dx);
    }

    /**
     * This allows real numbers to be used on the right-hand side of the  +, -, *, /, %, and ^^
     * operators. The real number is promoted to generalized dual number of the same degree as the
     * left-hand side with all derivatives being zero, a constant.
     *
     * Params:
     *   Op = the operator being applied
     *   val = the real number being combined on the right-hand side
     *
     * Returns:
     *   the GDN resulting from the of the operation
     */
    pure nothrow @nogc @safe GDN opBinary(string Op)(in real val) const
    {
        return mixin("this " ~ Op ~ " GDN(val, mkZeroDeriv())");
    }

    /**
     * This allows real numbers to be used on the left-hand side of the  +, -, *, /, %, and ^^
     * operators. The real number is promoted to generalized dual number of the same degree as the
     * right-hand side with all derivatives being zero, a constant.
     *
     * Params:
     *   Op = the operator being applied
     *   val = the real number being combined on the left-hand side
     *
     * Returns:
     *   the GDN resulting from the of the operation
     */
    pure nothrow @nogc @safe GDN opBinaryRight(string Op)(in real val) const
    {
        return mixin("GDN(val, mkZeroDeriv()) " ~ Op ~ " this");
    }

    /**
     * This generates a hash for a generalized dual number.
     *
     * Returns:
     *   the hash of the generalized dual number
     */
    pure nothrow @nogc @trusted hash_t toHash() const
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
     * This generates a string version of a generalized dual number. The form of the string will be
     * $(MATH f(x₀))` + `$(MATH f⁽¹⁾(x₀))`dx + `$(MATH f⁽²⁾(x₀))`(dx)² + `$(MATH …)` + `$(MATH f⁽ⁿ⁾(x₀))`(dx)`$(MATH ⁿ)
     * for a generalized dual number of degree $(MATH n) with value $(MATH f(x₀)), first derivative
     * $(MATH f⁽¹⁾(x₀)), second derivative $(MATH f⁽²⁾(x₀)), etc.
     *
     * Returns:
     *   the string representation of the generalized dual number
     *
     * Examples:
     *   ```
     *   const x = GDN!6(1, -1, +0., -0., real.infinity, -real.infinity, -real.nan);
     *   assert(x.toString == "1 + -1dx + +0(dx)² + -0(dx)³ + ∞(dx)⁴ + -∞(dx)⁵ + NaN(dx)⁶");
     *   ```
     */
    pure @safe string toString() const
    {
        return toString(0);
    }

    // Generates the string representation of the generalized dual number with a given derivative
    // order and higher.
    private pure @safe string toString(ulong derivOrd) const
    {
        static if (Degree == 1)
            const tail = format("%s%s", fmtNum(_dx), fmtDerivSuffix(derivOrd + 1));
        else
            const tail = _dx.toString(derivOrd + 1);

        return format("%s%s + %s", fmtNum(_x), fmtDerivSuffix(derivOrd), tail);
    }
}

// .init
unittest
{
    const GDN!2 w;
    assert(isNaN(w._x), "_x doesn't default init to NaN");
    assert(isNaN(w._dx._x), "_dx doesn't default init to NaN");
    assert(isNaN(w._dx._dx), "_dx._dx doesn't default init to NaN");
    assert(w.DEGREE == 2, "w has incorrect degree");
}

// DerivType
unittest
{
    assert(is(GDN!1.DerivType!1 == real), "The degree 1 derivative type should be real");
    assert(is(GDN!2.DerivType!1 == GDN!1), "The degree 2 derivative type should be a degree 1 GDN");
    assert(isNaN(GDN!1.mkNaNDeriv), "The degree 1 derivative NaN should be a real NaN");

    const x = GDN!2.mkNaNDeriv;
    assert(isNaN(x.val) && isNaN(x.d), "The degree 2 derivative NaN isn't correct");
    assert(GDN!2.mkZeroDeriv == GDN!1(0, 0), "The degree 2 derivative zero isn't correct");
}

// constructors
unittest
{
    const n = GDN!1();
    assert(isNaN(n._x), "The default GDN has value NaN");
    assert(isNaN(n._dx), "The default GDN has derivative NaN");

    const q = GDN!1(2);
    assert(q._x == 2.0L, "GDN!1(2)._x should be 2");
    assert(q._dx == 1.0L, "GDN!1(2)._dx should be 1");

    const w = GDN!2(2);
    assert(w._x == 2.0, "GDN!2(2)._x should be 2");
    assert(w._dx._x == 1 && w._dx._dx == 0, "GDN!2(2) should be GDN!1(1, 0)");

    const e = GDN!2(0, 1, 2);
    assert(e._x == 0.0L, "value is set incorrectly");
    assert(e._dx._x == 1 && e._dx._dx == 2, "derivatives are set incorrectly");

    const r = GDN!3(0, 1, 2, 0);
    const t = GDN!3(e);
    assert(r == t && r.d == t.d && r.d!2 == t.d!2 && r.d!3 == t.d!3, "GDN lifting not working");

    const y = GDN!1(0, 1);
    const u = GDN!1(e);
    assert(y == u && y.d == u.d, "GDN truncating not working");
}

// real properties
unittest
{
    assert(isInfinity(GDN!1.infinity._x), "infinity should be a variable of value infinity");
    assert(GDN!1.infinity._dx == 0, "infinity should have a derivative of 0");

    assert(isNaN(GDN!1.nan._x), "The value of NaN should be NaN");
    assert(isNaN(GDN!1.nan._dx), "The derivative of NaN should be NaN");

    const q = GDN!1();
    assert(q.re is q, "The real part of a GDN should be the GDN");
    assert(q.im is GDN!1.zero, "The imaginary part of a GDN should be zero");
}

// d
unittest
{
    const q = GDN!3(3, GDN!2(2, GDN!1(1, 0)));
    assert(q is q.d!0);

    const d2q = q.d!2;
    assert(is(typeof(d2q) == const(GDN!1)), "The type of GDN!3.d!2 should be GDN!1");
    assert(d2q.val == 1, "GDN!3.d!2.val has incorrect value");
}

// dirac
unittest
{
    const x = GDN!1(3, 4).dirac();
    assert(x.val == 0 && x.d == 0);

    const y = GDN!1(-0.).dirac();
    assert(y.val == real.infinity && y.d == real.infinity);

    const z = GDN!1(+0.).dirac();
    assert(z.val == real.infinity && z.d == -real.infinity);

    const w = GDN!2.one.dirac();
    assert(w.val == 0 && w.d == 0 && w.d!2 == 0, format("w != %s", w));

    const v = GDN!2(-0.).dirac();
    assert(v.val == real.infinity && v.d == real.infinity && isNaN(v.d!2), format("v != %s", v));

    const u = GDN!2(+0.).dirac();
    assert(u.val == real.infinity && u.d == -real.infinity && isNaN(u.d!2), format("u != %s", u));
}

// inv
unittest
{
    const x = GDN!1(3, 4);
    const res = x.inv;
    assert(res.val == 1. / 3 && res.d == -4. / 9, "x.inv is incorrect");

    const pz = GDN!1(+0.).inv;
    assert(isInfinity(pz._x) && pz._x > 0, "inverse of +0 should be inf");
    assert(isInfinity(pz._dx) && pz._dx < 0, "derivative of inverse of +0 should be -inf");

    const nz = GDN!1(-0.).inv;
    assert(isInfinity(nz._x) && nz._x < 0, "inverse of -0 should be -inf");
    assert(isInfinity(nz._dx) && nz._dx < 0, "inverse of -0 derivative should be -inf");
}

// log
unittest
{
    const p = GDN!2(1).log();
    assert(p.val == 0, "log(1) is incorrect");
    assert(p.d.val == 1, "derivative of log(1) is incorrect");
    assert(p.d!2 == -1, "second derivative of log(1) should be -1");

    const e = GDN!1(-0.);
    const q = e.log();
    assert(isNaN(q.val) && isNaN(q.d), "log(-0) should be NaN");

    assert(GDN!1(+0.).log() is GDN!1(-real.infinity, real.infinity), "log(+0) is incorrect");
    assert(GDN!1.infinity.log() is GDN!1(real.infinity, 0), "log(inf) incorrect");
    assert(GDN!1.nan.log() is GDN!1.nan, "log(nan) should be nan");
}

// comparison operations
unittest
{
    const q = GDN!1(4, 2);
    assert(q == 4, "GDN!1(4, 2) should equal 4");
    assert(q > 1, "GDN!1(4, 2) should be greater than 1");
    assert(q <= 4, "GDN!(4, 2) should be less than or equal to 4");

    const w = GDN!2(4, 3, 1);
    assert(w == q, "Two GDN objects with the same value should be equal");

    const e = GDN!2(3, 3, 1);
    assert(
        e < w,
        "If the value of 1 GDN object is less than another, the GDN object should be less than the"
            ~ " other.");

    const nz = GDN!1(-0.);
    const pz = GDN!1(+0.);
    assert(nz == pz, "-0 != +0");
}

// opUnary
unittest
{
    const q = GDN!1(2, 1);
    assert(q is +q, "+q should be the identical to q");
    assert(GDN!1(-2, -1) is -q, "-q should be the negation of q and all its derivatives");

    const nz = GDN!1(-0.);
    const pz = GDN!1(+0.);
    assert(!(nz is -pz), "-GDN(0) should not be the same as GDN(-0)");

    const ni = -GDN!1(real.infinity);
    assert(isInfinity(ni._x) && ni._x < 0, "-inf has incorrect value");
    assert(ni._dx == -1, "-inf has incorrect derivative");
}

// opBinary(+)
unittest
{
    const gdn1 = GDN!1(4, 5);
    const gdn2 = GDN!2(1, 2, 3);
    const sum = GDN!1(5, 7);
    assert(sum is gdn2 + gdn1, "GDN!2 + GDN!1 not working");
    assert(sum is gdn1 + gdn2, "GDN!1 + GDN!2 not working");

    const zi = GDN!1(0) + GDN!1(real.infinity);
    assert(isInfinity(zi._x) && zi._dx == 2, "0 + inf is incorrect");

    const zn = GDN!1.zero + GDN!1.nan;
    assert(isNaN(zn._x) && isNaN(zn._dx), "0 + nan is incorrect");

    const ni = GDN!1.nan + GDN!1.infinity;
    assert(isNaN(ni._x) && isNaN(ni._dx), "nan + inf is incorrect");
}

// opBinary(-)
unittest
{
    assert(GDN!1(1, -1) is GDN!1(2, 1) - GDN!1(1, 2), "GDN - GDN not working");
}

// opBinary(*)
unittest
{
    const q = GDN!1(2, -3);
    const w = GDN!1(5, 7);
    const nz = GDN!1(-0.);
    const e = q * w;

    assert(e._x == 10, "e.x should be 10");
    assert(e._dx == -1, "e.dx should be -1");

    const r = w * q;
    assert(r._x == e._x && r._dx == e._dx, "multiplication should be commutative");

    const nznz = nz * nz;
    assert(signbit(nznz._x) == 0 && signbit(nznz._dx) == 1, "-0 * -0 incorrect");

    const pz = GDN!1(+0.);
    const nzpz = nz * pz;
    assert(signbit(nzpz._x) == 1, "-0 * +0 should have -0 value");
    // NOTE: -0 + +0 is invalid, so nzpz._dx could be either -0 or +0.

    const zi = GDN!1.zero * GDN!1.infinity;
    assert(isNaN(zi._x) && isNaN(zi._dx), "0 * inf incorrect");

    const zn = GDN!1.zero * GDN!1.nan;
    assert(isNaN(zn._x) && isNaN(zn._dx), "0 * NaN incorrect");

    const ni = GDN!1.nan * GDN!1.infinity;
    assert(isNaN(ni.val) && isNaN(ni.d), "NaN * inf incorrect");

    const g1 = GDN!1(2, 3);
    const g2 = GDN!2(-1, 5, 7);
    const g1g2 = g1 * g2;
    const g2g1 = g2 * g1;
    assert(g1g2.DEGREE == 1, "g1 * g2 should have degree 1");
    assert(g2g1.DEGREE == 1, "g2 * g1 should have degree 1");
    assert(g1g2.val == -2 && g1g2.d == 7, "g1 * g2 is incorrect");

    const t = GDN!1(-1, 0) * GDN!1(1, -real.infinity);
    assert(t.val == -1 && t.d == real.infinity);
}

// opBinary(/)
unittest
{
    const q = GDN!1(6);
    const w = GDN!1(2);
    const r = q / w;

    assert(r._x == 3 && r._dx == -1, "division isn't correct");
}

// opBinary(%)
unittest
{
    const x = GDN!1(2, 3);
    const s = GDN!1(5, 4);
    const e = GDN!1(-0., 1);

    const d = s % x;
    assert(d.val == 1 && d.d == -2, "% positive divisor test failed");

    const f = GDN!1(-3, -1);
    const g = s % f;
    assert(g.val == 2 && g.d == 3, "% negative divisor test failed");

    const q = x % x;
    assert(q.val == 0 && isNaN(q.d), "x % x is incorrect");

    const w = e % x;
    assert(w.val == 0 && signbit(w.val) == 1, format("-0 %% 2 = %s is incorrect value", w.val));
    assert(isInfinity(w.d) && signbit(w.d) == 1, "-0 % 1 has incorrect derivative");

    const r = GDN!1(+0., 1);
    const t = r % x;
    assert(t.val == 0 && signbit(t.val) == 0, format("+0 %% 2 = %s is incorrect value", t.val));
    assert(isInfinity(t.d) && signbit(t.d) == 1, "+0 % 1 has incorrect derivative");

    const y = GDN!1.infinity % x;
    assert(isNaN(y.val) && isNaN(y.d), "inf % anything should be NaN");

    const u = x % e;
    assert(isNaN(u.val) && isNaN(u.d), "anything % -0 should be NaN");

    const i = GDN!1(+0., 1);
    const o = x % i;
    assert(isNaN(o.val) && isNaN(o.d), "anything % +0 should be NaN");

    const p = x % +GDN!1.infinity;
    assert(p.val == 2 && p.d == 3, format("x finite %% inf should be x, not %s", p));

    const h = GDN!2(-3, -2, -1);
    const j = GDN!2(2, 3, 4);
    const k = h % j;
    assert(k.val == -1 && k.d.val == 1 && k.d!2 == 3, format("%s %% %s is not %s", h, j, k));
}

// opBinary(^^)
unittest
{
    const a1 = GDN!1(2, -1);
    const a2 = GDN!1(-2, 3);
    const nz = GDN!1(-0.);
    const pz = GDN!1(+0.);
    const no = GDN!1(-1);
    const po = GDN!1(+1);
    const ni = GDN!1(-real.infinity);
    const pi = GDN!1(+real.infinity);
    const ns = GDN!1(-0.5);
    const ps = GDN!1(0.5);
    const pl = GDN!1(3);

    const q = a1 ^^ a2;
    // dq = 2^-2(-1*-2/2 + 3*ln2) = .25 + .75ln2)
    assert(q.val == 0.25 && q.d == 0.25 + 0.75 * LN2, format("%s ^^ %s is not %s", a1, a2, q));

    const w = a1 ^^ nz;
    // dw = 1(-1*-0/2 + 1ln2 = ln2
    assert(w.val == 1 && w.d == LN2, "anything ^ -0 is incorrect");

    const e = a1 ^^ pz;
    // de = 1(-1*+0/2 + 1ln2) = ln2
    assert(e.val == 1 && w.d == LN2, "anything ^ +0 is incorrect");

    const t = a2 ^^ pi;
    // dt = inf(3*inf/-2 + 1ln-2) = NaN
    assert(t.val == +real.infinity && isNaN(t.d), "x<-1 ^ +inf is incorrect");

    const r = a1 ^^ pi;
    // dr = inf(-1inf/2 + 1ln2) = -inf
    assert(r.val == +real.infinity && r.d == -real.infinity, "x>1 ^ +inf is incorrect");

    const o = ns ^^ pi;
    // do = +0(1*inf/-.5 + 1ln-.5) = NaN
    assert(o.val == +0. && isNaN(o.d), "-1<x<0 ^ +inf is incorrect");

    const u = ps ^^ pi;
    // du = +0(1*inf/.5 + 1ln.5) = NaN
    assert(u.val == +0. && isNaN(u.d), "0<x<1 ^ +inf is incorrect");

    const a = a2 ^^ ni;
    // da = +0(3*-inf/-2 + 1ln-2) = NaN
    assert(a.val == +0. && isNaN(a.d), "x<-1^ -inf is incorrect");

    const p = a1 ^^ ni;
    // dp = +0(-1inf/2 + 1ln2) = NaN
    assert(p.val == +0. && isNaN(p.d), "x>1 ^ -inf is incorrect");

    const d = ns ^^ ni;
    // dd = +inf(1*-inf/-.5 + 1ln-.5) = NaN
    assert(d.val == +real.infinity && isNaN(d.d), "-1<x<0 ^ -inf is incorrect");

    const s = ps ^^ ni;
    // ds = +inf(1*-inf/.5 + 1ln.5) = -inf
    assert(s.val == +real.infinity && s.d == -real.infinity, "0<x<1 ^ -inf is incorrect");

    const h = pi ^^ a1;
    // dh = +inf(1*2/inf + -1ln(inf)) = -inf
    assert(
        h.val == +real.infinity && h.d == -real.infinity, format("+inf ^ (%s) is not %s", a1, h));

    const j = pi ^^ a2;
    // dj = +0(1*-2/inf + 3ln(inf)) = NaN
    assert(j.val == +0. && isNaN(j.d), "+inf ^ y<0 is incorrect");

    const l = ni ^^ pl;
    // dl = -inf(1*3/-inf + 1ln(-inf)) = NaN
    assert(l.val == -real.infinity && isNaN(l.d), format("(-inf) ^ %s is not %s", pl, l));

    const z = ni ^^ a1;
    // dz = inf(1*2/-inf + -1ln(-inf) = NaN
    assert(z.val == +real.infinity && isNaN(z.d), "(-inf) ^ (y>0,not odd) is incorrect");

    const c = ni ^^ no;
    // dc = -0(1*-1/-inf + 1ln(-inf)) = NaN
    assert(c.val == -0. && isNaN(c.d), "(-inf) ^ -1 is incorrect");

    const v = ni ^^ a2;
    // dv = +0(1*-2/-inf + 3ln(-inf)) = NaN
    assert(v.val == +0. && isNaN(v.d), "(-inf) ^ (y<0, not odd) is incorrect");

    const b = no ^^ ni;
    // db = NaN
    assert(isNaN(b.val) && isNaN(b.d), "(-1) ^ -inf is incorrect");

    const n = no ^^ pi;
    // dn = NaN
    assert(isNaN(n.val) && isNaN(n.d), "(-1) ^ inf is incorrect");

    const m = po ^^ ni;
    // dm = NaN
    assert(isNaN(m.val) && isNaN(m.d), "1 ^ -inf is incorrect");

    const f = po ^^ pi;
    // df = NaN
    assert(isNaN(f.val) && isNaN(f.d), "1 ^ inf is incorrect");

    const g = a2 ^^ ns;
    assert(isNaN(g.val) && isNaN(g.d), "(x<0) ^ (finite, nonintegral) is incorrect");

    const i = nz ^^ no;
    // di = -inf(1*-1/-0 + 1ln(-0) = NaN
    assert(i.val == -real.infinity && isNaN(i.d), "(-0) ^ (y<0, odd) is incorrect");

    const k = pz ^^ no;
    // dk = inf(1*(-1)/+0 + 1ln(+0))
    assert(k.val == +real.infinity && k.d == -real.infinity, "(+0) ^ (y<0, odd) is incorrect");

    const y = nz ^^ a2;
    // dy = +inf(1*-2/-0 + 3ln(-0)) = NaN
    assert(y.val == +real.infinity && isNaN(y.d), "(-0) ^ (y<0, odd) is incorrect");

    const qq = pz ^^ a2;
    // dqq = inf(1*(-2)/+0 + 3ln(+0)) = -inf
    assert(qq.val == +real.infinity && qq.d == -real.infinity, "(+0) ^ (y<0. odd) is incorrect");

    const qw = nz ^^ pl;
    // dqw = -0(1*3/-0 + 1ln(-0)) = NaN
    assert(qw.val == -0. && isNaN(qw.d), "(-0) ^ (y>0, odd) is incorrect");

    const qe = pz ^^ pl;
    // dqe = +0*1*3 = +0
    assert(qe.val == +0. && qe.d == +0., "(+0) ^ (y>0, odd) is incorrect");

    const qr = nz ^^ a1;
    // dqr = +0(1*2/-0 + -1ln(-0)) = NaN
    assert(qr.val == +0. && isNaN(qr.d), "(-0) ^ (y>0, not odd) is incorrect");

    const qt = pz ^^ a1;
    // dqt = +0*1*2 = +0
    assert(qt.val == +0. && qt.d == +0., "(+0) ^ (y>0, not odd) is incorrect");

    const qy = GDN!1(2, real.infinity);
    const qu = po ^^ qy;
    // dqu = 1(1*2/1 + inf*ln1) = NaN
    assert(qu.val == 1 && isNaN(qu.d), "1 ^ <2, inf> is incorrect");

    const qi = GDN!1(real.infinity, -real.infinity);
    const qo = qi ^^ a1;
    // dqo = inf(-inf*2/inf + -3ln(inf) = NaN
    assert(qo.val == real.infinity && isNaN(qo.d), "<inf, -inf> ^ <2, -3> is incorrect");

    const qp = GDN!2(1, -1, 0);
    const qa = GDN!2(-2, 2, 3);
    const qs = qp ^^ qa;
    // qs = 1
    // dqs = 1(-1*-2/1 + 2ln1) = 2
    // <dqs,d2dqs> = <1,2>*(<-1,0>*<-2,2>/<1,-1> + <2,3>ln<1,-1>)
    //       <-1,0>*<-2,2>/<1,-1> + <2,3>ln<1,-1>
    //             <-1,0>*<-2,2>/<1,-1>
    //                   <-1,0>*<-2,2> = <2,0*2 + 2*-1> = <2,-2>
    //                = <2,-2>/<1,-1> = <2,(-2*1 - 2*-1)/1^2>
    //                = <2,0>
    //             <2,3>ln<1,-1>
    //                   ln<1,-1> = <0,-1>
    //                = <2,3><0,-1> = <0,3*0 + 2*-1>
    //                = <0,-2>
    //         = <2,0> + <0,-2>
    //         = <2,-2>
    //    = <1,2><2,-2> = <2,2*2 + 1*-2>
    //    = <2,2>
    assert(qs.val == 1 && qs.d == 2 && qs.d!2 == 2, format("%s ^ %s is not %s", qp, qa, qs));

    const qd = 2 ^^ GDN!1(0, -real.infinity);
    assert(qd.val == 1 && qd.d == -real.infinity);

    const qg = GDN!1(0, 1) ^^ GDN!1(2, 0);
    assert(qg == 0 && qg.d == 0);
    // f' = g^(h-1)g'h + g^h*h'ln(g)
    // <0,1> ^ <2,0>, <0,1>*<0,1> = <0,1*0 + 0*1> = <0,0>
    // f = 0
    // f' = 0^1*1*2 + 0^2*0*ln(0)
}

// toHash
unittest
{
    auto q = GDN!1(0.1L, 0.2L);
    auto w = GDN!1(0.1L, 0.0L);
    assert(q.toHash() == w.toHash(), "toHash isn't working correctly");
}

// toString
unittest
{
    const x = GDN!1(0, -1).toString;
    assert(x == "+0 + -1dx", format("GDN!1(0,-1).toString != '%s'", x));

    const y = GDN!2(-1, 1, -2).toString(0);
    assert(y == "-1 + 1dx + -2(dx)²", format("GDN!2(-1, 1, -2).toString(0) != '%s'", y));
}


private pure @safe
{
    // Formats a number, either the value or the value of a derivative, into a string
    string fmtNum(in real num)
    {
        if (isNaN(num))
            return "NaN";
        if (isInfinity(num))
            return (signbit(num) == 0 ? "" : "-") ~ "\u221E";
        if (num == 0)
            return signbit(num) == 0 ? "+0" : "-0";
        return format("%g", num);
    }

    unittest
    {
        string num;

        num = fmtNum(1);
        assert(num == "1", "formatted 1 incorrectly");

        num = fmtNum(-23.4);
        assert(num == "-23.4", "formatted -23.4 incorrectly");

        num = fmtNum(+0.);
        assert(num == "+0", "formatted +0 incorrectly");

        num = fmtNum(-0.);
        assert(num == "-0", "formatted -0 incorrectly");

        num = fmtNum(real.infinity);
        assert(num == "∞", "formatted +∞ incorrectly");

        num = fmtNum(-real.infinity);
        assert(num == "-∞", "formatted -∞ incorrectly");

        num = fmtNum(real.nan);
        assert(num == "NaN", "formatted NaN incorrectly");

        num = fmtNum(-real.nan);
        assert(num == "NaN", "formatted -NaN incorrectly");
    }

    // Formats the derivative identifier of a given order
    string fmtDerivSuffix(in ulong ord)
    {
        switch (ord)
        {
        case 0:
            return "";
        case 1:
            return "dx";
        default:
            return format("(dx)%s", fmtDerivOrd(ord));
        }
    }

    unittest
    {
        assert(fmtDerivSuffix(0) == "", "fmtDerivSuffix(0) is incorrect");
        assert(fmtDerivSuffix(1) == "dx", "fmtDerivSuffix(1) is incorrect");
        assert(fmtDerivSuffix(3) == "(dx)³", "fmtDerivSuffix(3) is incorrect");
    }

    // Formats the order of a derivative into a string
    nothrow string fmtDerivOrd(in ulong ord)
    {
        switch (ord)
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
            return fmtDerivOrd(ord / 10) ~ fmtDerivOrd(ord % 10);
        }
    }

    unittest
    {
        assert(fmtDerivOrd(0) == "⁰", "fmtDerivOrd(0) should be ⁰");
        assert(fmtDerivOrd(1) == "¹", "fmtDerivOrd(1) should be ¹");
        assert(fmtDerivOrd(2) == "²", "fmtDerivOrd(2) should be ²");
        assert(fmtDerivOrd(3) == "³", "fmtDerivOrd(3) should be ³");
        assert(fmtDerivOrd(4) == "⁴", "fmtDerivOrd(4) should be ⁴");
        assert(fmtDerivOrd(5) == "⁵", "fmtDerivOrd(5) should be ⁵");
        assert(fmtDerivOrd(6) == "⁶", "fmtDerivOrd(6) should be ⁶");
        assert(fmtDerivOrd(7) == "⁷", "fmtDerivOrd(7) should be ⁷");
        assert(fmtDerivOrd(8) == "⁸", "fmtDerivOrd(8) should be ⁸");
        assert(fmtDerivOrd(9) == "⁹", "fmtDerivOrd(9) should be ⁹");
        assert(fmtDerivOrd(10) == "¹⁰", "fmtDerivOrd(10) should be ¹⁰");
        assert(fmtDerivOrd(25) == "²⁵", "fmtDerivOrd(25) should be ²⁵");
        assert(fmtDerivOrd(3000) == "³⁰⁰⁰", "fmtDerivOrd(3000) should be ³⁰⁰⁰");
    }
}
