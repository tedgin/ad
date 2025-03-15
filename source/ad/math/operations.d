/// This module extends `std.math.operations` to support `GDN` objects.
module ad.math.operations;

static import std.math.operations;

import std.range: ElementType, empty, front, isInputRange, popFront;
import std.traits: isImplicitlyConvertible, Select;

import ad.core;
import ad.math.traits: asReal, isGDN, isGDNOrReal;


/// The default relative difference for operations
enum real DEFAULT_REL_DIFF = 10.0L ^^ -((real.dig + 1)/2 + 1);


/**
 * Defines a total order on all GDN objects. The one of the objects has a lower degree than the
 * other, it is promoted a GDN of the higher degree with all added derivatives being `0`. If one of
 * the arguments is a floating point number, it is promoted to a constant GDN with the same degree
 * as the other term. The result is `cmp(x.val, y.val)` unless it is `0`. Otherwise, the result is
 * `cmp(x.d, y.d)` unless it is also `0`, then the result is `cmp(x.d!2, y.d!2)` and so on. If all
 * derivatives are the same, the result is `0`.
 *
 * Params:
 *   F = a floating point type
 *   XDeg = the degree of `x`
 *   YDeg = the degree of `y`
 *   x = the `GDN` object being compared
 *   y = the `GDN` object being compared to
 *
 * Returns:
 *   a negative value if `x` precedes `y`, `0` if `x` and `y` are identical, and a positive value
 *   otherwise.
 */
pure nothrow @nogc @safe int cmp(ulong XDeg, ulong YDeg)(const GDN!XDeg x, const GDN!YDeg y)
{
    alias Deg = Select!(XDeg > YDeg, XDeg, YDeg);
    return cmp_impl(cast(GDN!Deg) x, cast(GDN!Deg) y);
}

/// ditto
pure nothrow @nogc @safe int cmp(F, ulong XDeg)(const GDN!XDeg x, const(F) y)
if (isImplicitlyConvertible!(F, real)) {
    return cmp_impl(x, GDN!XDeg.mkConst(y));
}

/// ditto
pure nothrow @nogc @safe int cmp(F, ulong YDeg)(const(F) x, const GDN!YDeg y)
if (isImplicitlyConvertible!(F, real)) {
    return cmp_impl(GDN!YDeg.mkConst(x), y);
}

///
unittest
{
    assert(cmp(GDN!1(0, 0), GDN!1(0, 1)) < 0);
    assert(cmp(GDN!2(0, 0, 1), GDN!1(0, 0)) > 0);
    assert(cmp(GDN!1(1, 0), 1.) == 0);
}

unittest
{
    assert(cmp(GDN!3(2), GDN!2(2)) == 0);
    assert(cmp(0., GDN!1(0)) < 0);
}

private pragma(inline, true) pure nothrow @nogc @safe
int cmp_impl(ulong Deg)(in GDN!Deg x, in GDN!Deg y)
{
    const res = std.math.operations.cmp(x.val, y.val);

    static if (Deg == 1)
        alias d_cmp = std.math.operations.cmp;
    else
        alias d_cmp = cmp_impl!(Deg - 1);

    return res != 0 ? res : d_cmp(x.d, y.d);
}

unittest
{
    assert(cmp_impl(GDN!1(0), GDN!1(0)) == 0);
    assert(cmp_impl(GDN!1(-1), GDN!1(0)) < 0);
    assert(cmp_impl(GDN!1(1), GDN!1(0)) > 0);
    assert(cmp_impl(GDN!1(0, -1), GDN!1(0, 0)) < 0);
    assert(cmp_impl(GDN!1(0, 1), GDN!1(0, 0)) > 0);
    assert(cmp_impl(GDN!2(0), GDN!2(0)) == 0);
    assert(cmp_impl(GDN!2(0, -1, 0), GDN!2(0, 0, 0)) < 0);
    assert(cmp_impl(GDN!2(0, 0, 0), GDN!2(0, 0, -1)) > 0);
}


/**
 * To what precision are the values of `x` and `y` equal?
 *
 * Params:
 *   X = GDN or real
 *   Y = GDN or real
 *   x = the `GDN` being compared
 *   y = a `GDN` being compared to
 *
 * Returns:
 *   the number of mantissa bits which are equal in `x` and `y`.
 */
pure nothrow @nogc @safe int feqrel(X, Y)(in X x, in Y y)
if (isGDNOrReal!X && isGDNOrReal!Y && (isGDN!X || isGDN!Y)) {
    return std.math.operations.feqrel(asReal(x), asReal(y));
}

///
unittest
{
    assert(feqrel(GDN!1(2),GDN!1(2)) == real.mant_dig);
}

unittest
{
    assert(feqrel(2., GDN!2(2)) == real.mant_dig);
    assert(feqrel(GDN!1(2), real.nan) == 0);
}


/**
 * Extract an integral payload from a NaN-valued `GDN` object.
 *
 * Params:
 *   Deg = the degree of `f`
 *   f = the `GDN` object carrying the NaN payload
 *
 * Returns:
 *   the payload
 */
pure nothrow @nogc @safe ulong getNaNPayload(ulong Deg)(in GDN!Deg f)
{
    return std.math.operations.getNaNPayload(f.val);
}

///
unittest
{
    import std.math: NaN;

    assert(getNaNPayload(GDN!1(NaN(1))) == 1);
}


/**
 * Determines whether the values of two `GDN` objects or a `GDN` object and a `real` are
 * approximately equal, i.e., with a given maximum relative difference and a maximum absolute
 * difference.
 *
 * Params:
 *   T = a GDN object or real or a range type of such
 *   U = a GDN object or real or a range type of such
 *   lhs = the value being compared
 *   rhs = the value being compared to
 *   maxRelDiff = maximum allowable relative difference, setting to 0 disables this check
 *   maxAbsDiff = maximum absolute difference, setting this to 0 disables this check
 *
 * Returns:
 *   `true` if the two values are approximately equal under either of the criteria.
 *
 *   If either object is a range, and the other is a single value, then `true` is returned if the
 *   single value satisfies either of the criteria for all of the elements in the range.
 *
 *   If both objects are ranges, then a `true` is returned only if they have the same length, and
 *   either of the criteria are met for each pair of elements.
 */
pure nothrow @nogc @safe
bool isClose(T, U)(in T lhs, in U rhs, in real maxRelDiff=DEFAULT_REL_DIFF, in real maxAbsDiff=0)
if (isGDNOrReal!T && isGDNOrReal!U && (isGDN!T || isGDN!U)) {
    return std.math.operations.isClose(asReal(lhs), asReal(rhs), maxRelDiff, maxAbsDiff);
}

/// ditto
pure nothrow @nogc @safe
bool isClose(T, U)(in T lhs, U rhs, in real maxRelDiff=DEFAULT_REL_DIFF, in real maxAbsDiff=0)
if (isGDNOrReal!T && isInputRange!U && isGDNOrReal!(ElementType!U)) {
    for(; !rhs.empty; rhs.popFront()) {
        if (!isClose(lhs, rhs.front, maxRelDiff, maxAbsDiff)) {
            return false;
        }
    }

    return true;
}

/// ditto
pure nothrow @nogc @safe
bool isClose(T, U)(T lhs, in U rhs, in real maxRelDiff=DEFAULT_REL_DIFF, in real maxAbsDiff=0)
if (isInputRange!T && isGDNOrReal!(ElementType!T) && isGDNOrReal!U) {
    for(; !lhs.empty; lhs.popFront()) {
        if (!isClose(lhs.front, rhs, maxRelDiff, maxAbsDiff)) {
            return false;
        }
    }

    return true;
}

/// ditto
pure nothrow @nogc @safe
bool isClose(T, U)(T lhs, U rhs, in real maxRelDiff=DEFAULT_REL_DIFF, in real maxAbsDiff=0)
if (isInputRange!T && isGDNOrReal!(ElementType!T) && isInputRange!U && isGDNOrReal!(ElementType!U))
{
    for(;; lhs.popFront(), rhs.popFront()) {
        if (lhs.empty) return rhs.empty;
        if (rhs.empty) return lhs.empty;
        if (!isClose(lhs.front, rhs.front, maxRelDiff, maxAbsDiff)) return false;
    }
}

///
unittest
{
    assert(isClose(GDN!1(1),GDN!2(0.999_999_999_9)));
    assert(isClose(GDN!1(2), [GDN!1(2), GDN!1(1.999_999_999_9), GDN!1(2.000_000_000_1)]));

    assert(isClose(
        [GDN!1(1), GDN!1(2), GDN!1(3)],
        [GDN!1(0.999_999_999_9), GDN!1(2.000_000_000_1), GDN!1(3)]));
}

unittest
{
    assert(isClose(GDN!1(17.123_456_789), GDN!1(17.123_45), 1e-6));
    assert(isClose(GDN!1(1e-100), GDN!1.zero, 0, 1e-90));

    assert(isClose(GDN!1(0.001), 0.000_999_999_999_99));
    assert(!isClose(GDN!1(17.123_456_789), 17.123_45, 1e-7));
    assert(isClose(GDN!1(1e-10), -1e-10, 0, 1e-9));

    assert(isClose(1_000_000_000.0, GDN!1(999_999_999.9)));
    assert(isClose(17.123_456_789, GDN!1(17.123_45), 1e-6));
    assert(isClose(1e-300, GDN!1(1e-298), 0, 1e-200));

    assert(isClose(GDN!1(2), [2, 1.999_999_999], 1e-8));
    assert(!isClose(0, [GDN!1(1e-100), GDN!1(1e-10)], 0, 1e-90));

    assert(isClose([GDN!1(2), GDN!1(1.999_999_999_9), GDN!1(2.000_000_000_1)], GDN!1(2)));
    assert(isClose([2, 1.999_999_999], GDN!1(2), 1e-8));
    assert(!isClose([1e-100, 1e-10], GDN!1(0), 0, 1e-90));

    assert(!isClose([GDN!1(1), GDN!1(2)], [0.999_999_999, 2.000_000_001, 3]));
    assert(!isClose([1, 2, 3], [GDN!1(0.999_999_999), GDN!1(2.000_000_001)]));
    assert(isClose([GDN!1(1), GDN!1(2)], [0.999_999_999, 2.000_000_001], 1e-9));
    assert(isClose([GDN!1(0), GDN!1(2)], [1e-100, 2.000_000_001], 1e-9, 1e-90));
}


// TODO: implement NaN
// TODO: implement nextafter
// TODO: implement nextDown
// TODO: implement nextUp

// TODO: implement fdim
// TODO: implement fma
// TODO: implement fmax
// TODO: implement fmin
