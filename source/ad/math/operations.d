/// This module extends `std.math.operations` to support `GDN` objects.
module ad.math.operations;

static import std.math.operations;

import std.algorithm: min;
import std.math.hardware: ieeeFlags, resetIeeeFlags;
import std.range: ElementType, empty, front, isInputRange, popFront;
import std.traits: isImplicitlyConvertible, Select;

import ad.core;
import ad.math.traits:
    areAll, asGDN, asReal, CommonGDN, isGDN, isGDNOrReal, isNaN, isOne, sgn, signbit;


/// The default relative difference for operations
enum real DEFAULT_REL_DIFF = 10.0L ^^ -((real.dig + 1)/2 + 1);


/**
 * Defines a total order on all GDN objects. Id one of the objects has a lower degree than the
 * other, it is promoted a GDN of the higher degree with all added derivatives being `0`. If one of
 * the arguments is a floating point number, it is promoted to a constant GDN with the same degree
 * as the other term. The result is `cmp(x.val, y.val)` unless it is `0`. Otherwise, the result is
 * `cmp(x.d, y.d)` unless it is also `0`, then the result is `cmp(x.d!2, y.d!2)` and so on. If all
 * derivatives are the same, the result is `0`.
 *
 * Params:
 *   X = GDN or real
 *   Y = GDN or real
 *   x = the object being compared
 *   y = the object being compared to
 *
 * Returns:
 *   a negative value if `x` precedes `y`, `0` if `x` and `y` are identical, and a positive value
 *   otherwise.
 */
pure nothrow @nogc @safe int cmp(X, Y)(in X x, in Y y)
if (isOne!(isGDN, X, Y) && areAll!(isGDNOrReal, X, Y)) {
    static if (isImplicitlyConvertible!(X, real))
        alias Deg = Y.DEGREE;
    else static if (isImplicitlyConvertible!(Y, real))
        alias Deg = X.DEGREE;
    else
        alias Deg = Select!(X.DEGREE > Y.DEGREE, X.DEGREE, Y.DEGREE);

    return cmp_impl(asGDN!Deg(x), asGDN!Deg(y));
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
if (areAll!(isGDNOrReal, T, U) && isOne!(isGDN, T, U)) {
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
if (areAll!(isInputRange, T, U) && areAll!(isGDNOrReal, ElementType!T, ElementType!U)) {
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


/**
 * Create a quiet NaN GDN, storing an integer as payload.
 *
 * Params:
 *   Deg = the degree of the `GDN` to create
 *   payload = the integral payload
 *
 * Returns:
 *   It returns a `GDN` object whose value is a NaN with the given payload.
 */
pure nothrow @nogc @safe GDN!Deg NaN(ulong Deg)(in ulong payload)
{
    return GDN!Deg(std.math.operations.NaN(payload));
}

///
unittest
{
    import std.math: getNaNPayload;

    assert(getNaNPayload(NaN!2(3).val) == 3);
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
    assert(getNaNPayload(NaN!1(1)) == 1);
}


/**
 * Returns the positive difference between `g` and `h`, i.e. $(MATH max{g-h, 0}). The degree of the
 * result with will be the lesser of `GDeg` and `HDeg`. If either input is a real number, it will be
 * converted to a constant `GDN` with the same degree as the other input.
 *
 * If $(MATH f(x) = 0) when $(MATH g(x) < h(x)), and $(MATH f(x) = g(x) - h(x)) when
 * $(MATH g(x) > h(x)), then $(MATH f' = 0) when $(MATH g < h), and $(MATH f' = g' - h') when
 * $(MATH g > h).
 *
 * Params:
 *   G = the type of `g` either a `GDN` or a `real`
 *   H = the type of `h` either a `GDN` or a `real`
 *   g = the minuend
 *   h = the subtrahend
 *
 * Returns:
 *   the resulting generalized dual number
 *
 * Bug:
 *   $(MATH f') may not be correct when $(MATH g = h).
 */
pure nothrow @nogc @safe
CommonGDN!(G, H) fdim(G, H)(in G g, in H h) if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H))
{
    alias Deg = typeof(return).DEGREE;

    const gmh = asGDN!Deg(g) - asGDN!Deg(h);
    if (signbit(gmh) == 1 && !isNaN(gmh)) return GDN!Deg.zero;
    return gmh;
}

///
unittest
{
    assert(fdim(GDN!1(2), GDN!2(2)) is GDN!1(0, 0));
}

unittest
{
    const g = fdim(GDN!1(0), -2);
    assert(g is GDN!1(2));

    const f = fdim(0, GDN!2(real.infinity));
    assert(f is GDN!2.zero,);

    assert(fdim(GDN!1(2) ,GDN!1(0)) is GDN!1(2, 0));
    assert(fdim(GDN!2(-2), GDN!2(0)) is GDN!2.zero);
    assert(isNaN(fdim(GDN!1.nan, GDN!1(2))));
}


/**
 * Computes $(MATH gh + i). The degree of the result will be the least of the degrees of `G`, `H`,
 * and `I`. If any of the inputs is a real number, it will be converted to a constant `GDN` with the
 * same degree as the common degree of the other inputs.
 *
 * If $(MATH f(x) = g(x)h(x) + i(x)), then $(MATH f' = g'h + gh' + i').
 *
 * Params:
 *   G = the type of `g`, either a `GDN` or a `real`
 *   H = the type of `h`, either a `GDN` or a `real`
 *   I = the type of `i`, either a `GDN` or a `real`
 *   g = the multiplicand
 *   h = the multiplier
 *   i = the addend
 *
 * Returns:
 *   the resulting generalized dual number
 */
pure nothrow @nogc @safe CommonGDN!(G, H, I) fma(G, H, I)(in G g, in H h, in I i)
if (isOne!(isGDN, G, H, I) && areAll!(isGDNOrReal, G, H, I)) {
    alias Deg = typeof(return).DEGREE;

    const gg = asGDN!Deg(g);
    const hh = asGDN!Deg(h);
    const ii = asGDN!Deg(i);

    return GDN!Deg(
        std.math.operations.fma(gg.val, hh.val, ii.val), gg.d*hh.val + gg.val*hh.d + ii.d);
}

///
unittest
{
    assert(fma(GDN!1(2), GDN!2(3), 4) is GDN!1(10, 5));
}


/**
 * Computes $(MATH max{g, h}). The degree of the result will be the least of the degrees of `G` and
 * `H`. If either input is a real number, it will be converted to a constant `GDN` with the same
 * degree as the other input.
 *
 * If $(MATH f(x) = max{g(x), h(x)}), then $(MATH f' = g')$ when $(MATH g > h)$, and
 * $(MATH f' = h') when $(MATH g < h).
 *
 * Params:
 *   G = the type of `g`, either a `GDN` or a `real`
 *   H = the type of `h`, either a `GDN` or a `real`
 *   g = the first argument
 *   h = the second argument
 *
 * Returns:
 *   the resulting generalized dual number
 *
 * Bug:
 *   $(MATH f') may not be correct when $(MATH g = h).
 */
pure nothrow @nogc @safe
CommonGDN!(G, H) fmax(G, H)(in G g, in H h) if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H))
{
    alias Deg = typeof(return).DEGREE;

    const gg = asGDN!Deg(g);
    const hh = asGDN!Deg(h);

    if (isNaN(gg) || hh > gg) {
        return hh;
    } else {
        return gg;
    }
}

///
unittest
{
    assert(fmax(GDN!1(2), GDN!2(3)) is GDN!1(3));
}

unittest
{
    assert(fmax(GDN!1(2), GDN!1(3)) is GDN!1(3));
    assert(fmax(GDN!1(2), 3) is GDN!1(3, 0));
    assert(fmax(2, GDN!1(3)) is GDN!1(3));
}


/**
 * Computes $(MATH min{g, h}). The degree of the result will be the least of the degrees of `G` and
 * `H`. If either input is a real number, it will be converted to a constant `GDN` with the same
 * degree as the other input.
 *
 * If $(MATH f(x) = min{g(x), h(x)}), then $(MATH f' = g')$ when $(MATH g < h)$, and
 * $(MATH f' = h') when $(MATH g > h).
 *
 * Params:
 *   G = the type of `g`, either a `GDN` or a `real`
 *   H = the type of `h`, either a `GDN` or a `real`
 *   g = the first argument
 *   h = the second argument
 *
 * Returns:
 *   the resulting generalized dual number
 *
 * Bug:
 *   $(MATH f') may not be correct when $(MATH g = h).
 */
pure nothrow @nogc @safe
CommonGDN!(G, H) fmin(G, H)(in G g, in H h) if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H))
{
    alias Deg = typeof(return).DEGREE;

    const gg = asGDN!Deg(g);
    const hh = asGDN!Deg(h);

    if (isNaN(gg) || hh < gg) {
        return hh;
    } else {
        return gg;
    }
}

///
unittest
{
    assert(fmin(GDN!1(2), GDN!2(3)) is GDN!1(2));
}

unittest
{
    assert(fmin(GDN!1(2), GDN!1(3)) is GDN!1(2));
    assert(fmin(GDN!1(2), 3) is GDN!1(2));
    assert(fmin(2, GDN!1(3)) is GDN!1(2, 0));
}


/**
 * Computes the next representable value after `g` in the direction of `h`.
 *
 * If $(MATH f(x) = g(x) + sgn(h(x))ε), where $(MATH ε) is an infinitesimal, then $(MATH f' = g').
 *
 * Params:
 *   Deg = the degree a `GDN`
 *   H = the type of `h`, either a `GDN` or a `real`
 *   g = the starting value
 *   h = a value in the target direction
 *
 * Returns:
 *   the value with the next representable value after `g` in the direction of `h`.
 */
pure nothrow @nogc @safe GDN!Deg nextafter(H, ulong Deg)(in GDN!Deg g, in H h) if (isGDNOrReal!H)
{
    const f = std.math.operations.nextafter(g.val, asGDN!Deg(h).val);
    return GDN!Deg(f, std.math.isNaN(f) ? GDN!Deg.mkNaNDeriv() : g.d);
}

/// ditto
pure nothrow @nogc @safe real nextafter(ulong Deg)(in real g, in GDN!Deg h) {
    return std.math.operations.nextafter(g, asReal(h));
}

///
unittest
{
    assert(nextafter(GDN!2(1), GDN!1(2)) is GDN!2(1+real.epsilon));
    assert(nextafter(1, GDN!1(0)) == 1 - real.epsilon/2);
}

unittest
{
    assert(nextafter(GDN!1(2), 2) is GDN!1(2));
    assert(nextafter(GDN!1(1), GDN!2(0)) is GDN!1(1 - real.epsilon/2));
    assert(isNaN(nextafter(GDN!1.nan, GDN!1(2))));
}


/**
 * Computes the largest representable GDN smaller than `g` having the same degree as `g`.
 *
 * If $(MATH f(x) = g(x) - ε), where $(MATH ε) is an infinitesimal, then $(MATH f' = g').
 *
 * Params:
 *   Deg = the degree a `g`
 *   g = the starting value
 *
 * Returns:
 *   a `GDN` object with the same degree as `g`.
 */
pure nothrow @nogc @safe GDN!Deg nextDown(ulong Deg)(in GDN!Deg g)
{
    const f = std.math.operations.nextDown(g.val);
    return GDN!Deg(f, std.math.isNaN(f) ? GDN!Deg.mkNaNDeriv() : g.d);
}

///
unittest
{
    assert(nextDown(GDN!2(1)) is GDN!2(1 - real.epsilon/2));
}

unittest
{
    assert(isNaN(nextDown(GDN!1.nan)));
}


/**
 * Computes the smallest representable GDN larger than `g` having the same degree as `g`.
 *
 * If $(MATH f(x) = g(x) + ε), where $(MATH ε) is an infinitesimal, then $(MATH f' = g').
 *
 * Params:
 *   Deg = the degree a `g`
 *   g = the starting value
 *
 * Returns:
 *   a `GDN` object with the same degree as `g`.
 */
pure nothrow @nogc @safe GDN!Deg nextUp(ulong Deg)(in GDN!Deg g)
{
    const f = std.math.operations.nextUp(g.val);
    return GDN!Deg(f, std.math.isNaN(f) ? GDN!Deg.mkNaNDeriv() : g.d);
}

///
unittest
{
    assert(nextUp(GDN!2(1)) is GDN!2(1 + real.epsilon));
}

unittest
{
    assert(isNaN(nextUp(GDN!1.nan)));
}
