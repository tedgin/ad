/// It extends `std.math.traits` to support `GDN` objects.
module ad.math.traits;

static import std.math.traits;

import std.traits:
    fullyQualifiedName, isFloatingPoint, isImplicitlyConvertible, isIntegral, TemplateOf;

import ad.core;
import ad.math.internal: dirac;


/**
 * Determines if a Given type is a `GDN`.
 *
 * Params:
 *   T = the type to test
 *
 * Returns:
 *   `true` if `T` is a `GDN`, otherwise `false`.
 */
enum bool isGDN(T) = fullyQualifiedName!(TemplateOf!T) == "ad.core.GDN";


///
unittest
{
    static assert(isGDN!(GDN!1));
    static assert(!isGDN!real);
}


/**
 * Determines if a Given type is a `GDN` or implicitly convertible to `real`.
 *
 * Params:
 *   T = the type to test
 *
 * Returns:
 *   `true` if `T` is a `GDN` or implicitly convertible to `real`, otherwise `false`.
 */
enum bool isGDNOrReal(T) = isGDN!T || isImplicitlyConvertible!(T, real);


///
unittest
{
    static assert(isGDNOrReal!(GDN!2));
    static assert(isGDNOrReal!double);
}


// Converts a GDN or something implicitly convertible to a real to a real.
package pragma(inline, true) pure nothrow @nogc @safe real asReal(F)(in F f) if (isGDNOrReal!F)
{
    static if (isGDN!F)
        return f.val;
    else
        return f;
}

unittest
{
    assert(asReal(GDN!2(3)) is 3.0L);
    assert(asReal(1) is 1.0L);
}


/**
 * Checks if the given `GDN` object has a finite value.
 *
 * Params:
 *   Deg = The degree of the `GDN` object.
 *   f = The `GDN` object to check.
 *
 * Returns:
 *   `true` if the value of the `GDN` object is finite, `false` otherwise.
 */
pure nothrow @nogc @safe bool isFinite(ulong Deg)(in GDN!Deg f)
{
    return std.math.traits.isFinite(f.val);
}

///
unittest
{
    assert(isFinite(GDN!1(1)));
}


/**
 * Checks if two `GDN` objects have the same binary representation.
 *
 * Params:
 *   FDeg = The degree of the first `GDN` object.
 *   GDeg = The degree of the second `GDN` object.
 *   f = The first `GDN` object to compare.
 *   g = The second `GDN` object to compare.
 *
 * Returns:
 *   `true` if the two `GDN` objects are identical, `false` otherwise.
 */
pure nothrow @nogc @safe bool isIdentical(ulong FDeg, ulong GDeg)(in GDN!FDeg f, in GDN!GDeg g)
{
    static if (FDeg == GDeg)
    {
        static if (FDeg == 1)
        {
            return std.math.traits.isIdentical(f.val, g.val)
                && std.math.traits.isIdentical(f.d, g.d);
        }
        else
        {
            return std.math.traits.isIdentical(f.val, g.val) && isIdentical(f.d, g.d);
        }
    }
    else
    {
        return false;
    }
}

///
unittest
{
    assert(isIdentical(GDN!1(1), GDN!1(1)));
    assert(!isIdentical(GDN!1(1), GDN!1(2)));
    assert(!isIdentical(GDN!1(1, 2), GDN!1(1, 1)));
    assert(!isIdentical(GDN!1(1), GDN!2(1)));
}


/**
 * Checks if the given `GDN` object has an infinite value.
 *
 * Params:
 *   Deg = The degree of the `GDN` object.
 *   f = The `GDN` object to check.
 *
 * Returns:
 *   `true` if the value of the `GDN` object is infinite, `false` otherwise.
 */
pure nothrow @nogc @safe bool isInfinity(ulong Deg)(in GDN!Deg f)
{
    return std.math.traits.isInfinity(f.val);
}

///
unittest
{
    assert(isInfinity(GDN!1(real.infinity)));
}


/**
 * This function determines whether the value of the given `GDN` object is `NaN`.
 *
 * Params:
 *   Deg = the degree of the `GDN` object
 *   f = the `GDN` object to check
 *
 * Returns:
 *   `true` if the value of the `GDN` object is `NaN`, `false` otherwise.
 */
pure nothrow @nogc @safe bool isNaN(ulong Deg)(in GDN!Deg f)
{
    return std.math.traits.isNaN(f.val);
}

///
unittest
{
    assert(isNaN(GDN!1.nan));
}


/**
 * Checks if the given `GDN` object has a normal value, i.e., it is finite, non-zero, and not
 * subnormal.
 *
 * Params:
 *   Deg = The degree of the `GDN` object.
 *   f = The `GDN` object to check.
 *
 * Returns:
 *   `true` if the value of the `GDN` object is normal, `false` otherwise.
 */
pure nothrow @nogc @safe bool isNormal(ulong Deg)(in GDN!Deg f)
{
    return std.math.traits.isNormal(f.val);
}

///
unittest
{
    assert(isNormal(GDN!1(1)));
}


/**
 * Checks if the given `GDN` object is a power of 2.
 *
 * Params:
 *   Deg = The degree of the `GDN` object.
 *   f = The `GDN` object to check.
 *
 * Returns:
 *   `true` if the value of the `GDN` object is a power of 2, `false` otherwise.
 */
pure nothrow @nogc @safe bool isPowerOf2(ulong Deg)(in GDN!Deg f)
{
    return std.math.traits.isPowerOf2(f.val);
}

///
unittest
{
    assert(isPowerOf2(GDN!1(1)));
}


/**
 * Checks if the given `GDN` object has a subnormal value.
 *
 * Params:
 *   Deg = The degree of the `GDN` object.
 *   f = The `GDN` object to check.
 *
 * Returns:
 *   `true` if the value of the `GDN` object is subnormal, `false` otherwise.
 */
pure nothrow @nogc @safe bool isSubnormal(ulong Deg)(in GDN!Deg f)
{
    return std.math.traits.isSubnormal(f.val);
}

///
unittest
{
    assert(isSubnormal(GDN!1(real.min_normal / 2)));
}


/**
 * Checks if the sign bit of the value of a given `GDN` object is set.
 *
 * Params:
 *   Deg = The degree of the `GDN` object.
 *   f = The `GDN` object to check.
 *
 * Returns:
 *   `1` if the sign bit of the `GDN` object's value is set, `0` otherwise.
 */
pure nothrow @nogc @safe int signbit(ulong Deg)(in GDN!Deg f)
{
    return std.math.traits.signbit(f.val);
}

///
unittest
{
    assert(signbit(GDN!1(-1.0)) == 1);
}


/**
 * This function makes `to` have the same sign as `from`.
 *
 * Params:
 *   G = GDN or implicitly convertible to real
 *   F = a floating-point type
 *   I = an integral type
 *   TDeg = the degree of the `GDN` object to change the sign of
 *   FDeg = the degree of the `GDN` object to copy the sign from
 *   to = the  value to change the sign of
 *   from = the value to copy the sign from
 *
 * Returns:
 *   `to` with the same sign as `from`
 */
pure nothrow @nogc @safe GDN!TDeg copysign(G, ulong TDeg)(in GDN!TDeg to, in G from)
if (isGDNOrReal!G) {
    return GDN!TDeg(std.math.traits.copysign(to.val, asReal(from)), to.d);
}

/// ditto
pure nothrow @nogc @safe F copysign(F, ulong FDeg)(in F to, in GDN!FDeg from) if (isFloatingPoint!F)
{
    return std.math.traits.copysign(to, from.val);
}

/// ditto
pure nothrow @nogc @safe real copysign(I, ulong FDeg)(in I to, in GDN!FDeg from) if (isIntegral!I)
{
    return std.math.traits.copysign(to, from.val);
}

///
unittest
{
    assert(isIdentical(copysign(GDN!1(-1), GDN!1(-2)), GDN!1(-1)));
    assert(isIdentical(copysign(GDN!1(-3), 4.), GDN!1(3)));
    assert(copysign(5., GDN!1(-6.)) is -5.);
    assert(copysign(7, GDN!1(8)) is 7);
}

unittest
{
    assert(isIdentical(copysign(GDN!3(1), GDN!1(-1)), GDN!3(-1)));
}


/**
 * This function computes the sign of a `GDN` Object.
 *
 * If $(MATH f(x) = sgn(g(x))), then $(MATH f' = 2𝛿(g)g'), where $(MATH 𝛿) is the Dirac delta
 * function.
 *
 * To be in agreement with `std.math.traits.sgn`, the sign of  $(MATH sgn(±0) = ±0).
 *
 * Params:
 *   Deg = the degree of the `GDN` object
 *   g = the `GDN` object to compute the sign of
 *
 * Returns:
 *   the sign of the `GDN` object
 */
pure nothrow @nogc @safe GDN!Deg sgn(ulong Deg)(in GDN!Deg g)
{
    return GDN!Deg(std.math.traits.sgn(g.val), 2*dirac(g.reduce())*g.d);
}

///
unittest
{
    assert(isIdentical(sgn(GDN!1(-2)), GDN!1(-1, 0)));
    assert(isIdentical(sgn(GDN!1(0)), GDN!1(0, real.infinity)));
}

unittest
{
    assert(isIdentical(sgn(GDN!1()), GDN!1.nan));
    assert(isIdentical(sgn(GDN!1(0, real.nan)), GDN!1(0, real.nan)));
    assert(isIdentical(sgn(GDN!2(1.0L, 2.0L, real.nan)), GDN!2(1.0L, 0.0L, real.nan)));
    assert(isIdentical(sgn(GDN!1(real.infinity, 0)), GDN!1(1, 0)));
    assert(isIdentical(sgn(GDN!1(-real.infinity, 0)), GDN!1(-1, 0)));
    assert(isIdentical(sgn(GDN!1(-1)), GDN!1(-1, 0)));
}
