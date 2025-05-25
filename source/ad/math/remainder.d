/// It extends the `std.math.remainder` to support `GDN` objects.
module ad.math.remainder;

static import std.math.remainder;

import ad.core;
import ad.math.internal: areAll, CommonGDN, isGDN, isGDNOrReal, isOne, trunc;


/**
 * Returns the remainder of `g` divided by `h`.
 *
 * If $(MATH f(x) = g(x) (mod h(x))), then $(MATH f' = g' - (g - f)h'/h - δ(f)(g'h - gh')/h), where
 * $(MATH δ) is the Dirac Delta function.
 *
 * Params:
 *   G = the type of g, either `GDN` or `real`.
 *   H = the type of h, either `GDN` or `real`.
 *   g = the dividend.
 *   h = the divisor.
 *
 * Returns:
 *   The remainder of `g` divided by `h`.
 */
nothrow @nogc @safe
CommonGDN!(G, H) fmod(G, H)(in G g, in H h) if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H))
{
	return g % h;
}

///
unittest
{
	const f = fmod(GDN!1(5), GDN!1(3));
	assert(f is GDN!1(2,0));
}


/**
 * Breaks g into an integer and a fraction, each with the same sign as g.
 *
 * `f = modf(g, i)` can be expressed mathematically as follows. When $MATH(g ≥ 0), $(MATH i = ⌊g⌋),
 * and when $(MATH g < 0), $(MATH i = ⌈g⌉). This is the mathematical definition of `trunc(g)`.
 * $(MATH f = g - i).
 *
 * Params:
 *   Deg = the degree of g
 *   g = the GDN object to break into an integer and a fraction.
 *   i = the integer part of g.
 *
 * Returns:
 *   The fractional part of g.
 */
nothrow @nogc @safe GDN!Deg modf(ulong Deg)(in GDN!Deg g, out GDN!Deg i)
{
	i = trunc(g);
	return g - i;
}

///
unittest
{
	import std.math: isClose;

	const g = GDN!1(3.14159);
	GDN!1 i;
	const f = modf(GDN!1(g), i);
	assert(i is GDN!1(3, 0));
	assert(isClose(f.val, 0.14159));
	assert(f.d == 1);
}


// TODO: implement
nothrow @nogc @safe
CommonGDN!(G, H) remainder(G, H)(in G g, in H h)
if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H));


// TODO: implement
nothrow @nogc @safe
CommonGDN!(G, H) remquo(G, H)(in G g, in H h, out int n)
if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H));
