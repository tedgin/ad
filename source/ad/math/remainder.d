/// It extends the `std.math.remainder` to support `GDN` objects.
module ad.math.remainder;

static import std.math.remainder;

import ad.core;
import ad.math.internal: areAll, CommonGDN, isGDN, isGDNOrReal, isOne;


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


// TODO: implement
nothrow @nogc @safe real modf(ulong Deg)(in GDN!Deg g, ref real i);


// TODO: implement
nothrow @nogc @safe
CommonGDN!(G, H) remainder(G, H)(in G g, in H h)
if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H));


// TODO: implement
nothrow @nogc @safe
CommonGDN!(G, H) remquo(G, H)(in G g, in H h, out int n)
if (isOne!(isGDN, G, H) && areAll!(isGDNOrReal, G, H));
