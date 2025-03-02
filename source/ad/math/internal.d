// This module implements functions used internally by the `ad.math` module.
module ad.math.internal;

import std.math: isNaN;

import ad.core;

package pure nothrow @nogc @safe
{
    // Exposes the GDN.dirac method as a function to help with overload resolution.
    pragma(inline, true) GDN!Deg dirac(ulong Deg)(in GDN!Deg g)
    {
        return g.dirac();
    }

    // The Dirac delta function for a `real`.
    real dirac(in real g)
    {
        if (isNaN(g)) {
            return real.nan;
        }

        return g == 0 ? real.infinity : 0;
    }

    unittest
    {
        assert(dirac(0) == real.infinity);
        assert(dirac(1) == 0);
        assert(isNaN(dirac(real.nan)));
    }
}