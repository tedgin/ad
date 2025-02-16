// This module implements the dirac delta function. It is intended for internal use only.
module ad.math.dirac;

import ad.core;

package nothrow pure @nogc @safe
{
    pragma(inline, true) GDN!Deg dirac(ulong Deg)(in GDN!Deg g)
    {
        return g.dirac();
    }

    real dirac(in real g)
    {
        return g == 0 ? real.infinity : 0;
    }
}