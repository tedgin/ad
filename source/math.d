/**
This module extends the std.math library to supporting GenDualNum objects.
*/
module ad.math;

import std.math;

public import ad.core;

/**
This function computes the absolute value of the argument. It is analogous to
std.math.abs().

The derivative of the abs(x) is x*x'/abs(x).
 
Params:
    x = the argument
*/
nothrow pure @safe real abs(in real x)
{
    return std.math.abs(x);
}

/// ditto
nothrow pure @safe GenDualNum!Order abs(ulong Order)(in GenDualNum!Order x)
{
    alias PN = GenDualNum!Order;

    const dx = x == 0 ? PN.DerivType!1.nan : x.d * sgn(x.reduce()); 
    return PN(abs(x.val), dx);
}

unittest 
{
    assert(abs(GenDualNum!1()).same(GenDualNum!1(real.nan, real.nan)));

    // XXX - derivSeq isn't implemented due to the following bug 
    // https: //issues.dlang.org/show_bug.cgi?id=22621 is fixed.
    // assert(abs(GenDualNum!1.zero).same(derivSeq(0.0L, real.nan)));
    // assert(abs(derivSeq(2.0L, 2.0L)).same(derivSeq(2.0L, 2.0L)));
    // assert(abs(GenDualNum!1.var(-1)).same(derivSeq(1.0L, -1.0L)));
    
    assert(abs(GenDualNum!1.infinity).same(GenDualNum!1.infinity));
    assert(abs(-GenDualNum!1.infinity).same(GenDualNum!1.infinity));
}

/** 
This function computes the cosine of the argument. It is analogous to 
std.math.cos().
 
Params:
    x = the argument
*/
nothrow pure @safe real cos(in real x)
{
    return std.math.cos(x);
}

/// ditto
nothrow pure @safe GenDualNum!Order cos(ulong Order)(in GenDualNum!Order x)
{
    return GenDualNum!Order(cos(x.val), -sin(x.reduce()) * x.d);
}

unittest 
{
    assert(cos(GenDualNum!1()).same(GenDualNum!1(real.nan, real.nan)));

    // XXX - derivSeq isn't implemented due to the following bug 
    // https: //issues.dlang.org/show_bug.cgi?id=22621 is fixed.
    // assert(cos(GenDualNum!1.zero).same(derivSeq(1.0L, 0.0L)));
    // assert(cos(GenDualNum!1.var(PI_2)).same(derivSeq(cos(PI_2), -1.0L)));
    
    assert(cos(GenDualNum!1.infinity).same(GenDualNum!1.nan));
    assert(cos(-GenDualNum!1.infinity).same(GenDualNum!1.nan));
}

/**
This function raises e to a given power. It is analogous to std.math.exp().

Params:
    x = the power e is raised to.
*/
nothrow pure @safe real exp(in real x)
{
    return std.math.exp(x);
}

/// ditto
nothrow pure @safe GenDualNum!Order exp(ulong Order)(in GenDualNum!Order x)
{
    return GenDualNum!Order(exp(x.val), x.d * exp(x.reduce()));
}

unittest 
{
    assert(exp(GenDualNum!1.nan).same(GenDualNum!1.nan));
    
    assert(exp(GenDualNum!1.var(0)).same(GenDualNum!1.var(1)));

    // XXX - derivSeq isn't implemented due to the following bug 
    // https: //issues.dlang.org/show_bug.cgi?id=22621 is fixed.
    // assert(exp(GenDualNum!1.infinity).same(derivSeq(real.infinity, real.infinity)));

    assert(exp(-GenDualNum!1.infinity).same(GenDualNum!1.zero));
}

/**
This function determines whether its argument is a NaN.  It is analogous to 
std.math.isNaN().
  
Params:
    x = the argument
*/
nothrow pure @safe bool isNaN(in real x)
{
    return std.math.isNaN(x);
}

/// ditto
nothrow pure @safe bool isNaN(ulong Order)(in GenDualNum!Order x) 
{
    return std.math.isNaN(x.val);
}

unittest 
{
    assert(isNaN(GenDualNum!1.nan));
}

/**
This function computes the natural logarithm of its argument. It is analogous
to std.math.log().
 
Params:
    x = the argument
*/
nothrow pure @safe real log(in real x)
{
    return std.math.log(x);
}

/// ditto
nothrow pure @safe GenDualNum!Order log(ulong Order)(in GenDualNum!Order x)
{
    return x.log();
}

/**
This function computes the raises a given number to a given power. It is 
analogous to std.math.pow().

Params:
    x = the base
    y = the exponent
*/
nothrow pure @safe real pow(in real x, real y)
{
    return std.math.pow(x, y);
}

/// ditto
nothrow pure @safe GenDualNum!Order pow(ulong Order)(in GenDualNum!Order x, in real y)
{
    return x ^^ y;
}

/// ditto
nothrow pure @safe GenDualNum!Order pow(ulong Order)(in real x, in GenDualNum!Order y)
{
    return x ^^ y;
}

/// ditto
nothrow pure @safe GenDualNum!Order pow(ulong Order)(in GenDualNum!Order x, in GenDualNum!Order y)
{
    return x ^^ y;
}

/**
This function computes the sign of the argument. It is analogous to 
std.math.sgn().
 
The derivate of sgn(x) evaluated at 0 is undefined or NaN. Otherwise it is 0 * 
x'.
 
Params:
    x = the argument
*/
nothrow pure @safe real sgn(in real x)
{ 
    return std.math.sgn(x); 
}

/// ditto
nothrow pure @safe GenDualNum!Order sgn(ulong Order)(in GenDualNum!Order x)
{
    alias PN = GenDualNum!Order;

    if (isNaN(x))
        return PN.nan;
    else if (x == 0)
        return PN(0, PN.DerivType!1.nan);
    else
        return PN(std.math.sgn(x.val), PN.DerivType!1.zero);
}

/// ditto
nothrow pure @safe GenDualNum!Order sgn(ulong Order : 1)(in GenDualNum!Order x)
{
    alias PN = GenDualNum!Order;

    if (isNaN(x))
        return PN.nan;
    else if (x == 0)
        return PN(0, real.nan);
    else
        return PN(std.math.sgn(x.val), 0);
}

unittest 
{
    import std.format;

    assert(
        sgn(GenDualNum!1()).same(GenDualNum!1.nan), 
        format("sgn(%s) is not %s", sgn(GenDualNum!1()), GenDualNum!1.nan));

    assert(sgn(GenDualNum!1(0, real.nan)).same(GenDualNum!1(0, real.nan)));

    // XXX - derivSeq isn't implemented due to the following bug 
    // https: //issues.dlang.org/show_bug.cgi?id=22621 is fixed.
    // assert(sgn(derivSeq(1.0L, 2.0L, real.nan)).same(derivSeq(1.0L, 0.0L, real.nan)));

    assert(sgn(GenDualNum!1(real.infinity, 0)).same(GenDualNum!1.one));
    assert(sgn(GenDualNum!1(-real.infinity, 0)).same(-GenDualNum!1.one));
}

/**
This function computes the square root of its argument. It is analogous to 
std.math.sqrt().
 
Params:
    x = the argument
*/
// XXX - This should be needed, because of implicit conversion
nothrow pure @safe real sqrt(in real x) 
{
    return std.math.sqrt(x);
}

/// ditto
nothrow pure @safe GenDualNum!Order sqrt(ulong Order)(in GenDualNum!Order x)
{
    alias PN = GenDualNum!Order;

    if (x == 0) 
        return PN(0, PN.DerivType!1.nan);
    else
        return PN(sqrt(x.val), x.d / (2 * sqrt(x.reduce())));        
}

unittest 
{
    // XXX - derivSeq isn't implemented due to the following bug 
    // https: //issues.dlang.org/show_bug.cgi?id=22621 is fixed.
    // assert(sqrt(GenDualNum!1.var(1)).same(derivSeq(1.0L, 0.5L)));
    // 
    // assert(sqrt(GenDualNum!1.var(0)).same(derivSeq(0.0L, real.nan)));
    assert(sqrt(GenDualNum!1.var(-1)).same(GenDualNum!1.nan));

    // XXX - derivSeq isn't implemented due to the following bug 
    // https: //issues.dlang.org/show_bug.cgi?id=22621 is fixed.
    // assert(sqrt(GenDualNum!1.infinity).same(derivSeq(real.infinity, 0.0L)));

    assert(sqrt(GenDualNum!1.nan).same(GenDualNum!1.nan));
}

/+
pure DerivSeqType sin(DerivSeqType)(const DerivSeqType u)
{
    static if (DerivSeqType.Order == 0) {
        return DerivSeqType(std.math.sin(u.val()));
    } else {
        return DerivSeqType(
            std.math.sin(u.val()), 
            u.d() * cos(u.reduce()));
    }
}
+/

/+
// TODO learn why remainder is not pure
DerivSeqType tan(DerivSeqType)(const DerivSeqType u)
in 
{
    assert(remainder(u.val() - PI_2, PI) != 0);
}
do 
{
    static if (DerivSeqType.Order == 0) {
        return DerivSeqType(std.math.tan(u.val()));
    } else {
        return sin(u) / cos(u);
    }
}

pure DerivSeqType acos(DerivSeqType)(const DerivSeqType u)
in 
{
    static if (DerivSeqType.Order == 0) {
        assert(-1 <= u.val() && u.val() <= 1);
    } else {
        assert(-1 < u.val() && u.val() < 1);
    }
}
do 
{
    static if (DerivSeqType.Order == 0) {
        return std.math.acos(u.val());
    } else {
        return DerivSeq( 
            std.math.acos(u.val()), 
            -u.d() / sqrt(1 - pow(u.reduce(), 2)));
    }
}
+/