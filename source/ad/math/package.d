/**
* This package extends the `std.math` Phobos package to support `GDN` objects.
*
* It is decomposed into modules in the same way that `std.math` is. It also exports all of the
* symbols from these modules just like `std.math` does.
*/
module ad.math;

public import ad.math.algebraic;
public import ad.math.constants;
public import ad.math.exponential;
public import ad.math.operations;
public import ad.math.remainder;
public import ad.math.rounding;
public import ad.math.traits;
public import ad.math.trigonometry;