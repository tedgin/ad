/** TODO: finish documentation after all of the submodules are completed.
 * This module extends the `core.math` and `std.math` libraries to support `GDN` objects. It is
 * decomposed into submodules in the same way that std.math is. It also exports all of the symbols
 * from `ad.math`, `core.math` and `std.math` to make it easier to work with real and generalized
 * dual numbers together.
 */
module ad.math;

public import ad.math.algebraic;
public import ad.math.constants;
public import ad.math.core;
public import ad.math.exponential;
public import ad.math.operations;
public import ad.math.remainder;
public import ad.math.rounding;
public import ad.math.special;
public import ad.math.traits;
public import ad.math.trigonometry;