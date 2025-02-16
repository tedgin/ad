# ad

This is an automatic differentiation library written in D supporting one-dimensional, real-valued derivatives of arbitrary order. It is not a high performance library. It could become one; I'm not against that by any means! It's just that I originally built this library after learning about automatic differentiation. The elegance of the concept struct me. I had to implement it.

## Features

* supports all of D's arithmetic operators
* support the same set of functions as `core.math`
* supports the same set of functions as `std.math`
  * `std.math.constants`
  * `std.math.algebraic`
    * `abs` (implicitly)
    * `fabs`
    * `sqrt`
    * `cbrt`
    * `hypot`
    * `poly`
    * `nextPow2` (not implemented yet)
    * `truncPow2` (not implemented yet)
  * `std.math.trigonometry` (not implemented yet)
  * `std.math.rounding`
    * `ceil`
    * `floor` (not implemented yet)
    * `lrint` (not implemented yet)
    * `lround` (not implemented yet)
    * `nearbyint` (not implemented yet)
    * `quantize` (not implemented yet)
    * `rint` (not implemented yet)
    * `rndtol` (not implemented yet)
    * `round` (not implemented yet)
    * `trunc` (not implemented yet)
  * `std.math.exponential` (not implemented yet)
  * `std.math.remainder` (not implemented yet)
  * `std.math.operations` (not implemented yet)
  * `std.math.traits`
    * `copysign` (not implemented yet)
    * `isFinite` (not implemented yet)
    * `isIdentical` (not implemented yet)
    * `isInfinity` (not implemented yet)
    * `isNaN` (not implemented yet)
    * `isNormal` (not implemented yet)
    * `isPowerOf2` (not implemented yet)
    * `isSubnormal` (not implemented yet)
    * `sgn`
    * `signbit` (not implemented yet)
  * `std.math.hardware` (not implemented yet)
* supports arbitrary order differentiation, must be fixed at compile time

## Examples

To be continued ...

## Overview

To be continued ...

## Building

To be continued ...
