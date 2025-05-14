# ad

This is an automatic differentiation library written in D supporting one-dimensional, real-valued derivatives of arbitrary order. It is not a high performance library. It could become one; I'm not against that by any means! It's just that I originally built this library after learning about automatic differentiation. The elegance of the concept struct me. I had to implement it.

## Features

* supports all of D's arithmetic operators
* supports the same set of functions as `core.math`
* supports the same set of functions as `std.math`
  * `std.math.algebraic`
  * `std.math.constants`
  * `std.math.exponential`
  * `std.math.operations`
  * `std.math.remainder` (not implemented yet)
  * `std.math.rounding`
    * `ceil`
    * `floor`
    * `lrint`
    * `lround`
    * `nearbyint`
    * `quantize`
    * `rint`
    * `rndtol`
    * `round` (not implemented yet)
    * `trunc` (not implemented yet)
  * `std.math.trigonometry`
    * `acos`
    * `acosh` (not implemented yet)
    * `asin`
    * `asinh` (not implemented yet)
    * `atan`
    * `atan2` (not implemented yet)
    * `atanh` (not implemented yet)
    * `cos`
    * `cosh` (not implemented yet)
    * `sin`
    * `sinh` (not implemented yet)
    * `tan`
    * `tanh` (not implemented yet)
  * `std.math.traits`
* supports the same set of functions as `std.mathspecial`
* supports arbitrary order differentiation, must be fixed at compile time

## Examples

To be continued ...

## Overview

To be continued ...

## Building

To be continued ...
