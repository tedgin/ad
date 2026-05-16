# ad

This is an automatic differentiation library written in D supporting one-dimensional, real-valued derivatives of arbitrary order. It is not a high performance library. It could become one; I'm not against that by any means! It's just that I originally built this library after learning about automatic differentiation. The elegance of the concept struct me. I had to implement it.

TODO: document modules and package

## Features

* supports all of D's arithmetic operators
* supports the same set of functions as `core.math`
* supports the same set of functions as `std.math`
* supports the same set of functions as `std.mathspecial`
* supports arbitrary order differentiation, must be fixed at compile time

## Overview

<!-- TODO: write this secion -->
To be continued ...

## Examples

Here are some examples of using this library.

### Basic Usage

To differentiate a function, create a `GDN` variable with the desired derivative degree and evaluate your function:

```d
import ad;

void main()
do {
   // Create a variable x = 3 that can track up to first order
   // derivatives
   const x = GDN!1(3);

   // Evaluate a function: f(x) = 2x + 1
   const f = 2 * x + 1;

   // The value and first derivative
   assert(f == 7);    // f(3)
   assert(f.d == 2);  // f'(3), f'(x) = 2
}
```

### Computing Higher-Order Derivatives

Specify a higher degree to compute multiple orders of derivatives:

```d
import ad;

void main()
do {
   // Create a variable x = 2 that tracks up to third order
   // derivatives
   const x = GDN!3(2);

   // Evaluate a function: f(x) = x³
   const f = x ^^ 3;

   // Access derivatives of different orders
   assert(f == 8);       // f(2)
   assert(f.d == 12);    // f'(2), f'(x) = 3x²
   assert(f.d!2 == 12);  // f''(2), f''(x) = 6x
   assert(f.d!3 == 6);   // f'''(2), f'''(x) = 6
}
```

### Using with Math Functions

The library supports standard math functions from `core.math`, `std.math`, and `std.mathspecial`:

```d
import std.math : E;
import std.mathspecial : digamma;

import ad;
import ad.math;

void main()
do {
   // Differentiate trigonometric functions
   const x = GDN!1(0);
   const y = sin(x);
   assert(y == 0);    // sin(0) = 0
   assert(y.d == 1);  // sin'(x) = cos(x)

   // Differentiate exponential functions
   const u = exp(GDN!1(1));
   assert(u == E);
   assert(u.d == E);  // (d/dx)eˣ = eˣ

   // Differentiate the Gamma function
   const q = gamma(GDN!1(2));
   assert(q == 1);             // Γ(2) = 1
   assert(q.d == digamma(2));  // Γ'(x) = Γ(x)Ψ(x)
}
```

### Combining Operations

GDN objects can be freely mixed with arithmetic operations and math functions:

```d
import std.math : E, cos;

import ad;
import ad.math;

void main()
do {
   // x is the result of evaluating a function whose value is 1
   // and derivative is 2.
   const x = GDN!1(1, 2);

   // Evaluate f(x) = eˣsin(x)
   const f = exp(x) * sin(x);
   assert(f == E * sin(1.0L));

   // f'(x) = eˣx'sin(x) + eˣcos(x)x'
   //       = x'eˣ[cos(x) + sin(x)]
   assert(f.d == 2*E*(cos(1.0L) + sin(1.0L)));
}
```

## Building

<!-- TODO: write this section -->
To be continued ...

## Future work

<!-- TODO: write this section -->
gammaIncomplete
betaIncomplete
