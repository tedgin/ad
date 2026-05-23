# ad

This is an automatic differentiation library written in D supporting one-dimensional, real-valued derivatives of arbitrary order. It is not a high performance library. It could become one; I'm not against that by any means! It's just that I originally built this library after learning about automatic differentiation. The elegance of the concept struct me. I had to implement it.

## Features

* supports all of D's arithmetic operators
* supports the same set of functions as `core.math`
* supports the same set of functions as `std.math`
* supports the same set of functions as `std.mathspecial`
* supports arbitrary order differentiation, must be fixed at compile time

## Overview

This library consists of a handful of modules. [`ad`](source/ad/package.d) provides the generalized dual number type `GDN` and its basic arithmetic operations that are available to the floating point type `real`. The remaining modules provide implementations of all of the mathematical operations defined in the phobos modules `core.math`, `std.math`, and `std.mathspecial` for real numbers to generalized dual numbers. The are organized similarly to the modules in the `std.math` package. Here is the mapping.

* [`ad.math.algebraic`](source/ad/math/algebraic.d) → `std.math.algebraic`
* [`ad.math.constants`](source/ad/math/constants.d) → `std.math.constants`
* [`ad.math.exponential`](source/ad/math/exponential.d) → `std.math.exponential`
* [`ad.math.operations`](source/ad/math/operations.d) → `std.math.operations`
* [`ad.math.remainder`](source/ad/math/remainder.d) → `std.math.remainder`
* [`ad.math.rounding`](source/ad/math/rounding.d) → `std.math.rounding`
* [`ad.math.special`](source/ad/math/special.d) → `std.mathspecial`
* [`ad.math.traits`](source/ad/math/traits.d) → `std.math.traits`
* [`ad.math.trigonometry`](source/ad/math/trigonometry.d) → `std.math.trigonometry`

The module [`ad.math`](source/ad/math/package.d) provides implementations of all of the operations defined in `core.math` for real numbers to the `GDN` type. It also aggregates and exposes all of the functions defined in its submodules in the same way the Phobos package `std.math` does.

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

This library is built with DUB.

To build the library, run:

```bash
dub build
```

To build the API documentation from the DDoc configuration, run:

```bash
dub build --config=docs
```

This command will generate `docs/` folder with HTML reference pages.

## Using the library in an application

Add `ad` github repository as a dependency in your application's `dub.json` or `dub.sdl`.

Example `dub.json` dependency:

```json
"dependencies": {
   "ad": {
      "repository": "git+https://github.com/tedgin/ad.git"
   }
}
```

## Future Work

The implementations of the regularized lower incomplete gamma function $P(s,x)$ and the regularized incomplete beta function $I_x(a,b)$ don't support differention of their parameters. I would like to extend the library to support this.

For the regularized lower incomplete gamma function $P(s,x)$, implementating $\frac{∂P}{∂s}$ looks possible. $\frac{∂P}{∂s} = -\frac{∂Q}{∂s}$, where $Q(s,x)$ is the regularized upper incomplete gamma function. $\frac{∂Q}{∂s} = \lbrack\frac{∂𝛤}{∂s}(s,x)𝛤(s) - 𝛤(s,x)\frac{d𝛤}{ds}(s)\rbrack/𝛤(s)^2$ where $𝛤(s,x)$ is the upper incomplete gamma function, and $𝛤(s)$ is the gamma function. $\frac{∂𝛤}{∂s}(s,x) = ln(x)𝛤(s,x) + xT(3,s,x)$, where $T(m,s,x) = G^{m,0}_{m-1,m}(0,0,…,0; s-1,-1,…,-1 |x)$, and $G$ is the Meijer G-function. $T$ has recurrent derivative formulas for both $s$ and $x$: $\frac{∂T}{∂s}(m,s,x) = ln(x)T(m,s,x) + (m-1)T(m+1,s,x)$, and $\frac{∂T}{∂x}(m,s,x) = -\frac{T(m-1,s,x) + T(m,s,x)}{x}$.

For the regularized incomplete beta function $I_x(a,b)$, I haven't found any recurrent derivative formulas for computing $\frac{∂^nI_x}{∂a^n}$ and $\frac{∂^nI_x}{∂b^n}$, so it may not be possible to implement these.
