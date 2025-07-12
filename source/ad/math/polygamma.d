module ad.math.polygamma;

import core.math: cos, fabs, rndtol, sin;
import std.math: abs, ceil, isInfinity, isNaN, log10, nextUp, PI, poly, sgn, signbit, trunc;
import std.traits: isIntegral;

private pure nothrow @nogc @safe bool isOdd(T)(in T n) if (isIntegral!T)
{
	return (n & 1) == 1;
}

unittest
{
	assert(!isOdd(0));
	assert(isOdd(1));
	assert(!isOdd(-2));
}

// n!
private pure nothrow @nogc @safe real factorial(in ulong n)
{
	real res = 1;

	for (ulong i = 2; i <= n; i++) {
		res *= i;
	}

	return res;
}

unittest
{
	assert(factorial(0) == 1);
	assert(factorial(1) == 1);
	assert(factorial(2) == 2);
	assert(factorial(10) == 3_628_800);
}


// Implemented using Seidel's algorithm
private nothrow pure @safe real[] initEulerZigZag()
{
	real[][] s;
	s.length = 1;
	s[0].length = 1;
	s[0][0] = 1;

	for (auto k = 1; true; k++) {
		s.length++;
		s[k].length = k + 1;

		if (isOdd(k)) {
			auto j = 0;
			s[k][j] = s[k-1][j];

			for (j++; j < k; j++) {
				s[k][j] = s[k][j-1] + s[k-1][j];
			}

			s[k][j] = s[k][j-1];
		} else {
			auto j = k;
			s[k][j] = s[k-1][j-1];

			for (j--; j > 0; j--) {
				s[k][j] = s[k][j+1] + s[k-1][j-1];
			}

			s[k][j] = s[k][j+1];
		}

		if (isInfinity(isOdd(k) ? s[k][0] : s[k][$-1])) break;
	}

	real [] t;
	t.length = cast(long)s.length - 1;
	t[0] = 1;
	for (auto k = 1; k < t.length; k++) t[k] = isOdd(k) ? s[k][0] : s[k][$-1];
	return t;
}

private static immutable real[] T = initEulerZigZag();

// Euler Zig Zag Numbers
private pure nothrow @nogc @safe real eulerZigZag(in ulong n)
{

	if (n < T.length) {
		return T[n];
	} else {
		return real.infinity;
	}
}

unittest
{
	assert(eulerZigZag(1) == 1);
	assert(eulerZigZag(1) == 1);
	assert(eulerZigZag(2) == 1);
	assert(eulerZigZag(3) == 2);
	assert(eulerZigZag(5) == 16);
	assert(eulerZigZag(11) == 353_792);
	assert(eulerZigZag(2000) == real.infinity);
}


// Bernoulli numbers
private pure nothrow @nogc @safe real bernoulli(int Sign)(in ulong n)
{
	if (n == 0) return 1;

	if (n == 1) {
		static if (signbit!float(Sign) == 0)
			return 1.0L / 2;
		else
			return -1.0L / 2;
	}

	if (isOdd(n)) return 0;
	return (-1.0L)^^(n / 2) * n * (eulerZigZag(n-1) / (2.0L^^n - 4.0L^^n));
}


unittest
{
	import std.format: format;

	assert(bernoulli!1(0) == 1);
	assert(bernoulli!1(1) == 1.0L / 2);
	assert(bernoulli!1(2) == 1.0L / 6, format("%s", bernoulli!1(2)));

	auto s = bernoulli!1(20);
	assert(s == -174_611.0L/330,  format("%s != %s", s, -174_611.0L/330));

	assert(bernoulli!(-1)(0) == 1);
	assert(bernoulli!(-1)(1) == -1.0L / 2);
	assert(bernoulli!(-1)(2) == 1.0L / 6);
}


// Polynomial coefficients for derivatives of cot(πx)
private pure nothrow @safe real[] polygammaReflectPolyCoef(in ulong n)
{
	if (n == 1) return [-1.0L, 0, 0];

	const prev = polygammaReflectPolyCoef(n - 1);

	real[] coef;
	coef.length = n + 2;  // polynomial has degree n+1
	coef[] = 0;

	// Apply recursion formula: Pₙ₊₁(x) = -{(n+1)xPₙ(x) + (1-x²)P'ₙ(x)}

	// Term from xPₙ(x)
	for (long k = 1; k <= n; k++) {
		coef[k-1] -= k * prev[k];
	}

	// Term from (1-x²)P'ₙ(x)
	for (long k = 0; k <= n; k++) {
		coef[k + 1] += (k - cast(long)n) * prev[k];
	}

	return coef;
}

unittest
{
	import std.format: format;

	// p1(x) = -[x^2 + (1-x^2)] = -1 + 0x + 0x^2
	const p1 = polygammaReflectPolyCoef(1);
	assert(p1 == [-1, 0, 0], format("p_1 = %s", p1));

	// p2(x) = -(0 + 0x) + (-n-1)(-1)x + 0x^2 + 0x^3
	const p2 = polygammaReflectPolyCoef(2);
	assert(p2 == [0, 2, 0, 0], format("p_2 = %s", p2));

	const p9 = polygammaReflectPolyCoef(9);
	assert(
		p9 == [-7_936, 0, -137_216, 0.0L, -185_856, 0, -31_616, 0, -256, 0, 0],
		format("p_9 = %s", p9));
}


private pure nothrow @nogc @safe real polygammaReflectDelta(ulong n)(in real x)
{
	static immutable coef = polygammaReflectPolyCoef(n);

	// Reduce the round-off error when PI*x = +/-PI/2
	const cs = (fabs(x)%1 == 0.5) ? 0.0L : cos(PI*x);

	// Compute πⁿ⁺¹Pₙ(cos(πx))/sinⁿ⁺¹(πx)
	return (PI / sin(PI*x))^^(n + 1) * poly(cs, coef);
}

unittest
{
	import std.format: format;
	import std.math: isClose;

	real act, exp;

	assert(polygammaReflectDelta!1(-0.5) == -PI^^2);

	exp = -(PI / sin(0.1*PI))^^2;
	act = polygammaReflectDelta!1(-0.1);
	assert(
		isClose(act, exp, real.epsilon),
		format("exp=%s, act=%s, err=%s > %s", exp, act, fabs(act - exp), fabs(real.epsilon*exp)));

	assert(polygammaReflectDelta!1(-1.5) == -PI^^2);
	assert(polygammaReflectDelta!1(-99.5) == -PI^^2);

	exp = 0;
	act = polygammaReflectDelta!2(-0.5);
	assert(
		isClose(act, exp, real.epsilon),
		format("exp=%s, act=%s, err=%s > %s", exp, act, fabs(act - exp), fabs(real.epsilon*exp)));

	exp = -2 * cos(1.2*PI) * (PI / sin(1.2*PI))^^3;
	act = polygammaReflectDelta!2(-1.2);
	assert(
		isClose(act, exp, real.epsilon),
		format("exp=%s, act=%s, err=%s > %s", exp, act, fabs(act - exp), fabs(real.epsilon*exp)));

	exp = -7_936 * PI^^10;
	act = polygammaReflectDelta!9(-0.5);
	assert(
		isClose(act, exp, real.epsilon),
		format("exp=%s, act=%s, err=%s > %s", exp, act, fabs(act - exp), fabs(real.epsilon*exp)));
}


private pure nothrow @nogc @safe real polygammaRecureDelta(ulong N)(in real x, in long shift)
{
	static const real scale = (-1.0L)^^N * factorial(N);
	static const pow = -1.0L*(N + 1);

	const x0 = x - (shift<0 ? 1 : 0);

	real delta = 0;
	for(auto i = 0; i < abs(shift); i++) delta += (x0 + i*sgn(shift))^^pow;
	return sgn(shift) * scale * delta;
}

unittest
{
	import std.format: format;

	real y;

	y = polygammaRecureDelta!1(1, 1);
	assert(y == -1, format("%s", y));

	assert(polygammaRecureDelta!1(2, 1) == -0.25);
	assert(polygammaRecureDelta!1(1, 3) == -1.25 - 1.0L/9);
	assert(polygammaRecureDelta!2(1, 1) == 2);
}


// Ψₙ(x) ~ (-1)ⁿ⁻¹{(n-1)!/xⁿ + n!/(2xⁿ⁺¹) + ∑ₖ₌₁∞[(-1)ᵏ|B⁻₂ₖ|/(2k)!](2k+n-1)!/x²ᵏ⁺ⁿ}, for n > 0,
// 2x > n + ⌊log(1/ε)⌋, and x > n.
private pure nothrow @nogc @safe real polygammaAsymptoticSeries(ulong N)(in real x) if (N > 0)
{
	static const real sign = (-1.) ^^ (N - 1);
	static const term_0 = sign * factorial(N - 1);
	static const term_1 = sign * factorial(N);
	static const max_x = (cast(real)long.max + N/2.0L) / PI;

	if (isNaN(x)) return real.nan;
	if (x == real.infinity) return sign * 0;

	if (x > max_x) {
		// (n-1)!/xⁿ + n!/(2xⁿ⁺¹) ≤ (-1)ⁿ⁺¹Ψₙ(x) ≤ (n-1)!/xⁿ + n!/xⁿ⁺¹
		// Ψₙ(x) ≈ (-1)ⁿ⁺¹[(n-1)!/xⁿ + 3n!/(4xⁿ⁺¹)]
		return term_0/x^^N + 0.75*term_1/x^^(N + 1);
	}

	real series = 0;

	for(ulong k = rndtol(PI*x - N/2.0L); k >= 1; k--) {
		real trunc_fac = 1;
		for(auto j = 2*k + 1; j < 2*k + N; j++) trunc_fac *= j;
		series += trunc_fac * (bernoulli!(-1)(2*k) / x^^(2*k + N));
	}

	return term_0 / x^^N + term_1 / (2 * x^^(N + 1)) + sign * series;
}

unittest
{
	import std.format: format;
	import std.math: isClose;

	real exp, act;
	real acc;

	// n = 1, x = 11
	acc = 0;
	for (auto k = 10; k > 0; k--) acc += 1.0L / k^^2;
	exp = PI^^2 / 6 - acc;
	act = polygammaAsymptoticSeries!1(11);
	assert(
		isClose(act, exp, 10*real.epsilon),
		format("exp=%s, act=%s, err=%s > %s", exp, act, fabs(act - exp), fabs(10*real.epsilon*exp)));

	// n = 1, x = 100
	acc = 0;
	for (auto k = 99; k > 0; k--) acc += 1.0L / k^^2;
	exp = PI^^2 / 6 - acc;
	act = polygammaAsymptoticSeries!1(100);
	assert(
		isClose(act, exp/*, 10*real.epsilon*/),
		format("exp=%s, act=%s, err=%s > %s", exp, act, fabs(act - exp), fabs(10*real.epsilon*exp)));
}


// polygamma
package pure nothrow @nogc @safe real polygamma(ulong N)(in real x) if (N > 0)
{
	static const odd_order = N%2 == 1;
	static const x_cut = (N + ceil(log10(1/real.epsilon))) / 2.0L;

	if (isNaN(x)) return real.nan;

	// If x is large enough, use the asymptotic expansion
	if (x > x_cut) {
		return polygammaAsymptoticSeries!N(x);
	}

	// x is too small. If x is positive, use the recurrence relation.
	if (x > 0) {
		const long x_shift = cast(long) ceil(nextUp(x_cut - x));
		return polygamma!N(x+x_shift) - polygammaRecureDelta!N(x, x_shift);
	}

	// x is not positive. If x is an integer, polygamma has a pole there. Handle
	// the pole. Also handle the case where x is -∞.
	if (x <= 0 && x == trunc(x)) {
		// If x is -0, Ψₙ approaches positive infinity.
		// If x is +0 and n is odd, Ψₙ approaches positive infinity.
		// If x is +0 and n is even Ψₙ approaches negative infinity.
		if (x == 0) {
			return (signbit(x) == 1 || odd_order) ? real.infinity : -real.infinity;
		}

		// x is a negative integer or -∞. If n is odd, Ψₙ approaches infinity,
		// otherwise it diverges
		return odd_order ? real.infinity : real.nan;
	}

	// x is nonintegral negative number. Use the reflection relation.
	return (odd_order ? -1 : 1)*polygamma!N(1-x) - polygammaReflectDelta!N(x);
}

unittest
{
	import std.format: format;
	import std.math: isClose, isNaN;

	void assert_close(ulong n)(in real x, in real expected_result) {
		const act = polygamma!n(x);

		const info = format("exp=%s, act=%s, err=%s > %s",
			expected_result, act, fabs(act - expected_result), fabs(real.epsilon*expected_result));

		assert(isClose(act, expected_result, 10*real.epsilon), info);
	}

	// test n = 1
	assert(isNaN(polygamma!1(real.nan)));
	assert(polygamma!1(-real.infinity) == real.infinity);
	assert(polygamma!1(+0.) == real.infinity);
	assert(polygamma!1(-0.) == real.infinity);
	assert(polygamma!1(-1) == real.infinity);
	assert_close!1(-12.75, 19.663_772_856_722_737_612_034_697_464_751_605L);
	assert_close!1(-0.25, 18.541_879_647_671_606_498_397_662_880_417_078L);
	assert_close!1(9.765_625e-4, 1.048_577_642_589_392_152_617_040_806_167_829_8e6L);
	assert_close!1(0.125, 65.388_133_444_988_034_473_142_999_334_395_961L);
	assert_close!1(10.75, 0.097_483_848_201_852_104_395_946_001_854_344_927L);
	assert_close!1(100, 0.010_050_166_663_333_571_395_245_668_465_701_423L);

	// test n = 2
	assert(isNaN(polygamma!2(-real.infinity)));
	assert(polygamma!2(+0.) == -real.infinity);
	assert(polygamma!2(-0.) == real.infinity);
	assert(isNaN(polygamma!2(-1)));
	assert_close!2(-12.75, -124.030_794_614_158_233_846_041_532_515_436_81L);
	assert_close!2(-0.25, 122.697_366_783_662_360_368_567_293_089_561_59L);
	assert_close!2(9.765_625e-4, -2.147_483_650_397_783_916_396_063_006_390_936_4e9L);
	assert_close!2(0.125, -1.025_753_338_118_135_682_595_668_930_056_517_4e3L);
	assert_close!2(10.75, -9.495_619_644_926_590_077_648_896_563_179_177_5e-3L);
	assert_close!2(100, -1.010_049_998_333_499_970_008_330_044_605_938_2e-4L);

	// test n >= 10
	assert_close!20(-9.5, -1.030_763_776_245_141_659_801_808_179_634_593_2e-3L);
	assert_close!10(-0.25, 1.522_020_442_846_279_111_956_134_706_549_276e13L);
	assert_close!15(9.765_625e-4, 1.911_168_229_927_653_680_431_359_258_893_451_1e60L);
	assert_close!12(0.125, -2.633_339_144_617_578_462_370_751_412_184_393_7e20L);
	assert_close!12(17.125, -8.747_949_394_125_802_795_334_672_654_399_567_7e-8L);
	assert_close!12(100, -4.236_368_168_960_810_441_389_986_390_777_533_3e-17L);
}
