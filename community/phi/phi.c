/* Make SNFS polynomials for cyclotomic numbers.

  Copyright 2005 Alexander Kruppa.
            2022 Jon Becker

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

/* Version: 1.0.0

History:	0.1 initial release
			0.1.1 Added GGNFS format, estimate of cost and time,
				changed x to mpz_t type.
			0.1.2 Makes better polynomials for numbers with composite
				base and no useful algebraic factors
			0.1.2.1 Allow reading N from stdin
			0.1.3 Allow polynomial for homogenous cyclotomic numbers by Polly T. Wang
				changed default output type to GGNFS
			0.1.4 Divide out the algebraic factor when making 6th-degree polynomials
				with exponents that are divisible by 3.
				Fix some bugs in the code from 0.1.2 dealing with composite bases.
			0.1.5 Add support for octic polynomials.  Fix handling of arguments
				for non-homogeneous forms.
			1.0.0 Add support for aurifeuillian factorizations.  These are
				auto-detected. Change default degree.

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <assert.h>

#define VERSION "1.0"

#define uint unsigned int

#define MAXDEGREE 8

/* Output formats */
#define FRANKE 0
#define CWI 1
#define GGNFS 2

#define MAXAURIFDEGREE 8

static int gAurifCoeffs2_3[] = { 1, 1, 1 };
static int gAurifCoeffs5_6[] = { 1, 3, 1, 1, 1 };
static int gAurifCoeffs7[] = { 1, 3, 3, 1, 1, 1, 1 };
static int gAurifCoeffs10[] = { 1, 5, 7, 5, 1, 1, 2, 2, 1 };
static int gAurifCoeffs11[] = { 1, 5, -1, -1, 5, 1, 1, 1, -1, 1, 1 };
static int gAurifCoeffs13[] = { 1, 7, 15, 19, 15, 7, 1, 1, 3, 5, 5, 3, 1 };
static int gAurifCoeffs14[] = { 1, 7, 3, -7, 3, 7, 1, 1, 2, -1, -1, 2,
				1 };
static int gAurifCoeffs15[] = { 1, 8, 13, 8, 1, 1, 3, 3, 1 };
static int gAurifCoeffs17[] = { 1, 9, 11, -5, -15, -5, 11, 9, 1, 1, 3, 1,
				-3, -3, 1, 3, 1 };
static int gAurifCoeffs21[] = { 1, 10, 13, 7, 13, 10, 1, 1, 3, 2, 2, 3,
				1 };
static int gAurifCoeffs30[] = { 1, 15, 38, 45, 43, 45, 38, 15, 1, 1, 5, 8,
				8, 8, 8, 5, 1 };


static int gcd(int a, int b)
{
	return b ? gcd(b, a % b) : a;
}


static int ipwr(int base, int exp)
{
	int i;
	int n = 1;
	
	for (i = 0; i < exp; i++) {
		n *= base;
	}
	return n;
}


/* Compute p[0]^v[0] * ... * p[k-1]^v[k-1] */

static void mpz_mul_list(mpz_t e, int *p, int *v, uint k)
{
	int i;
	uint j;
	
	mpz_set_ui(e, 1);
	for (j = 0; j < k; j++)
		for (i = 0; i < v[j]; i++)
			mpz_mul_ui(e, e, p[j]);
}

/* Find prime factors < maxp of R. Cofactor is put in R, 
 the i-th prime factor in p_i with exponent v_i.
 Returns number of distinct prime factors found. */

static uint trialdiv(mpz_t R, int *p, int *v, uint maxp, uint flen)
{
	uint i, k = 0;
	
#ifdef DEBUG
	printf("trialdiv: ");
#endif
	
	for (i = 2; i < maxp && k < flen; i = (i + 1) | 1)	/* 2 3 5 7 9... */
		if (mpz_divisible_ui_p(R, i)) {
			p[k] = i;
			v[k] = 0;
			do {
				v[k]++;
				mpz_divexact_ui(R, R, i);
			} while (mpz_divisible_ui_p(R, i));
#ifdef DEBUG
			printf("%d^%d ", i, v[k]);
#endif
			k++;
		}
#ifdef DEBUG
	printf("\n");
#endif
	
	return k;
}

static void usage()
{
	printf("Phi " VERSION "\n");
	printf("Usage: phi [options] [n x y] N\n");
	printf
	("Makes SNFS polynomial for factoring the primitive part Phi_n(x,y) of x^n-y^n\n");
	printf
	("Cofactor N must divide x^n-y^n if n is odd, or x^(n/2)+y^(n/2) if n is even\n");
	printf
	("If you do not specify x, y and n, the program will try to auto-detect them.\n");
	printf
	("If you omit y, it will default to a value of 1.\n");
	printf("Options:\n");
	printf("-deg4    Make a degree 4 polynomial\n");
	printf("-deg5    Make a degree 5 polynomial\n");
	printf("-deg6    Make a degree 6 polynomial\n");
	printf("-deg8	 Make a degree 8 polynomial\n");
	printf
	("         By default chooses the degree that minimizes the SNFS difficulty\n");
	printf("-franke  Use Franke file format for polynomial\n");
	printf("-ggnfs   Use GGNFS file format for polynomial\n");
	printf("-cwi     Use CWI file format for polynomial\n");
	printf("-short   Omit unneeded lines for Franke file format\n");
	printf("-        Read N from stdin\n");
}

static int autodet(int *x, int *y, uint * n, mpz_t N)
{
	const int flen = 32;
	const int maxp = 2000;	/* Limit for prime divisors */
	const int maxbase = 1000;	/* Limit for bases to try */
	int p[flen], v[flen];	/* p are primes, v exponents */
	int a, b, i, j, k;
	mpz_t R;			/* Remaining cofactor */
	mpz_t e;			/* exponent */
	
	mpz_init(R);
	mpz_sub_ui(R, N, 1);
	mpz_init(e);
	
	k = trialdiv(R, p, v, maxp, flen);
	
	mpz_mul_list(e, p, v, k);
#ifdef DEBUG
	gmp_printf("\ne = %Zd\n", e);
#endif
	
	
	/*
	 * Note that if we switch these two nested loops and put
	 * b on the outside, then we could just calculate b-inverse once
	 * each trip through the outer loop, instead of having to do
	 * it every time through the inner loop.  But then the inner
	 * loop would have to go up to maxbase every time through the
	 * outer loop, so we would find the numbers with small bases
	 * much later in the process.
	 */
	for (a = 2; a <= maxbase; a++)
		for (b = 1; b <= a - 1; b++) {
			if (gcd(a, b) > 1)
				continue;
			
#ifdef DEBUG
			printf("Trying base %d,%d\n", a, b);
#endif
			mpz_set_ui(R, b);
			mpz_invert(R, R, N);
			mpz_mul_ui(R, R, a);
			mpz_powm(R, R, e, N);
			if (mpz_cmp_ui(R, 1) == 0) {
#ifdef DEBUG
				printf("Discovered! base = %d,%d\n", a, b);
#endif
				goto foundbase;
			}
		}
foundbase:
	if (b > maxbase - 1) {
		/* No suitable base found */
#ifdef DEBUG
		printf("No suitable base found\n");
#endif
		mpz_clear(R);
		mpz_clear(e);
		return 0;
	}
	
	/* Find exponent */
	for (j = 0; j < k; j++) {
		while (v[j] > 0) {
			/* See if lowering this exponent by one still satisfies ord|e */
			v[j]--;
			mpz_mul_list(e, p, v, k);
#ifdef DEBUG
			gmp_printf("Trying if ord | e = %Zd\n", e);
#endif
			mpz_set_ui(R, b);
			mpz_invert(R, R, N);
			mpz_mul_ui(R, R, a);
			mpz_powm(R, R, e, N);
			if (mpz_cmp_ui(R, 1) != 0) {
				/* No, cannot decrease this exponent any further. */
				v[j]++;		/* Undo decrease */
#ifdef DEBUG
				printf("No, leaving exponent of %d at %d\n", p[j], v[j]);
#endif
				break;		/* Try next prime */
			}
#ifdef DEBUG
			else
				printf("Yes, decreasing exponent of %d to %d\n", p[j],
					   v[j]);
#endif
			
		}
	}
	
	mpz_mul_list(e, p, v, k);
	if (mpz_fits_uint_p(e)) {
		if (x != NULL)
			*x = a;
		if (y != NULL)
			*y = b;
		if (n != NULL)
			*n = mpz_get_ui(e);
		i = 1;
	} else
		i = 0;
	
	mpz_clear(R);
	mpz_clear(e);
	
	return i;
}

/* Difficulty must be in base 10 */
static double snfscost(const double difficulty)
{
	const double c = 1.5262856567;	/* (32/9)^(1/3) */
	const double logd = difficulty * 2.3025850930;	/* *log(10) */
	
	return exp(c * pow(logd, 1. / 3) * pow(log(logd), 2. / 3));
}



static void generate_std(mpz_t N, mpz_t x, mpz_t y, int n, int sign, int *degree,
				  mpz_t f[], mpz_t g[], double *difficulty, double *skewness)
{
	int k;
	int halved = 0;
	int i, jp, jq;
	double dx;			/* x as a double */
	double dy;			/* y as a double */
	mpz_t t;			/* Numerator of M */
	mpz_t u;			/* Denominator of M */
	mpz_t r;
	int d = *degree;

	dx = mpz_get_d(x);
	dy = mpz_get_d(y);
	mpz_init(t);
	mpz_init(u);
	mpz_init(r);

	if ((d == 0 || d == 4) && n % 15 == 0) {
		/*
		 * x^(15*k)+-1, remaining size 8/15 = 0.5333...
		 * We have to halve to make it a quartic
		 */
		d = 4;
		k = n / 15;
		mpz_set_si(f[4], 1);
		mpz_set_si(f[3], sign);
		mpz_set_si(f[2], -4);
		mpz_set_si(f[1], -sign * 4);
		mpz_set_si(f[0], 1);
		mpz_pow_ui(t, x, k);
		mpz_pow_ui(u, y, k);
		halved = 1;
		*skewness = 1.;
		*difficulty = log10(dx) * 8. * k;
	} else if ((d == 8) && n % 15 == 0) {
		/* A natural octic already*/
		k = n / 15;
		mpz_set_si(f[8], 1);
		mpz_set_si(f[7], sign);
		mpz_set_si(f[5], -sign);
		mpz_set_si(f[4], -1);
		mpz_set_si(f[3], -sign);
		mpz_set_si(f[1], sign);
		mpz_set_si(f[0], 1);
		mpz_pow_ui(t, x, k);
		mpz_pow_ui(u, y, k);
		*skewness = 1.;
		*difficulty = log10(dx) * 8. * k;
	} else if ((d == 0 || d == 6) && n % 21 == 0) {
		/* x^(21*k)+-1, remaining size 12/21 = 0.571428... */
		d = 6;
		k = n / 21;
		mpz_set_si(f[6], 1);
		mpz_set_si(f[5], sign);
		mpz_set_si(f[4], -6);
		mpz_set_si(f[3], -sign * 6);
		mpz_set_si(f[2], 8);
		mpz_set_si(f[1], sign * 8);
		mpz_set_si(f[0], 1);
		mpz_pow_ui(t, x, k);
		mpz_pow_ui(u, y, k);
		halved = 1;
		*skewness = 1.;
		*difficulty = log10(dx) * 12. * k;
	} else if (d == 4 && n % 3 == 0) {
		/* Exponent divisible by 3, making a quartic */
		if (n % 2 == 0) {
			/* x^(6*k)+1, remaining size 2/3 = 0.666... */
			k = n / 6;
			mpz_set_si(f[4], 1);
			mpz_set_si(f[2], -1);
			mpz_set_si(f[0], 1);
			mpz_pow_ui(t, x, k);
			mpz_pow_ui(u, y, k);
			*skewness = 1.;
			*difficulty = log10(dx) * 4. * k;
		} else {
			/* x^(6*k+3)+-1, remaining size 2/3 = 0.666... */
			k = (n - 3) / 6;
			mpz_mul(f[4], x, x);
			mpz_mul(f[2], x, y);
			mpz_mul_si(f[2], f[2], -sign);
			mpz_mul(f[0], y, y);
			mpz_pow_ui(t, x, k);
			mpz_pow_ui(u, y, k);
			*skewness = pow(dx / dy, -0.5);
			*difficulty = log10(dx) * (4. * k + 2);
		}
	} else if ((d == 0 || d == 6) && n % 3 == 0) {
		/* Exponent divisible by 3, making a sextic */
		d = 6;
		if (n % 9 == 0) {
			/* x^(9*k)+-1, remaining size 2/3 = 0.666... */
			k = n / 9;
			mpz_set_si(f[6], 1);
			mpz_set_si(f[3], -sign);
			mpz_set_si(f[0], 1);
			mpz_pow_ui(t, x, k);
			mpz_pow_ui(u, y, k);
			*skewness = 1.;
			*difficulty = log10(dx) * 6. * k;
		} else if (n % 9 == 3) {
			/* x^(9*k+3)+-1, remaining size 2/3 = 0.666... */
			k = (n - 3) / 9;
			mpz_mul(f[6], x, x);
			mpz_mul(f[3], x, y);
			mpz_mul_si(f[3], f[3], -sign);
			mpz_mul(f[0], y, y);
			mpz_pow_ui(t, x, k);
			mpz_pow_ui(u, y, k);
			*skewness = pow(dx / dy, -1./3);
			*difficulty = log10(dx) * (6. * k + 2);
		} else {
			/* x^(9*k-3)+-1, remaining size 2/3 = 0.666... */
			k = (n + 3) / 9;
			mpz_mul(f[6], y, y);
			mpz_mul(f[3], x, y);
			mpz_mul_si(f[3], f[3], -sign);
			mpz_mul(f[0], x, x);
			mpz_pow_ui(t, x, k);
			mpz_pow_ui(u, y, k);
			*skewness = pow(dx / dy, 1./3);
			*difficulty = log10(dy) * 2. + log10(dx) * 6. * k;
		}
	} else if (d == 8 && n % 3 == 0) {
		/* Exponent divisible by 3, making an octic */
		k = n / 3;
		if (k % 4 == 0) {
			k = k / 4;
			mpz_set_si(f[8], 1);
			mpz_set_si(f[4], -1);
			mpz_set_si(f[0], 1);
			mpz_pow_ui(t, x, k);
			mpz_pow_ui(u, y, k);
			*skewness = 1.;
			*difficulty = log10(dx) * 8. * k;
		} else if (k % 4 == 1) {
			k = (k - 1) / 4;
			mpz_mul(f[8], x, x);
			mpz_mul(f[4], x, y);
			mpz_mul_si(f[4], f[4], -sign);
			mpz_mul(f[0], y, y);
			mpz_pow_ui(t, x, k);
			mpz_pow_ui(u, y, k);
			*skewness = pow(dx / dy, -0.25);
			*difficulty = log10(dx) * (8. * k + 2);
		} else if (k % 4 == 2) {
			k = (k - 2) / 4;
			mpz_pow_ui(f[8], x, 4);
			mpz_mul(f[4], x, x);
			mpz_mul(f[4], f[4], y);
			mpz_mul(f[4], f[4], y);
			mpz_mul_si(f[4], f[4], -1);
			mpz_pow_ui(f[0], y, 4);
			mpz_pow_ui(t, x, k);
			mpz_pow_ui(u, y, k);
			*skewness = pow(dx / dy, -0.5);
			*difficulty = log10(dx) * (8. * k + 4);
		} else {
			k = (k + 1) / 4;
			mpz_mul(f[8], y, y);
			mpz_mul(f[4], x, y);
			mpz_mul_si(f[4], f[4], -sign);
			mpz_mul(f[0], x, x);
			mpz_pow_ui(t, x, k);
			mpz_pow_ui(u, y, k);
			*skewness = pow(dx / dy, 0.25);
			*difficulty = log10(dy) * 2. + log10(dx) * 8. * k;
		}
	} else if ((d == 0 || d == 4) && n % 5 == 0) {
		/* x^(5*k)+-1, remaining size 4/5 = 0.8 */
		d = 4;
		k = n / 5;
		mpz_set_si(f[4], 1);
		mpz_set_si(f[3], -sign);
		mpz_set_si(f[2], 1);
		mpz_set_si(f[1], -sign);
		mpz_set_si(f[0], 1);
		mpz_pow_ui(t, x, k);
		mpz_pow_ui(u, y, k);
		*skewness = 1.;
		*difficulty = log10(dx) * 4. * k;
	} else if (d == 8 && n % 5 == 0) {
		/* Exponent divisible by 5, making an octic */
		if (n % 2 == 0) {
			k = n / 10;
			mpz_set_si(f[8], 1);
			mpz_set_si(f[6], -1);
			mpz_set_si(f[4], 1);
			mpz_set_si(f[2], -1);
			mpz_set_si(f[0], 1);
			mpz_pow_ui(t, x, k);
			mpz_pow_ui(u, y, k);
			*skewness = 1.;
			*difficulty = log10(dx) * 8. * k;
		} else {
			k = (n - 5) / 10;
			mpz_pow_ui(f[8], x, 4);
			mpz_pow_ui(f[6], x, 3);
			mpz_mul(f[6], f[6], y);
			mpz_mul_si(f[6], f[6], -sign);
			mpz_mul(f[4], x, x);
			mpz_mul(f[4], f[4], y);
			mpz_mul(f[4], f[4], y);
			mpz_pow_ui(f[2], y, 3);
			mpz_mul(f[2], f[2], x);
			mpz_mul_si(f[2], f[2], -sign);
			mpz_pow_ui(f[0], y, 4);
			mpz_pow_ui(t, x, k);
			mpz_pow_ui(u, y, k);
			*skewness = pow(dx / dy, -0.5);
			*difficulty = log10(dx) * (8. * k + 4);
		}
	} else if ((d == 0 || d == 6) && n % 7 == 0) {
		/* x^(7*k)+-1, remaining size 6/7 = 0.857142... */
		d = 6;
		k = n / 7;
		mpz_set_si(f[6], 1);
		mpz_set_si(f[5], -sign);
		mpz_set_si(f[4], 1);
		mpz_set_si(f[3], -sign);
		mpz_set_si(f[2], 1);
		mpz_set_si(f[1], -sign);
		mpz_set_si(f[0], 1);
		mpz_pow_ui(t, x, k);
		mpz_pow_ui(u, y, k);
		*skewness = 1.;
		*difficulty = log10(dx) * 6. * k;
	} else if ((d == 0 || d == 5) && n % 11 == 0) {
		/* x^(11*k)+-1, remaining size 10/11 = 0.90... */
		d = 5;
		k = n / 11;
		mpz_set_si(f[5], 1);
		mpz_set_si(f[4], -sign * 1);
		mpz_set_si(f[3], -4);
		mpz_set_si(f[2], sign * 3);
		mpz_set_si(f[1], 3);
		mpz_set_si(f[0], -sign);
		mpz_pow_ui(t, x, k);
		mpz_pow_ui(u, y, k);
		halved = 1;
		*skewness = 1.;
		*difficulty = log10(dx) * 10. * k;
	} else if ((d == 0 || d == 6) && n % 13 == 0) {
		/* x^(13*k)+-1, remaining size 12/13 = 0.923076... */
		d = 6;
		k = n / 13;
		mpz_set_si(f[6], 1);
		mpz_set_si(f[5], -sign);
		mpz_set_si(f[4], -5);
		mpz_set_si(f[3], sign * 4);
		mpz_set_si(f[2], 6);
		mpz_set_si(f[1], -sign * 3);
		mpz_set_si(f[0], -1);
		mpz_pow_ui(t, x, k);
		mpz_pow_ui(u, y, k);
		halved = 1;
		*skewness = 1.;
		*difficulty = log10(dx) * 12. * k;
	} else if ((d == 0 || d == 8) && n % 17 == 0) {
		d = 8;
		k = n / 17;
		mpz_set_si(f[8], 1);
		mpz_set_si(f[7], -sign);
		mpz_set_si(f[6], -7);
		mpz_set_si(f[5], sign * 6);
		mpz_set_si(f[4], 15);
		mpz_set_si(f[3], -sign * 10);
		mpz_set_si(f[2], -10);
		mpz_set_si(f[1], sign * 4);
		mpz_set_si(f[0], 1);
		mpz_pow_ui(t, x, k);
		mpz_pow_ui(u, y, k);
		halved = 1;
		*skewness = 1.;
		*difficulty = log10(dx) * 16. * k;
	} else {
		/* Make polynomial without using any algebraic factor */
		const int maxp = 10000;
		const int flen = 15;	/* Only examine 16 different factors
								 (including cofactor). Should easily suffice */
		int p[flen], vp[flen], q[flen], vq[flen], bestp, bestq;
		int l, m;
		mpz_t R, S;
		double rating, bestrating;
		
		if (d == 0) d = 5;
		mpz_init_set(R, x);
		mpz_init_set(S, y);
		l = trialdiv(R, p, vp, maxp, flen);
		m = trialdiv(S, q, vq, maxp, flen);
		/* Now x = p[0] ^ v[0] * ... * p[l - 1] ^ v[l - 1] * R */
		/* Now y = q[0] ^ v[0] * ... * q[m - 1] ^ v[m - 1] * S */
		/* Here, l and m is the number of distinct prime factors */
		
		for (i = 0; i < l; i++)
			vp[i] *= n;
		for (i = 0; i < m; i++)
			vq[i] *= n;
		/* Now N = p[0] ^ v[0] * ... * p[l - 1] ^ v[l - 1] * R^n -
		 q[0] ^ v[0] * ... * q[m - 1] ^ v[m - 1] * S^n */
		
		/* Make all exponents multiples of $d$.
		 When decreasing the exponent of $p$ by $v$, we need to multiply
		 $f_d$ by $p^v$.
		 When increasing the exponent of $p$ by $v$, we need to multiply
		 $f_0$ by $p^v$ and $t$ by $p$.
		 When decreasing the exponent of $q$ by $v$, we need to multiply
		 $f_0$ by $q^v$.
		 When increasing the exponent of $q$ by $v$, we need to multiply
		 $f_d$ by $q^v$ and $u$ by $p$.
		 For each exponent $v$, we can choose between increasing by
		 $(-v) \bmod d$ or decreasing by $v \bmod d$.
		 So we have a search space of $2^{k+l+2}$ possibilities.
		 
		 For now, we simply use $f_d + f_0 + a + b$ as the function to minimise,
		 where a and b are what we multiplied t and u by.
		 This rating function could and should be improved, but it seems to
		 usually produce decent polynomials.
		 */
		
		for (jp = 0; jp < (1 << (l + (mpz_cmp_ui(R, 1) > 0))); jp++)
			for (jq = 0; jq < (1 << (m + (mpz_cmp_ui(S, 1) > 0))); jq++) {
				double fd = 1., f0 = 1., a = 1., b = 1.;
				for (i = 0; i < l; i++) {
					if (jp & (1 << i)) {
						/* Reduce exponent by v[i] % d */
						fd *= pow(p[i], vp[i] % d);
					} else if (vp[i] % d != 0) {
						/* Multiply by -v[i] % d */
						f0 *= pow(p[i], d - vp[i] % d);
						a *= p[i];
					}
				}
				if (jp & (1 << l))
					fd *= pow(mpz_get_d(R), n % d);
				else {
					f0 *= pow(mpz_get_d(R), (d - n % d) % d);
					a *= mpz_get_d(R);
				}
				for (i = 0; i < m; i++) {
					if (jq & (1 << i)) {
						/* Reduce exponent by v[i] % d */
						f0 *= pow(q[i], vq[i] % d);
					} else if (vq[i] % d != 0) {
						/* Multiply by -v[i] % d */
						fd *= pow(q[i], d - vq[i] % d);
						b *= q[i];
					}
				}
				if (jq & (1 << m))
					f0 *= pow(mpz_get_d(S), n % d);
				else {
					fd *= pow(mpz_get_d(S), (d - n % d) % d);
					b *= mpz_get_d(S);
				}
				
				/* rating = c * c * pow (fd, 1. + 1./d) * pow (f0, 1. - 1./d); */
				rating = fd + f0 + a + b;
#ifdef DEBUG
				printf
				("j = %d / %d, %.0f * x^%d %c %.0f, a = %.0f, b = %.0f, rating = %f\n",
				 jp, jq, fd, d, (sign > 0) ? '+' : '-', f0, a, b, rating);
#endif
				if ((jp == 0 && jq == 0) || rating < bestrating) {
					bestp = jp;
					bestq = jq;
					bestrating = rating;
#ifdef DEBUG
					printf("New best rating!\n");
#endif
				}
				
			}
		
		*difficulty = log10(dx) * n;
		mpz_set_ui(t, 1);
		mpz_set_ui(u, 1);
		
		/* Now make the actual coefficients f_d and f_0 */
		
		mpz_set_ui(f[d], 1);
		mpz_set_ui(f[0], 1);
		for (i = 0; i < l; i++) {
			if (bestp & (1 << i)) {
				uint v = vp[i] % d;
				mpz_ui_pow_ui(r, p[i], v);
				mpz_mul(f[d], f[d], r);
				mpz_ui_pow_ui(r, p[i], (vp[i] - v) / d);
				mpz_mul(t, t, r);
			} else {
				uint v = (d - vp[i] % d) % d;	/* v = (-vp[i]) % d */
				mpz_ui_pow_ui(r, p[i], v);
				mpz_mul(f[0], f[0], r);
				*difficulty += log10(p[i]) * v;
				mpz_ui_pow_ui(r, p[i], (vp[i] + v) / d);
				mpz_mul(t, t, r);
			}
		}
		if (bestp & (1 << l)) {
			uint v = n % d;
			mpz_pow_ui(r, R, v);
			mpz_mul(f[d], f[d], r);
			mpz_pow_ui(r, R, (n - v) / d);
			mpz_mul(t, t, r);
		} else {
			uint v = (d - n % d) % d;
			mpz_pow_ui(r, R, v);
			mpz_mul(f[0], f[0], r);
			mpz_pow_ui(r, R, (n + v) / d);
			mpz_mul(t, t, r);
			*difficulty += log10(mpz_get_d(R)) * v;
		}
		for (i = 0; i < m; i++) {
			if (bestq & (1 << i)) {
				uint v = vq[i] % d;
				mpz_ui_pow_ui(r, q[i], v);
				mpz_mul(f[0], f[0], r);
				mpz_ui_pow_ui(r, q[i], (vq[i] - v) / d);
				mpz_mul(u, u, r);
			} else {
				uint v = (d - vq[i] % d) % d;	/* v = (-vp[i]) % d */
				mpz_ui_pow_ui(r, q[i], v);
				mpz_mul(f[d], f[d], r);
				*difficulty += log10(q[i]) * v;
				mpz_ui_pow_ui(r, q[i], (vq[i] + v) / d);
				mpz_mul(u, u, r);
			}
		}
		if (bestq & (1 << m)) {
			uint v = n % d;
			mpz_pow_ui(r, S, v);
			mpz_mul(f[0], f[0], r);
			mpz_pow_ui(r, S, (n - v) / d);
			mpz_mul(u, u, r);
		} else {
			uint v = (d - n % d) % d;
			mpz_pow_ui(r, S, v);
			mpz_mul(f[d], f[d], r);
			mpz_pow_ui(r, S, (n + v) / d);
			mpz_mul(u, u, r);
			*difficulty += log10(mpz_get_d(S)) * v;
		}
		
		*skewness = pow(mpz_get_d(f[0]) / mpz_get_d(f[d]), 1. / d);
		
		if (sign == -1)
			mpz_neg(f[0], f[0]);
	}
	*degree = d;
	if (halved) {
		mpz_mul(g[1], t, u);	/* Y1 = -x^k y^k */
		mpz_neg(g[1], g[1]);
		mpz_mul(g[0], t, t);
		mpz_mul(r, u, u);
		mpz_add(g[0], g[0], r);	/* Y0 = x^(2k)+y^(2k) */
	} else {
		mpz_set(g[0], t);
		mpz_neg(g[1], u);
	}
	mpz_clear(t);
	mpz_clear(u);
	mpz_clear(r);
}


static int find_aurifeuillian(mpz_t N, mpz_t x, mpz_t y, int n, int sbp, int sp)
{
	int h;
	int *highCoeffs;
	int *lowCoeffs;
	int degree;
	mpz_t x_powers[MAXAURIFDEGREE+1];
	mpz_t y_powers[MAXAURIFDEGREE+1];
	int i;
	mpz_t high;
	mpz_t low;
	mpz_t prod;
	mpz_t num;
	int sign;
	
	h = n / sbp;
	switch (sbp) {
		case 2:
		case 3:
			highCoeffs = gAurifCoeffs2_3;
			degree = 1;
			break;
		case 5:
		case 6:
			highCoeffs = gAurifCoeffs5_6;
			degree = 2;
			break;
		case 7:
			highCoeffs = gAurifCoeffs7;
			degree = 3;
			break;
		case 10:
			highCoeffs = gAurifCoeffs10;
			degree = 4;
			break;
		case 11:
			highCoeffs = gAurifCoeffs11;
			degree = 5;
			break;
		case 13:
			highCoeffs = gAurifCoeffs13;
			degree = 6;
			break;
		case 14:
			highCoeffs = gAurifCoeffs14;
			degree = 6;
			break;
		case 15:
			highCoeffs = gAurifCoeffs15;
			degree = 4;
			break;
		case 17:
			highCoeffs = gAurifCoeffs17;
			degree = 8;
			break;
		case 21:
			highCoeffs = gAurifCoeffs21;
			degree = 6;
			break;
		case 30:
			highCoeffs = gAurifCoeffs30;
			degree = 8;
			break;
	}
	lowCoeffs = highCoeffs + (degree + 1);
	mpz_init_set_ui(x_powers[0], 1);
	mpz_init_set_ui(y_powers[0], 1);
	mpz_init(x_powers[1]);
	mpz_pow_ui(x_powers[1], x, h);
	mpz_init(y_powers[1]);
	mpz_pow_ui(y_powers[1], y, h);
	for (i = 2; i <= degree; i++) {
		mpz_init(x_powers[i]);
		mpz_mul(x_powers[i], x_powers[i-1], x_powers[1]);
		mpz_init(y_powers[i]);
		mpz_mul(y_powers[i], y_powers[i-1], y_powers[1]);
	}
	mpz_init(prod);
	
	// Calculate the high half
	mpz_init_set_ui(high, 0);
	for (i = 0; i <= degree; i++) {
		mpz_set(prod, x_powers[degree-i]);
		mpz_mul(prod, prod, y_powers[i]);
		mpz_mul_si(prod, prod, highCoeffs[i]);
		mpz_add(high, high, prod);
	}
	
	// Now do the low half
	mpz_init_set_ui(low, 0);
	for (i = 0; i < degree; i++) {
		mpz_set(prod, x_powers[degree-1-i]);
		mpz_mul(prod, prod, y_powers[i]);
		mpz_mul_si(prod, prod, lowCoeffs[i]);
		mpz_add(low, low, prod);
	}
	mpz_pow_ui(prod, x, (h-1)/2);
	mpz_mul(low, low, prod);
	mpz_pow_ui(prod, y, (h-1)/2);
	mpz_mul(low, low, prod);
	mpz_mul_ui(low, low, sbp * sp);

	mpz_init(num);
	mpz_sub(num, high, low);
	if (mpz_divisible_p(num, N)) {
		sign = -1;
	} else {
		mpz_add(num, high, low);
		if (mpz_divisible_p(num, N)) {
			sign = 1;
		} else {
			sign = 0;
		}
	}
	for (i = 0; i <= degree; i++) {
		mpz_clear(x_powers[i]);
		mpz_clear(y_powers[i]);
	}
	mpz_clear(prod);
	mpz_clear(high);
	mpz_clear(low);
	mpz_clear(num);
	return sign;
}


static int generate_aurif(mpz_t N, mpz_t x, mpz_t y, int n, int sign, int *degree,
					mpz_t f[], mpz_t g[], double *difficulty, double *skewness, int *lm)
{
	int sbp;
	int a, b;
	int xsqr;
	int ysqr;
	int i;
	int d = *degree;
	mpz_t r;
	int primes[6];
	int exps[6];
	int nprimes;
	int success = 0;
	int h;
	int k;
	int LM;  /* -1 for L, 1 for M, 0 for error */
	double dx;			/* x as a double */
	double dy;			/* y as a double */
	enum { MODE_HALF, MODE_NATURAL, MODE_EXPAND };
	int mode;
	int exp_factor = 2;
	int exp_offset = 1;

	dx = mpz_get_d(x);
	dy = mpz_get_d(y);
	mpz_init(r);
	
	/* Find the squarefree base product and extract the square parts */
	mpz_set(r, x);
	nprimes = trialdiv(r, primes, exps, 31, 6);
	if (mpz_cmp_ui(r, 1) != 0) {
		/* Too many prime factors of bases.  Bailing out. */
		goto exit;
	}
	a = b = xsqr = ysqr = 1;
	for (i = 0; i < nprimes; i++) {
		if (exps[i] % 2 == 1) {
			a *= primes[i];
		}
		xsqr *= ipwr(primes[i], exps[i]/2);	/* May truncate */
	}
	mpz_set(r, y);
	nprimes = trialdiv(r, primes, exps, 31, 6);
	if (mpz_cmp_ui(r, 1) != 0) {
		/* Too many prime factors of bases.  Bailing out. */
		goto exit;
	}
	for (i = 0; i < nprimes; i++) {
		if (exps[i] % 2 == 1) {
			b *= primes[i];
		}
		ysqr *= ipwr(primes[i], exps[i]/2);	/* May truncate */
	}
	sbp = a * b;
	
	/* Got the SBP.  Is there an Aurifeuillian factorization? */
	if ((sbp % 4 == 1 && sign > 0) || (sbp % 4 != 1 && sign < 0)) {
		/* We're on the wrong side. */
		goto exit;
	}
	if (n % sbp != 0 || (n / sbp) % 2 != 1) {
		/* Exponent isn't an odd multiple of SBP. */
		goto exit;
	}
	if (sbp == 19 || sbp == 22 || sbp == 23 || sbp == 26 || sbp == 29 || sbp > 30) {
		/* We can't make a polynomial with a reasonable degree. */
		goto exit;
	}
	
	/*
	 * We now know there should be an Aurifeuillian factorization.
	 * Let's see if we can use it.
	 */
	LM = find_aurifeuillian(N, x, y, n, sbp, xsqr * ysqr);
	if (LM == 0) {
		/* N didn't divide either form, bail out. */
		goto exit;
	}
	h = n / sbp;
	k = h;
	/* There are 52 distinct cases, separated by SBP */
	switch (sbp) {
		case 2:
			/* 17 cases */
			if (h % 15 == 0 && (d == 0 || d == 8)) {
				/* Case 1: Use algebraic factors for both 3 and 5 */
				d = 8;
				k = h / 15;
				mpz_set_si(f[8], 1);
				mpz_set_si(f[7], LM * 2);
				mpz_set_si(f[6], -14);
				mpz_set_si(f[5], -LM * 24);
				mpz_set_si(f[4], 64);
				mpz_set_si(f[3], LM * 88);
				mpz_set_si(f[2], -96);
				mpz_set_si(f[1], -LM * 96);
				mpz_set_si(f[0], 16);
				mode = MODE_HALF;
			} else if (h % 9 == 0 && (d == 0 || d == 6)) {
				/*
				 * Case 2: Use algebraic factor of 3, and play
				 * tricks to make it a sextic.
				 */
				d = 6;
				k = h / 9;
				mpz_set_si(f[6], 1);
				mpz_set_si(f[4], -12);
				mpz_set_si(f[3], LM * 4);
				mpz_set_si(f[2], 36);
				mpz_set_si(f[1], -LM * 24);
				mpz_set_si(f[0], -8);
				mode = MODE_HALF;
			} else if (h % 3 == 0 && (d == 0 || d == 4 || d == 8)) {
				/* Algebraic factor for 3 */
				k = h / 3;
				if (d == 0 || d == 4) {
					/* Case 3: quartic */
					d = 4;
					mpz_set_si(f[4], a * a);
					mpz_set_si(f[3], LM * 2 * a);
					mpz_set_si(f[2], 2);
					mpz_set_si(f[1], LM * 2 * b);
					mpz_set_si(f[0], b * b);
					mode = MODE_NATURAL;
				} else {
					/* Cases 4, 5: Need to expand it to an octic */
					if (k % 4 == 1) {
						mpz_set(f[8], x);
						mpz_mul(f[8], f[8], x);
						mpz_set(f[6], x);
						mpz_set(f[2], y);
						mpz_set(f[0], y);
						mpz_mul(f[0], f[0], y);
						exp_offset = -1;
					} else {
						mpz_set(f[8], y);
						mpz_mul(f[8], f[8], y);
						mpz_set(f[6], y);
						mpz_set(f[2], x);
						mpz_set(f[0], x);
						mpz_mul(f[0], f[0], x);
						exp_offset = 1;
					}
					mpz_mul_si(f[6], f[6], LM * 2 * xsqr * ysqr);
					mpz_set(f[4], x);
					mpz_mul(f[4], f[4], y);
					mpz_mul_si(f[2], f[2], LM * 2 * xsqr * ysqr);
					exp_factor = 2;
					mode = MODE_EXPAND;
				}
			} else if (h % 5 == 0 && (d == 0 || d == 4 || d == 8)) {
				/* Algebraic factor for 5 */
				k = h / 5;
				if (d == 0 || d == 4) {
					/* Case 6: quartic */
					d = 4;
					mpz_set_si(f[4], 1);
					mpz_set_si(f[3], LM * 2);
					mpz_set_si(f[2], -6);
					mpz_set_si(f[1], -LM * 12);
					mpz_set_si(f[0], -4);
					mode = MODE_HALF;
				} else {
					/* Case 7: octic */
					mpz_set_si(f[8], a * a * a * a);
					mpz_set_si(f[7], LM * 2 * a * a * a);
					mpz_set_si(f[6], 2 * a * a);
					mpz_set_si(f[4], -4);
					mpz_set_si(f[2], 2 * b * b);
					mpz_set_si(f[1], LM * 2 * b * b * b);
					mpz_set_si(f[0], b * b * b * b);
					mode = MODE_NATURAL;
				}
			} else if (h % 7 == 0 && (d == 0 || d == 6)) {
				/* Case 8: algebraic factor for 7 */
				d = 6;
				k = h / 7;
				mpz_set_si(f[6], 1);
				mpz_set_si(f[5], -LM * 2);
				mpz_set_si(f[4], -10);
				mpz_set_si(f[3], LM * 20);
				mpz_set_si(f[2], 16);
				mpz_set_si(f[1], -LM * 32);
				mpz_set_si(f[0], 8);
				mode = MODE_HALF;
			} else {
				/* Can't use algebraic factors.  What degree shall we use? */
				if (d == 4) {
					/* Cases 9, 10: quartic */
					if (h % 4 == 1) {
						mpz_set(f[4], x);
						mpz_set(f[0], y);
						exp_offset = -1;
					} else {
						mpz_set(f[4], y);
						mpz_set(f[0], x);
						exp_offset = 1;
					}
					mpz_set_si(f[2], LM * 2 * xsqr * ysqr);
					exp_factor = 2;
					mode = MODE_EXPAND;
				} else if (d == 0 || d == 6) {
					d = 6;
					if (h % 6 == 1 || h % 6 == 5) {
						if (h % 6 == 1) {
							/*
							 * Cases 11, 12: sextic. These two
							 * cases have a lot in common, unlike when
							 * h == 3 (mod 6).
							 */
							mpz_set(f[6], x);
							mpz_set(f[0], y);
							exp_offset = -1;
						} else {
							mpz_set(f[6], y);
							mpz_set(f[0], x);
							exp_offset = 1;
						}
						mpz_set_si(f[3], LM * 2 * xsqr * ysqr);
						exp_factor = 3;
						mode = MODE_EXPAND;
					} else {
						/*
						 * Case 13: sextic.
						 * This really shouldn't happen. It would mean that
						 * h is a multiple of 3, so we should be using
						 * the algebraic factor here, which we want to catch
						 * above.  But if h is not a multiple of 9 and the
						 * user is demanding a sextic anyway, then we will
						 * end up here.  It also has unusual algebraic
						 * properties: we're expanding, but this behaves more
						 * like a natural sextic with k = h/3, so we use
						 * that mode.
						 */
						k = h / 3;
						mpz_set_si(f[6], a * a * a);
						mpz_set_si(f[3], LM * 4);
						mpz_set_si(f[0], b * b * b);
						mode = MODE_NATURAL;
					}
				} else if (d == 8) {
					/* Cases 14-17: octic */
					if (h % 8 == 1) {
						mpz_set(f[8], x);
						mpz_set_si(f[4], 1);
						mpz_set(f[0], y);
						exp_offset = -1;
					} else if (h % 8 == 3) {
						mpz_set(f[8], x);
						mpz_mul(f[8], f[8], x);
						mpz_mul(f[8], f[8], x);
						mpz_set(f[4], x);
						mpz_mul(f[4], f[4], y);
						mpz_set(f[0], y);
						mpz_mul(f[0], f[0], y);
						mpz_mul(f[0], f[0], y);
						exp_offset = -3;
					} else if (h % 8 == 5) {
						mpz_set(f[8], y);
						mpz_mul(f[8], f[8], y);
						mpz_mul(f[8], f[8], y);
						mpz_set(f[4], x);
						mpz_mul(f[4], f[4], y);
						mpz_set(f[0], x);
						mpz_mul(f[0], f[0], x);
						mpz_mul(f[0], f[0], x);
						exp_offset = 3;
					} else {
						mpz_set(f[8], y);
						mpz_set_si(f[4], 1);
						mpz_set(f[0], x);
						exp_offset = 1;
					}
					mpz_mul_si(f[4], f[4], LM * 2 * xsqr * ysqr);
					exp_factor = 4;
					mode = MODE_EXPAND;
				} else {
					goto exit;
				}
			}
			break;
		case 3:
			/* 12 cases */
			if (h % 5 == 0 && (d == 0 || d == 4)) {
				/* Case 1: quartic with algebraic factor of 5 */
				d = 4;
				k = h / 5;
				mpz_set_si(f[4], 1);
				mpz_set_si(f[3], LM * 3);
				mpz_set_si(f[2], -6);
				mpz_set_si(f[1], -LM * 18);
				mpz_set_si(f[0], -9);
				mode = MODE_HALF;
			} else if (h % 5 == 0 && d == 8) {
				/* Case 2: octic with algebraic factor of 5 */
				k = h / 5;
				mpz_set_si(f[8], a * a * a * a);
				mpz_set_si(f[7], LM * 3 * a * a * a);
				mpz_set_si(f[6], 6 * a * a);
				mpz_set_si(f[5], LM * 9 * a);
				mpz_set_si(f[4], 9);
				mpz_set_si(f[3], LM * 9 * b);
				mpz_set_si(f[2], 6 * b * b);
				mpz_set_si(f[1], LM * 3 * b * b * b);
				mpz_set_si(f[0], b * b * b * b);
				mode = MODE_NATURAL;
			} else if (h % 7 == 0 && (d == 0 || d == 6)) {
				/* Case 3: sextic with algebraic factor of 7 */
				d = 6;
				k = h / 7;
				mpz_set_si(f[6], 1);
				mpz_set_si(f[5], LM * 3);
				mpz_set_si(f[4], -12);
				mpz_set_si(f[3], -LM * 36);
				mpz_set_si(f[2], 18);
				mpz_set_si(f[1], LM * 54);
				mpz_set_si(f[0], -27);
				mode = MODE_HALF;
			} else {
				/* Can't use algebraic factor.  What degree? */
				if (d == 4) {
					/* Cases 4, 5: quartic */
					if (h % 4 == 1) {
						mpz_set(f[4], x);
						mpz_set(f[0], y);
						exp_offset = -1;
					} else {
						mpz_set(f[4], y);
						mpz_set(f[0], x);
						exp_offset = 1;
					}
					mpz_set_si(f[2], LM * 3 * xsqr * ysqr);
					exp_factor = 2;
					mode = MODE_EXPAND;
				} else if (d == 0 || d == 6) {
					d = 6;
					if (h % 6 == 1 || h % 6 == 5) {
						/*
						 * Cases 6, 7: sextic.  These cases share some
						 * features.  See comments above for cases 11, 12
						 * with SBP = 2.
						 */
						if (h % 6 == 1) {
							mpz_set(f[6], x);
							mpz_set(f[0], y);
							exp_offset = -1;
						} else {
							mpz_set(f[6], y);
							mpz_set(f[0], x);
							exp_offset = 1;
						}
						mpz_set_si(f[3], LM * 3 * xsqr * ysqr);
						exp_factor = 3;
						mode = MODE_EXPAND;
					} else {
						/*
						 * Case 8: sextic.  Unlike with SBP = 2, this
						 * is a perfectly reasonable case, since there can't
						 * be an algebraic factor of 3 when SBP = 3.  However
						 * the comment for that case about the unusual
						 * properties that make it more like a natural
						 * sextic do apply here as well.
						 */
						k = h / 3;
						mpz_set_si(f[6], a * a * a);
						mpz_set_si(f[3], LM * 9);
						mpz_set_si(f[0], b * b * b);
						mode = MODE_NATURAL;
					}
				} else if (d == 8) {
					/* Cases 9-12: octic */
					if (h % 8 == 1) {
						mpz_set(f[8], x);
						mpz_set_si(f[4], 1);
						mpz_set(f[0], y);
						exp_offset = -1;
					} else if (h % 8 == 3) {
						mpz_set(f[8], x);
						mpz_mul(f[8], f[8], x);
						mpz_mul(f[8], f[8], x);
						mpz_set(f[4], x);
						mpz_mul(f[4], f[4], y);
						mpz_set(f[0], y);
						mpz_mul(f[0], f[0], y);
						mpz_mul(f[0], f[0], y);
						exp_offset = -3;
					} else if (h % 8 == 5) {
						mpz_set(f[8], y);
						mpz_mul(f[8], f[8], y);
						mpz_mul(f[8], f[8], y);
						mpz_set(f[4], x);
						mpz_mul(f[4], f[4], y);
						mpz_set(f[0], x);
						mpz_mul(f[0], f[0], x);
						mpz_mul(f[0], f[0], x);
						exp_offset = 3;
					} else {
						mpz_set(f[8], y);
						mpz_set_si(f[4], 1);
						mpz_set(f[0], x);
						exp_offset = 1;
					}
					mpz_mul_si(f[4], f[4], LM * 3 * xsqr * ysqr);
					exp_factor = 4;
					mode = MODE_EXPAND;
				} else {
					goto exit;
				}
			}
			break;
		case 5:
			/* 5 cases */
			if (d != 0 && d != 4 && d != 8) goto exit;
			if (h % 3 == 0) {
				/* Use the algebraic factor.  Octic or quartic? */
				k = h / 3;
				if (d == 0 || d == 4) {
					/* Case 1: quartic with algebraic factor of 3 */
					d = 4;
					mpz_set_si(f[4], 1);
					mpz_set_si(f[3], LM * 5);
					mpz_set_si(f[2], -10);
					mpz_set_si(f[1], -LM * 50);
					mpz_set_si(f[0], 25);
					mode = MODE_HALF;
				} else {
					/* Case 2: octic with algebraic factor of 3 */
					mpz_set_si(f[8], a * a * a * a);
					mpz_set_si(f[7], LM * 5 * a * a * a);
					mpz_set_si(f[6], 10 * a * a);
					mpz_set_si(f[5], LM * 25 * a);
					mpz_set_si(f[4], 75);
					mpz_set_si(f[3], LM * 25 * b);
					mpz_set_si(f[2], 10 * b * b);
					mpz_set_si(f[1], LM * 5 * b * b * b);
					mpz_set_si(f[0], b * b * b * b);
					mode = MODE_NATURAL;
				}
			} else {
				/* No algebraic factor.  Octic or quartic? */
				if (d == 0 || d == 4) {
					/* Case 3: quartic */
					d = 4;
					mpz_set_si(f[4], a * a);
					mpz_set_si(f[3], LM * 5 * a);
					mpz_set_si(f[2], 15);
					mpz_set_si(f[1], LM * 5 * b);
					mpz_set_si(f[0], b * b);
					mode = MODE_NATURAL;
				} else {
					/* Cases 4, 5: octic */
					if (h % 4 == 1) {
						mpz_set(f[8], x);
						mpz_mul(f[8], f[8], x);
						mpz_set(f[6], x);
						mpz_set(f[2], y);
						mpz_set(f[0], y);
						mpz_mul(f[0], f[0], y);
						exp_offset = -1;
					} else {
						mpz_set(f[8], y);
						mpz_mul(f[8], f[8], y);
						mpz_set(f[6], y);
						mpz_set(f[2], x);
						mpz_set(f[0], x);
						mpz_mul(f[0], f[0], x);
						exp_offset = 1;
					}
					mpz_mul_si(f[6], f[6], LM * 5 * xsqr * ysqr);
					mpz_set_si(f[4], 3);
					mpz_mul(f[4], f[4], x);
					mpz_mul(f[4], f[4], y);
					mpz_mul_si(f[2], f[2], LM * 5 * xsqr * ysqr);
					exp_factor = 2;
					mode = MODE_EXPAND;
				}
			}
			break;
		case 6:
			/* 5 cases */
			if (h % 5 == 0 && (d == 0 || d == 8)) {
				/* Case 1: octic with algebraic factor of 5 */
				d = 8;
				k = h / 5;
				mpz_set_si(f[8], 1);
				mpz_set_si(f[7], -LM * 6);
				mpz_set_si(f[6], -30);
				mpz_set_si(f[5], LM * 216);
				mpz_set_si(f[4], 144);
				mpz_set_si(f[3], -LM * 1944);
				mpz_set_si(f[2], 0);
				mpz_set_si(f[1], LM * 5184);
				mpz_set_si(f[0], 1296);
				mode = MODE_HALF;
			} else if (h % 3 == 0 && (d == 0 || d == 6)) {
				/*
				 * Case 2: no algebraic factor here, but when h is a multiple
				 * of 3 there's a trick we can use to make a sextic.
				 */
				d = 6;
				k = h / 3;
				mpz_set_si(f[6], 1);
				mpz_set_si(f[4], -36);
				mpz_set_si(f[3], LM * 36);
				mpz_set_si(f[2], 324);
				mpz_set_si(f[1], -LM * 648);
				mpz_set_si(f[0], 216);
				mode = MODE_HALF;
			} else {
				/* No algebraic, or sextic was not allowed. Octic or quartic? */
				if (d == 0 || d == 4) {
					/* Case 3: quartic */
					d = 4;
					mpz_set_si(f[4], a * a);
					mpz_set_si(f[3], LM * 6 * a);
					mpz_set_si(f[2], 18);
					mpz_set_si(f[1], LM * 6 * b);
					mpz_set_si(f[0], b * b);
					mode = MODE_NATURAL;
				} else if (d == 8) {
					/* Cases 4, 5: octic, expanded from quartic */
					if (h % 4 == 1) {
						mpz_set(f[8], x);
						mpz_mul(f[8], f[8], x);
						mpz_set(f[6], x);
						mpz_set_si(f[4], 3);
						mpz_mul(f[4], f[4], x);
						mpz_mul(f[4], f[4], y);
						mpz_set(f[2], y);
						mpz_set(f[0], y);
						mpz_mul(f[0], f[0], y);
						exp_offset = -1;
					} else {
						mpz_set(f[8], y);
						mpz_mul(f[8], f[8], y);
						mpz_set(f[6], y);
						mpz_set_si(f[4], 3);
						mpz_mul(f[4], f[4], x);
						mpz_mul(f[4], f[4], y);
						mpz_set(f[2], x);
						mpz_set(f[0], x);
						mpz_mul(f[0], f[0], x);
						exp_offset = 1;
					}
					mpz_mul_si(f[6], f[6], LM * 6 * xsqr * ysqr);
					mpz_mul_si(f[2], f[2], LM * 6 * xsqr * ysqr);
					exp_factor = 2;
					mode = MODE_EXPAND;
				} else {
					goto exit;
				}
			}
			break;
		case 7:
			/* 2 cases */
			if (d != 0 && d != 6) goto exit;
			d = 6;
			if (h % 3 == 0) {
				k = h / 3;
				mpz_set_si(f[6], 1);
				mpz_set_si(f[5], -LM * 7);
				mpz_set_si(f[4], -14);
				mpz_set_si(f[3], LM * 196);
				mpz_set_si(f[2], -392);
				mpz_set_si(f[0], 343);
				mode = MODE_HALF;
			} else {
				mpz_set_si(f[6], a * a * a);
				mpz_set_si(f[5], LM * 7 * a * a);
				mpz_set_si(f[4], 21 * a);
				mpz_set_si(f[3], LM * 49);
				mpz_set_si(f[2], 21 * b);
				mpz_set_si(f[1], LM * 7 * b * b);
				mpz_set_si(f[0], b * b * b);
				mode = MODE_NATURAL;
			}
			break;
		case 10:
			/* 3 cases */
			if ((d == 0 || d == 8) && h % 3 == 0) {
				/* Use the algebraic factor */
				if (d != 0 && d != 8) goto exit;
				d = 8;
				k = h / 3;
				mpz_set_si(f[8], 1);
				mpz_set_si(f[7], -LM * 10);
				mpz_set_si(f[6], -30);
				mpz_set_si(f[5], LM * 600);
				mpz_set_si(f[4], -1200);
				mpz_set_si(f[3], -LM * 7000);
				mpz_set_si(f[2], 32000);
				mpz_set_si(f[1], -LM * 40000);
				mpz_set_si(f[0], 10000);
				mode = MODE_HALF;
			} else {
				/*
				 * No algebraic factors, or a non-octic was requested.
				 * Octic or quartic?
				 */
				if (d == 0 || d == 4) {
					d = 4;
					mpz_set_si(f[4], 1);
					mpz_set_si(f[3], LM * 10);
					mpz_set_si(f[2], 10);
					mpz_set_si(f[1], -LM * 100);
					mpz_set_si(f[0], -100);
					mode = MODE_HALF;
				} else if (d == 8) {
					d = 8;
					mpz_set_si(f[8], a * a * a * a);
					mpz_set_si(f[7], LM * 10 * a * a * a);
					mpz_set_si(f[6], 50 * a * a);
					mpz_set_si(f[5], LM * 200 * a);
					mpz_set_si(f[4], 700);
					mpz_set_si(f[3], LM * 200 * b);
					mpz_set_si(f[2], 50 * b * b);
					mpz_set_si(f[1], LM * 10 * b * b * b);
					mpz_set_si(f[0], b * b * b * b);
					mode = MODE_NATURAL;
				} else {
					goto exit;
				}
			}
			break;
		case 11:
			/* Single case, no algebraic factors */
			if (d != 0 && d != 5) goto exit;
			d = 5;
			mpz_set_si(f[5], 1);
			mpz_set_si(f[4], LM * 11);
			mpz_set_si(f[2], -LM * 363);
			mpz_set_si(f[1], -1331);
			mpz_set_si(f[0], -LM * 1331);
			mode = MODE_HALF;
			break;
		case 13:
			/* Single case, no algebraic factors */
			if (d != 0 && d != 6) goto exit;
			d = 6;
			mpz_set_si(f[6], 1);
			mpz_set_si(f[5], LM * 13);
			mpz_set_si(f[4], 13);
			mpz_set_si(f[3], -LM * 338);
			mpz_set_si(f[2], -676);
			mpz_set_si(f[1], LM * 2197);
			mpz_set_si(f[0], 2197);
			mode = MODE_HALF;
			break;
		case 14:
			/* Single case, no algebraic factors */
			if (d != 0 && d != 6) goto exit;
			d = 6;
			mpz_set_si(f[6], 1);
			mpz_set_si(f[5], LM * 14);
			mpz_set_si(f[4], 14);
			mpz_set_si(f[3], -LM * 588);
			mpz_set_si(f[2], -3136);
			mpz_set_si(f[1], -LM * 5488);
			mpz_set_si(f[0], -2744);
			mode = MODE_HALF;
			break;
		case 15:
			/* 2 cases.  Octic or quartic (default to quartic) */
			if (d == 0 || d == 4) {
				d = 4;
				mpz_set_si(f[4], 1);
				mpz_set_si(f[3], LM * 15);
				mpz_set_si(f[2], 60);
				mpz_set_si(f[0], -225);
				mode = MODE_HALF;
			} else if (d == 8) {
				mpz_set_si(f[8], a * a * a * a);
				mpz_set_si(f[7], LM * 15 * a * a * a);
				mpz_set_si(f[6], 120 * a * a);
				mpz_set_si(f[5], LM * 675 * a);
				mpz_set_si(f[4], 2925);
				mpz_set_si(f[3], LM * 675 * b);
				mpz_set_si(f[2], 120 * b * b);
				mpz_set_si(f[1], LM * 15 * b * b * b);
				mpz_set_si(f[0], b * b * b * b);
				mode = MODE_NATURAL;
			} else {
				goto exit;
			}
			break;
		case 17:
			/* Single case, no algebraic factors */
			if (d != 0 && d != 8) goto exit;
			d = 8;
			mpz_set_si(f[8], 1);
			mpz_set_si(f[7], LM * 17);
			mpz_set_si(f[6], 17);
			mpz_set_si(f[5], -LM * 1156);
			mpz_set_si(f[4], -6647);
			mpz_set_si(f[2], 78608);
			mpz_set_si(f[1], LM * 167042);
			mpz_set_si(f[0], 83521);
			mode = MODE_HALF;
			break;
		case 21:
			/* Single case, no algebraic factors */
			if (d != 0 && d != 6) goto exit;
			d = 6;
			mpz_set_si(f[6], 1);
			mpz_set_si(f[5], LM * 21);
			mpz_set_si(f[4], 84);
			mpz_set_si(f[3], -LM * 882);
			mpz_set_si(f[2], -7938);
			mpz_set_si(f[1], -LM * 18522);
			mpz_set_si(f[0], -9261);
			mode = MODE_HALF;
			break;
		case 30:
			/* Single case, no algebraic factors */
			if (d != 0 && d != 8) goto exit;
			d = 8;
			mpz_set_si(f[8], 1);
			mpz_set_si(f[7], LM * 30);
			mpz_set_si(f[6], 210);
			mpz_set_si(f[5], -LM * 1800);
			mpz_set_si(f[4], -28800);
			mpz_set_si(f[3], -LM * 81000);
			mpz_set_si(f[2], 324000);
			mpz_set_si(f[1], LM * 1620000);
			mpz_set_si(f[0], 810000);
			mode = MODE_HALF;
			break;
		default:
			goto exit;
	}
	
	switch (mode) {
		case MODE_HALF:
			mpz_pow_ui(g[0], x, k);
			mpz_pow_ui(r, y, k);
			mpz_add(g[0], g[0], r);
			mpz_pow_ui(g[1], x, (k-1)/2);
			mpz_pow_ui(r, y, (k-1)/2);
			mpz_mul(g[1], g[1], r);
			mpz_mul_ui(g[1], g[1], xsqr * ysqr);
			/*
			 * Note the integer division here.  Usually d is even,
			 * in which case (d+1)/2 is just d/2, so we're just
			 * taking the square root of sbp.  But if d is odd
			 * (currently only when sbp is 11), then the correct
			 * exponent depends on d.
			 */
			*skewness = pow((double)sbp, (double)((d+1)/2) / d);
			*difficulty = log10(dx) * d * k;
			break;
		case MODE_NATURAL:
			mpz_pow_ui(g[0], x, (k-1)/2);
			mpz_mul_ui(g[0], g[0], xsqr);
			mpz_pow_ui(g[1], y, (k-1)/2);
			mpz_mul_ui(g[1], g[1], ysqr);
			*skewness = pow((double)b / a, 0.5);
			*difficulty = log10(dx) * d/2 * k;
			break;
		case MODE_EXPAND:
			mpz_pow_ui(g[0], x, (k+exp_offset)/(2*exp_factor));
			mpz_pow_ui(g[1], y, (k+exp_offset)/(2*exp_factor));
			*skewness = pow(dy / dx, (double)(-exp_offset) / (2*exp_factor));
			if (exp_offset < 0) {
				*difficulty = log10(dx) * (double)(d / (2*exp_factor)) * k;
			} else {
				*difficulty = (double)(d / (2*exp_factor)) *
					(log10(dx) * (k + exp_offset) + log10(dy) * exp_offset);
			}
			break;
	}
	mpz_neg(g[1], g[1]);
	*degree = d;
	*lm = LM;
	success = 1;
	
exit:
	mpz_clear(r);
	return success;
}


int main(int argc, char **argv)
{
	int format = GGNFS;		/* Output format */
	mpz_t x;
	mpz_t y;
	int n = 0;			/* We factor phi_n(x,y) */
	int degree = 0;
	int sign, i;
	int shortform = 0;		/* Short output for Franke format? */
	int read_stdin = 0;		/* Read number from stdin? */
	int check_aurif = 1;
	int lm;
	int aurif = 0;
	int use_auto = 0;
	mpz_t N;			/* The cofactor of the number to factor */
	mpz_t f[MAXDEGREE+1];/* Coefficients on algebraic side, max degree = 8 */
	mpz_t g[2];			/* Coefficients on rational side, max degree = 1 */
	mpz_t M;			/* Common root */
	mpz_t t;			/* Numerator of M */
	mpz_t u;			/* Denominator of M */
	mpz_t r;
	mpz_t s;			/* Temporaries */
	double difficulty;
	double skewness;
	const double timefactor = 1. / 2.1e15;
	char aurif_str[] = " Aurifeuillian L";
	
	/* Parse command line parameters */
	while (argc > 1 && argv[1][0] == '-') {
		if (strcmp(argv[1], "-cwi") == 0) {
			format = CWI;
			argc--;
			argv++;
		} else if (strcmp(argv[1], "-franke") == 0) {
			format = FRANKE;
			argc--;
			argv++;
		} else if (strcmp(argv[1], "-ggnfs") == 0) {
			format = GGNFS;
			argc--;
			argv++;
		} else if (strcmp(argv[1], "-short") == 0) {
			shortform = 1;
			argc--;
			argv++;
		} else if (strcmp(argv[1], "-deg4") == 0) {
			degree = 4;
			argc--;
			argv++;
		} else if (strcmp(argv[1], "-deg5") == 0) {
			degree = 5;
			argc--;
			argv++;
		} else if (strcmp(argv[1], "-deg6") == 0) {
			degree = 6;
			argc--;
			argv++;
		} else if (strcmp(argv[1], "-deg8") == 0) {
			degree = 8;
			argc--;
			argv++;
		} else if (strcmp(argv[1], "-noaurif") == 0) {
			check_aurif = 0;
			argc--;
			argv++;
		} else if (strcmp(argv[1], "-") == 0) {
			read_stdin = 1;
			argc--;
			argv++;
		} else {
			fprintf(stderr, "Error, unknown command line option: %s\n",
					argv[1]);
			usage();
			exit(EXIT_FAILURE);
		}
	}
	
	mpz_init(N);
	mpz_init(x);
	mpz_init(y);
	mpz_init(s);
	mpz_init(r);
	
	if (argc + read_stdin == 2) {
		int xi = 0;
		int yi = 0;
		if (read_stdin)
			mpz_inp_str(N, stdin, 10);
		else
			mpz_set_str(N, argv[1], 10);
		n = 0;
		if (autodet(&xi, &yi, (uint *)&n, N)) {
			mpz_set_ui(x, xi);
			mpz_set_ui(y, yi);
#ifdef DEBUG
			printf("Auto-detected Phi_%d(%d,%d)\n", n, xi, yi);
		} else {
			printf("Auto-detection failed\n");
#endif
		}
		use_auto = 1;
	} else if (argc + read_stdin > 3) {
		n = atoi(argv[1]);
		mpz_set_str(x, argv[2], 10);
		if (argc + read_stdin == 4) {
			mpz_set_ui(y, 1);
		} else {
			mpz_set_str(y, argv[3], 10);
		}
		if (read_stdin)
			mpz_inp_str(N, stdin, 10);
		else {
			if (argc == 4)
				mpz_set_str(N, argv[3], 10);
			else
				mpz_set_str(N, argv[4], 10);
		}
	}
	if (n <= 0 || mpz_sgn(x) <= 0 || mpz_sgn(y) <= 0 || mpz_sgn(N) <= 0) {
		if (use_auto) {
			printf("Unable to auto-detect cyclotomic form for this number.\n");
			printf("Try providing values for n, x and y.\n");
		} else {
			gmp_printf("Error, invalid parameters: n = %d, x = %Zd, y = %Zd\n",
					   n, x, y);
			usage();
		}
		exit(EXIT_FAILURE);
	}
	mpz_gcd(r, x, y);
	if (mpz_cmp_si(r, 1) > 0) {
		printf("Error, GCD(x,y) != 1\n");
		exit(EXIT_FAILURE);
	}
	if (mpz_cmp(x, y) <= 0)
		mpz_swap(x, y);
	mpz_init(t);
	mpz_init(u);
	mpz_init(g[0]);
	mpz_init(g[1]);
	for (i = 0; i < MAXDEGREE + 1; i++)
		mpz_init(f[i]);
	mpz_init(M);
	
	sign = -1;
	/* Simplify even n, that is x^(n/2)+y^(n/2) numbers */
	if (n % 2 == 0) {
		sign = +1;
		n /= 2;
	}
	/* Thus the number we factor is x^n+sign*y^n */
	
	/* Check that N | x^n+sign*y^n */
	mpz_powm_ui(r, x, n, N);
	mpz_powm_ui(s, y, n, N);
	(sign > 0 ? mpz_add : mpz_sub) (r, r, s);
	mpz_mod(r, r, N);
	if (mpz_sgn(r) != 0) {
		gmp_fprintf(stderr, "Error, N does not divide %Zd^%d%c%Zd^%d\n",
					x, n, sign < 0 ? '-' : '+', y, n);
		exit(EXIT_FAILURE);
	}
	
	/* Find out what kind of number we have (algebraic factors etc.) and
	 choose polynomial */
	
	if (check_aurif) {
		aurif = generate_aurif(N, x, y, n, sign, &degree, f, g, &difficulty,
							   &skewness, &lm);
		if (aurif && lm == 1) {
			aurif_str[strlen(aurif_str)-1] = 'M';
		}
	}
	if (!aurif) {
		aurif_str[0] = '\0';
		generate_std(N, x, y, n, sign, &degree, f, g, &difficulty, &skewness);
	}
	mpz_neg(M, g[1]);
	mpz_invert(M, M, N);
	mpz_mul(M, M, g[0]);
	mpz_mod(M, M, N);
	
	
	/* Test polynomials (check that values vanish (mod N) at root M) */
	mpz_set_ui(r, 0);
	for (i = degree; i >= 0; i--) {
		mpz_mul(r, r, M);
		mpz_add(r, r, f[i]);
		mpz_mod(r, r, N);
	}
	if (mpz_sgn(r) != 0) {
		gmp_fprintf(stderr, "Error: M=%Zd is not a root of f(x) % N\n", M);
		fprintf(stderr, "f(x) = ");
		for (i = degree; i >= 0; i--)
			gmp_fprintf(stderr, "%Zd*x^%d +", f[i], i);
		gmp_fprintf(stderr, "\n" "Remainder is %Zd\n", r);
		gmp_fprintf(stderr, "%Zd^%d%c%Zd^%d\n", x, n, sign < 0 ? '-' : '+',
					y, n);
		fprintf(stderr, "Please report this bug.\n");
		exit(EXIT_FAILURE);
	}
	
	mpz_mul(r, g[1], M);
	mpz_add(r, r, g[0]);
	mpz_mod(r, r, N);
	if (mpz_sgn(r) != 0) {
		gmp_fprintf(stderr, "Error: M=%Zd is not a root of g(x) % N\n", M);
		gmp_fprintf(stderr, "Remainder is %Zd\n", r);
		fprintf(stderr, "Please report this bug.\n");
		exit(EXIT_FAILURE);
	}
	
	/* Output polynomials */
	if (format == FRANKE) {
		const double cost =
		snfscost(difficulty);
		
		gmp_printf("%Zd\n", N);
		gmp_printf("# %Zd^%d%c%Zd", x, n, (sign > 0 ? '+' : '-'), y);
		if (mpz_cmp_ui(y, 1) > 0) {
			printf("^%d", n);
		}
		gmp_printf
		("%s, difficulty: %.2f, skewness: %.2f\n"
		 "# Cost: %g, est. time: %.2f GHz days (not accurate yet!)\n",
		 aurif_str, difficulty, skewness, cost, cost * timefactor);
		
		for (i = degree; i >= 0; i--)
			if (!shortform || mpz_sgn(f[i]) != 0)
				gmp_printf("X%d %Zd\n", i, f[i]);
		
		gmp_printf("Y1 %Zd\n", g[1]);
		gmp_printf("Y0 %Zd\n", g[0]);
		gmp_printf("M %Zd\n", M);
	} else if (format == GGNFS) {
		const double cost = snfscost(difficulty);
		
		gmp_printf("n: %Zd\n", N);
		gmp_printf("# %Zd^%d%c%Zd", x, n, (sign > 0 ? '+' : '-'), y);
		if (mpz_cmp_ui(y, 1) > 0) {
			printf("^%d", n);
		}
		gmp_printf
		("%s, difficulty: %.2f, skewness: %.2f\n"
		 "# cost: %g, est. time: %.2f GHz days (not accurate yet!)\n",
		 aurif_str, difficulty, skewness, cost, cost * timefactor);
		
		printf("skew: %.5f\n", skewness);
		for (i = degree; i >= 0; i--)
			if (mpz_sgn(f[i]) != 0)
				gmp_printf("c%d: %Zd\n", i, f[i]);
		
		gmp_printf("Y1: %Zd\n", g[1]);
		gmp_printf("Y0: %Zd\n", g[0]);
		gmp_printf("m: %Zd\n", M);
		printf("type: snfs\n");
	} else {
		gmp_printf("%Zd\n", N);
		gmp_printf("%Zd\n\n2\n\n1\n", M);
		gmp_printf("%Zd %Zd\n", g[0], g[1]);
		printf("\n\n%d\n", degree);
		for (i = 0; i <= degree; i++)
			gmp_printf("%Zd ", f[i]);
		printf("\n");
	}
	
	mpz_clear(N);
	mpz_clear(x);
	mpz_clear(y);
	mpz_clear(s);
	mpz_clear(r);
	mpz_clear(t);
	mpz_clear(u);
	mpz_clear(g[0]);
	mpz_clear(g[1]);
	for (i = 0; i < MAXDEGREE + 1; i++)
		mpz_clear(f[i]);
	mpz_clear(M);
	return 0;
}
