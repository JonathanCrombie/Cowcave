// all source from msieve svn 839.  type "svn co -r 839 http://svn.code.sf.net/p/msieve/code/trunk svn839" to get it.

#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <gmp.h>
#include <string.h>
#include <values.h>

// msieve stuff
#define MAX_POLY_DEGREE 8
#define POLY_HEAP_SIZE 1
#define PRECOMPUTED_NUM_PRIMES 9592
#define SMALL_PRIME_BOUND 100
#define ROOT_BUF_SIZE 200
#define NUM_POLY_COEFFS (MAX_POLY_DEGREE+1)
#define MAX_ROOTFINDER_DEGREE 10
#define NEWTON_ITER 10
#define MOD(i) buf[NUM_POLY_COEFFS+i]
#define OP1(i) buf[i]
#define SIZE_EPS 1e-6
#define MAX_MP_WORDS 32
#define POSITIVE 0
#define NEGATIVE 1
#define PRIME_BOUND 2000
#define mp_is_zero(a) ((a)->nwords == 0)
#define MP_RADIX 4294967296.0
#define MAX_TRAP_LEVELS 7
#define MAX_COEFFS 55
#define MAX_LINE 30
#define DICKMAN_ACCURACY 1e-8
#define INTEGRATE_LIMIT 1e12
#define STAGE1_ITER 5
#define STAGE3_ITER 10
#define MAX_LEVELS 30
#define MAX(a,b) ((a) > (b)? (a) : (b))
#define MIN(a,b) ((a) < (b)? (a) : (b))
#define HEAP_SWAP(a,b) { tmp = a; a = b; b = tmp; }
#define HEAP_PARENT(i)  (((i)-1) >> 1)
#define HEAP_LEFT(i)    (2 * (i) + 1)
#define HEAP_RIGHT(i)   (2 * (i) + 2)
#define FUNC_EVAL(base, off) func(base, off, params)
#define NUM_RECIP_FACTORIALS (sizeof(recip_factorial) / sizeof(double))

#define QUICK_TWO_SUM(a, b, sum, err) {		\
	double u = (a) + (b);			\
	(err) = (b) - (u - (a));		\
	(sum) = u;				\
}

#define TWO_SUM(a, b, sum, err) { 		\
	double u = (a) + (b);			\
	double bb = u - (a);			\
	(sum) = u;				\
	(err) = ((a) - (u - bb)) + ((b) - bb);	\
}

#define SPLIT_CONST 134217729.0               /* = 2^27 + 1 */

#define SPLIT(a, hi, lo) {		\
	double u = SPLIT_CONST * (a);\
	(hi) = u - (u - (a));	\
	(lo) = (a) - (hi);		\
}

#define TWO_PROD(a, b, prod, err) {		\
	double a_hi, a_lo, b_hi, b_lo;		\
	double u = (a) * (b);			\
	SPLIT((a), a_hi, a_lo);			\
	SPLIT((b), b_hi, b_lo);			\
	(prod) = u;				\
	(err) = ((a_hi * b_hi - u) + 		\
		a_hi * b_lo + a_lo * b_hi) + 	\
		a_lo * b_lo;			\
}

#define TWO_DIFF(a, b, sum, err) { 		\
	double u = (a) - (b);			\
	double bb = u - (a);			\
	(sum) = u;				\
	(err) = ((a) - (u - bb)) - ((b) + bb);	\
}

typedef struct {
	double hi;
	double lo;
} dd_t;

typedef struct {
	uint32_t nwords;		/* number of nonzero words in val[] */
	uint32_t val[MAX_MP_WORDS];
} mp_t;

typedef struct {
	uint32_t sign;	/* POSITIVE or NEGATIVE */
	mp_t num;
} signed_mp_t;

typedef struct {
	uint32_t degree;
	double coeff[MAX_POLY_DEGREE + 1];
} dpoly_t;

typedef struct {
	uint32_t degree;
	dd_t coeff[MAX_POLY_DEGREE + 1];
} ddpoly_t;

typedef uint16_t dd_precision_t;

typedef struct {
	dd_t r, i;
} dd_complex_t;

typedef struct {
	double r, i;
} complex_t;

typedef struct {
	complex_t poly[MAX_ROOTFINDER_DEGREE + 1]; 
	complex_t poly_aux[MAX_ROOTFINDER_DEGREE + 1];
	complex_t hpoly[MAX_ROOTFINDER_DEGREE + 1]; 
	complex_t hpoly_aux[MAX_ROOTFINDER_DEGREE + 1];

	complex_t angle, angle_inc;
	uint32_t hpoly_root_found;
	uint32_t degree;
} jt_t;

typedef struct {
	double power;
	ddpoly_t *rpoly;
	ddpoly_t *apoly;
} dparam_t;

typedef struct {
	uint32_t num_coeffs;
	uint32_t coeff_offset;
} dickman_line_t;

typedef struct {
	dickman_line_t *lines;
	double *coeffs;
} dickman_t;

typedef struct {
	dickman_t *dickman_aux;
	double root_score_r;
	double root_score_a;
	double skew_x;
	double skew_y;
	double rfb_limit;
	double afb_limit;
	double log_rfb_limit;
	double log_afb_limit;
	double sieve_area;
	dpoly_t rpoly;
	dpoly_t apoly;
} murphy_param_t;

typedef double (*integrand_t)(double, double, void *);

enum integrator_type {
	double_exponential
};

typedef struct {
	double result;		/* value of last integration */
	double error;		/* absolute error in last integration */

	enum integrator_type type;/* type of integrator */
	void *internal;		/* integrator-specific data */
} integrate_t;


/* representation of polynomials with multiple-
   precision coefficients. For polynomial p(x),
   element i of coeff[] gives the coefficient
   of x^i */

typedef struct {
	uint32_t degree;
	mpz_t coeff[MAX_POLY_DEGREE + 1];

	/* scratch quantities for evaluating the homogeneous
	   form of poly */
	mpz_t tmp1, tmp2, tmp3;
} mpz_poly_t;

typedef struct {
	mpz_poly_t rpoly;
	mpz_poly_t apoly;
	double size_score;
	double root_score;
	double combined_score;
	double skewness;
	uint32_t num_real_roots;
} poly_select_t;


typedef struct {
	poly_select_t *heap[POLY_HEAP_SIZE];
	uint32_t heap_num_filled;

	integrate_t integ_aux;
	dickman_t dickman_aux;
} poly_config_t;

typedef struct {
	uint32_t coef[2 * MAX_POLY_DEGREE + 1];
	uint32_t degree;
} _poly_t;

typedef _poly_t poly_t[1];


static void mp_clear(mp_t *a) {
	memset(a, 0, sizeof(mp_t));
}

/*---------- double exponential integrator, definite interval ------------*/

/* describes the integration of one subinterval */

typedef struct {
	double left;          /* left endpoint */
	double right;         /* right endpoint */
	double result;        /* numerical value of the integral */
	double error;         /* estimated error */
	uint32_t level;
} subinterval_t;

/* describes one point in the interval to evaluate */

typedef struct {
	double abscissa;        /* distance from the left or right endpoint
				   where the integrand is evaluated */
	double weight;          /* weight applied to the integrand at
				   'abscissa' before adding to integral value */
	double residual_weight; /* weight applied to the integrand at
				   'abscissa' before adding to residual */
} de_point_t;

/* structure controlling numerical integration. Note that introductory
   numerical methods describe the trapezoid rule as having equal-size
   panels, with the rule at half the step size containing twice as many
   half-size panels. The DE transform is different: we compute the
   integral from -1 to 1 via a transformation to an integral from -infinity
   to +infinity, then evaluate the latter integral using a series of
   trapezoid panels that increase super-exponentially in width. The 'step-size'
   refers to the distance of the first sample point from the origin, and
   halving the step-size here means doubling the number of trapezoid panels. 
   While there are really an infinite number of trapezoid panels to sum, 
   we truncate the sum when the size of the abscissas become too large to 
   be represented accurately. Because the abscissas increase in size so 
   quickly, it only takes a few samples to reach this point.  The de_point_t 
   structures store the value of the abscissa after conversion back to the 
   (-1,1) interval, as offsets from +-1 */

typedef struct {
	uint32_t num_points;         /* maximum number of trapezoid panels
				      allowed, beyond which the abscissas
				      are of size close to the limit of double 
				      precision */
	double relative_error;     /* the allowed relative error */
	de_point_t center_point;   /* initial sample point of the interval */
	de_point_t *initial_rule;  /* further trapezoid sample points at the 
				      coarsest step size (1.0). There are 
				      num_points of these, and since the rule 
				      is symmetric about the origin we have
				      up to 2*num_points samples at this
				      resolution */
	de_point_t *refined_rules; /* list of groups of num_points samples
				      that implement the trapezoid rule with
				      successively halved step sizes. For
				      level i, 0 <= i < num_trap_levels, the
				      step size (i.e. the distance of the
				      first abscissa to the origin) is
				      2^-(i+1) and there are 2^i groups of
				      up to num_points samples to add. Sample
				      groups are stored in order of increasing 
				      i, then increasing value of initial
				      abscissa. Each group of samples is also
				      symmetric about the origin.  */
	subinterval_t *heap;       /* used to store integrals of subintervals */
	uint32_t num_heap;           /* number of subintervals */
	uint32_t num_heap_alloc;     /* memory allocated for subintervals */
} de_t;


/*--------------------------------------------------------------------*/
static void * xmalloc(size_t len) {
	void *ptr = malloc(len);
	if (ptr == NULL) {
		printf("failed to allocate %u bytes\n", (uint32_t)len);
		exit(-1);
	}
	return ptr;
}

static void * xcalloc(size_t num, size_t len) {
	void *ptr = calloc(num, len);
	if (ptr == NULL) {
		printf("failed to calloc %u bytes\n", (uint32_t)(num * len));
		exit(-1);
	}
	return ptr;
}

static void * xrealloc(void *iptr, size_t len) {
	void *ptr = realloc(iptr, len);
	if (ptr == NULL) {
		printf("failed to reallocate %u bytes\n", (uint32_t)len);
		exit(-1);
	}
	return ptr;
}

static dd_t dd_set_dd(double hi, double lo) {
	dd_t dd;
	dd.hi = hi;
	dd.lo = lo;
	return dd;
}

static dd_t dd_set_d(double x) {
	dd_t dd;
	dd.hi = x;
	dd.lo = 0.0;
	return dd;
}

static dd_t dd_neg(dd_t a) {

	return dd_set_dd(-a.hi, -a.lo);
}

static int32_t dd_cmp_d(dd_t a, double b) {
	if (a.hi < b)
		return -1;
	if (a.hi > b)
		return 1;

	if (a.lo < 0)
		return -1;
	if (a.lo > 0)
		return 1;
	return 0;
}

static int32_t dd_cmp_dd(dd_t a, dd_t b) {
	if (a.hi < b.hi)
		return -1;
	if (a.hi > b.hi)
		return 1;

	if (a.lo < b.lo)
		return -1;
	if (a.lo > b.lo)
		return 1;
	return 0;
}

static dd_t dd_fabs(dd_t a) {
	if (a.hi < 0 || (a.hi == 0 && a.lo < 0))
		return dd_set_dd(-a.hi, -a.lo);
	else
		return a;
}

static dd_t dd_add_d(dd_t a, double b) {
	
	double s1, s2;

	TWO_SUM(a.hi, b, s1, s2);
	s2 += a.lo;
	QUICK_TWO_SUM(s1, s2, s1, s2);
	return dd_set_dd(s1, s2);
}

static dd_t dd_add_dd(dd_t a, dd_t b) {

	double s1, s2, t1, t2;

	TWO_SUM(a.hi, b.hi, s1, s2);
	TWO_SUM(a.lo, b.lo, t1, t2);
	s2 += t1;
	QUICK_TWO_SUM(s1, s2, s1, s2);
	s2 += t2;
	QUICK_TWO_SUM(s1, s2, s1, s2);
	return dd_set_dd(s1, s2);
}

static dd_t dd_mul_dpow2(dd_t a, double b) {

	return dd_set_dd(a.hi * b, a.lo * b);
}

static dd_t dd_mul_dd(dd_t a, dd_t b) {

	double p1, p2;

	TWO_PROD(a.hi, b.hi, p1, p2);
	p2 += (a.hi * b.lo);
	p2 += (a.lo * b.hi);
	QUICK_TWO_SUM(p1, p2, p1, p2);
	return dd_set_dd(p1, p2);
}

static dd_t dd_mp2dd(mp_t *x) {

	/* convert a multiple-precision number to 
	   a double-double */

	int32_t i;
	dd_t d;

	if (mp_is_zero(x))
		return dd_set_d(0.0);
	
	i = x->nwords - 1;
	d = dd_set_d((double)x->val[i]);
	for (i--; i >= 0; i--) {
		d = dd_mul_dpow2(d, MP_RADIX);
		d = dd_add_d(d, (double)x->val[i]);
	}

	return d;
}

static dd_t dd_signed_mp2dd(signed_mp_t *x) {

	if (x->sign == NEGATIVE)
		return dd_neg(dd_mp2dd(&x->num));
	else
		return dd_mp2dd(&x->num);
}

static dd_t dd_mul_d(dd_t a, double b) {

	double p1, p2;

	TWO_PROD(a.hi, b, p1, p2);
	p2 += (a.lo * b);
	QUICK_TWO_SUM(p1, p2, p1, p2);
	return dd_set_dd(p1, p2);
}

static dd_t dd_sub_dd(dd_t a, dd_t b) {

	double s1, s2, t1, t2;

	TWO_DIFF(a.hi, b.hi, s1, s2);
	TWO_DIFF(a.lo, b.lo, t1, t2);
	s2 += t1;
	QUICK_TWO_SUM(s1, s2, s1, s2);
	s2 += t2;
	QUICK_TWO_SUM(s1, s2, s1, s2);
	return dd_set_dd(s1, s2);
}

static dd_t dd_div_dd(dd_t a, dd_t b) {

	double q1, q2, q3;
	dd_t r;

	q1 = a.hi / b.hi;
	r = dd_sub_dd(a, dd_mul_d(b, q1));

	q2 = r.hi / b.hi;
	r = dd_sub_dd(r, dd_mul_d(b, q2));

	q3 = r.hi / b.hi;
	QUICK_TWO_SUM(q1, q2, q1, q2);
	return dd_add_d(dd_set_dd(q1, q2), q3);
}

static dd_complex_t cplx_set(dd_t re, dd_t im) {
	dd_complex_t c;
	c.r = re;
	c.i = im;
	return c;
}

static dd_complex_t cplx_set_d(double re, double im) {
	return cplx_set(dd_set_d(re), dd_set_d(im));
}

static dd_complex_t cplx_add(dd_complex_t a, dd_complex_t b) {
	return cplx_set(dd_add_dd(a.r, b.r),
			dd_add_dd(a.i, b.i));
}

static dd_complex_t cplx_sub(dd_complex_t a, dd_complex_t b) {
	return cplx_set(dd_sub_dd(a.r, b.r),
			dd_sub_dd(a.i, b.i));
}

static dd_complex_t cplx_mul(dd_complex_t a, dd_complex_t b) {
	return cplx_set(dd_sub_dd(dd_mul_dd(a.r, b.r),
				  dd_mul_dd(a.i, b.i)),
			dd_add_dd(dd_mul_dd(a.r, b.i),
				  dd_mul_dd(a.i, b.r)));
}

static dd_complex_t cplx_div(dd_complex_t a, dd_complex_t b) {

	dd_complex_t ans;
	if (dd_cmp_dd(dd_fabs(b.r), dd_fabs(b.i)) >= 0) {
		dd_t q = dd_div_dd(b.i, b.r);
		dd_t den = dd_add_dd(b.r, dd_mul_dd(q, b.i));
		ans = cplx_set(dd_div_dd(dd_add_dd(a.r, 
						    dd_mul_dd(q, a.i)),
					  den),
				dd_div_dd(dd_sub_dd(a.i,
						    dd_mul_dd(q, a.r)),
					  den));
	}
	else {
		dd_t q = dd_div_dd(b.r, b.i);
		dd_t den = dd_add_dd(b.i, dd_mul_dd(q, b.r));
		ans = cplx_set(dd_div_dd(dd_add_dd(dd_mul_dd(q, a.r),
						    a.i), 
					  den),
				dd_div_dd(dd_sub_dd(dd_mul_dd(q, a.i),
						    a.r),
					  den));
	}

	return ans;
}

static uint32_t mp_mod64(uint64_t a, uint32_t n) {

	uint32_t hi = (uint32_t)(a >> 32);
	uint32_t lo = (uint32_t)a;

	__asm__("divl %2"
	     : "+d"(hi), "+a"(lo)
	     : "rm"(n) : "cc");
	return hi;
}

static uint32_t mp_modmul_1(uint32_t a, uint32_t b, uint32_t n) {
	uint64_t acc = (uint64_t)a * (uint64_t)b;
	return mp_mod64(acc, n);
}

static uint32_t mp_modsub_1(uint32_t a, uint32_t b, uint32_t p) {

	uint32_t t = 0, tr;
	tr = a - b;
	if (tr > a)
		t = p;
	return tr + t;
}

static uint32_t mp_modadd_1(uint32_t a, uint32_t b, uint32_t p) {

	return mp_modsub_1(a, p - b, p);
}

static uint32_t mp_modinv_1(uint32_t a, uint32_t p) {

	uint32_t ps1, ps2, dividend, divisor, rem, q, t;
	uint32_t parity;

	q = 1; rem = a; dividend = p; divisor = a;
	ps1 = 1; ps2 = 0; parity = 0;

	while (divisor > 1) {
		rem = dividend - divisor;
		t = rem - divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t;
		if (rem >= divisor) {
			q = dividend / divisor;
			rem = dividend % divisor;
			q *= ps1;
		} } } } } } } } }

		q += ps2;
		parity = ~parity;
		dividend = divisor;
		divisor = rem;
		ps2 = ps1;
		ps1 = q;
	}
	
	if (parity == 0)
		return ps1;
	else
		return p - ps1;
}

static uint32_t mp_expo_1(uint32_t a, uint32_t b, uint32_t n) {

	uint32_t res = 1;
	while (b) {
		if (b & 1)
			res = mp_modmul_1(res, a, n);
		a = mp_modmul_1(a, a, n);
		b = b >> 1;
	}
	return res;
}

/*------------------------------------------------------------------*/
static uint64_t sqr_mac(uint32_t a, uint32_t b, uint64_t c, 
				uint32_t q, uint32_t _mod, uint64_t psq) {

	/* For 32-bit a, b, q, _mod and 64-bit c, psq compute

	   (a * b + c - q * _mod) % psq

	   psq is the square of a prime, c is less than psq 
	   and the 32-bit quantities are all assumed less 
	   than the prime */
	
	uint64_t ans;
	uint64_t tmp;
	ans = (uint64_t)a * (uint64_t)b;
	tmp = psq - c;
	ans = ans - tmp + (tmp > ans ? psq : 0);
	tmp = (uint64_t)q * (uint64_t)_mod;
	ans = ans - tmp + (tmp > ans ? psq : 0);

	return ans;
}

/*------------------------------------------------------------------*/
static uint64_t sqr_mac0(uint32_t a, uint32_t b,
			uint32_t q, uint32_t _mod, uint64_t psq) {

	/* For 32-bit a, b, q, _mod, compute

	   (a * b - q * _mod) % psq

	   psq is the square of a prime and the 32-bit quantities 
	   are all assumed less than the prime */

	uint64_t ans;
	uint64_t tmp;
	ans = (uint64_t)a * (uint64_t)b;
	tmp = (uint64_t)q * (uint64_t)_mod;
	ans = ans - tmp + (tmp > ans ? psq : 0);

	return ans;
}

/*------------------------------------------------------------------*/
static uint32_t mul_mac(uint32_t a, uint32_t b, uint32_t c, 
				uint32_t q, uint32_t _mod, 
				uint32_t p, uint64_t psq) {

	/* For 32-bit a, b, q, c, _mod and 64-bit psq compute

	   (a * b + c - q * _mod) % p

	   psq is the square of p, and the 32-bit quantities 
	   are all assumed less than p
	
	   Note that this is intended for a high-performance
	   processor with conditional move instructions. If the compiler 
	   does not inline, or use these instructions, or realize that
	   many inputs are cached from previous calls then the overall
	   performance will be dismal. */

	uint64_t ans, tmp;

	/* for a,b,c < p the following will never exceed psq,
	   so no correction is needed */
	ans = (uint64_t)a * (uint64_t)b + (uint64_t)c;

	tmp = (uint64_t)q * (uint64_t)_mod;
	ans = ans - tmp + (tmp > ans ? psq : 0);
	return mp_mod64(ans, p);
}


static void gmp2mp(mpz_t src, mp_t *dest) {

	size_t count;

	mp_clear(dest);
	mpz_export(dest->val, &count, -1, sizeof(uint32_t),
			0, (size_t)0, src);
	dest->nwords = count;
}

static dd_t dd_gmp2dd(mpz_t x) {

	signed_mp_t mpx;

	gmp2mp(x, &mpx.num);

	mpx.sign = POSITIVE;
	if (mpz_sgn(x) < 0)
		mpx.sign = NEGATIVE;

	return dd_signed_mp2dd(&mpx);
}

static dd_precision_t dd_set_precision_ieee(void) {
	dd_precision_t old_prec, new_prec;
	__asm__ volatile ("fnstcw %0":"=m"(old_prec));
	new_prec = (old_prec & ~0x0300) | 0x0200;
	__asm__ volatile ("fldcw %0": :"m"(new_prec));
	return old_prec;
}

static void dd_clear_precision(dd_precision_t old_prec) {
	__asm__ volatile ("fldcw %0": :"m"(old_prec));
}

static uint32_t dd_precision_is_ieee(void) {
	dd_precision_t prec;
	__asm__ volatile ("fnstcw %0":"=m"(prec));
	return ((prec & ~0x0300) == 0x0200) ? 1 : 0;
}

static void mpz_poly_init(mpz_poly_t * poly) {
	uint32_t i;

	memset(poly, 0, sizeof(mpz_poly_t));

	mpz_init(poly->tmp1);
	mpz_init(poly->tmp2);
	mpz_init(poly->tmp3);
	for (i = 0; i <= MAX_POLY_DEGREE; i++)
		mpz_init_set_ui(poly->coeff[i], 0);
}

static void mpz_poly_free(mpz_poly_t * poly) {
	uint32_t i;

	mpz_clear(poly->tmp1);
	mpz_clear(poly->tmp2);
	mpz_clear(poly->tmp3);
	for (i = 0; i <= MAX_POLY_DEGREE; i++)
		mpz_clear(poly->coeff[i]);
}

void analyze_one_poly_hook( long, mpz_t*, mpz_t*, double, double*, double*, double*, unsigned int* );
void analyze_poly(poly_config_t *, poly_select_t *);
uint32_t analyze_poly_roots(mpz_poly_t *, uint32_t, double *);
uint32_t analyze_poly_size(integrate_t *, ddpoly_t *, ddpoly_t *, double *);
static double eval_dpoly(dpoly_t *, double, double, double, double);
static void get_bernstein_combined_score(poly_select_t *, double);
uint32_t analyze_poly_murphy(integrate_t *, dickman_t *, ddpoly_t *, double,
			ddpoly_t *, double, double, double *, uint32_t *);
static double murphy_integrand(double, double, void *);
void poly_config_init(poly_config_t *);
void poly_config_free(poly_config_t *);
void poly_select_init(poly_select_t *);
void poly_select_free(poly_select_t *);
uint32_t poly_get_zeros(uint32_t *, mpz_poly_t *, uint32_t, uint32_t *, uint32_t);
static void poly_gcd(poly_t g_in, poly_t h_in, uint32_t p);
static void get_zeros_rec(uint32_t *, uint32_t, uint32_t *, poly_t, uint32_t);
static void poly_xpow(poly_t, uint32_t, uint32_t, poly_t, uint32_t);
static void poly_expo_modmul(uint32_t *, uint32_t, uint32_t, uint32_t, uint64_t);
static void poly_expo_square(uint32_t *, uint32_t, uint32_t, uint64_t);
uint32_t poly_get_zeros_and_mult(uint32_t *, uint32_t *, mpz_poly_t *, uint32_t, uint32_t *);
static void poly_reduce_mod_p(poly_t, mpz_poly_t *, uint32_t);
static void poly_fix_degree(poly_t);
static void poly_mod(poly_t, poly_t, poly_t, uint32_t);
static void poly_cp(poly_t, poly_t);
static void poly_make_monic(poly_t, poly_t, uint32_t);
static uint32_t jenkins_traub(complex_t [], uint32_t, complex_t []);
static double get_root_freq(mpz_poly_t *, uint32_t, uint32_t);
uint32_t find_poly_roots(dd_t *, uint32_t, dd_complex_t *);
static uint32_t polish_root(dd_complex_t *, uint32_t, dd_complex_t, dd_complex_t *, double);
static complex_t complex(double, double);
static int find_one_root(complex_t *, jt_t *);
static void stage1(uint32_t, complex_t [], complex_t []);
static uint32_t stage2(uint32_t, complex_t *, jt_t *);
static complex_t cadd(complex_t, complex_t);
static complex_t cneg(complex_t);
static complex_t cscale(double, complex_t);
static complex_t cmul(complex_t, complex_t);
static complex_t cmac(complex_t, complex_t, complex_t);
static double cmod(complex_t);
static complex_t cdiv(complex_t, complex_t);
static double cauchy_bound(uint32_t, const complex_t []);
static complex_t poly_val(uint32_t, complex_t, const complex_t [], complex_t []);
static void next_hpoly(complex_t, jt_t *);
static complex_t next_correction(complex_t, complex_t, jt_t *);
static uint32_t stage3(complex_t *, jt_t *);
uint32_t mp_modsqrt_1(uint32_t, uint32_t);
int32_t mp_legendre_1(uint32_t, uint32_t);
static int compare_double(const void *, const void *);
static double get_polyval(ddpoly_t *, double, double);
static double integrand(double, double, void *);
static void de_run_core(integrate_t *, integrand_t, void *, subinterval_t *);
static void heapify(subinterval_t *, uint32_t, uint32_t);
static void heap_insert(subinterval_t *, uint32_t);
static void make_heap(subinterval_t *, uint32_t);
static uint32_t de_run(integrate_t *, integrand_t, void *, double *, uint32_t);
void integrate_init(integrate_t *, double, enum integrator_type);
void integrate_free(integrate_t *);
uint32_t integrate_run(integrate_t *, integrand_t, void *, double *, uint32_t);
static void de_init(integrate_t *, double);
static void de_free(integrate_t *);
static void de_fill_sample_group(de_point_t *, uint32_t, double, double, double, double);
void dickman_init(dickman_t *);
void dickman_free(dickman_t *);
double dickman(dickman_t *, double);

static const complex_t czero = { 0, 0 };
static const double mre = 2.0 * M_SQRT2 * DBL_EPSILON;

static const double recip_factorial[] = {
  5.00000000000000000e-01,
  1.66666666666666657e-01,
  4.16666666666666644e-02,
  8.33333333333333322e-03,
  1.38888888888888894e-03,
  1.98412698412698413e-04,
  2.48015873015873016e-05,
  2.75573192239858925e-06,
  2.75573192239858883e-07,
  2.50521083854417202e-08,
  2.08767569878681002e-09,
  1.60590438368216133e-10,
  1.14707455977297245e-11,
  7.64716373181981641e-13,
  4.77947733238738525e-14,
  2.81145725434552060e-15,
};

const uint8_t prime_delta[PRECOMPUTED_NUM_PRIMES] = {
 2, 1, 2, 2, 4, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8,
 4, 2, 4, 2, 4,14, 4, 6, 2,10, 2, 6, 6, 4, 6, 6, 2,10, 2, 4, 2,12,12, 4, 2,
 4, 6, 2,10, 6, 6, 6, 2, 6, 4, 2,10,14, 4, 2, 4,14, 6,10, 2, 4, 6, 8, 6, 6,
 4, 6, 8, 4, 8,10, 2,10, 2, 6, 4, 6, 8, 4, 2, 4,12, 8, 4, 8, 4, 6,12, 2,18,
 6,10, 6, 6, 2, 6,10, 6, 6, 2, 6, 6, 4, 2,12,10, 2, 4, 6, 6, 2,12, 4, 6, 8,
10, 8,10, 8, 6, 6, 4, 8, 6, 4, 8, 4,14,10,12, 2,10, 2, 4, 2,10,14, 4, 2, 4,
14, 4, 2, 4,20, 4, 8,10, 8, 4, 6, 6,14, 4, 6, 6, 8, 6,12, 4, 6, 2,10, 2, 6,
10, 2,10, 2, 6,18, 4, 2, 4, 6, 6, 8, 6, 6,22, 2,10, 8,10, 6, 6, 8,12, 4, 6,
 6, 2, 6,12,10,18, 2, 4, 6, 2, 6, 4, 2, 4,12, 2, 6,34, 6, 6, 8,18,10,14, 4,
 2, 4, 6, 8, 4, 2, 6,12,10, 2, 4, 2, 4, 6,12,12, 8,12, 6, 4, 6, 8, 4, 8, 4,
14, 4, 6, 2, 4, 6, 2, 6,10,20, 6, 4, 2,24, 4, 2,10,12, 2,10, 8, 6, 6, 6,18,
 6, 4, 2,12,10,12, 8,16,14, 6, 4, 2, 4, 2,10,12, 6, 6,18, 2,16, 2,22, 6, 8,
 6, 4, 2, 4, 8, 6,10, 2,10,14,10, 6,12, 2, 4, 2,10,12, 2,16, 2, 6, 4, 2,10,
 8,18,24, 4, 6, 8,16, 2, 4, 8,16, 2, 4, 8, 6, 6, 4,12, 2,22, 6, 2, 6, 4, 6,
14, 6, 4, 2, 6, 4, 6,12, 6, 6,14, 4, 6,12, 8, 6, 4,26,18,10, 8, 4, 6, 2, 6,
22,12, 2,16, 8, 4,12,14,10, 2, 4, 8, 6, 6, 4, 2, 4, 6, 8, 4, 2, 6,10, 2,10,
 8, 4,14,10,12, 2, 6, 4, 2,16,14, 4, 6, 8, 6, 4,18, 8,10, 6, 6, 8,10,12,14,
 4, 6, 6, 2,28, 2,10, 8, 4,14, 4, 8,12, 6,12, 4, 6,20,10, 2,16,26, 4, 2,12,
 6, 4,12, 6, 8, 4, 8,22, 2, 4, 2,12,28, 2, 6, 6, 6, 4, 6, 2,12, 4,12, 2,10,
 2,16, 2,16, 6,20,16, 8, 4, 2, 4, 2,22, 8,12, 6,10, 2, 4, 6, 2, 6,10, 2,12,
10, 2,10,14, 6, 4, 6, 8, 6, 6,16,12, 2, 4,14, 6, 4, 8,10, 8, 6, 6,22, 6, 2,
10,14, 4, 6,18, 2,10,14, 4, 2,10,14, 4, 8,18, 4, 6, 2, 4, 6, 2,12, 4,20,22,
12, 2, 4, 6, 6, 2, 6,22, 2, 6,16, 6,12, 2, 6,12,16, 2, 4, 6,14, 4, 2,18,24,
10, 6, 2,10, 2,10, 2,10, 6, 2,10, 2,10, 6, 8,30,10, 2,10, 8, 6,10,18, 6,12,
12, 2,18, 6, 4, 6, 6,18, 2,10,14, 6, 4, 2, 4,24, 2,12, 6,16, 8, 6, 6,18,16,
 2, 4, 6, 2, 6, 6,10, 6,12,12,18, 2, 6, 4,18, 8,24, 4, 2, 4, 6, 2,12, 4,14,
30,10, 6,12,14, 6,10,12, 2, 4, 6, 8, 6,10, 2, 4,14, 6, 6, 4, 6, 2,10, 2,16,
12, 8,18, 4, 6,12, 2, 6, 6, 6,28, 6,14, 4, 8,10, 8,12,18, 4, 2, 4,24,12, 6,
 2,16, 6, 6,14,10,14, 4,30, 6, 6, 6, 8, 6, 4, 2,12, 6, 4, 2, 6,22, 6, 2, 4,
18, 2, 4,12, 2, 6, 4,26, 6, 6, 4, 8,10,32,16, 2, 6, 4, 2, 4, 2,10,14, 6, 4,
 8,10, 6,20, 4, 2, 6,30, 4, 8,10, 6, 6, 8, 6,12, 4, 6, 2, 6, 4, 6, 2,10, 2,
16, 6,20, 4,12,14,28, 6,20, 4,18, 8, 6, 4, 6,14, 6, 6,10, 2,10,12, 8,10, 2,
10, 8,12,10,24, 2, 4, 8, 6, 4, 8,18,10, 6, 6, 2, 6,10,12, 2,10, 6, 6, 6, 8,
 6,10, 6, 2, 6, 6, 6,10, 8,24, 6,22, 2,18, 4, 8,10,30, 8,18, 4, 2,10, 6, 2,
 6, 4,18, 8,12,18,16, 6, 2,12, 6,10, 2,10, 2, 6,10,14, 4,24, 2,16, 2,10, 2,
10,20, 4, 2, 4, 8,16, 6, 6, 2,12,16, 8, 4, 6,30, 2,10, 2, 6, 4, 6, 6, 8, 6,
 4,12, 6, 8,12, 4,14,12,10,24, 6,12, 6, 2,22, 8,18,10, 6,14, 4, 2, 6,10, 8,
 6, 4, 6,30,14,10, 2,12,10, 2,16, 2,18,24,18, 6,16,18, 6, 2,18, 4, 6, 2,10,
 8,10, 6, 6, 8, 4, 6, 2,10, 2,12, 4, 6, 6, 2,12, 4,14,18, 4, 6,20, 4, 8, 6,
 4, 8, 4,14, 6, 4,14,12, 4, 2,30, 4,24, 6, 6,12,12,14, 6, 4, 2, 4,18, 6,12,
 8, 6, 4,12, 2,12,30,16, 2, 6,22,14, 6,10,12, 6, 2, 4, 8,10, 6, 6,24,14, 6,
 4, 8,12,18,10, 2,10, 2, 4, 6,20, 6, 4,14, 4, 2, 4,14, 6,12,24,10, 6, 8,10,
 2,30, 4, 6, 2,12, 4,14, 6,34,12, 8, 6,10, 2, 4,20,10, 8,16, 2,10,14, 4, 2,
12, 6,16, 6, 8, 4, 8, 4, 6, 8, 6, 6,12, 6, 4, 6, 6, 8,18, 4,20, 4,12, 2,10,
 6, 2,10,12, 2, 4,20, 6,30, 6, 4, 8,10,12, 6, 2,28, 2, 6, 4, 2,16,12, 2, 6,
10, 8,24,12, 6,18, 6, 4,14, 6, 4,12, 8, 6,12, 4, 6,12, 6,12, 2,16,20, 4, 2,
10,18, 8, 4,14, 4, 2, 6,22, 6,14, 6, 6,10, 6, 2,10, 2, 4, 2,22, 2, 4, 6, 6,
12, 6,14,10,12, 6, 8, 4,36,14,12, 6, 4, 6, 2,12, 6,12,16, 2,10, 8,22, 2,12,
 6, 4, 6,18, 2,12, 6, 4,12, 8, 6,12, 4, 6,12, 6, 2,12,12, 4,14, 6,16, 6, 2,
10, 8,18, 6,34, 2,28, 2,22, 6, 2,10,12, 2, 6, 4, 8,22, 6, 2,10, 8, 4, 6, 8,
 4,12,18,12,20, 4, 6, 6, 8, 4, 2,16,12, 2,10, 8,10, 2, 4, 6,14,12,22, 8,28,
 2, 4,20, 4, 2, 4,14,10,12, 2,12,16, 2,28, 8,22, 8, 4, 6, 6,14, 4, 8,12, 6,
 6, 4,20, 4,18, 2,12, 6, 4, 6,14,18,10, 8,10,32, 6,10, 6, 6, 2, 6,16, 6, 2,
12, 6,28, 2,10, 8,16, 6, 8, 6,10,24,20,10, 2,10, 2,12, 4, 6,20, 4, 2,12,18,
10, 2,10, 2, 4,20,16,26, 4, 8, 6, 4,12, 6, 8,12,12, 6, 4, 8,22, 2,16,14,10,
 6,12,12,14, 6, 4,20, 4,12, 6, 2, 6, 6,16, 8,22, 2,28, 8, 6, 4,20, 4,12,24,
20, 4, 8,10, 2,16, 2,12,12,34, 2, 4, 6,12, 6, 6, 8, 6, 4, 2, 6,24, 4,20,10,
 6, 6,14, 4, 6, 6, 2,12, 6,10, 2,10, 6,20, 4,26, 4, 2, 6,22, 2,24, 4, 6, 2,
 4, 6,24, 6, 8, 4, 2,34, 6, 8,16,12, 2,10, 2,10, 6, 8, 4, 8,12,22, 6,14, 4,
26, 4, 2,12,10, 8, 4, 8,12, 4,14, 6,16, 6, 8, 4, 6, 6, 8, 6,10,12, 2, 6, 6,
16, 8, 6, 6,12,10, 2, 6,18, 4, 6, 6, 6,12,18, 8, 6,10, 8,18, 4,14, 6,18,10,
 8,10,12, 2, 6,12,12,36, 4, 6, 8, 4, 6, 2, 4,18,12, 6, 8, 6, 6, 4,18, 2, 4,
 2,24, 4, 6, 6,14,30, 6, 4, 6,12, 6,20, 4, 8, 4, 8, 6, 6, 4,30, 2,10,12, 8,
10, 8,24, 6,12, 4,14, 4, 6, 2,28,14,16, 2,12, 6, 4,20,10, 6, 6, 6, 8,10,12,
14,10,14,16,14,10,14, 6,16, 6, 8, 6,16,20,10, 2, 6, 4, 2, 4,12, 2,10, 2, 6,
22, 6, 2, 4,18, 8,10, 8,22, 2,10,18,14, 4, 2, 4,18, 2, 4, 6, 8,10, 2,30, 4,
30, 2,10, 2,18, 4,18, 6,14,10, 2, 4,20,36, 6, 4, 6,14, 4,20,10,14,22, 6, 2,
30,12,10,18, 2, 4,14, 6,22,18, 2,12, 6, 4, 8, 4, 8, 6,10, 2,12,18,10,14,16,
14, 4, 6, 6, 2, 6, 4, 2,28, 2,28, 6, 2, 4, 6,14, 4,12,14,16,14, 4, 6, 8, 6,
 4, 6, 6, 6, 8, 4, 8, 4,14,16, 8, 6, 4,12, 8,16, 2,10, 8, 4, 6,26, 6,10, 8,
 4, 6,12,14,30, 4,14,22, 8,12, 4, 6, 8,10, 6,14,10, 6, 2,10,12,12,14, 6, 6,
18,10, 6, 8,18, 4, 6, 2, 6,10, 2,10, 8, 6, 6,10, 2,18,10, 2,12, 4, 6, 8,10,
12,14,12, 4, 8,10, 6, 6,20, 4,14,16,14,10, 8,10,12, 2,18, 6,12,10,12, 2, 4,
 2,12, 6, 4, 8, 4,44, 4, 2, 4, 2,10,12, 6, 6,14, 4, 6, 6, 6, 8, 6,36,18, 4,
 6, 2,12, 6, 6, 6, 4,14,22,12, 2,18,10, 6,26,24, 4, 2, 4, 2, 4,14, 4, 6, 6,
 8,16,12, 2,42, 4, 2, 4,24, 6, 6, 2,18, 4,14, 6,28,18,14, 6,10,12, 2, 6,12,
30, 6, 4, 6, 6,14, 4, 2,24, 4, 6, 6,26,10,18, 6, 8, 6, 6,30, 4,12,12, 2,16,
 2, 6, 4,12,18, 2, 6, 4,26,12, 6,12, 4,24,24,12, 6, 2,12,28, 8, 4, 6,12, 2,
18, 6, 4, 6, 6,20,16, 2, 6, 6,18,10, 6, 2, 4, 8, 6, 6,24,16, 6, 8,10, 6,14,
22, 8,16, 6, 2,12, 4, 2,22, 8,18,34, 2, 6,18, 4, 6, 6, 8,10, 8,18, 6, 4, 2,
 4, 8,16, 2,12,12, 6,18, 4, 6, 6, 6, 2, 6,12,10,20,12,18, 4, 6, 2,16, 2,10,
14, 4,30, 2,10,12, 2,24, 6,16, 8,10, 2,12,22, 6, 2,16,20,10, 2,12,12,18,10,
12, 6, 2,10, 2, 6,10,18, 2,12, 6, 4, 6, 2,24,28, 2, 4, 2,10, 2,16,12, 8,22,
 2, 6, 4, 2,10, 6,20,12,10, 8,12, 6, 6, 6, 4,18, 2, 4,12,18, 2,12, 6, 4, 2,
16,12,12,14, 4, 8,18, 4,12,14, 6, 6, 4, 8, 6, 4,20,12,10,14, 4, 2,16, 2,12,
30, 4, 6,24,20,24,10, 8,12,10,12, 6,12,12, 6, 8,16,14, 6, 4, 6,36,20,10,30,
12, 2, 4, 2,28,12,14, 6,22, 8, 4,18, 6,14,18, 4, 6, 2, 6,34,18, 2,16, 6,18,
 2,24, 4, 2, 6,12, 6,12,10, 8, 6,16,12, 8,10,14,40, 6, 2, 6, 4,12,14, 4, 2,
 4, 2, 4, 8, 6,10, 6, 6, 2, 6, 6, 6,12, 6,24,10, 2,10, 6,12, 6, 6,14, 6, 6,
52,20, 6,10, 2,10, 8,10,12,12, 2, 6, 4,14,16, 8,12, 6,22, 2,10, 8, 6,22, 2,
22, 6, 8,10,12,12, 2,10, 6,12, 2, 4,14,10, 2, 6,18, 4,12, 8,18,12, 6, 6, 4,
 6, 6,14, 4, 2,12,12, 4, 6,18,18,12, 2,16,12, 8,18,10,26, 4, 6, 8, 6, 6, 4,
 2,10,20, 4, 6, 8, 4,20,10, 2,34, 2, 4,24, 2,12,12,10, 6, 2,12,30, 6,12,16,
12, 2,22,18,12,14,10, 2,12,12, 4, 2, 4, 6,12, 2,16,18, 2,40, 8,16, 6, 8,10,
 2, 4,18, 8,10, 8,12, 4,18, 2,18,10, 2, 4, 2, 4, 8,28, 2, 6,22,12, 6,14,18,
 4, 6, 8, 6, 6,10, 8, 4, 2,18,10, 6,20,22, 8, 6,30, 4, 2, 4,18, 6,30, 2, 4,
 8, 6, 4, 6,12,14,34,14, 6, 4, 2, 6, 4,14, 4, 2, 6,28, 2, 4, 6, 8,10, 2,10,
 2,10, 2, 4,30, 2,12,12,10,18,12,14,10, 2,12, 6,10, 6,14,12, 4,14, 4,18, 2,
10, 8, 4, 8,10,12,18,18, 8, 6,18,16,14, 6, 6,10,14, 4, 6, 2,12,12, 4, 6, 6,
12, 2,16, 2,12, 6, 4,14, 6, 4, 2,12,18, 4,36,18,12,12, 2, 4, 2, 4, 8,12, 4,
36, 6,18, 2,12,10, 6,12,24, 8, 6, 6,16,12, 2,18,10,20,10, 2, 6,18, 4, 2,40,
 6, 2,16, 2, 4, 8,18,10,12, 6, 2,10, 8, 4, 6,12, 2,10,18, 8, 6, 4,20, 4, 6,
36, 6, 2,10, 6,24, 6,14,16, 6,18, 2,10,20,10, 8, 6, 4, 6, 2,10, 2,12, 4, 2,
 4, 8,10, 6,12,18,14,12,16, 8, 6,16, 8, 4, 2, 6,18,24,18,10,12, 2, 4,14,10,
 6, 6, 6,18,12, 2,28,18,14,16,12,14,24,12,22, 6, 2,10, 8, 4, 2, 4,14,12, 6,
 4, 6,14, 4, 2, 4,30, 6, 2, 6,10, 2,30,22, 2, 4, 6, 8, 6, 6,16,12,12, 6, 8,
 4, 2,24,12, 4, 6, 8, 6, 6,10, 2, 6,12,28,14, 6, 4,12, 8, 6,12, 4, 6,14, 6,
12,10, 6, 6, 8, 6, 6, 4, 2, 4, 8,12, 4,14,18,10, 2,16, 6,20, 6,10, 8, 4,30,
36,12, 8,22,12, 2, 6,12,16, 6, 6, 2,18, 4,26, 4, 8,18,10, 8,10, 6,14, 4,20,
22,18,12, 8,28,12, 6, 6, 8, 6,12,24,16,14, 4,14,12, 6,10,12,20, 6, 4, 8,18,
12,18,10, 2, 4,20,10,14, 4, 6, 2,10,24,18, 2, 4,20,16,14,10,14, 6, 4, 6,20,
 6,10, 6, 2,12, 6,30,10, 8, 6, 4, 6, 8,40, 2, 4, 2,12,18, 4, 6, 8,10, 6,18,
18, 2,12,16, 8, 6, 4, 6, 6, 2,52,14, 4,20,16, 2, 4, 6,12, 2, 6,12,12, 6, 4,
14,10, 6, 6,14,10,14,16, 8, 6,12, 4, 8,22, 6, 2,18,22, 6, 2,18, 6,16,14,10,
 6,12, 2, 6, 4, 8,18,12,16, 2, 4,14, 4, 8,12,12,30,16, 8, 4, 2, 6,22,12, 8,
10, 6, 6, 6,14, 6,18,10,12, 2,10, 2, 4,26, 4,12, 8, 4,18, 8,10,14,16, 6, 6,
 8,10, 6, 8, 6,12,10,20,10, 8, 4,12,26,18, 4,12,18, 6,30, 6, 8, 6,22,12, 2,
 4, 6, 6, 2,10, 2, 4, 6, 6, 2, 6,22,18, 6,18,12, 8,12, 6,10,12, 2,16, 2,10,
 2,10,18, 6,20, 4, 2, 6,22, 6, 6,18, 6,14,12,16, 2, 6, 6, 4,14,12, 4, 2,18,
16,36,12, 6,14,28, 2,12, 6,12, 6, 4, 2,16,30, 8,24, 6,30,10, 2,18, 4, 6,12,
 8,22, 2, 6,22,18, 2,10, 2,10,30, 2,28, 6,14,16, 6,20,16, 2, 6, 4,32, 4, 2,
 4, 6, 2,12, 4, 6, 6,12, 2, 6, 4, 6, 8, 6, 4,20, 4,32,10, 8,16, 2,22, 2, 4,
 6, 8, 6,16,14, 4,18, 8, 4,20, 6,12,12, 6,10, 2,10, 2,12,28,12,18, 2,18,10,
 8,10,48, 2, 4, 6, 8,10, 2,10,30, 2,36, 6,10, 6, 2,18, 4, 6, 8,16,14,16, 6,
14, 4,20, 4, 6, 2,10,12, 2, 6,12, 6, 6, 4,12, 2, 6, 4,12, 6, 8, 4, 2, 6,18,
10, 6, 8,12, 6,22, 2, 6,12,18, 4,14, 6, 4,20, 6,16, 8, 4, 8,22, 8,12, 6, 6,
16,12,18,30, 8, 4, 2, 4, 6,26, 4,14,24,22, 6, 2, 6,10, 6,14, 6, 6,12,10, 6,
 2,12,10,12, 8,18,18,10, 6, 8,16, 6, 6, 8,16,20, 4, 2,10, 2,10,12, 6, 8, 6,
10,20,10,18,26, 4, 6,30, 2, 4, 8, 6,12,12,18, 4, 8,22, 6, 2,12,34, 6,18,12,
 6, 2,28,14,16,14, 4,14,12, 4, 6, 6, 2,36, 4, 6,20,12,24, 6,22, 2,16,18,12,
12,18, 2, 6, 6, 6, 4, 6,14, 4, 2,22, 8,12, 6,10, 6, 8,12,18,12, 6,10, 2,22,
14, 6, 6, 4,18, 6,20,22, 2,12,24, 4,18,18, 2,22, 2, 4,12, 8,12,10,14, 4, 2,
18,16,38, 6, 6, 6,12,10, 6,12, 8, 6, 4, 6,14,30, 6,10, 8,22, 6, 8,12,10, 2,
10, 2, 6,10, 2,10,12,18,20, 6, 4, 8,22, 6, 6,30, 6,14, 6,12,12, 6,10, 2,10,
30, 2,16, 8, 4, 2, 6,18, 4, 2, 6, 4,26, 4, 8, 6,10, 2, 4, 6, 8, 4, 6,30,12,
 2, 6, 6, 4,20,22, 8, 4, 2, 4,72, 8, 4, 8,22, 2, 4,14,10, 2, 4,20, 6,10,18,
 6,20,16, 6, 8, 6, 4,20,12,22, 2, 4, 2,12,10,18, 2,22, 6,18,30, 2,10,14,10,
 8,16,50, 6,10, 8,10,12, 6,18, 2,22, 6, 2, 4, 6, 8, 6, 6,10,18, 2,22, 2,16,
14,10, 6, 2,12,10,20, 4,14, 6, 4,36, 2, 4, 6,12, 2, 4,14,12, 6, 4, 6, 2, 6,
 4,20,10, 2,10, 6,12, 2,24,12,12, 6, 6, 4,24, 2, 4,24, 2, 6, 4, 6, 8,16, 6,
 2,10,12,14, 6,34, 6,14, 6, 4, 2,30,22, 8, 4, 6, 8, 4, 2,28, 2, 6, 4,26,18,
22, 2, 6,16, 6, 2,16,12, 2,12, 4, 6, 6,14,10, 6, 8,12, 4,18, 2,10, 8,16, 6,
 6,30, 2,10,18, 2,10, 8, 4, 8,12,24,40, 2,12,10, 6,12, 2,12, 4, 2, 4, 6,18,
14,12, 6, 4,14,30, 4, 8,10, 8, 6,10,18, 8, 4,14,16, 6, 8, 4, 6, 2,10, 2,12,
 4, 2, 4, 6, 8, 4, 6,32,24,10, 8,18,10, 2, 6,10, 2, 4,18, 6,12, 2,16, 2,22,
 6, 6, 8,18, 4,18,12, 8, 6, 4,20, 6,30,22,12, 2, 6,18, 4,62, 4, 2,12, 6,10,
 2,12,12,28, 2, 4,14,22, 6, 2, 6, 6,10,14, 4, 2,10, 6, 8,10,14,10, 6, 2,12,
22,18, 8,10,18,12, 2,12, 4,12, 2,10, 2, 6,18, 6, 6,34, 6, 2,12, 4, 6,18,18,
 2,16, 6, 6, 8, 6,10,18, 8,10, 8,10, 2, 4,18,26,12,22, 2, 4, 2,22, 6, 6,14,
16, 6,20,10,12, 2,18,42, 4,24, 2, 6,10,12, 2, 6,10, 8, 4, 6,12,12, 8, 4, 6,
12,30,20, 6,24, 6,10,12, 2,10,20, 6, 6, 4,12,14,10,18,12, 8, 6,12, 4,14,10,
 2,12,30,16, 2,12, 6, 4, 2, 4, 6,26, 4,18, 2, 4, 6,14,54, 6,52, 2,16, 6, 6,
12,26, 4, 2, 6,22, 6, 2,12,12, 6,10,18, 2,12,12,10,18,12, 6, 8, 6,10, 6, 8,
 4, 2, 4,20,24, 6, 6,10,14,10, 2,22, 6,14,10,26, 4,18, 8,12,12,10,12, 6, 8,
16, 6, 8, 6, 6,22, 2,10,20,10, 6,44,18, 6,10, 2, 4, 6,14, 4,26, 4, 2,12,10,
 8, 4, 8,12, 4,12, 8,22, 8, 6,10,18, 6, 6, 8, 6,12, 4, 8,18,10,12, 6,12, 2,
 6, 4, 2,16,12,12,14,10,14, 6,10,12, 2,12, 6, 4, 6, 2,12, 4,26, 6,18, 6,10,
 6, 2,18,10, 8, 4,26,10,20, 6,16,20,12,10, 8,10, 2,16, 6,20,10,20, 4,30, 2,
 4, 8,16, 2,18, 4, 2, 6,10,18,12,14,18, 6,16,20, 6, 4, 8, 6, 4, 6,12, 8,10,
 2,12, 6, 4, 2, 6,10, 2,16,12,14,10, 6, 8, 6,28, 2, 6,18,30,34, 2,16,12, 2,
18,16, 6, 8,10, 8,10, 8,10,44, 6, 6, 4,20, 4, 2, 4,14,28, 8, 6,16,14,30, 6,
30, 4,14,10, 6, 6, 8, 4,18,12, 6, 2,22,12, 8, 6,12, 4,14, 4, 6, 2, 4,18,20,
 6,16,38,16, 2, 4, 6, 2,40,42,14, 4, 6, 2,24,10, 6, 2,18,10,12, 2,16, 2, 6,
16, 6, 8, 4, 2,10, 6, 8,10, 2,18,16, 8,12,18,12, 6,12,10, 6, 6,18,12,14, 4,
 2,10,20, 6,12, 6,16,26, 4,18, 2, 4,32,10, 8, 6, 4, 6, 6,14, 6,18, 4, 2,18,
10, 8,10, 8,10, 2, 4, 6, 2,10,42, 8,12, 4, 6,18, 2,16, 8, 4, 2,10,14,12,10,
20, 4, 8,10,38, 4, 6, 2,10,20,10,12, 6,12,26,12, 4, 8,28, 8, 4, 8,24, 6,10,
 8, 6,16,12, 8,10,12, 8,22, 6, 2,10, 2, 6,10, 6, 6, 8, 6, 4,14,28, 8,16,18,
 8, 4, 6,20, 4,18, 6, 2,24,24, 6, 6,12,12, 4, 2,22, 2,10, 6, 8,12, 4,20,18,
 6, 4,12,24, 6, 6,54, 8, 6, 4,26,36, 4, 2, 4,26,12,12, 4, 6, 6, 8,12,10, 2,
12,16,18, 6, 8, 6,12,18,10, 2,54, 4, 2,10,30,12, 8, 4, 8,16,14,12, 6, 4, 6,
12, 6, 2, 4,14,12, 4,14, 6,24, 6, 6,10,12,12,20,18, 6, 6,16, 8, 4, 6,20, 4,
32, 4,14,10, 2, 6,12,16, 2, 4, 6,12, 2,10, 8, 6, 4, 2,10,14, 6, 6,12,18,34,
 8,10, 6,24, 6, 2,10,12, 2,30,10,14,12,12,16, 6, 6, 2,18, 4, 6,30,14, 4, 6,
 6, 2, 6, 4, 6,14, 6, 4, 8,10,12, 6,32,10, 8,22, 2,10, 6,24, 8, 4,30, 6, 2,
12,16, 8, 6, 4, 6, 8,16,14, 6, 6, 4, 2,10,12, 2,16,14, 4, 2, 4,20,18,10, 2,
10, 6,12,30, 8,18,12,10, 2, 6, 6, 4,12,12, 2, 4,12,18,24, 2,10, 6, 8,16, 8,
 6,12,10,14, 6,12, 6, 6, 4, 2,24, 4, 6, 8, 6, 4, 2, 4, 6,14, 4, 8,10,24,24,
12, 2, 6,12,22,30, 2, 6,18,10, 6, 6, 8, 4, 2, 6,10, 8,10, 6, 8,16, 6,14, 6,
 4,24, 8,10, 2,12, 6, 4,36, 2,22, 6, 8, 6,10, 8, 6,12,10,14,10, 6,18,12, 2,
12, 4,26,10,14,16,18, 8,18,12,12, 6,16,14,24,10,12, 8,22, 6, 2,10,60, 6, 2,
 4, 8,16,14,10, 6,24, 6,12,18,24, 2,30, 4, 2,12, 6,10, 2, 4,14, 6,16, 2,10,
 8,22,20, 6, 4,32, 6,18, 4, 2, 4, 2, 4, 8,52,14,22, 2,22,20,10, 8,10, 2, 6,
 4,14, 4, 6,20, 4, 6, 2,12,12, 6,12,16, 2,12,10, 8, 4, 6, 2,28,12, 8,10,12,
 2, 4,14,28, 8, 6, 4, 2, 4, 6, 2,12,58, 6,14,10, 2, 6,28,32, 4,30, 8, 6, 4,
 6,12,12, 2, 4, 6, 6,14,16, 8,30, 4, 2,10, 8, 6, 4, 6,26, 4,12, 2,10,18,12,
12,18, 2, 4,12, 8,12,10,20, 4, 8,16,12, 8, 6,16, 8,10,12,14, 6, 4, 8,12, 4,
20, 6,40, 8,16, 6,36, 2, 6, 4, 6, 2,22,18, 2,10, 6,36,14,12, 4,18, 8, 4,14,
10, 2,10, 8, 4, 2,18,16,12,14,10,14, 6, 6,42,10, 6, 6,20,10, 8,12, 4,12,18,
 2,10,14,18,10,18, 8, 6, 4,14, 6,10,30,14, 6, 6, 4,12,38, 4, 2, 4, 6, 8,12,
10, 6,18, 6,50, 6, 4, 6,12, 8,10,32, 6,22, 2,10,12,18, 2, 6, 4,30, 8, 6, 6,
18,10, 2, 4,12,20,10, 8,24,10, 2, 6,22, 6, 2,18,10,12, 2,30,18,12,28, 2, 6,
 4, 6,14, 6,12,10, 8, 4,12,26,10, 8, 6,16, 2,10,18,14, 6, 4, 6,14,16, 2, 6,
 4,12,20, 4,20, 4, 6,12, 2,36, 4, 6, 2,10, 2,22, 8, 6,10,12,12,18,14,24,36,
 4,20,24,10, 6, 2,28, 6,18, 8, 4, 6, 8, 6, 4, 2,12,28,18,14,16,14,18,10, 8,
 6, 4, 6, 6, 8,22,12, 2,10,18, 6, 2,18,10, 2,12,10,18,32, 6, 4, 6, 6, 8, 6,
 6,10,20, 6,12,10, 8,10,14, 6,10,14, 4, 2,22,18, 2,10, 2, 4,20, 4, 2,34, 2,
12, 6,10, 2,10,18, 6,14,12,12,22, 8, 6,16, 6, 8, 4,12, 6, 8, 4,36, 6, 6,20,
24, 6,12,18,10, 2,10,26, 6,16, 8, 6, 4,24,18, 8,12,12,10,18,12, 2,24, 4,12,
18,12,14,10, 2, 4,24,12,14,10, 6, 2, 6, 4, 6,26, 4, 6, 6, 2,22, 8,18, 4,18,
 8, 4,24, 2,12,12, 4, 2,52, 2,18, 6, 4, 6,12, 2, 6,12,10, 8, 4, 2,24,10, 2,
10, 2,12, 6,18,40, 6,20,16, 2,12, 6,10,12, 2, 4, 6,14,12,12,22, 6, 8, 4, 2,
16,18,12, 2, 6,16, 6, 2, 6, 4,12,30, 8,16, 2,18,10,24, 2, 6,24, 4, 2,22, 2,
16, 2, 6,12, 4,18, 8, 4,14, 4,18,24, 6, 2, 6,10, 2,10,38, 6,10,14, 6, 6,24,
 4, 2,12,16,14,16,12, 2, 6,10,26, 4, 2,12, 6, 4,12, 8,12,10,18, 6,14,28, 2,
 6,10, 2, 4,14,34, 2, 6,22, 2,10,14, 4, 2,16, 8,10, 6, 8,10, 8, 4, 6, 2,16,
 6, 6,18,30,14, 6, 4,30, 2,10,14, 4,20,10, 8, 4, 8,18, 4,14, 6, 4,24, 6, 6,
18,18, 2,36, 6,10,14,12, 4, 6, 2,30, 6, 4, 2, 6,28,20, 4,20,12,24,16,18,12,
14, 6, 4,12,32,12, 6,10, 8,10, 6,18, 2,16,14, 6,22, 6,12, 2,18, 4, 8,30,12,
 4,12, 2,10,38,22, 2, 4,14, 6,12,24, 4, 2, 4,14,12,10, 2,16, 6,20, 4,20,22,
12, 2, 4, 2,12,22,24, 6, 6, 2, 6, 4, 6, 2,10,12,12, 6, 2, 6,16, 8, 6, 4,18,
12,12,14, 4,12, 6, 8, 6,18, 6,10,12,14, 6, 4, 8,22, 6, 2,28,18, 2,18,10, 6,
14,10, 2,10,14, 6,10, 2,22, 6, 8, 6,16,12, 8,22, 2, 4,14,18,12, 6,24, 6,10,
 2,12,22,18, 6,20, 6,10,14, 4, 2, 6,12,22,14,12, 4, 6, 8,22, 2,10,12, 8,40,
 2, 6,10, 8, 4,42,20, 4,32,12,10, 6,12,12, 2,10, 8, 6, 4, 8, 4,26,18, 4, 8,
28, 6,18, 6,12, 2,10, 6, 6,14,10,12,14,24, 6, 4,20,22, 2,18, 4, 6,12, 2,16,
18,14, 6, 6, 4, 6, 8,18, 4,14,30, 4,18, 8,10, 2, 4, 8,12, 4,12,18, 2,12,10,
 2,16, 8, 4,30, 2, 6,28, 2,10, 2,18,10,14, 4,26, 6,18, 4,20, 6, 4, 8,18, 4,
12,26,24, 4,20,22, 2,18,22, 2, 4,12, 2, 6, 6, 6, 4, 6,14, 4,24,12, 6,18, 2,
12,28,14, 4, 6, 8,22, 6,12,18, 8, 4,20, 6, 4, 6, 2,18, 6, 4,12,12, 8,28, 6,
 8,10, 2,24,12,10,24, 8,10,20,12, 6,12,12, 4,14,12,24,34,18, 8,10, 6,18, 8,
 4, 8,16,14, 6, 4, 6,24, 2, 6, 4, 6, 2,16, 6, 6,20,24, 4, 2, 4,14, 4,18, 2,
 6,12, 4,14, 4, 2,18,16, 6, 6, 2,16,20, 6, 6,30, 4, 8, 6,24,16, 6, 6, 8,12,
30, 4,18,18, 8, 4,26,10, 2,22, 8,10,14, 6, 4,18, 8,12,28, 2, 6, 4,12, 6,24,
 6, 8,10,20,16, 8,30, 6, 6, 4, 2,10,14, 6,10,32,22,18, 2, 4, 2, 4, 8,22, 8,
18,12,28, 2,16,12,18,14,10,18,12, 6,32,10,14, 6,10, 2,10, 2, 6,22, 2, 4, 6,
 8,10, 6,14, 6, 4,12,30,24, 6, 6, 8, 6, 4, 2, 4, 6, 8, 6, 6,22,18, 8, 4, 2,
18, 6, 4, 2,16,18,20,10, 6, 6,30, 2,12,28, 6, 6, 6, 2,12,10, 8,18,18, 4, 8,
18,10, 2,28, 2,10,14, 4, 2,30,12,22,26,10, 8, 6,10, 8,16,14, 6, 6,10,14, 6,
 4, 2,10,12, 2, 6,10, 8, 4, 2,10,26,22, 6, 2,12,18, 4,26, 4, 8,10, 6,14,10,
 2,18, 6,10,20, 6, 6, 4,24, 2, 4, 8, 6,16,14,16,18, 2, 4,12, 2,10, 2, 6,12,
10, 6, 6,20, 6, 4, 6,38, 4, 6,12,14, 4,12, 8,10,12,12, 8, 4, 6,14,10, 6,12,
 2,10,18, 2,18,10, 8,10, 2,12, 4,14,28, 2,16, 2,18, 6,10, 6, 8,16,14,30,10,
20, 6,10,24, 2,28, 2,12,16, 6, 8,36, 4, 8, 4,14,12,10, 8,12, 4, 6, 8, 4, 6,
14,22, 8, 6, 4, 2,10, 6,20,10, 8, 6, 6,22,18, 2,16, 6,20, 4,26, 4,14,22,14,
 4,12, 6, 8, 4, 6, 6,26,10, 2,18,18, 4, 2,16, 2,18, 4, 6, 8, 4, 6,12, 2, 6,
 6,28,38, 4, 8,16,26, 4, 2,10,12, 2,10, 8, 6,10,12, 2,10, 2,24, 4,30,26, 6,
 6,18, 6, 6,22, 2,10,18,26, 4,18, 8, 6, 6,12,16, 6, 8,16, 6, 8,16, 2,42,58,
 8, 4, 6, 2, 4, 8,16, 6,20, 4,12,12, 6,12, 2,10, 2, 6,22, 2,10, 6, 8, 6,10,
14, 6, 6, 4,18, 8,10, 8,16,14,10, 2,10, 2,12, 6, 4,20,10, 8,52, 8,10, 6, 2,
10, 8,10, 6, 6, 8,10, 2,22, 2, 4, 6,14, 4, 2,24,12, 4,26,18, 4, 6,14,30, 6,
 4, 6, 2,22, 8, 4, 6, 2,22, 6, 8,16, 6,14, 4, 6,18, 8,12, 6,12,24,30,16, 8,
34, 8,22, 6,14,10,18,14, 4,12, 8, 4,36, 6, 6, 2,10, 2, 4,20, 6, 6,10,12, 6,
 2,40, 8, 6,28, 6, 2,12,18, 4,24,14, 6, 6,10,20,10,14,16,14,16, 6, 8,36, 4,
12,12, 6,12,50,12, 6, 4, 6, 6, 8, 6,10, 2,10, 2,18,10,14,16, 8, 6, 4,20, 4,
 2,10, 6,14,18,10,38,10,18, 2,10, 2,12, 4, 2, 4,14, 6,10, 8,40, 6,20, 4,12,
 8, 6,34, 8,22, 8,12,10, 2,16,42,12, 8,22, 8,22, 8, 6,34, 2, 6, 4,14, 6,16,
 2,22, 6, 8,24,22, 6, 2,12, 4, 6,14, 4, 8,24, 4, 6, 6, 2,22,20, 6, 4,14, 4,
 6, 6, 8, 6,10, 6, 8, 6,16,14, 6, 6,22, 6,24,32, 6,18, 6,18,10, 8,30,18, 6,
16,12, 6,12, 2, 6, 4,12, 8, 6,22, 8, 6, 4,14,10,18,20,10, 2, 6, 4, 2,28,18,
 2,10, 6, 6, 6,14,40,24, 2, 4, 8,12, 4,20, 4,32,18,16, 6,36, 8, 6, 4, 6,14,
 4, 6,26, 6,10,14,18,10, 6, 6,14,10, 6, 6,14, 6,24, 4,14,22, 8,12,10, 8,12,
18,10,18, 8,24,10, 8, 4,24, 6,18, 6, 2,10,30, 2,10, 2, 4, 2,40, 2,28, 8, 6,
 6,18, 6,10,14, 4,18,30,18, 2,12,30, 6,30, 4,18,12, 2, 4,14, 6,10, 6, 8, 6,
10,12, 2, 6,12,10, 2,18, 4,20, 4, 6,14, 6, 6,22, 6, 6, 8,18,18,10, 2,10, 2,
 6, 4, 6,12,18, 2,10, 8, 4,18, 2, 6, 6, 6,10, 8,10, 6,18,12, 8,12, 6, 4, 6,
14,16, 2,12, 4, 6,38, 6, 6,16,20,28,20,10, 6, 6,14, 4,26, 4,14,10,18,14,28,
 2, 4,14,16, 2,28, 6, 8, 6,34, 8, 4,18, 2,16, 8, 6,40, 8,18, 4,30, 6,12, 2,
30, 6,10,14,40,14,10, 2,12,10, 8, 4, 8, 6, 6,28, 2, 4,12,14,16, 8,30,16,18,
 2,10,18, 6,32, 4,18, 6, 2,12,10,18, 2, 6,10,14,18,28, 6, 8,16, 2, 4,20,10,
 8,18,10, 2,10, 8, 4, 6,12, 6,20, 4, 2, 6, 4,20,10,26,18,10, 2,18, 6,16,14,
 4,26, 4,14,10,12,14, 6, 6, 4,14,10, 2,30,18,22, 2,16, 2, 4, 8, 6, 6,16, 2,
 6,12,10, 8,12, 4,14, 4, 6,20,10,12, 2, 6, 6, 4, 2,10, 2,30,16,12,20,18, 4,
 6, 2, 4, 8,16,14,18,22, 6, 2,22, 6, 6,18, 2,10,36, 8, 4, 6,20, 4,12, 6,14,
 4, 2,28,24, 8, 4, 6,12,30,18,32,22, 8,36, 6, 4,12, 2,12, 4, 6,20,10,18,18,
 8, 6, 4,24, 8,10,14, 6, 4, 8,12,16, 2,16, 6, 8,16,12,14,10,30,14, 4,12, 8,
12, 6,10, 2,12,28, 6,12,12,20,10, 2,10,14, 6, 6,30, 4, 8,12, 4, 2,10,14, 4,
26,18,12,10, 6, 8, 4,12, 6,24,18, 8,10, 2,12, 4,12,12, 6, 2,22, 2, 4, 2,12,
16,14,10, 2,16,18,32, 4, 6,20,22, 8,10, 2,10, 6, 2, 4,14, 6,24, 4, 8, 4, 6,
12,12, 8, 6,10,12, 8,10, 2,10,12, 6,12,12,20,28,20,10,14,10, 8,10, 6, 2, 4,
14, 6, 6,12, 6,12,10,14,10,14,16, 8,10,26, 4, 2, 6, 4,14, 4, 6,12, 8, 6,30,
18,12, 6,12,16,12,12, 2,28, 6,14,10,36, 2, 4, 6, 8,12,22,18, 2,30,18,22,20,
18,10,38, 6, 4, 2,24, 4, 6, 6, 2,10, 6,14,10, 8, 4,24,14,16,14,22, 6,20,10,
14, 4,12,12, 2,16, 8, 6, 6,18, 4, 6,14,22, 6, 2,42,16, 2,10, 6, 2, 4, 6, 8,
10,20,16,30, 8,10, 8,10, 2,30, 6, 6,36,10, 8,16, 6, 2,12,28, 2, 4, 6,18,12,
 6, 8,10, 2, 4,50, 4,20, 4,30, 8, 4, 6,12, 2,24, 4, 8,18, 6, 4, 6, 8,10, 2,
 4, 2,40,18,36,30,30, 8,16,14, 6,12,28, 2,22, 2, 4,12,30,12, 6, 2, 4,14,10,
 2,18,22,12,18, 2,10,18,32, 6, 4, 2, 6,10,20,12,10, 6,12,20,12, 6, 4, 2,16,
 2,16, 6,14, 4, 2,16, 2, 6,16, 6, 8, 4, 8,22,18, 8,12, 4, 8, 6,24,22, 6, 2,
12,30, 6,10,12, 6, 2,22, 6, 2,12, 6,22, 8,12,22, 2,10, 6,18,12, 2, 6,12,18,
 6, 4,20,22, 8,12,24,16,14,10,30,18, 2, 6, 4,14,10, 2,12,10,12, 6, 2,16,12,
 2, 6,12,10, 2,10, 6, 2,12,12,16,20,10,12, 8,30,10,14, 4, 6, 8, 6, 4,20,18,
24, 4,12, 8, 4, 2,24, 6,24,10, 2, 4, 6, 2, 6, 6, 6, 4,24, 2,10,12, 2, 6,10,
 8, 6,10,18, 2, 6, 4,20,24,10,12, 2,12, 6,24, 4,36,14,16, 8,22, 6, 8, 4, 2,
 6,22,20,16,12,18, 2,12,16, 6, 6,12, 6,12, 2, 6,12,10, 8,16, 8, 6,16, 8,12,
 4, 6, 6,20,12,12, 4, 6,20, 4,12, 2,10, 2, 6,30,22, 6, 2, 4,38,10, 2, 4, 2,
22, 2,16, 2, 6,10,20, 6,24, 4,12,14,12, 4,38,10,30, 6, 2,12,12, 4, 6,30,14,
 4, 8,18,36, 4, 6,20, 4, 2,12,10, 2, 6,10,12, 6,12, 8, 6, 6,24, 4,30,20, 6,
36,10, 2,12, 6, 4, 8, 6, 4,12, 8, 6,12, 4, 6,14, 4,20,12, 4, 6,18, 2, 4,18,
 2,16,12,30, 6, 6, 8,40, 8,48, 6,16,18,14,12, 6,18, 4,20,10, 2, 6,10, 8,30,
 4,12,20, 6,12, 6, 6,34, 6, 6,18, 6, 8,10,12, 6, 8,10, 2, 4,24, 6, 8,22, 6,
 2,12, 6,10,12, 6,24, 6,14,12,36, 4,24, 2,10, 8,10, 6,14,10,32, 4, 8,10,12,
26,18, 4, 6,20, 4,20, 6,16, 6, 2,30,12, 6,10, 2, 6,10,12, 8, 4, 2, 6,10,12,
26,22, 8, 6, 4,14, 6, 6,30, 4, 6,14, 4, 2,28, 2, 6,22, 8, 4,18,18,18, 2,12,
 6, 4,20,10, 6, 6,14,10,12, 2,12,30,34,12, 8, 6, 4, 2,10, 2,16,12, 2,10, 8,
18,24, 6, 4,12,14, 4, 8, 4,14, 4, 6, 6,20, 6, 4, 8,18,52, 2, 4,12, 8, 4,38,
 4,26,24,16,12, 6, 2,12,12,16, 2, 6, 6, 4,12,14,16, 8,12,18,16, 6, 8,10, 6,
14,10,12, 2,10, 2, 4,24, 6,42,24, 8,10, 6, 6, 6, 2,12, 4,14, 6, 6,28, 6, 2,
10,12,12, 6,20, 4, 6,14, 4, 2,12,10,12,24, 6, 8, 6, 6, 4,24,12,20,16,14,30,
18, 6, 4,26,12, 4, 6, 2, 6, 4, 2,28, 8,40, 2,10, 8, 4,20, 6,18,10, 2, 4,44,
 6,18,12, 6, 4, 6, 2,22, 6,14,30,10,24, 2,10, 8,16,18, 2,18,22, 8,10, 6, 6,
14, 4, 8,18, 4, 2,18,18,18, 6, 4,24,18, 2,16, 6, 6,18,20,16,20, 4,14, 6, 4,
20,18,10, 2, 6,10,24, 2,10,24, 6, 6,24, 6,12, 2,28,12,14, 6, 6,12, 6,22,12,
12, 8,36, 4,12,14, 4,20,10,12,24, 2, 4, 6,12, 2, 4, 2,10,12,26, 6,16, 8, 4,
 8,10, 8, 6,34, 2,12,16,24, 6, 2,10, 2,18, 4, 8, 6,16, 6, 2, 6, 6, 6, 4,14,
 4,20, 6, 4,20, 6,12,22, 6, 2,10,12, 2, 6, 4, 8,12, 4,14,12,10,14, 4,12,26,
10,14, 4,26, 6,30, 4,18,18, 8, 6,16, 8,10,14,10, 8,10,20,22,20,16, 2,18, 6,
 4, 6, 6,12, 2,10,26, 4, 8,18,18, 6,18, 6, 4, 6,24, 6,20,34,26,10, 2,28,12,
 8,10,12, 2, 6,22, 2,12,16, 2, 6, 6,10,14,16,20, 6, 4,38, 6,10, 6, 8,16,42,
 2, 6, 4, 6, 6, 6,14,16,14, 4,20,10, 2, 4, 8,18,10,12,36, 2,10,42, 8, 4,20,
24,16, 8,22, 6, 8, 4, 2, 6,22, 6, 6, 8,28, 2,10,18,14, 6, 4,18, 8,10,14, 4,
12, 8,10,12,14, 4, 2,12,12, 4, 6,18,30,12,38, 6,12,10, 2,18,10,12, 8, 4, 8,
 6, 4, 2,24,12,18, 4, 2, 4, 2,58,12, 8,24,10, 2, 4, 6, 6,12, 2, 4,14, 6, 6,
16,12, 2, 4,32, 4,24, 6, 6, 8,10, 2,22,18,12,20, 6,30, 4,30, 6, 2, 4,14, 6,
 4,14,16, 2,12,10, 2, 6,12,12,10, 6, 8,22, 8,12,12, 6,16, 6,18,20,22,18, 2,
22, 2,16, 2,22,14,10,20,10,32, 4, 8,10, 6, 2,22, 6,12, 2, 6, 4, 2, 4,14,12,
24,10, 2,12,16, 2, 4, 6,14, 6,10,12, 2,16,14,34,12, 2, 6, 6, 6, 4,20,10,26,
12,12, 4, 2, 4, 8,10, 2, 4, 2,22, 6, 6,14, 4,18,12,26, 6,10, 8,16, 2, 4,20,
10, 6,42, 2,10, 6, 8,24,12, 6, 4, 6,12, 2,28, 8,12,18,18, 6,46, 8,10, 6,14,
 4, 2, 6, 4, 6,42, 8,10, 8,10, 2,18, 4, 6,12,12, 2, 4,20,10,12,12, 8, 4,26,
18,22, 8, 6,16,14,16, 2,18,10, 2, 6, 6,10,14, 4, 2,30, 4, 2, 4, 8,10, 6, 2,
12,16, 6,56,10, 2,12,10, 8,12, 6, 4,14,10, 2, 4, 8, 6, 4,20, 6,12,22, 6,32,
10, 2,10,12,14, 6,28,36, 6, 6, 2,12, 4, 6, 6, 8,22, 2,18,10, 2, 6, 4,20,10,
 8, 4, 6,14,18, 6,42,22, 2, 4, 2,28, 2, 4,18, 6, 6, 6,12, 2,24,10,36, 6, 2,
12,10,26,24,18,16, 6, 6,14,24,12, 4, 8, 6,12, 4, 8,16,20,40,26, 4,12, 2, 6,
 4, 2,10,14,10, 2, 4,26,12,28, 2,16,26, 6,10, 2, 6,10, 6, 8, 6, 6, 6,10,12,
 6,20,40,20, 4, 2,16,12, 6,12, 8, 4,18, 2,12,10,26,12,16, 2,18,24,12, 4,14,
22,20,10,14,12, 4,18,12, 8,10,12, 6,30,14, 4,24, 6,30, 6, 6, 2, 6,22,32, 6,
 4, 6, 6,20,16, 2,10, 8,12,10, 2, 6,10, 8,16,36, 8, 6, 4, 2,28, 2,28,12, 2,
10, 6,14,10, 6, 6, 6, 8, 6, 4,14,18, 4, 6,12, 2,10,18, 8,30,40, 2,18, 4, 6,
14,18, 6, 4,12, 6,12, 6,14,10,26, 6,16, 2,16,30, 2,10, 2,42, 6,28,14, 6,10,
 2,12,18,12, 6,10,12,12,20, 6, 4, 2,10, 6,12,12,14,12,34, 6, 2,12,10, 6, 8,
 6, 4,12,38, 6,10,18, 2,28, 2, 6,12,30,16, 2,10, 8, 4, 2,16,18,26, 4, 6, 8,
18,22, 6,20, 4, 6,12, 2, 6,12, 4,18, 6, 2,22,12, 8, 6,16,18,30,12,24, 2,10,
 2, 6, 6, 4, 6,36,14, 6,22, 2,58, 8,12, 6,10, 2,40, 8, 6,28, 2, 4,14, 6, 6,
18,10, 8, 4,14, 4, 8,30, 4, 6, 8, 6, 6,18, 4, 2, 4,14,12,18,10, 2, 4,12, 2,
10, 8,10,14,10,18,12, 8, 6,10,14,10, 8,22, 2, 6,22,12, 6, 8,12,28, 2,48,12,
 4,18, 8,10,14,10,14, 4,12,30,24, 6, 8, 6, 4, 8,54, 4, 2,10,12, 8,10,12,12,
18, 2,24, 4, 8,22,12,20, 4,12, 2,12,16, 2,28, 2, 6,24,10, 2,28, 2, 4,20, 4,
12, 6,14, 4, 6,14,22,24,20, 4,14, 6, 6,10,30, 8,10,18, 2, 6, 6,16, 2, 6, 6,
 4, 2,24, 4, 2,24,10, 6, 2,10, 2, 6,22, 8, 4, 8, 6, 4,18, 2,18, 4, 8,16,26,
 4, 6, 8,22,20,16, 8, 4, 6,24, 6,14,12,16, 2,12, 4,14,10, 2, 4,12,18,32,10,
14,24,12,40, 8,34,12,14, 4,18, 2,28,12,20, 6,10, 2,40,18,14,12, 4,36, 6, 2,
22, 6,14,10,24,42, 2,16, 2,34, 8, 6, 4, 2, 4,14,40, 8,12, 6,24,18, 4, 6, 2,
 6, 4, 2, 4, 2,24,10, 8, 6, 6,10,14, 6,16,18,14,18,24, 4, 6, 6, 8, 4,20,10,
 6,12, 2,12, 4,14, 6, 6, 6, 4,14,16,36,14, 6, 4,14, 4, 6,24, 8, 4,20,10,14,
12,34, 8,10, 6, 6, 6,14, 4,14,12, 6,10,18,14,10,12, 6, 2, 6, 6,28, 2, 4,24,
 6, 2, 4, 8,16, 6,20, 4, 2,10, 2,10, 8,64, 6, 8,12, 4,14,12,10, 2,12, 6,10,
18,24, 6, 2,10, 8, 6,16,20, 4,14, 6, 6,12, 6, 4, 6, 2, 4, 8,22, 6, 8, 4, 2,
16,18,14, 6,22,14,10,14, 4, 6, 2, 4,14,10,12, 8,16, 8,10, 8,24,40, 6,12, 2,
 6,18, 4, 2, 4,30, 2,30, 4, 8,18,12,12, 4, 2, 4,14,36,16,18, 2,12,10, 6,12,
18, 2,18, 6, 6,22,18,38, 6,10,18, 2,10, 8, 6,16,24,14, 6, 4, 6,14,16,24, 6,
12, 8,12,10,14,46, 2,16, 2,22, 6, 2,10, 2,10, 2, 6, 4,20,10, 6,30, 8, 6, 6,
 4,30, 8, 6, 6, 6,22,36, 2, 4, 8, 6, 6, 4,14,12,10,20, 4, 2, 4,30, 6,14,16,
12,30, 2, 4, 6, 8,30,10, 8,34,18,12, 8,22,20, 4,14,10,20, 6, 4, 2,10,14, 4,
26, 6,36,12,18, 4, 8, 6, 4, 6, 2,28, 6, 6,24, 8,10,26, 6,24, 4, 8,24,10,20,
 4, 2,10,14,16, 2, 6, 6, 4, 6, 8,18,28,14, 6,16,14, 6, 4, 6, 6, 8, 4, 2, 4,
12, 2,12, 6,12,28, 2, 6,12,10,14, 4,44, 6,10, 2,12,12,30, 4,12, 2, 6,10,12,
 2,10, 2,10, 6, 8,10, 6,14,16, 8, 6,12,10, 2,10, 8,12,10,18, 8, 4, 2, 4,26,
 6,22, 6,14,10, 6, 2,28, 6, 8,46, 6, 6,18, 6, 6, 8, 6,10,18, 2, 6,12,18,10,
 8,12,30,10, 2,10, 2, 4, 6,18, 2, 4,20,12, 4, 6, 8,34, 6, 6,24,12, 8,36,16,
 2, 6, 4, 2, 4, 6,20, 6,24, 4, 2, 4,18,20, 6,22, 8,46,18, 2,16,20,22, 2,24,
22, 2,16,24,20,16, 2, 4, 8,10, 2,10,14, 4, 8,18, 4, 8, 4,14,10, 2,24,16, 8,
 6,16,20,10, 2, 6, 4,30, 2,16,32, 6,12,10,24, 8,12,18,16, 2,12, 6, 4,12, 6,
 2,28,18, 2,22, 6, 6, 6, 2, 6,16,14, 6,30,16, 2,10, 2, 4,12, 2,12,10,14, 6,
10, 8,28, 2,36, 6,16,14, 4,20,24, 6, 4, 8, 4,18, 8, 4,14, 4, 6, 2,24,16,14,
 4,26,16, 2,10,32, 6, 4, 6,12, 6,36, 8,12, 4, 2, 4, 8, 6, 4,20,12,10,24,12,
 2,12,10, 6,12, 2, 6,18, 4, 6, 6, 6, 8,24, 6,10,12,30,14,10, 8,12, 6,10,12,
 2,18, 6, 4, 8, 4,24,20, 4, 8,10,12, 8,12,16, 6,14, 4, 8, 4,18,50, 6, 6, 4,
 6, 8, 6,10,26,10, 6, 2,10, 2,10, 6,38,12, 4, 8,10,20, 6, 6, 6,18,10, 2,12,
16, 2,12,12, 4,26,10, 6,20,18,40,12, 8,10,12, 2,18,12,10, 2,10,26, 4, 6,12,
 8, 4,30, 6, 2, 6,16,24,24,18,12,12, 8, 6, 4, 8,10, 8, 6, 4,20,10,26, 4,24,
 6, 2,12,42,18, 6, 4,26, 6,28, 6, 2,10, 8, 6, 6,10, 8,10, 2,22, 2, 4,20, 4,
 6,36,14, 4,20,22, 6,14, 6,10, 8, 4, 2, 4,14,18,34, 8,22,14,10,24, 6, 2,10,
 2, 6,10,26,18,10,18,24,18, 2,24,40, 2, 4, 6, 2, 6,10,26, 6,12,12, 6, 4,36,
 2,10,12,24, 2, 4, 8,10, 6, 2, 4,24, 2, 4,36, 2,22,14,24,18,42, 6,10, 2,24,
16,12, 2, 4, 2,10, 2,10, 8, 4,36, 8, 4,12,18, 6, 6,14,22, 2, 6,24, 6,10,24,
20,22, 6,14,36,28, 6, 8, 6,24, 6,12,28, 2,18, 4, 2, 4,20,22, 8,10, 2,18, 4,
 8,10,14,10, 6, 8, 6, 6,12,16,12,14,10,18, 2,10,24,24, 6,12, 2,22, 6,20,22,
 2, 4,12, 2, 6,36, 6,22, 6, 2,28,12,18, 2, 4,14, 6, 4, 2,10, 2,16, 2,10, 8,
 6,10,18,12, 6,14, 4, 6,18,12,26, 4, 6,14, 6,10,12, 2, 4, 2,10,24, 8,10,32,
10, 8,10, 6, 2,18,12,28,30, 2,18, 4, 6,14, 6, 4, 8,22, 8,30,18,10,26, 4, 2,
22, 8, 4, 8, 6, 4,26, 4,12,20,18, 6,12,10,18, 2, 4, 6, 2,12,28, 6,20, 6,16,
 8, 6, 6, 4, 6,20,12, 6, 4,20, 6,16, 6,32,10,18, 2
};

// everything below is from the msieve library

void analyze_one_poly_hook( long degree, mpz_t* acoeffs, mpz_t* rcoeffs, double skewness, double* size_score,
       double* root_score, double* combined_score, unsigned int* num_real_roots ) {

	uint32_t i;
	poly_config_t config;
	poly_select_t s;
	dd_precision_t prec = 0;
	uint32_t prec_changed = 0;

	if (!dd_precision_is_ieee()) {
		prec_changed = 1;
		prec = dd_set_precision_ieee();
	}

	poly_config_init(&config);
	poly_select_init(&s);

	s.skewness = skewness;

	s.rpoly.degree = 1;
	for (i = 0; i <= s.rpoly.degree; i++)
		mpz_set(s.rpoly.coeff[i], rcoeffs[i]);

	s.apoly.degree = degree;
	for (i = 0; i <= s.apoly.degree; i++)
		mpz_set(s.apoly.coeff[i], acoeffs[i]);

	analyze_poly(&config, &s);

        *size_score = s.size_score;
        *root_score = s.root_score;
        *combined_score = s.combined_score;
        *num_real_roots = s.num_real_roots;

	poly_config_free(&config);
	poly_select_free(&s);

	if (prec_changed)
		dd_clear_precision(prec);
}

/*------------------------------------------------------------------*/
void analyze_poly(poly_config_t *config, poly_select_t *poly) {

	/* analyze a polynomial for sieving goodness
	  
	   The analysis routines are general enough so that
	   any polynomials can be tested, independent of
	   degree and skewness. The score assigned is
	   directly comparable to that of any other polynomials 
	   given to this routine */

	uint32_t i;
	double root_score_r, root_score_a;
	ddpoly_t ddr, dda;
	mpz_poly_t *rpoly = &poly->rpoly;
	mpz_poly_t *apoly = &poly->apoly;

	poly->size_score = 0.0;
	poly->root_score = 0.0;
	poly->combined_score = 0.0;

	/* convert the polynomial coefficients from arbitrary
	   precision to dd_t floating point */

	ddr.degree = rpoly->degree;
	for (i = 0; i <= rpoly->degree; i++)
		ddr.coeff[i] = dd_gmp2dd(rpoly->coeff[i]);

	dda.degree = apoly->degree;
	for (i = 0; i <= apoly->degree; i++)
		dda.coeff[i] = dd_gmp2dd(apoly->coeff[i]);

	if (analyze_poly_roots(&poly->rpoly, PRIME_BOUND, &root_score_r))
		return;
	if (analyze_poly_roots(&poly->apoly, PRIME_BOUND, &root_score_a))
		return;

	if (analyze_poly_size(&config->integ_aux, 
				&ddr, &dda, &poly->size_score))
		return;

	get_bernstein_combined_score(poly, root_score_r + root_score_a);

	poly->root_score = root_score_a;

	analyze_poly_murphy(&config->integ_aux, &config->dickman_aux,
				&ddr, root_score_r, &dda, root_score_a, 
				poly->skewness, &poly->combined_score,
				&poly->num_real_roots);
}

uint32_t analyze_poly_roots(mpz_poly_t *poly, uint32_t prime_bound,
				double *result) {

	/* analyze a polynomial for root properties (i.e.
	   compute Murphy's 'alpha' value) */

	uint32_t i, j;
	double root_score;
	uint32_t prime;
	mpz_poly_t rev_poly;

	/* all linear polynomials have the same root properties;
	   note that we still must treat them as generating 
	   homogeneous polynomial values and *not* as random
	   numbers. A linear NFS polynomial will always have one
	   root modulo each prime p, leading to a fixed but
	   nonzero alpha value */

	if (poly->degree == 1) {
		*result = 0.569959993064325;
		return 0;
	}

	/* handling projective roots requires a reversed version
	   of poly (we find roots of poly(1/x) instead of poly(x)) */

	mpz_poly_init(&rev_poly);
	j = poly->degree;
	for (i = 0; i <= j; i++) {
		mpz_set(rev_poly.coeff[i], poly->coeff[j-i]);
	}
	while (j && mpz_cmp_ui(rev_poly.coeff[j], 0) == 0) {
		j--;
	}
	rev_poly.degree = j;

	/* roots of (poly mod p) can be either ordinary or special.
	   An ordinary root mod p corresponds to exactly one root
	   mod p^k for all k > 1. Murphy shows that for each ordinary 
	   root of (poly mod p) the contribution to the size of a
	   random sieve value due to all powers of p is log(p) * 
	   p/(p^2-1). 
	   
	   We compare this to the contribution of powers of p to 
	   the size of a random number. 1 in p random numbers are 
	   divisible by p, 1 in p^2 random numbers are divisible by
	   p^2, etc. Adding up all of these contributions and summing
	   the resulting geometric series yields an average contribution
	   of log(p) / (p-1). 
	   
	   Hence the root score for poly is the sum of the difference 
	   between the contribution of powers of p in sieve values 
	   compared to the contribution to a random number. When there 
	   are k ordinary roots mod p, root_score is the sum of 

	               log(p) * (1/(p-1) - k*p/(p^2-1))

	   for the first few p, meaning that dividing the small primes 
	   out of sieve values makes their size exp(root_score) times 
	   smaller than dividing the small primes out of a random number.

	   When the number of roots mod p^k changes with k, Murphy's 
	   formula for the contribution is inaccurate (usually it's too
	   small), so we have to take special measures to compute the 
	   contribution of the roots accurately. Roots become special 
	   if p divides the discriminant of poly, and basically the only 
	   effective way to compute the size of the contribution from 
	   special roots is to find the number of roots mod p^k for the
	   first few k and then explicitly sum the first few terms in 
	   the (formerly geometric) series. Many thanks to Kleinjung 
	   and Franke whose pol5 code handles special roots efficiently */

	root_score = 0;
	for (i = prime = 0; i < PRECOMPUTED_NUM_PRIMES; i++) {

		uint32_t roots[MAX_POLY_DEGREE];
		uint32_t mult[MAX_POLY_DEGREE];
		uint32_t high_coeff;
		uint32_t num_roots, num_ordinary_roots;
		double contrib;

		/* test if next prime will exceed the bound */

		if (prime + prime_delta[i] > SMALL_PRIME_BOUND || 
		    prime + prime_delta[i] > prime_bound)
			break;

		/* find the roots of poly mod prime and determine if
		   any are multiple roots (this is the telltale sign
		   that the root is special) */

		prime += prime_delta[i];
		num_roots = poly_get_zeros_and_mult(roots, mult, poly, 
					prime, &high_coeff);

		contrib = 1.0 / (prime - 1);

		/* test each root for specialness */

		for (j = num_ordinary_roots = 0; j < num_roots; j++) {
			if (mult[j] == 0)
				num_ordinary_roots++;
			else
				contrib -= get_root_freq(poly, roots[j], prime);
		}
		
		/* projective roots are always handled specially */

		if (high_coeff == 0)
			contrib -= get_root_freq(&rev_poly, 0, prime);

		/* add in the root score contribution for prime */

		root_score += log((double)prime) * (contrib - 
						(double)num_ordinary_roots * 
						prime / (prime * prime - 1));
	}

	/* for the larger primes, assume that all roots are
	   ordinary, so that only the number of roots modulo
	   prime i matter */

	for (; i < PRECOMPUTED_NUM_PRIMES; i++) {
		uint32_t num_roots;
		uint32_t dummy_roots[MAX_POLY_DEGREE];
		uint32_t high_coeff;

		prime += prime_delta[i];
		if (prime > prime_bound)
			break;

		num_roots = poly_get_zeros(dummy_roots, poly, 
					prime, &high_coeff, 1);
		if (high_coeff == 0)
			num_roots++;

		root_score += (1.0 - (double)num_roots * prime / 
					(prime + 1)) * log((double)prime) / 
					(prime - 1);
	}

	*result = root_score;
	mpz_poly_free(&rev_poly);
	return 0;
}

uint32_t analyze_poly_size(integrate_t *integ_aux,
			ddpoly_t *rpoly, ddpoly_t *apoly, 
			double *result) {

	uint32_t i, j;
	uint32_t rdeg, adeg;
	dparam_t params;
	dd_complex_t roots[2*MAX_POLY_DEGREE];
	uint32_t num_roots;
	double endpoints[2*MAX_POLY_DEGREE+2];
	uint32_t num_endpoints;

	rdeg = rpoly->degree;
	adeg = apoly->degree;
	params.power = -2.0 / (rdeg + adeg);
	params.rpoly = rpoly;
	params.apoly = apoly;

	/* Rather than use Murphy's method to rate the average
	   size of polynomial values, we use a result of D.J.
	   Bernstein: For rational poly R(x) and algebraic poly
	   A(x), the number of polynomial values where x is rational
	   and the size of A(x)*R(x) achieves a given value is 
	   proportional to the following superelliptic integral:

	   	dx / ((R(x) * A(x))^2) ^ (1/(deg(R(x))+deg(A(x))))

	   for x from -infinity to +infinity. Larger values of
	   this integral imply more values that get found by
	   a siever. Bernstein writes that this estimate is 
	   "extremely accurate", but has not elaborated to date.
	   The integration variable x refers to the a/b values used
	   by a siever, and it is both more realistic and simpler to
	   make the integration interval finite but large 
	   
	   The integrand takes on a wide range of values depending
	   on how close x is to a root of R(x)*A(x), so we have to
	   partition the integral */


	/* find the roots of R(x)*A(x) */

	if (find_poly_roots(params.rpoly->coeff, rdeg, roots)) {
		printf("rational poly rootfinder failed\n");
		return 1;
	}
	num_roots = rdeg;
	if (find_poly_roots(params.apoly->coeff, adeg, roots + num_roots)) {
		printf("algebraic poly rootfinder failed\n");
		return 1;
	}
	num_roots += adeg;

	/* squeeze out roots whose real part would not lie
	   in the integration interval, along with all complex 
	   conjugates of complex roots.

	   We cannot remove all complex roots because if the
	   imaginary part of the root is small then the value
	   of A(x)*R(x) is 'almost' a singularity at this point,
	   and ignoring it would throw off the numerical integrator */
	   
	for (i = j = 0; i < num_roots; i++) {
		if (roots[i].i.hi < 0.0 || 
		    roots[i].r.hi <= -INTEGRATE_LIMIT ||
		    roots[i].r.hi >= INTEGRATE_LIMIT)
			continue;
		endpoints[j++] = roots[i].r.hi;
	}
	endpoints[j] = -INTEGRATE_LIMIT;
	endpoints[j+1] = INTEGRATE_LIMIT;
	num_endpoints = j + 2;

	/* sort the endpoints into increasing order */

	qsort(endpoints, (size_t)num_endpoints,
			sizeof(double), compare_double);

	/* integrate */

	*result = 0.0;

	integrate_run(integ_aux, integrand,
			&params, endpoints, num_endpoints);

	/* test for convergence */

	if (integ_aux->result > 1 ||
	    integ_aux->error > SIZE_EPS * fabs(integ_aux->result)) {
		printf("integrator failed %e %e\n", 
				integ_aux->error,
				SIZE_EPS * fabs(integ_aux->result));
		return 2;
	}
	*result = integ_aux->result;
	return 0;
}

/*------------------------------------------------------------------*/
static double eval_dpoly(dpoly_t *poly,
			double x0, double xh, 
			double y0, double yh) {

	/* evaluate poly(x0+xh, y0+yh) using a Taylor series
	   centered at poly(x0,y0). We have to do this because
	   the numerical integrator will sample the integrand 
	   extremely close to poly(x,y)=0, and in that case
	   xh and yh will be of such different magnitude compared
	   to x0 and y0 that we will run into numerical trouble */

	uint32_t i, j, k;
	uint32_t deg = poly->degree;
	double *coeff = poly->coeff;
	double xbinom[MAX_POLY_DEGREE+1][MAX_POLY_DEGREE+1];
	double ybinom[MAX_POLY_DEGREE+1][MAX_POLY_DEGREE+1];
	double sum[MAX_POLY_DEGREE+1];
	double *xrow, *yrow;
	double res;

	xbinom[1][0] = xh;
	xbinom[1][1] = x0;
	for (i = 2; i <= deg; i++) {

		double *polyj = xbinom[i-1];
		double *polyk = xbinom[1];
		double *row = xbinom[i];

		for (j = 0; j <= i; j++)
			row[j] = 0;

		for (j = 0; j < i; j++) {
			row[j+0] += polyj[j] * polyk[0];
			row[j+1] += polyj[j] * polyk[1];
		}
	}

	ybinom[1][0] = yh;
	ybinom[1][1] = y0;
	for (i = 2; i <= deg; i++) {

		double *polyj = ybinom[i-1];
		double *polyk = ybinom[1];
		double *row = ybinom[i];

		for (j = 0; j <= i; j++)
			row[j] = 0;

		for (j = 0; j < i; j++) {
			row[j+0] += polyj[j] * polyk[0];
			row[j+1] += polyj[j] * polyk[1];
		}
	}

	xrow = xbinom[deg];
	yrow = ybinom[deg];
	for (i = 0; i <= deg; i++)
		sum[i] = xrow[i] * coeff[deg] + yrow[i] * coeff[0];

	for (i = 1; i < deg; i++) {

		uint32_t xdeg = i;
		uint32_t ydeg = deg - i;

		xrow = xbinom[xdeg];
		yrow = ybinom[ydeg];
		for (j = 0; j <= xdeg; j++) {
			for (k = 0; k <= ydeg; k++) {
				sum[j+k] += coeff[i] * xrow[j] * yrow[k];
			}
		}
	}

	res = sum[0];
	for (i = 1; i <= deg; i++)
		res += sum[i];

	return res;
}

static void get_bernstein_combined_score(poly_select_t *poly,
					double root_score) {

	/* An interesting problem: you have a root score and
	   a size score for each polynomial; how do you combine
	   them to get an overall rating? The value of root_score
	   implies that values of R(x)*A(x) effectively have size
	   R(x)*A(x)*exp(root_score) after dividing out the small
	   primes from sieve values. Technically this will
	   overestimate the size, because root_score includes a 
	   bias that accounts for random numbers; but this bias
	   is the same for all polynomials. The exponential term 
	   is constant, so pulling this expression out of the 
	   integrand used to calculate the size above gives: */

	poly->size_score *= pow(exp(root_score), -2.0 / 
				(poly->rpoly.degree + poly->apoly.degree));
}

/*------------------------------------------------------------------*/
uint32_t analyze_poly_murphy(integrate_t *integ_aux, dickman_t *dickman_aux,
			ddpoly_t *rpoly, double root_score_r,
			ddpoly_t *apoly, double root_score_a,
			double skewness, double *result,
			uint32_t *num_real_roots) {

	/* Given the skewness and root score for an NFS polynomial
	   pair, calculate the probability that an average sieve 
	   value in the sieving region has all rational (resp. algebraic)
	   factors less than rfb_limit (resp. afb_limit) 
	 
	   Ideally the sieving area and factor base limits should vary
	   with the size of the NFS input, but we fix them here to
	   be compatible with the code in pol51. That code uses a
	   trapezoidal approximation to the integral and computes
	   Dickman's function via linear interpolation from pre-tabulated
	   values. The following uses a full numerical integrator and 
	   the classical Dickman series instead, but the integrand is
	   so smooth most of the time that the end effect is the same.
	   This code uses about 90% fewer integrand evaluations though */

	uint32_t i, j;
	murphy_param_t params;
	dd_complex_t roots[2*MAX_POLY_DEGREE];
	double angles[2*MAX_POLY_DEGREE+2];
	uint32_t num_roots, num_angles;
	uint32_t rdeg = rpoly->degree;
	uint32_t adeg = apoly->degree;

	const double rfb_limit = 5000000;
	const double afb_limit = 10000000;
	const double sieve_area = 1e16;

	params.dickman_aux = dickman_aux;
	params.root_score_r = root_score_r;
	params.root_score_a = root_score_a;
	params.rfb_limit = rfb_limit;
	params.afb_limit = afb_limit;
	params.log_rfb_limit = log(rfb_limit);
	params.log_afb_limit = log(afb_limit);
	params.skew_x = sqrt(sieve_area * skewness);
	params.skew_y = params.skew_x / skewness;

	params.rpoly.degree = rdeg;
	for (i = 0; i <= rdeg; i++)
		params.rpoly.coeff[i] = rpoly->coeff[i].hi;

	params.apoly.degree = adeg;
	for (i = 0; i <= adeg; i++)
		params.apoly.coeff[i] = apoly->coeff[i].hi;

	/* find the roots of rpoly * apoly */

	if (find_poly_roots(rpoly->coeff, rdeg, roots)) {
		printf("rational poly rootfinder failed\n");
		return 1;
	}
	num_roots = rdeg;
	if (find_poly_roots(apoly->coeff, adeg, roots + num_roots)) {
		printf("algebraic poly rootfinder failed\n");
		return 1;
	}
	num_roots += adeg;

	*num_real_roots = 0;
	for (i = 0; i < adeg; i++) {
		if (dd_cmp_d(roots[rdeg+i].i, 0.0) == 0)
			(*num_real_roots)++;
	}

	/* convert the roots to angles between 0 and pi. Since the
	   integrator will skew values of x and y derived from these
	   angles, skew the roots the other way to compensate */

	for (i = j = 0; i < num_roots; i++) {

		if (roots[i].i.hi < 0.0)
			continue;

		angles[j++] = atan2(1.0, roots[i].r.hi / skewness);
	}
	angles[j] = 0.0;
	angles[j+1] = M_PI;
	num_angles = j + 2;

	/* sort in order of increasing angle */

	qsort(angles, (size_t)num_angles, sizeof(double), compare_double);

	integrate_run(integ_aux, murphy_integrand, 
			&params, angles, num_angles);

	*result = integ_aux->result / M_PI;
	return 0;
}

/*------------------------------------------------------------------*/
static double murphy_integrand(double r, double h, void *params) {

	murphy_param_t *p = (murphy_param_t *)params;
	double x0, xh, y0, yh;
	double polyval_r, polyval_a;

	/* we need to convert the angular measure (r+h), with
	   r and h of possibly very different magnitudes, into
	   cartesian measures (x0+xh, y0+yh) so that we can use
	   the homogeneous form of the NFS polynomials */

	if (fabs(h) > 0.1) {
		x0 = cos(r + h);
		y0 = sin(r + h);
		xh = yh = 0;
	}
	else {
		uint32_t i;
		double term, hpow;

		/* sum the Taylor series for cos(r+h) and sin(r+h)
		   simultaneously */

		x0 = cos(r);
		y0 = sin(r);
		hpow = h;
		xh = -y0 * h;
		yh = x0 * h;

		for (i = 0; i < NUM_RECIP_FACTORIALS; i += 2) {

			hpow *= -h;
			term = hpow * recip_factorial[i];

			xh += x0 * term;
			yh += y0 * term;

			hpow *= h;
			term = hpow * recip_factorial[i+1];

			xh -= y0 * term;
			yh += x0 * term;

			if (fabs(term) <= 1e-16 * (fabs(xh) + fabs(yh)))
				break;
		}
	}

	/* skew the coordinates */

	x0 *= p->skew_x;
	xh *= p->skew_x;
	y0 *= p->skew_y;
	yh *= p->skew_y;

	/* evaluate the NFS polynomials at the skewed coordinates */

	polyval_r = eval_dpoly(&p->rpoly, x0, xh, y0, yh);
	polyval_a = eval_dpoly(&p->apoly, x0, xh, y0, yh);

	/* polyval_[ra] give a measure of the size of sieve values
	   contained in a ray emanating from the origin at an angle 
	   of (r+h) radians. The integrand is the probability that
	   sieve values on this ray are [ra]fb-limit-smooth, after 
	   small primes are removed */

	return dickman(p->dickman_aux,
			(log(fabs(polyval_r)) + p->root_score_r) / 
				p->log_rfb_limit) *
	       dickman(p->dickman_aux,
			(log(fabs(polyval_a)) + p->root_score_a) / 
				p->log_afb_limit);
}

/*------------------------------------------------------------------*/
void poly_config_init(poly_config_t *config) {

	/* one-time initialization for polynomial search */

	uint32_t i;

	config->heap_num_filled = 0;
	for (i = 0; i < POLY_HEAP_SIZE; i++) {
		config->heap[i] = (poly_select_t *)xcalloc((size_t)1, 
					sizeof(poly_select_t));
		poly_select_init(config->heap[i]);
	}

	integrate_init(&config->integ_aux, SIZE_EPS,
			double_exponential);

	dickman_init(&config->dickman_aux);
}

/*------------------------------------------------------------------*/
void poly_config_free(poly_config_t *config) {

	uint32_t i;

	for (i = 0; i < POLY_HEAP_SIZE; i++) {
		poly_select_free(config->heap[i]);
		free(config->heap[i]);
	}
	integrate_free(&config->integ_aux);
	dickman_free(&config->dickman_aux);
}

/*------------------------------------------------------------------*/
void poly_select_init(poly_select_t *s) {

	memset(s, 0, sizeof(poly_select_t));
	mpz_poly_init(&s->rpoly);
	mpz_poly_init(&s->apoly);
}

void poly_select_free(poly_select_t *s) {

	mpz_poly_free(&s->rpoly);
	mpz_poly_free(&s->apoly);
}

/*------------------------------------------------------------------*/
uint32_t poly_get_zeros(uint32_t *zeros, mpz_poly_t *_f, 
			uint32_t p, uint32_t *high_coeff,
			uint32_t count_only) { 

        /* Find all roots of multiplicity 1 for polynomial _f,
	   when the coefficients of _f are reduced mod p. 
	   The leading coefficient of _f mod p is returned
	   
	   Make count_only nonzero if only the number of roots
	   and not their identity matters; this is much faster */

	poly_t g, f;
	uint32_t i, j, num_zeros;

	/* reduce the coefficients mod p */

	poly_reduce_mod_p(f, _f, p);
	*high_coeff = f->coef[_f->degree];

	/* bail out if the polynomial is zero */

	if (f->degree == 0)
		return 0;

	/* pull out roots of zero. We do this early to
	   avoid having to handle degree-1 polynomials
	   in later code */

	num_zeros = 0;
	if (f->coef[0] == 0) {
		for (i = 1; i <= f->degree; i++) {
			if (f->coef[i])
				break;
		}
		for (j = i; i <= f->degree; i++) {
			f->coef[i - j] = f->coef[i];
		}
		f->degree = i - j - 1;
		zeros[num_zeros++] = 0;
	}

	/* handle trivial cases */

	if (f->degree == 0) {
		return num_zeros;
	}
	else if (f->degree == 1) {
		uint32_t w = f->coef[1];

		if (count_only)
			return num_zeros + 1;

		if (w != 1) {
			w = mp_modinv_1(w, p);
			zeros[num_zeros++] = mp_modmul_1(p - f->coef[0], 
						w, p);
		}
		else {
			zeros[num_zeros++] = (f->coef[0] == 0 ? 
						0 : p - f->coef[0]);
		}
		return num_zeros;
	}

	/* the rest of the algorithm assumes p is odd, which
	   will not work for p=2. Fortunately, in that case
	   there are only two possible roots, 0 and 1. The above
	   already tried 0, so try 1 here */

	if (p == 2) {
		uint32_t parity = 0;
		for (i = 0; i <= f->degree; i++)
			parity ^= f->coef[i];
		if (parity == 0)
			zeros[num_zeros++] = 1;
		return num_zeros;
	}
	 
	/* Compute g = gcd(f, x^(p-1) - 1). The result is
	   a polynomial that is the product of all the linear
	   factors of f. A given factor only occurs once in
	   this polynomial */

	poly_xpow(g, 0, p-1, f, p);
	g->coef[0] = mp_modsub_1(g->coef[0], 1, p);
	poly_fix_degree(g);
	poly_gcd(g, f, p);

	/* no linear factors, no service */

	if (g->degree < 1 || count_only)
		return num_zeros + g->degree;

	/* isolate the linear factors */

	get_zeros_rec(zeros, 0, &num_zeros, g, p);
	return num_zeros;
}

/*------------------------------------------------------------------*/
static void poly_gcd(poly_t g_in, poly_t h_in, uint32_t p) { 

	poly_t g, h;

	/* make sure the first GCD iteration actually
	   performs useful work */

	if (g_in->degree > h_in->degree) {
		poly_cp(g, g_in);
		poly_cp(h, h_in);
	}
	else {
		poly_cp(h, g_in);
		poly_cp(g, h_in);
	}

	while ((h->degree > 0) || (h->coef[h->degree])) {
		poly_t r;
		poly_mod(r, g, h, p);
		poly_cp(g, h);
		poly_cp(h, r);
	}
	if (g->degree == 0)
		g->coef[0] = 1;
	poly_cp(g_in, g);
}

/*------------------------------------------------------------------*/
static void get_zeros_rec(uint32_t *zeros, uint32_t shift, 
			uint32_t *num_zeros, poly_t f, uint32_t p) {

	/* get the zeros of a poly, f, that is known to split
	   completely over Z/pZ. Many thanks to Bob Silverman 
	   for a neat implementation of Cantor-Zassenhaus splitting */

	poly_t g, xpow;
	uint32_t degree1, degree2;

	/* base cases of the recursion: we can find the roots
	   of linear and quadratic polynomials immediately */

	if (f->degree == 1) {
		uint32_t w = f->coef[1];
		if (w != 1) {
			w = mp_modinv_1(w, p);
			zeros[(*num_zeros)++] = mp_modmul_1(p - f->coef[0],w,p);
		}
		else {
			zeros[(*num_zeros)++] = (f->coef[0] == 0 ? 0 : 
							p - f->coef[0]);
		}
		return;
	}
	else if (f->degree == 2) {

		/* if f is a quadratic polynomial, then it will 
		   always have two distinct nonzero roots or else
		   we wouldn't have gotten to this point. The two 
		   roots are the solution of a general quadratic 
		   equation, mod p */

		uint32_t d = mp_modmul_1(f->coef[0], f->coef[2], p);
		uint32_t root1 = p - f->coef[1];
		uint32_t root2 = root1;
		uint32_t ainv = mp_modinv_1(
				mp_modadd_1(f->coef[2], f->coef[2], p),
				p);

		d = mp_modsub_1(mp_modmul_1(f->coef[1], f->coef[1], p),
				mp_modmul_1(4, d, p),
				p);
		d = mp_modsqrt_1(d, p);

		root1 = mp_modadd_1(root1, d, p);
		root2 = mp_modsub_1(root2, d, p);
		zeros[(*num_zeros)++] = mp_modmul_1(root1, ainv, p);
		zeros[(*num_zeros)++] = mp_modmul_1(root2, ainv, p);
		return;
	}

	/* For an increasing sequence of integers 's', compute 
	   the polynomial gcd((x-s)^(p-1)/2 - 1, f). If the result is
	   not g = 1 or g = f, this is a nontrivial splitting 
	   of f. References require choosing s randomly, but however
	   s is chosen there is a 50% chance that it will split f.
	   Since only 0 <= s < p is valid, we choose each s in turn;
	   choosing random s allows the possibility that the same
	   s gets chosen twice (mod p), which would waste time */

	while (shift < p) {
		poly_xpow(xpow, shift, (p-1)/2, f, p);

		poly_cp(g, xpow);
		g->coef[0] = mp_modsub_1(g->coef[0], 1, p);
		poly_fix_degree(g);

		poly_gcd(g, f, p);

		if (g->degree > 0)
			break;
		shift++;
	}

	/* f was split; repeat the splitting process on
	   the two halves of f. The linear factors of f are
	   either somewhere in x^((p-1)/2) - 1, in 
	   x^((p-1)/2) + 1, or 'shift' itself is a linear
	   factor. Test each of these possibilities in turn.
	   In the first two cases, begin trying values of s
	   strictly greater than have been tried thus far */

	degree1 = g->degree;
	get_zeros_rec(zeros, shift + 1, num_zeros, g, p);

	poly_cp(g, xpow);
	g->coef[0] = mp_modadd_1(g->coef[0], 1, p);
	poly_fix_degree(g);
	poly_gcd(g, f, p);
	degree2 = g->degree;

	if (degree2 > 0)
		get_zeros_rec(zeros, shift + 1, num_zeros, g, p);

	if (degree1 + degree2 < f->degree)
		zeros[(*num_zeros)++] = (shift == 0 ? 0 : p - shift);
}

/*------------------------------------------------------------------*/
static void poly_xpow(poly_t res, uint32_t shift, uint32_t n, 
			poly_t mod, uint32_t p) { 

	/* Modular exponentiation of polynomials with 
	   finite-field coefficients, i.e. res = (x-shift) ^ n % mod
	   with all polynomial coefficients reduced modulo p.
	   n is assumed nonzero */

	poly_t modnorm;
	uint32_t msw;
	uint32_t i, d;
	uint32_t buf[2*NUM_POLY_COEFFS] = {0};
	uint64_t psq;

	poly_make_monic(modnorm, mod, p);
	d = modnorm->degree;

	OP1(0) = shift;
	OP1(1) = 1;
	for (i = 0; i <= d; i++)
		MOD(i) = modnorm->coef[i];

	msw = 0x80000000;
	while (!(n & msw)) {
		msw >>= 1;
	}
	msw >>= 1;
	 
	psq = (uint64_t)p * (uint64_t)p;

	/* use left-to-right binary exponentiation, not
	   the right-to-left variety. For factor base generation
	   the base always has degree less than modnorm, and the 
	   left-to-right method preserves that, saving time during
	   modular multiplication */

	while (msw) {
		poly_expo_square(buf, d, p, psq);
		if (n & msw) {
			poly_expo_modmul(buf, d, shift, p, psq);
		}
		msw >>= 1;
	}

	res->degree = d;
	for (i = 0; i <= d; i++)
		res->coef[i] = OP1(i);
	poly_fix_degree(res);
}

/*------------------------------------------------------------------*/
static void poly_expo_modmul(uint32_t *buf, uint32_t dm, uint32_t shift,
			   uint32_t p, uint64_t psq) { 

	/* OP1 = OP1 * (x - shift) mod MOD
	   OP1 and MOD are of degree dm */

	uint32_t q;
	uint32_t zero = 0;

	q = OP1(dm-1);
	switch(dm-1) {
	case 7: OP1(7) = mul_mac(OP1(7), shift, OP1(6), q, MOD(7), p, psq);
	case 6: OP1(6) = mul_mac(OP1(6), shift, OP1(5), q, MOD(6), p, psq);
	case 5: OP1(5) = mul_mac(OP1(5), shift, OP1(4), q, MOD(5), p, psq);
	case 4: OP1(4) = mul_mac(OP1(4), shift, OP1(3), q, MOD(4), p, psq);
	case 3: OP1(3) = mul_mac(OP1(3), shift, OP1(2), q, MOD(3), p, psq);
	case 2: OP1(2) = mul_mac(OP1(2), shift, OP1(1), q, MOD(2), p, psq);
	case 1: OP1(1) = mul_mac(OP1(1), shift, OP1(0), q, MOD(1), p, psq);
	case 0: OP1(0) = mul_mac(OP1(0), shift,   zero, q, MOD(0), p, psq);
		break;
	}
}

/*------------------------------------------------------------------*/
static void poly_expo_square(uint32_t *buf, uint32_t dm, uint32_t p, uint64_t psq) { 

	/* OP1 = OP1 * OP1 mod MOD
	   OP1 and MOD are both of degree dm */

	uint32_t i;
	uint32_t q;
	uint64_t acc[NUM_POLY_COEFFS];

	for (i = 0; i < dm; i++)
		acc[i] = (uint64_t)(OP1(i)) * (uint64_t)(OP1(dm-1));

	for (i = dm - 2; (int32_t)i >= 0; i--) {
		q = mp_mod64(acc[dm-1], p);
		switch(dm-1) {
  		case 7: acc[7] = sqr_mac(OP1(7), OP1(i), acc[6], 
  							q, MOD(7), psq);
		case 6: acc[6] = sqr_mac(OP1(6), OP1(i), acc[5], 
							q, MOD(6), psq);
		case 5: acc[5] = sqr_mac(OP1(5), OP1(i), acc[4], 
							q, MOD(5), psq);
		case 4: acc[4] = sqr_mac(OP1(4), OP1(i), acc[3], 
							q, MOD(4), psq);
		case 3: acc[3] = sqr_mac(OP1(3), OP1(i), acc[2], 
							q, MOD(3), psq);
		case 2: acc[2] = sqr_mac(OP1(2), OP1(i), acc[1], 
							q, MOD(2), psq);
		case 1: acc[1] = sqr_mac(OP1(1), OP1(i), acc[0], 
							q, MOD(1), psq);
		case 0: acc[0] = sqr_mac0(OP1(0), OP1(i), q, MOD(0), psq);
			break;
		}
	}

	for (i = 0; i < dm; i++)
		OP1(i) = mp_mod64(acc[i], p);
}

/*------------------------------------------------------------------*/
uint32_t poly_get_zeros_and_mult(uint32_t *zeros, uint32_t *mult,
				mpz_poly_t *_f, uint32_t p,
				uint32_t *high_coeff) {

	uint32_t i;
	uint32_t num_roots;
	poly_t f;

	num_roots = poly_get_zeros(zeros, _f, p, high_coeff, 0);
	if (num_roots == 0)
		return num_roots;

	poly_reduce_mod_p(f, _f, p);
	for (i = 0; i < num_roots; i++)
		mult[i] = 0;
	if (f->degree == num_roots)
		return num_roots;

	for (i = 0; i < num_roots; i++) {

		poly_t g, r;
		uint32_t root = zeros[i];

		g->degree = 2;
		g->coef[0] = mp_modmul_1(root, root, p);
		g->coef[1] = p - mp_modadd_1(root, root, p);
		g->coef[2] = 1;

		poly_mod(r, f, g, p);
		if (r->degree == 0)
			mult[i] = 1;
	}
	return num_roots;
}

/*------------------------------------------------------------------*/
static void poly_reduce_mod_p(poly_t res, mpz_poly_t *_f, uint32_t p) {

	uint32_t i;

	res->degree = _f->degree;
	for (i = 0; i <= _f->degree; i++)
		res->coef[i] = mpz_fdiv_ui(_f->coeff[i], p);
	poly_fix_degree(res);
}

/*------------------------------------------------------------------*/
static void poly_fix_degree(poly_t op) { 

	int32_t i = op->degree;

	while ((i > 0) && (op->coef[i] == 0))
		i--;
	op->degree = i;
}

/*------------------------------------------------------------------*/
static void poly_mod(poly_t res, poly_t op, poly_t _mod, uint32_t p) { 

	/* divide the polynomial 'op' by the polynomial '_mod'
	   and write the remainder to 'res'. All polynomial
	   coefficients are reduced modulo 'p' */

	int32_t i;
	uint32_t msw;
	poly_t tmp, mod;

	if (_mod->degree == 0) {
		memset(res, 0, sizeof(res[0]));
		return;
	}
	poly_cp(tmp, op);
	poly_make_monic(mod, _mod, p);

	while (tmp->degree >= mod->degree) {

		/* tmp <-- tmp - msw * mod * x^{deg(tmp)- deg(mod)} */

		msw = tmp->coef[tmp->degree];

		tmp->coef[tmp->degree] = 0;
		for (i = mod->degree-1; i >= 0; i--) {
			uint32_t c = mp_modmul_1(msw, mod->coef[i], p);
			uint32_t j = tmp->degree - (mod->degree - i);
			tmp->coef[j] = mp_modsub_1(tmp->coef[j], c, p);
		}
		poly_fix_degree(tmp);
	}
	poly_cp(res, tmp);
}

/*------------------------------------------------------------------*/
static void poly_cp(poly_t dest, poly_t src) {

	dest[0] = src[0];
}

/*------------------------------------------------------------------*/
static void poly_make_monic(poly_t res, poly_t a, uint32_t p) {

	uint32_t i;
	uint32_t d = a->degree;
	uint32_t msw = a->coef[d];

	if (msw != 1) {
		msw = mp_modinv_1(msw, p);
		res->degree = d;
		res->coef[d] = 1;
		for (i = 0; i < d; i++)
			res->coef[i] = mp_modmul_1(msw, a->coef[i], p);
	}
	else {
		poly_cp(res, a);
	}
}

/*-----------------------------------------------------------------------*/
static uint32_t jenkins_traub(complex_t poly[], 
			uint32_t degree, complex_t roots[]) {

	/* main Jenkins-Traub driver; returns number 
	   of roots found */

	uint32_t i; 
	uint32_t roots_found;
	jt_t w;

	/* remove any zeros at the origin */

	for (i = degree, roots_found = 0; i; i--, roots_found++) {
		if (poly[i].r != 0.0 || poly[i].i != 0.0)
			break;

		roots[roots_found] = czero;
	}
	w.degree = i;

	/* initialize */

	for (i = 0; i <= w.degree; i++)
		w.poly[i] = poly[i];

	w.angle = complex(M_SQRT1_2, -M_SQRT1_2);
	w.angle_inc = complex(cos(94 * M_PI / 180), 
			   sin(94 * M_PI / 180));

	/* loop to find roots */

	for (; roots_found < degree; roots_found++) {

		if (find_one_root(roots + roots_found, &w) == 0)
			break;

		/* deflate the polynomial */

		(w.degree)--;
		for (i = 0; i <= w.degree; i++)
			w.poly[i] = w.poly_aux[i];
	}

	return roots_found;
}

static double get_root_freq(mpz_poly_t *poly, 
			uint32_t initial_root, uint32_t p) {
	
	/* find the number of roots of poly mod p^k for
	   all k where p^k <= PRIME_BOUND, then add up
	   the resulting contribution to random values of
	   poly */

	uint32_t i, j;
	uint32_t coeffs[MAX_POLY_DEGREE + 1];
	uint32_t roots[ROOT_BUF_SIZE];
	uint32_t new_roots[ROOT_BUF_SIZE];
	uint32_t max_power, power, degree;
	uint32_t num_roots, num_new_roots;
	double contrib;

	/* find the largest power of p that does not exceed 2^16 */

	max_power = 1;
	power = p;
	while (1) {
		uint32_t next_power = power * p;
		if (next_power > 65000)
			break;

		power = next_power;
		max_power++;
	}

	/* compute poly mod power */

	degree = poly->degree;
	for (i = 0; i <= degree; i++)
		coeffs[i] = mpz_fdiv_ui(poly->coeff[i], power);
	while (degree > 0 && coeffs[degree] == 0)
		degree--;

	/* there is only one root mod p */

	num_roots = 1;
	roots[0] = initial_root;
	power = p;
	contrib = 1.0 / power;

	for (i = 2; i <= max_power; i++) {

		uint32_t next_power = power * p;
		uint32_t x;

		/* search for roots mod p^i using the roots
		   mod p^(i-1) */

		for (j = num_new_roots = 0; j < num_roots; j++) {

			/* all roots mod p^i must have a known value
			   mod p^(i-1). One root mod p^(i-1) can be a
			   multiple root, and can give rise to several 
			   roots mod p^i, or to no roots at all */

			for (x = roots[j]; x < next_power; x += power) {

				uint32_t k;
				uint32_t polyval = coeffs[degree] % next_power;
				for (k = degree; k; k--) {
					polyval = (polyval * x + 
						coeffs[k-1]) % next_power;
				}
				if (polyval == 0 && 
				    num_new_roots < ROOT_BUF_SIZE) {
					new_roots[num_new_roots++] = x;
				}
			}
		}

		if (num_new_roots == 0)
			break;

		/* add in the contribution due to roots mod p^i 
		   and make the new roots into the old roots */

		contrib += (double)num_new_roots / next_power;
		memcpy(roots, new_roots, num_new_roots * sizeof(uint32_t));
		power = next_power;
		num_roots = num_new_roots;
	}

	/* for ordinary roots, contrib would be 1/(p-1)
	   asymptotically */

	return contrib * p / (p + 1);
}

/*------------------------------------------------------------------*/
uint32_t find_poly_roots(dd_t *poly, uint32_t degree, dd_complex_t *roots) {

	uint32_t i;
	dd_complex_t ddcoeffs[MAX_ROOTFINDER_DEGREE + 1];
	complex_t dcoeffs[MAX_ROOTFINDER_DEGREE + 1];
	complex_t droots[MAX_ROOTFINDER_DEGREE + 1];

	if (degree == 1) {
		roots[0].r = dd_div_dd(dd_neg(poly[0]), poly[1]);
		roots[0].i = dd_set_d(0.0);
		return 0;
	}

	for (i = 0; i <= degree; i++) {
		ddcoeffs[i].r = poly[i];
		ddcoeffs[i].i = dd_set_d(0.0);
		dcoeffs[degree - i] = complex(poly[i].hi, 0.0);
	}

	/* find the roots to a relative error close to the
	   double-precision limit */

	if (jenkins_traub(dcoeffs, degree, droots) != degree)
		return 1;

	/* polish each root */

	for (i = 0; i < degree; i++) {

		if (polish_root(ddcoeffs, degree,
				cplx_set_d(droots[i].r, droots[i].i),
				roots + i, 1e-30) != 0)
			return 2;

		/* change roots with very small imaginary part to
		   be explicitly real roots */

		if (dd_cmp_dd(dd_fabs(roots[i].i),
			      dd_mul_d(dd_fabs(roots[i].r), 1e-30)) <= 0) {
			roots[i].i = dd_set_d(0.0);
		}
	}

	return 0;
}

static uint32_t polish_root(dd_complex_t *poly, uint32_t degree,
			dd_complex_t x, dd_complex_t *root,
			double eps) {

	uint32_t i = 0;
	double eps2 = eps * eps;

	for (i = 0; i < NEWTON_ITER; i++) {

		uint32_t j = degree;
		dd_complex_t f = poly[j];
		dd_complex_t df = cplx_set_d(0.0, 0.0);
		dd_complex_t dx;
		dd_t abs_x, abs_dx;

		for (j--; (int32_t)j >= 0; j--) {
			df = cplx_add(cplx_mul(df, x), f);
			f = cplx_add(cplx_mul(f, x), poly[j]);
		}
		dx = cplx_div(f, df);
		x = cplx_sub(x, dx);

		abs_x = dd_add_dd(dd_mul_dd(x.r, x.r),
				  dd_mul_dd(x.i, x.i));
		abs_dx = dd_add_dd(dd_mul_dd(dx.r, dx.r),
				  dd_mul_dd(dx.i, dx.i));

		if (dd_cmp_dd(abs_dx, dd_mul_d(abs_x, eps2)) <= 0)
			break;
	}

	*root = x;
	return 0;
}

/*-----------------------------------------------------------------------*/
static complex_t complex(double re, double im) {

	complex_t res;

	res.r = re;
	res.i = im;
	return res;
}

/*-----------------------------------------------------------------------*/
static int find_one_root(complex_t *root, jt_t *w) {

	uint32_t i, j, k;
	double bound;
	complex_t hpoly_start[MAX_ROOTFINDER_DEGREE + 1];

	/* find linear roots immediately */

	if (w->degree <= 1) {
		*root = cdiv(cneg(w->poly[1]), w->poly[0]);
		return 1;
	}

	/* calculate a lower bound on the modulus of the zeros */

	bound = cauchy_bound(w->degree, w->poly);

	/* stage 1 sets up the initial h polynomial only */

	stage1(w->degree, w->poly, hpoly_start);

	/* try the fixed-shift sequence twice */

	for (i = 0; i < 2; i++) {

		/* inner loop to select a shift */

		for (j = 1; j < 10; j++) {

			/* start point is chosen with modulus 'bound'
			   and a pseudo-random angle. In practice
			   we don't want to repeat previous work,
			   so the starting angle is rotated a fixed 
			   amount (94 degrees) from the previous 
			   start point */

			w->angle = cmul(w->angle, w->angle_inc);
			*root = cscale(bound, w->angle);

			/* do the second stage, with a varying
			   number of iterations.
			   
			   Note that every starting point uses the same
			   h polynomial. This is a change from all other
			   cpoly() versions, including the original 1972 
			   fortran, which uses a global h array that is
			   not reinitialized when a new start point is
			   chosen (I think erroneously) */

			for (k = 0; k < w->degree; k++)
				w->hpoly[k] = hpoly_start[k];

			if (stage2(10 * j, root, w) == 1)
				return 1;
		}
	}

	return 0;
}

static void stage1(uint32_t n, complex_t p[], complex_t h[]) {

	uint32_t i, j;

	/* the initial h polynomial is a scaled version of the
	   derivative of the input polynomial p(x) */

	for (i = 0; i < n; i++)
		h[i] = cscale((double) (n - i) / n, p[i]);

	/* compute a series of no-shift h polynomials */

	for (i = 0; i < STAGE1_ITER; i++) {

		/* if the constant term is essentially zero, 
		   shift the h coefficients */

		if (cmod(h[n-1]) <= 10.0 * DBL_EPSILON * cmod(p[n-1])) {

			for (j = n - 1; j; j--)
				h[j] = h[j-1];
			h[j] = czero;
		}
		else {

			complex_t tmp = cdiv(cneg(p[n]), h[n-1]);
			for (j = n - 1; j; j--)
				h[j] = cmac(h[j-1], tmp, p[j]);
			h[j] = p[0];
		}
	}
}

/*-----------------------------------------------------------------------*/
static uint32_t stage2(uint32_t stage2_iter, complex_t *root, jt_t *w) {

	uint32_t i;
	complex_t curr_root;
	complex_t correction;
	complex_t pval;

	/* calculate first correction */

	curr_root = *root;
	pval = poly_val(w->degree, curr_root, w->poly, w->poly_aux);
	correction = next_correction(pval, curr_root, w);

	for (i = 0; i < stage2_iter; i++) {

		/* compute next h polynomial and new correction;
		   note that the fixed-shift iteration never changes
		   the value of curr_root, only the h polynomial */

		next_hpoly(correction, w);
		correction = next_correction(pval, curr_root, w);

		if (w->hpoly_root_found == 1)
			break;
	}

	/* attempt stage 3 with the final h polynomial and
	   final correction */

	*root = cadd(curr_root, correction);
	return stage3(root, w);
}

static complex_t cadd(complex_t x, complex_t y) {
	return complex(x.r + y.r, x.i + y.i);
}

static complex_t cneg(complex_t x) {
	return complex(-x.r, -x.i);
}

static complex_t cscale(double a, complex_t b) {
	return complex(a * b.r, a * b.i);
}

static complex_t cmul(complex_t a, complex_t b) {
	return complex(a.r * b.r - a.i * b.i, 
			a.r * b.i + a.i * b.r);
}

static complex_t cmac(complex_t a, complex_t x, complex_t b) {
	return complex(a.r * x.r - a.i * x.i + b.r,
		       a.r * x.i + a.i * x.r + b.i);
}

static double cmod(complex_t x) {

	/* modulus of a complex number avoiding overflow */

	double re = fabs(x.r); 
	double im = fabs(x.i); 

	if (re == im)
		return re * M_SQRT2;

	if (re < im) {
		double t = re / im;
		return im * sqrt(1.0 + t * t);
	}
	else {
		double t = im / re;
		return re * sqrt(1.0 + t * t);
	}
}

static complex_t cdiv(complex_t x, complex_t y) {

	/* complex division avoiding overflow */

	double zr, zi;
	double xr = x.r;
	double xi = x.i;
	double yr = y.r;
	double yi = y.i;

	if (yr == 0 && yi == 0) {
		zr = zi = DBL_MAX;
	}
	else if (fabs(yr) < fabs(yi)) {
		double u = yr / yi;
		double v = yi + u * yr;

		zr = (xr * u + xi) / v;
		zi = (xi * u - xr) / v;
	}
	else {
		double u = yi / yr;
		double v = yr + u * yi;

		zr = (xr + xi * u) / v;
		zi = (xi - xr * u) / v;
	}

	return complex(zr, zi);
}

/*-----------------------------------------------------------------------*/
static double cauchy_bound(uint32_t n, const complex_t p[]) {

	/* computes a lower bound on the moduli of the zeros of a
	   polynomial p(x). norms(x) is a polynomial whose i_th coefficient
	   is the modulus of the i_th coeffcient of p(x) but
	   whose constant term is negative. The lower bound is the 
	   (unique) positive root x of norms(x) */

	double x, xmax, f, dx, df;
	double norms[MAX_ROOTFINDER_DEGREE + 1];
	uint32_t i;

	for (i = 0; i < n; i++)
		norms[i] = cmod(p[i]);
	norms[i] = -cmod(p[i]);

	/* compute upper estimate of bound: assume all the
	   middle terms of norms(x) are zero */

	xmax = exp((log(-norms[n]) - log(norms[0])) / (double) n);

	/* if ignoring the nonlinear terms of norms(x) produces
	   a smaller root, use that instead */

	if (norms[n - 1] != 0.0) {
		x = -norms[n] / norms[n - 1];
		xmax = MIN(x, xmax);
	}

	/* chop the interval (0, x) until until x is about
	   to make norms(x) change sign */

	do {
		x = xmax;
		xmax = 0.1 * x;

		f = norms[0];
		for (i = 1; i <= n; i++)
			f = f * xmax + norms[i];
	} while (f > 0.0);

	/* do newton iteration until x converges to two decimal places */

	dx = x;
	while (fabs(dx / x) > 0.005) {
		df = 0;
		f = norms[0];
		for (i = 1; i <= n; i++) {
			df = df * x + f;
			f = f * x + norms[i];
		}
		dx = f / df;
		x -= dx;
	}

	return x;
}

/*-----------------------------------------------------------------------*/
static complex_t poly_val(uint32_t n, complex_t s, 
			const complex_t p[], complex_t q[]) {

	/* evaluates a polynomial p at s by the horner 
	   recurrence, placing the partial sums in q and 
	   returning the computed value */

	complex_t pv;
	uint32_t i;

	pv = q[0] = p[0];
	for (i = 1; i <= n; i++)
		pv = q[i] = cmac(pv, s, p[i]);

	return pv;
}

/*-----------------------------------------------------------------------*/
static void next_hpoly(complex_t correction, jt_t * w) {

	/* calculates the next shifted h polynomial */

	uint32_t i;
	complex_t *poly_aux = w->poly_aux;
	complex_t *hpoly_aux = w->hpoly_aux;
	complex_t *hpoly = w->hpoly;

	if (w->hpoly_root_found == 0) {
		hpoly[0] = poly_aux[0];
		for (i = 1; i < w->degree; i++) {
			hpoly[i] = cmac(hpoly_aux[i - 1], correction, 
						poly_aux[i]);
		}
	}
	else {
		/* we are essentially at a root of h(x); remove 
		   it by deflating the polynomial. Calling code 
		   always expects h(x) to have the same degree, 
		   so the high-order coefficient becomes zero */

		hpoly[0] = czero;
		for (i = 1; i < w->degree; i++)
			hpoly[i] = hpoly_aux[i - 1];
	}
}

/*-----------------------------------------------------------------------*/
static complex_t next_correction(complex_t pval, complex_t curr_root,
				jt_t * w) {

	/* computes -pval / hpoly(curr_root)
	   sets flag to true if hval is essentially zero. */

	complex_t *hpoly = w->hpoly;
	complex_t *hpoly_aux = w->hpoly_aux;
	complex_t hval = poly_val(w->degree - 1, curr_root, hpoly, hpoly_aux);

	if (cmod(hval) <= 10.0 * DBL_EPSILON * cmod(hpoly[w->degree - 1])) {
		w->hpoly_root_found = 1;
		return czero;
	}
	else {
		w->hpoly_root_found = 0;
		return cdiv(cneg(pval), hval);
	}
}

static uint32_t stage3(complex_t *root, jt_t *w) {

	/* carries out the third stage iteration,
	   returns 1 if iteration converges */

	double mp, ms, tp;
	uint32_t i, j;
	complex_t pval;
	complex_t correction;
	complex_t curr_root = *root;
	complex_t *poly = w->poly;
	complex_t *poly_aux = w->poly_aux;

	for (i = 0; i < STAGE3_ITER; i++) {

		/* evaluate poly at current root value */

		pval = poly_val(w->degree, curr_root, poly, poly_aux);

		/* calculate bound on the error in evaluating the polynomial 
		   by the horner recurrence */

		mp = cmod(pval);
		ms = cmod(curr_root);
		tp = cmod(poly_aux[0]) * mre / (DBL_EPSILON + mre);
		for (j = 0; j <= w->degree; j++)
			tp = tp * ms + cmod(poly_aux[j]);
		tp = tp * (DBL_EPSILON + mre) - mp * mre;

		if (mp <= 20.0 * tp) {
			/* polynomial value is smaller in value 
			   than a bound on the error in evaluating p, 
			   terminate the iteration */
			*root = curr_root;
			return 1;
		}

		/* calculate next h polynomial */

		correction = next_correction(pval, curr_root, w);
		next_hpoly(correction, w);

		/* use the next h polynomial to calculate the next
		   root estimate, using the current root estimate */

		correction = next_correction(pval, curr_root, w);
		curr_root = cadd(curr_root, correction);
	}

	return 0;
}

/*---------------------------------------------------------------*/
uint32_t mp_modsqrt_1(uint32_t a, uint32_t p) {

	uint32_t a0 = a;

	if ( (p & 7) == 3 || (p & 7) == 7 ) {
		return mp_expo_1(a0, (p+1)/4, p);
	}
	else if ( (p & 7) == 5 ) {
		uint32_t x, y;
		
		if (a0 >= p)
			a0 = a0 % p;
		x = mp_expo_1(a0, (p+3)/8, p);

		if (mp_modmul_1(x, x, p) == a0)
			return x;

		y = mp_expo_1(2, (p-1)/4, p);

		return mp_modmul_1(x, y, p);
	}
	else {
		uint32_t d0, d1, a1, s, t, m;
		uint32_t i;

		if (a0 == 1)
			return 1;

		for (d0 = 2; d0 < p; d0++) {
			if (mp_legendre_1(d0, p) != -1)
				continue;
	
			t = p - 1;
			s = 0;
			while (!(t & 1)) {
				s++;
				t = t / 2;
			}
	
			a1 = mp_expo_1(a0, t, p);
			d1 = mp_expo_1(d0, t, p);
	
			for (i = 0, m = 0; i < s; i++) {
				uint32_t ad;
	
				ad = mp_expo_1(d1, m, p);
				ad = mp_modmul_1(ad, a1, p);
				ad = mp_expo_1(ad, (uint32_t)(1) << (s-1-i), p);
				if (ad == (p - 1))
					m += (1 << i);
			}
	
			a1 = mp_expo_1(a0, (t+1)/2, p);
			d1 = mp_expo_1(d1, m/2, p);
			return mp_modmul_1(a1, d1, p);
		}
	}

	printf("modsqrt_1 failed\n");
	exit(-1);
}

/*---------------------------------------------------------------*/
int32_t mp_legendre_1(uint32_t a, uint32_t p) {

	uint32_t x, y, tmp;
	int32_t out = 1;

	x = a;
	y = p;
	while (x) {
		while ((x & 1) == 0) {
			x = x / 2;
			if ( (y & 7) == 3 || (y & 7) == 5 )
				out = -out;
		}

		tmp = x;
		x = y;
		y = tmp;

		if ( (x & 3) == 3 && (y & 3) == 3 )
			out = -out;

		x = x % y;
	}
	if (y == 1)
		return out;
	return 0;
}

/*------------------------------------------------------------------*/
static int compare_double(const void *x, const void *y) {
	double *xx = (double *)x;
	double *yy = (double *)y;

	if (*xx > *yy)
		return 1;
	if (*xx < *yy)
		return -1;
	return 0;
}

/*------------------------------------------------------------------*/
static double get_polyval(ddpoly_t *poly, double x, double h) {

	/* evaluate poly(x+h) using a taylor series centered
	   at poly(x). We compute poly(x) and the correction
	   from poly(x) separately, because in practice x and h
	   can have such different magnitudes that x+h could
	   equal x to machine precision. */

	double base; 
	double off;
	double hpow;
	dd_t *p = poly->coeff;
	double p0, p1, p2, p3, p4, p5, p6, p7, p8;

	switch (poly->degree) {
	case 0:
		p0 = p[0].hi;

		base = p0;
		off = 0;
		break;
	case 1:
		p0 = p[0].hi;
		p1 = p[1].hi;

		base = p1*x+p0;
		off = p1*h;
		break;
	case 2:
		p0 = p[0].hi;
		p1 = p[1].hi;
		p2 = p[2].hi;

		base = (p2*x+p1)*x+p0;
		hpow = h;
		off = hpow * (2*p2*x+p1);
		hpow *= h;
		off += hpow * p2;
		break;
	case 3:
		p0 = p[0].hi;
		p1 = p[1].hi;
		p2 = p[2].hi;
		p3 = p[3].hi;

		base = ((p3*x+p2)*x+p1)*x+p0;
		hpow = h;
		off = hpow * ((3*p3*x+2*p2)*x+p1);
		hpow *= h;
		off += hpow * (3*p3*x+p2);
		hpow *= h;
		off += hpow * p3;
		break;
	case 4:
		p0 = p[0].hi;
		p1 = p[1].hi;
		p2 = p[2].hi;
		p3 = p[3].hi;
		p4 = p[4].hi;

		base = (((p4*x+p3)*x+p2)*x+p1)*x+p0;
		hpow = h;
		off = hpow * (((4*p4*x+3*p3)*x+2*p2)*x+p1);
		hpow *= h;
		off += hpow * ((6*p4*x+3*p3)*x+p2);
		hpow *= h;
		off += hpow * (4*p4*x+p3);
		hpow *= h;
		off += hpow * p4;
		break;
	case 5:
		p0 = p[0].hi;
		p1 = p[1].hi;
		p2 = p[2].hi;
		p3 = p[3].hi;
		p4 = p[4].hi;
		p5 = p[5].hi;

		base = ((((p5*x+p4)*x+p3)*x+p2)*x+p1)*x+p0;
		hpow = h;
		off = hpow * ((((5*p5*x+4*p4)*x+3*p3)*x+2*p2)*x+p1);
		hpow *= h;
		off += hpow * (((10*p5*x+6*p4)*x+3*p3)*x+p2);
		hpow *= h;
		off += hpow * ((10*p5*x+4*p4)*x+p3);
		hpow *= h;
		off += hpow * (5*p5*x+p4);
		hpow *= h;
		off += hpow * p5;
		break;
	case 6:
		p0 = p[0].hi;
		p1 = p[1].hi;
		p2 = p[2].hi;
		p3 = p[3].hi;
		p4 = p[4].hi;
		p5 = p[5].hi;
		p6 = p[6].hi;

		base = (((((p6*x+p5)*x+p4)*x+p3)*x+p2)*x+p1)*x+p0;
		hpow = h;
		off = hpow * (((((6*p6*x+5*p5)*x+4*p4)*x+3*p3)*x+2*p2)*x+p1);
		hpow *= h;
		off += hpow * ((((15*p6*x+10*p5)*x+6*p4)*x+3*p3)*x+p2);
		hpow *= h;
		off += hpow * (((20*p6*x+10*p5)*x+4*p4)*x+p3);
		hpow *= h;
		off += hpow * ((15*p6*x+5*p5)*x+p4);
		hpow *= h;
		off += hpow * (6*p6*x+p5);
		hpow *= h;
		off += hpow * p6;
		break;
 	case 7:
		p0 = p[0].hi;
		p1 = p[1].hi;
		p2 = p[2].hi;
		p3 = p[3].hi;
		p4 = p[4].hi;
		p5 = p[5].hi;
		p6 = p[6].hi;
		p7 = p[7].hi;

 		base = ((((((p7*x+p6)*x+p5)*x+p4)*x+p3)*x+p2)*x+p1)*x+p0;
 		hpow = h;
 		off = hpow * ((((((7*p7*x+6*p6)*x+5*p5)*x+4*p4)*x+
 						3*p3)*x+2*p2)*x+p1);
 		hpow *= h;
 		off += hpow * (((((21*p7*x+15*p6)*x+10*p5)*x+
 						6*p4)*x+3*p3)*x+p2);
 		hpow *= h;
 		off += hpow * ((((35*p7*x+20*p6)*x+10*p5)*x+4*p4)*x+p3);
 		hpow *= h;
 		off += hpow * (((35*p7*x+15*p6)*x+5*p5)*x+p4);
 		hpow *= h;
 		off += hpow * ((21*p7*x+6*p6)*x+p5);
 		hpow *= h;
 		off += hpow * (7*p7*x+p6);
 		hpow *= h;
 		off += hpow * p7;
 		break;
	case 8:
		p0 = p[0].hi;
		p1 = p[1].hi;
		p2 = p[2].hi;
		p3 = p[3].hi;
		p4 = p[4].hi;
		p5 = p[5].hi;
		p6 = p[6].hi;
		p7 = p[7].hi;
		p8 = p[8].hi;

		base = (((((((p8*x+p7)*x+p6)*x+p5)*x+p4)*x+p3)*x+p2)*x+p1)*x+p0;
		hpow = h;
		off = hpow * (((((((8*p8*x+7*p7)*x+6*p6)*x+5*p5)*x+4*p4)*x+
						3*p3)*x+2*p2)*x+p1);
		hpow *= h;
		off += hpow * ((((((28*p8*x+21*p7)*x+15*p6)*x+10*p5)*x+
						6*p4)*x+3*p3)*x+p2);
		hpow *= h;
		off += hpow * (((((56*p8*x+35*p7)*x+20*p6)*x+10*p5)*x+4*p4)*x+p3);
		hpow *= h;
		off += hpow * ((((70*p8*x+35*p7)*x+15*p6)*x+5*p5)*x+p4);
		hpow *= h;
		off += hpow * (((56*p8*x+21*p7)*x+6*p6)*x+p5);
		hpow *= h;
		off += hpow * ((28*p8*x+7*p7)*x+p6);
		hpow *= h;
		off += hpow * (8*p8*x+p7);
		hpow *= h;
		off += hpow * p8;
		break;
	default:
		base = off = 0;
		break;
	}
	return base + off;
}

/*------------------------------------------------------------------*/
static double integrand(double x, double h, void *params) {

	/* callback for numerical integration */

	dparam_t *aux = (dparam_t *)params;
	double polyval;

	polyval = get_polyval(aux->rpoly, x, h) *
		  get_polyval(aux->apoly, x, h);
	return pow(fabs(polyval), aux->power);
}

/*-----------------------------------------------------------------------*/
static void de_run_core(integrate_t *aux, integrand_t func, 
			void *params, subinterval_t *range) {

	/* compute the integral of func(x, params) in the
	   range defined by 'range' */

	uint32_t i, j, k;
	uint32_t num_groups;
	de_t *integ = (de_t *)aux->internal;
	de_point_t *points;
	uint32_t num_points = integ->num_points;
	double lower_limit = range->left;
	double upper_limit = range->right;
	double interval = upper_limit - lower_limit;
	double step_size = 1.0;
	double result;
	double residual;
	double curr_error;
	double target_error;
	double left_val, right_val;

	/* evaluate func at the middle of the interval */

	result = FUNC_EVAL(0.5 * (lower_limit + upper_limit), 0.0);
	residual = result * integ->center_point.residual_weight;
	result *= integ->center_point.weight;

	/* compute the trapezoid rule approximation at the coarsest
	   step size. Sample at abscissa values that are symmetric 
	   about the origin */

	points = integ->initial_rule;

	for (i = 0; i < num_points; i++) {
		de_point_t *p = points + i;
		double abscissa = interval * p->abscissa;
		left_val = FUNC_EVAL(lower_limit, abscissa);
		right_val = FUNC_EVAL(upper_limit, -abscissa);
		result += (left_val + right_val) * p->weight;
		residual += (left_val + right_val) * p->residual_weight;
	}

	/* now compute trapezoid rules corresponding to successively
	   halved step sizes, until the estimate of the error indicates
	   convergence to the desired error tolerance */

	points = integ->refined_rules;
	i = 0;
	num_groups = 1;

	do {
		double old_result = result;
		double old_residual = residual;

		/* compute the trapezoid rule with step size 2^-(i+1) by
		   adding in 2^i groups of additional sample points.
		   Each sample falls somewhere in between two previous
	           samples, and all groups march from near the middle
		   of the interval out to the integration endpoints */

		for (j = 0; j < num_groups; j++, points += num_points) {

			/* update the integral value  and residual */

			for (k = 0; k < num_points; k++) {
				de_point_t *p = points + k;
				double abscissa = interval * p->abscissa;
				left_val = FUNC_EVAL(lower_limit, abscissa);
				right_val = FUNC_EVAL(upper_limit, -abscissa);
				result += (left_val + right_val) * p->weight;
				residual += (left_val + right_val) * 
							p->residual_weight;
			}
		}

		/* use the computed integral value and residual, at
		   the two different step sizes, to update the error 
		   estimate */

		curr_error = fabs(result - 2.0 * old_result) +
			     fabs(residual - 2.0 * old_residual);

		target_error = 0.1 * integ->relative_error * fabs(result);
		step_size *= 0.5;
		num_groups *= 2;

	} while (++i < MAX_TRAP_LEVELS && curr_error > target_error);

	/* compute the final integral value */

	range->result = result * interval * step_size;
	range->error = 2 * curr_error * interval;

#if 0
	printf("done: %le %le res %le err %le lvl %u\n", 
			lower_limit, upper_limit, 
			range->result, range->error, i-1);
#endif
}

static void heapify(subinterval_t *h, uint32_t index, uint32_t size) {

	uint32_t c;
	subinterval_t tmp;
	for (c = HEAP_LEFT(index); c < (size-1); 
			index = c, c = HEAP_LEFT(index)) {

		if (h[c].error < h[c+1].error)
			c++;

		if (h[index].error < h[c].error) {
			HEAP_SWAP(h[index], h[c]);
		}
		else {
			return;
		}
	}
	if (c == (size-1) && h[index].error < h[c].error) {
		HEAP_SWAP(h[index], h[c]);
	}
}

static void heap_insert(subinterval_t *h, uint32_t size) {
	uint32_t i = size - 1;
	subinterval_t tmp = h[i];

	while ((i > 0) && (h[HEAP_PARENT(i)].error < tmp.error) ) {
		h[i] = h[HEAP_PARENT(i)];
		i = HEAP_PARENT(i);
	}

	h[i] = tmp;
}

static void make_heap(subinterval_t *h, uint32_t size) {

	int32_t i;
	for (i = HEAP_PARENT(size); i >= 0; i--)
		heapify(h, (uint32_t)i, size);
}

static uint32_t de_run(integrate_t *aux, integrand_t func, 
			void *params, double *endpoints,
			uint32_t num_endpoints) {

	uint32_t i;
	de_t *integ = (de_t *)aux->internal;
	double result = 0;
	double error = 0;
	double rel_err = integ->relative_error;

	if (num_endpoints < 2)
		return 1;

	if (num_endpoints >= integ->num_heap_alloc) {
		integ->num_heap_alloc = MAX(num_endpoints + 100,
						2 * integ->num_heap_alloc);
		integ->heap = (subinterval_t *)xrealloc(integ->heap,
						integ->num_heap_alloc *
						sizeof(subinterval_t));
	}

	for (i = 0; i < num_endpoints - 1; i++) {
		subinterval_t *s = integ->heap + i;

		s->left = endpoints[i];
		s->right = endpoints[i+1];
		s->level = 0;
		de_run_core(aux, func, params, s);
		result += s->result;
		error += s->error;
	}
	integ->num_heap = i;

	if (error < rel_err * fabs(result)) {
		aux->result = result;
		aux->error = error;
		return 0;
	}

	make_heap(integ->heap, integ->num_heap);

	while (error > rel_err * fabs(result)) {

		subinterval_t *s = integ->heap;
		double left = s->left;
		double right = s->right;
		uint32_t level = s->level;

		if (level > MAX_LEVELS)
			break;

		result -= s->result;
		error -= s->error;
		s->left = left;
		s->right = 0.5 * (left + right);
		s->level = level + 1;
		de_run_core(aux, func, params, s);
		result += s->result;
		error += s->error;
		heapify(integ->heap, 0, integ->num_heap);

		if (integ->num_heap == integ->num_heap_alloc) {
			integ->num_heap_alloc *= 2;
			integ->heap = (subinterval_t *)xrealloc(integ->heap,
						integ->num_heap_alloc *
						sizeof(subinterval_t));
		}

		s = integ->heap + integ->num_heap++;
		s->left = 0.5 * (left + right);
		s->right = right;
		s->level = level + 1;
		de_run_core(aux, func, params, s);
		result += s->result;
		error += s->error;
		heap_insert(integ->heap, integ->num_heap);
	}

	aux->result = result;
	aux->error = error;
	return 0;
}

void integrate_init(integrate_t *aux, double relative_error,
			enum integrator_type type) {

	aux->type = type;
	if (type == double_exponential)
		de_init(aux, relative_error);
}

void integrate_free(integrate_t *aux) {

	if (aux->type == double_exponential)
		de_free(aux);
}

uint32_t integrate_run(integrate_t *aux, 
			integrand_t func, void *params,
			double *endpoints, uint32_t num_endpoints) {

	aux->result = 0.0;
	aux->error = 0.0;
	if (aux->type == double_exponential)
		return de_run(aux, func, params, 
				endpoints, num_endpoints);
	return 1;
}

/*-----------------------------------------------------------------------*/
static void de_init(integrate_t *aux, double relative_error) {

	/* initialize for all future integrations using aux */

	/* parameters of the algorithm: 
	   
	   min_val is the smallest function value we expect to
	   encounter

	   abscissa_scale expands the DE region from (-1,1) to 
	   (-abscissa_scale,abscissa_scale), presumably to avoid 
	   truncation error */

	const double min_val = 1e-30;
	const double abscissa_scale = 8.5;

	uint32_t i, j;
	de_t *integ;
	de_point_t *curr_samples;
	uint32_t num_points;
	double log_min_val = -log(min_val);
	double log_relative_error = 1.0 - log(relative_error);
	double abscissa_range = abscissa_scale / log_relative_error;
	double expo_large = exp(abscissa_range);
	double expo_small = 1.0 / expo_large;

	/* fill error bounds */

	integ = (de_t *)xmalloc(sizeof(de_t));
	integ->relative_error = relative_error;

	/* describe the initial sample point, in the middle 
	   of the interval */

	integ->center_point.weight = M_PI_2 * abscissa_range * 0.5;
	integ->center_point.residual_weight = abscissa_range;

	/* figure out the maximum number of samples needed */

	num_points = integ->num_points = (uint32_t)(log(log_min_val / M_PI_2) /
						 abscissa_range) + 1;

	curr_samples = integ->initial_rule = (de_point_t *)xmalloc(
					(num_points << MAX_TRAP_LEVELS) *
					sizeof(de_point_t));

	/* create the initial set of trapezoid sample points, with
	   step size of 1.0 * abscissa_range */

	de_fill_sample_group(curr_samples, num_points, 1.0,
				abscissa_range, expo_large, expo_small);

	/* now create trapezoid sample points in order of increasing
	   level. Level i has 2^i groups of samples; the j_th group
	   contains samples that start a distance of 
	   abscissa_range*j*2^-(i+1) from the origin and march away 
	   to inifinity in transformed coordinates at a given rate
	   (identical for all groups). Combining all the samples 
	   at level i, along with initial_rule and all samples at 
	   levels < i, produces the complete set of samples needed 
	   to implement the trapezoid rule with step size 2^-(i+1) */

	curr_samples += num_points;
	integ->refined_rules = curr_samples;

	for (i = 0; i < MAX_TRAP_LEVELS; i++) {
		uint32_t num_groups = 1 << i;
		double step_size = 1.0 / num_groups;
		double sample_val = 0.5 * step_size;

		for (j = 0; j < num_groups; j++) {
			de_fill_sample_group(curr_samples, num_points, 
					sample_val, abscissa_range, 
					expo_large, expo_small);
			sample_val += step_size;
			curr_samples += num_points;
		}
	}

	/* set up the heap of subintervals */

	integ->num_heap = 0;
	integ->num_heap_alloc = 100;
	integ->heap = (subinterval_t *)xmalloc(integ->num_heap_alloc *
						sizeof(subinterval_t));

	aux->internal = (void *)integ;
}

/*-----------------------------------------------------------------------*/
static void de_free(integrate_t *aux) {

	de_t *i = (de_t *)(aux->internal);
	free(i->initial_rule);
	free(i->heap);
	free(i);
	aux->internal = NULL;
}

static void de_fill_sample_group(de_point_t *rule, 
			uint32_t num_points, double first_sample, 
			double abscissa_range,
			double expo_large, double expo_small) {

	/* fill the information for the num_points trapezoid
	   panels, with the first abscissa at (first_sample * 
	   abscissa_range), 0 < first_sample <= 1. The sample 
	   points march off to infinity in transformed coordinates;
	   they are stored in untransformed coordinates, as offsets 
	   from +-abscissa_range, converging super-exponentially 
	   quickly to zero */

	uint32_t i;
	double curr_expo = exp(first_sample * abscissa_range);
	double curr_expo_large = M_PI_2 * curr_expo;
	double curr_expo_small = M_PI_2 / curr_expo;

	for (i = 0; i < num_points; i++) {
		de_point_t *curr_point = rule + i;
		double abscissa = 1.0 / (1.0 + exp(curr_expo_large -
						   curr_expo_small));
		double weight = abscissa * (1.0 - abscissa) * abscissa_range;

		curr_point->abscissa = abscissa;
		curr_point->weight = weight * (curr_expo_large +
						curr_expo_small);
		curr_point->residual_weight = 4.0 * weight;

		curr_expo_large *= expo_large;
		curr_expo_small *= expo_small;
	}
}


/*------------------------------------------------------------------*/
void dickman_init(dickman_t *aux) {

	uint32_t i, j, k;
	uint32_t num_coeffs;
	uint32_t num_coeffs_alloc;
	double last_coeffs[MAX_COEFFS];
	double *coeffs;
	double sum, s, t;

	/* initialize */

	num_coeffs = 0;
	num_coeffs_alloc = 1000;
	aux->coeffs = (double *)xmalloc(num_coeffs_alloc *
					sizeof(double));
	aux->lines = (dickman_line_t *)xmalloc((MAX_LINE + 1) *
					sizeof(dickman_line_t));

	/* dickman(x) is 1.0 for x <= 1; for x <= 2 the
	   value of dickman(x) is known explicitly */

	last_coeffs[0] = 1.0 - M_LN2;
	t = 0.5;
	for (i = 1; i < MAX_COEFFS; i++, t *= 0.5) {
		last_coeffs[i] = t / i;
	}

	/* save only enough coefficients for x=2 to 
	   achieve the desired accuracy, but use all the
	   rest for below to avoid roundoff error for
	   larger x 

	   We know the largest argument to the series
	   version of dickman(x) will be 1.0 and all terms
	   are positive, so this bounds the work needed */

	coeffs = aux->coeffs;
	sum = coeffs[0] = last_coeffs[0];
	num_coeffs = MAX_COEFFS;

	for (i = 1; i < MAX_COEFFS; i++) {
		coeffs[i] = last_coeffs[i];
		sum += coeffs[i];

		if (coeffs[i] < sum * 0.5 * DICKMAN_ACCURACY) {
			num_coeffs = i + 1;
			break;
		}
	}
	aux->lines[2].num_coeffs = num_coeffs;
	aux->lines[2].coeff_offset = 0;

	/* proceed with the rest of the integer aguments
	   to dickman(x) */

	for (i = 3; i <= MAX_LINE; i++) {
		dickman_line_t *line = aux->lines + i;
		double recip_di = 1.0 / i;

		if (num_coeffs + MAX_COEFFS >= num_coeffs_alloc) {
			num_coeffs_alloc *= 2;
			aux->coeffs = (double *)xrealloc(aux->coeffs,
							num_coeffs_alloc *
							sizeof(double));
		}

		/* derive the coefficients for dickman(x) from
		   those used in dickman(x-1) */

		sum = 0;
		for (j = MAX_COEFFS - 1; j; j--) {

			s = 0;
			t = recip_di / j;
			for (k = j - 1; (int32_t)k >= 0; k--) {
				s += t * last_coeffs[k];
				t *= recip_di;
			}

			last_coeffs[j] = s;
			sum += s / (j + 1);
		}
		last_coeffs[j] = sum / (i - 1);

		/* save enough coefficients to achieve the
		   desired accuracy for fractional arguments */

		coeffs = aux->coeffs + num_coeffs;
		sum = coeffs[0] = last_coeffs[0];
		line->coeff_offset = num_coeffs;
		line->num_coeffs = MAX_COEFFS;

		for (j = 1; j < MAX_COEFFS; j++) {
			coeffs[j] = last_coeffs[j];
			sum += coeffs[j];

			if (coeffs[j] < sum * 0.5 * DICKMAN_ACCURACY) {
				line->num_coeffs = j + 1;
				break;
			}
		}
		num_coeffs += line->num_coeffs;
	}
}

/*------------------------------------------------------------------*/
void dickman_free(dickman_t *aux) {

	free(aux->lines);
	free(aux->coeffs);
}

/*------------------------------------------------------------------*/
double dickman(dickman_t *aux, double arg) {

	uint32_t i;
	double int_arg, frac_arg;
	double sum, term;
	dickman_line_t *line;
	uint32_t num_coeffs;
	double *coeffs;

	if (arg <= 1.0)
		return 1.0;

	if (arg > MAX_LINE)
		return 0.0;

	int_arg = ceil(arg);
	frac_arg = int_arg - arg;
	line = aux->lines + (int32_t)int_arg;
	num_coeffs = line->num_coeffs;
	coeffs = aux->coeffs + line->coeff_offset;

	sum = coeffs[0];
	term = frac_arg;

	for (i = 1; i < num_coeffs; i++) {

		double iterm = term * coeffs[i];

		sum += iterm;

		if (iterm < sum * 0.5 * DICKMAN_ACCURACY)
			break;

		term *= frac_arg;
	}

	return sum;
}



