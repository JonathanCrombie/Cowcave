/*
 * This source distribution is placed into the public domain by its author, Jonathan Crombie.
 * You may use it for whatever purpose, without charge, and without having to notify anyone.
 * I disclaim any responsibility for any errors or lack of fitness for any purpose.
 *
 * cyclo Version 1.3.5  Dec 19, 2023
 *
 * (Dec 19, 2023. Add a check early on to filter out non-LMs.)
 * (Oct 13, 2022. Smarter ComputeLM. Eliminate mpz_pow() call from within for loops.)
 * (Jun 30, 2016. Removed the pari dependency.)
 * (Sep 22, 2014. Improved EulerTotient() -- skip computing all composite factors.)
 * (Jul 16, 2013. Changes version 1.03 --> 1.02:
 *  Added some error messages and did some refactoring.)
 * (Jul 1, 2013.  Changes version 1.01 --> 1.02:
 *  some optimizations.  Instead of trial factoring up to sqrt(n) for RemoveSquares(),
 *  ComputeAllFactors() etc., now only compute the prime factorization and then permute
 *  factors to compute the remaining composite factors.)
 *
 * (Dec 9, 2012.  Changes version 1.00 --> 1.01:
 *  minor tweak -- modify pari_init() to not build a large number of prime differences.
 *    this speeds up initialization quite a bit.)
 *
 * To build try:
 * g++ cyclo.cpp -lgmp -o cyclo
 *
 * (You will need to have g++ and the gmp library already installed)
 * (For linux bash shell for Windows 10, try:
    sudo apt-get install g++ libgmp3-dev                           )
 *
 * This source file is currently located at:
 * http://myfactors.mooo.com/source/cyclo.cpp
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <values.h>
#include <gmp.h>
#include <string.h>
#include <ctype.h>

struct factor_info {
mpz_t           the_factor;
long            occurrences;
long            digits;
char            factor_status;
};

struct factor_infos {
long                 count;
struct factor_info*  the_factors;
};


void ConvertPerfectPowers( mpz_t, mpz_t, mpz_t, mpz_t, mpz_t );
int IsLM( mpz_t, mpz_t, char );
int ComputeLM( mpz_t, mpz_t, mpz_t, char, mpz_t, mpz_t );
int ComputeLucasCD( mpz_t, mpz_t**, long*, mpz_t**, long* );
char quickprimecheck( mpz_t );
void TFDivideOut( mpz_t, long, char*, mpz_t, struct factor_infos* );
void Factorize( mpz_t, struct factor_infos* );
void PermuteFactors( struct factor_infos*, struct factor_infos*, long, mpz_t );
void ComputeAllFactors( mpz_t, struct factor_infos* );
long ComputeOccurrences( mpz_t, mpz_t );
int Mobius( mpz_t );
void EulerTotient( mpz_t, mpz_t );
void Init_Factor_Infos( struct factor_infos* );
void AddFactorInfo( struct factor_infos*, mpz_t, long, char, long );
void Cleanup_Factor_Infos( struct factor_infos* );
long NumberOfDigits( mpz_t );
void mpz_pow( mpz_t, mpz_t, mpz_t );
int factor_info_cmpfunc( const void*, const void* );


unsigned char wheel7[48] = { 2, 4, 2, 4, 6, 2, 6, 4, 2, 4,
                             6, 6, 2, 6, 4, 2, 6, 4, 6, 8,
                             4, 2, 4, 2, 4, 8, 6, 4, 6, 2,
                             4, 6, 2, 6, 6, 4, 2, 4, 6, 2,
                             6, 4, 2, 4, 2, 10, 2, 10 };


int main(int argc, char* argv[] )
{
  if ( argc != 4 ) {
    printf("\n");
    printf("cyclo Version 1.3.5\n\n");
    printf("Usage: \"cyclo base exponent +|-\"");
    printf("\n\n");
    return -1;
  }

  mpz_t base;
  mpz_init( base );
  mpz_set_str( base, argv[1], 10 );

  if ( mpz_cmp_ui( base, 2 ) < 0 ) {
    printf( "\n\nbase must be >= 2.  Aborting.\n\n" );
    return -1;
  }

  mpz_t exponent;
  mpz_init( exponent );
  mpz_set_str( exponent, argv[2], 10 );

  char c0 = argv[3][0];

  if ( c0 != '-' && c0 != '+' ) {
    printf( "\n\nError: Third argument must be either \'-\' or \'+\'.  Aborting.\n\n" );
    mpz_clear( exponent );
    mpz_clear( base );
    exit( -1 );
  }

  { // Do a quick test to weed out non-LMs early on
    mpz_t gcd;
    mpz_init( gcd );
    mpz_gcd( gcd, base, exponent );
    if ( mpz_cmp_ui( gcd, 1 ) == 0 ) {
      gmp_printf("\n\n%Zd^%Zd %c 1 is not an LM number.  Aborting.\n\n", base, exponent, c0 );
      mpz_clear( gcd );
      mpz_clear( exponent );
      mpz_clear( base );
      return 0;
    }
    mpz_clear( gcd );
  }

  mpz_t ConvertedBase;
  mpz_init( ConvertedBase );

  mpz_t ConvertedExponent;
  mpz_init( ConvertedExponent );

  mpz_t SqrFreeBase;
  mpz_init( SqrFreeBase );

  ConvertPerfectPowers( base, exponent, ConvertedBase, ConvertedExponent, SqrFreeBase );

  if ( IsLM( SqrFreeBase, ConvertedExponent, c0 ) ) {
    mpz_t L;
    mpz_init( L );
    mpz_t M;
    mpz_init( M );

    ComputeLM( ConvertedBase, SqrFreeBase, ConvertedExponent, c0, L, M );

    gmp_printf( "L=%Zd\n", L );
    gmp_printf( "M=%Zd\n", M );

    mpz_clear( M );
    mpz_clear( L );
  }
  else {
    gmp_printf("\n\n%Zd^%Zd %c 1 is not an LM number.  Aborting.\n\n", base, exponent, c0 );
  }

  mpz_clear( ConvertedBase );
  mpz_clear( ConvertedExponent );
  mpz_clear( SqrFreeBase );
  mpz_clear( exponent );
  mpz_clear( base );

  return 0;
}

// Need to fix base, exponent if base is a perfect power
// Compute the square-free part of base as well, since it's really easy to do.
void ConvertPerfectPowers( mpz_t base, mpz_t exponent, mpz_t ConvertedBase, mpz_t ConvertedExponent, mpz_t SqrFreeBase ) {
  struct factor_infos PrimeFactorization;
  Init_Factor_Infos( &PrimeFactorization );

  Factorize( base, &PrimeFactorization );

  mpz_t gcd;
  mpz_init_set_ui( gcd, PrimeFactorization.the_factors[0].occurrences );

  long i;
  for ( i = 1; i < PrimeFactorization.count && mpz_cmp_ui( gcd, 1 ) != 0; i++ )
    mpz_gcd_ui( gcd, gcd, PrimeFactorization.the_factors[i].occurrences );

  long power = mpz_get_d( gcd );

  mpz_set_ui( ConvertedBase, 1 );
  mpz_mul_ui( ConvertedExponent, exponent, power );

  mpz_t tempz1;
  mpz_init( tempz1 );

  mpz_set_ui( SqrFreeBase, 1 );
  for ( i = 0; i < PrimeFactorization.count; i++ ) {
    long reducedpower = PrimeFactorization.the_factors[i].occurrences / power;
    mpz_pow_ui( tempz1, PrimeFactorization.the_factors[i].the_factor, reducedpower );
    mpz_mul( ConvertedBase, ConvertedBase, tempz1 );

    if ( reducedpower % 2 )
      mpz_mul( SqrFreeBase, SqrFreeBase, PrimeFactorization.the_factors[i].the_factor );
  }

  mpz_clear( tempz1 );
  mpz_clear( gcd );
  Cleanup_Factor_Infos( &PrimeFactorization );
}

// Returns 1 if base^exponent + c0 has an Aurifeuillian factorization else 0.
// Assumption: base is not a perfect power and is square-free
int IsLM( mpz_t base, mpz_t exponent, char c0 ) {

  // Not handled values
  if ( mpz_cmp_ui( base, 2 ) < 0 ||
       mpz_cmp_ui( exponent, 2 ) < 0 ||
       (c0 != '-' && c0 != '+') ) {
    return 0;
  }

  {  // check that we have the correct c0.  c0 should be -1 if base == 1 mod 4, else c0 should be 1
    mpz_t r;
    mpz_init( r );
    mpz_mod_ui( r, base, 4 );
    if (  ( (c0 == '+') && mpz_cmp_ui( r, 1 ) == 0) ||
          ( (c0 == '-') && mpz_cmp_ui( r, 1 ) != 0) ) {
      mpz_clear( r );
      return 0;
    }
    mpz_clear( r );
  }

  if ( !mpz_divisible_p( exponent, base ) )  // not divisible.  L,M do not exist
    return 0;

  mpz_t h;
  mpz_init( h );
  mpz_divexact( h, exponent, base );
  if ( mpz_even_p( h ) ) { // h is even.  L,M do not exist
    mpz_clear( h );
    return 0;
  }

  mpz_clear( h );

  return 1;
}

// Compute C and D polynomial coefficients.  Code implemented based on
// Richard Brent's "On computing factors of cyclotomic polynomials" Algorithm L.
int ComputeLucasCD( mpz_t n, mpz_t** C, long* C_count, mpz_t** D, long* D_count ) {

  mpz_t tempz1;
  mpz_init( tempz1 );
  mpz_t tempz2;
  mpz_init( tempz2 );

  mpz_t nmod4;
  mpz_init( nmod4 );
  mpz_mod_ui( nmod4, n, 4 );

  mpz_t nmod8;
  mpz_init( nmod8 );
  mpz_mod_ui( nmod8, n, 8 );

  mpz_t nprime;
  mpz_init( nprime );
  if ( mpz_cmp_ui( nmod4, 1 ) == 0 )
    mpz_set( nprime, n );
  else
    mpz_mul_ui( nprime, n, 2 );

  mpz_t s;
  mpz_init( s );
  if ( mpz_cmp_ui( nmod4, 3 ) == 0 )
    mpz_set_si( s, -1 );
  else
    mpz_set_ui( s, 1 );

  mpz_t sprime;
  mpz_init( sprime );
  if ( mpz_cmp_ui( nmod8, 5 ) == 0 )
    mpz_set_si( sprime, -1 );
  else
    mpz_set_ui( sprime, 1 );

  mpz_t d;
  mpz_init( d );
  EulerTotient(nprime, tempz1 );
  mpz_divexact_ui( d, tempz1, 2 );

  // Step 1
  long ld = mpz_get_ui( d );
  mpz_t* q = (mpz_t*) calloc( ld + 1, sizeof( mpz_t ) );
  long k;
  for ( k = 0; k <= ld; k++ )
    mpz_init_set_ui( q[k], 0 );

  mpz_t mpz_k;
  mpz_init( mpz_k );

  for ( k = 1; k <= ld; k++ ) {
    mpz_set_ui( mpz_k, k );

    if ( k % 2 )
      mpz_set_si( q[k], mpz_jacobi( n, mpz_k ) );
    else {
      mpz_t gprime_k;
      mpz_init( gprime_k );
      mpz_gcd( gprime_k, mpz_k, nprime );

      mpz_t tempz1;
      mpz_init( tempz1 );
      mpz_divexact( tempz1, nprime, gprime_k );
      int mobius_value = Mobius( tempz1 );

      mpz_t mpz_totient_value;
      mpz_init( mpz_totient_value );
      EulerTotient( gprime_k, mpz_totient_value );

      int cosine_value = 0;
      {
        mpz_sub_ui( tempz1, n, 1 );
        mpz_mul( tempz1, tempz1, mpz_k );

        mpz_t mod8;
        mpz_init( mod8 );
        mpz_mod_ui( mod8, tempz1, 8 );

        // we won't check for odd values since k is even and therefore mod8 is also even
        if ( mpz_cmp_ui( mod8, 0 ) == 0 )
          cosine_value = 1;
        else if ( mpz_cmp_ui( mod8, 2 ) == 0 )
          cosine_value = 0;
        else if ( mpz_cmp_ui( mod8, 4 ) == 0 )
          cosine_value = -1;
        else // mod8 == 6
          cosine_value = 0;

        mpz_clear( mod8 );        
      }

      mpz_mul_si( tempz1, mpz_totient_value, mobius_value );
      mpz_mul_si( tempz1, tempz1, cosine_value );
      mpz_set( q[k], tempz1 );

      mpz_clear( mpz_totient_value );
      mpz_clear( tempz1 ); 
      mpz_clear( gprime_k );
    }
  }

  mpz_t* gamma = (mpz_t*) calloc( ld + 1, sizeof( mpz_t ) );
  for ( k = 0; k <= ld; k++ )
    mpz_init_set_ui( gamma[k], 0 );

  mpz_t* delta = (mpz_t*) calloc( ld, sizeof( mpz_t ) );
  for ( k = 0; k < ld; k++ )
    mpz_init_set_ui( delta[k], 0 );

  // Step 2
  mpz_set_ui( gamma[0], 1 );
  mpz_set_ui( delta[0], 1 );

  // Step 3
  int d_is_odd = ld % 2;
  long max_gamma = d_is_odd ? ( ld - 1 ) / 2 : ld / 2;
  long max_delta = d_is_odd ? ( ld - 1 ) / 2 : ( ld - 2 ) / 2;
  for ( k = 1; k <= max_gamma; k++ ) {
    // compute gamma[k]
    long j = 0;
    for ( j = 0; j <= k - 1; j++ ) {
      long q_index = 2 * k - 2 * j - 1;
      mpz_mul( tempz1, n, q[q_index] );
      mpz_mul( tempz1, tempz1, delta[j] );
      mpz_add( gamma[k], gamma[k], tempz1 );

      q_index++;
      mpz_mul( tempz1, q[q_index], gamma[j] );
      mpz_sub( gamma[k], gamma[k], tempz1 );
    }
    mpz_divexact_ui( gamma[k], gamma[k], 2 * k );

    if ( k > max_delta )
      continue;

    // compute delta[k]
    for ( j = 0; j <= k - 1; j++ ) {
      long q_index = 2 * k + 1 - 2 * j;
      mpz_mul( tempz1, q[q_index], gamma[j] );
      mpz_add( delta[k], delta[k], tempz1 );

      q_index--;
      mpz_mul( tempz1, q[q_index], delta[j] );
      mpz_sub( delta[k], delta[k], tempz1 );
    }
    mpz_add( delta[k], delta[k], gamma[k] );
    mpz_divexact_ui( delta[k], delta[k], 2 * k + 1 );
  }

  // Step 4
  for ( k = max_gamma + 1; k <= ld; k++ )
    mpz_set( gamma[k], gamma[ld-k] );

  // Step 5
  for ( k = max_delta + 1; k < ld; k++ )
    mpz_set( delta[k], delta[ld-1-k] );

  // Pass back
  *C = gamma;
  *C_count = ld + 1;
  *D = delta;
  *D_count = ld;

  // clean up
  mpz_clear( mpz_k );

  if ( q != NULL ) {
    for ( k = 0; k <= ld; k++ )
      mpz_clear( q[k] );
    free( q );
    q = NULL;
  }

  mpz_clear( d );
  mpz_clear( sprime );
  mpz_clear( s );
  mpz_clear( nprime );
  mpz_clear( nmod8 );
  mpz_clear( nmod4 );
  mpz_clear( tempz2 );
  mpz_clear( tempz1 );

  return 0;
}


// Computes the Aurifeuillian factors L,M for a^x +/- 1.  Returns 0 on success.
// Assumption: base is not a perfect power
int ComputeLM( mpz_t base, mpz_t sqrfreebase, mpz_t exponent, char c0, mpz_t L, mpz_t M ) {

  // Not handled values
  if ( mpz_cmp_ui( base, 2 ) < 0 ||
       mpz_cmp_ui( exponent, 2 ) < 0 ||
       (c0 != '-' && c0 != '+') ) {
    return -1;
  }

  {  // check that we have the correct c0.  c0 should be -1 if sqrfreebase == 1 mod 4, else c0 should be 1
    mpz_t r;
    mpz_init( r );
    mpz_mod_ui( r, sqrfreebase, 4 );
    if ( (c0 == '+' && mpz_cmp_ui( r, 1 ) == 0) ||
         (c0 == '-' && mpz_cmp_ui( r, 1 ) != 0) ) {
      mpz_clear( r );
      return -1;
    }
    mpz_clear( r );
  }

  if ( !mpz_divisible_p( exponent, sqrfreebase ) )  // not divisible.  L,M do not exist
    return -1;

  unsigned long h = 0;
  {
    mpz_t temph;
    mpz_init( temph );
    mpz_divexact( temph, exponent, sqrfreebase );
    h = mpz_get_ui( temph );
    mpz_clear( temph );
    if ( (h % 2) == 0 ) // h is even.  L,M do not exist
      return -1;
  }

  mpz_t* C = NULL;
  long C_count = 0;
  mpz_t* D = NULL;
  long D_count = 0;
  if ( ComputeLucasCD( sqrfreebase, &C, &C_count, &D, &D_count ) != 0 ) {
    gmp_printf("Error: Could not compute Lucas C,D polynomial coefficients for base (%Zd). Aborting.\n", sqrfreebase );
    return -1;
  }

  mpz_t tempz1;
  mpz_init( tempz1 );

  unsigned long k = (h + 1) / 2;

  // Both Polynomial A & B are initially computed as follows:
  //
  //   coef[0] * base^(h*0) + coef[1] * base^(h*1) + coef[2] * base^(h*2) + ...
  // = coef[0] * (base^h)^0 + coef[1] * (base^h)^1 + coef[2] * (base^h)^2 + ...
  // = coef[0] * 1 + coef[1] * base^h + coef[2] * base^h * base^h + ...
  //
  // B polynomial is then multiplied by additional terms
  //
  // Then,
  //   L = A - B
  //   M = A + B

  // NOTE: Reading coefficients from beginning to end is exactly the same as reading from the end to the beginning.
  // eg. for base 31, C = { 1,15,43,83,125,151,169,173,173,169,151,125,83,43,15,1 } and D = { 1,5,11,19,25,29,31,31,31,29,25,19,11,5,1 }
  // Also, D_count = C_count - 1


  // Compute Polynomial A with the Lucas 'C' coefficients
  mpz_t A;
  mpz_init( A );
  long i;
  {
    mpz_t base_h; // base^h
    mpz_init( base_h );
    mpz_pow_ui( base_h, base, h );

    mpz_t base_raised; // base_h^i
    mpz_init( base_raised );

    for ( i = 0; i < C_count; i++ ) { // each iteration computes one term and then adds to A
      if ( i == 0 )
        mpz_set_ui( base_raised, 1 );
      else if ( i == 1 )
        mpz_set( base_raised, base_h );
      else
        mpz_mul( base_raised, base_raised, base_h );

      mpz_mul( tempz1, base_raised, C[i] );
      mpz_add( A, A, tempz1 );
    }

    mpz_clear( base_raised );
    mpz_clear( base_h );
  }

  // Compute Polynomial B with Lucas 'D' coefficients
  mpz_t B;
  mpz_init( B );

  {
    mpz_t base_h; // base^h
    mpz_init( base_h );
    mpz_pow_ui( base_h, base, h );

    mpz_t base_raised; // base_h^i
    mpz_init( base_raised );

    for ( i = 0; i < D_count; i++ ) { // each iteration computes one term and then adds to B
      if ( i == 0 )
        mpz_set_ui( base_raised, 1 );
      else if ( i == 1 )
        mpz_set( base_raised, base_h );
      else
        mpz_mul( base_raised, base_raised, base_h );

      mpz_mul( tempz1, base_raised, D[i] );
      mpz_add( B, B, tempz1 );
    }

    mpz_clear( base_raised );
    mpz_clear( base_h );
  }

  if ( mpz_cmp( base, sqrfreebase ) != 0 ) {
    mpz_divexact( tempz1, base, sqrfreebase );
    mpz_sqrt( tempz1, tempz1 );

    mpz_pow_ui( tempz1, tempz1, h );
    mpz_mul( B, B, tempz1 );
  }

  mpz_pow_ui( tempz1, sqrfreebase, k );
  mpz_mul( B, B, tempz1 );

  mpz_sub( L, A, B );
  mpz_add( M, A, B );

  mpz_clear( B );
  mpz_clear( A );
  mpz_clear( tempz1 );

  for ( i = 0; i < C_count; i++ )
    mpz_clear( C[i] );
  free( C );
  C = NULL;

  for ( i = 0; i < D_count; i++ )
    mpz_clear( D[i] );
  free( D );
  D = NULL;

  return 0;
}


char quickprimecheck( mpz_t the_number ) {

  char retval = 'C';
  if ( mpz_cmp_ui( the_number, 1 ) == 0 ) {
    retval = 'N';
    return retval;
  }

  int factor_type = mpz_probab_prime_p( the_number, 30 );
  if ( factor_type != 0 )
    retval = 'P';

  return retval;
}


// divide out denominator from running_N
void TFDivideOut( mpz_t running_N, long ldenominator, char* running_N_status, mpz_t square_root, struct factor_infos* Factor_Infos ) {

  mpz_t denominator;
  mpz_init_set_ui( denominator, ldenominator );

  long occurrences = ComputeOccurrences( running_N, denominator );

  long i = 0;
  for ( i = occurrences; i > 0; i-- )
    mpz_divexact_ui( running_N, running_N, ldenominator );

  *running_N_status = quickprimecheck( running_N );
  if ( *running_N_status == 'N' )
    mpz_set_ui( square_root, 1 );
  else
    mpz_sqrt( square_root, running_N );

  AddFactorInfo( Factor_Infos, denominator, occurrences, 'P', NumberOfDigits( denominator ) );

  mpz_clear( denominator );
}


// Compute the prime factorization of a general number where square_root( the_number ) < MAXLONG - 255
void Factorize( mpz_t the_number, struct factor_infos* Factor_Infos ) {
  if ( Factor_Infos == NULL )
    return;

  // Wheel Factorization. eg. http://programmingpraxis.com/2009/05/08/wheel-factorization/
  mpz_t running_N;
  mpz_init_set( running_N, the_number );
  char running_N_status = quickprimecheck( running_N );
  if ( running_N_status == 'P' ) {
    AddFactorInfo( Factor_Infos, running_N, 1, 'P', NumberOfDigits( running_N ) );
    mpz_clear( running_N );
    return;
  }

  mpz_t square_root;
  mpz_init( square_root );
  mpz_sqrt( square_root, running_N );

  if ( mpz_cmp_ui( square_root, MAXLONG - 255 ) > 0 ) {
    printf( "Internal Error:  Trying to factorize too big a number.\n" );
    mpz_clear( square_root );
    mpz_clear( running_N );
    return;
  }

  if ( mpz_divisible_ui_p( running_N, 2 ) != 0 ) {
    TFDivideOut( running_N, 2, &running_N_status, square_root, Factor_Infos );
    if ( running_N_status == 'P' ) {
      AddFactorInfo( Factor_Infos, running_N, 1, 'P', NumberOfDigits( running_N ) );
      mpz_clear( square_root );
      mpz_clear( running_N );
      return;
    }
  }

  if ( mpz_divisible_ui_p( running_N, 3 ) != 0 ) {
    TFDivideOut( running_N, 3, &running_N_status, square_root, Factor_Infos );
    if ( running_N_status == 'P' ) {
      AddFactorInfo( Factor_Infos, running_N, 1, 'P', NumberOfDigits( running_N ) );
      mpz_clear( square_root );
      mpz_clear( running_N );
      return;
    }
  }

  if ( mpz_divisible_ui_p( running_N, 5 ) != 0 ) {
    TFDivideOut( running_N, 5, &running_N_status, square_root, Factor_Infos );
    if ( running_N_status == 'P' ) {
      AddFactorInfo( Factor_Infos, running_N, 1, 'P', NumberOfDigits( running_N ) );
      mpz_clear( square_root );
      mpz_clear( running_N );
      return;
    }
  }

  if ( mpz_divisible_ui_p( running_N, 7 ) != 0 ) {
    TFDivideOut( running_N, 7, &running_N_status, square_root, Factor_Infos );
    if ( running_N_status == 'P' ) {
      AddFactorInfo( Factor_Infos, running_N, 1, 'P', NumberOfDigits( running_N ) );
      mpz_clear( square_root );
      mpz_clear( running_N );
      return;
    }
  }

  unsigned long i = 0;
  unsigned long denominator = 11;
  for ( ;
        running_N_status == 'C';
        denominator += wheel7[i], i = ( i == 47 ? 0 : i + 1 ) ) {

      if ( mpz_divisible_ui_p( running_N, denominator ) != 0 )
        TFDivideOut( running_N, denominator, &running_N_status, square_root, Factor_Infos );
  }

  if ( running_N_status == 'P' )
    AddFactorInfo( Factor_Infos, running_N, 1, 'P', NumberOfDigits( running_N ) );

  mpz_clear( square_root );
  mpz_clear( running_N );
}


void PermuteFactors( struct factor_infos* permutefactors, struct factor_infos* allfactors, long lefttumbler, mpz_t leftproduct, mpz_t the_number ) {

  long tumbler = lefttumbler + 1;

  mpz_t currfactor;
  mpz_init( currfactor );

  mpz_t currproduct;
  mpz_init( currproduct );

  long exponent;
  for ( exponent = 0; exponent <= permutefactors->the_factors[tumbler].occurrences; exponent++ ) {
    // compute primefactor^exponent
    if ( exponent == 0 )
      mpz_set_ui( currfactor, 1 );
    else if ( exponent == 1 )
      mpz_set( currfactor, permutefactors->the_factors[tumbler].the_factor );
    else
      mpz_mul( currfactor, currfactor, permutefactors->the_factors[tumbler].the_factor );

    mpz_mul( currproduct, leftproduct, currfactor );

    if ( tumbler < permutefactors->count - 1 )
      PermuteFactors( permutefactors, allfactors, tumbler, currproduct, the_number );
    else if ( mpz_cmp_ui( currproduct, 1 ) != 0 && mpz_cmp( currproduct, the_number ) != 0 ) {
      char status = quickprimecheck( currproduct );
      long occurrences = ComputeOccurrences( the_number, currproduct );
      AddFactorInfo( allfactors, currproduct, occurrences, status, NumberOfDigits( currproduct ) );
    }
  }

  mpz_clear( currproduct );
  mpz_clear( currfactor );
}


// Compute all factors of the_number except for 1 and the number itself
void ComputeAllFactors( mpz_t the_number, struct factor_infos* Factor_Infos ) {
  if ( Factor_Infos == NULL )
    return;

  if ( quickprimecheck( the_number ) != 'C' )
    return;

  struct factor_infos primes_only;
  Init_Factor_Infos( &primes_only );

  Factorize( the_number, &primes_only );

  mpz_t leftproduct;
  mpz_init_set_ui( leftproduct, 1 );
  PermuteFactors( &primes_only, Factor_Infos, -1, leftproduct, the_number );
  mpz_clear( leftproduct );

  Cleanup_Factor_Infos( &primes_only );

  qsort( Factor_Infos->the_factors, Factor_Infos->count, sizeof(struct factor_info), factor_info_cmpfunc );
}


// Compute the number of times denominator divides evenly into numerator
long ComputeOccurrences( mpz_t numerator, mpz_t denominator ) {

  // not handled values
  if ( mpz_cmp_ui( denominator, 1 ) <= 0 )
    return -1;

  long occurrences = 0;

  mpz_t quotient;
  mpz_init_set( quotient, numerator );

  while ( mpz_divisible_p( quotient, denominator ) ) {
    mpz_divexact( quotient, quotient, denominator );
    occurrences++;
  }

  mpz_clear( quotient );

  return occurrences;
}


// Compute the Mobius function
int Mobius( mpz_t n ) {

  if ( mpz_cmp_ui( n, 1 ) == 0 )
    return 1;

  long total_occurrences = 0;
  int has_duplicate_occurrences = 0;

  if ( quickprimecheck( n ) == 'P' ) {
    total_occurrences++;
  }
  else {
    struct factor_infos Factor_Infos;
    Init_Factor_Infos( &Factor_Infos );

    ComputeAllFactors( n, &Factor_Infos );

    long i = 0;
    for ( ; i < Factor_Infos.count; i++ ) {
      if ( Factor_Infos.the_factors[i].factor_status == 'P' ) {
        total_occurrences += Factor_Infos.the_factors[i].occurrences;
        if ( Factor_Infos.the_factors[i].occurrences > 1 )
          has_duplicate_occurrences = 1;
      }
    }

    Cleanup_Factor_Infos( &Factor_Infos );
  }

  if ( has_duplicate_occurrences )
    return 0;

  if ( total_occurrences % 2 )
    return -1;

  return 1;
}


// Compute Euler's Totient function
// tot(mn) = tot(m)*tot(n)*d/tot(d)   where d = gcd(m,n)
// see: http://en.wikipedia.org/wiki/Euler's_totient_function
void EulerTotient( mpz_t the_number, mpz_t tot ) {

  if ( mpz_cmp_si( the_number, 0 ) < 0 ) {
    mpz_set_si( tot, -1 );
    return;
  }

  mpz_set_ui( tot, 1 );
  // by convention, tot(0) is 1
  if ( mpz_cmp_si( the_number, 0 ) == 0 || mpz_cmp_ui( the_number, 1 ) == 0 )
    return;

  if ( quickprimecheck( the_number ) == 'P' ) {
    mpz_sub_ui( tot, the_number, 1 );
    return;
  }

  struct factor_infos primefactorization;
  Init_Factor_Infos( &primefactorization );
  Factorize( the_number, &primefactorization );

  mpz_t tempZ;
  mpz_init( tempZ );

  long i,j;
  // tot(n^k) = n^(k-1) * tot(n) see: http://mathworld.wolfram.com/TotientFunction.html  (14)
  // and if n is prime, then this is just: tot(p^k) = p^(k-1) * (p - 1)
  for ( i = 0; i < primefactorization.count; i++ ) {
    for ( j = 1; j < primefactorization.the_factors[i].occurrences; j++ )
      mpz_mul( tot, tot, primefactorization.the_factors[i].the_factor );

    mpz_sub_ui( tempZ, primefactorization.the_factors[i].the_factor, 1 );
    mpz_mul( tot, tot, tempZ );
  }

  mpz_clear( tempZ );

  Cleanup_Factor_Infos( &primefactorization );
}

// Initialize factor_infos
void Init_Factor_Infos( struct factor_infos* Factor_Infos ) {
  if ( Factor_Infos == NULL )
    return;
  Factor_Infos->count = 0;
  Factor_Infos->the_factors = NULL;
}

// Add a factor to factor_infos structure
void AddFactorInfo( struct factor_infos* Factor_Infos, mpz_t the_factor, long occurrences, char factor_status, long numdigits ) {

  long index = 0;

  // allocate memory
  if ( Factor_Infos->count == 0 ) {
    Factor_Infos->count = 1;
    Factor_Infos->the_factors = (struct factor_info*) calloc( 1, sizeof(struct factor_info) );
  }
  else {
    Factor_Infos->count++;
    Factor_Infos->the_factors = (struct factor_info*) realloc( Factor_Infos->the_factors, sizeof(struct factor_info) * Factor_Infos->count );
    index = Factor_Infos->count - 1;
    memset( &Factor_Infos->the_factors[index], 0, sizeof(struct factor_info) );
  }

  mpz_init( Factor_Infos->the_factors[index].the_factor );
  mpz_set( Factor_Infos->the_factors[index].the_factor, the_factor );
  Factor_Infos->the_factors[index].digits = numdigits;
  Factor_Infos->the_factors[index].occurrences    = occurrences;
  Factor_Infos->the_factors[index].factor_status  = factor_status;
}

// Free the memory allocated
void Cleanup_Factor_Infos( struct factor_infos* Factor_Infos ) {
  if ( Factor_Infos == NULL )
    return;

  long i;
  for ( i = 0; i < Factor_Infos->count; i++ )
    mpz_clear( Factor_Infos->the_factors[i].the_factor );

  if ( Factor_Infos->the_factors != NULL ) {
    free( Factor_Infos->the_factors );
    Factor_Infos->the_factors = NULL;
  }
  Factor_Infos->count = 0;
}


// Return the number of digits (base 10) of "the_number".  Positive numbers only.
long NumberOfDigits( mpz_t the_number ) {

  if ( mpz_cmp_si( the_number, 0 ) < 0 )
    return -1;

  char* p = NULL;
  gmp_asprintf( &p, "%Zd", the_number );

  long len = strlen( p );

  if ( p != NULL ) {
    free( p );
    p = NULL;
  }

  return len;
}

// Just call mpz_pow_ui() with some error checking first
void mpz_pow( mpz_t result, mpz_t b, mpz_t p ) {

  if ( mpz_cmp_ui( p, ULONG_MAX ) > 0 ) {
    printf( "Error: exponent too large in power function.\n" );
    return;
  }

  mpz_pow_ui( result, b, mpz_get_ui( p ) );
}


int factor_info_cmpfunc( const void* p1, const void* p2 ) {
  return mpz_cmp( (*((struct factor_info*)p1)).the_factor, (*((struct factor_info*)p2)).the_factor );
}

