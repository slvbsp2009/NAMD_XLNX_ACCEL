/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
 * Copyright (c) 1993 Martin Birgmeier
 * All rights reserved.
 *
 * You may redistribute unmodified or modified versions of this source
 * code provided that the above copyright notice and this and the
 * following conditions are retained.
 *
 * This software is provided ``as is'', and comes with no warranties
 * of any kind. I shall in no event be liable for anything that happens
 * to anyone/anything when using this software.
 */

#ifndef RANDOM_H
#define RANDOM_H

#include "common.h"
#include "Vector.h"

#ifdef _MSC_VER
#define UINT64_LITERAL(X) X ## ui64
#else
#define UINT64_LITERAL(X) X ## ULL
#endif

#define	RAND48_SEED   UINT64_LITERAL(0x00001234abcd330e)
#define	RAND48_MULT   UINT64_LITERAL(0x00000005deece66d)
#define	RAND48_ADD    UINT64_LITERAL(0x000000000000000b)
#define RAND48_MASK   UINT64_LITERAL(0x0000ffffffffffff)

class Random {

private:

  double second_gaussian;
  uint64_t second_gaussian_waiting;
  uint64_t rand48_seed;
  uint64_t rand48_mult;
  uint64_t rand48_add;

public:

  // default constructor
  Random(void) {
    init(0);
    rand48_seed = RAND48_SEED;
  }

  // constructor with seed
  Random(unsigned long seed) {
    init(seed);
  }

  // reinitialize with seed
  void init(unsigned long seed) {
    second_gaussian = 0;
    second_gaussian_waiting = 0;
    rand48_seed = seed & UINT64_LITERAL(0x00000000ffffffff);
    rand48_seed = rand48_seed << 16;
    rand48_seed |= RAND48_SEED & UINT64_LITERAL(0x0000ffff);
    rand48_mult = RAND48_MULT;
    rand48_add = RAND48_ADD;
  }

  // advance generator by one (seed = seed * mult + add, to 48 bits)
  void skip(void) {
    rand48_seed = ( rand48_seed * rand48_mult + rand48_add ) & RAND48_MASK;
  }

  // split into numStreams different steams and take stream iStream
  void split(int iStream, int numStreams) {

    int i;

    // make sure that numStreams is odd to ensure maximum period
    numStreams |= 1;

    // iterate to get to the correct stream
    for ( i = 0; i < iStream; ++i ) skip();

    // save seed and add so we can use skip() for our calculations
    uint64_t save_seed = rand48_seed;

    // calculate c *= ( 1 + a + ... + a^(numStreams-1) )
    rand48_seed = rand48_add;
    for ( i = 1; i < numStreams; ++i ) skip();
    uint64_t new_add = rand48_seed;

    // calculate a = a^numStreams
    rand48_seed = rand48_mult;
    rand48_add  = 0;
    for ( i = 1; i < numStreams; ++i ) skip();
    rand48_mult = rand48_seed;

    rand48_add  = new_add;
    rand48_seed = save_seed;

    second_gaussian = 0;
    second_gaussian_waiting = 0;
  }

  // return a number uniformly distributed between 0 and 1
  BigReal uniform(void) {
    skip();
    const double exp48 = ( 1.0 / (double)(UINT64_LITERAL(1) << 48) );
    return ( (double) rand48_seed * exp48 );
  }

  // return a number from a standard gaussian distribution
  BigReal gaussian(void) {
    BigReal fac, r, v1, v2;

    if (second_gaussian_waiting) {
      second_gaussian_waiting = 0;
      return second_gaussian;
    } else {
      r = 2.;                 // r >= 1.523e-8 ensures abs result < 6
      while (r >=1. || r < 1.523e-8) { // make sure we are within unit circle
        v1 = 2.0 * uniform() - 1.0;
        v2 = 2.0 * uniform() - 1.0;
        r = v1*v1 + v2*v2;
      }
      fac = sqrt(-2.0 * log(r)/r);
      // now make the Box-Muller transformation to get two normally
      // distributed random numbers. Save one and return the other.
      second_gaussian_waiting = 1;
      second_gaussian = v1 * fac;
      return v2 * fac;
    }
  }

 /** \brief Return the sum of i.i.d. squared Gaussians.
  *
  *  The sum of k squared standard normal random variables is equivalent to
  *  drawing a value from a chi-squared distribution with the shape given as
  *  the number of Gaussians. That is,
  *
  *  X ~ sum_n=1^k N(0, 1)^2 ~ chi^2(k).
  *
  *  This is in turn a special case of the Gamma distribution with shape k/2
  *  and scale 2 (we could also use rate = 1/scale)
  *
  *  X ~ chi^2(k) ~ Gamma(k/2, 2) ~ 2*Gamma(k/2, 1).
  *
  *  The second relation follows from the scaling property of the Gamma
  *  distribution. Furthermore, when a Gamma distribution has unit scale and a
  *  large shape, it can be well-approximated by a normal distribution with
  *  mean and variance equal to the shape. Thus,
  *
  *  X ~ 2*Gamma(k/2, 1) ~= 2*N(k/2, sqrt(k/2)).
  *
  *  A quick numerical test shows that the mean integrated square error for
  *  this approximation is <10^-5 for shape >30 and <10^-6 for shape >100.
  *  We'll be conservative and use the latter cutoff.
  *
  *  We thus have three cases for k Gaussians:
  *
  *   0 < k <=   2 - just brute force generate and sum the Gaussians
  *   2 < k <= 200 - use a (slightly modified) version of the Gamma
  *                  distribution algorithm from Marsaglia and Tsang that is
  *                  implemented in the GNU Science Library (GSL)
  *   else         - use a single Gaussian distribution
  *
  *   The brute force method is almost certainly the slowest method, even for
  *   k = 3. The rigorous method takes about 150% as long as the approximate
  *   method (but this is invariant to the value of k).
  *
  *  \param num_gaussians A positive integer number of Gaussians
  *  \return a random variable equal to the sum of num_gaussians squared
  *    standard normal variables
  */
  BigReal sum_of_squared_gaussians(int64_t num_gaussians) {
    BigReal z, u, v;

    if (num_gaussians <= 2) {
        v = 0.0;
        for(int i=0; i<num_gaussians; ++i) {
          z = gaussian();
          v += z*z;
        }
        return v;
    } else if (2 < num_gaussians && num_gaussians <= 200) {
      const BigReal d = 0.5*num_gaussians - 1./3;
      const BigReal c = 1. / (3*sqrt(d));
      const BigReal zmin = -1. / c;
      do {
        do {
          z = gaussian();
        }
        while(z <= zmin);
        u = uniform();
        v = (1 + c*z); v = v*v*v; // v = (1 + c*z)^3
      }
      while(log(u) >= (0.5*z*z + d*(1 - v + log(v))));
      return 2*d*v;
    } else {
      return num_gaussians + sqrt(2*num_gaussians)*gaussian();
    }
  }

  // return a vector of gaussian random numbers
  Vector gaussian_vector(void) {
    return Vector( gaussian(), gaussian(), gaussian() );
  }

  // return a random long
  // signed int32 ranges from 0 to (2^31)-1
  long integer(void) {
    skip();
    return ( ( rand48_seed >> 17 ) & UINT64_LITERAL(0x000000007fffffff) );
  }

  // randomly order an array of whatever
  template <class Elem> void reorder(Elem *a, int n) {
    for ( int i = 0; i < (n-1); ++i ) {
      int ie = i + ( integer() % (n-i) );
      if ( ie == i ) continue;
      const Elem e = a[ie];
      a[ie] = a[i];
      a[i] = e;
    }
  }

  ///
  /// Fill length n array from uniform distribution.
  ///
  void uniform_array_f(float *a, int n) {
    // for now just call uniform()
    // ultimately we will provide faster implementation
    for (int i=0;  i < n;  i++) {
      a[i] = (float)uniform();
    }
  }

  ///
  /// Fill length n array from standard Gaussian distribution.
  ///
  void gaussian_array_f(float *a, int n) {
    // for now just call gaussian()
    // ultimately we will provide faster implementation
    for (int i=0;  i < n;  i++) {
      a[i] = (float)gaussian();
    }
  }

};

#endif  // RANDOM_H

