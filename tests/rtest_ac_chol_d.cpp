/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 3.5                                                 *
 *                                                                        *
 *  Release Date    : Thu Feb  8 17:36:42 PST 2024                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.5.0                                               *
 *                                                                        *
 *  Copyright 2018 Siemens                                                *
 *                                                                        *
 **************************************************************************
 *  Licensed under the Apache License, Version 2.0 (the "License");       *
 *  you may not use this file except in compliance with the License.      * 
 *  You may obtain a copy of the License at                               *
 *                                                                        *
 *      http://www.apache.org/licenses/LICENSE-2.0                        *
 *                                                                        *
 *  Unless required by applicable law or agreed to in writing, software   * 
 *  distributed under the License is distributed on an "AS IS" BASIS,     * 
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or       *
 *  implied.                                                              * 
 *  See the License for the specific language governing permissions and   * 
 *  limitations under the License.                                        *
 **************************************************************************
 *                                                                        *
 *  The most recent version of this package is available at github.       *
 *                                                                        *
 *************************************************************************/
// =========================TESTBENCH=======================================
// This testbench file contains a stand-alone testbench that exercises the
// ac_chol_d() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_chol_d.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_chol_d.h>
using namespace ac_math;

// ==============================================================================
// Test Designs
//   These simple functions allow executing the ac_chol_d() function
//   using multiple data types at the same time. Template parameters are
//   used to configure the bit-widths of the types.

// Test Design for real and complex fixed point values.
template <unsigned M, bool use_pwl, int Wfi, int Ifi, int outWfi, int outIfi, bool outSfi>
void test_ac_chol_d_fixed(
  const ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP> A1[M][M],
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> L1[M][M],
  const ac_complex<ac_fixed<Wfi + 1, Ifi + 1, true, AC_TRN, AC_WRAP> > A2[M][M],
  ac_complex<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> > L2[M][M],
  const ac_matrix<ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP>, M, M> &A3,
  ac_matrix<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP>, M, M> &L3,
  const ac_matrix<ac_complex<ac_fixed<Wfi + 1, Ifi + 1, true, AC_TRN, AC_WRAP> >, M, M> &A4,
  ac_matrix<ac_complex<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> >, M, M> &L4
)
{
  #ifdef _WIN32 
    ac_chol_d<M, use_pwl>(A1, L1);
    ac_chol_d<M, use_pwl>(A2, L2);
  #else
    ac_chol_d<use_pwl>(A1, L1);
    ac_chol_d<use_pwl>(A2, L2);
  #endif
  ac_chol_d<use_pwl>(A3, L3);
  ac_chol_d<use_pwl>(A4, L4);
}

// Test Design for real ac_float values.
template <unsigned M, bool use_pwl,int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
void test_ac_chol_d_float(
  const ac_float<Wfl, Ifl, Efl, AC_TRN> A1[M][M],
  ac_float<outWfl, outIfl, outEfl, AC_TRN> L1[M][M],
  const ac_matrix<ac_float<Wfl, Ifl, Efl, AC_TRN>, M, M> &A2,
  ac_matrix<ac_float<outWfl, outIfl, outEfl, AC_TRN>, M, M> &L2
)
{
  #ifdef _WIN32 
  ac_chol_d<M, use_pwl>(A1, L1);
  #else
  ac_chol_d<use_pwl>(A1, L1);
  #endif
  ac_chol_d<use_pwl>(A2, L2);
}

// Test Design for real ac_std_float values.
template <unsigned M, bool use_pwl,int Wstfl, int Estfl, int outWstfl, int outEstfl>
void test_ac_chol_d_stfloat(
  const ac_std_float<Wstfl, Estfl> A1[M][M],
  ac_std_float<outWstfl, outEstfl> L1[M][M],
  const ac_matrix<ac_std_float<Wstfl, Estfl>, M, M> &A2,
  ac_matrix<ac_std_float<outWstfl, outEstfl>, M, M> &L2
)
{
  #ifdef _WIN32 
  ac_chol_d<M, use_pwl>(A1, L1);
  #else
  ac_chol_d<use_pwl>(A1, L1);
  #endif
  ac_chol_d<use_pwl>(A2, L2);
}

// Test Design for real ac_ieee_float values.
template <unsigned M, bool use_pwl,ac_ieee_float_format in_format, ac_ieee_float_format out_format>
void test_ac_chol_d_ifloat(
  const ac_ieee_float<in_format> A1[M][M],
  ac_ieee_float<out_format> L1[M][M],
  const ac_matrix<ac_ieee_float<in_format>, M, M> &A2,
  ac_matrix<ac_ieee_float<out_format>, M, M> &L2
)
{
  #ifdef _WIN32 
  ac_chol_d<M, use_pwl>(A1, L1);
  #else
  ac_chol_d<use_pwl>(A1, L1);
  #endif
  ac_chol_d<use_pwl>(A2, L2);
}

// ==============================================================================

#include <ac_math/ac_normalize.h>
using namespace ac_math;
#include <math.h>
#include <string>
#include <fstream>
#include <limits>
#include <random>
#include <iostream>
using namespace std;

// ------------------------------------------------------------------------------
// Helper functions

// Helper structs for printing out type info for ac_std_float and ac_ieee_float

// Generic struct, enables template specialization
template <typename T>
struct type_string_st { };

// Specialized struct, handles ac_std_floats
template <int W, int E>
struct type_string_st<ac_std_float<W, E> > {
  static string type_string() {
    string format_string = "ac_std_float<";
    format_string += ac_int<32,true>(W).to_string(AC_DEC);
    format_string += ",";
    format_string += ac_int<32,true>(E).to_string(AC_DEC);
    format_string += ">";

    return format_string;
  }
};

// Specialized struct, handles ac_ieee_floats
template <ac_ieee_float_format Format>
struct type_string_st<ac_ieee_float<Format> > {
  static string type_string() {
    string format_string = "ac_ieee_float<";
    if (Format == binary16)  { format_string += "binary16"; }
    if (Format == binary32)  { format_string += "binary32"; }
    if (Format == binary64)  { format_string += "binary64"; }
    if (Format == binary128) { format_string += "binary128"; }
    if (Format == binary256) { format_string += "binary256"; }
    format_string += ">";

    return format_string;
  }
};

#ifdef DEBUG
// print_matrix functions: Print 2D C-style matrix for debugging purposes.

// print_matrix for ac_fixed/ac_complex<ac_fixed>/ac_float matrix.
template<unsigned M, class T>
void print_matrix(const T mat[M][M])
{
  cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl;
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {cout << "mat[" << i << "][" << j << "] = " << mat[i][j] << endl;}
  }
}

// print_matrix for ac_std_float matrix.
template<unsigned M, int W, int E>
void print_matrix(const ac_std_float<W, E> mat[M][M])
{
  cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl;
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {cout << "mat[" << i << "][" << j << "] = " << mat[i][j].to_ac_float().to_double() << endl;}
  }
}

// print_matrix for ac_ieee_float matrix.
template<unsigned M, ac_ieee_float_format Format>
void print_matrix(const ac_ieee_float<Format> mat[M][M])
{
  cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl;
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {cout << "mat[" << i << "][" << j << "] = " << mat[i][j].to_ac_float().to_double() << endl;}
  }
}
#endif

// Generate positive definite matrix of ac_fixed values
template<unsigned M, int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
void gen_matrix(ac_fixed<W, I, S, Q, O> A[M][M])
{
  static_assert(I - int(S) >= ac::nbits<M - 1>::val, "Not enough integer bits in input type.");
  static_assert(W - I >= 2, "Input type must have at least 2 fractional bits.");
  
  // Declare two MxM matrices, once of which (tbmatT) is the transpose of the other.
  ac_fixed<(W - I)/2, 0, false> tbmat[M][M];
  ac_fixed<(W - I)/2, 0, false> tbmatT[M][M];
  
  // Make sure the minimum limit for the random number generator is the quantum double value.
  double min_rand_limit = std::numeric_limits<double>::min();
  // default_random_engine and uniform_real_distribution are libraries from the <random> header,
  // packaged with C++11 and later standards.
  default_random_engine generator;
  // Use a uniform distribution to maximize the chance of obtaining an invertible tbmat/tbmatT matrix.
  // The output produced by this is similar to that of matlab's "rand" function, i.e. a randomly
  // selected value picked out of a uniform probability density in the range of (0, 1).
  uniform_real_distribution<double> distribution(min_rand_limit, 1.0);

  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      tbmat[i][j] = distribution(generator);
      tbmatT[j][i] = tbmat[i][j];
    }
  }

  #ifdef DEBUG
  cout << "tbmat is : " << endl;
  print_matrix<M>(tbmat);
  cout << "tbmatT is : " << endl;
  print_matrix<M>(tbmatT);
  #endif

  // Multiply tbmat by its transpose to get the positive definite input matrix
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      A[i][j] = 0;
      for (int k = 0; k < (int)M; k++) {
        A[i][j] += tbmat[i][k] * tbmatT[k][j];
      }
    }
  }
  
  #ifdef DEBUG
  cout << "A in gen_matrix function is : " << endl;
  print_matrix<M>(A);
  #endif
}

// Generate positive definite matrix of ac_complex<ac_fixed> values
template<unsigned M, int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
void gen_matrix(ac_complex<ac_fixed<W, I, S, Q, O> > A[M][M])
{
  static_assert(S, "Input type must be signed");
  static_assert(I - 1 >= ac::nbits<M - 1>::val, "Not enough integer bits in input type.");
  static_assert(W - I >= 4, "Input type must have at least 4 fractional bits.");
  
  // Declare two MxM matrices, once of which (tbmatT) is the conjugate transpose of the other.
  ac_complex<ac_fixed<(W - I)/2, 0, true> > tbmat[M][M];
  ac_complex<ac_fixed<(W - I)/2, 0, true> > tbmatT[M][M];
  
  // Make sure the minimum limit for the random number generator is the quantum double value.
  double min_rand_limit = std::numeric_limits<double>::min();
  // default_random_engine and uniform_real_distribution are libraries from the <random> header,
  // packaged with C++11 and later standards.
  default_random_engine generator;
  // Use a uniform distribution to maximize the chance of obtaining an invertible tbmat/tbmatT matrix.
  uniform_real_distribution<double> distribution(min_rand_limit, 0.5);

  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      ac_fixed<(W - I)/2 - 1, -1, false> rand_val;
      tbmat[i][j].r() = distribution(generator);
      tbmat[i][j].i() = distribution(generator);
      tbmatT[j][i] = tbmat[i][j].conj();
    }
  }

  #ifdef DEBUG
  cout << "tbmat is : " << endl;
  print_matrix<M>(tbmat);
  cout << "tbmatT is : " << endl;
  print_matrix<M>(tbmatT);
  #endif

  // Multiply tbmat by its transpose to get the positive definite input matrix
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      A[i][j] = 0;
      for (int k = 0; k < (int)M; k++) {
        A[i][j] += tbmat[i][k] * tbmatT[k][j];
      }
    }
  }

  #ifdef DEBUG
  cout << "A in gen_matrix function is : " << endl;
  print_matrix<M>(A);
  #endif
}

// Generate positive definite matrix of ac_float values
template<unsigned M, int W, int I, int E, ac_q_mode Q>
void gen_matrix(ac_float<W, I, E, Q> A[M][M]) {
  static_assert(I >= 1, "Input mantissa must have at least 1 integer bit.");
  static_assert(W >= ac::nbits<M - 1>::val + 2 + 1, "Not enough bitwidth for mantissa.");
  
  enum {
    max_int_val = ac::nbits<M - 1>::val,
    E_val_1 = ac::nbits<W - 1>::val,
    val_2 = (max_int_val - (I - 1) < 0) ? -(max_int_val - (I - 1)) : max_int_val - (I - 1),
    E_val_2 = ac::nbits<AC_MAX(val_2, 1)>::val,
    E_min = AC_MAX(int(E_val_1), int(E_val_2)) + 1,
  };
  
  static_assert(E >= E_min, "Not enough exponent bits.");
  
  // Make sure the minimum limit for the random number generator is the quantum double value.
  double min_rand_limit = std::numeric_limits<double>::min();
  // default_random_engine and uniform_real_distribution are libraries from the <random> header,
  // packaged with C++11 and later standards.
  default_random_engine generator;
  // Use a uniform distribution to maximize the chance of obtaining an invertible tbmat/tbmatT matrix.
  // The output produced by this is similar to that of matlab's "rand" function, i.e. a randomly
  // selected value picked out of a uniform probability density in the range of (0, 1).
  uniform_real_distribution<double> distribution(min_rand_limit, 1.0);
  
  ac_fixed<(W - max_int_val - 1)/2, 0, false> tbmat[M][M], tbmatT[M][M];
  
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      tbmat[i][j] = distribution(generator);
      tbmatT[j][i] = tbmat[i][j];
    }
  }
  
  ac_fixed<W - 1, max_int_val, false> mac_val;

  // Multiply tbmat by its transpose to get the positive definite input matrix
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      mac_val = 0;
      for (int k = 0; k < (int)M; k++) {
        mac_val += tbmat[i][k] * tbmatT[k][j];
      }
      A[i][j] = mac_val;
    }
  }

  #ifdef DEBUG
  cout << "tbmat:" << endl;
  print_matrix<M>(tbmat);
  cout << "tbmatT:" << endl;
  print_matrix<M>(tbmatT);
  cout << "A:" << endl;
  print_matrix<M>(A);
  #endif
}

// Generate positive definite matrix of ac_std_float values
template<unsigned M, int W, int E>
void gen_matrix(ac_std_float<W, E> A[M][M])
{
  static_assert(W - E - 1 >= ac::nbits<M - 1>::val + 2 + 1, "Not enough bitwidth for mantissa.");
  
  enum {
    max_int_val = ac::nbits<M - 1>::val,
    E_val_1 = ac::nbits<W - E - 1>::val,
    E_val_2 = ac::nbits<AC_MAX(max_int_val - 1, 1)>::val,
    E_min = AC_MAX(int(E_val_1), int(E_val_2)) + 1,
  };
  
  static_assert(E >= E_min, "Not enough exponent bits.");
  
  // Make sure the minimum limit for the random number generator is the quantum double value.
  double min_rand_limit = std::numeric_limits<double>::min();
  // default_random_engine and uniform_real_distribution are libraries from the <random> header,
  // packaged with C++11 and later standards.
  default_random_engine generator;
  // Use a uniform distribution to maximize the chance of obtaining an invertible tbmat/tbmatT matrix.
  // The output produced by this is similar to that of matlab's "rand" function, i.e. a randomly
  // selected value picked out of a uniform probability density in the range of (0, 1).
  uniform_real_distribution<double> distribution(min_rand_limit, 1.0);
  
  ac_fixed<(W - E - max_int_val - 1)/2, 0, false> tbmat[M][M], tbmatT[M][M];
  
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      tbmat[i][j] = distribution(generator);
      tbmatT[j][i] = tbmat[i][j];
    }
  }
  
  ac_fixed<W - E - 1, max_int_val, false> mac_val;

  // Multiply tbmat by its transpose to get the positive definite input matrix
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      mac_val = 0;
      for (int k = 0; k < (int)M; k++) {
        mac_val += tbmat[i][k] * tbmatT[k][j];
      }
      A[i][j] = ac_std_float<W, E>(mac_val);
    }
  }

  #ifdef DEBUG
  cout << "tbmat:" << endl;
  print_matrix<M>(tbmat);
  cout << "tbmatT:" << endl;
  print_matrix<M>(tbmatT);
  cout << "A_ac_fl:" << endl;
  print_matrix<M>(A);
  #endif
}

// Generate positive definite matrix of ac_ieee_float values
template<unsigned M, ac_ieee_float_format Format>
void gen_matrix(ac_ieee_float<Format> A[M][M])
{
  typedef ac_ieee_float<Format> T_in;
  
  enum {
    W = T_in::width,
    E = T_in::e_width,
    max_int_val = ac::nbits<M - 1>::val,
  };
  
  // Make sure the minimum limit for the random number generator is the quantum double value.
  double min_rand_limit = std::numeric_limits<double>::min();
  // default_random_engine and uniform_real_distribution are libraries from the <random> header,
  // packaged with C++11 and later standards.
  default_random_engine generator;
  // Use a uniform distribution to maximize the chance of obtaining an invertible tbmat/tbmatT matrix.
  // The output produced by this is similar to that of matlab's "rand" function, i.e. a randomly
  // selected value picked out of a uniform probability density in the range of (0, 1).
  uniform_real_distribution<double> distribution(min_rand_limit, 1.0);
  
  ac_fixed<(W - E - max_int_val - 1)/2, 0, false> tbmat[M][M], tbmatT[M][M];
  
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      tbmat[i][j] = distribution(generator);
      tbmatT[j][i] = tbmat[i][j];
    }
  }
  
  ac_fixed<W - E - 1, max_int_val, false> mac_val;

  // Multiply tbmat by its transpose to get the positive definite input matrix
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      mac_val = 0;
      for (int k = 0; k < (int)M; k++) {
        mac_val += tbmat[i][k] * tbmatT[k][j];
      }
      A[i][j] = ac_ieee_float<Format>(mac_val);
    }
  }

  #ifdef DEBUG
  cout << "tbmat:" << endl;
  print_matrix<M>(tbmat);
  cout << "tbmatT:" << endl;
  print_matrix<M>(tbmatT);
  cout << "A_ac_fl:" << endl;
  print_matrix<M>(A);
  #endif
}

// Testbench for cholesky decomposition for ac_fixed matrix
// The testbench uses the cholesky-crout algorithm
template<unsigned M, int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
void chol_d_tb(
  const ac_fixed<W, I, S, Q, O> A[M][M],
  double L_tb[M][M]
)
{
  double sum_Ajj_Ljk_sq, sum_Aij_Lik_Ljk;

  // All elements of output initialized to zero by default
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {L_tb[i][j] = 0;}
  }

  // L_tb = 0;

  for (int j = 0; j < (int)M; j++) {
    sum_Ajj_Ljk_sq = A[j][j].to_double();

    for (int k = 0; k < j; k++) {
      sum_Ajj_Ljk_sq -= L_tb[j][k] * L_tb[j][k];
    }

    // Check to make sure that the matrix is positive definite. If "sum_Ajj_Ljk_sq" is negative/zero, then the diagonal
    // element, i.e. L_tb(j, j) will be complex/zero, which is not valid. This condition will not be encountered if the
    // input matrix is positive definite
    assert(sum_Ajj_Ljk_sq > 0);
    // Assign value to diagonal elements.
    L_tb[j][j] = sqrt(sum_Ajj_Ljk_sq);

    for (int i = (j+1); i < (int)M; i++) {
      sum_Aij_Lik_Ljk = A[i][j].to_double();
      for (int k = 0; k < j; k++) {
        sum_Aij_Lik_Ljk -= L_tb[i][k] * L_tb[j][k];
      }
      // Assign value to non-diagonal elements below the diagonal.
      L_tb[i][j] = sum_Aij_Lik_Ljk / L_tb[j][j];
    }
  }
}

// Testbench for cholesky decomposition for ac_complex<ac_fixed> matrices
// The testbench uses the cholesky-crout algorithm
template<unsigned M, int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
void chol_d_tb(
  ac_complex<ac_fixed<W, I, S, Q, O> > A[M][M],
  ac_complex<double> L_tb[M][M]
)
{
  typedef ac_complex<double> output_type;
  output_type zero_complex(0, 0), sum_Ajj_Ljk_sq, sum_Aij_Lik_Ljk;

  // All elements of output initialized to zero by default
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {L_tb[i][j] = zero_complex;}
  }

  // L_tb = zero_complex;

  for (int j = 0; j < (int)M; j++) {
    sum_Ajj_Ljk_sq.r() = A[j][j].r().to_double();
    sum_Ajj_Ljk_sq.i() = A[j][j].i().to_double();

    for (int k = 0; k < j; k++) {
      sum_Ajj_Ljk_sq -= L_tb[j][k] * L_tb[j][k].conj();
    }

    // Check to make sure that the matrix is positive definite. If "sum_Ajj_Ljk_sq" is negative/zero, then the diagonal
    // element, i.e. L_tb(j, j) will be complex/zero, which is not valid. This condition will not be encountered if the
    // input matrix is positive definite
    assert(sum_Ajj_Ljk_sq.r() > 0);
    // Assign value to diagonal elements. Since the diagonal elements are real, only initialize the real part.
    L_tb[j][j].r() = sqrt(sum_Ajj_Ljk_sq.r());

    for (int i = (j+1); i < (int)M; i++) {
      sum_Aij_Lik_Ljk.r() = A[i][j].r().to_double();
      sum_Aij_Lik_Ljk.i() = A[i][j].i().to_double();
      for (int k = 0; k < j; k++) {
        sum_Aij_Lik_Ljk -= L_tb[i][k] * L_tb[j][k].conj();
      }
      // Assign value to non-diagonal elements below the diagonal.
      L_tb[i][j].r() = (1 / L_tb[j][j].r())*sum_Aij_Lik_Ljk.r();
      L_tb[i][j].i() = (1 / L_tb[j][j].r())*sum_Aij_Lik_Ljk.i();
    }
  }
}

// Testbench for cholesky decomposition for ac_float matrix
// The testbench uses the cholesky-crout algorithm
template<unsigned M, int W, int I, int E, ac_q_mode Q>
void chol_d_tb(
  const ac_float<W, I, E, Q> A[M][M],
  double L_tb[M][M]
)
{
  double sum_Ajj_Ljk_sq, sum_Aij_Lik_Ljk;

  // All elements of output initialized to zero by default
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {L_tb[i][j] = 0.0;}
  }

  for (int j = 0; j < (int)M; j++) {
    sum_Ajj_Ljk_sq = A[j][j].to_double();

    for (int k = 0; k < j; k++) {
      sum_Ajj_Ljk_sq -= L_tb[j][k] * L_tb[j][k];
    }

    // Check to make sure that the matrix is positive definite. If "sum_Ajj_Ljk_sq" is negative/zero, then the diagonal
    // element, i.e. L_tb(j, j) will be complex/zero, which is not valid. This condition will not be encountered if the
    // input matrix is positive definite
    assert(sum_Ajj_Ljk_sq > 0);
    // Assign value to diagonal elements.
    L_tb[j][j] = sqrt(sum_Ajj_Ljk_sq);

    for (int i = (j+1); i < (int)M; i++) {
      sum_Aij_Lik_Ljk = A[i][j].to_double();
      for (int k = 0; k < j; k++) {
        sum_Aij_Lik_Ljk -= L_tb[i][k] * L_tb[j][k];
      }
      // Assign value to non-diagonal elements below the diagonal.
      L_tb[i][j] = sum_Aij_Lik_Ljk / L_tb[j][j];
    }
  }
}

// Testbench for cholesky decomposition for ac_std_float matrix
// The testbench uses the cholesky-crout algorithm
template<unsigned M, int W, int E>
void chol_d_tb(
  const ac_std_float<W, E> A[M][M],
  double L_tb[M][M]
)
{
  double sum_Ajj_Ljk_sq, sum_Aij_Lik_Ljk;

  // All elements of output initialized to zero by default
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {L_tb[i][j] = 0.0;}
  }

  for (int j = 0; j < (int)M; j++) {
    sum_Ajj_Ljk_sq = A[j][j].to_ac_float().to_double();

    for (int k = 0; k < j; k++) {
      sum_Ajj_Ljk_sq -= L_tb[j][k] * L_tb[j][k];
    }

    // Check to make sure that the matrix is positive definite. If "sum_Ajj_Ljk_sq" is negative/zero, then the diagonal
    // element, i.e. L_tb(j, j) will be complex/zero, which is not valid. This condition will not be encountered if the
    // input matrix is positive definite
    assert(sum_Ajj_Ljk_sq > 0);
    // Assign value to diagonal elements.
    L_tb[j][j] = sqrt(sum_Ajj_Ljk_sq);

    for (int i = (j+1); i < (int)M; i++) {
      sum_Aij_Lik_Ljk = A[i][j].to_ac_float().to_double();
      for (int k = 0; k < j; k++) {
        sum_Aij_Lik_Ljk -= L_tb[i][k] * L_tb[j][k];
      }
      // Assign value to non-diagonal elements below the diagonal.
      L_tb[i][j] = sum_Aij_Lik_Ljk / L_tb[j][j];
    }
  }
}

// Testbench for cholesky decomposition for ac_ieee_float matrix
// The testbench uses the cholesky-crout algorithm
template<unsigned M, ac_ieee_float_format Format>
void chol_d_tb(
  const ac_ieee_float<Format> A[M][M],
  double L_tb[M][M]
) {
  double sum_Ajj_Ljk_sq, sum_Aij_Lik_Ljk;

  // All elements of output initialized to zero by default
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {L_tb[i][j] = 0.0;}
  }

  for (int j = 0; j < (int)M; j++) {
    sum_Ajj_Ljk_sq = A[j][j].to_ac_float().to_double();

    for (int k = 0; k < j; k++) {
      sum_Ajj_Ljk_sq -= L_tb[j][k] * L_tb[j][k];
    }

    // Check to make sure that the matrix is positive definite. If "sum_Ajj_Ljk_sq" is negative/zero, then the diagonal
    // element, i.e. L_tb(j, j) will be complex/zero, which is not valid. This condition will not be encountered if the
    // input matrix is positive definite
    assert(sum_Ajj_Ljk_sq > 0);
    // Assign value to diagonal elements.
    L_tb[j][j] = sqrt(sum_Ajj_Ljk_sq);

    for (int i = (j+1); i < (int)M; i++) {
      sum_Aij_Lik_Ljk = A[i][j].to_ac_float().to_double();
      for (int k = 0; k < j; k++) {
        sum_Aij_Lik_Ljk -= L_tb[i][k] * L_tb[j][k];
      }
      // Assign value to non-diagonal elements below the diagonal.
      L_tb[i][j] = sum_Aij_Lik_Ljk / L_tb[j][j];
    }
  }
}

// Return the absolute value of the matrix element that has the
// maximum absolute value.
template<unsigned M>
double abs_mat_max(
  const double L_tb[M][M]
)
{
  double max_val = 0;

  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      if (abs(L_tb[i][j]) > max_val) {max_val = abs(L_tb[i][j]);}
    }
  }

  return max_val;
}

// Return the absolute value of the matrix element that has the
// maximum absolute value.
template<unsigned M>
double abs_mat_max(
  const ac_complex<double> L_tb[M][M]
)
{
  double max_val = 0;

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      if (L_tb[i][j].mag_sqr() > max_val) {max_val = L_tb[i][j].mag_sqr();}
    }
  }

  return max_val;
}

// Keep real, double element as it is (just kept in order to ensure that the error checking
// function is compatible even for real, double values)
double conv_val(double x)
{
  return x;
}

// Convert real, ac_fixed element to double.
template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
double conv_val(ac_fixed<W, I, S, Q, O> x)
{
  return x.to_double();
}

// Convert complex double element to mag_sqr.
double conv_val(ac_complex<double> x)
{
  return x.mag_sqr();
}


// Convert complex ac_fixed element to the double of it's
// mag_sqr
template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
double conv_val(ac_complex<ac_fixed<W, I, S, Q, O> > x)
{
  return x.mag_sqr().to_double();
}

// Convert real, ac_float element to double.
template <int W, int I, int E, ac_q_mode Q>
double conv_val(ac_float<W, I, E, Q> x) {
  return x.to_double();
}

// Convert real, ac_std_float element to double.
template <int W, int E>
double conv_val(ac_std_float<W, E> x) {
  return x.to_ac_float().to_double();
}

// Convert real, ac_ieee_float element to double.
template <ac_ieee_float_format Format>
double conv_val(ac_ieee_float<Format> x) {
  return x.to_ac_float().to_double();
}

// Compare DUT output matrix to testbench computed matrix and find error.
template<unsigned M, class T, class T_tb>
double compare_matrices(
  const T L[M][M],
  const T_tb L_tb[M][M],
  const double allowed_error
)
{
  double this_error, max_error = 0, max_val;
  // Find the max. abs. value in the matrix. In case of complex matrices, this is the max.
  // mag_sqr value. For real matrices, it's the max. absolute value stored in the matrix.
  max_val = abs_mat_max<M>(L_tb);

  #ifdef DEBUG
  cout << "max_val = " << max_val << endl;
  #endif

  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      // The error is calculated as the difference between the value of the expected
      // vs. the actual value, normalized w.r.t. the max_value in the matrix. For complex
      // numbers, the expected vs actual values are first converted to the their
      // mag_sqr() representations. For real numbers, the values are passed as they are
      // for the error calculation.
      this_error = 100.0 * abs( conv_val(L[i][j]) - conv_val(L_tb[i][j]) ) / (max_val);
      if (this_error > max_error) { max_error = this_error;}
      #ifdef DEBUG
      cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl;
      cout << "L[" << i << "][" << j << "]    = " << L[i][j] << endl;
      cout << "L_tb[" << i << "][" << j << "] = " << L_tb[i][j] << endl;
      cout << "this_error = " << this_error << endl;
      assert(this_error < allowed_error);
      #endif
    }
  }

  return max_error;

}

// Check if real ac_fixed/ac_float matrix is zero matrix.
template<unsigned M, class T>
bool check_if_zero_matrix(
  const T L[M][M]
)
{
  bool is_zero_matrix = true;

  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      if (L[i][j] != 0) {
        is_zero_matrix = false;
        #ifdef DEBUG
        cout << "Matrix was not zero for a non positive definite input. Non-zero element found:" << endl;
        cout << "L[" << i << "][" << j << "] = " << L[i][j] << endl;
        assert(false);
        #endif
      }
    }
  }

  return is_zero_matrix;
}

// Check if complex matrix is zero matrix.
template<unsigned M, class T>
bool check_if_zero_matrix(
  const ac_complex<T> L[M][M]
)
{
  bool is_zero_matrix = true;

  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      if (L[i][j].r() != 0 || L[i][j].i() != 0) {
        is_zero_matrix = false;
        #ifdef DEBUG
        cout << "Matrix was not zero for a non positive definite input. Non-zero element found:" << endl;
        cout << "L[" << i << "][" << j << "] = " << L[i][j] << endl;
        assert(false);
        #endif
      }
    }
  }

  return is_zero_matrix;
}

// Check if ac_std_float matrix is zero matrix.
template <unsigned M, int W, int E>
bool check_if_zero_matrix(
  const ac_std_float<W, E> L[M][M]
) {
  typedef ac_std_float<W, E> T;
  bool is_zero_matrix = true;

  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      if (L[i][j] != T::zero()) {
        is_zero_matrix = false;
        #ifdef DEBUG
        cout << "Matrix was not zero for a non positive definite input. Non-zero element found:" << endl;
        cout << "L[" << i << "][" << j << "] = " << L[i][j] << endl;
        assert(false);
        #endif
      }
    }
  }

  return is_zero_matrix;
}

// Check if ac_ieee_float matrix is zero matrix.
template <unsigned M, ac_ieee_float_format Format>
bool check_if_zero_matrix(
  const ac_ieee_float<Format> L[M][M]
) {
  typedef ac_ieee_float<Format> T;
  bool is_zero_matrix = true;

  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      if (L[i][j] != T::zero()) {
        is_zero_matrix = false;
        #ifdef DEBUG
        cout << "Matrix was not zero for a non positive definite input. Non-zero element found:" << endl;
        cout << "L[" << i << "][" << j << "] = " << L[i][j] << endl;
        assert(false);
        #endif
      }
    }
  }

  return is_zero_matrix;
}

// Copy a C-style array's contents over to an ac_matrix.
template<unsigned M, class T_matrix, class T_ac_matrix>
void copy_to_ac_matrix(
  const T_matrix array_2D[M][M],
  T_ac_matrix &output
)
{
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {output(i, j) = array_2D[i][j];}
  }
}

// Copy an ac_matrix's contents over to a C-style array.
template<unsigned M, class T_matrix, class T_ac_matrix>
void copy_to_array_2D(
  const T_ac_matrix &input,
  T_matrix array_2D[M][M]
)
{
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {array_2D[i][j] = input(i, j);}
  }
}

// ==============================================================================
// Functions: test_driver functions
// Description: Templatized functions that can be configured for certain bit-
//   widths of AC datatypes. They use the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   DUT cholesky decomposition with the computed cholesky decomposition using a
//   standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

// ==============================================================================
// Function: test_driver_fixed()
// Description: test_driver function for ac_fixed and ac_complex<ac_fixed> inputs
//   and outputs.

// FBfi = number of fractional bits in input type.
template <bool use_pwl, unsigned M, int FBfi, int outWfi, int outIfi, bool outSfi>
int test_driver_fixed(
  double &cumulative_max_error,
  double &cumulative_max_error_cmplx,
  const double allowed_error
)
{
  enum {
    Ifi = ac::nbits<M - 1>::val,
    Wfi = FBfi + Ifi,
    Ifi_c = Ifi + 1,
    Wfi_c = Wfi + 1
  };

  bool passed = true;

  ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP> A_C_array[M][M];
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> L_C_array[M][M];
  ac_complex<ac_fixed<Wfi_c, Ifi_c, true, AC_TRN, AC_WRAP> > cmplx_A_C_array[M][M];
  ac_complex<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> > cmplx_L_C_array[M][M];
  ac_matrix<ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP>, M, M> A_ac_matrix;
  ac_matrix<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP>, M, M> L_ac_matrix;
  ac_matrix<ac_complex<ac_fixed<Wfi_c, Ifi_c, true, AC_TRN, AC_WRAP> >, M, M> cmplx_A_ac_matrix;
  ac_matrix<ac_complex<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> >, M, M> cmplx_L_ac_matrix;

  if (use_pwl) {
    cout << "TEST: ac_chol_d() with PWL fns. M = "; // ac_chol_d uses PWL functions.
  } else {
    cout << "TEST: ac_chol_d() with acc fns. M = "; // ac_chol_d uses accurate functions.
  }
  cout.width(2);
  cout << left << M;
  cout << ", INPUT: ";
  cout.width(37);
  cout << left << A_C_array[0][0].type_name();
  cout << "OUTPUT: ";
  cout.width(37);
  cout << left << L_C_array[0][0].type_name();
  cout << "RESULT: ";

  double L_tb[M][M];
  ac_complex<double> cmplx_L_tb[M][M];

  // The gen_matrix function takes an MxM matrix, and multiplies it by its
  // conjugate transpose to obtain a positive definite input matrix
  gen_matrix<M>(A_C_array);
  gen_matrix<M>(cmplx_A_C_array);

  copy_to_ac_matrix<M>(A_C_array, A_ac_matrix);
  copy_to_ac_matrix<M>(cmplx_A_C_array, cmplx_A_ac_matrix);

  test_ac_chol_d_fixed<M, use_pwl>(A_C_array, L_C_array, cmplx_A_C_array, cmplx_L_C_array, A_ac_matrix, L_ac_matrix, cmplx_A_ac_matrix, cmplx_L_ac_matrix);

  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> L_ac_matrix_converted[M][M];
  ac_complex<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> > L_cmplx_ac_matrix_converted[M][M];

  copy_to_array_2D<M>(L_ac_matrix, L_ac_matrix_converted);
  copy_to_array_2D<M>(cmplx_L_ac_matrix, L_cmplx_ac_matrix_converted);

  // Get output of testbench function for cholesky decomposition.
  chol_d_tb<M>(A_C_array, L_tb);
  chol_d_tb<M>(cmplx_A_C_array, cmplx_L_tb);

  #ifdef DEBUG
  cout << "A_C_array = " << endl;
  print_matrix<M>(A_C_array);
  cout << "L_C_array = " << endl;
  print_matrix<M>(L_C_array);
  cout << "L_tb = " << endl;
  print_matrix<M>(L_tb);
  cout << "cmplx_A_C_array = " << endl;
  print_matrix<M>(cmplx_A_C_array);
  cout << "cmplx_L_C_array = " << endl;
  print_matrix<M>(cmplx_L_C_array);
  cout << "cmplx_L_tb = " << endl;
  print_matrix<M>(cmplx_L_tb);
  cout << "A_ac_matrix = " << endl;
  cout << A_ac_matrix << endl;
  cout << "L_ac_matrix = " << endl;
  cout << L_ac_matrix << endl;
  cout << "cmplx_A_ac_matrix = " << endl;
  cout << cmplx_A_ac_matrix << endl;
  cout << "cmplx_L_ac_matrix = " << endl;
  cout << cmplx_L_ac_matrix << endl;
  #endif

  // Compare matrices and get the max error
  double max_error = compare_matrices<M>(L_C_array, L_tb, allowed_error);
  double max_error_cmplx = compare_matrices<M>(cmplx_L_C_array, cmplx_L_tb, allowed_error);
  double max_error_ac_matrix = compare_matrices<M>(L_ac_matrix_converted, L_tb, allowed_error);
  double max_error_cmplx_ac_matrix = compare_matrices<M>(L_cmplx_ac_matrix_converted, cmplx_L_tb, allowed_error);

  // Put max overall error in a separate variable.
  double max_error_overall = max_error > max_error_ac_matrix ? max_error : max_error_ac_matrix;
  double max_error_cmplx_overall = max_error_cmplx > max_error_cmplx_ac_matrix ? max_error_cmplx : max_error_cmplx_ac_matrix;

  passed = (max_error_overall < allowed_error) && (max_error_cmplx_overall < allowed_error);

  // Also, we must make sure that the output on passing a non-positive definite matrix is a zero matrix. To do this, we pass a matrix
  // with all the values set to the quantum values of the ac_fixed type as the input.
  ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP> ac_fixed_quantum_value;
  ac_fixed_quantum_value.template set_val<AC_VAL_QUANTUM>();
  ac_complex<ac_fixed<Wfi_c, Ifi_c, false, AC_TRN, AC_WRAP> > ac_complex_quantum_value(ac_fixed_quantum_value, ac_fixed_quantum_value);
  A_ac_matrix = ac_fixed_quantum_value;
  cmplx_A_ac_matrix = ac_complex_quantum_value;

  // Copy over a non-positive definite matrix to the standard C array inputs.
  copy_to_array_2D<M>(A_ac_matrix, A_C_array);
  copy_to_array_2D<M>(cmplx_A_ac_matrix, cmplx_A_C_array);

  test_ac_chol_d_fixed<M, use_pwl>(A_C_array, L_C_array, cmplx_A_C_array, cmplx_L_C_array, A_ac_matrix, L_ac_matrix, cmplx_A_ac_matrix, cmplx_L_ac_matrix);

  copy_to_array_2D<M>(L_ac_matrix, L_ac_matrix_converted);
  copy_to_array_2D<M>(cmplx_L_ac_matrix, L_cmplx_ac_matrix_converted);

  // Make sure that a zero matrix is returned at the output.
  passed = passed && check_if_zero_matrix<M>(L_C_array) && check_if_zero_matrix<M>(cmplx_L_C_array) && check_if_zero_matrix<M>(L_ac_matrix_converted) && check_if_zero_matrix<M>(L_cmplx_ac_matrix_converted);

  if (passed) { printf("PASSED , max err (%f) (%f complex)\n", max_error_overall, max_error_cmplx_overall); }
  else        { printf("FAILED , max err (%f) (%f complex)\n", max_error_overall, max_error_cmplx_overall); } // LCOV_EXCL_LINE

  if (max_error_overall > cumulative_max_error) { cumulative_max_error = max_error_overall; }
  if (max_error_cmplx_overall > cumulative_max_error_cmplx) { cumulative_max_error_cmplx = max_error_cmplx_overall; }

  return 0;
}

// ==============================================================================
// Function: test_driver_float()
// Description: test_driver function for ac_float inputs and outputs.

// FBfl = number of fractional bits in input mantissa type.
template <bool use_pwl, unsigned M, int FBfl, int outWfl, int outIfl, int outEfl>
int test_driver_float(
  double &cumulative_max_error,
  const double allowed_error
)
{
  bool passed = true;

  enum {
    Ifl = 1,
    Wfl = Ifl + FBfl,
    E_val_1 = ac::nbits<Wfl - 1>::val,
    max_int_val = ac::nbits<M - 1>::val,
    val_2 = (max_int_val - (Ifl - 1) < 0) ? -(max_int_val - (Ifl - 1)) : max_int_val - (Ifl - 1),
    E_val_2 = ac::nbits<AC_MAX(int(val_2), 1)>::val,
    Efl = AC_MAX(int(E_val_1), int(E_val_2)) + 1,
  };

  ac_float<Wfl, Ifl, Efl, AC_TRN> A_C_array[M][M];
  ac_float<outWfl, outIfl, outEfl, AC_TRN> L_C_array[M][M];
  ac_matrix<ac_float<Wfl, Ifl, Efl, AC_TRN>, M, M> A_ac_matrix;
  ac_matrix<ac_float<outWfl, outIfl, outEfl, AC_TRN>, M, M> L_ac_matrix;

  if (use_pwl) {
    cout << "TEST: ac_chol_d() with PWL fns. M = "; // ac_chol_d uses PWL functions.
  } else {
    cout << "TEST: ac_chol_d() with acc fns. M = "; // ac_chol_d uses accurate functions.
  }
  cout.width(2);
  cout << left << M;
  cout << ", INPUT: ";
  cout.width(37);
  cout << left << A_C_array[0][0].type_name();
  cout << "OUTPUT: ";
  cout.width(37);
  cout << left << L_C_array[0][0].type_name();
  cout << "RESULT: ";

  double L_tb[M][M];

  // The gen_matrix function takes an MxM matrix, and multiplies it by its
  // transpose to obtain a positive definite input matrix
  gen_matrix<M>(A_C_array);
  copy_to_ac_matrix<M>(A_C_array, A_ac_matrix);

  test_ac_chol_d_float<M, use_pwl>(A_C_array, L_C_array, A_ac_matrix, L_ac_matrix);

  ac_float<outWfl, outIfl, outEfl, AC_TRN> L_ac_matrix_converted[M][M];

  copy_to_array_2D<M>(L_ac_matrix, L_ac_matrix_converted);

  // Get output of testbench function for cholesky decomposition.
  chol_d_tb<M>(A_C_array, L_tb);

  #ifdef DEBUG
  cout << "A_C_array = " << endl;
  print_matrix<M>(A_C_array);
  cout << "L_C_array = " << endl;
  print_matrix<M>(L_C_array);
  cout << "L_tb = " << endl;
  print_matrix<M>(L_tb);
  cout << "A_ac_matrix = " << endl;
  cout << A_ac_matrix << endl;
  cout << "L_ac_matrix = " << endl;
  cout << L_ac_matrix << endl;
  #endif

  // Compare matrices and get the max error
  double max_error = compare_matrices<M>(L_C_array, L_tb, allowed_error);
  double max_error_ac_matrix = compare_matrices<M>(L_ac_matrix_converted, L_tb, allowed_error);

  // Put max overall error in a separate variable.
  double max_error_overall = max_error > max_error_ac_matrix ? max_error : max_error_ac_matrix;

  passed = (max_error_overall < allowed_error);

  // Also, we must make sure that the output on passing a non-positive definite matrix is a zero matrix. To do this, we pass a matrix
  // with all the values set to the quantum values of the ac_float type as the input.
  ac_float<Wfl, Ifl, Efl, AC_TRN> ac_float_quantum_value;
  ac_float_quantum_value.template set_val<AC_VAL_QUANTUM>();
  A_ac_matrix = ac_float_quantum_value;

  // Copy over a non-positive definite matrix to the standard C array inputs.
  copy_to_array_2D<M>(A_ac_matrix, A_C_array);

  test_ac_chol_d_float<M, use_pwl>(A_C_array, L_C_array, A_ac_matrix, L_ac_matrix);

  copy_to_array_2D<M>(L_ac_matrix, L_ac_matrix_converted);

  // Make sure that a zero matrix is returned at the output.
  passed = passed && check_if_zero_matrix<M>(L_C_array) && check_if_zero_matrix<M>(L_ac_matrix_converted);

  if (passed) { printf("PASSED , max err (%f)\n", max_error_overall); }
  else        { printf("FAILED , max err (%f)\n", max_error_overall); } // LCOV_EXCL_LINE

  if (max_error_overall > cumulative_max_error) { cumulative_max_error = max_error_overall; }

  return 0;
}

// ==============================================================================
// Function: test_driver_stfloat()
// Description: test_driver function for ac_std_float inputs and outputs.

// FBstfl: Number of fractional bits in input mantissa, i.e. number of bits in
// significand field of ac_std_float datatype.
template <bool use_pwl, unsigned M, int FBstfl, int outWstfl, int outEstfl>
int test_driver_stfloat(
  double &cumulative_max_error,
  const double allowed_error
)
{
  bool passed = true;

  enum {
    E_val_1 = ac::nbits<FBstfl>::val,
    E_val_2 = ac::nbits<AC_MAX(ac::nbits<M - 1>::val - 1, 1)>::val,
    Estfl = AC_MAX(int(E_val_1), int(E_val_2)) + 1,
    Wstfl = 1 + Estfl + FBstfl,
  };

  typedef ac_std_float<Wstfl, Estfl> T_in;
  typedef ac_std_float<outWstfl, outEstfl> T_out;

  T_in A_C_array[M][M];
  T_out L_C_array[M][M];
  ac_matrix<T_in, M, M> A_ac_matrix;
  ac_matrix<T_out, M, M> L_ac_matrix;

  if (use_pwl) {
    cout << "TEST: ac_chol_d() with PWL fns. M = "; // ac_chol_d uses PWL functions.
  } else {
    cout << "TEST: ac_chol_d() with acc fns. M = "; // ac_chol_d uses accurate functions.
  }
  cout.width(2);
  cout << left << M;
  cout << ", INPUT: ";
  cout.width(37);
  cout << left << type_string_st<T_in>::type_string();
  cout << "OUTPUT: ";
  cout.width(37);
  cout << left << type_string_st<T_out>::type_string();
  cout << "RESULT: ";

  double L_tb[M][M];

  // The gen_matrix function takes an MxN matrix, and multiplies it by its
  // transpose to obtain a positive definite input matrix
  gen_matrix<M>(A_C_array);
  copy_to_ac_matrix<M>(A_C_array, A_ac_matrix);

  test_ac_chol_d_stfloat<M, use_pwl>(A_C_array, L_C_array, A_ac_matrix, L_ac_matrix);

  T_out L_ac_matrix_converted[M][M];

  copy_to_array_2D<M>(L_ac_matrix, L_ac_matrix_converted);

  // Get output of testbench function for cholesky decomposition.
  chol_d_tb<M>(A_C_array, L_tb);

  #ifdef DEBUG
  cout << "A_C_array = " << endl;
  print_matrix<M>(A_C_array);
  cout << "L_C_array = " << endl;
  print_matrix<M>(L_C_array);
  cout << "L_tb = " << endl;
  print_matrix<M>(L_tb);
  cout << "A_ac_matrix = " << endl;
  cout << A_ac_matrix << endl;
  cout << "L_ac_matrix = " << endl;
  cout << L_ac_matrix << endl;
  #endif

  // Compare matrices and get the max error
  double max_error = compare_matrices<M>(L_C_array, L_tb, allowed_error);
  double max_error_ac_matrix = compare_matrices<M>(L_ac_matrix_converted, L_tb, allowed_error);

  // Put max overall error in a separate variable.
  double max_error_overall = max_error > max_error_ac_matrix ? max_error : max_error_ac_matrix;

  passed = (max_error_overall < allowed_error);

  // Also, we must make sure that the output on passing a non-positive definite matrix is a zero matrix. To do this, we pass a matrix
  // with all the values set to unity.
  A_ac_matrix = T_in::one();

  // Copy over a non-positive definite matrix to the standard C array inputs.
  copy_to_array_2D<M>(A_ac_matrix, A_C_array);

  test_ac_chol_d_stfloat<M, use_pwl>(A_C_array, L_C_array, A_ac_matrix, L_ac_matrix);

  copy_to_array_2D<M>(L_ac_matrix, L_ac_matrix_converted);

  // Make sure that a zero matrix is returned at the output.
  passed = passed && check_if_zero_matrix<M>(L_C_array) && check_if_zero_matrix<M>(L_ac_matrix_converted);

  if (passed) { printf("PASSED , max err (%f)\n", max_error_overall); }
  else        { printf("FAILED , max err (%f)\n", max_error_overall); } // LCOV_EXCL_LINE

  if (max_error_overall > cumulative_max_error) { cumulative_max_error = max_error_overall; }

  return 0;
}

// ==============================================================================
// Function: test_driver_ifloat()
// Description: test_driver function for ac_ieee_float inputs and outputs.

template <bool use_pwl, unsigned M, ac_ieee_float_format in_format, ac_ieee_float_format out_format>
int test_driver_ifloat(
  double &cumulative_max_error,
  const double allowed_error
)
{
  bool passed = true;

  typedef ac_ieee_float<in_format> T_in;
  typedef ac_ieee_float<out_format> T_out;

  T_in A_C_array[M][M];
  T_out L_C_array[M][M];
  ac_matrix<T_in, M, M> A_ac_matrix;
  ac_matrix<T_out, M, M> L_ac_matrix;

  if (use_pwl) {
    cout << "TEST: ac_chol_d() with PWL fns. M = "; // ac_chol_d uses PWL functions.
  } else {
    cout << "TEST: ac_chol_d() with acc fns. M = "; // ac_chol_d uses accurate functions.
  }
  cout.width(2);
  cout << left << M;
  cout << ", INPUT: ";
  cout.width(37);
  cout << left << type_string_st<T_in>::type_string();
  cout << "OUTPUT: ";
  cout.width(37);
  cout << left << type_string_st<T_out>::type_string();
  cout << "RESULT: ";

  double L_tb[M][M];

  // The gen_matrix function takes an MxN matrix, and multiplies it by its
  // transpose to obtain a positive definite input matrix
  gen_matrix<M>(A_C_array);
  copy_to_ac_matrix<M>(A_C_array, A_ac_matrix);

  test_ac_chol_d_ifloat<M, use_pwl>(A_C_array, L_C_array, A_ac_matrix, L_ac_matrix);

  T_out L_ac_matrix_converted[M][M];

  copy_to_array_2D<M>(L_ac_matrix, L_ac_matrix_converted);

  // Get output of testbench function for cholesky decomposition.
  chol_d_tb<M>(A_C_array, L_tb);

  #ifdef DEBUG
  cout << "A_C_array = " << endl;
  print_matrix<M>(A_C_array);
  cout << "L_C_array = " << endl;
  print_matrix<M>(L_C_array);
  cout << "L_tb = " << endl;
  print_matrix<M>(L_tb);
  cout << "A_ac_matrix = " << endl;
  cout << A_ac_matrix << endl;
  cout << "L_ac_matrix = " << endl;
  cout << L_ac_matrix << endl;
  #endif

  // Compare matrices and get the max error
  double max_error = compare_matrices<M>(L_C_array, L_tb, allowed_error);
  double max_error_ac_matrix = compare_matrices<M>(L_ac_matrix_converted, L_tb, allowed_error);

  // Put max overall error in a separate variable.
  double max_error_overall = max_error > max_error_ac_matrix ? max_error : max_error_ac_matrix;

  passed = (max_error_overall < allowed_error);

  // Also, we must make sure that the output on passing a non-positive definite matrix is a zero matrix. To do this, we pass a matrix
  // with all the values set to unity.
  A_ac_matrix = T_in::one();

  // Copy over a non-positive definite matrix to the standard C array inputs.
  copy_to_array_2D<M>(A_ac_matrix, A_C_array);

  test_ac_chol_d_ifloat<M, use_pwl>(A_C_array, L_C_array, A_ac_matrix, L_ac_matrix);

  copy_to_array_2D<M>(L_ac_matrix, L_ac_matrix_converted);

  // Make sure that a zero matrix is returned at the output.
  passed = passed && check_if_zero_matrix<M>(L_C_array) && check_if_zero_matrix<M>(L_ac_matrix_converted);

  if (passed) { printf("PASSED , max err (%f)\n", max_error_overall); }
  else        { printf("FAILED , max err (%f)\n", max_error_overall); } // LCOV_EXCL_LINE

  if (max_error_overall > cumulative_max_error) { cumulative_max_error = max_error_overall; }

  return 0;
}

int main(int argc, char *argv[])
{
  double max_error_pwl = 0, cmplx_max_error_pwl = 0, max_error_acc = 0, cmplx_max_error_acc = 0;
  double allowed_error_pwl = 4;
  double allowed_error_acc = 0.005;

  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_chol_d(), for scalar and complex datatypes - allowed_error_pwl = " << allowed_error_pwl << ", allowed_error_acc = " << allowed_error_acc << endl;

  // template <bool use_pwl, unsigned M, int FBfi, int outWfi, int outIfi, bool outSfi>
  test_driver_fixed< true,  7, 16, 64, 32, true>(max_error_pwl, cmplx_max_error_pwl, allowed_error_pwl);
  test_driver_fixed< true,  8, 16, 64, 32, true>(max_error_pwl, cmplx_max_error_pwl, allowed_error_pwl);
  test_driver_fixed< true, 10, 16, 64, 32, true>(max_error_pwl, cmplx_max_error_pwl, allowed_error_pwl);
  test_driver_fixed< true, 12, 16, 64, 32, true>(max_error_pwl, cmplx_max_error_pwl, allowed_error_pwl);

  test_driver_fixed<false,  7, 16, 64, 32, true>(max_error_acc, cmplx_max_error_acc, allowed_error_acc);
  test_driver_fixed<false,  8, 16, 64, 32, true>(max_error_acc, cmplx_max_error_acc, allowed_error_acc);
  test_driver_fixed<false,  9, 16, 64, 32, true>(max_error_acc, cmplx_max_error_acc, allowed_error_acc);
  test_driver_fixed<false, 10, 16, 64, 32, true>(max_error_acc, cmplx_max_error_acc, allowed_error_acc);
  test_driver_fixed<false, 11, 16, 64, 32, true>(max_error_acc, cmplx_max_error_acc, allowed_error_acc);
  test_driver_fixed<false, 12, 16, 64, 32, true>(max_error_acc, cmplx_max_error_acc, allowed_error_acc);
  test_driver_fixed<false, 13, 16, 64, 32, true>(max_error_acc, cmplx_max_error_acc, allowed_error_acc);
  test_driver_fixed<false, 14, 16, 64, 32, true>(max_error_acc, cmplx_max_error_acc, allowed_error_acc);
  
  // template <bool use_pwl, unsigned M, int FBfl, int outWfl, int outIfl, int outEfl>
  test_driver_float< true,  7, 16, 32, 2, 10>(max_error_pwl, allowed_error_pwl);
  test_driver_float< true,  8, 16, 32, 2, 10>(max_error_pwl, allowed_error_pwl);
  test_driver_float< true, 10, 16, 32, 2, 10>(max_error_pwl, allowed_error_pwl);
  test_driver_float< true, 11, 16, 32, 2, 10>(max_error_pwl, allowed_error_pwl);
  test_driver_float< true, 12, 16, 32, 2, 10>(max_error_pwl, allowed_error_pwl);
  test_driver_float< true, 13, 16, 32, 2, 10>(max_error_pwl, allowed_error_pwl);
  test_driver_float< true, 14, 16, 32, 2, 10>(max_error_pwl, allowed_error_pwl);
  
  test_driver_float<false,  7, 16, 32, 2, 10>(max_error_acc, allowed_error_acc);
  test_driver_float<false,  8, 16, 32, 2, 10>(max_error_acc, allowed_error_acc);
  test_driver_float<false,  9, 16, 32, 2, 10>(max_error_acc, allowed_error_acc);
  test_driver_float<false, 10, 16, 32, 2, 10>(max_error_acc, allowed_error_acc);
  test_driver_float<false, 11, 16, 32, 2, 10>(max_error_acc, allowed_error_acc);
  test_driver_float<false, 12, 16, 32, 2, 10>(max_error_acc, allowed_error_acc);
  test_driver_float<false, 13, 16, 32, 2, 10>(max_error_acc, allowed_error_acc);
  test_driver_float<false, 14, 16, 32, 2, 10>(max_error_acc, allowed_error_acc);

  // template <bool use_pwl, unsigned M, int FBstfl, int outWstfl, int outEstfl>
  test_driver_stfloat< true,  7, 23, 64, 11>(max_error_pwl, allowed_error_pwl);
  test_driver_stfloat< true,  8, 23, 64, 11>(max_error_pwl, allowed_error_pwl);
  test_driver_stfloat< true, 10, 23, 64, 11>(max_error_pwl, allowed_error_pwl);
  test_driver_stfloat< true, 11, 23, 64, 11>(max_error_pwl, allowed_error_pwl);
  test_driver_stfloat< true, 12, 23, 64, 11>(max_error_pwl, allowed_error_pwl);
  test_driver_stfloat< true, 13, 23, 64, 11>(max_error_pwl, allowed_error_pwl);
  
  test_driver_stfloat<false,  7, 23, 64, 11>(max_error_acc, allowed_error_acc);
  test_driver_stfloat<false,  8, 23, 64, 11>(max_error_acc, allowed_error_acc);
  test_driver_stfloat<false,  9, 23, 64, 11>(max_error_acc, allowed_error_acc);
  test_driver_stfloat<false, 10, 23, 64, 11>(max_error_acc, allowed_error_acc);
  test_driver_stfloat<false, 11, 23, 64, 11>(max_error_acc, allowed_error_acc);
  test_driver_stfloat<false, 12, 23, 64, 11>(max_error_acc, allowed_error_acc);
  test_driver_stfloat<false, 13, 23, 64, 11>(max_error_acc, allowed_error_acc);
  test_driver_stfloat<false, 14, 23, 64, 11>(max_error_acc, allowed_error_acc);
  
  // template <bool use_pwl, unsigned M, ac_ieee_float_format in_format, ac_ieee_float_format out_format>
  test_driver_ifloat< true,  7, binary32, binary64>(max_error_pwl, allowed_error_pwl);
  test_driver_ifloat< true,  8, binary32, binary64>(max_error_pwl, allowed_error_pwl);
  test_driver_ifloat< true, 10, binary32, binary64>(max_error_pwl, allowed_error_pwl);
  test_driver_ifloat< true, 11, binary32, binary64>(max_error_pwl, allowed_error_pwl);
  test_driver_ifloat< true, 12, binary32, binary64>(max_error_pwl, allowed_error_pwl);
  test_driver_ifloat< true, 13, binary32, binary64>(max_error_pwl, allowed_error_pwl);
  
  test_driver_ifloat<false,  7, binary32, binary64>(max_error_acc, allowed_error_acc);
  test_driver_ifloat<false,  8, binary32, binary64>(max_error_acc, allowed_error_acc);
  test_driver_ifloat<false,  9, binary32, binary64>(max_error_acc, allowed_error_acc);
  test_driver_ifloat<false, 10, binary32, binary64>(max_error_acc, allowed_error_acc);
  test_driver_ifloat<false, 11, binary32, binary64>(max_error_acc, allowed_error_acc);
  test_driver_ifloat<false, 12, binary32, binary64>(max_error_acc, allowed_error_acc);
  test_driver_ifloat<false, 13, binary32, binary64>(max_error_acc, allowed_error_acc);
  test_driver_ifloat<false, 14, binary32, binary64>(max_error_acc, allowed_error_acc);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum errors observed across all data type / bit-width variations:" << endl;
  cout << "    max_error_pwl       = " << max_error_pwl << endl;
  cout << "    cmplx_max_error_pwl = " << cmplx_max_error_pwl << endl;
  cout << "    max_error_acc       = " << max_error_acc << endl;
  cout << "    cmplx_max_error_acc = " << cmplx_max_error_acc << endl;

  // If error limits on any tested datatype have been crossed, the test has failed
  bool test_fail = (max_error_pwl > allowed_error_pwl) || (cmplx_max_error_pwl > allowed_error_pwl) || (max_error_acc > allowed_error_acc) || (cmplx_max_error_acc > allowed_error_acc);

  // Notify the user whether or not the test was a failure.
  if (test_fail) {
    cout << "  ac_chol_d - FAILED - Error tolerance(s) exceeded" << endl; // LCOV_EXCL_LINE
    cout << "=============================================================================" << endl; // LCOV_EXCL_LINE
    return -1; // LCOV_EXCL_LINE
  } else {
    cout << "  ac_chol_d - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }
  return 0;

}
