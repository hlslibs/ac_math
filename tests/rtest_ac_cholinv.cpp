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
// ac_cholinv() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_cholinv.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_cholinv.h>
using namespace ac_math;

// ==============================================================================
// Test Designs
//   These simple functions allow executing the ac_cholinv() function
//   using multiple data types at the same time. Template parameters are
//   used to configure the bit-widths of the types.

// Test Design for real and complex fixed point values.
template <unsigned M, bool use_pwl, int Wfi, int Ifi, int outWfi, int outIfi, bool outSfi>
void test_ac_cholinv(
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
    ac_cholinv<M, use_pwl, use_pwl>(A1, L1);
    ac_cholinv<M, use_pwl, use_pwl>(A2, L2);
  #else
    ac_cholinv<use_pwl, use_pwl>(A1, L1);
    ac_cholinv<use_pwl, use_pwl>(A2, L2);
  #endif
  ac_cholinv<use_pwl, use_pwl>(A3, L3);
  ac_cholinv<use_pwl, use_pwl>(A4, L4);
}

// ==============================================================================

#include <math.h>
#include <random>
#include <limits>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

// ------------------------------------------------------------------------------
// Helper functions for positive definite matrix generation and error calculation

// Print matrix, for debugging purposes.

#ifdef DEBUG
template<unsigned M, class T>
void print_matrix(T mat[M][M])
{
  cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl;
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {cout << "mat[" << i << "][" << j << "] = " << mat[i][j] << endl;}
  }
}
#endif

// Generate positive definite matrix of ac_fixed values
template<unsigned M, int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
void gen_matrix(ac_matrix<ac_fixed<W, I, S, Q, O>, M, M> &A)
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
  print_matrix(tbmat);
  cout << "tbmatT is : " << endl;
  print_matrix(tbmatT);
  #endif

  // Multiply tbmat by its transpose to get the positive definite input matrix
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      A(i, j) = 0;
      for (int k = 0; k < (int)M; k++) {
        A(i, j) += tbmat[i][k] * tbmatT[k][j];
      }
    }
  }

  #ifdef DEBUG
  cout << "A in gen_matrix function is : " << endl;
  cout << A << endl;
  #endif
}

// Generate positive definite matrix of ac_complex<ac_fixed> values
template<unsigned M, int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
void gen_matrix(ac_matrix<ac_complex<ac_fixed<W, I, S, Q, O> >, M, M> &A)
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
  print_matrix(tbmat);
  cout << "tbmatT is : " << endl;
  print_matrix(tbmatT);
  #endif

  // Multiply tbmat by its transpose to get the positive definite input matrix
  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      A(i, j) = 0;
      for (int k = 0; k < (int)M; k++) {
        A(i, j) += tbmat[i][k] * tbmatT[k][j];
      }
    }
  }

  #ifdef DEBUG
  cout << "A in gen_matrix function is : " << endl;
  cout << A << endl;
  #endif
}

template<class T_tb>
void getCofactor(T_tb &A_tb, T_tb &temp, int p, int q, int n)
{
  int i = 0, j = 0;

  // Looping for each element of the matrix
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      // Copying into temporary matrix only those element
      // which are not in given row and column
      if (row != p && col != q) {
        temp(i,j) = A_tb(row,col);
        j++;
        // Row is filled, so increase row index and
        // reset col index
        if (j == n - 1) {
          j = 0;
          i++;
        }
      }
    }
  }
}

/* Recursive function for finding determinant of matrix.
n is current dimension of A[][]. */
template<unsigned M, class T_tb>
double determinant(ac_matrix<T_tb, M, M> &A_tb, unsigned n)
{
  double D = 0; // Initialize result

  // Base case : if matrix contains single element
  if (n == 1)
  { return A_tb(0,0); }

  ac_matrix<T_tb, M, M> temp; // To store cofactors

  int sign = 1; // To store sign multiplier

  // Iterate for each element of first row
  for (int f = 0; f < n; f++) {
    // Getting Cofactor of A[0][f]
    getCofactor(A_tb, temp, 0, f, n);
    D += sign * A_tb(0,f) * determinant(temp, n-1);

    // terms are to be added with alternate sign
    sign = -sign;
  }

  return D;
}

/* Recursive function for finding determinant of matrix.
n is current dimension of A[][]. */
template<unsigned M, class T_tb>
ac_complex<double> determinant(ac_matrix<ac_complex<T_tb>, M, M> &A_tb, unsigned n)
{
  ac_complex<double> D = 0; // Initialize result

  // Base case : if matrix contains single element
  if (n == 1)
  { return A_tb(0,0); }

  ac_matrix<ac_complex<T_tb>, M, M> temp; // To store cofactors

  int sign = 1; // To store sign multiplier

  // Iterate for each element of first row
  for (int f = 0; f < n; f++) {
    // Getting Cofactor of A[0][f]
    getCofactor(A_tb, temp, 0, f, n);
    D += sign * A_tb(0,f) * determinant(temp, n-1);

    // terms are to be added with alternate sign
    sign = -sign;
  }

  return D;
}

// Function to get adjoint of A[N][N] in adj[N][N].
template<unsigned M, class T_tb>
void adjoint(ac_matrix<T_tb, M, M> &A_tb, ac_matrix<T_tb, M, M> &adj)
{
  if (M == 1) {
    adj(0,0) = 1;
    return;
  }

  // temp is used to store cofactors of A[][]
  int sign = 1;
  ac_matrix<T_tb, M, M> temp;

  for (int i=0; i<M; i++) {
    for (int j=0; j<M; j++) {
      // Get cofactor of A[i][j]
      getCofactor(A_tb, temp, i, j, M);

      // sign of adj[j][i] positive if sum of row
      // and column indexes is even.
      sign = ((i+j)%2==0)? 1: -1;

      // Interchanging rows and columns to get the
      // transpose of the cofactor matrix
      adj(j,i) = (sign)*(determinant(temp, M-1));
    }
  }
}
// Function to get adjoint of A[N][N] in adj[N][N].
template<unsigned M, class T_tb>
void adjoint(ac_matrix<ac_complex<T_tb>, M, M> &A_tb, ac_matrix<ac_complex<T_tb>, M, M> &adj)
{
  if (M == 1) {
    adj(0,0).r() = 1;
    adj(0,0).i() = 1;
    return;
  }

  // temp is used to store cofactors of A[][]
  int sign = 1;
  ac_matrix<ac_complex<T_tb>, M, M> temp;

  for (int i=0; i<M; i++) {
    for (int j=0; j<M; j++) {
      // Get cofactor of A[i][j]
      getCofactor(A_tb, temp, i, j, M);

      // sign of adj[j][i] positive if sum of row
      // and column indexes is even.
      sign = ((i+j)%2==0)? 1: -1;

      // Interchanging rows and columns to get the
      // transpose of the cofactor matrix
      adj(j,i) = (sign)*(determinant(temp, M-1));
    }
  }
}

template<unsigned M, class T_tb>
void chol_inv_tb(
  ac_matrix<T_tb, M, M> &A_tb,
  ac_matrix<T_tb, M, M> &Ainv_tb
)
{
  // C++ program to find adjoint and inverse of a matrix

  // Function to get cofactor of A[p][q] in temp[][]. n is current
  // dimension of A[][]
  // Find determinant of A[][]
  //typedef ac_matrix<double, M, M> input_tb_datatype;

  double det = determinant(A_tb, M);

  if (det == 0) {
    cout << "Singular matrix, can't find its inverse" << endl;
  }

  // Find adjoint
  ac_matrix<T_tb, M, M> adj;
  adjoint(A_tb, adj);

  // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
  for (int i=0; i<M; i++)
    for (int j=0; j<M; j++)
    { Ainv_tb(i,j) = adj(i,j)/det; }
}

template<unsigned M, class T_tb>
void chol_inv_tb(
  ac_matrix<ac_complex<T_tb>, M, M> &A_tb,
  ac_matrix<ac_complex<T_tb>, M, M> &Ainv_tb
)
{
  // C++ program to find adjoint and inverse of a matrix

  // Function to get cofactor of A[p][q] in temp[][]. n is current
  // dimension of A[][]
  // Find determinant of A[][]
  //typedef ac_matrix<double, M, M> input_tb_datatype;
  ac_complex<double> det = determinant(A_tb, M);
  if (det == 0) {
    cout << "Singular matrix, can't find its inverse" << endl;
  }

  // Find adjoint
  ac_matrix<ac_complex<T_tb>, M, M> adj;
  adjoint(A_tb, adj);

  // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
  for (int i=0; i<M; i++)
    for (int j=0; j<M; j++)
    { Ainv_tb(i,j) = adj(i,j)/det; }

}

template<unsigned M>
double abs_mat_max(
  const ac_matrix<double, M, M> &Ainv_tb
)
{
  double max_val = 0;

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      if (abs(Ainv_tb(i, j)) > max_val) {max_val = abs(Ainv_tb(i, j));}
    }
  }

  return max_val;
}

template<unsigned M>
double abs_mat_max(
  const ac_matrix<ac_complex<double>, M, M> &Ainv_tb
)
{
  double max_val = 0;

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      if (Ainv_tb(i, j).mag_sqr() > max_val) {max_val = Ainv_tb(i, j).mag_sqr();}
    }
  }

  return max_val;
}

double conv_val(double x)
{
  return (x);
}

template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
double conv_val(ac_fixed<W, I, S, Q, O> x)
{
  return (x.to_double());
}

double conv_val(ac_complex<double> x)
{
  return x.mag_sqr();
}

template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
double conv_val(ac_complex<ac_fixed<W, I, S, Q, O> > x)
{
  return x.mag_sqr().to_double();
}

template<unsigned M, class T, class T_tb>
double compare_matrices(
  const ac_matrix<T, M, M> &Ainv,
  const ac_matrix<T_tb, M, M> &Ainv_tb,
  const double allowed_error
)
{
  double this_error, max_error = 0, max_val;
  //Find the max. abs. value in the matrix. In case of complex matrices, this is the max.
  //mag_sqr value. For real matrices, it's the max. absolute value stored in the matrix.
  max_val = abs_mat_max(Ainv_tb);

  for (int i = 0; i < (int)M; i++) {
    for (int j = 0; j < (int)M; j++) {
      //The error is calculated as the difference between the value of the expected
      //vs. the actual value, normalized w.r.t. the max_value in the matrix. For complex
      //numbers, the expected vs actual values are first converted to the their
      //mag_sqr() representations. For real numbers, the values are passed as they are
      //for the error calculation.
      this_error = 100 * abs( conv_val(Ainv(i, j)) - conv_val(Ainv_tb(i, j)) ) / (max_val);
      if (this_error > max_error) { max_error = this_error;}
      #ifdef DEBUG
      cout << "Ainv(i, j)    = " << Ainv(i, j) << endl;
      cout << "Ainv_tb(i, j) = " << Ainv_tb(i, j) << endl;
      cout << "this_error = " << this_error << endl;
      assert(this_error < allowed_error);
      #endif
    }
  }

  #ifdef DEBUG
  cout << "max_val = " << max_val << endl;
  #endif

  return max_error;

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
//   DUT cholesky inverse with the computed cholesky inverse using a
//   standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

// ==============================================================================
// Function: test_driver()
// Description: test_driver function for ac_fixed and ac_complex<ac_fixed> inputs
//   and outputs.

// FBfi = number of fractional bits in input type.
template <bool use_pwl, unsigned M, int FBfi, int outWfi, int outIfi, bool outSfi>
int test_driver(
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
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> Ainv_C_array[M][M];
  ac_complex<ac_fixed<Wfi_c, Ifi_c, true, AC_TRN, AC_WRAP> > cmplx_A_C_array[M][M];
  ac_complex<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> > cmplx_Ainv_C_array[M][M];
  ac_matrix<ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP>, M, M> A_ac_matrix;
  ac_matrix<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP>, M, M> Ainv_ac_matrix;
  ac_matrix<ac_complex<ac_fixed<Wfi_c, Ifi_c, true, AC_TRN, AC_WRAP> >, M, M> cmplx_A_ac_matrix;
  ac_matrix<ac_complex<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> >, M, M> cmplx_Ainv_ac_matrix;

  cout << "TEST: ac_cholinv(), with " << (use_pwl ? "PWL" : "accurate") << " functions. M = ";
  cout << M << ",";
  cout << "INPUT: ";
  cout.width(38);
  cout << left << A_C_array[0][0].type_name();
  cout << "OUTPUT: ";
  cout.width(38);
  cout << left << Ainv_C_array[0][0].type_name();
  cout << "RESULT: ";

  ac_matrix<double, M, M> Ainv_tb, A_tb;
  ac_matrix<ac_complex<double>, M, M> cmplx_Ainv_tb, cmplx_A_tb;

  gen_matrix(A_ac_matrix);
  gen_matrix(cmplx_A_ac_matrix);
  //CCS_DESIGN(project)(A, Ainv);

  copy_to_array_2D<M>(A_ac_matrix, A_C_array);
  copy_to_array_2D<M>(cmplx_A_ac_matrix, cmplx_A_C_array);

  test_ac_cholinv<M, use_pwl>(A_C_array, Ainv_C_array, cmplx_A_C_array, cmplx_Ainv_C_array, A_ac_matrix, Ainv_ac_matrix, cmplx_A_ac_matrix, cmplx_Ainv_ac_matrix);

  ac_matrix<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP>, M, M> Ainv_C_array_converted;
  ac_matrix<ac_complex<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> >, M, M> cmplx_Ainv_C_array_converted;

  copy_to_ac_matrix<M>(Ainv_C_array, Ainv_C_array_converted);
  copy_to_ac_matrix<M>(cmplx_Ainv_C_array, cmplx_Ainv_C_array_converted);

  //ac_matrix<double, M, M> A_tb;
  for (unsigned i=0; i<M; i++) {
    for (unsigned j=0; j<M; j++) {
      A_tb(i, j) = A_ac_matrix(i,j).to_double();
      cmplx_A_tb(i, j).r() = cmplx_A_ac_matrix(i,j).r().to_double();
      cmplx_A_tb(i, j).i() = cmplx_A_ac_matrix(i,j).i().to_double();
    }
  }

  chol_inv_tb(A_tb, Ainv_tb);
  chol_inv_tb(cmplx_A_tb, cmplx_Ainv_tb);

  #ifdef DEBUG
  cout << "A_ac_matrix = " << endl;
  cout << A_ac_matrix << endl;
  cout << "A_tb = " << endl;
  cout << A_tb << endl;
  cout << "Ainv_tb = " << endl;
  cout << Ainv_tb << endl;
  cout << "cmplx_A_ac_matrix = " << endl;
  cout << cmplx_A_ac_matrix << endl;
  cout << "cmplx_A_tb = " << endl;
  cout << cmplx_A_tb << endl;
  cout << "cmplx_Ainv_tb = " << endl;
  cout << cmplx_Ainv_tb << endl;
  cout << "A_C_array = " << endl;
  print_matrix(A_C_array);
  cout << "Ainv_C_array = " << endl;
  print_matrix(Ainv_C_array);
  cout << "Ainv_tb = " << endl;
  print_matrix(Ainv_tb);
  cout << "cmplx_A_C_array = " << endl;
  print_matrix(cmplx_A_C_array);
  cout << "cmplx_Ainv_C_array = " << endl;
  print_matrix(cmplx_Ainv_C_array);
  cout << "cmplx_Ainv_tb = " << endl;
  print_matrix(cmplx_Ainv_tb);
  cout << "A_ac_matrix = " << endl;
  cout << A_ac_matrix << endl;
  cout << "Ainv_ac_matrix = " << endl;
  cout << Ainv_ac_matrix << endl;
  cout << "cmplx_A_ac_matrix = " << endl;
  cout << cmplx_A_ac_matrix << endl;
  cout << "cmplx_Ainv_ac_matrix = " << endl;
  cout << cmplx_Ainv_ac_matrix << endl;
  #endif

  // Compare matrices and get the max error
  double max_error = compare_matrices(Ainv_C_array_converted, Ainv_tb, allowed_error);
  double max_error_cmplx = compare_matrices(cmplx_Ainv_C_array_converted, cmplx_Ainv_tb, allowed_error);
  double max_error_ac_matrix = compare_matrices(Ainv_ac_matrix, Ainv_tb, allowed_error);
  double max_error_cmplx_ac_matrix = compare_matrices(cmplx_Ainv_ac_matrix, cmplx_Ainv_tb, allowed_error);

  // Put max overall error in a separate variable.
  double max_error_overall = max_error > max_error_ac_matrix ? max_error : max_error_ac_matrix;
  double max_error_cmplx_overall = max_error_cmplx > max_error_cmplx_ac_matrix ? max_error_cmplx : max_error_cmplx_ac_matrix;

  passed = (max_error_overall < allowed_error) && (max_error_cmplx_overall < allowed_error);

  if (passed) { printf("PASSED , max err (%f) (%f complex)\n", max_error_overall, max_error_cmplx_overall); }
  else        { printf("FAILED , max err (%f) (%f complex)\n", max_error_overall, max_error_cmplx_overall); } // LCOV_EXCL_LINE

  if (max_error_overall > cumulative_max_error) { cumulative_max_error = max_error_overall; }
  if (max_error_cmplx_overall > cumulative_max_error_cmplx) { cumulative_max_error_cmplx = max_error_cmplx_overall; }

  return 0;
}

int main(int argc, char *argv[])
{
  double max_error_pwl = 0, cmplx_max_error_pwl = 0, max_error_accurate = 0, cmplx_max_error_accurate = 0;
  double allowed_error_pwl = 4;
  double allowed_error_accurate = 0.005;

  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_cholinv(), for scalar and complex datatypes - allowed_error_pwl = " << allowed_error_pwl << ", allowed_error_accurate = " << allowed_error_accurate << endl;

  // template <bool use_pwl, unsigned M, int FBfi, int outWfi, int outIfi, bool outSfi>
  test_driver<true, 5, 16, 64, 32, true> (max_error_pwl, cmplx_max_error_pwl, allowed_error_pwl);

  test_driver<false, 4, 16, 64, 32, true> (max_error_accurate, cmplx_max_error_accurate, allowed_error_accurate);
  test_driver<false, 5, 16, 64, 32, true> (max_error_accurate, cmplx_max_error_accurate, allowed_error_accurate);
  test_driver<false, 6, 16, 64, 32, true> (max_error_accurate, cmplx_max_error_accurate, allowed_error_accurate);
  test_driver<false, 7, 16, 64, 32, true> (max_error_accurate, cmplx_max_error_accurate, allowed_error_accurate);
  test_driver<false, 8, 16, 64, 32, true> (max_error_accurate, cmplx_max_error_accurate, allowed_error_accurate);
  test_driver<false, 9, 16, 64, 32, true> (max_error_accurate, cmplx_max_error_accurate, allowed_error_accurate);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum errors observed across all data type / bit-width variations:" << endl;
  cout << "    max_error_pwl            = " << max_error_pwl << endl;
  cout << "    cmplx_max_error_pwl      = " << cmplx_max_error_pwl << endl;
  cout << "    max_error_accurate       = " << max_error_accurate << endl;
  cout << "    cmplx_max_error_accurate = " << cmplx_max_error_accurate << endl;

  // If error limits on any tested datatype have been crossed, the test has failed
  bool test_fail = (max_error_pwl > allowed_error_pwl) || (cmplx_max_error_pwl > allowed_error_pwl) || (max_error_accurate > allowed_error_accurate) || (cmplx_max_error_accurate > allowed_error_accurate);

  // Notify the user whether or not the test was a failure.
  if (test_fail) {
    cout << "  ac_cholinv - FAILED - Error tolerance(s) exceeded" << endl; // LCOV_EXCL_LINE
    cout << "=============================================================================" << endl; // LCOV_EXCL_LINE
    return -1; // LCOV_EXCL_LINE
  } else {
    cout << "  ac_cholinv - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }
  return 0;

}
