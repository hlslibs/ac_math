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
//********************************************************************************************************************************
// File: ac_qrd.h
//
// Created on: Jul 11, 2017
//
// Author: Sachchidanand Deo
//
// Description: This function computes Q and R matrices which are decompositions of an
//    input matrix (given by A), such that:
//    A is an input square matrix with dimension MxM
//    Q matrix is an orthogonal matrix (such that, Q = inverse (Q)
//    and R is an upper triangular matrix.
//    Supported datatypes are: ac_fixed and ac_complex<ac_fixed>
//
// Revision History:
//    3.3.0  - [CAT-25797] Added CDesignChecker fixes/waivers for code check violations in ac_math PWL and Linear Algebra IPs.
//             Waivers added for CNS and ABW violations.
//             Fixes added for FXD and MXS violations.
//               - FXD violations fixed by changing integer literals to floating point literals or typecasting to ac_fixed values.
//               - MXS violations fixed by typecasting unsigned variables to int.
//    3.2.3 - CAT-24362, CAT-24269.
//
//********************************************************************************************************************************

#ifndef _INCLUDED_AC_QRD_H_
#define _INCLUDED_AC_QRD_H_

// Include headers for data types supported by these implementations
#include <ac_fixed.h>
#include <ac_complex.h>
#include <ac_matrix.h>

// Include headers for required functions
#include <ac_math/ac_div.h>
#include <ac_math/ac_sqrt.h>
#include <ac_math/ac_sqrt_pwl.h>
#include <ac_math/ac_inverse_sqrt_pwl.h>

#ifndef __SYNTHESIS__
#include <iostream>
#endif

namespace ac_math
{
  //=========================================================================
  // Inline functions to perform required matrix operations:
  //=========================================================================

  // Assign unity or zero value to output argument for ac_fixed representation.
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
  void assign_unity_or_zero(ac_fixed<W, I, S, Q, O> &output, bool assign_unity)
  {
    typedef ac_fixed<1, 1, false> oneBitType;

    output = assign_unity ? oneBitType(1) : oneBitType(0);
  }

  // Assign unity value to output argument for ac_complex<ac_fixed> representation
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
  void assign_unity_or_zero(ac_complex<ac_fixed<W, I, S, Q, O> > &output, bool assign_unity)
  {
    typedef ac_fixed<1, 1, false> oneBitType;

    #pragma hls_waive ABW
    output.r() = assign_unity ? oneBitType(1) : oneBitType(0);
    #pragma hls_waive ABW
    output.i() = 0.0;
  }

  // Initializing the larger matrix as well as the smaller matrices. Formation of A1 from A and identity matrix.
  template <unsigned M, typename T1, typename T2>
  void initialize_matrix (ac_matrix<T1, M, M> &A, ac_matrix <T2,M,2*M> &A1)
  {
    // For matrices that are mapped to register banks, IDENT_ROW and and IDENT_COLUMN can be fully unrolled.
    IDENT_ROW:
    for (unsigned i = 0; i < M; i++) {
      IDENT_COLUMN:
      for (unsigned j = 0; j < M; j++) {
        A1(i,j) = A(i,j); // formation of A1 using A
        // setting remaining part of A to identity matrix
        assign_unity_or_zero(A1(i,j+M), j+M-i == M);
        // Assigment to real and imaginary values for ac_complex variables is handled separately in the ac_complex version of the
        // assign_unity_or_zero function. This is done in order to avoid using the ac_complex constructor used to assign real
        // values, due to the fact that doing so results in FXD violations.
      }
    }
  }

  // This function is used to get back final Q and R from large matrix A1 on which Givens rotations are performed.
  // It also functions to modify the last column if real_diag is set to true for complex matrices.
  template <bool real_diag, unsigned M, typename T1, typename T2, typename T3>
  void qr_separate (ac_matrix <T1,M,2*M> &A1, ac_matrix <T2,M,M> &Q, ac_matrix <T2,M,M> &R, T3 e_coeff)
  {
    // For matrices that are mapped to register banks, SEPARATE_ROW and and SEPARATE_COLUMN can be fully unrolled.
    SEPARATE_ROW:
    for (unsigned i = 0; i < M; i++) {
      SEPARATE_COLUMN:
      for (unsigned j = 0; j < M; j++) {
        #pragma hls_waive CNS
        if (real_diag && i == M - 1) {
          Q(j, i) = A1(i, j + M)*e_coeff;
        } else {
          Q(j, i) = A1(i, j + M);
        }
        R(i,j) = A1(i,j); // R is a part where identity matrix was present
      }
    }
  }

  // This function calculates the sum of squares for real diagonal_PE calculations.
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O, int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO>
  void sqr_calc(ac_fixed<W, I, S, Q, O> in1, ac_fixed<W, I, S, Q, O> in2, ac_fixed<outW, outI, outS, outQ, outO> &out)
  {
    out = in1*in1 + in2*in2;
  }

  // This function calculates the sum of squares for complex diagonal_PE calculations.
  template<class T, int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO>
  void sqr_calc(ac_complex<T> in1, ac_complex<T> in2, ac_fixed<outW, outI, outS, outQ, outO> &out)
  {
    out = in1.mag_sqr() + in2.mag_sqr();
  }

  // This function is used to produce sine and cosine coefficients for the given's rotation transform, for real matrices.
  template <bool ispwl, class T1, class T2>
  void diagonal_PE (T1 b, T1 a, T2 &c, T2 &s)
  {
    typedef typename T1::rt_unary::mag_sqr s_type;
    enum {
      sqr_I = s_type::i_width,
      n_f_b = 32,
    };
    typedef ac_fixed<n_f_b + sqr_I, sqr_I, false> sqr_type;

    sqr_type sqr;
    sqr_calc(b, a, sqr);

    enum {
      root_I   = sqr_I%2 == 0 ? sqr_I/2 : (sqr_I + 1)/2,
      i_root_I = n_f_b%2 == 0 ? n_f_b/2 : (n_f_b + 1)/2,
    };

    ac_fixed<n_f_b + i_root_I, i_root_I, false> i_root;

    #pragma hls_waive CNS
    if (ispwl) {
      // Use PWL functions.
      ac_math::ac_inverse_sqrt_pwl(sqr, i_root);
      c = a*i_root;
      s = -b*i_root;
    } else {
      // Use accurate ac_math functions.
      ac_fixed<n_f_b + root_I, root_I, false> root;
      ac_math::ac_sqrt(sqr, root);
      // Instead of dividing a and -b by root, we find the inverse of root and then
      // multiply a and -b with it, so as to minimize the usage of dividers.
      ac_math::ac_div(ac_fixed<1, 1, false>(1), root, i_root);
      c = a*i_root;
      s = -b*i_root;
    }
  }

  // Perform rotation for real matrices.
  template<unsigned M, int W1, int I1, bool S1, ac_q_mode Q1, ac_o_mode O1, int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
  void rotate(
    unsigned i, unsigned j,
    ac_fixed<W1, I1, S1, Q1, O1> row1_0, ac_fixed<W1, I1, S1, Q1, O1> row1_1,
    ac_fixed<W1, I1, S1, Q1, O1> &row2_0, ac_fixed<W1, I1, S1, Q1, O1> &row2_1,
    ac_fixed<W2, I2, S2, Q2, O2> c, ac_fixed<W2, I2, S2, Q2, O2> s
  )
  {
    // Set a temporary variable to store the negative of s. Doing so and using a multiply-add
    // instead of a multiply-subtract later reduces area.
    ac_fixed<W2, I2, S2, Q2, O2> s_neg = -s;
    row2_0 = c*row1_0 + s_neg*row1_1;
    row2_1 = c*row1_1 + s*row1_0;
  }

  // Perform rotation for complex matrices.
  template<unsigned M, class T1, class T2>
  void rotate(
    unsigned i, unsigned j,
    ac_complex<T1> row1_0, ac_complex<T1> row1_1,
    ac_complex<T1> &row2_0, ac_complex<T1> &row2_1,
    ac_complex<T2> c, ac_complex<T2> s
  )
  {
    // *_R => Coefficient for R matrix multiplication
    // *_Q => Coefficient for Q matrix multiplication
    ac_complex<T2> c_R = c.conj(), c2_R = c, s_R = s, s2_R = -s.conj(), c_Q = c, s2_Q = -s, c2_Q = c.conj(), s_Q = s.conj();

    typedef typename T1::rt_unary::mag_sqr s_type; // Find type of mag_sqr value.
    enum {
      sqr_I = s_type::i_width, // Find integer width of magnitude type
      n_f_b = 32, // Arbitrarily assign 32 fractional bits for intermediate types.
      // Calculate number of integer bits for root and inverse root variables.
      root_I   = sqr_I%2 == 0 ? sqr_I/2 : (sqr_I + 1)/2,
      i_root_I = n_f_b%2 == 0 ? n_f_b/2 : (n_f_b + 1)/2,
    };
    typedef ac_fixed<n_f_b + root_I, root_I, false> T_mag;
    typedef ac_fixed<n_f_b + i_root_I, i_root_I, false> T_i_mag;

    if (i < M) {
      // Calculations for R matrix.
      row2_0 = c_R*row1_0  + s2_R*row1_1;
      row2_1 = c2_R*row1_1 + s_R*row1_0;
    } else {
      // Calculations for Q matrix.
      row2_0 = c_Q*row1_0  + s2_Q*row1_1;
      row2_1 = c2_Q*row1_1 + s_Q*row1_0;
    }
  }

  //=========================================================================
  // Off-Diagonal Processing Elements:
  // Description: Off-diagonal processing elements take the broadcasted
  // cosine and sine parameters from the diagonal processing elements, larger
  // matrix on which Givens rotation is to be performed and then systolically
  // update two rows of the matrices involved in one rotation. After the
  // processing is done, new values of row for that perticular iteration are
  // assigned back to the input matrix.
  //-------------------------------------------------------------------------

  template<unsigned M, typename T1, typename T2>
  void offdiagonal_PE (ac_matrix<T1, M, 2*M> (&A1), unsigned pivot, T2 c, T2 s, unsigned j)
  {
    T1 row1[2];
    T1 row2[2];

    OFFDIAG_PROC:
    for (unsigned i = 0; i < 2 * M; i++) {
      row1[0] = A1(pivot-1, i);
      row1[1] = A1(pivot, i);
      // Perform givens rotation through the rotate function.
      rotate<M>(i, j, row1[0], row1[1], row2[0], row2[1], c, s);
      A1(pivot - 1, i) = row2[0];
      A1(pivot, i)     = row2[1];
    }

    #if !defined(__SYNTHESIS__) && defined(AC_QRD_H_DEBUG)
    std::cout << "pivot = " << pivot << std::endl;
    std::cout << "row1 = {" << row1[0] << ", " << row1[1] << std::endl;
    std::cout << "row2 = {" << row2[0] << ", " << row2[1] << std::endl;
    #endif
  }

  //=========================================================================
  // Function : ac_qrd (ac_fixed and ac_complex<ac_fixed> implementation)
  // Description :
  // The implementation is necessary similar to the ac_fixed implementation.
  // The function takes three matrices as input, which are passed by reference.
  // First one (A) is input matrix, where as Q and R are results of decomposition,
  // which are written at memory location for them.
  // First A1 is initialized using A and a unitary matrix.
  // Then for each and every element below diagonal of input_matrix, c and s are
  // computed using diagonal processing elements.
  // These values are then propagated to off-diagonal processing elements which then
  // systolically performs Givens rotation in every iteration of algorithm and
  // modifies the input matrix.
  // Then finally, function qr_separate is called which divides the modified matrix into
  // two matrices which are final output matrices Q and R.
  //----------------------------------------------------------------------------------------------------------------

  template <bool ispwl = true, unsigned M, int W1, int I1, ac_q_mode q1, ac_o_mode o1, int W2, int I2, ac_q_mode q2, ac_o_mode o2>
  void ac_qrd (ac_matrix <ac_fixed <W1, I1, true, q1, o1>, M, M> &A, ac_matrix <ac_fixed <W2, I2, true, q2, o2>, M, M> &Q, ac_matrix <ac_fixed <W2, I2, true, q2, o2>, M, M> &R)
  {
    enum {
      n_gr = (M*(M - 1))/2, // number of given's rotations = number of elements below matrix diagonal
      add_bits = ac::log2_ceil<n_gr + 1>::val, // Each given's rotation = one extra addition per element. Number of extra bits for n_gr rotations/additions = ac_log2_ceil<n_gr + 1>
      I_imd = I1 + add_bits,
      n_f_b = 32, // Arbitrarily choose 32 fractional bits for intermediate type.
      W_imd = I_imd + n_f_b,
    };

    typedef ac_fixed <W_imd, I_imd, true, AC_RND, AC_SAT> intermediate_type;
    typedef ac_fixed <n_f_b + 2, 2, true, AC_RND, AC_SAT> sin_cosine_type;

    // Defining intermediate variables
    sin_cosine_type c, s;
    ac_matrix <intermediate_type, M, 2*M> A1;
    initialize_matrix(A, A1);

    REAL_PROC_COLUMN:
    for (unsigned column = 0; column < M; column++) {
      REAL_PROC_ROW:
      for (int row = int(M) - 1; row > 0; row--) {
        if (row == int(column)) { break; }
        diagonal_PE<ispwl>(A1(row, column), A1(row - 1, column), c, s);
        offdiagonal_PE(A1, row, c, s, column);
      }
    }
    qr_separate<true>(A1, Q, R, ac_fixed<1, 1, false>(1));
  }

  // real_diag: Make sure that all diagonal matrix elements are real, including bottom right element.
  template <bool real_diag = false, bool ispwl = true, unsigned M, int W1, int I1, ac_q_mode q1, ac_o_mode o1, int W2, int I2, ac_q_mode q2, ac_o_mode o2>
  void ac_qrd (ac_matrix < ac_complex <ac_fixed <W1, I1, true, q1, o1> >, M, M> &A, ac_matrix <ac_complex <ac_fixed <W2, I2, true, q2, o2> >, M, M> &Q, ac_matrix <ac_complex <ac_fixed <W2, I2, true, q2, o2> >, M, M> &R)
  {
    enum {
      n_gr = (M*(M - 1))/2, // number of given's rotations = number of elements below matrix diagonal
      add_bits = ac::log2_ceil<2*n_gr + 1>::val, // number of real additions = 2*number of givens rotations
      I_imd = I1 + add_bits,
      n_f_b = 32, // Arbitrarily choose 32 fractional bits for intermediate type
      W_imd = I_imd + n_f_b,
    };

    typedef ac_complex<ac_fixed <W_imd, I_imd, true, AC_RND, AC_SAT> > intermediate_type;
    typedef ac_complex<ac_fixed <n_f_b + 2, 2, true, AC_RND, AC_SAT> > sin_cosine_type;

    // Defining intermediate variables
    sin_cosine_type c, s;
    ac_matrix <intermediate_type, M, 2*M> A1;
    initialize_matrix(A, A1);

    COMPLEX_PROC_COLUMN:
    for (unsigned column = 0; column < M; column++) {
      COMPLEX_PROC_ROW:
      for (int row = int(M) - 1; row > 0; row--) {
        if (row == int(column)) { break; }
        diagonal_PE<ispwl>(A1(row,column), A1(row-1,column), c, s);
        offdiagonal_PE(A1, row, c, s, column);
      }
    }

    typedef typename intermediate_type::rt_unary::mag_sqr s_type; // Find type of mag_sqr value.
    enum {
      sqr_I = s_type::i_width, // Find integer width of magnitude type.
      // Calculate number of integer bits for root and inverse root variables.
      root_I   = sqr_I%2 == 0 ? sqr_I/2 : (sqr_I + 1)/2,
      i_root_I = n_f_b%2 == 0 ? n_f_b/2 : (n_f_b + 1)/2,
    };
    typedef ac_fixed<n_f_b + root_I, root_I, false> T_mag;
    typedef ac_fixed<n_f_b + i_root_I, i_root_I, false> T_i_mag;

    ac_complex<ac_fixed<n_f_b + 2, 2, true> > exp_arg_br_conj;

    #pragma hls_waive CNS
    if (real_diag) {
      // If all diagonals need to be real, an extra stages are required to make the bottom
      // right element real and adjust the last column of the Q matrix accordingly.
      ac_complex<ac_fixed<W_imd, I_imd, true> > br_elem = A1(M - 1, M - 1);
      T_mag mag_br;
      T_i_mag i_mag_br;
      if (ispwl) {
        // Use PWL functions
        mag_br = ac_math::ac_sqrt_pwl<T_mag>(br_elem.mag_sqr());
        i_mag_br = ac_math::ac_inverse_sqrt_pwl<T_i_mag>(br_elem.mag_sqr());
      } else {
        // Use accurate ac_math functions.
        ac_math::ac_sqrt(br_elem.mag_sqr(), mag_br);
        // Only one divider is used, and it finds the inverse of the magnitude of the
        // bottom right element. We use the ac_div function instead of the "/" operator,
        // to avoid any issues that may occur if the synthesis library doesn't have a
        // component for division.
        ac_math::ac_div(ac_fixed<1, 1, false>(1), mag_br, i_mag_br);
      }
      exp_arg_br_conj.r() = br_elem.r()*i_mag_br;
      exp_arg_br_conj.i() = br_elem.i()*i_mag_br;
      A1(M - 1, M - 1).r() = mag_br;
      A1(M - 1, M - 1).i() = 0.0;
    }

    // If all diagonal elements are required to be real, qr_separate modifies the last column
    // of the Q matrix.
    qr_separate<real_diag>(A1, Q, R, exp_arg_br_conj);
  }

  #ifdef _WIN32
  template <unsigned M, bool ispwl = true, int W1, int I1, ac_q_mode q1, ac_o_mode o1, int W2, int I2, ac_q_mode q2, ac_o_mode o2>
  #else
  template <bool ispwl = true, unsigned M, int W1, int I1, ac_q_mode q1, ac_o_mode o1, int W2, int I2, ac_q_mode q2, ac_o_mode o2>
  #endif 
  void ac_qrd (ac_fixed <W1, I1, true, q1, o1> A[M][M], ac_fixed <W2, I2, true, q2, o2> Q[M][M], ac_fixed <W2, I2, true, q2, o2> R[M][M])
  {
    ac_matrix <ac_fixed <W1, I1, true, q1, o1>, M, M> input_mat;
    ac_matrix <ac_fixed <W2, I2, true, q2, o2>, M, M> Q_mat;
    ac_matrix <ac_fixed <W2, I2, true, q2, o2>, M, M> R_mat;

    COPY_REAL_C_ARRAY_INPUT_ROW:
    for (unsigned i = 0; i < M; i++) {
      COPY_REAL_C_ARRAY_INPUT_COL:
      for (unsigned j = 0; j < M; j++) {
        input_mat (i,j) = A[i][j];
      }
    }

    ac_qrd <ispwl> (input_mat, Q_mat, R_mat);

    COPY_REAL_C_ARRAY_OUTPUT_ROW:
    for (unsigned i = 0; i < M; i++) {
      COPY_REAL_C_ARRAY_OUTPUT_COL:
      for (unsigned j = 0; j < M; j++) {
        Q [i][j] = Q_mat (i,j);
        R [i][j] = R_mat (i,j);
      }
    }

  }

  // real_diag: Make sure that all diagonal matrix elements are real, including bottom right element.
  #ifdef _WIN32
  template <unsigned M, bool real_diag = false, bool ispwl = true, int W1, int I1, ac_q_mode q1, ac_o_mode o1, int W2, int I2, ac_q_mode q2, ac_o_mode o2>
  #else
  template <bool real_diag = false, bool ispwl = true, unsigned M, int W1, int I1, ac_q_mode q1, ac_o_mode o1, int W2, int I2, ac_q_mode q2, ac_o_mode o2>
  #endif
  void ac_qrd (ac_complex <ac_fixed <W1, I1, true, q1, o1> > A[M][M], ac_complex <ac_fixed <W2, I2, true, q2, o2> > Q[M][M], ac_complex <ac_fixed <W2, I2, true, q2, o2> > R[M][M])
  {
    ac_matrix <ac_complex <ac_fixed <W1, I1, true, q1, o1> >, M, M> input_mat;
    ac_matrix <ac_complex <ac_fixed <W2, I2, true, q2, o2> >, M, M> Q_mat;
    ac_matrix <ac_complex <ac_fixed <W2, I2, true, q2, o2> >, M, M> R_mat;

    COPY_COMPLEX_C_ARRAY_INPUT_ROW:
    for (unsigned i = 0; i < M; i++) {
      COPY_COMPLEX_C_ARRAY_INPUT_COL:
      for (unsigned j = 0; j < M; j++) {
        input_mat(i,j) = A[i][j];
      }
    }

    ac_qrd <real_diag, ispwl> (input_mat, Q_mat, R_mat);

    COPY_COMPLEX_C_ARRAY_OUTPUT_ROW:
    for (unsigned i = 0; i < M; i++) {
      COPY_COMPLEX_C_ARRAY_OUTPUT_COL:
      for (unsigned j = 0; j < M; j++) {
        Q[i][j] = Q_mat(i,j);
        R[i][j] = R_mat(i,j);
      }
    }

  }
}

#endif
