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
// *******************************************************************************************
// File: ac_chol_d.h
//
// Description: Provides an implementation for the Cholesky Decomposition of ac_fixed
//    and ac_complex<ac_fixed> matrices.
//
// Usage:
//    (see the type-specific examples below)
//
// Notes:
//    This file uses C++ function overloading to target implementations specific to each type
//    of data. Attempting to call the function with a type that is not implemented will result
//    in a compile error.
//
//    The user may choose to use the PWL functions for the internal calculations or they may
//    choose to use the accurate div and sqrt functions provided in the ac_math library,
//    using an optional template parameter. Please refer to the documentation for more details.
//
//    This library uses the ac_inverse_sqrt_pwl(), ac_sqrt_pwl(), ac_div() and ac_sqrt()
//    functions from the ac_math header files.
//
// Revision History:
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//    3.3.0  - [CAT-25797] Added CDesignChecker fixes/waivers for code check violations in ac_math PWL and Linear Algebra IPs.
//             Waivers added for CNS, CCC, ABR and ABW violations.
//             FXD violations fixed by changing integer literals to floating point literals.
//    2.0.10 - Official open-source release as part of the ac_math library.
//    Niramay Sanghvi : Nov 24 2017 : Used friend function to handle ac_matrix inputs/outputs.
//    Niramay Sanghvi : Nov 12 2017 : Added overloaded functions to handle standard C arrays.
//    Niramay Sanghvi : Oct 02 2017 : Incorporated the use of the inverse_sqrt PWL function.
//    Niramay Sanghvi : Aug 10 2017 : Added template parameters for precision configuration.
//    Niramay Sanghvi : Aug 09 2017 : Made output go to zero for non-positive def. matrix.
//    Niramay Sanghvi : Jul 28 2017 : Added choice between PWL and accurate functions.
//
// *******************************************************************************************

#ifndef _INCLUDED_AC_CHOL_D_H_
#define _INCLUDED_AC_CHOL_D_H_

#if (defined(__GNUC__) && (__cplusplus < 201103L))
#error Please use C++11 or a later standard for compilation.
#endif
#if (defined(_MSC_VER) && (_MSC_VER < 1920) && !defined(__EDG__))
#error Please use Microsoft VS 2019 or a later standard for compilation.
#endif


// Include headers for data types supported by these implementations
#include <ac_fixed.h>
#include <ac_float.h>
#include <ac_std_float.h>
#include <ac_complex.h>
#include <ac_matrix.h>

// Include headers for required functions
#include <ac_math/ac_sqrt_pwl.h>
#include <ac_math/ac_inverse_sqrt_pwl.h>
#include <ac_math/ac_div.h>
#include <ac_math/ac_sqrt.h>

#if !defined(__SYNTHESIS__) && defined(AC_CHOL_D_H_DEBUG)
#include <iostream>
#endif

// =========================================================================
// Function: ac_chol_d (for ac_fixed using native C-style arrays)
//
// Description:
//    Calculation of the Cholesky Decomposition of real-valued matrices of
//    ac_fixed variables.
//
//    The Cholesky-Crout algorithm is used for Cholesky Decomposition.
//    By default, all temporary variables use the same precision as the
//    output. The user can change that by specifying number of bits to be
//    added or taken away from this default value. The user can also
//    add extra template parameters to change the rounding and saturation
//    modes of the intermediate variables (by default, the rounding mode is
//    AC_RND and saturation mode is AC_SAT).
//
//    The user can also choose between using PWL vs. using accurate
//    division/sqrt functions, as mentioned earlier.
//
// Usage:
//    A sample implementation and its testbench looks like this:
//
//    #include <ac_fixed.h>
//    #include <ac_math/ac_chol_d.h>
//    using namespace ac_math;
//
//    // Define data types for input and output matrices
//    typedef ac_fixed<41, 21, true, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<64, 32, true, AC_RND, AC_SAT> output_type;
//    const unsigned M = 7;
//
//    #pragma hls_design top
//    void project(
//      const input_type input[M][M],
//      output_type output[M][M]
//    )
//    {
//      ac_chol_d(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input[7][7];
//      output_type output[7][7];
//      // Hypothetical function that generates a pos. def. matrix:
//      gen_pos_def_matrix(input);
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    A runtime error is thrown by AC_ASSERT if the input matrix is not
//    positive definite. Additionally the function will return a matrix
//    of all zeros if the input matrix is not positive definite, in case
//    the AC_ASSERT cannot provide an indication of that.
//
// -------------------------------------------------------------------------

namespace ac_math
{
  template<class T_out, int W1>
  T_out ind_conv(const ac_int<W1, false> &i, const ac_int<W1, false> &j)
  {
    // Modified index = (i*(i + 1)/2) + j
    T_out index = ((i*(i + 1))>>1) + j;
    return index;
  }
  #ifdef _WIN32
  template<unsigned M, bool use_pwl = false,
           int delta_w = 0, int delta_i = 0, ac_q_mode imd_Q = AC_RND, ac_o_mode imd_O = AC_SAT,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO>  
  #else
  template<bool use_pwl = false,
         int delta_w = 0, int delta_i = 0, ac_q_mode imd_Q = AC_RND, ac_o_mode imd_O = AC_SAT,
         int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
         int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO,
         unsigned M>
  #endif
  void ac_chol_d(
    const ac_fixed<W, I, S, Q, O> A[M][M],
    ac_fixed<outW, outI, outS, outQ, outO> L[M][M]
  )
  {
    // Make this go false if the input matrix is not positive definite (checking is done later)
    bool pos_def = true;

    typedef ac_fixed<outW, outI, outS, outQ, outO> T_out;
    // The number of non-zero elements in a square lower triangular matrix is (M*(M + 1)/2).
    // Hence, the array for storing the results of intermediate calculations should also be of that size.
    // The size of the resulting array is roughly half of what it would've been if the area was 2D square (M^2)
    // In order to map the indexing operations from a square matrix to this lower triangular matrix, the
    // intermediate matrix is one-dimensional with the above size, and any indexing operation for location
    // (i, j) is modified for the 1D array access by using the formula in the ind_conv() function.
    enum { L_IMD_SIZE = (M*(M + 1))/2 };
    typedef ac_int<ac::nbits<M>::val, false> index_type;
    typedef ac_int<ac::nbits<L_IMD_SIZE - 1>::val, false> t_1D;
    T_out L_imd_1D[L_IMD_SIZE];
    // Define type for the intermediate variables
    // Add an extra bit to W and I for intermediate variable type if the output is unsigned, and make sure that i_s_t is signed.
    typedef class ac_fixed<outW + delta_w, outI + delta_i, true, imd_Q, imd_O> i_s_t;
    // Unsigned versions of i_s_t and T_out
    typedef ac_fixed<outW, outI, false, outQ, outO> T_out_u;
    typedef ac_fixed<i_s_t::width, i_s_t::i_width, false, i_s_t::q_mode, i_s_t::o_mode> i_s_t_u;

    ARRAY_AC_FIXED_L_COL:
    for (index_type j = 0; j < M; j++) {
      i_s_t sum_Ajj_Ljk_sq = A[j][j];
      ARRAY_AC_FIXED_LJJ_K:
      for (index_type k = 0; k < M; k++) {
        // Break statement is used to avoid incorrect accesses to L_imd_1D.
        if (k == j) { break; }
        T_out L_imd_1D_temp = L_imd_1D[ind_conv<t_1D>(j, k)];
        sum_Ajj_Ljk_sq -= L_imd_1D_temp * L_imd_1D_temp;
      }

      // Use a macro to activate the AC_ASSERT
      #ifdef ASSERT_ON_INVALID_INPUT
      // Check to make sure that the matrix is positive definite. If "sum_Ajj_Ljk_sq" is negative/zero, then the diagonal
      // element will be complex/infinite, which is not valid. This condition will not be encountered if the
      // input matrix is positive definite
      AC_ASSERT(sum_Ajj_Ljk_sq > 0, "Input matrix is not positive definite");
      #endif
      if (sum_Ajj_Ljk_sq <= 0) {pos_def = false;}
      i_s_t_u recip_Ljj;
      T_out_u imd_Ljj;
      // Compute values for and initialize diagonal elements using PWL/accurate math functions, as may be the case.
      #pragma hls_waive CNS
      if (use_pwl) {
        // Use the PWL functions
        ac_math::ac_sqrt_pwl((i_s_t_u) sum_Ajj_Ljk_sq, imd_Ljj);
        L_imd_1D[ind_conv<t_1D>(j, j)] = imd_Ljj;
        // Store inverse of diagonal element in separate variable (i.e. "recip_Ljj") for later calculations.
        ac_math::ac_inverse_sqrt_pwl((i_s_t_u) sum_Ajj_Ljk_sq, recip_Ljj);
      } else {
        // Use accurate math functions.
        ac_math::ac_sqrt((i_s_t_u)sum_Ajj_Ljk_sq, imd_Ljj);
        L_imd_1D[ind_conv<t_1D>(j, j)] = imd_Ljj;
        // Make sure that every variable to be passed to the div function has the same sign.
        const ac_fixed<1, 1, false> unity = 1.0;
        // Store inverse of diagonal element in separate variable (i.e. "recip_Ljj") for later calculations.
        #pragma hls_waive DBZ
        ac_math::ac_div(unity, imd_Ljj, recip_Ljj);
      }
      #if !defined(__SYNTHESIS__) && defined(AC_CHOL_D_H_DEBUG)
      std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
      std::cout << "sum_Ajj_Ljk_sq = " << sum_Ajj_Ljk_sq << std::endl;
      std::cout << "imd_Ljj        = " << imd_Ljj << std::endl;
      std::cout << "recip_Ljj      = " << recip_Ljj << std::endl;
      #endif

      // Initializing non-diagonal elements.
      ARRAY_AC_FIXED_L_ROW:
      for (index_type i = M - 1; i > 0; i--) {
        i_s_t sum_Aij_Lik_Ljk;
        sum_Aij_Lik_Ljk = A[i][j];
        ARRAY_AC_FIXED_LIJ_K:
        for (index_type k = 0; k < M; k++) {
          // Break statement is used to avoid incorrect accesses to L_imd_1D.
          if (k == j) { break; }
          sum_Aij_Lik_Ljk -= L_imd_1D[ind_conv<t_1D>(i, k)] * L_imd_1D[ind_conv<t_1D>(j, k)];
        }
        // Break statement is used to avoid incorrect accesses to L_imd_1D.
        if (i == j) { break; }
        L_imd_1D[ind_conv<t_1D>(i, j)] = (T_out)(recip_Ljj * sum_Aij_Lik_Ljk);
      }
    }

    ARRAY_AC_FIXED_OUTPUT_COPY_ROW:
    for (index_type i = 0; i < M; i++) {
      ARRAY_AC_FIXED_SET_OUTPUT_COPY_COL:
      for (index_type j = 0; j < M; j++) {
        // If L is not positive definite, all elements should be initialized to zero.
        if (j > i || !pos_def) { L[i][j] = 0.0; }
        // If L is positive definite, only initialize elements above the diagonal to zero.
        else {
          L[i][j] = L_imd_1D[ind_conv<t_1D>(i, j)];
        }
      }
    }
  }

// =========================================================================
// Function: ac_chol_d (for ac_float using native C-style arrays)
//
// Description:
//    Calculation of the Cholesky Decomposition of real-valued matrices of
//    ac_float variables.
//
//    The Cholesky-Crout algorithm is used for Cholesky Decomposition.
//    By default, all temporary variables use the same precision as the
//    output. The user can change that by specifying number of bits to be
//    added or taken away from this default value. The user can also
//    add extra template parameters to change the rounding mode of the
//    intermediate variables (AC_RND by default).
//
//    The user can also choose between using PWL vs. using accurate
//    division/sqrt functions, as mentioned earlier.
//
// Usage:
//    A sample implementation and its testbench looks like this:
//
//    #include <ac_float.h>
//    #include <ac_math/ac_chol_d.h>
//    using namespace ac_math;
//
//    // Define data types for input and output matrices
//    typedef ac_float<18, 2, 10, AC_TRN> input_type;
//    typedef ac_float<18, 2, 10, AC_TRN> output_type;
//    const unsigned M = 7;
//
//    #pragma hls_design top
//    void project(
//      const input_type input[M][M],
//      output_type output[M][M]
//    )
//    {
//      ac_chol_d(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input[7][7];
//      output_type output[7][7];
//      // Hypothetical function that generates a pos. def. matrix:
//      gen_pos_def_matrix(input);
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    A runtime error is thrown by AC_ASSERT if the input matrix is not
//    positive definite. Additionally the function will return a matrix
//    of all zeros if the input matrix is not positive definite, in case
//    the AC_ASSERT cannot provide an indication of that.
//
// -------------------------------------------------------------------------
  #ifdef _WIN32
  template<unsigned M, bool use_pwl = false,
           int delta_w = 0, int delta_i = 0, int delta_e = 0, ac_q_mode imd_Q = AC_RND,
           int W, int I, int E, ac_q_mode Q,
           int outW, int outI, int outE, ac_q_mode outQ>  
  #else
  template<bool use_pwl = false,
           int delta_w = 0, int delta_i = 0, int delta_e = 0, ac_q_mode imd_Q = AC_RND,
           int W, int I, int E, ac_q_mode Q,
           int outW, int outI, int outE, ac_q_mode outQ,
           unsigned M>
  #endif
  void ac_chol_d(
    const ac_float<W, I, E, Q> A[M][M],
    ac_float<outW, outI, outE, outQ> L[M][M]
  )
  {
    // Make this go false if the input matrix is not positive definite (checking is done later)
    bool pos_def = true;

    typedef ac_float<outW, outI, outE, outQ> T_out;
    // The number of non-zero elements in a square lower triangular matrix is (M*(M + 1)/2).
    // Hence, the array for storing the results of intermediate calculations should also be of that size.
    // The size of the resulting array is roughly half of what it would've been if the area was 2D square (M^2)
    // In order to map the indexing operations from a square matrix to this lower triangular matrix, the
    // intermediate matrix is one-dimensional with the above size, and any indexing operation for location
    // (i, j) is modified for the 1D array access by using the formula in the ind_conv() function.
    enum { L_IMD_SIZE = (M*(M + 1))/2 };
    typedef ac_int<ac::nbits<M>::val, false> index_type;
    typedef ac_int<ac::nbits<L_IMD_SIZE - 1>::val, false> t_1D;
    T_out L_imd_1D[L_IMD_SIZE];
    // Define type for the intermediate variables
    typedef class ac_float<outW + delta_w, outI + delta_i, outE + delta_e, imd_Q> T_imd;

    ARRAY_AC_FLOAT_L_COL:
    for (index_type j = 0; j < M; j++) {
      T_imd sum_Ajj_Ljk_sq = A[j][j];
      ARRAY_AC_FLOAT_LJJ_K:
      for (index_type k = 0; k < M; k++) {
        // Break statement is used to avoid incorrect accesses to L_imd_1D.
        if (k == j) { break; }
        T_out L_imd_1D_temp = L_imd_1D[ind_conv<t_1D>(j, k)];
        sum_Ajj_Ljk_sq -= L_imd_1D_temp * L_imd_1D_temp;
      }

      // Use a macro to activate the AC_ASSERT
      #ifdef ASSERT_ON_INVALID_INPUT
      // Check to make sure that the matrix is positive definite. If "sum_Ajj_Ljk_sq" is negative/zero, then the diagonal
      // element will be complex/infinite, which is not valid. This condition will not be encountered if the
      // input matrix is positive definite
      AC_ASSERT(sum_Ajj_Ljk_sq > 0.0, "Input matrix is not positive definite");
      #endif
      if (sum_Ajj_Ljk_sq <= 0) { pos_def = false; }
      T_imd recip_Ljj;
      T_out imd_Ljj;
      // Compute values for and initialize diagonal elements using PWL/accurate math functions, as may be the case.
      #pragma hls_waive CNS
      if (use_pwl) {
        // Use the PWL functions
        ac_math::ac_sqrt_pwl(sum_Ajj_Ljk_sq, imd_Ljj);
        L_imd_1D[ind_conv<t_1D>(j, j)] = imd_Ljj;
        // Store inverse of diagonal element in separate variable (i.e. "recip_Ljj") for later calculations.
        ac_math::ac_inverse_sqrt_pwl(sum_Ajj_Ljk_sq, recip_Ljj);
      } else {
        // Use accurate math functions.
        ac_math::ac_sqrt(sum_Ajj_Ljk_sq, imd_Ljj);
        L_imd_1D[ind_conv<t_1D>(j, j)] = imd_Ljj;
        const ac_float<2, 2, 2> unity = 1.0;
        // Store inverse of diagonal element in separate variable (i.e. "recip_Ljj") for later calculations.
        #pragma hls_waive DBZ
        ac_math::ac_div(unity, imd_Ljj, recip_Ljj);
      }
      #if !defined(__SYNTHESIS__) && defined(AC_CHOL_D_H_DEBUG)
      std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
      std::cout << "sum_Ajj_Ljk_sq = " << sum_Ajj_Ljk_sq << std::endl;
      std::cout << "imd_Ljj        = " << imd_Ljj << std::endl;
      std::cout << "recip_Ljj      = " << recip_Ljj << std::endl;
      #endif

      // Initializing non-diagonal elements.
      ARRAY_AC_FLOAT_L_ROW:
      for (index_type i = M - 1; i > 0; i--) {
        T_imd sum_Aij_Lik_Ljk = A[i][j];
        ARRAY_AC_FLOAT_LIJ_K:
        for (index_type k = 0; k < M; k++) {
          // Break statement is used to avoid incorrect accesses to L_imd_1D.
          if (k == j) { break; }
          sum_Aij_Lik_Ljk -= L_imd_1D[ind_conv<t_1D>(i, k)] * L_imd_1D[ind_conv<t_1D>(j, k)];
        }
        // Break statement is used to avoid incorrect accesses to L_imd_1D.
        if (i == j) { break; }
        L_imd_1D[ind_conv<t_1D>(i, j)] = (T_out)(recip_Ljj * sum_Aij_Lik_Ljk);
      }
    }

    ARRAY_AC_FLOAT_OUTPUT_COPY_ROW:
    for (index_type i = 0; i < M; i++) {
      ARRAY_AC_FLOAT_SET_OUTPUT_COPY_COL:
      for (index_type j = 0; j < M; j++) {
        // If L is not positive definite, all elements should be initialized to zero.
        if (j > i || !pos_def) { L[i][j] = 0.0; }
        // If L is positive definite, only initialize elements above the diagonal to zero.
        else {
          L[i][j] = L_imd_1D[ind_conv<t_1D>(i, j)];
        }
      }
    }
  }

// =========================================================================
// Function: ac_chol_d (for ac_std_float using native C-style arrays)
//
// Description:
//    Calculation of the Cholesky Decomposition of real-valued matrices of
//    ac_std_float variables.
//
//    The Cholesky-Crout algorithm is used for Cholesky Decomposition.
//    By default, all temporary variables use the same precision as the
//    output. The user can change that by specifying number of bits to be
//    added or taken away from this default value. The user can also
//    add extra template parameters to change the rounding mode of the
//    intermediate variables (AC_RND by default).
//
//    The user can also choose between using PWL vs. using accurate
//    division/sqrt functions, as mentioned earlier.
//
//    This version depends on the ac_float implementation for computation.
//
// Usage:
//    A sample implementation and its testbench looks like this:
//
//    #include <ac_std_float.h>
//    #include <ac_math/ac_chol_d.h>
//    using namespace ac_math;
//
//    // Define data types for input and output matrices
//    typedef ac_std_float<32, 8> input_type;
//    typedef ac_std_float<32, 8> output_type;
//    const unsigned M = 7;
//
//    #pragma hls_design top
//    void project(
//      const input_type input[M][M],
//      output_type output[M][M]
//    )
//    {
//      ac_chol_d(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input[7][7];
//      output_type output[7][7];
//      // Hypothetical function that generates a pos. def. matrix:
//      gen_pos_def_matrix(input);
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    A runtime error is thrown by AC_ASSERT if the input matrix is not
//    positive definite. Additionally the function will return a matrix
//    of all zeros if the input matrix is not positive definite, in case
//    the AC_ASSERT cannot provide an indication of that.
//
// -------------------------------------------------------------------------
  #ifdef _WIN32
  template<unsigned M, bool use_pwl = false,
           int delta_w = 0, int delta_i = 0, int delta_e = 0, ac_q_mode imd_Q = AC_RND,
           int W, int E,
           int outW, int outE>  
  #else
  template<bool use_pwl = false,
           int delta_w = 0, int delta_i = 0, int delta_e = 0, ac_q_mode imd_Q = AC_RND,
           int W, int E,
           int outW, int outE,
           unsigned M>
  #endif
  void ac_chol_d(
    const ac_std_float<W, E> A[M][M],
    ac_std_float<outW, outE> L[M][M]
  )
  {
    // Intermediate variables that enable interfacing with ac_float version.
    ac_float<W - E + 1, 2, E> A_ac_fl[M][M];
    ac_float<outW - outE + 1, 2, outE> L_ac_fl[M][M];

    ARRAY_AC_STD_FLOAT_INPUT_COPY_ROW:
    for (int i = 0; i < M; i++) {
      ARRAY_AC_STD_FLOAT_INPUT_COPY_COL:
      for (int j = 0; j < M; j++) {
        A_ac_fl[i][j] = A[i][j].to_ac_float();
      }
    }


  #ifdef _WIN32 
    ac_chol_d<M, use_pwl, delta_w, delta_i, delta_e, imd_Q>(A_ac_fl, L_ac_fl);
  #else
    ac_chol_d<use_pwl, delta_w, delta_i, delta_e, imd_Q>(A_ac_fl, L_ac_fl);
  #endif
    ARRAY_AC_STD_FLOAT_OUTPUT_COPY_ROW:
    for (int i = 0; i < M; i++) {
      ARRAY_AC_STD_FLOAT_OUTPUT_COPY_COL:
      for (int j = 0; j < M; j++) {
        ac_std_float<outW, outE> L_val_temp(L_ac_fl[i][j]);
        L[i][j] = L_val_temp;
      }
    }
  }

// =========================================================================
// Function: ac_chol_d (for ac_ieee_float using native C-style arrays)
//
// Description:
//    Calculation of the Cholesky Decomposition of real-valued matrices of
//    ac_ieee_float variables.
//
//    The Cholesky-Crout algorithm is used for Cholesky Decomposition.
//    By default, all temporary variables use the same precision as the
//    output. The user can change that by specifying number of bits to be
//    added or taken away from this default value. The user can also
//    add extra template parameters to change the rounding mode of the
//    intermediate variables (AC_RND by default).
//
//    The user can also choose between using PWL vs. using accurate
//    division/sqrt functions, as mentioned earlier.
//
//    This version depends on the ac_float implementation for computation.
//
// Usage:
//    A sample implementation and its testbench looks like this:
//
//    #include <ac_ieee_float.h>
//    #include <ac_math/ac_chol_d.h>
//    using namespace ac_math;
//
//    // Define data types for input and output matrices
//    typedef ac_ieee_float<binary32> input_type;
//    typedef ac_ieee_float<binary32> output_type;
//    const unsigned M = 7;
//
//    #pragma hls_design top
//    void project(
//      const input_type input[M][M],
//      output_type output[M][M]
//    )
//    {
//      ac_chol_d(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input[7][7];
//      output_type output[7][7];
//      // Hypothetical function that generates a pos. def. matrix:
//      gen_pos_def_matrix(input);
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    A runtime error is thrown by AC_ASSERT if the input matrix is not
//    positive definite. Additionally the function will return a matrix
//    of all zeros if the input matrix is not positive definite, in case
//    the AC_ASSERT cannot provide an indication of that.
//
// -------------------------------------------------------------------------

  #ifdef _WIN32
  template<unsigned M, bool use_pwl = false,
           int delta_w = 0, int delta_i = 0, int delta_e = 0, ac_q_mode imd_Q = AC_RND,
           ac_ieee_float_format Format,
           ac_ieee_float_format outFormat>  
  #else
  template<bool use_pwl = false,
           int delta_w = 0, int delta_i = 0, int delta_e = 0, ac_q_mode imd_Q = AC_RND,
           ac_ieee_float_format Format,
           ac_ieee_float_format outFormat,
           unsigned M>
  #endif
  void ac_chol_d(
    const ac_ieee_float<Format> A[M][M],
    ac_ieee_float<outFormat> L[M][M]
  )
  {
    typedef ac_ieee_float<Format> T_in;
    const int W = T_in::width;
    const int E = T_in::e_width;
    typedef ac_ieee_float<outFormat> T_out;
    const int outW = T_out::width;
    const int outE = T_out::e_width;
    // Intermediate variables that enable interfacing with ac_float version.
    ac_float<W - E + 1, 2, E> A_ac_fl[M][M];
    ac_float<outW - outE + 1, 2, outE> L_ac_fl[M][M];

    ARRAY_AC_IEEE_FLOAT_INPUT_COPY_ROW:
    for (int i = 0; i < M; i++) {
      ARRAY_AC_IEEE_FLOAT_INPUT_COPY_COL:
      for (int j = 0; j < M; j++) {
        A_ac_fl[i][j] = A[i][j].to_ac_float();
      }
    }

#ifdef _WIN32
    ac_chol_d<M, use_pwl, delta_w, delta_i, delta_e, imd_Q>(A_ac_fl, L_ac_fl);
#else
    ac_chol_d<use_pwl, delta_w, delta_i, delta_e, imd_Q>(A_ac_fl, L_ac_fl);
#endif

    ARRAY_AC_IEEE_FLOAT_OUTPUT_COPY_ROW:
    for (int i = 0; i < M; i++) {
      ARRAY_AC_IEEE_FLOAT_OUTPUT_COPY_COL:
      for (int j = 0; j < M; j++) {
        T_out L_val_temp(L_ac_fl[i][j]);
        L[i][j] = L_val_temp;
      }
    }
  }

// ==============================================================================================
// Function: ac_chol_d (for ac_complex<ac_fixed> using native C-style arrays)
//
// Description:
//    Calculation of the Cholesky Decomposition of complex-valued matrices
//    with ac_complex<ac_fixed> variables.
//
//    This function can also be configured to change the type for the
//    intermediate variables, and to use the pwl/accurate math functions for
//    intermediate computations, just like the ac_fixed version.
//
// Usage:
//    A sample implementation and its testbench looks like this:
//
//    #include <ac_fixed.h>
//    #include <ac_math/ac_chol_d.h>
//    using namespace ac_math;
//
//    // Define data types for input and output matrices
//    typedef ac_complex<ac_fixed<41, 21, true, AC_RND, AC_SAT> > input_type;
//    typedef ac_complex<ac_fixed<64, 32, true, AC_RND, AC_SAT> > output_type;
//    const unsigned M = 7;
//
//    #pragma hls_design top
//    void project(
//      const input_type input[M][M],
//      output_type output[M][M]
//    )
//    {
//      ac_chol_d(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input[M][M];
//      output_type output[M][M];
//      gen_pos_def_matrix(input); // Hypothetical function that generates a pos. def. matrix
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    Similar to the ac_fixed version, the ac_complex version also has an AC_ASSERT and a backup
//    functionality to return a matrix of zeros if the input matrix is not positive definite.
//
//    As the diagonal elements of the output matrix are always real, the calculations involved
//    are optimized to ensure that only the real part of diagonal elements is ever calculated,
//    and only the real part is used for future calculations. The imaginary part is always set
//    to zero.
//
// ----------------------------------------------------------------------------------------------

  #ifdef _WIN32
  template<unsigned M, bool use_pwl = false,
           int delta_w = 0, int delta_i = 0, ac_q_mode imd_Q = AC_RND, ac_o_mode imd_O = AC_SAT,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO>  
  #else
  template<bool use_pwl = false,
         int delta_w = 0, int delta_i = 0, ac_q_mode imd_Q = AC_RND, ac_o_mode imd_O = AC_SAT,
         int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
         int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO,
         unsigned M>
  #endif
  void ac_chol_d(
    const ac_complex<ac_fixed<W, I, S, Q, O> > A[M][M],
    ac_complex<ac_fixed<outW, outI, outS, outQ, outO> > L[M][M]
  )
  {
    // Make this go false if the input matrix is not positive definite (checking is done later)
    bool pos_def = true;

    typedef ac_fixed<outW, outI, outS, outQ, outO> T_out;
    // The number of non-zero elements in a square lower triangular matrix is (M*(M + 1)/2).
    // Hence, the array for storing the results of intermediate calculations should also be of that size.
    // The size of the resulting array is roughly half of what it would've been if the area was 2D square (M^2)
    // In order to map the indexing operations from a square matrix to this lower triangular matrix, the
    // intermediate matrix is one-dimensional with the above size, and any indexing operation for location
    // (i, j) is modified for the 1D array access by using the formula in the ind_conv() function.
    enum { L_IMD_SIZE = (M*(M + 1))/2 };
    typedef ac_int<ac::nbits<M>::val, false> index_type;
    typedef ac_int<ac::nbits<L_IMD_SIZE - 1>::val, false> t_1D;
    ac_complex<T_out> L_imd_1D[L_IMD_SIZE];
    // Define type for the intermediate variables
    // Add an extra bit to W and I for intermediate variable type if the output is unsigned, and make sure that i_s_t is signed.
    typedef ac_fixed<outW + delta_w, outI + delta_i, true, imd_Q, imd_O> i_s_t;
    typedef ac_fixed<outW, outI, false, outQ, outO> T_out_fixed_u;
    typedef ac_fixed<i_s_t::width, i_s_t::i_width, false, i_s_t::q_mode, i_s_t::o_mode> i_s_t_u;

    ARRAY_AC_COMP_AC_FIXED_L_COL:
    for (index_type j = 0; j < M; j++) {
      i_s_t sum_Ajj_Ljk_sq;
      #pragma hls_waive ABR
      sum_Ajj_Ljk_sq = A[j][j].r();
      ARRAY_AC_COMP_AC_FIXED_LJJ_K:
      for (index_type k = 0; k < M; k++) {
        // Break statement is used to avoid incorrect accesses to L_imd_1D.
        if (k == j) { break; }
        sum_Ajj_Ljk_sq -= L_imd_1D[ind_conv<t_1D>(j, k)].mag_sqr();
      }

      // Use a macro to activate the AC_ASSERT
      #ifdef ASSERT_ON_INVALID_INPUT
      // Check to make sure that the input matrix is positive definite. If "sum_Ajj_Ljk_sq" is negative/zero, then
      // the diagonal element will be complex/infinite, which is not valid. This condition will not be encountered if the
      // input matrix is positive definite
      AC_ASSERT(sum_Ajj_Ljk_sq > 0, "Input matrix is not positive definite");
      #endif

      if (sum_Ajj_Ljk_sq <= 0) { pos_def = false; }

      i_s_t_u recip_Ljj;
      T_out_fixed_u imd_Ljj;

      // Compute values for and initialize diagonal elements using PWL/accurate math functions, as may be the case.
      #pragma hls_waive CNS
      if (use_pwl) {
        // Use the PWL functions
        // Initialize diagonal elements. Since the diagonal elements are real, initialize the imaginary part to 0.
        // Only bother with the real part, as the diagonal elements of the decomposed matrix are always real.
        ac_math::ac_sqrt_pwl((i_s_t_u) sum_Ajj_Ljk_sq, imd_Ljj);
        t_1D act_ind = ind_conv<t_1D>(j, j);
        #pragma hls_waive ABW
        L_imd_1D[act_ind].r() = imd_Ljj;
        #pragma hls_waive ABW
        L_imd_1D[act_ind].i() = 0.0;
        // Store inverse of real part of diagonal element in separate variable (i.e. "recip_Ljj") for later calculations.
        ac_math::ac_inverse_sqrt_pwl((i_s_t_u) sum_Ajj_Ljk_sq, recip_Ljj);
      } else {
        // Use accurate math functions.
        // Only bother with the real part, as the diagonal elements of the decomposed matrix are always real.
        ac_math::ac_sqrt((i_s_t_u)sum_Ajj_Ljk_sq, imd_Ljj);
        t_1D act_ind = ind_conv<t_1D>(j, j);
        #pragma hls_waive ABW
        L_imd_1D[act_ind].r() = imd_Ljj;
        #pragma hls_waive ABW
        L_imd_1D[act_ind].i() = 0.0;
        // Make sure that every variable to be passed to the div function has the same sign.
        const ac_fixed<1, 1, false> unity = 1.0;
        // Store inverse of diagonal element in separate variable (i.e. "recip_Ljj") for later calculations.
        #pragma hls_waive DBZ
        ac_math::ac_div(unity, imd_Ljj, recip_Ljj);
      }
      #if !defined(__SYNTHESIS__) && defined(AC_CHOL_D_H_DEBUG)
      std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
      std::cout << "j = " << j << std::endl;
      std::cout << "sum_Ajj_Ljk_sq = " << sum_Ajj_Ljk_sq << std::endl;
      std::cout << "imd_Ljj        = " << imd_Ljj << std::endl;
      std::cout << "recip_Ljj      = " << recip_Ljj << std::endl;
      #endif

      // Initializing non-diagonal elements.
      ARRAY_AC_COMP_AC_FIXED_L_ROW:
      for (index_type i = M - 1; i > 0; i--) {
        ac_complex<i_s_t> sum_Aij_Lik_Ljk = A[i][j];
        ARRAY_AC_COMP_AC_FIXED_LIJ_K:
        for (index_type k = 0; k < M; k++) {
          // Break statement is used to avoid incorrect accesses to L_imd_1D.
          if (k == j) { break; }
          sum_Aij_Lik_Ljk -= L_imd_1D[ind_conv<t_1D>(i, k)] * L_imd_1D[ind_conv<t_1D>(j, k)].conj();
        }
        // Break statement is used to avoid incorrect accesses to L_imd_1D.
        if (i == j) { break; }
        t_1D act_ind = ind_conv<t_1D>(i, j);
        #pragma hls_waive ABW
        L_imd_1D[act_ind].r() = (T_out)(recip_Ljj * sum_Aij_Lik_Ljk.r());
        #pragma hls_waive ABW
        L_imd_1D[act_ind].i() = (T_out)(recip_Ljj * sum_Aij_Lik_Ljk.i());
      }
    }

    ARRAY_AC_COMP_AC_FIXED_SET_OUTPUT_COPY_ROW:
    for (index_type i = 0; i < M; i++) {
      ARRAY_AC_COMP_AC_FIXED_SET_OUTPUT_COPY_COL:
      for (index_type j = 0; j < M; j++) {
        // If input matrix is not positive definite, all elements should be initialized to zero.
        // If input matrix is positive definite, only initialize elements above the diagonal to zero.
        if (j > i || !pos_def) {
          #pragma hls_waive ABW
          L[i][j].r() = 0.0;
          #pragma hls_waive ABW
          L[i][j].i() = 0.0;
        } else { L[i][j] = L_imd_1D[ind_conv<t_1D>(i, j)]; }
      }
    }

  }
} // namespace ac_math

// =============================================================================================
// Function: indirect_chol_d
// Helper function for using ac_chol_d on ac_matrix<ac_fixed/ac_complex<ac_fixed>> objects

template<bool use_pwl = false,
         int delta_w = 0, int delta_i = 0, ac_q_mode imd_Q = AC_RND, ac_o_mode imd_O = AC_SAT,
         class T1, unsigned M1, class T2>
void indirect_chol_d(const ac_matrix<T1, M1, M1> &input, ac_matrix<T2, M1, M1> &output)
{
  // Extract 2D array member data, and pass it over to the 2D array implementation.
#ifdef _WIN32
  ac_math::ac_chol_d<M1, use_pwl, delta_w, delta_i, imd_Q, imd_O>(input.m_data, output.m_data);
#else
  ac_math::ac_chol_d<use_pwl, delta_w, delta_i, imd_Q, imd_O>(input.m_data, output.m_data);
#endif

}

// =============================================================================================
// Function: indirect_chol_d_float
// Helper function for using ac_chol_d on ac_matrix<ac_float/ac_std_float/ac_ieee_float> objects

template<bool use_pwl = false,
         int delta_w = 0, int delta_i = 0, ac_q_mode imd_Q = AC_RND,
         class T1, unsigned M1, class T2>
void indirect_chol_d_float(const ac_matrix<T1, M1, M1> &input, ac_matrix<T2, M1, M1> &output)
{
  // Extract 2D array member data, and pass it over to the 2D array implementation.
#ifdef _WIN32
  ac_math::ac_chol_d<M1, use_pwl, delta_w, delta_i, imd_Q>(input.m_data, output.m_data);
#else
  ac_math::ac_chol_d<use_pwl, delta_w, delta_i, 0, imd_Q>(input.m_data, output.m_data);
#endif
}

namespace ac_math
{
// ===============================================================================
// Function: ac_chol_d (for ac_fixed using ac_matrix 2-D storage class)
//
// Description:
//    Calculation of the Cholesky Decomposition of real-valued matrices of
//    ac_fixed variables.
//
//    The Cholesky-Crout algorithm is used for Cholesky Decomposition.
//    By default, all temporary variables use the same precision as the
//    output. The user can change that by specifying number of bits to be
//    added or taken away from this default value. The user can also
//    add extra template parameters to turn off rounding and saturation
//    and rounding for temp. variables (these are turned on by default).
//
//    The user can also choose between using PWL vs. using accurate
//    division/sqrt functions, as mentioned earlier.
//
// Usage:
//    A sample implementation and its testbench looks like this:
//
//    #include <ac_fixed.h>
//    #include <ac_math/ac_chol_d.h>
//    #include <ac_matrix.h>
//
//    // Define data types for input and output matrices
//    typedef ac_matrix<ac_fixed<41, 21, true, AC_RND, AC_SAT>, 7, 7> input_type;
//    typedef ac_matrix<ac_fixed<64, 32, true, AC_RND, AC_SAT>, 7, 7> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_chol_d(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input;
//      output_type output;
//      // Hypothetical function that generates a pos. def. matrix:
//      gen_pos_def_matrix(input);
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    This version uses the C++ 2D array implementation for it's functioning.
//    It does this by using indirect_chol_d friend function from the
//    ac_matrix class to pass the 2D array member data to the relevant
//    implementation.
//
// -------------------------------------------------------------------------------

  template<bool use_pwl = false,
           int delta_w = 0, int delta_i = 0, ac_q_mode imd_Q = AC_RND, ac_o_mode imd_O = AC_SAT,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO,
           unsigned M>
  void ac_chol_d(
    const ac_matrix<ac_fixed<W, I, S, Q, O>, M, M> &A,
    ac_matrix<ac_fixed<outW, outI, outS, outQ, outO>, M, M> &L
  )
  {
    indirect_chol_d<use_pwl, delta_w, delta_i, imd_Q, imd_O>(A, L);
  }

// ===============================================================================
// Function: ac_chol_d (for ac_float using ac_matrix 2-D storage class)
//
// Description:
//    Calculation of the Cholesky Decomposition of real-valued matrices of
//    ac_float variables.
//
//    The Cholesky-Crout algorithm is used for Cholesky Decomposition.
//    By default, all temporary variables use the same precision as the
//    output. The user can change that by specifying number of bits to be
//    added or taken away from this default value. The user can also
//    add extra template parameters to turn off rounding for temp. variables
//    (turned on by default).
//
//    The user can also choose between using PWL vs. using accurate
//    division/sqrt functions, as mentioned earlier.
//
// Usage:
//    A sample implementation and its testbench looks like this:
//
//    #include <ac_float.h>
//    #include <ac_math/ac_chol_d.h>
//    #include <ac_matrix.h>
//
//    // Define data types for input and output matrices
//    typedef ac_matrix<ac_float<18, 2, 10, AC_RND>, 7, 7> input_type;
//    typedef ac_matrix<ac_float<18, 2, 10, AC_RND>, 7, 7> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_chol_d(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input;
//      output_type output;
//      // Hypothetical function that generates a pos. def. matrix:
//      gen_pos_def_matrix(input);
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    This version uses the C++ 2D array implementation for it's functioning.
//    It does this by using indirect_chol_d_float friend function from the
//    ac_matrix class to pass the 2D array member data to the relevant
//    implementation.
//
// -------------------------------------------------------------------------------

  template<bool use_pwl = false,
           int delta_w = 0, int delta_i = 0, int delta_e = 0, ac_q_mode imd_Q = AC_RND,
           int W, int I, int E, ac_q_mode Q,
           int outW, int outI, int outE, ac_q_mode outQ,
           unsigned M>
  void ac_chol_d(
    const ac_matrix<ac_float<W, I, E, Q>, M, M> &A,
    ac_matrix<ac_float<outW, outI, outE, outQ>, M, M> &L
  )
  {
    indirect_chol_d_float<use_pwl, delta_w, delta_i, imd_Q>(A, L);
  }

// ===============================================================================
// Function: ac_chol_d (for ac_std_float using ac_matrix 2-D storage class)
//
// Description:
//    Calculation of the Cholesky Decomposition of real-valued matrices of
//    ac_std_float variables.
//
//    The Cholesky-Crout algorithm is used for Cholesky Decomposition.
//    By default, all temporary variables use the same precision as the
//    output. The user can change that by specifying number of bits to be
//    added or taken away from this default value. The user can also
//    add extra template parameters to turn off rounding for temp. variables
//    (turned on by default).
//
//    The user can also choose between using PWL vs. using accurate
//    division/sqrt functions, as mentioned earlier.
//
//    This version depends on the ac_float implementation for computation.
//
// Usage:
//    A sample implementation and its testbench looks like this:
//
//    #include <ac_std_float.h>
//    #include <ac_math/ac_chol_d.h>
//    #include <ac_matrix.h>
//
//    // Define data types for input and output matrices
//    typedef ac_matrix<ac_std_float<32, 8>, 7, 7> input_type;
//    typedef ac_matrix<ac_std_float<32, 8>, 7, 7> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_chol_d(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input;
//      output_type output;
//      // Hypothetical function that generates a pos. def. matrix:
//      gen_pos_def_matrix(input);
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    This version uses the C++ 2D array implementation for it's functioning.
//    It does this by using indirect_chol_d_float friend function from the
//    ac_matrix class to pass the 2D array member data to the relevant
//    implementation.
//
// -------------------------------------------------------------------------------

  template<bool use_pwl = false,
           int delta_w = 0, int delta_i = 0, int delta_e = 0, ac_q_mode imd_Q = AC_RND,
           int W, int E,
           int outW, int outE,
           unsigned M>
  void ac_chol_d(
    const ac_matrix<ac_std_float<W, E>, M, M> &A,
    ac_matrix<ac_std_float<outW, outE>, M, M> &L
  )
  {
    indirect_chol_d_float<use_pwl, delta_w, delta_i, imd_Q>(A, L);
  }

// ===============================================================================
// Function: ac_chol_d (for ac_ieee_float using ac_matrix 2-D storage class)
//
// Description:
//    Calculation of the Cholesky Decomposition of real-valued matrices of
//    ac_ieee_float variables.
//
//    The Cholesky-Crout algorithm is used for Cholesky Decomposition.
//    By default, all temporary variables use the same precision as the
//    output. The user can change that by specifying number of bits to be
//    added or taken away from this default value. The user can also
//    add extra template parameters to turn off rounding for temp. variables
//    (turned on by default).
//
//    The user can also choose between using PWL vs. using accurate
//    division/sqrt functions, as mentioned earlier.
//
//    This version depends on the ac_float implementation for computation.
//
// Usage:
//    A sample implementation and its testbench looks like this:
//
//    #include <ac_ieee_float.h>
//    #include <ac_math/ac_chol_d.h>
//    #include <ac_matrix.h>
//
//    // Define data types for input and output matrices
//    typedef ac_matrix<ac_ieee_float<binary32>, 7, 7> input_type;
//    typedef ac_matrix<ac_ieee_float<binary32>, 7, 7> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_chol_d(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input;
//      output_type output;
//      // Hypothetical function that generates a pos. def. matrix:
//      gen_pos_def_matrix(input);
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    This version uses the C++ 2D array implementation for it's functioning.
//    It does this by using indirect_chol_d_float friend function from the
//    ac_matrix class to pass the 2D array member data to the relevant
//    implementation.
//
// -------------------------------------------------------------------------------

  template<bool use_pwl = false,
           int delta_w = 0, int delta_i = 0, int delta_e = 0, ac_q_mode imd_Q = AC_RND,
           ac_ieee_float_format Format,
           ac_ieee_float_format outFormat,
           unsigned M>
  void ac_chol_d(
    const ac_matrix<ac_ieee_float<Format>, M, M> &A,
    ac_matrix<ac_ieee_float<outFormat>, M, M> &L
  )
  {
    indirect_chol_d_float<use_pwl, delta_w, delta_i, imd_Q>(A, L);
  }

// ==============================================================================================
// Function: ac_chol_d (for ac_complex<ac_fixed> using ac_matrix 2-D storage class)
//
// Description:
//    Calculation of the Cholesky Decomposition of complex-valued matrices
//    with ac_complex<ac_fixed> variables.
//
//    This function can also be configured to change the type for the
//    temporary variables, and to use the pwl/accurate math functions for
//    intermediate computations, just like the ac_fixed version.
//
// Usage:
//    A sample implementation and its testbench looks like this:
//
//    #include <ac_fixed.h>
//    #include <ac_math/ac_chol_d.h>
//    using namespace ac_math;
//    #include <ac_matrix.h>
//
//    // Define data types for input and output matrices
//    typedef ac_matrix<ac_complex<ac_fixed<41, 21, true, AC_RND, AC_SAT> >, 7, 7> input_type;
//    typedef ac_matrix<ac_complex<ac_fixed<64, 32, true, AC_RND, AC_SAT> >, 7, 7> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_chol_d(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input;
//      output_type output;
//      // Hypothetical function that generates a pos. def. matrix:
//      gen_pos_def_matrix(input);
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    This version uses the C++ 2D array implementation for its functioning.
//    It does this by using indirect_chol_d friend function from the
//    ac_matrix class to pass the 2D array member data to the relevant
//    implementation.
//
// ----------------------------------------------------------------------------------------------

  template<bool use_pwl = false,
           int delta_w = 0, int delta_i = 0, ac_q_mode imd_Q = AC_RND, ac_o_mode imd_O = AC_SAT,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO,
           unsigned M>
  void ac_chol_d(
    const ac_matrix<ac_complex<ac_fixed<W, I, S, Q, O> >, M, M> &A,
    ac_matrix<ac_complex<ac_fixed<outW, outI, outS, outQ, outO> >, M, M> &L
  )
  {
    indirect_chol_d<use_pwl, delta_w, delta_i, imd_Q, imd_O>(A, L);
  }

} // namespace ac_math

#endif
