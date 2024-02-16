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
//*****************************************************************************************
// File: ac_softmax_pwl.h (for ac_fixed)
//
// Description: Synthesizable softmax function for AC fixed point datatypes.
// Usage:
//    Calculation of softmax of an array of real inputs, passed as ac_fixed variables.
//
// Notes:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_softmax_pwl.h>
//    using namespace ac_math;
//
//    const int num_logits_tb = 20;
//
//    typedef ac_fixed<16, 8, true, AC_TRN, AC_SAT> input_type;
//    typedef ac_fixed<16, 8, false, AC_TRN, AC_SAT> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type (&input)[num_logits_tb],
//      output_type (&output)[num_logits_tb]
//    )
//    {
//      ac_softmax_pwl(input,output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input[num_logits_tb];
//      for (int i = 0; i < num_logits_tb; i++) { input[i] = 9 - i; }
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    Attempting to call this function with a type that is not implemented will
//    result in a compile-time error.
//
//    This file uses the ac_exp_pwl() function from the file ac_pow_pwl.h and the
//    ac_reciprocal_pwl() function from the file ac_reciprocal_pwl.h.
//
// Revision History:
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//    3.3.0  - [CAT-25798] Added CDesignChecker fixes/waivers for code check and Synthesis-simulation mismatch/violations in ac_math PWL and Linear Algebra IPs.
//    3.2.0 - Initial version
//
//*****************************************************************************************

#ifndef _INCLUDED_AC_SOFTMAX_PWL_H_
#define _INCLUDED_AC_SOFTMAX_PWL_H_

// The function uses default template parameters, which are only supported by C++11 or later
// compiler standards. Hence, the user should be informed if they are not using those standards.
#if (defined(__GNUC__) && (__cplusplus < 201103L))
#error Please use C++11 or a later standard for compilation.
#endif
#if (defined(_MSC_VER) && (_MSC_VER < 1920) && !defined(__EDG__))
#error Please use Microsoft VS 2019 or a later standard for compilation.
#endif

// Include header for supported datatype.
#include <ac_fixed.h>

// Include headers for required functions
#include <ac_math/ac_pow_pwl.h>
#include <ac_math/ac_reciprocal_pwl.h>

// Encapsulate within ac_math namespace.
namespace ac_math
{
  // or_e: Override the default type assigned to the exponent variable. Use user-supplied type information instead.
  // or_r: Same as above, but for the reciprocal variable.
  // K: number of input elements.
  template<ac_q_mode pwl_Q = AC_TRN,
           bool or_e = false, int iW_e = 0, int iI_e = 0, ac_q_mode iQ_e = AC_TRN, ac_o_mode iO_e = AC_WRAP,
           bool or_r = false, int iW_r = 0, int iI_r = 0, ac_q_mode iQ_r = AC_TRN, ac_o_mode iO_r = AC_WRAP,
           unsigned K,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  // By declaring the function parameters as references to arrays (i.e. "(&input)[K]" and "(&output)[K]"), we ensure
  // that template parameter deduction infers the value of K.
  void ac_softmax_pwl(
    const ac_fixed<W, I, S, Q, O> (&input)[K],
    ac_fixed<outW, outI, false, outQ, outO> (&output)[K]
  )
  {
    // The default number of fractional bits assigned to the exponent variables assumes that the PWL implementation
    // has four segments and uses 10 fractional bits to store the slope and intercept values. For any other implementation
    // the default number may have to be changed to ensure no loss of precision.
    const int exp_frac_bits = or_e ? iW_e - iI_e : AC_MAX(1, AC_MIN(10, W - I - 2)) + 10;
    const int exp_int_bits  = or_e ? iI_e : int(1.443*double(1 << (I - S))) + 1;
    // The default rounding and saturation modes are AC_TRN and AC_WRAP, respectively.
    const ac_q_mode exp_Q = or_e ? ac_q_mode(iQ_e) : ac_q_mode(AC_TRN);
    const ac_o_mode exp_O = or_e ? ac_o_mode(iO_e) : ac_o_mode(AC_WRAP);
    // The below static_assert limits the size the integer width of the exponent variable can reach.
    static_assert(exp_int_bits <= 64, "Intermediate bitwidth calculation gives a very large value for integer bits. Consider reducing the number of input integer bits.");
    typedef ac_fixed<exp_frac_bits + exp_int_bits, exp_int_bits, false, exp_Q, exp_O> T_exp;
    typedef ac_fixed<T_exp::width + ac::log2_ceil<K>::val, T_exp::i_width + ac::log2_ceil<K>::val, false> T_sum;

    #pragma hls_waive APT
    T_exp exp_arr[K];
    // The reciprocal variable is also assigned a default bitwidth based on the default PWL bitwidths and segments (10 fractional bits and 8 segments).
    // For any PWL implementation other than the default, these bitwidths may change.
    const int recip_W = or_r ? iW_r : T_sum::width + 20;
    const int recip_I = or_r ? iI_r : T_sum::width - T_sum::i_width + 1;
    const ac_q_mode recip_Q = or_r ? ac_q_mode(iQ_r) : ac_q_mode(AC_TRN);
    const ac_o_mode recip_O = or_r ? ac_o_mode(iO_r) : ac_o_mode(AC_WRAP);
    typedef ac_fixed<recip_W, recip_I, false, recip_Q, recip_O> T_recip;
    T_recip sum_exp_recip;

    // All the loops used in this function can be pipelined/unrolled to give the desired area/throughput score.
    // 1. Pipelining all the loops and the main function call with an II of 1 gives (2K) number of throughput cycles.
    // 2. Unrolling all the loops and pipelining the main function call with an II of 1 ensures a throughput of 1.
    // The second option can also result in a very large area score.

    // Calculate exponential of all inputs.
    CALC_EXP_LOOP: for (unsigned i = 0; i < K; i++) { ac_exp_pwl<pwl_Q>(input[i], exp_arr[i]); }

    // Perform a MAC operation to add all the exponential values.
    T_sum sum_exp = 0.0;
    SUM_EXP_LOOP: for (unsigned i = 0; i < K; i++) { sum_exp += exp_arr[i]; }

    // Find the reciprocal of the sum of exponentials.
    ac_reciprocal_pwl<pwl_Q>(sum_exp, sum_exp_recip);

    // The types for the exponent and reciprocal variable are configurable primarily to ensure that the size of the multiplier used for the multiplication below
    // does not become too large. A large multiplier can become an issue if the loop below is being unrolled and the number of iterations (K) is large,
    // hence resulting in many large multipliers and a very large area.
    CALC_SOFTMAX_LOOP: for (unsigned i = 0; i < K; i++) { output[i] = sum_exp_recip*exp_arr[i]; }
  }

};


#endif // #ifndef _INCLUDED_AC_SOFTMAX_PWL_H_















