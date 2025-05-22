/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 3.8                                                 *
 *                                                                        *
 *  Release Date    : Tue May 13 15:34:32 PDT 2025                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.8.1                                               *
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
// File: ac_softmax_pwl_new.h (for ac_fixed)
//
// Description:
//    Synthesizable softmax function for AC fixed point datatypes.
//    New design normalizes the numerator and denominator terms in the softmax equation
//    to minimize bit growth and improve QofR.
//
// Usage:
//    Calculation of softmax of an array of real inputs, passed as ac_fixed variables.
//
// Notes:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_softmax_pwl_new.h>
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
//      ac_math::ac_softmax_pwl_new(input,output);
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

#ifndef _INCLUDED_AC_SOFTMAX_PWL_NEW_H_
#define _INCLUDED_AC_SOFTMAX_PWL_NEW_H_

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
  template<unsigned K, // Must always be specified as inputs are passed as pointers and size info is unavailable by default.
           int exp2f_fw_ = 0, // Fractional width for exp2_frac variables. If this is non-positive, the default value will be used.
           int l2e_fw_ = 0, // Fractional width for log2e constant. If this is non-positive, the default value will be used.
           ac_q_mode pwl_Q = AC_TRN,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  void ac_softmax_pwl_new_ptr(
    const ac_fixed<W, I, S, Q, O> input[K],
    ac_fixed<outW, outI, false, outQ, outO> output[K]
  )
  {
    typedef ac_fixed<W, I, S, Q, O> in_type;

    const ac_fixed<33, 1, false> full_prec_log2e = 1.4426950407214462757110595703125;
    static_assert(l2e_fw_ <= full_prec_log2e.width - full_prec_log2e.i_width, "l2e_fw_ must not exceed the number of fractional bits in full_prec_log2e.");
    constexpr int l2e_fw = l2e_fw_ > 0 ? l2e_fw_ : 16;
    typedef ac_fixed<l2e_fw + 1, 1, false> log2e_type;
    const log2e_type log2e = full_prec_log2e;
    typedef typename ac::rt_2T<log2e_type, in_type>::mult b2_pow_type;
    typedef ac_int<b2_pow_type::i_width, b2_pow_type::sign> b2p_I_type;
    b2p_I_type b2p_I[K];
    b2p_I_type b2p_Im;
    constexpr int exp2f_fw = exp2f_fw_ > 0 ? exp2f_fw_ : AC_MAX(1, AC_MIN(10, W - I - 2)) + 10;
    ac_fixed<exp2f_fw + 1, 1, false> exp2_frac[K];

    SEPARATE_INTO_INT_AND_FRAC: for (unsigned i = 0; i < K; i++) {
      b2_pow_type b2_pow = input[i]*log2e;
      b2p_I[i] = b2_pow.to_ac_int();
      if (b2_pow_type::width > b2_pow_type::i_width) {
        // fwidth must always be positive to prevent compiler errors.
        constexpr int fwidth = AC_MAX(b2_pow_type::width - b2_pow_type::i_width, 1);
        ac_fixed<fwidth, 0, false> frac_part;
        frac_part.set_slc(0, b2_pow.template slc<fwidth>(0));
        ac_math::ac_pow2_pwl<pwl_Q>(frac_part, exp2_frac[i]);
      } else {
        exp2_frac[i] = 1.0; // No fractional part => exp2(frac_part) = exp2(0) = 1
      }
      if (i == 0) { b2p_Im = b2p_I[i]; }
      else { b2p_Im = AC_MAX(b2p_I[i], b2p_Im); }
    }

    ac_fixed<exp2f_fw + 1, 1, false> prod[K];
    constexpr int sprod_Iw = ac::nbits<2*K - 1>::val;
    ac_fixed<exp2f_fw + sprod_Iw, sprod_Iw, false> prod_sum = 0.0;

    CALCULATE_NORM_NUM_AND_DEN: for (unsigned i = 0; i < K; i++) {
      ac_int<b2_pow_type::i_width, false> del_b2p_I = b2p_Im - b2p_I[i];
      prod[i] = exp2_frac[i] >> del_b2p_I;
      prod_sum += prod[i];
    }

    constexpr int recip_W = prod_sum.width + 20;
    constexpr int recip_iW = prod_sum.width - prod_sum.i_width + 1;

    ac_fixed<recip_W, recip_iW, false> recip_sum;
    ac_reciprocal_pwl<pwl_Q>(prod_sum, recip_sum);

    CALCULATE_SOFTMAX_OUT: for (unsigned i = 0; i < K; i++) { output[i] = prod[i]*recip_sum; }
  }

  template<int exp2f_fw_ = 0, // Fractional width for exp2_frac variables. If this is non-positive, the default value will be used.
           int l2e_fw_ = 0, // Fractional width for log2e constant. If this is non-positive, the default value will be used.
           ac_q_mode pwl_Q = AC_TRN,
           unsigned K,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  // By declaring the function parameters as references to arrays (i.e. "(&input)[K]" and "(&output)[K]"), we ensure
  // that template parameter deduction infers the value of K.
  void ac_softmax_pwl_new(
    const ac_fixed<W, I, S, Q, O> (&input)[K],
    ac_fixed<outW, outI, false, outQ, outO> (&output)[K]
  )
  {
    ac_softmax_pwl_new_ptr<K, exp2f_fw_, l2e_fw_, pwl_Q>(input, output);
  }
};


#endif // #ifndef _INCLUDED_AC_SOFTMAX_PWL_NEW_H_
