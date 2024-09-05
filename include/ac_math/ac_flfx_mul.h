/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 3.6                                                 *
 *                                                                        *
 *  Release Date    : Sun Aug 25 18:24:45 PDT 2024                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.6.0                                               *
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
// File: ac_flfx_mul.h
//
// Description: Optimized functions for the following mixed ac_float, ac_fixed multiplication.
//    i.   ac_float*ac_fixed=ac_float
//    ii.  ac_float*ac_fixed=ac_fixed
//    iii. ac_float*ac_float=ac_fixed
//    The functions accept eithe an ac_fixed and an ac_float input or two ac_float inputs 
//    and produce an ac_float output for case (i) and ac_fixed output for cases (ii and iii).
//
// Usage:
//    A sample testbench and its implementation look like this:
//
//    #include <ac_math/ac_flfx_mul.h>
//
//    typedef ac_float<14, 2, 6> fl_in_type;
//    typedef ac_float<14, 2, 6> fl_out_type;
//    typedef ac_fixed<10, 2, true> fx_in_type;
//    
//    typedef fl_in_type inA_type;
//    typedef fx_in_type inB_type;
//    typedef fl_out_type out_type;
//
//    #pragma hls_design top
//    void project(
//      const inA_type &fl_in,
//      const inB_type &fx_in,
//      out_type &out
//    ) {
//      ac_math::ac_flfx_mul(fl_in, fx_in, out);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    
//    CCS_MAIN(int, char**) {
//      inA_type fl_in = 1.5;
//      inB_type fx_in = 0.5;
//      out_type out;
//      
//      CCS_DESIGN(project)(fl_in, fx_in, out);
//
//      CCS_RETURN(0);
//    }
//    #endif
//
// Notes:
//    This file uses C++ function overloading where the order of ac_fixed and ac_float
//    inputs in the parameter lists of the functions is interchanged. Hence, the user
//    does not have to worry about the order in which they supply the function inputs while
//    calling it, when the two inputs are ac_float and ac_fixed.
//
//    Attempting to call the function with a type that is not implemented will result in a
//    compile error.
//
//    Activating subnormal support can be done using the subn_support template parameter.
//    More details on this parameter are given in the AC Math reference manual.
//
//    Output truncation mode must always be AC_TRN for ac_float outputs.
//
// Revision History:
//    3.6.0  - [CAT-36087] Added library to ac_math subproject.
//    3.7.1  - Added support for fixed*float=fixed and float*float=fixed
//
//*****************************************************************************************

#ifndef _INCLUDED_AC_FLFX_MUL_H_
#define _INCLUDED_AC_FLFX_MUL_H_

#include <ac_fixed.h>
#include <ac_float.h>

namespace ac_math
{
  // Subnormals supported.
  template <class in1_type, class in2_type, class out_type, bool subn_support>
  struct ac_flfx_mul_helper_st {
    static const ac_q_mode Q1 = in1_type::q_mode;
    static const ac_q_mode outQ = out_type::q_mode;

    enum {
      W1 = in1_type::width,
      I1 = in1_type::i_width,
    };

    typedef typename ac::rt_2T<ac_fixed<W1, I1, true>, in2_type>::mult mprod_type;

    enum {
      mprod_W = mprod_type::width,
      mprod_I = mprod_type::i_width,
      outW = out_type::width,
      outI = out_type::i_width,
      rshift = mprod_I - outI,
      max_oexp_t = in1_type::MAX_EXP + rshift,
      min_oexp_t = in1_type::MIN_EXP - (mprod_W - 1) + rshift,
      min_oexp_t_abs = (min_oexp_t < 0) ? -min_oexp_t : min_oexp_t,
      temp_oexp_w = ac::nbits<AC_MAX(max_oexp_t, min_oexp_t_abs)>::val + 1,
      // Max temp. exponent value possible if output is zero. This is assuming that a zero ac_float input
      // has the most extreme exponent value, i.e. MAX_EXP.
      max_oexp_t_oz = in1_type::MAX_EXP - (mprod_W - 1) + rshift,
      // Is extra masking on oexp_temp to set it to zero for zero inputs needed?
      extra_exp_masking = max_oexp_t_oz > out_type::MAX_EXP,
    };
    
    typedef ac_int<temp_oexp_w, true> oet_type; // Type of temporary output exponent
    // Type for in2_ls output. Since we don't normalize in2, in2_ls is always 0 and we only need one
    // bit to represent that.
    typedef ac_int<1, false> in2_ls_type;
    // Type of mprod_ls output. Since we normalize all of mprod_ls, we use the leading_sign typedef
    // internal to mprod_type.
    typedef typename mprod_type::rt_unary::leading_sign mprod_ls_type;
    // type for sum of ls values.
    typedef mprod_ls_type ls_sum_type;

    in2_ls_type calc_in2_ls(const in2_type &in2) {
      return 0;
    }

    mprod_ls_type calc_mprod_ls(const mprod_type &mprod) {
      // We're normalize all of mprod => leading_sign() will also traverse all mprod bits.
      return mprod.leading_sign();
    }

    // Perform some final calculations on output mantissa to account for exponents that are greater or
    // less than out_type::MAX_EXP or out_type::MIN_EXP, respectively.
    void final_omant_calc(const oet_type &oexp_temp, const in1_type &in1, ac_fixed<outW, outI, true, outQ> &out_mant) {
      constexpr int max_oexp = out_type::MAX_EXP;
      constexpr int min_oexp = out_type::MIN_EXP;

      // Saturate as close to +inf or -inf as possible if the exponent exceeds the maximum. This is done
      // via inexpensive bit-masking logic.
      bool exp_gt_max = (oexp_temp > max_oexp);
      ac_int<outW - 1, false> out_mant_lsbs = out_mant.template slc<outW - 1>(0);
      bool amask_bit_neg = (out_mant[outW - 1] && exp_gt_max);
      ac_int<outW - 1, false> out_mant_lsbs_amask = -int(!amask_bit_neg);
      out_mant_lsbs &= out_mant_lsbs_amask;
      ac_int<outW - 1, false> out_mant_lsbs_omask = -int(!out_mant[outW - 1] && exp_gt_max);
      out_mant_lsbs |= out_mant_lsbs_omask;
      out_mant.set_slc(0, out_mant_lsbs);
      
      // Take care of subnormal outputs, where the temporary output exponent is less than the minimum.
      int exp_dif = oexp_temp.to_int() - min_oexp;

      if (exp_dif < 0) {
        out_mant <<= exp_dif;
      }
    }
  };

  // Subnormals not supported.
  template <class in1_type, class in2_type, class out_type>
  struct ac_flfx_mul_helper_st<in1_type, in2_type, out_type, false> {
    static const ac_q_mode Q1 = in1_type::q_mode;
    static const ac_q_mode outQ = out_type::q_mode;

    enum {
      W1 = in1_type::width,
      I1 = in1_type::i_width,
    };

    typedef typename ac::rt_2T<ac_fixed<W1, I1, true>, in2_type>::mult mprod_type;

    enum {
      mprod_W = mprod_type::width,
      mprod_I = mprod_type::i_width,
      W2 = in2_type::width,
      S = in2_type::sign,
      outW = out_type::width,
      outI = out_type::i_width,
      rshift = mprod_I - outI,
      max_oexp_t = in1_type::MAX_EXP + rshift,
      min_oexp_t = in1_type::MIN_EXP - (W2 - 1) - (1 + int(S)) + rshift,
      min_oexp_t_abs = (min_oexp_t < 0) ? -min_oexp_t : min_oexp_t,
      temp_oexp_w = ac::nbits<AC_MAX(max_oexp_t, min_oexp_t_abs)>::val + 1,
      // Max temp. exponent value possible if output is zero. This is assuming that a zero ac_float input
      // has the most extreme exponent value, i.e. MAX_EXP.
      max_oexp_t_oz = in1_type::MAX_EXP - (1 + int(S)) + rshift,
      // Is extra masking on oexp_temp to set it to zero for zero inputs needed?
      extra_exp_masking = max_oexp_t_oz > out_type::MAX_EXP,
    };
    
    typedef ac_int<temp_oexp_w, true> oet_type; // Type of temporary output exponent
    // Type for in2_ls output. Since we normalize all of in2, we use the leading_sign typedef internal
    // to in2_type.
    typedef typename in2_type::rt_unary::leading_sign in2_ls_type;
    // The maximum leading_sign output for unsigned fixed point inputs is 1, and that for signed fixed
    // point inputs is 2. The number of bits needed to store either maximum output is given by 1 + int(S).
    typedef ac_int<1 + int(S), false> mprod_ls_type;
    // type for sum of ls values.
    typedef ac_int<ac::nbits<(W2 - 1) + 1 + int(S)>::val, false> ls_sum_type;

    in2_ls_type calc_in2_ls(const in2_type &in2) {
      return in2.leading_sign();
    }

    mprod_ls_type calc_mprod_ls(const mprod_type &mprod) {
      // Since both the fixed and floating point inputs to the multipliers are considered normal
      // (subnormal values are treated as zero later, in final_omant_calc, and the final ac_float
      // constructor takes care of zero values), mprod is mostly normalized. The leading_sign calculation
      // will only need to traverse 3 MSBs of mprod, if the fixed point input is signed, or 2 MSBs, if the
      // fixed point input is unsigned. This number of MSBs is given by the 2 + int(S) calculation below.
      constexpr int num_msbs = 2 + int(S);
      ac_int<num_msbs, true> mprod_msbs = mprod.template slc<num_msbs>(mprod_W - num_msbs);
      return mprod_msbs.leading_sign();
    }

    // Perform some final calculations on output mantissa to account for exponents that are greater or
    // less than out_type::MAX_EXP or out_type::MIN_EXP, respectively.
    void final_omant_calc(const oet_type &oexp_temp, const in1_type &in1, ac_fixed<outW, outI, true, outQ> &out_mant) {
      constexpr int max_oexp = out_type::MAX_EXP;
      constexpr int min_oexp = out_type::MIN_EXP;

      bool exp_gt_max = (oexp_temp > max_oexp);
      bool exp_lt_min = (oexp_temp < min_oexp);

      // The handling below does the following to the output mantissa:
      // 1. Saturate as close to +inf or -inf as possible if the exponent exceeds the maximum.
      // 2. Set to zero if the exponent is less than the minimum. In other words, outputs that are too
      //    small to be normalized are set to 0. Subnormals will not be represented.
      // 3. Set to zero if the floating point input is subnormal.
      //
      // All of this is done via inexpensive bit-masking.
      ac_int<outW - 1, false> out_mant_lsbs = out_mant.template slc<outW - 1>(0);
      bool amask_bit_neg = ((out_mant[outW - 1] && exp_gt_max) || exp_lt_min || !in1.is_norm());
      ac_int<outW - 1, false> out_mant_lsbs_amask = -int(!amask_bit_neg);
      out_mant_lsbs &= out_mant_lsbs_amask;
      ac_int<outW - 1, false> out_mant_lsbs_omask = -int(!out_mant[outW - 1] && exp_gt_max);
      out_mant_lsbs |= out_mant_lsbs_omask;
      out_mant[outW - 1] = out_mant[outW - 1] && !exp_lt_min && in1.is_norm();

      out_mant.set_slc(0, out_mant_lsbs);
    }
  };
  
  //==========================================================================
  // Function: ac_flfx_mul (ac_float and ac_fixed inputs, ac_float output)
  //
  // Description:
  //    Optimized ac_float*ac_fixed multiplication where the ac_float input
  //    comes first in the parameter list and the ac_fixed input comes second.
  //
  // Usage:
  //    See above example code for usage.
  //
  //--------------------------------------------------------------------------

  // Declaring subn_support as a template parameter allows us to use class templates to differentiate
  // between designs that have subnormal support and those that don't. This in turn means better coverage,
  // as we don't have to use if-else branches for different versions where the condition is a constant.
  template <
    bool subn_support = false,
    int W1, int I1, int E1, ac_q_mode Q1,
    int W2, int I2, bool S, ac_q_mode Q2, ac_o_mode O,
    int outW, int outI, int outE, ac_q_mode outQ
    >
  void ac_flfx_mul(
    const ac_float<W1, I1, E1, Q1> &in1,
    const ac_fixed<W2, I2, S, Q2, O> &in2,
    ac_float<outW, outI, outE, outQ> &out
  )
  {
    static_assert(outQ == AC_TRN, "Output quantization mode can only be AC_TRN.");

    AC_ASSERT(in1.is_norm() || in1.is_subn() || in1.mantissa() == 0, "Floating point inputs must either be (a) normalized (b) subnormal or (c) zero. Make sure you at least attempted to normalize in1.");

    typedef ac_float<W1, I1, E1, Q1> in1_type;
    typedef ac_fixed<W2, I2, S, Q2, O> in2_type;
    typedef ac_float<outW, outI, outE, outQ> out_type;
    typedef ac_flfx_mul_helper_st<in1_type, in2_type, out_type, subn_support> hst_type;
    hst_type hst_obj;

    // If subn_support = false, we normalize in2 by finding the position of the leading sign and
    // left-shifting accordingly.
    // If subn_support = true, in2 is not normalized and in2_ls is always equal to 0.
    typename hst_type::in2_ls_type in2_ls = hst_obj.calc_in2_ls(in2);
    in2_type in2_mul = (in2 << in2_ls);

    // Multiply in1 mantissa with in2_mul. The result is called the "mantissa product" and is stored in
    // the variable "mprod".
    typedef typename hst_type::mprod_type mprod_type;
    constexpr int mprod_W = mprod_type::width;
    constexpr int mprod_I = mprod_type::i_width;
    mprod_type mprod = in1.mantissa()*in2_mul;

    // Normalize mantissa product.
    // If subn_support = false, we need to traverse fewer mprod bits for the leading_sign operation. Refer
    //   to the code in the calc_mprod_ls function for more details. The shifter size will also be smaller
    //   because the leading_sign output is smaller.
    // If subn_support = true, we need to traverse all of mprod for the leading_sign operation.
    typename hst_type::mprod_ls_type mprod_ls = hst_obj.calc_mprod_ls(mprod);
    mprod_type mprod_norm = (mprod << mprod_ls);

    // To prevent extra assign_from logic, the final mantissa and exponent arguments to the ac_float
    // constructor call below must match the precision of the output mantissa and exponent.
    // The first step towards this is slicing all the mantissa product bits into a variable with outI
    // integer bits, i.e. out_mant_all.
    ac_fixed<mprod_W, outI, true> out_mant_all;
    out_mant_all.set_slc(0, mprod_norm.template slc<mprod_W>(0));
    // Finally, assign out_mant_all to a variable with outW bits. Since outQ is only allowed to be AC_TRN
    // for now, all extra fractional bits are truncated.
    ac_fixed<outW, outI, true, outQ> out_mant = out_mant_all;
    // The slicing operation is equivalent to right-shifting by (mprod_I - outI).
    constexpr int rshift = mprod_I - outI;
    
    typename hst_type::ls_sum_type ls_sum = in2_ls + mprod_ls;
    // Account for all the normalization and shifting we've done in the calculation of the temporary
    // output exponent value (oexp_temp).
    typename hst_type::oet_type oexp_temp = (in1.exp() - ls_sum).to_int() + rshift;
    
    if (hst_type::extra_exp_masking) {
      ac_int<oexp_temp.width, false> oexp_t_bits = oexp_temp;
      ac_int<oexp_temp.width, false> oexp_t_mask = -int(!!out_mant);
      oexp_t_bits &= oexp_t_mask;
      oexp_temp = oexp_t_bits;
    }
    
    // Bring the temporary output exponent into the correct range. The only reason an intermediate fixed
    // point variable is used is so that we can turn on saturation.
    ac_fixed<outE, outE, true, AC_TRN, AC_SAT> out_exp_fx = oexp_temp;
    ac_int<outE, true> out_exp = out_exp_fx.to_int();
    
    // Account for large positive or negative exponents.
    hst_obj.final_omant_calc(oexp_temp, in1, out_mant);
    
    // Calculate temporary output with ac_float constructor call.
    ac_float<outW, outI, outE, outQ> out_temp(out_mant, out_exp, false);
    out = out_temp;
  }
  
  //==========================================================================
  // Function: ac_flfx_mul (ac_float and ac_fixed inputs, ac_float output)
  //
  // Description:
  //    Optimized ac_float*ac_fixed multiplication where the ac_fixed input
  //    comes first in the parameter list and the ac_float input comes second.
  //
  // Usage:
  //    #include <ac_math/ac_flfx_mul.h>
  //    
  //    typedef ac_fixed<10, 2, true> fx_in_type;
  //    typedef ac_float<14, 2, 6> fl_in_type, out_type;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const fx_in_type &fx_in,
  //      const fl_in_type &fl_in,
  //      out_type &out
  //    )
  //    {
  //      ac_math::ac_flfx_mul(fx_in, fl_in, out);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //    
  //    CCS_MAIN(int, char **)
  //    {
  //      fx_in_type fx_in = 0.5;
  //      fl_in_type fl_in = 1.5;
  //      out_type out;
  //    
  //      CCS_DESIGN(project)(fx_in, fl_in, out);
  //    
  //      CCS_RETURN(0);
  //    }
  //    #endif
  //    
  // Notes:
  //    This function uses the ac_flfx_mul() version above for the actual
  //    computation.
  //
  //--------------------------------------------------------------------------

  template <
    bool subn_support = false,
    int W1, int I1, bool S, ac_q_mode Q1, ac_o_mode O,
    int W2, int I2, int E2, ac_q_mode Q2,
    int outW, int outI, int outE, ac_q_mode outQ
    >
  void ac_flfx_mul(
    const ac_fixed<W1, I1, S, Q1, O> &in1,
    const ac_float<W2, I2, E2, Q2> &in2,
    ac_float<outW, outI, outE, outQ> &out
  )
  {
    ac_flfx_mul<subn_support>(in2, in1, out);
  }


  //==========================================================================
  // Function: ac_flfx_mul (ac_float and ac_fixed inputs, ac_fixed output)
  //
  // Description:
  //    Optimized ac_float*ac_fixed multiplication where the ac_float input
  //    comes first in the parameter list and the ac_fixed input comes second.
  //    The result is an ac_fixed.
  //
  // Usage:
  //    See above example code for usage.
  //
  //--------------------------------------------------------------------------

  // Declaring subn_support as a template parameter allows us to use class templates to differentiate
  // between designs that have subnormal support and those that don't. This in turn means better coverage,
  // as we don't have to use if-else branches for different versions where the condition is a constant.
  // When subn_support is true, the operation has increased area because it tries to identify subnormal inputs
  // and use a zero value instead.
  template <
    bool subn_support = false,
    int W1, int I1, int E1, ac_q_mode Q1,
    int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2,
    int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO
    >
  void ac_flfx_mul(
    const ac_float<W1, I1, E1, Q1> &in1,
    const ac_fixed<W2, I2, S2, Q2, O2> &in2,
    ac_fixed<outW, outI, outS, outQ, outO> &out
  )
  { 
    AC_ASSERT(in1.is_norm() || in1.is_subn() || in1.mantissa() == 0, "Floating point inputs must either be (a) normalized (b) subnormal or (c) zero. Make sure you at least attempted to normalize in1.");

    typedef ac_float<W1, I1, E1, Q1> in1_type;
    typedef ac_fixed<W2, I2, S2, Q2, O2> in2_type;
    typedef ac_fixed<outW, outI, outS, outQ, outO> out_type;

    typedef typename ac::rt_2T<ac_fixed<W1, I1, true>, ac_fixed<W2, I2, S2>>::mult mprod_type;
    typedef mprod_type mprod_t;

    // Check if the float input is subnormal. When subn_support=false the subnormal input will result in zero output.
    // Check only if the input is not normalized. This is enough because it means that the input is either zero or subnormal.
    // In either case the result will be zero, if the subn_support is false.
    bool in1_subn= !in1.is_norm();
    bool r_zero = !subn_support & in1_subn;

    // Multiply in1 mantissa with in2. The result is called the "mantissa product" and is stored in
    // the variable "mprod".
    mprod_t mprod = in1.mantissa()*in2;

    // Shift the product according to the exponent to create the ac_fixed value.
    out_type fxpt_temp;
    
    typedef typename ac_float<mprod_t::width, mprod_t::i_width, E1>::rt_unary::to_ac_fixed_t full_prec_t;
    full_prec_t to_shift = full_prec_t(mprod);
    to_shift <<= in1.exp();
    fxpt_temp = out_type(to_shift);

    // Set the output to the result of the product or to zero, depending on the subnormal check.
    out = (r_zero) ? out_type(0) : out_type(fxpt_temp);
  }

  //==========================================================================
  // Function: ac_flfx_mul (ac_float and ac_fixed inputs, ac_fixed output)
  //
  // Description:
  //    Optimized ac_float*ac_fixed multiplication where the ac_fixed input
  //    comes first in the parameter list and the ac_float input comes second.
  //    The result is retuned as an ac_fixed.
  //
  // Usage:
  //    #include <ac_math/ac_flfx_mul.h>
  //    
  //    typedef ac_fixed<7, 3, true> fx_in_type;
  //    typedef ac_float<9, 2, 3> fl_in_type;
  //    typedef ac_fixed<25, 7, true> out_type;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const fl_in_type &fi_in,
  //      const fx_in_type &fx_in,
  //      out_type &out
  //    )
  //    {
  //      ac_math::ac_flfx_mul(fx_in, fl_in, out);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //    
  //    CCS_MAIN(int, char **)
  //    {
  //      fx_in_type fx_in = 0.5;
  //      fl_in_type fl_in = 1.5;
  //      out_type out;
  //    
  //      CCS_DESIGN(project)(fx_in, fl_in, out);
  //    
  //      CCS_RETURN(0);
  //    }
  //    #endif
  //    
  // Notes:
  //    This function uses the ac_flfx_mul() function defined above to compute
  //    the result.
  //
  //--------------------------------------------------------------------------
  template <
    bool subn_support = false,
    int W1, int I1, bool S1, ac_q_mode Q1, ac_o_mode O1,
    int W2, int I2, int E2, ac_q_mode Q2,
    int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO
    >
  void ac_flfx_mul(
    const ac_fixed<W1, I1, S1, Q1, O1> &in1,
    const ac_float<W2, I2, E2, Q2> &in2,
    ac_fixed<outW, outI, outS, outQ, outO> &out
  ) {
    ac_flfx_mul<subn_support>(in2, in1, out);
  }
  //==========================================================================
  // Function: ac_flfx_mul (two ac_float inputs, ac_fixed output)
  //
  // Description:
  //    Optimized ac_float*ac_float multiplication where the result is ac_fixed.
  //
  // Usage:
  //    #include <ac_math/ac_flfx_mul.h>
  //
  //    typedef ac_float<5, 2, 3> fl_in_type;
  //    typedef ac_fixed<18, 8, true> out_type;
  //
  //    #pragma hls_design top
  //    void project(
  //      const fl_in_type &fl1_in,
  //      const fl_in_type &fl2_in,
  //      out_type &out
  //    ) {
  //      ac_math::ac_flfx_mul(fl1_in, fl2_in, out);
  //    }
  //
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //    
  //    CCS_MAIN(int, char**) {
  //      fl_in_type fl1_in = 1.5;
  //      fl_in_type fl2_in = 0.5;
  //      out_type out;
  //      
  //      CCS_DESIGN(project)(fl1_in, fl2_in, out);
  //
  //      CCS_RETURN(0);
  //    }
  //    #endif
  //
  //--------------------------------------------------------------------------
  template <
    bool subn_support = false,
    int W1, int I1, int E1, ac_q_mode Q1,
    int W2, int I2, int E2, ac_q_mode Q2, 
    int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO
    >
  void ac_flfx_mul(
    const ac_float<W1, I1, E1, Q1> &in1,
    const ac_float<W2, I2, E2, Q2> &in2,
    ac_fixed<outW, outI, outS, outQ, outO> &out
  )
  {
    AC_ASSERT(in1.is_norm() || in1.is_subn() || in1.mantissa() == 0, "Floating point inputs must either be (a) normalized (b) subnormal or (c) zero. Make sure you at least attempted to normalize in1.");
    AC_ASSERT(in2.is_norm() || in2.is_subn() || in2.mantissa() == 0, "Floating point inputs must either be (a) normalized (b) subnormal or (c) zero. Make sure you at least attempted to normalize in2.");
    
    typedef ac_float<W1, I1, E1, Q1> in1_type;
    typedef ac_float<W2, I2, E2, Q2> in2_type;
    typedef ac_fixed<outW, outI, outS, outQ, outO> out_type;

    typedef typename ac::rt_2T<ac_fixed<W1, I1, true>, ac_fixed<W2, I2, true>>::mult mprod_type;
    typedef typename ac::rt_2T<ac_int<E1, true>, ac_int<E2, true>>::plus esum_type;
    typedef mprod_type mprod_t;
    typedef esum_type esum_t;
    
    // Check if either input is subnormal. When subn_support=false a subnormal input will result in zero output.
    // Check only if the input is not normalized. This is enough because it means that the input is either zero or subnormal.
    // In either case the result will be zero, if the subn_support is false.
    bool in1_subn= !in1.is_norm();
    bool in2_subn= !in2.is_norm();
    bool r_zero = !subn_support & (in1_subn | in2_subn);

    // Add in1 exponent and in2 exponent. The result is called the "exponent sum" and is stored in the 
    // variable "esum"
    esum_t esum = in1.exp() + in2.exp();

    // Multiply in1 mantissa with in2 mantissa. The result is called the "mantissa product" and is stored in
    // the variable "mprod".
    mprod_t mprod = in1.mantissa()*in2.mantissa();

    // Shift the product according to the exponent to create the ac_fixed value.
    out_type fxpt_temp;

    typedef typename ac_float<mprod_t::width, mprod_t::i_width, esum_t::width>::rt_unary::to_ac_fixed_t full_prec_t;
    full_prec_t to_shift = full_prec_t(mprod);
    to_shift <<= esum;
    fxpt_temp = out_type(to_shift);

    // Set the output to the result of the product or to zero, depending on the subnormal check.
    out = (r_zero) ? out_type(0) : out_type(fxpt_temp);
  }
}

#endif // #ifndef _INCLUDED_AC_FLFX_MUL_H_
