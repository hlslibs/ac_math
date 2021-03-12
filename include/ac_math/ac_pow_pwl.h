/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 3.4                                                 *
 *                                                                        *
 *  Release Date    : Sat Jan 23 14:58:27 PST 2021                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.4.0                                               *
 *                                                                        *
 *  Copyright , Mentor Graphics Corporation,                     *
 *                                                                        *
 *  All Rights Reserved.                                                  *
 *  
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
//********************************************************************************************
// File: ac_pow_pwl.h
//
// Description:
//    Provides piece-wise linear implementations of the
//    base 2 and base e exponential functions for ac_fixed inputs.
//
// Usage:
//    A sample testbench and its implementation look like this:
//
//    #include <ac_math/ac_pow_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_fixed<20, 11, true, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<24, 14, false, AC_RND, AC_SAT> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_pow2_pwl(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 2.5;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN(0);
//    }
//    #endif
//
// Notes:
//    Attempting to call the function with a type that is not implemented will result
//    in a compile error.
//    The file uses the ac_log2_pwl function from the ac_log_pwl.h file and the
//    ac_shift_left function from the ac_shift.h header file.
//
// Revision History:
//    3.3.0  - [CAT-25797] Added CDesignChecker fixes/waivers for code check violations in ac_math PWL and Linear Algebra IPs.
//             Waivers added for CNS and CCC violations.
//             Fixes added for FXD, STF and MXS violations.
//               - FXD violations fixed by changing integer literals to floating point literals or typecasting to ac_fixed values.
//               - STF violations fixed by using "const" instead of "static const" parameters. LUT generator files also print out "const" LUTs instead of "static const" LUTs.
//               - MXS violations fixed by typecasting unsigned variables to int.
//    Niramay Sanghvi : Aug 15 2017 : Default template parameters for configurability added
//    Niramay Sanghvi : Aug 07 2017 : Used ac_shift_left function for RND and SAT support.
//    Niramay Sanghvi : Jul 12 2017 : Added header style format.
//    Niramay Sanghvi : Jul 11 2017 : Added support for all possible values of integer widths.
//    Niramay Sanghvi : Jul 05 2017 : Passed output by reference.
//    Niramay Sanghvi : Jun 29 2017 : Renamed header files and functions.
//
//********************************************************************************************

#ifndef _INCLUDED_AC_POW_PWL_H_
#define _INCLUDED_AC_POW_PWL_H_

// The functions use default template parameters, which are only supported by C++11 or later
// compiler standards. Hence, the user should be informed if they are not using those standards.

#if __cplusplus < 201103L
#error Please use C++11 or a later standard for compilation.
#endif

#include <ac_int.h>
#include <ac_fixed.h>

// Include headers for required functions
#include <ac_math/ac_shift.h>
#include <ac_math/ac_log_pwl.h>

#if !defined(__SYNTHESIS__)
#include <iostream>
#endif

//=========================================================================
// Function: ac_pow2_pwl (for ac_fixed)
//
// Description:
//    Calculation of base 2 exponential of real inputs, passed as ac_fixed
//    variables.
//
//    Separates input into integer and fractional part, the fractional part
//    is passed to the PWL approximation. The output is then left-shifted
//    by the value of the integer part, in order to de-normalize.
//
// Usage:
//    See above example code for usage.
//
// Notes:
//    The PWL implementation utilizes 3 elements, which has a small impact
//    on accuracy.
//    Function only supports unsigned output types.
//
//-------------------------------------------------------------------------

namespace ac_math
{
  template<ac_q_mode pwl_Q = AC_TRN,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  void ac_pow2_pwl(
    const ac_fixed<W, I, S, Q, O> &input,
    ac_fixed<outW, outI, false, outQ, outO> &output
  )
  {
    // Stores the fractional part of the input. By default it is set to 0
    ac_fixed<AC_MAX(W - I, 1), 0, false> input_frac_part = 0.0;

    // Take out the fractional part of the input
    // This serves as a sort of normalization, with the fractional part being
    // the normalized data (can only vary from 0 to 0.9999...)

    // Only carry out slicing if the input has a fractional component.
    // If the input doesn't have a fractional part, the default value of input_frac_part, i.e. 0,
    // is suitable to be used in later calculations.
#pragma hls_waive CNS
    if (W > I) {input_frac_part.set_slc(0, input.template slc<AC_MAX(W - I, 1)>(0));}

    // Start of code outputted by ac_pow_pwl_lutgen.cpp
    // Note that the LUT generator file also outputs values for x_min_lut (lower limit of PWL domain), x_max_lut (upper limit of PWL domain)
    // and sc_constant_lut (scaling factor used to scale the input from 0 to n_segments_lut). However, these values aren't considered in the header
    // file because it has been optimized to work with a 4-segment PWL model that covers the domain of [0, 1). For other PWL implementations, the user will probably have
    // to take these values into account explicitly. Guidelines for doing so are given in the comments.
    // In addition, some of the slope values here are modified slightly in order to ensure monotonicity of the PWL function as the input crosses segment boundaries.
    // The user might want to take care to ensure that for their own PWL versions.

    // Initialization for PWL LUT
    const unsigned n_segments_lut = 4;
    const int int_bits = ac::nbits<n_segments_lut - 1>::val;
    // The number of fractional bits for the LUT values is chosen by first finding the maximum absolute error over the domain of the PWL
    // when double-precision values are used for LUT values. This error will correspond to a number of fractional bits that are always
    // guaranteed to be error-free, for fixed-point PWL outputs.
    // This number of fractional bits is found out by the formula:
    // nbits = abs(ceil(log2(abs_error_max)).
    // The number of fractional bits hereafter used to store the LUT values is nbits + 2.
    // For this particular PWL implementation, the number of fractional bits is 9.
    // Initializing the LUT arrays
    const int n_frac_bits = 10;
    // The intermediate values for the scaled input will have the number of fractional bits set to n_frac_bits or the number of fractional bits in the
    // input minus 2, whichever is lower. The subtraction of 2 is done in order to take the scaling-related left-shift by 2 into account.
    // This value for the number of fractional bits will changed if the PWL implementation changes; the value provided by default works with an 4-segment
    // PWL that covers the domain of [0, 1)
    // In order to ensure that the intermediate variables always have a length >= 1 and thereby prevent a compile-time error, AC_MAX is used.
    const int sc_input_frac_bits = AC_MAX(1, AC_MIN(n_frac_bits, W - I - 2));
    const ac_fixed<n_frac_bits, 0, false> m_lut[n_segments_lut] = {.189453125, .224609375, 0.2666015625, .3173828125};
    const ac_fixed<n_frac_bits + 1, 1, false> c_lut[n_segments_lut] = {.998046875, 1.1875, 1.412109375, 1.6787109375};

    // End of code outputted by ac_pow_pwl_lutgen.cpp

    // Compute power of two using pwl
    // Scale the normalized input from 0 to n_segments_lut. Any other PWL implementation
    // with a different number of segments/domain should be scaled according to the formula: x_in_sc = (input_frac_part - x_min_lut) * sc_constant_lut
    // where sc_constant_lut = n_segments_lut / (x_max_lut - x_min_lut)
    // (x_min_lut and and x_max_lut are the lower and upper limits of the domain)
    ac_fixed<sc_input_frac_bits + int_bits, int_bits, false> x_in_sc = ((ac_fixed<sc_input_frac_bits + int_bits + 2, int_bits, false>)input_frac_part) << 2;
    ac_fixed<sc_input_frac_bits, 0, false> x_in_sc_frac;
    // Slice out the fractional part from the scaled input, store it in another variable.
    x_in_sc_frac.set_slc(0, x_in_sc.template slc<sc_input_frac_bits>(0));
    // The integer part of the scaled input is the index of the LUT table
    ac_int<int_bits, false> index = x_in_sc.to_int();
    typedef ac_fixed<sc_input_frac_bits + n_frac_bits + 1, 1, false, pwl_Q> output_pwl_type;
    output_pwl_type output_pwl = m_lut[index] * x_in_sc_frac + c_lut[index];

    // Shift left by the integer part of the input to cancel out the previous normalization.
    ac_math::ac_shift_left(output_pwl, input.to_int(), output);

#if !defined(__SYNTHESIS__) && defined(AC_POW_PWL_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "Actual input              = " << input << std::endl;
    std::cout << "normalized input          = " << input_frac_part << std::endl;
    std::cout << "output up-scaled by exp   = " << output << std::endl;
    std::cout << "index                     = " << index  << std::endl;
#endif
  }

//=============================================================================
// Version that allows the return of values.
  template<class T_out, ac_q_mode pwl_Q = AC_TRN, class T_in>
  T_out ac_pow2_pwl(
    const T_in &input
  )
  {
    // Create a temporary variable for output and use the pass-by-reference version
    // to evaluate it. This temporary variable is returned as the output.
    T_out output;
    ac_pow2_pwl<pwl_Q>(input, output);
    return output;
  }

//=============================================================================
// Function: ac_exp_pwl (for ac_fixed)
//
// Description:
//    Calculation of base e exponential of real inputs, passed as ac_fixed
//    variables.
//
// Usage:
//    A sample testbench and its implementation look like this:
//
//    #include <ac_math/ac_pow_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_fixed<20, 11, true, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<24, 14, false, AC_RND, AC_SAT> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_exp_pwl(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 2.5;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN(0);
//    }
//    #endif
//
// Notes:
//    This function relies on the ac_pow2_pwl function for its computation. It
//    does this by multiplying the input with log2(e), then passing it to
//    the ac_pow2_pwl function. In doing so, we also make sure that the
//    product variable has enough precision to store the result of
//    input*log2(e).
//    Function only supports unsigned output types.
//
//-----------------------------------------------------------------------------

  // This struct computes precision of input_inter variable ("pii") for base e exponent. It also makes sure
  // that there are a set minimum no. of fractional bits to represent the multiplication of x with log2(e)
  // (this is decided by the n_f_b variable).
  template <int W, int I, bool S, ac_q_mode Q, ac_o_mode O, int n_f_b>
  struct comp_pii_exp {
    enum {
      pit_i       = I + 1,
      pit_w_inter = W + 1,
      pit_w       = (W - I) > n_f_b ? pit_w_inter : pit_i + n_f_b
    };
    typedef ac_fixed<pit_w, pit_i, S, Q, O> pit_t;
  };

  //n_f_b = minimum no of fractional bits used in storing the result of multiplication by log2(e)
  template<int n_f_b = 11, ac_q_mode pwl_Q = AC_TRN,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  void ac_exp_pwl(
    const ac_fixed<W, I, S, Q, O> &input,
    ac_fixed<outW, outI, false, outQ, outO> &output
  )
  {
    const ac_fixed<17, 3, true> log2e = 1.44269504089;
    // Find type of intermediate variable used to store output of x*log2(e)
    typedef typename comp_pii_exp<W, I, S, Q, O, n_f_b>::pit_t input_inter_type;
    input_inter_type input_inter;
    // e^x = 2^(x*log2(e))
    input_inter = input*log2e;
    ac_pow2_pwl<pwl_Q>(input_inter, output);

#if !defined(__SYNTHESIS__) && defined(AC_POW_PWL_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "input_inter.width       = " << input_inter.width << std::endl;
    std::cout << "input_inter.i_width     = " << input_inter.i_width << std::endl;
    std::cout << "input (power_exp)       = " << input << std::endl;
    std::cout << "log2e (power_exp)       = " << log2e << std::endl;
    std::cout << "input_inter (power_exp) = " << input_inter << std::endl;
    std::cout << "output (power_exp)      = " << output << std::endl;
#endif
  }

//=============================================================================
// Version that allows the return of values.
  template<class T_out, int n_f_b = 11, ac_q_mode pwl_Q = AC_TRN, class T_in>
  T_out ac_exp_pwl(
    const T_in &input
  )
  {
    // Create a temporary variable for output and use the pass-by-reference version
    // to evaluate it. This temporary variable is returned as the output.
    T_out output;
    ac_exp_pwl<n_f_b, pwl_Q>(input, output);
    return output;
  }

//=============================================================================
// Function: ac_pow_pwl (for ac_fixed)
//
// Description:
//    Calculation of exponentials with any base, for ac_fixed variables.
//
// Usage:
//    A sample testbench and its implementation look like this:
//
//    #include <ac_math/ac_pow_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_fixed<20, 11, false, AC_RND, AC_SAT> base_type;
//    typedef ac_fixed<21, 12, true, AC_RND, AC_SAT> expon_type;
//    typedef ac_fixed<24, 14, false, AC_RND, AC_SAT> output_type;
//
//    #pragma hls_design top
//    void project(
//      const base_type &base,
//      const expon_type  &expon,
//      output_type &output
//    )
//    {
//      ac_pow_pwl(base, expon, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      base_type base = 2.5;
//      expon_type expon = 2;
//      output_type output;
//      CCS_DESIGN(project)(base, expon, output);
//      CCS_RETURN(0);
//    }
//    #endif
//
// Notes:
//    This function relies on the ac_pow2_pwl and ac_log_pwl functions for its
//    computation. It does this by multiplying expon with log2(base), then
//    passing it to the ac_pow2_pwl function. In doing so, we also make sure
//    that the product variable has enough precision to store the result of
//    expon*log2(base).
//    Input for the base and output for the exponential value have to be
//    unsigned ac_fixed variables.
//
//-----------------------------------------------------------------------------

  template<ac_q_mode pwl_Q = AC_TRN,
           int baseW, int baseI, ac_q_mode baseQ, ac_o_mode baseO,
           int exponW, int exponI, bool exponS, ac_q_mode exponQ, ac_o_mode exponO,
           int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  void ac_pow_pwl(
    const ac_fixed<baseW, baseI, false, baseQ, baseO> &base,
    const ac_fixed<exponW, exponI, exponS, exponQ, exponO> &expon,
    ac_fixed<outW, outI, false, outQ, outO> &output
  )
  {
    // Find the number of integer bits required to represent the minimum and maximum values expressable for log2 of the base.
    // The number of integer bits used for the temporary variable that stores log2(base) is whichever is larger + 1.
    const int t_I_frac = ac::nbits<AC_MAX(baseW - baseI, 0)>::val;
    const int t_I_int  = ac::nbits<AC_MAX(baseI, 0)>::val;
    const int t_I      = (t_I_frac > t_I_int ? t_I_frac : t_I_int) + 1;
    // Store the number of fractional bits for the temp output that ensures losslessness. This can change based on the log2
    // PWL implementation, hence, the user must handle these changes appropriately.
    const int n_f_b_pwl_out = 22;
    ac_fixed <n_f_b_pwl_out + t_I, t_I, true, pwl_Q> log2_base;
    // Find log2(base)
    ac_math::ac_log2_pwl(base, log2_base);
    // Multiply expon by log2(base) and pass it to the ac_pow2_pwl function
    ac_pow2_pwl(expon*log2_base, output);

#if !defined(__SYNTHESIS__) && defined(AC_POW_PWL_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "log2_base          = " << log2_base << std::endl;
    std::cout << "output(ac_pow_pwl) = " << output << std::endl;
#endif
  }

//=============================================================================
// Version that allows the return of values.
  template<class T_out, ac_q_mode pwl_Q = AC_TRN, class T_in_base, class T_in_expon>
  T_out ac_pow_pwl(
    const T_in_base  &base,
    const T_in_expon &expon
  )
  {
    // Create a temporary variable for output and use the pass-by-reference version
    // to evaluate it. This temporary variable is returned as the output.
    T_out output;
    ac_pow_pwl<pwl_Q>(base, expon, output);
    return output;
  }
}

#endif

