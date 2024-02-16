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
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//    3.4.2 - [CAT-29828] Added floating point support for ac_pow2_pwl() and ac_exp_pwl().
//    3.3.0 - [CAT-25797] Added CDesignChecker fixes/waivers for code check violations in ac_math PWL and Linear Algebra IPs.
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
#if (defined(__GNUC__) && (__cplusplus < 201103L))
#error Please use C++11 or a later standard for compilation.
#endif
#if (defined(_MSC_VER) && (_MSC_VER < 1920) && !defined(__EDG__))
#error Please use Microsoft VS 2019 or a later standard for compilation.
#endif

#include <ac_int.h>
// Include headers for data types supported by these implementations
#include <ac_fixed.h>
#include <ac_float.h>
#include <ac_std_float.h>

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
    if (W > I) {
      input_frac_part.set_slc(0, input.template slc<AC_MAX(W - I, 1)>(0));
    }

    // Start of code outputted by ac_pow_pwl_lutgen.cpp

    const unsigned n_segments_lut = 4; // Number of PWL segments.
    const int n_frac_bits = 10; // Number of fractional bits
    // Since scaling constant is a positive power-of-two, multiplication with it is the same as left-shifting by 2.
    // Accordingly, the scaled normalized input will have 2 less fractional bits than the normalized input, provided that this
    // number of fractional bits is lesser than n_frac_bits. If not, the number of fractional bits in the scaled input is set to n_frac_bits.
    const int sc_input_frac_bits = AC_MAX(1, AC_MIN(n_frac_bits, W - I - 2));
    // Slope and intercept LUT values.
    const ac_fixed<0 + n_frac_bits, 0, false> m_lut[n_segments_lut] = {.189453125, .224609375, .2666015625, .3173828125};
    const ac_fixed<2 + n_frac_bits, 2, false> c_lut[n_segments_lut] = {.998046875, 1.1875, 1.412109375, 1.6787109375};
    const ac_fixed<1 + n_frac_bits, 1, false> x_min_lut = .0; // Minimum limit of PWL domain
    const ac_fixed<3 + n_frac_bits, 3, false> sc_constant_lut = 4.0; // Scaling constant

    // End of code outputted by ac_pow_pwl_lutgen.cpp

    const int int_bits = ac::nbits<n_segments_lut - 1>::val;
    // Compute power of two using pwl
    // Scale the normalized input to a value that ranges from 0 to n_segments_lut.
    ac_fixed<sc_input_frac_bits + int_bits, int_bits, false> x_in_sc = (input_frac_part - x_min_lut)*sc_constant_lut;
    ac_fixed<sc_input_frac_bits, 0, false> x_in_sc_frac;
    // Slice out the fractional part from the scaled input, store it in another variable.
		#pragma hls_waive UMR
    x_in_sc_frac.set_slc(0, x_in_sc.template slc<sc_input_frac_bits>(0));
    // The integer part of the scaled input is the index of the LUT table
    ac_int<int_bits, false> index = x_in_sc.to_int();
    typedef ac_fixed<sc_input_frac_bits + n_frac_bits + 1, 1, false, pwl_Q> output_pwl_type;
    output_pwl_type output_pwl = m_lut[index] * x_in_sc_frac + c_lut[index];

    #pragma hls_waive CNS
    if (!S && I <= 0) {
      // No de-normalization is required if the input is unsigned and has no integer bits.
      output = output_pwl;
    } else {
      // Shift left by the integer part of the input to cancel out the previous normalization.
      ac_math::ac_shift_left(output_pwl, input.to_int(), output);
    }

    #if !defined(__SYNTHESIS__) && defined(AC_POW_PWL_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "Actual input              = " << input << std::endl;
    std::cout << "normalized input          = " << input_frac_part << std::endl;
    std::cout << "output up-scaled by exp   = " << output << std::endl;
    std::cout << "index                     = " << index  << std::endl;
    #endif
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
  template <int W, int I, bool S, int n_f_b>
  struct comp_pii_exp {
    enum {
      pit_i       = I + 1,
      pit_w_inter = W + 1,
      pit_w       = (W - I) > n_f_b ? pit_w_inter : pit_i + n_f_b
    };
    typedef ac_fixed<pit_w, pit_i, S> pit_t;
  };

  // n_f_b = If the "AC_POW_PWL_CHANGE_FRAC_BITS" macro is defined, this becomes the minimum no of 
  // fractional bits used in storing the result of multiplication by log2(e).
  // If the macro is not defined, this template parameter is not used in the code.
  template<int n_f_b = 11, ac_q_mode pwl_Q = AC_TRN,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  void ac_exp_pwl(
    const ac_fixed<W, I, S, Q, O> &input,
    ac_fixed<outW, outI, false, outQ, outO> &output
  )
  {
    const ac_fixed<17, 3, true> log2e = 1.44269504089;
    // n_segments_lut and n_frac_bits should be the same as the corresponding values in the ac_fixed
    // implementation above.
    // The bitwidth calculations which use these constants are designed with a PWL domain of [0, 1) in
    // mind. They will still work for other PWL domains, but at a possibly lower-than-optimal accuracy.
    //
    // NOTE: Change these constants if you change the either of the corresponding values in the
    // ac_fixed implementation above.
    const unsigned n_segments_lut = 4;
    const int n_frac_bits = 10;
    
    #ifdef AC_POW_PWL_CHANGE_FRAC_BITS
    // Find type of intermediate variable used to store output of x*log2(e)
    typedef class comp_pii_exp<W, I, S, n_f_b>::pit_t input_inter_type;
    #else
    const bool is_n_seg_po2 = !bool(n_segments_lut & (n_segments_lut - 1));
    const int extra_f_bits = is_n_seg_po2 ? ac::nbits<n_segments_lut - 1>::val : 0;
    // Find type of intermediate variable used to store output of x*log2(e)
    typedef class comp_pii_exp<W, I, S, n_frac_bits + extra_f_bits>::pit_t input_inter_type;
    #endif
    
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

//=================================================================================
// Helper class: ac_generic_float_pow_pwl
//
// Description:
//   Helper class that contains interace() method that lets us calculate the pow
//   outputs with different bases. The ac_pow2_pwl and ac_exp_pwl functions for
//   floating point inputs are friend functions to this class.
//
//   The selection between base 2 and base e exponential outputs is done via a
//   template argument which accepts a "base" enum from the "float_base" struct
//   in ac_log_pwl.h. This template parameter must be supplied in the while
//   declaring an object of the class.
//
// Notes:
//   The actual PWL calculation is done by the ac_fixed ac_pow2_pwl implementation.
//
//   All the class functions, including the constructor, are private. Objects of
//   the class can only be declared and accessed through the friend functions.
//
//---------------------------------------------------------------------------------

  // float_base: struct defined in ac_log_pwl.h
  template <ac_q_mode pwl_Q_, enum float_base::base base_val>
  class ac_generic_float_pow_pwl
  {
  public:
    template <ac_q_mode pwl_Q, int W, int I, int E, ac_q_mode Q, int outW, int outI, int outE, ac_q_mode outQ>
    friend void ac_pow2_pwl(const ac_float<W, I, E, Q> input, ac_float<outW, outI, outE, outQ> &output);

    template <ac_q_mode pwl_Q, int W, int I, int E, ac_q_mode Q, int outW, int outI, int outE, ac_q_mode outQ>
    friend void ac_exp_pwl(const ac_float<W, I, E, Q> input, ac_float<outW, outI, outE, outQ> &output);

  private:
    template <int W, int I, int E, ac_q_mode Q, int outW, int outI, int outE, ac_q_mode outQ>
    void interface(const ac_float<W, I, E, Q> input, ac_float<outW, outI, outE, outQ> &output) {
      const int out_exp_min = 1 << (outE - 1);
      const int I_fxpt = ac::nbits<out_exp_min + AC_MAX(outI - 1, outW - outI)>::val + 1;
      // n_segments_lut and n_frac_bits should be the same as the corresponding values in the ac_fixed
      // implementation above.
      // The bitwidth calculations which use these constants are designed with a PWL domain of [0, 1) in
      // mind. They will still work for other PWL domains, but at a possibly lower-than-optimal accuracy.
      //
      // NOTE: Change these constants if you change the either of the corresponding values in the
      // ac_fixed implementation above.
      const unsigned n_segments_lut = 4;
      const int n_frac_bits = 10;
      
      const bool is_n_seg_po2 = !bool(n_segments_lut & (n_segments_lut - 1));
      const int extra_f_bits = is_n_seg_po2 ? ac::nbits<n_segments_lut - 1>::val : 0;
      const int W_fxpt = I_fxpt + n_frac_bits + extra_f_bits;

      ac_fixed<W, I, true> input_mant = input.mantissa();
      int input_exp = input.exp().to_int();
      // Rounding mode of input_to_pow2_fxpt is defined by pwl_Q_ template parameter.
      ac_fixed<W_fxpt, I_fxpt, true, pwl_Q_, AC_SAT> input_to_pow2_fxpt;

      #pragma hls_waive CNS
      if (base_val == float_base::base_2) {
        // Convert ac_float to ac_fixed variable, through the following left-shift operation:
        // mantissa << exp. The precision of input_to_pow2_fxpt is such that it will saturate a little
        // above the input which would cause the floating point output to saturate if we were using an
        // accurate math function to calculate the pow2 output.
        ac_math::ac_shift_left(input_mant, input_exp, input_to_pow2_fxpt);
      } else {
        // This step utilizes almost the same left-shifting operation as the one for base 2 exponentials,
        // except the input mantissa is also multiplied by log2(e).
        const ac_fixed<17, 3, true> log2e = 1.44269504089;
        ac_math::ac_shift_left(input_mant*log2e, input_exp, input_to_pow2_fxpt);
      }

      ac_fixed<n_frac_bits + extra_f_bits, 0, false> input_to_pow2_fxpt_frac;
      input_to_pow2_fxpt_frac.set_slc(0, input_to_pow2_fxpt.template slc<n_frac_bits + extra_f_bits>(0));

      ac_fixed<2*n_frac_bits + 1, 1, false> output_fxpt_frac;

      ac_pow2_pwl<pwl_Q_>(input_to_pow2_fxpt_frac, output_fxpt_frac);
      ac_int<I_fxpt, true> input_to_pow2_fxpt_int = input_to_pow2_fxpt.to_int();

      ac_float<outW, outI, outE, outQ> output_temp(output_fxpt_frac, input_to_pow2_fxpt_int, true);

      output = output_temp;

      #if !defined(__SYNTHESIS__) && defined(AC_POW_PWL_H_DEBUG)
      std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
      std::cout << "Floating point input = " << input << std::endl;
      std::cout << "input_to_pow2_fxpt.type_name() = " << input_to_pow2_fxpt.type_name() << std::endl;
      std::cout << "input_to_pow2_fxpt = " << input_to_pow2_fxpt << std::endl;
      std::cout << "input_to_pow2_fxpt_int = " << input_to_pow2_fxpt_int << std::endl;
      std::cout << "input_to_pow2_fxpt_frac = " << input_to_pow2_fxpt_frac << std::endl;
      std::cout << "output_fxpt_frac = " << output_fxpt_frac << std::endl;
      std::cout << "output_temp = " << output_temp << std::endl;
      #endif
    }

    ac_generic_float_pow_pwl() {}
  };

//=========================================================================
// Function: ac_pow2_pwl (for ac_float)
//
// Description:
//    Calculation of base 2 exponential of real inputs, passed as ac_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_math/ac_pow_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_float<20, 5, 6> input_type;
//    typedef ac_float<32, 5, 10> output_type;
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
//
//      CCS_DESIGN(project)(input, output);
//
//      CCS_RETURN(0);
//    }
//    #endif
//
// Note:
//    The ac_float version of ac_pow2_pwl is a friend function to the
//    ac_generic_float_pow_pwl class, and depends on the interface() function in
//    said class to call the ac_fixed version of the ac_pow2_pwl function and
//    calculate the final output.
//
//-------------------------------------------------------------------------

  template <ac_q_mode pwl_Q = AC_TRN,
            int W, int I, int E, ac_q_mode Q,
            int outW, int outI, int outE, ac_q_mode outQ>
  void ac_pow2_pwl(
    const ac_float<W, I, E, Q> input,
    ac_float<outW, outI, outE, outQ> &output
  )
  {
    ac_generic_float_pow_pwl<pwl_Q, float_base::base_2> gen_pow_obj;
    gen_pow_obj.interface(input, output);
  }

//=============================================================================
// Function: ac_exp_pwl (for ac_float)
//
// Description:
//    Calculation of base e exponential of real inputs, passed as ac_float
//    variables.
//
// Usage:
//    #include <ac_math/ac_pow_pwl.h>
//    using namespace ac_math;
//    
//    typedef ac_float<16, 5, 6> input_type;
//    typedef ac_float<20, 5, 10> output_type;
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
//    The ac_float version of ac_exp_pwl is a friend function to the
//    ac_generic_float_pow_pwl class, and depends on the interface() function in
//    said class to call the ac_fixed version of the ac_pow2_pwl function and
//    calculate the final output.
//
//-----------------------------------------------------------------------------

  template <ac_q_mode pwl_Q = AC_TRN,
            int W, int I, int E, ac_q_mode Q,
            int outW, int outI, int outE, ac_q_mode outQ>
  void ac_exp_pwl(
    const ac_float<W, I, E, Q> input,
    ac_float<outW, outI, outE, outQ> &output
  ) {
    ac_generic_float_pow_pwl<pwl_Q, float_base::base_e> gen_pow_obj;
    gen_pow_obj.interface(input, output);
  }

//=========================================================================
// Function: ac_pow2_pwl (for ac_std_float)
//
// Description:
//    Calculation of base 2 exponential of real inputs, passed as
//    ac_std_float variables.
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_math/ac_pow_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_std_float<32, 8> input_type;
//    typedef ac_std_float<32, 8> output_type;
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
//      input_type input(2.5);
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN(0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template<ac_q_mode pwl_Q = AC_TRN,
           int W, int E, int outW, int outE>
  void ac_pow2_pwl(
    const ac_std_float<W, E> &input,
    ac_std_float<outW, outE> &output
  )
  {
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_pow2_pwl<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_std_float<outW, outE> output_temp(output_ac_fl); // Convert output ac_float to ac_std_float.
    
    #ifdef AC_POW_PWL_NAN_SUPPORTED
    // If the input is +/-nan, the output will be +nan. This mirrors the behavior of the
    // pow library to produce the same kind of output as ac_pow2_pwl, i.e. pow(2.0, input).
    if (input.isnan()) {
      output_temp = output_temp.nan();
    }
    #endif
    
    output = output_temp;
  }

//=========================================================================
// Function: ac_exp_pwl (for ac_std_float)
//
// Description:
//    Calculation of base e exponential of real inputs, passed as
//    ac_std_float variables.
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_math/ac_pow_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_std_float<32, 8> input_type;
//    typedef ac_std_float<32, 8> output_type;
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
//      input_type input(2.5);
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN(0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template<ac_q_mode pwl_Q = AC_TRN,
           int W, int E, int outW, int outE>
  void ac_exp_pwl(
    const ac_std_float<W, E> &input,
    ac_std_float<outW, outE> &output
  )
  {
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_exp_pwl<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_std_float<outW, outE> output_temp(output_ac_fl); // Convert output ac_float to ac_std_float.
    
    #ifdef AC_POW_PWL_NAN_SUPPORTED
    // If the input is +/-nan, the output is set to +/-nan, depending on the signedness of the input.
    // This mirrors the functioning of the exp(input) function in math.h.
    if (input.isnan()) {
      output_temp = output_temp.nan();
      output_temp.set_signbit(input.signbit());
    }
    #endif
    
    output = output_temp;
  }

//=========================================================================
// Function: ac_pow2_pwl (for ac_ieee_float)
//
// Description:
//    Calculation of base 2 exponential of real inputs, passed as
//    ac_ieee_float variables.
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_math/ac_pow_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_ieee_float<binary32> input_type;
//    typedef ac_ieee_float<binary32> output_type;
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
//      input_type input(2.5);
//      output_type output;
//
//      CCS_DESIGN(project)(input, output);
//
//      CCS_RETURN(0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template<ac_q_mode pwl_Q = AC_TRN,
           ac_ieee_float_format Format,
           ac_ieee_float_format outFormat>
  void ac_pow2_pwl(
    const ac_ieee_float<Format> &input,
    ac_ieee_float<outFormat> &output
  )
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int outW = T_out::width;
    const int outE = T_out::e_width;
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_pow2_pwl<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_ieee_float<outFormat> output_temp(output_ac_fl); // Convert output ac_float to ac_ieee_float.
    
    #ifdef AC_POW_PWL_NAN_SUPPORTED
    // If the input is +/-nan, the output will be +nan. This mirrors the behavior of the
    // pow library to produce the same kind of output as ac_pow2_pwl, i.e. pow(2.0, input).
    if (input.isnan()) {
      output_temp = output_temp.nan();
    }
    #endif
    
    output = output_temp;
  }

//=========================================================================
// Function: ac_exp_pwl (for ac_ieee_float)
//
// Description:
//    Calculation of base e exponential of real inputs, passed as
//    ac_ieee_float variables.
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_math/ac_pow_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_ieee_float<binary32> input_type;
//    typedef ac_ieee_float<binary32> output_type;
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
//      input_type input(2.5);
//      output_type output;
//
//      CCS_DESIGN(project)(input, output);
    // If the input is +/-nan, the output is set to +/-nan, depending on the signedness of the input.
    // This mirrors the functioning of the exp(input) function in math.h.
//
//      CCS_RETURN(0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template<ac_q_mode pwl_Q = AC_TRN,
           ac_ieee_float_format Format,
           ac_ieee_float_format outFormat>
  void ac_exp_pwl(
    const ac_ieee_float<Format> &input,
    ac_ieee_float<outFormat> &output
  )
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int outW = T_out::width;
    const int outE = T_out::e_width;
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_exp_pwl<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_ieee_float<outFormat> output_temp(output_ac_fl); // Convert output ac_float to ac_ieee_float.
    
    #ifdef AC_POW_PWL_NAN_SUPPORTED
    // If the input is +/-nan, the output is set to +/-nan, depending on the signedness of the input.
    // This mirrors the functioning of the exp(input) function in math.h.
    if (input.isnan()) {
      output_temp = output_temp.nan();
      output_temp.set_signbit(input.signbit());
    }
    #endif
    
    output = output_temp;
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
// Version that allows the return of values for variables that aren't ac_fixed
  template<class T_out, ac_q_mode pwl_Q = AC_TRN, class T_in>
  T_out ac_exp_pwl(
    const T_in &input
  )
  {
    // Create a temporary variable for output and use the pass-by-reference version
    // to evaluate it. This temporary variable is returned as the output.
    T_out output;
    ac_exp_pwl<pwl_Q>(input, output);
    return output;
  }

//=============================================================================
// Version that allows the return of values for variables that are ac_fixed
  template <class T_out, int n_f_b = 11, ac_q_mode pwl_Q = AC_TRN,
            int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
  T_out ac_exp_pwl(
    const ac_fixed<W, I, S, Q, O> &input
  ) {
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
//      expon_type expon = 2.75;
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
    // n_frac_bits should be the same as the n_frac_bits constant in the ac_fixed ac_log2_pwl
    // implementation.
    // NOTE: Change this quantity if you change the n_frac_bits value in said implementation.
    const int n_frac_bits = 11;
    const int n_f_b_pwl_out = 2*n_frac_bits;
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
  // Version that allows the return of values, for fixed point types.
  template<class T_out, ac_q_mode pwl_Q = AC_TRN,
           class T_base, class T_expon>
  T_out ac_pow_pwl(
    const T_base &base,
    const T_expon &expon
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

