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
//***********************************************************************************
// File: ac_log_pwl.h
//
// Description: Provides piece-wise linear implementations of the
//    log function for the AC (tm) Datatypes: ac_fixed.
//    Two different functions compute values with bases as 2 and e.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_log_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_fixed<16, 8, false, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<16, 8, true, AC_RND, AC_SAT> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_log2_pwl(input,output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 1.25;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    This file uses the ac_normalize() function from ac_normalize.h
//
// Revision History:
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//    3.3.0  - [CAT-25797] Added CDesignChecker fixes/waivers for code check violations in ac_math PWL and Linear Algebra IPs.
//             Waivers added for CNS and CCC violations.
//             Fixes added for FXD, STF and MXS violations.
//               - FXD violations fixed by changing integer literals to floating point literals or typecasting to ac_fixed values.
//               - STF violations fixed by using "const" instead of "static const" parameters. LUT generator files also print out "const" LUTs instead of "static const" LUTs.
//               - MXS violations fixed by typecasting unsigned variables to int.
//    3.2.4  - CAT-24697, Added floating point support.
//    3.1.2  - Improved bitwidth calculations for ac_fixed.
//    2.0.10 - Official open-source release as part of the ac_math library.
//
//***********************************************************************************

#ifndef _INCLUDED_AC_LOG_PWL_H_
#define _INCLUDED_AC_LOG_PWL_H_

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

// Include headers for required functions
#include <ac_math/ac_normalize.h>

#if !defined(__SYNTHESIS__) && defined(AC_LOG_PWL_H_DEBUG)
#include <iostream>
#endif

namespace ac_math
{
//=================================================================================
// Function: ac_log2_pwl (for ac_fixed, returns log2(input) )
//
// Description:
//    Calculation of log base 2 of real inputs, passed as ac_fixed
//    variables.
//
// Usage:
//    See above example code for usage.
//
// Notes:
//    The PWL approximation uses 8 segments, which has a small impact on accuracy.
//    If the user is sure that the input is already normalized, they can set the
//    call_normalize function parameter to false (the parameter is set to
//    true by default)
//
//---------------------------------------------------------------------------------

  template <ac_q_mode pwlQ = AC_TRN, int W, int I, ac_q_mode Q, ac_o_mode O, int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO>
  void ac_log2_pwl (const ac_fixed <W, I, false, Q, O> input, ac_fixed <outW, outI, outS, outQ, outO> &output, const bool call_normalize = true)
  {
    // Use a macro to activate the AC_ASSERT
    // If AC_ASSERT is activated: the program will stop running as soon as a zero input
    // is encountered.
    // If AC_ASSERT is not activated: the output will saturate when a zero input is encountered.
    // The functionality behind this is taken care of by other sections of the code.
    #ifdef ASSERT_ON_INVALID_INPUT
    AC_ASSERT(!!input, "Log of zero not supported.");
    #endif

    ac_fixed<W, 0, false, Q, O> input_normalized;

    // exp is used to store the final exponent value
    int exp;

    // If call_normalize is set to true, the design does not assume that the input is already normalized and performs normalization by calling ac_normalize().
    // If call_normalize is set to false, the design assumes that the input is already normalized and does not call ac_normalize, thereby saving on hardware.
    #pragma hls_waive CNS
    if (call_normalize) {
      exp = ac_math::ac_normalize(input, input_normalized);
    } else {
      exp = I;
      input_normalized.set_slc(0, input.template slc<W>(0));
    }

    // Start of code outputted by ac_log_pwl_lutgen.cpp

    const unsigned n_segments_lut = 8;
    const int n_frac_bits = 11;
    const ac_fixed<-1 + n_frac_bits, -1, false> m_lut[n_segments_lut] = {.169921875, .1513671875, .13720703125, .125, .115234375, .1064453125, .099609375, .09326171875};
    const ac_fixed<2 + n_frac_bits, 2, true> c_lut[n_segments_lut] = {-.99853515625, -.82861328125, -.67724609375, -.53955078125, -.41455078125, -.298828125, -.1923828125, -.0927734375};
    const ac_fixed<0 + n_frac_bits, 0, false> x_min_lut = .5;
    const ac_fixed<1 + n_frac_bits, 1, false> x_max_lut = 1.0;
    const ac_fixed<5 + n_frac_bits, 5, false> sc_constant_lut = 16.0;

    // End of code outputted by ac_log_pwl_lutgen.cpp

    const int int_bits = ac::nbits<n_segments_lut - 1>::val;

    // Scale the normalized input from 0 to n_segments_lut
    ac_fixed <n_frac_bits + int_bits, int_bits, false> input_sc = (input_normalized - x_min_lut)*sc_constant_lut;
    // Take out the fractional bits of the scaled input
    ac_fixed<n_frac_bits, 0, false> input_sc_frac;
    input_sc_frac.set_slc(0, input_sc.template slc<n_frac_bits>(0));
    // Integer part of scaled input is index
    ac_int <int_bits, false> index = input_sc.to_int();
    // computation of the pwl output
    // The precision given below will ensure that there is no precision lost in the assignment to t, hence rounding for the variable is switched off by default.
    // However, if the user uses less fractional bits and turn rounding on instead, they are welcome to do so by changing giving a different value for pwlQ.
    ac_fixed <2*n_frac_bits + 1, 1, true, pwlQ> t = m_lut[index]*input_sc_frac + c_lut[index];
    // Add the exponent to get the final function output
    ac_fixed <outW, outI, outS, outQ, outO> t2 = t + exp;

    ac_fixed <outW, outI, outS, outQ, outO> output_min;
    // If 0 is supplied as the function input, maximum negative value is returned at the output
    output_min.template set_val<AC_VAL_MIN>();
    // assignment to the final output
    output = (input == 0) ? output_min : t2;

    #if !defined(__SYNTHESIS__) && defined(AC_LOG_PWL_H_DEBUG)
    std::cout << "input.type_name(): " << input.type_name() << std::endl;
    std::cout << "input = " << input << std::endl;
    std::cout << "input to normalization function = " << input << std::endl;
    std::cout << "input_normalized = " << input_normalized << std::endl;
    std::cout << "normalized exp = " << exp << std::endl;
    std::cout << "index of element chosen from ROM = " << index << std::endl;
    std::cout << "final output = " << output << std::endl;
    #endif
  }

//===========================================================================
// Function: ac_log_pwl (for ac_fixed, returns ln(input) )
//
// Description:
//    Calculation of ln (natural log) of real inputs, passed as ac_fixed
//    variables.
//
//    This implementation uses change of base method to compute log base e of
//    input. The value of ln of input is given by,
//    ln(x) = log2(x)/log2(e)
//    which makes, ln(x) = log2(x)*constant
//    log2(x) is computed using above piecewise linear implementation.
//    Note that, 1/log2(e) = loge(2), which is declared as a constant.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_log_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_fixed<16, 8, false, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<16, 8, true, AC_RND, AC_SAT> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_log_pwl(input,output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 1.25;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
//---------------------------------------------------------------------------

  template <ac_q_mode pwlQ = AC_TRN, int W, int I, ac_q_mode Q, ac_o_mode O, int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO>
  void ac_log_pwl (const ac_fixed <W, I, false, Q, O> input, ac_fixed <outW, outI, outS, outQ, outO> &output)
  {
    // Store ln(2) as a constant ac_fixed value.
    const ac_fixed <12, 0, false> ln2 = 0.693145751953125;
    // Find the number of integer bits required to represent the minimum and maximum values expressable for the input type. The number of integer bits
    // used for the temporary variable is whichever is larger.
    const int t_I_frac = ac::nbits<AC_MAX(W - I, 0)>::val;
    const int t_I_int  = ac::nbits<AC_MAX(I, 0)>::val;
    const int t_I      = (t_I_frac > t_I_int ? t_I_frac : t_I_int) + 1;
    const int n_frac_bits = 11; // Number of fractional bits used in ac_fixed PWL implementation. Change this if the ac_fixed PWL implementation changes.
    const int n_f_b_pwl_out = 2*n_frac_bits;
    // The above precision will ensure that the assignment to output_temp is lossless, hence rounding is turned off by default. However, the user is free to use less
    // bits and turn rounding on instead if they wish, by changing the pwlQ variable.
    ac_fixed <n_f_b_pwl_out + t_I, t_I, true, pwlQ> output_temp;
    // call to the log base 2 pwl function
    ac_log2_pwl<pwlQ> (input, output_temp);

    ac_fixed <outW, outI, outS, outQ, outO> output_min;
    // If 0 is supplied as the function input, maximum negative value is returned at the output
    output_min.template set_val<AC_VAL_MIN>();// If input is non-zero, multiply the output of log2(input) by ln(2) to get the final result
    // If input is zero, then assign minimum value to output.
    output = input != 0 ? (ac_fixed <outW, outI, outS, outQ, outO>)(output_temp*ln2) : output_min;

    #if !defined(__SYNTHESIS__) && defined (AC_LOG_PWL_H_DEBUG)
    std::cout << "Input to the log base e function = " << input << std::endl;
    std::cout << "output_temp = " << output_temp << std::endl;
    std::cout << "output = " << output << std::endl;
    #endif
  }

//=================================================================================
// Helper struct: float_base
//
// Description:
//   Helper struct with enum parameters that let the ac_generic_float_log_pwl class
//   below switch between base 2 and base e implementations.
//
//---------------------------------------------------------------------------------

  struct float_base {
    enum base {
      base_2,
      base_e,
    };
  };

//=================================================================================
// Helper class: ac_generic_float_log_pwl
//
// Description:
//   Helper class that contains interace() method that lets us calculate the log
//   outputs with different bases. The ac_log_pwl and ac_log2_pwl functions for
//   floating point inputs are friend functions to this class.
//
//   The selection between log base 2 and log base e outputs is done via a template
//   argument which accepts a "base" enum from the "float_base" struct above. This
//   template parameter must be supplied in the while declaring an object of the
//   class.
//
// Notes:
//   The actual PWL calculation is done by the ac_fixed ac_log2_pwl implementation.
//
//   All the class functions, including the constructor, are private. Objects of
//   the class can only be declared and accessed through the friend functions.
//
//---------------------------------------------------------------------------------

  template <ac_q_mode pwlQ_, enum float_base::base base_val>
  class ac_generic_float_log_pwl
  {
  public:
    template <ac_q_mode pwlQ, int W, int I, int E, ac_q_mode Q, int outW, int outI, int outE, ac_q_mode outQ>
    friend void ac_log2_pwl(const ac_float<W, I, E, Q> input, ac_float<outW, outI, outE, outQ> &output);

    template <ac_q_mode pwlQ, int W, int I, int E, ac_q_mode Q, int outW, int outI, int outE, ac_q_mode outQ>
    friend void ac_log_pwl(const ac_float<W, I, E, Q> input, ac_float<outW, outI, outE, outQ> &output);

  private:
    template <int W, int I, int E, ac_q_mode Q, int outW, int outI, int outE, ac_q_mode outQ>
    void interface(const ac_float<W, I, E, Q> input, ac_float<outW, outI, outE, outQ> &output) {
      // Use a macro to activate the AC_ASSERT
      // If AC_ASSERT is activated, the program will stop running as soon as a negative input
      // is encountered.
      #ifdef ASSERT_ON_INVALID_INPUT
      AC_ASSERT(input >= 0, "Negative input not supported.");
      #endif
      // Store ln(2) as a constant ac_fixed value.
      const ac_fixed <15, 0, false> ln2 =  0.693145751953125;

      const int n_frac_bits = 11; // Number of fractional bits used in ac_fixed PWL implementation. Change this if the ac_fixed PWL implementation changes.

      // Extract input mantissa, discard sign bit and normalize it to the domain of [0.5, 1)
      ac_fixed<W - 1, I - 1, false> input_mant = input.mantissa();
      ac_fixed<W - 1, 0, false> input_mant_norm;
      input_mant_norm.set_slc(0, input_mant.template slc<W - 1>(0));

      ac_fixed<2*n_frac_bits + 1, 1, true> output_mant_norm;
      // Call ac_fixed version to find log of normalized mantissa.
      // ac_float mantissa is already normalized -> ac_fixed version doesn't have to call ac_normalize.
      const bool call_normalize = false;
      ac_log2_pwl<pwlQ_>(input_mant_norm, output_mant_norm, call_normalize);
      ac_float<outW, outI, outE, outQ> output_temp;
      enum {
        I_m_1_neg = bool(I - 1 < 0),
        I_m_1_abs = I_m_1_neg ? -(I - 1) : I - 1,
      };
      typedef ac_int<ac::nbits<I_m_1_abs>::val + int(I_m_1_neg), I_m_1_neg> I_m_1_type;
      // There are two factors combined that result in the ac_fixed input to ac_log2_pwl being normalized: One is the normalization inherent in having a
      // floating point mantissa, and the other is the normalization done while slicing the bit contents of input_mant into input_mant_norm. The
      // denorm_sum variable stores the sum total of these two factors, and is used to denormalize the output.
      typename ac::rt_2T<I_m_1_type, ac_int<E, true> >::plus denorm_sum = I_m_1_type(I - 1) + input.exp(); // Convert (I - 1) to an ac_int(I_m_1_type) variable to reduce area.
      // De-normalize the output.
      #pragma hls_waive CNS
      if (base_val == float_base::base_2) {
        ac_float<outW, outI, outE, outQ> output_temp_b2(output_mant_norm + denorm_sum);
        output_temp = output_temp_b2;
      } else {
        ac_float<outW, outI, outE, outQ> output_temp_be((output_mant_norm + denorm_sum)*ln2); // Multiply the output by ln(2) if base is e.
        output_temp = output_temp_be;
      }

      if (input_mant == 0) {
        output_temp.template set_val<AC_VAL_MIN>(); // If input is 0, set output to the minimum possible value.
      }
      output = output_temp;

      #if !defined(__SYNTHESIS__) && defined(AC_LOG_PWL_H_DEBUG)
      std::cout << "input = " << input << std::endl;
      std::cout << "input_mant_norm = " << input_mant_norm << std::endl;
      std::cout << "output_mant_norm = " << output_mant_norm << std::endl;
      std::cout << "I_m_1_type::type_name() : " << I_m_1_type::type_name() << std::endl;
      std::cout << "input.exp().type_name() : " << input.exp().type_name() << std::endl;
      std::cout << "denorm_sum.type_name() : " << denorm_sum.type_name() << std::endl;
      #endif
    }

    ac_generic_float_log_pwl() {}
  };

//=================================================================================
// Function: ac_log2_pwl (for ac_float, returns log2(input) )
//
// Description:
//    Calculation of log base 2 of real inputs, passed as ac_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_log_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_float<16, 2, 6, AC_TRN> input_type;
//    typedef ac_float<16, 2, 9, AC_TRN> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type input,
//      output_type &output
//    )
//    {
//      ac_log2_pwl(input,output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 1.25;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    The ac_float version of ac_log2_pwl is a friend function to the
//    ac_generic_float_log_pwl class, and depends on the interface() function in
//    said class to call the ac_fixed version of the ac_log2_pwl function and
//    calculate the final output.
//
//---------------------------------------------------------------------------------

  template <ac_q_mode pwlQ = AC_TRN, int W, int I, int E, ac_q_mode Q, int outW, int outI, int outE, ac_q_mode outQ>
  void ac_log2_pwl(const ac_float<W, I, E, Q> input, ac_float<outW, outI, outE, outQ> &output)
  {
    ac_generic_float_log_pwl<pwlQ, float_base::base_2> gen_log_obj;
    gen_log_obj.interface(input, output);
  }

//=================================================================================
// Function: ac_log_pwl (for ac_float, returns log(input) )
//
// Description:
//    Calculation of log base e of real inputs, passed as ac_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_log_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_float<16, 8, 6, AC_TRN> input_type;
//    typedef ac_float<16, 8, 9, AC_TRN> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type input,
//      output_type &output
//    )
//    {
//      ac_log_pwl(input,output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 1.25;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    The ac_float version of ac_log_pwl is a friend function to the
//    ac_generic_float_log_pwl class, and depends on the interface() function in
//    said class to call the ac_fixed version of the ac_log2_pwl function and
//    calculate the final output.
//
//---------------------------------------------------------------------------------

  template <ac_q_mode pwlQ = AC_TRN, int W, int I, int E, ac_q_mode Q, int outW, int outI, int outE, ac_q_mode outQ>
  void ac_log_pwl(const ac_float<W, I, E, Q> input, ac_float<outW, outI, outE, outQ> &output)
  {
    ac_generic_float_log_pwl<pwlQ, float_base::base_e> gen_log_obj;
    gen_log_obj.interface(input, output);
  }

// For this section of the code to work, the user must include ac_std_float.h in their testbench before including the log header,
// so as to have the code import the ac_std_float and ac_ieee_float datatypes and define the __AC_STD_FLOAT_H macro.
  #ifdef __AC_STD_FLOAT_H
//=========================================================================
// Function: ac_log2_pwl (for ac_std_float, returns 1/sqrt(input) )
//
// Description:
//    Calculation of base 2 logarithm of real inputs, passed as
//    ac_std_float variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_log_pwl.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_log_pwl.h>
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
//      ac_log2_pwl(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input(0.25);
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template <ac_q_mode pwl_Q = AC_TRN, int W, int E, int outW, int outE>
  void ac_log2_pwl(
    const ac_std_float<W, E> &input,
    ac_std_float<outW, outE> &output
  )
  {
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_log2_pwl<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_std_float<outW, outE> output_temp(output_ac_fl); // Convert output ac_float to ac_std_float.
    output = output_temp;
  }

//=========================================================================
// Function: ac_log_pwl (for ac_std_float, returns 1/sqrt(input) )
//
// Description:
//    Calculation of logarithm of real inputs, passed as ac_std_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_log_pwl.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_log_pwl.h>
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
//      ac_log_pwl(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input(0.25);
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template <ac_q_mode pwl_Q = AC_TRN, int W, int E, int outW, int outE>
  void ac_log_pwl(
    const ac_std_float<W, E> &input,
    ac_std_float<outW, outE> &output
  )
  {
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_log_pwl<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_std_float<outW, outE> output_temp(output_ac_fl); // Convert output ac_float to ac_std_float.
    output = output_temp;
  }

//=================================================================================
// Function: ac_log2_pwl (for ac_ieee_float, returns log2(input) )
//
// Description:
//    Calculation of log base 2 of real inputs, passed as ac_ieee_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_log_pwl.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_log_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_ieee_float<binary32> input_type;
//    typedef ac_ieee_float<binary32> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type input,
//      output_type &output
//    )
//    {
//      ac_log2_pwl(input,output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input(1.25);
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    The ac_ieee_float version of ac_log2_pwl relies on the ac_float version to
//    perform the actual PWL computation, which in turn relies on the ac_fixed
//    implementation.
//
//---------------------------------------------------------------------------------

  template<ac_q_mode pwl_Q = AC_TRN, ac_ieee_float_format Format, ac_ieee_float_format outFormat>
  void ac_log2_pwl(const ac_ieee_float<Format> input, ac_ieee_float<outFormat> &output)
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int outW = T_out::width;
    const int outE = T_out::e_width;
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_log2_pwl<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_ieee_float<outFormat> output_temp(output_ac_fl); // Convert output ac_float to ac_ieee_float.
    output = output_temp;

    #if !defined(__SYNTHESIS__) && defined(AC_LOG_PWL_H_DEBUG)
    std::cout << "input.to_ac_float().type_name() : " << input.to_ac_float().type_name() << std::endl;
    std::cout << "output_ac_fl.type_name() : " << output_ac_fl.type_name() << std::endl;
    #endif
  }

//=================================================================================
// Function: ac_log_pwl (for ac_ieee_float, returns log(input) )
//
// Description:
//    Calculation of log base e of real inputs, passed as ac_ieee_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_log_pwl.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_log_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_ieee_float<binary32> input_type;
//    typedef ac_ieee_float<binary32> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type input,
//      output_type &output
//    )
//    {
//      ac_log_pwl(input,output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input(1.25);
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    The ac_ieee_float version of ac_log_pwl relies on the ac_float version to
//    perform the actual PWL computation, which in turn relies on the ac_fixed
//    implementation.
//
//---------------------------------------------------------------------------------

  template<ac_q_mode pwl_Q = AC_TRN, ac_ieee_float_format Format, ac_ieee_float_format outFormat>
  void ac_log_pwl(const ac_ieee_float<Format> input, ac_ieee_float<outFormat> &output)
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int outW = T_out::width;
    const int outE = T_out::e_width;
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_log_pwl<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_ieee_float<outFormat> output_temp(output_ac_fl); // Convert output ac_float to ac_ieee_float.
    output = output_temp;

    #if !defined(__SYNTHESIS__) && defined(AC_LOG_PWL_H_DEBUG)
    std::cout << "input.to_ac_float().type_name() : " << input.to_ac_float().type_name() << std::endl;
    std::cout << "output_ac_fl.type_name() : " << output_ac_fl.type_name() << std::endl;
    #endif
  }

  #endif

//=========================================================================
// Version that allows returning of values for log2.
  template<class T_out, ac_q_mode pwlQ = AC_TRN, class T_in>
  T_out ac_log2_pwl(
    const T_in input
  )
  {
    // create a variable that is to be returned
    T_out output;
    // call above implementation of log base 2
    ac_log2_pwl<pwlQ>(input, output);
    // return the final output
    return output;
  }

//=========================================================================
// Version that allows returning of values for ln.
  template<class T_out, ac_q_mode pwlQ = AC_TRN, class T_in>
  T_out ac_log_pwl(
    const T_in input
  )
  {
    // create a variable that is to be returned
    T_out output;
    // call above implementation of log base e
    ac_log_pwl<pwlQ>(input, output);
    // return the final output
    return output;
  }
}
#endif
