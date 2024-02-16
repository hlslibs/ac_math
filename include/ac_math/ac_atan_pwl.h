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
//******************************************************************************************
// Function: ac_atan_pwl
//
// Description:
//    Calculation of arctangent of real inputs.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_atan_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_fixed<16, 8, false, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<16, 8, false, AC_RND, AC_SAT> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_atan_pwl(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 0.25;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    This file uses C++ function overloading and templates for various datatype
//    implementations. Attempting to call it with a type that is not implemented will
//    result in a compile-time error.
//
//    This file uses the ac_reciprocal_pwl() function from ac_reciprocal_pwl.h as well as the
//    ac_shift_left() function from ac_shift.h.
//
// Revision History:
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//    3.3.0  - [CAT-25797] Added CDesignChecker fixes/waivers for code check violations in ac_math PWL and Linear Algebra IPs.
//             Waivers added for CNS and CCC violations.
//             Fixes added for FXD, STF and MXS violations.
//               - FXD violations fixed by changing integer literals to floating point literals or typecasting to ac_fixed values.
//               - STF violations fixed by using "const" instead of "static const" parameters. LUT generator files also print out "const" LUTs instead of "static const" LUTs.
//               - MXS violations fixed by typecasting unsigned variables to int.
//    2.0.10 - Official open-source release as part of the ac_math library.
//
//******************************************************************************************

#ifndef _INCLUDED_AC_ATAN_PWL_H_
#define _INCLUDED_AC_ATAN_PWL_H_

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

#if !defined(__SYNTHESIS__) && defined(AC_ATAN_PWL_H_DEBUG)
#include <iostream>
#endif

#include <ac_math/ac_reciprocal_pwl.h>
#include <ac_math/ac_shift.h>

//=========================================================================
// Function: ac_atan_pwl (for ac_fixed)
//
// Description:
//    Calculation of arctangent of real, positive inputs, passed as
//    ac_fixed variables.
//
// Usage:
//    See above example for usage.
//
//-------------------------------------------------------------------------

namespace ac_math
{
  template<ac_q_mode pwl_Q = AC_TRN,
           int W, int I, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  void ac_atan_pwl(
    const ac_fixed<W, I, false, Q, O> &input,
    ac_fixed<outW, outI, false, outQ, outO> &output
  )
  {
    // Store the approximate value of pi by 2
    const ac_fixed<11, 1, false> pi_by_2 = 1.5703125;

    // Give the number of fraction bits to be assigned for the normalized input. 14 is chosen for this implementation, because it results in
    // a very low error.
    const int f_b_n_i = 14;
    ac_fixed<f_b_n_i, 0, false, pwl_Q, AC_SAT> normalized_input;

    // Start of code outputted by ac_atan_pwl_lutgen.cpp
    // Note that the LUT generator file also outputs values for x_min_lut (lower limit of PWL domain), x_max_lut (upper limit of PWL domain)
    // and sc_constant_lut (scaling factor used to scale the input from 0 to n_segments_lut). However, these values aren't considered in the header
    // file because it has been optimized to work with a 4-segment PWL model that covers the domain of [0, 1). For other PWL implementations, the user will probably have
    // to take these values into account explicitly. Guidelines for doing so are given in the comments.
    // In addition, some of the slope values here are modified slightly in order to ensure monotonicity of the PWL function as the input crosses segment boundaries.
    // The user might want to take care to ensure that for their own PWL versions.

    // Initialization for PWL LUT
    const unsigned n_segments_lut = 4;
    // The number of fractional bits for the LUT values is chosen by first finding the maximum absolute error over the domain of the PWL
    // when double-precision values are used for LUT values. This error will correspond to a number of fractional bits that are always
    // guaranteed to be error-free, for fixed-point PWL outputs.
    // This number of fractional bits is found out by the formula:
    // nbits = abs(ceil(log2(abs_error_max)).
    // The number of fractional bits hereafter used to store the LUT values is nbits + 2.
    // For this particular PWL implementation, the number of fractional bits is 10.
    const int n_frac_bits = 10;
    // slope and intercept value array.
    const ac_fixed<n_frac_bits, 0, false> m_lut[n_segments_lut] = {.2451171875, .2197265625, .1796875, 0.138671875};
    const ac_fixed<n_frac_bits, 0, false> c_lut[n_segments_lut] = {.0009765625, .24609375, .4658203125, .646484375};

    // End of code outputted by ac_atan_pwl_lutgen.cpp

    // If the input exceeds or equals 1, we take the reciprocal of the input and find the arctangent of that reciprocal. We then use the formula
    // atan(x) = pi/2 - atan(1 / x) to find the arctangent of the original input.
    // Also, keep in mind that the input can only exceed 1 if the number of integer bits are greater than or equal to 1 and we won't need a reciprocal operation if that's
    // not the case. The "if (I >= 1)" condition will then ensure that the reciprocal block is optimized away.
    bool input_exceeds_1 = false;
    #pragma hls_waive CNS
    if (I >= 1) { input_exceeds_1 = (input > 1); }
    if ((I >= 1) && input_exceeds_1) { ac_math::ac_reciprocal_pwl<pwl_Q>(input, normalized_input); }
    // If input is lesser than 1, then it is within the domain of the PWL function. Hence, no reciprocal operation is required.
    else {normalized_input = input;}

    // Compute atan using pwl.
    const int int_bits = ac::nbits<n_segments_lut - 1>::val;
    // Scale the input from 0 to 4. Note that the below expression is simplified and optimized for 4 segments and a domain of [0, 1). Any other PWL implementation
    // with a different number of segments/domain should be scaled according to the formula: x_in_sc = (normalized_input - x_min_lut) * sc_constant_lut
    // where sc_constant_lut = n_segments_lut / (x_max_lut - x_min_lut)
    // where x_min_lut and x_max_lut are the lower and upper limits of the PWL Domain.
    ac_fixed<f_b_n_i, int_bits, false> x_in_sc = (ac_fixed<f_b_n_i + int_bits, int_bits, false>)normalized_input << 2;
    // Take out the fractional bits of the scaled input
    ac_fixed<f_b_n_i - int_bits, 0, false> x_in_sc_frac;
    x_in_sc_frac.set_slc(0, x_in_sc.template slc<f_b_n_i - int_bits>(0));
    ac_int<int_bits, false> index;
    // The integer part of the scaled input is the index of the LUT table
    index = x_in_sc.to_int();
    // The precision given below will ensure that there is no precision lost in the assignment to output_pwl, hence rounding for the variable is switched off by default.
    // However, if the user uses less fractional bits and wishes to turn rounding on instead, they are welcome to do so by giving a different value for pwl_Q.
    ac_fixed<f_b_n_i - int_bits + n_frac_bits + 1, 1, false, pwl_Q> output_pwl = m_lut[index]*x_in_sc_frac + c_lut[index];

    // If the input exceeds 1, apply the previously mentioned formula.
    if ((I >= 1) && input_exceeds_1) { output_pwl = pi_by_2 - output_pwl; }

    output = output_pwl;

    #if !defined(__SYNTHESIS__) && defined(AC_ATAN_PWL_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "input            = " << input << std::endl;
    std::cout << "normalized_input = " << normalized_input << std::endl;
    std::cout << "x_in_sc          = " << x_in_sc << std::endl;
    std::cout << "output_pwl       = " << output_pwl << std::endl;
    std::cout << "output           = " << output << std::endl;
    #endif
  }

//=========================================================================
// Function: ac_atan_pwl (for ac_float)
//
// Description:
//    Calculation of arctangent of real, positive inputs, passed as
//    ac_float variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_atan_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_float<16, 8, 8, AC_TRN> input_type;
//    typedef ac_float<16, 8, 8, AC_TRN> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_atan_pwl(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 0.25;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template<ac_q_mode pwl_Q = AC_TRN,
           int W, int I, int E, ac_q_mode Q,
           int outW, int outI, int outE, ac_q_mode outQ>
  void ac_atan_pwl(
    const ac_float<W, I, E, Q> &input,
    ac_float<outW, outI, outE, outQ> &output
  )
  {
    #ifdef ASSERT_ON_INVALID_INPUT
    // ac_atan_pwl only works for first quadrant angles and hence can only accept positive
    // inputs. The AC_ASSERT below will ensure that that is always the case.
    AC_ASSERT(input.mantissa() >= 0, "Input must be positive.");
    #endif
    const ac_fixed<11, 1, false> pi_by_2 = 1.5703125; // Store the approximate value of pi by 2

    // Calculate intermediate bitwidth of reciprocal output. Consider changing the formula if PWL implementation changes from the default.
    const int recipW = AC_MAX(1, AC_MIN(10, W - 4)) + 10 + 2;

    ac_float<recipW, I, E> atan_input_fl; // Temp variable to store the reciprocal of the ac_float input, if necessary.

    // If input exceeds one, pass the reciprocal of the input to the atan function and use
    // the formula atan(x) = pi/2 - atan(1/x) to calculate the final output.
    bool input_exceeds_1 = (input >= 1);
    #pragma hls_waive CNS
    if (input_exceeds_1) { ac_math::ac_reciprocal_pwl<pwl_Q>(input, atan_input_fl); }
    else { atan_input_fl = input; }

    // Since inputs are always supposed to be positive, we don't need the signed bit of the
    // mantissa as that will always be zero. Reduce the bitwidth of mantVal accordingly.
    ac_fixed<recipW - 1, I - 1, false> mantVal = atan_input_fl.mantissa();
    int exp_val = atan_input_fl.exp().to_int();
    // Ultimately, we will have to convert the input to an ac_fixed value and pass it to the ac_fixed
    // implementation of ac_atan_pwl. The purpose of performing the reciprocal operation before passing
    // the input to the ac_fixed implementation is to restrict the input domain to [0, 1) and hence restrict the
    // bitwidth of the ac_fixed input. We need zero integer bits (since the ac_fixed input is always fractional),
    // which will enable further optimizations in the ac_fixed version code.
    ac_fixed<14, 0, false, AC_RND, AC_SAT> atan_input_fi;
    // atan_input_fi = mantVal*2^(exp_val) = mantVal << exp_val
    // Use ac_shift_left, as well as output rounding/saturation, to prevent overflow and ensure rounding.
    ac_math::ac_shift_left(mantVal, exp_val, atan_input_fi);
    ac_fixed<23, 1, false> atan_output_fi;
    ac_atan_pwl<pwl_Q>(atan_input_fi, atan_output_fi); // Call ac_fixed version.
    if (input_exceeds_1) { atan_output_fi = pi_by_2 - atan_output_fi; } // atan(x) = pi/2 - atan(1/x)
    // Convert ac_fixed output to ac_float by using a constructor.
    ac_float<outW, outI, outE, outQ> output_temp(atan_output_fi);
    output = output_temp;
  }

// For this section of the code to work, the user must include ac_std_float.h in their testbench before including the arctangent header,
// so as to have the code import the ac_std_float and ac_ieee_float datatypes and define the __AC_STD_FLOAT_H macro.
  #ifdef __AC_STD_FLOAT_H
//=========================================================================
// Function: ac_atan_pwl (for ac_std_float)
//
// Description:
//    Calculation of arctangent of real inputs, passed as ac_std_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_atan_pwl.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_atan_pwl.h>
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
//      ac_atan_pwl(input, output);
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
  void ac_atan_pwl(
    const ac_std_float<W, E> &input,
    ac_std_float<outW, outE> &output
  )
  {
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_atan_pwl<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_std_float<outW, outE> output_temp(output_ac_fl); // Convert output ac_float to ac_std_float.
    output = output_temp;
  }

//=========================================================================
// Function: ac_atan_pwl (for ac_ieee_float)
//
// Description:
//    Calculation of arctangent of real, positive inputs, passed as
//    ac_ieee_float variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_atan_pwl.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_atan_pwl.h>
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
//      ac_atan_pwl(input, output);
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

  template<ac_q_mode pwl_Q = AC_TRN,
           ac_ieee_float_format Format,
           ac_ieee_float_format outFormat>
  void ac_atan_pwl(
    const ac_ieee_float<Format> &input,
    ac_ieee_float<outFormat> &output
  )
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int outW = T_out::width;
    const int outE = T_out::e_width;
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_atan_pwl<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_ieee_float<outFormat> output_temp(output_ac_fl); // Convert output ac_float to ac_ieee_float.
    output = output_temp;
  }
  #endif

  // The following version enables a return-by-value.
  template<class T_out,
           ac_q_mode pwl_Q = AC_TRN,
           class T_in>
  T_out ac_atan_pwl(const T_in &input)
  {
    T_out output;
    ac_atan_pwl<pwl_Q>(input, output);
    return output;
  }

}

#endif // _INCLUDED_AC_ATAN_PWL_H_
