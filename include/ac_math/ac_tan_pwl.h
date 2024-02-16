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
//**************************************************************************************************
// File: ac_tan_pwl.h
//
// Description: Provides piece-wise linear implementations of the
//    tangent function for the AC (tm) Datatypes: ac_fixed, ac_float,
//    ac_ieee_float.
//
//    The inputs to each function in this header must be positive and
//    must lie in the domain of [0, pi/2)
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_tan_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_fixed<18, 2, false, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<64, 32, false, AC_RND, AC_SAT> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_tan_pwl(input,output);
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
//    This file uses the ac_reciprocal_pwl() function from ac_reciprocal_pwl.h as well as
//    the ac_shift_left() function from ac_shift.h.
//
// Revision History:
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//    3.3.0  - CAT-25797, waived CNS violations and fixed FXD violations reported by CDesignChecker.
//    3.2.4  - CAT-24698, Added floating point support.
//    3.1.2  - Improved bitwidth calculations for PWL.
//    2.0.10 - Official open-source release as part of the ac_math library.
//
//**************************************************************************************************

#ifndef _INCLUDED_AC_TAN_PWL_H_
#define _INCLUDED_AC_TAN_PWL_H_

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

#if !defined(__SYNTHESIS__) && defined(AC_TAN_PWL_H_DEBUG)
#include <iostream>
#endif

#include <ac_math/ac_reciprocal_pwl.h>
#include <ac_math/ac_shift.h>

//******************************************************************************************
// Function: ac_tan_pwl (for ac_fixed)
//
// Description:
//    Calculation of tangent of real inputs, passed as ac_fixed
//    variables. Inputs and outputs must both be unsigned.
//
// Usage:
//    See above example code for usage.
//
// Notes:
//    The PWL implementation utilizes 8 segments, which has a small impact
//    on accuracy.
//    The ac_fixed implementation also serves as the PWL-based "backbone"
//    for the ac_float and ac_ieee_float versions.
//
//******************************************************************************************

namespace ac_math
{
  template<ac_q_mode pwl_Q = AC_TRN,
           int W, int I, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  void ac_tan_pwl(
    const ac_fixed<W, I, false, Q, O> &input,
    ac_fixed<outW, outI, false, outQ, outO> &output
  )
  {
    #ifdef ASSERT_ON_INVALID_INPUT
    AC_ASSERT(input < 1.5704021453857421875, "Input must not exceed pi/2");
    #endif

    // Store the approximate value of 89.17 degrees in radian
    const ac_fixed<11, 1, false> sat_limit = 1.556640625;
    ac_fixed<outW, outI, false, outQ, outO> output_temp;

    // Start of code outputted by ac_tan_pwl_lutgen.cpp
    // NOTE: Any change from the default PWL implementation using values generated by the lutgen file is not guaranteed to give a monotonically increasing PWL output,
    // owing to errors incurred in quantizing the double slope/intercept values to fixed point. The user must change the slope-intercept values themself to ensure monotonicity.

    // The number of fractional bits for the LUT values is chosen by first finding the maximum absolute error over the domain of the PWL
    // when double-precision values are used for LUT values. This error will correspond to a number of fractional bits that are always
    // guaranteed to be error-free, for fixed-point PWL outputs.
    // This number of fractional bits is found out by the formula:
    // nbits = abs(ceil(log2(abs_error_max)) - 1
    // The number of fractional bits hereafter used to store the LUT values is nbits + 2.
    // For this particular PWL implementation, the number of fractional bits is 10.
    const int n_frac_bits = 10;
    // Initialization for PWL LUT
    const unsigned n_segments_lut = 8;
    const ac_fixed<n_frac_bits, 0, false> m_lut[n_segments_lut] = {.0986328125, .099609375, .1044921875, .1103515625, .1201171875, .1328125, .15234375, .1787109375};
    const ac_fixed<n_frac_bits, 0, false> c_lut[n_segments_lut] = {.0, .0986328125, .1982421875, .3037109375, .4140625, .5341796875, .6669921875, .8193359375};
    // Domain of PWL
    const ac_fixed<n_frac_bits, 0, false> x_max_lut = .78515625;
    const ac_fixed<1, 1, false> x_min_lut = 0.0;
    const ac_fixed<14, 4, false> sc_constant_lut = 10.1884765625;

    // End of code outputted by ac_tan_pwl_lutgen.cpp

    // If the input equals or exceeds pi/4, we halve it and use the formula tan(2*x) = 2*tan(x) / (1 - tan(x)^2) to get the tan value.
    // You will only need to add an extra fractional bit if the number of integer bits is greater than or equal to zero, because only then
    // will the input have a chance of exceeding pi/4, and only then will halving be required.

    // The intermediate variable for the input can have a maximum of 20 fractional bits. This maximum is chosen empirically, to limit the
    // area used while also making sure that the error values are the same, up to 3 decimal places, as they would be if we were to not use the
    // 20-bit maximum.
    const int n_f_b_int = AC_MIN(19, W - I) + int(I >= 0);
    const int I_int = AC_MIN(I, 1);
    const int W_int = I_int + n_f_b_int;
    ac_fixed<W_int, I_int, false, Q, O> input_int;
    bool input_exceeds_pi_by_4;
    // Keep in mind that the input can only exceed pi/4 if the number of integer bits is greater than or equal to zero.
    #pragma hls_waive CNS
    if (I >= 0) { input_exceeds_pi_by_4 = (input >= x_max_lut) ? true : false; }

    #pragma hls_waive CNS
    if ((I >= 0) && input_exceeds_pi_by_4) { input_int = (ac_fixed<W_int + 1, I_int, false, Q, O>)input >> 1; }
    else { input_int = input; }

    const int int_bits = ac::nbits<n_segments_lut - 1>::val;
    // Compute tan using pwl.
    // Scale the input from 0 to n_segments_lut
    ac_fixed<n_frac_bits + int_bits, int_bits, false> x_in_sc = (input_int - x_min_lut) * sc_constant_lut;
    // Take out the fractional bits of the scaled input
    ac_fixed<n_frac_bits, 0, false> x_in_sc_frac;
    x_in_sc_frac.set_slc(0, x_in_sc.template slc<n_frac_bits>(0));
    ac_int<int_bits, false> index;
    // The integer part of the input is the index of the LUT table
    index = x_in_sc.to_int();
    // The precision given below will ensure that there is no precision lost in the assignment to output_pwl, hence rounding for the variable is switched off by default.
    // However, if the user wishes to use less fractional bits and turn rounding on instead, they are welcome to do so by giving a different value for pwl_Q.
    typedef ac_fixed<2*n_frac_bits, 0, false, pwl_Q> output_pwl_type;
    output_pwl_type output_pwl = m_lut[index]*x_in_sc_frac + c_lut[index];

    // As mentioned earlier, if the input equals or exceeds pi/4, we use the formula tan(2*x) = 2*tan(x) / (1 - tan(x)^2) to get the tan value.
    #pragma hls_waive CNS
    if ((I >= 0) && input_exceeds_pi_by_4) {
      // The bitwidths for the following two declarations were chosen after careful, exhaustive testing to make sure that there were negligible error penalties while using
      // the lowest possible bitwidths, so as to limit the area used, preserve monotonicity and optimize QofR. These bitwidths are optimized for a domain of [0, pi/4) and
      // 8 segments. The user may have to change them in order to suit their own design, in case they change the number of segments and/or the domain.
      ac_fixed<16, 0, false, pwl_Q> one_minus_tan_theta_by_2_sqr = 1 - output_pwl * output_pwl;
      ac_fixed<16, 6, false, pwl_Q> recip_value;
      // Use the reciprocal_pwl function to calculate 1 / (1 - tan(x)^2)
      ac_math::ac_reciprocal_pwl<pwl_Q>(one_minus_tan_theta_by_2_sqr, recip_value);
      output_temp = 2 * output_pwl * recip_value;
    } else {
      output_temp = output_pwl;
    }

    // If input crosses or equals 89.17 degrees (roughly), set the output to saturate
    #pragma hls_waive CNS
    if (input >= sat_limit) { output_temp.template set_val<AC_VAL_MAX>(); }

    output = output_temp;

    #if !defined(__SYNTHESIS__) && defined(AC_TAN_PWL_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "input           = " << input << std::endl;
    std::cout << "sc_constant_lut = " << sc_constant_lut << std::endl;
    std::cout << "input_int       = " << input_int << std::endl;
    std::cout << "x_in_sc         = " << x_in_sc << std::endl;
    std::cout << "output_pwl      = " << output_pwl << std::endl;
    std::cout << "output_temp     = " << output_temp << std::endl;
    std::cout << "output          = " << output << std::endl;
    #endif
  }

//******************************************************************************************
// Function: ac_tan_pwl (for ac_float)
//
// Description:
//    Calculation of tangent of real inputs, passed as ac_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_tan_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_float<10, 1, 6, AC_TRN> input_type;
//    typedef ac_float<16, 1, 9, AC_TRN> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_tan_pwl(input,output);
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
//    This implementation relies on the ac_fixed implementation of ac_tan_pwl.
//
//******************************************************************************************

  template<ac_q_mode pwl_Q = AC_TRN,
           int W, int I, int E, ac_q_mode Q,
           int outW, int outI, int outE, ac_q_mode outQ>
  void ac_tan_pwl(
    const ac_float<W, I, E, Q> &input,
    ac_float<outW, outI, outE, outQ> &output
  )
  {
    #ifdef ASSERT_ON_INVALID_INPUT
    AC_ASSERT(input >= 0, "Input must be positive.");
    AC_ASSERT(input < 1.57040202617645263671875, "Input must not exceed pi/2");
    #endif

    // Input is always assumed to be positive -> sign bit is unnecessary.
    ac_fixed<W - 1, I - 1, false> mantVal = input.mantissa();
    int exp_val = input.exp().to_int();
    // Intermediate ac_fixed variables to store the value of inputs and outputs
    // and enable compatibility with ac_fixed implementation.
    ac_fixed<22, 1, false, AC_RND> tan_input_fi;
    // tan_input_fi = mantVal*2^(exp_val) = mantVal << exp_val
    // Use ac_shift_left instead of "<<" operator to ensure rounding.
    // Saturation is not required because the input is assumed to never exceed pi/2.
    ac_math::ac_shift_left(mantVal, exp_val, tan_input_fi);
    ac_fixed<36, 7, false> tan_output_fi;
    ac_tan_pwl<pwl_Q>(tan_input_fi, tan_output_fi);
    // Convert ac_fixed output to ac_float by using a constructor.
    ac_float<outW, outI, outE, outQ> output_temp(tan_output_fi);

    output = output_temp;
  }

// For this section of the code to work, the user must include ac_std_float.h in their testbench before including the tan header,
// so as to have the code import the ac_std_float and ac_ieee_float datatypes, and define the __AC_STD_FLOAT_H macro.
  #ifdef __AC_STD_FLOAT_H
//=========================================================================
// Function: ac_tan_pwl (for ac_std_float)
//
// Description:
//    Calculation of tangent of real inputs, passed as ac_std_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_tan_pwl.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_tan_pwl.h>
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
//      ac_tan_pwl(input, output);
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
  void ac_tan_pwl(
    const ac_std_float<W, E> &input,
    ac_std_float<outW, outE> &output
  )
  {
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_tan_pwl<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_std_float<outW, outE> output_temp(output_ac_fl); // Convert output ac_float to ac_std_float.
    output = output_temp;
  }

//=========================================================================
// Function: ac_tan_pwl (for ac_ieee_float)
//
// Description:
//    Calculation of tangent of real inputs, passed as ac_ieee_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_tan_pwl.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_tan_pwl.h>
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
//      ac_tan_pwl(input, output);
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
  void ac_tan_pwl(
    const ac_ieee_float<Format> &input,
    ac_ieee_float<outFormat> &output
  )
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int outW = T_out::width;
    const int outE = T_out::e_width;
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_tan_pwl<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_ieee_float<outFormat> output_temp(output_ac_fl); // Convert output ac_float to ac_ieee_float.
    output = output_temp;
  }
  #endif

  // The following version enables a return-by-value.
  template<class T_out, ac_q_mode pwl_Q = AC_TRN, class T_in>
  T_out ac_tan_pwl(const T_in &input)
  {
    T_out output;
    ac_tan_pwl(input, output);
    return output;
  }

}

#endif // _INCLUDED_AC_TAN_PWL_H_
