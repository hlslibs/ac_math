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
// Function: ac_atan_pwl_ha
//
// Description:
//    High-accuracy tangent calculation of real inputs.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_atan_pwl_ha.h>
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
//      ac_atan_pwl_ha(input, output);
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
//    This file uses the ac_reciprocal_pwl_ha() function from ac_reciprocal_pwl_ha.h as
//    well as the ac_shift_left() function from ac_shift.h.
//
// Revision History:
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//    3.3.0  - [CAT-25798] Added CDesignChecker fixes/waivers for code check and Synthesis-simulation mismatch/violations in ac_math PWL and Linear Algebra IPs.
//    3.2.4 - Improved bitwidth calculations for floating point datatypes.
//    3.2.3 - Created header file.
//
//******************************************************************************************

#ifndef _INCLUDED_AC_ATAN_PWL_HA_H_
#define _INCLUDED_AC_ATAN_PWL_HA_H_

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

#if !defined(__SYNTHESIS__) && defined(AC_ATAN_PWL_HA_H_DEBUG)
#include <iostream>
#endif

#include <ac_math/ac_reciprocal_pwl_ha.h>
#include <ac_math/ac_shift.h>

//=========================================================================
// Function: ac_atan_pwl_ha (for ac_fixed)
//
// Description:
//    High accuracy calculation of arctangent of real, positive inputs,
//    passed as ac_fixed variables.
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
  void ac_atan_pwl_ha(
    const ac_fixed<W, I, false, Q, O> &input,
    ac_fixed<outW, outI, false, outQ, outO> &output
  )
  {
    // Number of fraction bits to be assigned for the normalized input.
    const int f_b_n_i = 32;
    ac_fixed<f_b_n_i, 0, false, pwl_Q, AC_SAT> normalized_input;

    // Start of code outputted by ac_atan_pwl_ha_lutgen.cpp

    // Initialization for PWL LUT.
    const unsigned n_segments_lut = 32;
    const int n_frac_bits = 16;
    const ac_fixed<32, 1, false> pi_by_2 = 1.570796326734125614166259765625;
    const ac_fixed<-4 + n_frac_bits, -4, false> m_lut[n_segments_lut] = {.0312347412109375, .03118896484375, .03106689453125, .0308837890625, .0306396484375, .0303497314453125, .0300140380859375, .0296173095703125, .0291900634765625, .028717041015625, .0282135009765625, .027679443359375, .0271148681640625, .026519775390625, .0259246826171875, .0253143310546875, .024688720703125, .0240631103515625, .0234222412109375, .0227813720703125, .02215576171875, .0215301513671875, .020904541015625, .020294189453125, .0196990966796875, .0191192626953125, .0185394287109375, .017974853515625, .017425537109375, .0168914794921875, .0163726806640625, .0158233642578125};
    const ac_fixed<1 + n_frac_bits, 1, false> c_lut[n_segments_lut] = {0.0, .0312347412109375, .0624237060546875, .0934906005859375, .1243743896484375, .1550140380859375, .18536376953125, .2153778076171875, .2450103759765625, .274200439453125, .30291748046875, .3311309814453125, .3588104248046875, .38592529296875, .412445068359375, .4383697509765625, .46368408203125, .488372802734375, .5124359130859375, .535858154296875, .5586395263671875, .5807952880859375, .602325439453125, .6232452392578125, .6435394287109375, .663238525390625, .6823577880859375, .700897216796875, .7188720703125, .736297607421875, .7531890869140625, .769561767578125};
    const ac_fixed<1 + n_frac_bits, 1, false> x_min_lut = .0;
    const ac_fixed<1 + n_frac_bits, 1, false> x_max_lut = 1.0;
    const ac_fixed<6 + n_frac_bits, 6, false> sc_constant_lut = 32.0;
    // End of code outputted by ac_atan_pwl_ha_lutgen.cpp

    // If the input exceeds or equals 1, we take the reciprocal of the input and find the arctangent of that reciprocal. We then use the formula
    // atan(x) = pi/2 - atan(1 / x) to find the arctangent of the original input.
    // Also, keep in mind that the input can only exceed 1 if the number of integer bits are greater than or equal to 1 and we won't need a reciprocal operation if that's
    // not the case. The "if (I >= 1)" condition will then ensure that the reciprocal block is optimized away.
    bool input_exceeds_1 = false;
    if (I >= 1) { input_exceeds_1 = input > 1 ? true : false; }
    if ((I >= 1) && input_exceeds_1) { ac_math::ac_reciprocal_pwl_ha<pwl_Q>(input, normalized_input); }
    // If input is lesser than 1, then it is within the domain of the PWL function. Hence, no reciprocal operation is required.
    else { normalized_input = input; }

    // Compute atan using pwl.
    const int int_bits = ac::nbits<n_segments_lut - 1>::val;
    ac_fixed<f_b_n_i, int_bits, false> x_in_sc = (ac_fixed<f_b_n_i + int_bits, int_bits, false>)(normalized_input - x_min_lut)*sc_constant_lut;
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

    #if !defined(__SYNTHESIS__) && defined(AC_ATAN_PWL_HA_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "input            = " << input << std::endl;
    std::cout << "normalized_input = " << normalized_input << std::endl;
    std::cout << "x_in_sc          = " << x_in_sc << std::endl;
    std::cout << "output_pwl       = " << output_pwl << std::endl;
    std::cout << "output           = " << output << std::endl;
    #endif
  }

//=========================================================================
// Function: ac_atan_pwl_ha (for ac_float)
//
// Description:
//    High-accuracy alculation of arctangent of real, positive inputs,
//    passed as ac_float variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_atan_pwl_ha.h>
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
//      ac_atan_pwl_ha(input, output);
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
  void ac_atan_pwl_ha(
    const ac_float<W, I, E, Q> &input,
    ac_float<outW, outI, outE, outQ> &output
  )
  {
    #ifdef ASSERT_ON_INVALID_INPUT
    // ac_atan_pwl_ha only works for first quadrant angles and hence can only accept positive
    // inputs. The AC_ASSERT below will ensure that that is always the case.
    AC_ASSERT(input.mantissa() >= 0, "Input must be positive.");
    #endif
    const ac_fixed<32, 1, false> pi_by_2 = 1.570796326734125614166259765625;

    // Calculate intermediate bitwidth of reciprocal output.
    const int n_frac_bits_fixed_recip = 16; // Enter the number of fractional bits specified by the "n_frac_bits" variable in the fixed point reciprocal PWL.
    // The temporary floating point variable will either store the reciprocal of the input or the input itself. Make sure it has enough bits to store either.
    const int recipW = AC_MAX(2*n_frac_bits_fixed_recip + 2, W);
    // Make sure that the temporary floating point variable has enough bits in the exponent part to account for the minimum possible value that might be stored.
    const int min_exp_val = AC_MAX(2*I - 3, 1) + (1 << (E - 1)) - 1;
    const int recipE = ac::nbits<min_exp_val>::val + 1;

    ac_float<recipW, I, recipE> atan_input_fl; // Temp variable to store the reciprocal of the ac_float input, if necessary.

    // If input exceeds one, pass the reciprocal of the input to the atan function and use
    // the formula atan(x) = pi/2 - atan(1/x) to calculate the final output.
    bool input_exceeds_1 = (input >= 1);
    if (input_exceeds_1) { ac_math::ac_reciprocal_pwl_ha<pwl_Q>(input, atan_input_fl); }
    else { atan_input_fl = input; }

    // Since inputs are always supposed to be positive, we don't need the signed bit of the
    // mantissa as that will always be zero. Reduce the bitwidth of mantVal accordingly.
    ac_fixed<recipW - 1, I - 1, false> mantVal = atan_input_fl.mantissa();
    int exp_val = atan_input_fl.exp().to_int();
    // Ultimately, we will have to convert the input to an ac_fixed value and pass it to the ac_fixed
    // implementation of ac_atan_pwl_ha. The purpose of performing the reciprocal operation before passing
    // the input to the ac_fixed implementation is to restrict the input domain to [0, 1) and hence restrict the
    // bitwidth of the ac_fixed input. We need zero integer bits (since the ac_fixed input is always fractional),
    // which will enable further optimizations in the ac_fixed version code.
    ac_fixed<32, 0, false, AC_RND, AC_SAT> atan_input_fi;
    // atan_input_fi = mantVal*2^(exp_val) = mantVal << exp_val
    // Use ac_shift_left, as well as output rounding/saturation, to prevent overflow and ensure rounding.
    ac_math::ac_shift_left(mantVal, exp_val, atan_input_fi);
    ac_fixed<43, 1, false> atan_output_fi;
    ac_atan_pwl_ha<pwl_Q>(atan_input_fi, atan_output_fi); // Call ac_fixed version.
    if (input_exceeds_1) { atan_output_fi = pi_by_2 - atan_output_fi; } // atan(x) = pi/2 - atan(1/x)
    // Convert ac_fixed output to ac_float by using a constructor.
    ac_float<outW, outI, outE, outQ> output_temp(atan_output_fi);
    output = output_temp;
  }

// For this section of the code to work, the user must include ac_std_float.h in their testbench before including the arctangent header,
// so as to have the code import the ac_std_float and ac_ieee_float datatypes and define the __AC_STD_FLOAT_H macro.
  #ifdef __AC_STD_FLOAT_H
//=========================================================================
// Function: ac_atan_pwl_ha (for ac_std_float)
//
// Description:
//    Calculation of arctangent of real inputs, passed as ac_std_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_atan_pwl_ha.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_atan_pwl_ha.h>
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
//      ac_atan_pwl_ha(input, output);
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
  void ac_atan_pwl_ha(
    const ac_std_float<W, E> &input,
    ac_std_float<outW, outE> &output
  )
  {
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_atan_pwl_ha<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_std_float<outW, outE> output_temp(output_ac_fl); // Convert output ac_float to ac_std_float.
    output = output_temp;
  }

//=========================================================================
// Function: ac_atan_pwl_ha (for ac_ieee_float)
//
// Description:
//    Calculation of arctangent of real, positive inputs, passed as
//    ac_ieee_float variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_atan_pwl_ha.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_atan_pwl_ha.h>
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
//      ac_atan_pwl_ha(input, output);
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
  void ac_atan_pwl_ha(
    const ac_ieee_float<Format> &input,
    ac_ieee_float<outFormat> &output
  )
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int outW = T_out::width;
    const int outE = T_out::e_width;
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_atan_pwl_ha<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_ieee_float<outFormat> output_temp(output_ac_fl); // Convert output ac_float to ac_ieee_float.
    output = output_temp;
  }
  #endif

  // The following version enables a return-by-value.
  template<class T_out,
           ac_q_mode pwl_Q = AC_TRN,
           class T_in>
  T_out ac_atan_pwl_ha(const T_in &input)
  {
    T_out output;
    ac_atan_pwl_ha<pwl_Q>(input, output);
    return output;
  }

}

#endif // _INCLUDED_AC_ATAN_PWL_HA_H_
