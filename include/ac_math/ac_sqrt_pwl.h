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
// *****************************************************************************************
// File : ac_sqrt_pwl.h
//
// Created on: Jun 14, 2017
//
// Author: Sachchidanand Deo
//
// Description: Provides piece-wise linear implementations of the
//  square root function for the AC (tm) Datatypes: ac_fixed, ac_float,
//  ac_complex<ac_fixed> and ac_complex<ac_float>.
//  The ac_fixed version must be provided with unsigned input/output types,
//  while the ac_complex<ac_fixed> version must have signed real/imaginary parts.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_sqrt_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_fixed<16, 8, false, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<16, 8, false, AC_RND, AC_SAT> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type input,
//      output_type &output
//    )
//    {
//      ac_sqrt_pwl(input,output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
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
//    This file uses C++ function overloading and templates for various datatype
//    implementations. Attempting to call it with a type that is not implemented will
//    result in a compile-time error.
//
//    This file uses the ac_normalize() function from ac_normalize.h and the ac_shift_left()
//    function from ac_shift.h
//
// Revision History:
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//    3.3.0  - [CAT-25797] Added CDesignChecker fixes/waivers for code check violations in ac_math PWL and Linear Algebra IPs.
//             Waivers added for CNS and CCC violations.
//             Fixes added for FXD, STF and MXS violations.
//               - FXD violations fixed by changing integer literals to floating point literals or typecasting to ac_fixed values.
//               - STF violations fixed by using "const" instead of "static const" parameters. LUT generator files also print out "const" LUTs instead of "static const" LUTs.
//               - MXS violations fixed by typecasting unsigned variables to int.
//    3.1.2  - Improved bitwidth calculations for PWL. Removed direct access to ac_float
//             data members. Fixed bug in output near normalized 1.
//    3.1.0  - bug51145 - Improved ac_float outputting.
//    2.0.10 - Official open-source release as part of the ac_math library.
//
// *****************************************************************************************

#ifndef _INCLUDED_AC_SQRT_PWL_H_
#define _INCLUDED_AC_SQRT_PWL_H_

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
#include <ac_complex.h>

// Include headers for required functions
#include <ac_math/ac_shift.h>
#include <ac_math/ac_normalize.h>

#if !defined(__SYNTHESIS__) && defined(AC_SQRT_PWL_H_DEBUG)
#include <iostream>
#endif

//=========================================================================
// Function: ac_sqrt_pwl (for ac_fixed)
//
// Description:
//    Calculation of square root of positive real inputs, passed as ac_fixed
//    variables.
//
//    Passes inputs for normalization to the function in ac_normalize.h,
//    which gives the exponent and output normalized between 0.5 and 1.
//    The normalized value is then subject to the piecewise linear
//    implementation to calculate the square root.
//
//    For simplification square root of exponent is simply computed by
//    dividing the exponent value by 2 and square root of normalized value is
//    computed using piecewise linear mechanism.
//
// Usage:
//    Please check code snippet from above for usage.
//
//-------------------------------------------------------------------------

namespace ac_math
{
  // Only unsigned inputs/outputs are accepted.
  template <ac_q_mode pwlQ = AC_TRN, int W, int I, ac_q_mode Q, ac_o_mode O, int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  void ac_sqrt_pwl(const ac_fixed <W, I, false, Q, O> input, ac_fixed <outW, outI, false, outQ, outO> &output, const bool call_normalize = true)
  {
    // Temporary variable to store output
    ac_fixed <outW, outI, false, outQ, outO> output_temp;

    // Declaring square root of 2 as constant
    const ac_fixed <13, 1, false> root2 = 1.414306640625;

    // normalized_input is basically output of the normalization function
    ac_fixed <W, 0, false, Q, O> normalized_input;
    int normalized_exp;

    // If call_normalize is set to true, the design does not assume that the input is already normalized and performs normalization by calling ac_normalize().
    // If call_normalize is set to false, the design assumes that the input is already normalized and does not call ac_normalize, thereby saving on hardware.
    #pragma hls_waive CNS
    if (call_normalize) {
      normalized_exp = ac_math::ac_normalize(input, normalized_input);
    } else {
      normalized_exp = I;
      normalized_input.set_slc(0, input.template slc<W>(0));
    }

    // Start of code outputted by ac_sqrt_pwl_lutgen.cpp

    const unsigned n_segments_lut = 4; // Number of PWL segments.
    const int n_frac_bits = 12; // Number of fractional bits
    // Since scaling constant is a positive power-of-two, multiplication with it is the same as left-shifting by 3.
    // Accordingly, the scaled normalized input will have 3 less fractional bits than the normalized input, provided that this
    // number of fractional bits is lesser than n_frac_bits. If not, the number of fractional bits in the scaled input is set to n_frac_bits.
    const int sc_input_frac_bits = AC_MAX(1, AC_MIN(n_frac_bits, W - 3));
    // Slope and intercept LUT values.
    const ac_fixed<-2 + n_frac_bits, -2, false> m_lut[n_segments_lut] = {.08349609375, .0751953125, .0693359375, .064453125};
    const ac_fixed<1 + n_frac_bits, 1, false> c_lut[n_segments_lut] = {.707763671875, .791259765625, .866455078125, .935791015625};
    const ac_fixed<0 + n_frac_bits, 0, false> x_min_lut = .5; // Minimum limit of PWL domain
    const ac_fixed<4 + n_frac_bits, 4, false> sc_constant_lut = 8.0; // Scaling constant

    // End of code outputted by ac_sqrt_pwl_lutgen.cpp

    const int int_bits = ac::nbits<n_segments_lut - 1>::val;
    // Scale input to the range [0, n_segments_lut)
    ac_fixed<int_bits + sc_input_frac_bits, int_bits, false> input_sc = (normalized_input - x_min_lut)*sc_constant_lut;
    // Take out the fractional bits of the scaled input
    ac_fixed<sc_input_frac_bits, 0, false> input_sc_frac;
    input_sc_frac.set_slc(0, input_sc.template slc<sc_input_frac_bits>(0));
    // index is taken as integer part of scaled value and used for selection of m and c values
    ac_int <int_bits, false> index = input_sc.to_int();

    // normalized output provides square root of normalized value
    ac_fixed <sc_input_frac_bits + n_frac_bits + 1, 1, false, pwlQ> normalized_output = m_lut[index]*input_sc_frac + c_lut[index];

    // Handling of odd exponents
    ac_fixed <sc_input_frac_bits + n_frac_bits + 1, 1, false, pwlQ> normalized_output_temp = normalized_output * root2;
    // The precision given below will ensure that there is no precision lost in the assignment to m1, hence rounding for the variable is switched off by default.
    // However, if the user uses less fractional bits and turn rounding on instead, they are welcome to do so by giving a different value for pwlQ.
    #pragma hls_waive CNS
    ac_fixed <sc_input_frac_bits + n_frac_bits + 1, 1, false, pwlQ> m1 = (normalized_exp % 2 == 0) ? normalized_output : normalized_output_temp;

    // exponent and normalized output are combined to get the final ac_fixed value, which is written at memory location of output
    ac_math::ac_shift_left(m1, normalized_exp >> 1, output_temp);
    output = (input == 0) ? 0.0 : output_temp;

    #if !defined(__SYNTHESIS__) && defined(AC_SQRT_PWL_H_DEBUG)
    std::cout << "W = " << W << std::endl;
    std::cout << "I = " << I << std::endl;
    std::cout << "outW = " << outW << std::endl;
    std::cout << "outI = " << outI << std::endl;
    std::cout << "input to normalization function = " << input << std::endl;
    std::cout << "output (fractional of normalization function = " << normalized_input << std::endl;
    std::cout << "input_sc = " << input_sc << std::endl;
    std::cout << "index of element chosen from ROM = " << index << std::endl;
    std::cout << "normalized_output = " << normalized_output << std::endl;
    std::cout << "normalized_output_temp = " << normalized_output_temp << std::endl;
    std::cout << "m1 = " << m1 << std::endl;
    std::cout << "normalized_exp = " << normalized_exp << std::endl;
    std::cout << "final output = " << output << std::endl;
    #endif
  }

  // This struct provides parameterized bitwidths to ensure a lossless return type for the monotonous PWL function provided by default,
  // that operates with 4 segments and uses 12 fractional bits to store slope and intercept values.
  // n_f_b is the number of fractional bits and I is the number of integer bits in the input. The input and output are assumed to be
  // unsigned. Other PWL implementations might require different calculations for the parameterized bitwidths.
  template <int n_f_b, int I>
  struct find_rt_sqrt_pwl {
    enum {
      // An extra bit is added for I1 for even integer widths because as we approach 1 for the normalized input value, the output for the sqrt function can slightly exceed
      // the expected value, taking into account the upward shifting of PWL segments against the direction of concavity for the sqrt function.
      // This addition of an extra bit is only required for even integer widths because the maximum value fixed point configurations with even integer widths can store
      // approaches an even power of two. This maximum value would correspond to an input to the PWL that would approach the maximum limit of the PWL domain after normalization.
      // For the default implementation, this maximum limit is 1. The corresponding output value slightly exceeds 1 owing to the upward shifting of PWL segments as discussed earlier.
      // If we don't add the extra bit, this slight increase can result in an overflow at the output.
      I1 = I % 2 == 0 ? (I/2) + 1 : (I + 1)/2,
      n_f_b_floor = n_f_b % 2 == 0 ? n_f_b / 2 : (n_f_b - 1) / 2,
      W1 = I1 + n_f_b_floor + 22
    };
    typedef ac_fixed<W1, I1, false> rt_sqrt_pwl;
  };

//=========================================================================
// Function: ac_sqrt_pwl (for ac_float)
//
// Description:
//    Calculation of square root of positive real inputs, passed as ac_float
//    variables.
//
//    This function uses the piecewise linear implementation by separating the
//    mantissa and exponent. Exponent is simply divided by two, if it is even or
//    is made even and then multiplied by square root of 2, which is stored as
//    as a constant, where as mantissa undergoes piecewise linear implementation
//    using helper function defined above.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_sqrt_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_float<16, 10, 8, AC_RND> input_type;
//    typedef ac_float<18, 20, 9, AC_RND> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type input,
//      output_type &output
//    )
//    {
//      ac_sqrt_pwl(input,output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 1.25;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template <ac_q_mode pwlQ = AC_TRN, int W, int I, int E, ac_q_mode Q, int outW, int outI, int outE, ac_q_mode outQ>
  void ac_sqrt_pwl (const ac_float <W, I, E, Q> input, ac_float <outW, outI, outE, outQ> &output)
  {
    // Use a macro to activate the AC_ASSERT
    // If AC_ASSERT is activated, the program will stop running as soon as a negative input
    // is encountered.
    #ifdef ASSERT_ON_INVALID_INPUT
    AC_ASSERT(input >= 0, "Negative input not supported.");
    #endif

    const ac_fixed <16, 1, false> root2 = 1.414215087890625;

    // mantissa and exponent are separately considered
    ac_fixed <W - 1, I - 1, false, Q> mantissa = input.mantissa();

    const int W1 = find_rt_sqrt_pwl<W - I, I - 1>::W1;
    const int I1 = find_rt_sqrt_pwl<W - I, I - 1>::I1;
    // declaring variable to store square root of mantissa
    ac_fixed <W1, I1, false> m2;

    const bool call_normalize = false; // ac_float input -> Mantissa is already normalized and fixed point function doesn't have to call ac_normalize.
    ac_sqrt_pwl<pwlQ>(mantissa, m2, call_normalize); // Call the ac_fixed implementation to get square root of mantissa

    // Multiplication by root 2 for odd exponent
    // Note that an extra bit does not need to be added to the word and integer widths for the product, for odd integer widths at the input. This is because the addition of the
    // extra bit in such a case is already done in the find_rt_sqrt_pwl struct, to take into account the fact that the PWL segments for the default (4 segments
    // and 12 fractional bits) are shifted slightly upward. For any other PWL implementation, this might change and the user might have to add the extra bit.
    typedef ac_fixed<W1 + int(I%2 == 0), I1 + int(I%2 == 0), false> type_prod;
    type_prod m3 = (input.exp() % 2 == 0) ? (type_prod)m2 : (type_prod)(m2 * root2);

    // The mantissa without normalization will be either the square root of the original mantissa, or that square root multiplied by sqrt(2),
    // depending upon whether the input exponent is even or not.
    // The exponent without normalization will be the input exponent right-shifted by 1/divided by 2. This follows the formula:
    // sqrt(mant * (2^exp)) = sqrt(mant) * (2^(exp/2))
    // These two values are passed to an ac_float constructor that takes care of normalization.
    ac_float <outW, outI, outE, outQ> output_temp(m3, input.exp() >> 1, true);

    output = output_temp;

    #if !defined(__SYNTHESIS__) && defined(AC_SQRT_PWL_FL_H_DEBUG)
    std::cout << "input = " << input << std::endl;
    std::cout << "W = " << W << std::endl;
    std::cout << "I = " << I << std::endl;
    std::cout << "W1 = " << W1 << std::endl;
    std::cout << "I1 = " << I1 << std::endl;
    std::cout << "m2.type_name() : " << m2.type_name() << std::endl;
    std::cout << "output of call to ac_fixed version of sqrt_pwl = " << m2 << std::endl;
    std::cout << "m3.type_name() : " << m3.type_name() << std::endl;
    std::cout << "m3 = " << m3 << std::endl;
    std::cout << "output_temp.m = " << output_temp.m << std::endl;
    std::cout << "output_temp.e = " << output_temp.e << std::endl;
    std::cout << "final output = " << output << std::endl;
    #endif
  }

// For this section of the code to work, the user must include ac_std_float.h in their testbench before including the square root header,
// so as to have the code import the ac_std_float and ac_ieee_float datatypes and define the __AC_STD_FLOAT_H macro.
  #ifdef __AC_STD_FLOAT_H
//=========================================================================
// Function: ac_sqrt_pwl (for ac_std_float)
//
// Description:
//    Calculation of square root of real inputs, passed as ac_std_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_sqrt_pwl.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_sqrt_pwl.h>
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
//      ac_sqrt_pwl(input, output);
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
  void ac_sqrt_pwl(
    const ac_std_float<W, E> &input,
    ac_std_float<outW, outE> &output
  )
  {
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_sqrt_pwl<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_std_float<outW, outE> output_temp(output_ac_fl); // Convert output ac_float to ac_std_float.
    output = output_temp;
  }

//=================================================================================
// Function: ac_sqrt_pwl (for ac_ieee_float, returns sqrt(input) )
//
// Description:
//    Calculation of square root of real inputs, passed as ac_ieee_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_sqrt_pwl.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_sqrt_pwl.h>
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
//      ac_sqrt_pwl(input,output);
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
//    The ac_ieee_float version of ac_sqrt_pwl relies on the ac_float version to
//    perform the actual PWL computation, which in turn relies on the ac_fixed
//    implementation.
//
//---------------------------------------------------------------------------------

  template<ac_q_mode pwl_Q = AC_TRN, ac_ieee_float_format Format, ac_ieee_float_format outFormat>
  void ac_sqrt_pwl(const ac_ieee_float<Format> input, ac_ieee_float<outFormat> &output)
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int outW = T_out::width;
    const int outE = T_out::e_width;
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_sqrt_pwl<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_ieee_float<outFormat> output_temp(output_ac_fl); // Convert output ac_float to ac_ieee_float.
    output = output_temp;

    #if !defined(__SYNTHESIS__) && defined(AC_SQRT_PWL_H_DEBUG)
    std::cout << "input.to_ac_float().type_name() : " << input.to_ac_float().type_name() << std::endl;
    std::cout << "output_ac_fl.type_name() : " << output_ac_fl.type_name() << std::endl;
    #endif
  }
  #endif

//=====================================================================================
// Function: ac_sqrt_pwl (for ac_complex <ac_fixed>)
//
// Description:
//    Calculation of square root of positive complex inputs, passed as ac_complex
//    data with real and imaginary part as ac_fixed type of data.
//
//    This function uses following mathematical formula for computation of square root:
//
//    output_real_part = sqrt((input_real_part + sqrt(input.mag_sqr()))/2.0)
//    output_imaginary_part = sqrt((-input_real_part + sqrt(input.mag_sqr()))/2.0)
//    Then square root is calculated by using the PWL model for ac_fixed values.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_sqrt_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_complex<ac_fixed<16, 8, true, AC_RND, AC_SAT> > input_type;
//    typedef ac_complex<ac_fixed<18, 8, true, AC_RND, AC_SAT> > output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type input,
//      output_type &output
//    )
//    {
//      ac_sqrt_pwl(input,output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input(1.0, 2.0);
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------------------

  // Real/imaginary part of input and output must be signed.
  template <ac_q_mode pwlQ = AC_TRN, int W, int I, ac_q_mode Q, ac_o_mode O, int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  void ac_sqrt_pwl (const ac_complex <ac_fixed <W, I, true, Q, O> > input,  ac_complex <ac_fixed <outW, outI, true, outQ, outO> > &output)
  {
    // Calculate parameterized bitwidths for all intermediate types.
    typedef class find_rt_sqrt_pwl<(2*(W - I)), (2*I - 1)>::rt_sqrt_pwl sqrt_mod_type;
    const int W1 = sqrt_mod_type::width;
    const int I1 = sqrt_mod_type::i_width;
    const int n_f_b_1 = W1 - I1;
    const int t_I = I + 1;
    const int t_n_f_b = (W - I) > n_f_b_1 ? W - I : n_f_b_1;
    const int t_W = t_I + t_n_f_b;
    typedef class find_rt_sqrt_pwl<t_W, t_I - 1>::rt_sqrt_pwl x_y_type;
    const int W2 = x_y_type::width;
    const int I2 = x_y_type::i_width;

    sqrt_mod_type sqrt_mod;
    ac_sqrt_pwl <pwlQ> (input.mag_sqr(), sqrt_mod);     // computation of square root of mag_sqr
    // Since the output of the PWL function is not exact, temp_real and temp_imag can possibly go negative.
    // Since negative inputs aren't supported for the real sqrt function, we must ensure that such a thing
    // never happens. To do so, we use a type that supports saturation, as shown below. By doing
    // so, we can ensure that the temporary variables saturate at 0 instead of wrapping around to
    // undesirable values.
    ac_fixed <t_W, t_I, false, AC_TRN, AC_SAT> temp_real = sqrt_mod + input.r();
    ac_fixed <t_W, t_I, false, AC_TRN, AC_SAT> temp_imag = sqrt_mod - input.r();
    ac_fixed <t_W, t_I - 1, false> sqr_real = ((ac_fixed <t_W + 1, t_I, false>)temp_real) >> 1; // calculating square of the output's real part
    ac_fixed <t_W, t_I - 1, false> sqr_imag = ((ac_fixed <t_W + 1, t_I, false>)temp_imag) >> 1; // calculating square of the output's imaginary part
    x_y_type x, y;
    ac_sqrt_pwl <pwlQ> (sqr_real, x); // calculating output's real part
    ac_sqrt_pwl <pwlQ> (sqr_imag, y); // calculating output's imaginary part
    output.r() = x;
    output.i() = (input.i() < 0) ? -y : (ac_fixed <W2 + 1, I2 + 1, true>)y; // if imaginary part is less than zero, assign output value as negative otherwise positive

    #if !defined(__SYNTHESIS__) && defined(AC_SQRT_PWL_H_DEBUG)
    std::cout << "initial input = " << input << std::endl;
    std::cout << "W1 = " << W1 << std::endl;
    std::cout << "I1 = " << I1 << std::endl;
    std::cout << "Value of square root of mag_sqr = " << sqrt_mod << std::endl;
    std::cout << "Type of square root of mag_sqr = " << sqrt_mod.type_name() << std::endl;
    std::cout << "Result of addition = " << temp_real << std::endl;
    std::cout << "Result of subtraction = " << temp_imag << std::endl;
    std::cout << "Type of sqr_real and sqr_imag = " << sqr_real.type_name() << std::endl;
    std::cout << "Result of square of output real part = " << sqr_real << std::endl;
    std::cout << "Result of square of output imaginary part = " << sqr_imag << std::endl;
    std::cout << "Type of real and imaginary part of answer = " << x.type_name() << std::endl;
    std::cout << "Absolute value of real part = " << x << std::endl;
    std::cout << "Absolute value of imaginary part = " << y << std::endl;
    std::cout << "Final value of square root of complex number = " << output << std::endl;
    #endif
  }

//=========================================================================
// Version that allows returning of values
  template<class T_out, ac_q_mode pwlQ = AC_TRN, class T_in>
  T_out ac_sqrt_pwl(
    const T_in input
  )
  {
    // Initializing the final output value that is to be returned
    T_out output;
    // Call the function by referencing the output variable. This is call to one of above implementations
    ac_sqrt_pwl<pwlQ>(input, output);
    // Return the final computed output
    return output;
  }
}

#endif

