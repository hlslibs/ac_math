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
// File: ac_reciprocal_pwl_ha.h
//
// Description: Provides high-accuracy piece-wise linear implementations of the
//    reciprocal function for the AC (tm) Datatypes: ac_fixed, ac_float,
//    ac_complex<ac_fixed>, ac_complex<ac_float>, ac_ieee_float.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_reciprocal_pwl_ha.h>
//    using namespace ac_math;
//
//    typedef ac_fixed<20, 11, true, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<24, 14, true, AC_RND, AC_SAT> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_reciprocal_pwl_ha(input, output);
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
//    This file uses C++ function overloading to target implementations
//    specific to each type of data. Attempting to call the function
//    with a type that is not implemented will result in a compile error.
//
//    This library uses the ac_normalize() and ac_shift_right() function
//    from the other ac_math header files.
//
// Revision History:
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//    3.2.3 - Added to help resolve CAT-24022. Further revised in accordance to CAT-24139.
//
//*****************************************************************************************

#ifndef _INCLUDED_AC_RECIPROCAL_PWL_HA_H_
#define _INCLUDED_AC_RECIPROCAL_PWL_HA_H_

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
#include <ac_complex.h>

// Include headers for required functions
#include <ac_math/ac_normalize.h>
#include <ac_math/ac_shift.h>

#if !defined(__SYNTHESIS__) && defined(AC_RECIPROCAL_PWL_HA_H_DEBUG)
#include <iostream>
#endif

//=========================================================================
// Function: ac_reciprocal_pwl_ha (for ac_fixed)
//
// Description:
//    Calculation of reciprocal of real inputs, passed as ac_fixed
//    variables.
//
//    Passes inputs for normalization to the function in ac_normalize.h,
//    which gives the exponent and output normalized between 0.5 and 1.
//    The normalized value is then subject to the piecewise linear
//    implementation to calculate the reciprocal.
//
//    This reciprocal is then up-scaled according to the exponent value
//    returned by the normalization function, and returned to the calling
//    function.
//
// Usage:
//    See above example code for usage.
//
// Notes:
//    The PWL implementation utilizes 7 elements, which has a small impact
//    on accuracy.
//
//-------------------------------------------------------------------------

namespace ac_math
{
  template<ac_q_mode pwl_Q = AC_TRN,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO>
  void ac_reciprocal_pwl_ha(
    const ac_fixed<W, I, S, Q, O> &input,
    ac_fixed<outW, outI, outS, outQ, outO> &output
  )
  {
    // Use a macro to activate the AC_ASSERT
    // If AC_ASSERT is activated: the program will stop running as soon as a zero input
    // is encountered.
    // If AC_ASSERT is not activated: the output will saturate when a zero input is encountered.
    // The functionality behind this is taken care of by other sections of the code.
    #ifdef ASSERT_ON_INVALID_INPUT
    AC_ASSERT(!!input, "Reciprocal of zero not supported.");
    #endif

    ac_fixed<outW, outI, outS, outQ, outO> output_temp;

    // Start of code outputted by ac_reciprocal_pwl_ha_lutgen.cpp
    const unsigned n_segments_lut = 32; // Number of PWL segments.
    const int n_frac_bits = 16; // Number of fractional bits
    // Since scaling constant is a positive power-of-two, multiplication with it is the same as left-shifting by 6.
    // Accordingly, the scaled normalized input will have 6 less fractional bits than the normalized input, provided that this
    // number of fractional bits is lesser than n_frac_bits. If not, the number of fractional bits in the scaled input is set to n_frac_bits.
    const int sc_input_frac_bits = AC_MAX(1, AC_MIN(n_frac_bits, W - 6 - int(S))); // One less bit is used if the function input is signed, due to how ac_normalize works.
    static const ac_fixed<-2 + n_frac_bits, -2, true> m_lut[n_segments_lut] = {
      -.06060791015625, -.0570220947265625, -.0537567138671875, -.0507659912109375, -.04803466796875, -.045501708984375, -.0431671142578125, -.041015625,
        -.0390167236328125, -.0371551513671875, -.035430908203125, -.0338134765625, -.032318115234375, -.030914306640625, -.02960205078125, -.0283660888671875,
        -.0272064208984375, -.0261077880859375, -.0251007080078125, -.0241241455078125, -.023223876953125, -.0223541259765625, -.02154541015625, -.020782470703125,
        -.020050048828125, -.01934814453125, -.018707275390625, -.01806640625, -.017486572265625, -.0169219970703125, -.016387939453125, -.015869140625
      };
    static const ac_fixed<2 + n_frac_bits, 2, false> c_lut[n_segments_lut] = {
      1.9997711181640625, 1.9391632080078125, 1.88214111328125, 1.828369140625, 1.7776031494140625, 1.7295684814453125, 1.684051513671875, 1.6408843994140625,
      1.5998687744140625, 1.56085205078125, 1.5236968994140625, 1.4882659912109375, 1.4544525146484375, 1.4221343994140625, 1.3912200927734375, 1.3616180419921875,
      1.333251953125, 1.3060455322265625, 1.279937744140625, 1.2548370361328125, 1.230712890625, 1.207489013671875, 1.1851348876953125, 1.1635894775390625,
      1.1428070068359375, 1.1227569580078125, 1.1034088134765625, 1.0847015380859375, 1.0666351318359375, 1.0491485595703125, 1.0322265625, 1.015838623046875
    };
    static const ac_fixed<0 + n_frac_bits, 0, false> x_min_lut = .5; // Minimum limit of PWL domain
    static const ac_fixed<7 + n_frac_bits, 7, false> sc_constant_lut = 64.0; // Scaling constant
    // End of code outputted by ac_reciprocal_pwl_ha_lutgen.cpp

    // The absolute value of the input is taken and passed to the normalization function. Initialize variables for the same.
    ac_fixed<W, I, false> input_abs_value;
    ac_fixed<W, 0, false> normalized_fixed;

    // If input is signed, take absolute value and assign to intermediate variable.
    if (S) {
      input_abs_value = ((input >= 0) ? (ac_fixed <W, I, false>)input : (ac_fixed <W, I, false>)(-input));
    }
    // If input is unsigned, assign value of input to intermediate variable.
    else {
      input_abs_value = input;
    }

    // Normalize the absolute value. expret stores the value of the returned base 2 exponential.
    int expret_temp = ac_math::ac_normalize(input_abs_value, normalized_fixed);

    const int int_bits = ac::nbits<n_segments_lut - 1>::val;
    // Compute reciprocal using pwl.
    // Scale the normalized input from 0 to n_segments_lut.
    ac_fixed<int_bits + sc_input_frac_bits, int_bits, false> x_in_sc = (normalized_fixed - x_min_lut)*sc_constant_lut;
    // Take out the fractional bits of the scaled input
    ac_fixed<sc_input_frac_bits, 0, false> x_in_sc_frac;
    x_in_sc_frac.set_slc(0, x_in_sc.template slc<sc_input_frac_bits>(0));
    // The integer part of the input is the index of the LUT table
    ac_int<int_bits, false> index = x_in_sc.to_int();
    // The output of the PWL approximation should have the same signedness as the output of the function.
    // The precision given below will ensure that there is no precision lost in the assignment to output_pwl, hence rounding for the variable is switched off by default.
    // However, if the user uses less fractional bits and turn rounding on instead, they are welcome to do so by giving a different value for pwl_Q.
    typedef ac_fixed<sc_input_frac_bits + n_frac_bits + 1 + int(outS), 1 + int(outS), outS, pwl_Q> output_pwl_type;
    output_pwl_type output_pwl = m_lut[index]*x_in_sc_frac + c_lut[index];

    if (input != 0) { // If input is non-zero, De-normalize output by shifting right by expret_temp
      // If input and output are signed, change sign of output_pwl based on whether input is positive or negative.
      if (S && outS) { output_pwl = (input < 0) ? (output_pwl_type)(-output_pwl) : output_pwl; }
      // ac_shift_right function used for denormalization so as to ensure saturation and rounding.
      ac_math::ac_shift_right(output_pwl, expret_temp, output_temp);
    } else {
      // If zero input is encountered, set output to the max possible value.
      output_temp.template set_val<AC_VAL_MAX>();
    }

    output = output_temp;

    #if !defined(__SYNTHESIS__) && defined(AC_RECIPROCAL_PWL_HA_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "input                   = " << input << std::endl;
    std::cout << "input_abs_value         = " << input_abs_value << std::endl;
    std::cout << "normalized input        = " << normalized_fixed << std::endl;
    std::cout << "expret_temp             = " << expret_temp << std::endl;
    std::cout << "x_in_sc                 = " << x_in_sc << std::endl;
    std::cout << "x_in_sc_frac            = " << x_in_sc_frac << std::endl;
    std::cout << "output_pwl              = " << output_pwl << std::endl;
    std::cout << "output_temp             = " << output_temp << std::endl;
    std::cout << "output up-scaled by exp = " << output << std::endl;
    #endif
  }

//=========================================================================
// Function: ac_reciprocal_pwl_ha (for ac_float)
//
// Description:
//    Calculation of reciprocal of real inputs, passed as ac_float
//    variables.
//
//    The mantissa of the ac_float number is passed as an ac_fixed variable
//    to the ac_reciprocal_pwl_ha function for ac_fixed numbers. The mantissa of
//    the output is then set to the resultant reciprocal. The exponent
//    of the output ac_float variable is set to the negative of the
//    exponent of the input.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_reciprocal_pwl_ha.h>
//    using namespace ac_math;
//
//    typedef ac_float<20, 11, 7, AC_RND> input_type;
//    typedef ac_float<24, 14, 7, AC_RND> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_reciprocal_pwl_ha(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 1.2;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Note:
//    This implementation relies on the ac_fixed implementation of the
//    reciprocal function.
//
//-------------------------------------------------------------------------

  template<ac_q_mode pwl_Q = AC_TRN,
           int W, int I, int E, ac_q_mode Q,
           int outW, int outI, int outE, ac_q_mode outQ>
  void ac_reciprocal_pwl_ha(
    const ac_float<W, I, E, Q> &input,
    ac_float<outW, outI, outE, outQ> &output
  )
  {
    // Start of code to be changed if the ac_fixed PWL implementation changes.

    const bool is_sc_constant_po2 = true; // is sc_constant_lut used in ac_fixed implementation a power of two?
    const int log2_sc_constant = 6; // log2(sc_constant_lut), if sc_constant_lut is a power-of-two. If not, this value is ignored.
    const int n_frac_bits = 16; // Number of fractional bits used in ac_fixed PWL implementation.

    // End of code to be changed if the ac_fixed PWL implementation changes.

    const int sc_input_frac_bits = is_sc_constant_po2 ? AC_MAX(1, AC_MIN(n_frac_bits, W - log2_sc_constant - 1)) : n_frac_bits;

    const int W1 = W + n_frac_bits + sc_input_frac_bits + 1;
    const int I1 = W - I + 2;

    // Find the reciprocal of the mantissa using the ac_fixed implementation.
    ac_fixed<W1, I1, true, outQ> recip_mantissa;
    ac_reciprocal_pwl_ha<pwl_Q>(input.mantissa(), recip_mantissa);

    // Find the additive inverse of the input's exponent.
    // Pass it and recip_mantissa to an ac_float constructor that takes care of
    // normalization.
    ac_float<outW, outI, outE, outQ> output_temp(recip_mantissa, -input.exp(), true);

    // If the input is zero, set the temp output to the max. possible value.
    if (input.mantissa() == 0) { output_temp.template set_val<AC_VAL_MAX>(); }

    output = output_temp;

    #if !defined(__SYNTHESIS__) && defined(AC_RECIPROCAL_PWL_HA_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "input                   = " << input << std::endl;
    std::cout << "input.mantissa()        = " << input.mantissa() << std::endl;
    std::cout << "input.exp()             = " << input.exp() << std::endl;
    std::cout << "recip_mantissa          = " << recip_mantissa << std::endl;
    std::cout << "output.type_name()      = " << output.type_name() << std::endl;
    std::cout << "output_temp.type_name() = " << output_temp.type_name() << std::endl;
    std::cout << "output_temp             = " << output_temp << std::endl;
    std::cout << "output                  = " << output << std::endl;
    std::cout << "output.mantissa()       = " << output.mantissa() << std::endl;
    std::cout << "output.exp()            = " << output.exp() << std::endl;
    #endif
  }

//=========================================================================
// Function: ac_reciprocal_pwl_ha (for ac_complex<ac_fixed>)
//
// Description:
//    Calculation of reciprocal of complex inputs, passed as ac_complex
//    variables with ac_fixed real and imaginary parts.
//
//    The reciprocal is calculated using the formula:
//
//    1 / (a + bi) = (a - bi) / (a^2 + b^2)
//
//    Where a and b are the real and imaginary parts, respectively.In
//    order to do this, the value of a^2 + b^2 is first calculated, and the
//    reciprocal of this real number is calculated by passing it to the
//    ac_reciprocal_pwl_ha function for ac_fixed numbers.
//
//    The resultant reciprocal is multiplied by the real and imaginary
//    parts of the number passed, the results of which are then passed to
//    the output.
//
// Usage:
//
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_reciprocal_pwl_ha.h>
//    using namespace ac_math;
//
//    typedef ac_complex<ac_fixed<20, 11, true, AC_RND, AC_SAT> > input_type;
//    typedef ac_complex<ac_fixed<40, 18, true, AC_RND, AC_SAT> > output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_reciprocal_pwl_ha(input, output);
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
// Note:
//    This implementation relies on the ac_fixed implementation of the
//    reciprocal function.
//
//-------------------------------------------------------------------------

  template<ac_q_mode pwl_Q = AC_TRN,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO>
  void ac_reciprocal_pwl_ha(
    const ac_complex<ac_fixed<W, I, S, Q, O> > &input,
    ac_complex<ac_fixed<outW, outI, outS, outQ, outO> > &output
  )
  {
    const int W1 = S ? 2*W - 1 : 2*W + 1;
    const int I1 = S ? 2*I - 1 : 2*I + 1;

    // Start of code to be changed if the ac_fixed PWL implementation changes.

    const bool is_sc_constant_po2 = true; // is sc_constant_lut used in ac_fixed implementation a power of two?
    const int log2_sc_constant = 6; // log2(sc_constant_lut), if sc_constant_lut is a power-of-two. If not, this value is ignored.
    const int n_frac_bits = 16; // Number of fractional bits used in ac_fixed PWL implementation.

    // End of code to be changed if the ac_fixed PWL implementation changes.

    const int sc_input_frac_bits = is_sc_constant_po2 ? AC_MAX(1, AC_MIN(n_frac_bits, W1 - log2_sc_constant - int(S))) : n_frac_bits;

    // The derived type for the reciprocal of the mag_sqr() of the input has its
    // bitwidths calculated to ensure a lossless return type for said reciprocal.
    const int W2 = W1 + sc_input_frac_bits + n_frac_bits + int(!S);
    const int I2 = W1 - I1 + 1;

    ac_fixed<W2, I2, false, outQ, outO> recip_mag_sqr;
    ac_reciprocal_pwl_ha<pwl_Q>(input.mag_sqr(), recip_mag_sqr);

    ac_complex<ac_fixed<outW, outI, outS, outQ, outO> > output_temp;

    if (input.r() != 0 || input.i() != 0) {
      // Use the formula "1/(a+bi) = (a-bi)/(a^2+b^2)" to assign values to the output.
      output_temp.r() =  input.r() * recip_mag_sqr;
      output_temp.i() = -input.i() * recip_mag_sqr;
    } else {
      // If zero input is passed, then assign the maximum possible value for output's real and imaginary part
      output_temp.r().template set_val<AC_VAL_MAX>();
      output_temp.i().template set_val<AC_VAL_MAX>();
    }

    output = output_temp;

    #if !defined(__SYNTHESIS__) && defined(AC_RECIPROCAL_PWL_HA_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "input.mag_sqr() = " << input.mag_sqr() << std::endl;
    std::cout << "recip_mag_sqr   = " << recip_mag_sqr << std::endl;
    std::cout << "output_temp     = " << output_temp << std::endl;
    std::cout << "output          = " << output << std::endl;
    #endif
  }

//=========================================================================
// Function: ac_reciprocal_pwl_ha (for ac_complex<ac_float>)
//
// Description:
//    Calculation of reciprocal of complex inputs, passed as ac_complex
//    variables with ac_float real and imaginary parts.
//
//    The reciprocal is calculated using the formula:
//
//    (a + bi) = (a - bi) / (a^2 + b^2)
//
//    Where a and b are the real and imaginary parts, respectively. In
//    order to do this, the value of a^2 + b^2 is first calculated, and the
//    reciprocal of this real number is calculated using the
//    ac_reciprocal_pwl_ha function for ac_float numbers.
//
//    The resultant reciprocal is multiplied by the real and imaginary
//    parts of the number passed, the results of which are then passed to
//    the output.
//
// Usage:
//
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_reciprocal_pwl_ha.h>
//    using namespace ac_math;
//
//    typedef ac_complex<ac_float<20, 11, 7, AC_RND> > input_type;
//    typedef ac_complex<ac_float<40, 18, 7, AC_RND> > output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_reciprocal_pwl_ha(input, output);
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
//      CCS_RETURN(0);
//    }
//    #endif
//
// Note:
//    This implementation relies on both the ac_float and ac_fixed
//    implementations of the reciprocal function.
//
//-------------------------------------------------------------------------

  template<ac_q_mode pwl_Q = AC_TRN,
           int W, int I, int E, ac_q_mode Q,
           int outW, int outI, int outE, ac_q_mode outQ>
  void ac_reciprocal_pwl_ha(
    const ac_complex<ac_float<W, I, E, Q> > &input,
    ac_complex<ac_float<outW, outI, outE, outQ> > &output
  )
  {
    // Start of code to be changed if the ac_fixed PWL implementation changes.

    const bool is_sc_constant_po2 = true; // is sc_constant_lut used in ac_fixed implementation a power of two?
    const int log2_sc_constant = 6; // log2(sc_constant_lut), if sc_constant_lut is a power-of-two. If not, this value is ignored.
    const int n_frac_bits = 16; // Number of fractional bits used in ac_fixed PWL implementation.

    // End of code to be changed if the ac_fixed PWL implementation changes.

    const int sc_input_frac_bits = is_sc_constant_po2 ? AC_MAX(1, AC_MIN(n_frac_bits, 2*W + 1 - log2_sc_constant)) : n_frac_bits;

    // Calculate real^2 + imag^2
    ac_float<W, I, E, Q> input_real = input.r();
    ac_float<W, I, E, Q> input_imag = input.i();

    // Define type for input_mag_sqr.
    // E + 2 is chosen for the exponent width of input_mag_sqr
    // Because it works well with all input exponential widths,
    // as well as corner cases with very low exponential widths (as low as 1 or 2)
    typedef ac_float<2*W + 1, 2*I + 1, E + 2, Q> i_m_s_type;
    i_m_s_type input_mag_sqr;

    // Calculate bitwidths for the return type of reciprocal of the input magnitude.
    const int W_r_m_s = 2*W + 1 + sc_input_frac_bits + n_frac_bits + 1;
    const int I_r_m_s = 2*W + 1 - (2*I + 1) + 2;

    ac_float<W_r_m_s, I_r_m_s, E + 3, outQ> recip_mag_sqr;

    // Store value of input_mag_sqr to the variable.
    input_mag_sqr.add(input_real*input_real, input_imag*input_imag);
    ac_reciprocal_pwl_ha<pwl_Q>(input_mag_sqr, recip_mag_sqr);

    ac_complex<ac_float<outW, outI, outE, outQ> > output_temp;

    if (input_mag_sqr.mantissa() != 0) {
      // 1/(a+bi) = (a-bi)/(a^2+b^2)
      output_temp.r() =  input.r() * recip_mag_sqr;
      output_temp.i() = -input.i() * recip_mag_sqr;
    } else {
      // If zero input is passed, then assign the maximum possible value for output's real part,
      // and a zero value for imaginary.
      output_temp.r().template set_val<AC_VAL_MAX>();
      output_temp.i().template set_val<AC_VAL_MAX>();
    }

    output = output_temp;

    #if !defined(__SYNTHESIS__) && defined(AC_RECIPROCAL_PWL_HA_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "input                     = " << input << std::endl;
    std::cout << "input_mag_sqr             = " << input_mag_sqr << std::endl;
    std::cout << "input_mag_sqr.to_double() = " << input_mag_sqr.to_double() << std::endl;
    std::cout << "recip_mag_sqr             = " << recip_mag_sqr << std::endl;
    std::cout << "recip_mag_sqr.to_double() = " << recip_mag_sqr.to_double() << std::endl;
    std::cout << "output_temp               = " << output_temp << std::endl;
    std::cout << "output                    = " << output << std::endl;
    std::cout << "output.r().to_double()    = " << output.r().to_double() << std::endl;
    std::cout << "output.i().to_double()    = " << output.i().to_double() << std::endl;
    #endif
  }

// For this section of the code to work, the user must include ac_std_float.h in their testbench before including the reciprocal header,
// so as to have the code import the ac_std_float and ac_ieee_float datatypes and define the __AC_STD_FLOAT_H macro.
  #ifdef __AC_STD_FLOAT_H
//=========================================================================
// Function: ac_reciprocal_pwl_ha (for ac_std_float)
//
// Description:
//    Calculation of reciprocal of real inputs, passed as ac_std_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_reciprocal_pwl_ha.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_reciprocal_pwl_ha.h>
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
//      ac_reciprocal_pwl_ha(input, output);
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
  void ac_reciprocal_pwl_ha(
    const ac_std_float<W, E> &input,
    ac_std_float<outW, outE> &output
  )
  {
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_reciprocal_pwl_ha<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_std_float<outW, outE> output_temp(output_ac_fl); // Convert output ac_float to ac_std_float.
    output = output_temp;
  }

//=========================================================================
// Function: ac_reciprocal_pwl_ha (for ac_ieee_float)
//
// Description:
//    Calculation of reciprocal of real inputs, passed as ac_ieee_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_reciprocal_pwl_ha.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_reciprocal_pwl_ha.h>
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
//      ac_reciprocal_pwl_ha(input, output);
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
  void ac_reciprocal_pwl_ha(
    const ac_ieee_float<Format> &input,
    ac_ieee_float<outFormat> &output
  )
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int outW = T_out::width;
    const int outE = T_out::e_width;
    ac_float<outW - outE + 1, 2, outE> output_ac_fl; // Equivalent ac_float representation for output.
    ac_reciprocal_pwl_ha<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_ieee_float<outFormat> output_temp(output_ac_fl); // Convert output ac_float to ac_ieee_float.
    output = output_temp;
  }
  #endif

//=========================================================================
// Version that allows returning of values.
  template<class T_out, ac_q_mode pwl_Q = AC_TRN, class T_in>
  T_out ac_reciprocal_pwl_ha(
    const T_in &input
  )
  {
    // Create an intermediate variable for output and use the pass-by-reference version
    // to evaluate it. This intermediate variable is returned as the output.
    T_out output;
    ac_reciprocal_pwl_ha<pwl_Q>(input, output);
    return output;
  }
}

#endif
