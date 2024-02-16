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
// File: ac_reciprocal_pwl_vha.h
//
// Description: Provides very high-accuracy piece-wise linear implementations of the
//    reciprocal function for the AC (tm) Datatypes: ac_fixed, ac_float,
//    ac_complex<ac_fixed>, ac_complex<ac_float>, ac_ieee_float.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_reciprocal_pwl_vha.h>
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
//      ac_reciprocal_pwl_vha(input, output);
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
//    3.4.0 - Added to help resolve CAT-27696.
//
//*****************************************************************************************

#ifndef _INCLUDED_AC_RECIPROCAL_PWL_VHA_H_
#define _INCLUDED_AC_RECIPROCAL_PWL_VHA_H_

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

#ifdef __SYNTHESIS__
#include <iostream>
#endif

//=========================================================================
// Function: ac_reciprocal_pwl_vha (for ac_fixed)
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
  void ac_reciprocal_pwl_vha(
    const ac_fixed<W, I, S, Q, O> &input,
    ac_fixed<outW, outI, outS, outQ, outO> &output,
    const bool call_normalize = true
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

    // Start of code outputted by ac_reciprocal_pwl_vha_lutgen.cpp

    const unsigned n_segments_lut = 64; // Number of PWL segments.
    const int n_frac_bits = 16; // Number of fractional bits
    // Since scaling constant is a positive power-of-two, multiplication with it is the same as left-shifting by 7.
    // Accordingly, the scaled normalized input will have 7 less fractional bits than the normalized input, provided that this
    // number of fractional bits is lesser than n_frac_bits. If not, the number of fractional bits in the scaled input is set to n_frac_bits.
    const int sc_input_frac_bits = AC_MAX(1, AC_MIN(n_frac_bits, W - 7 - int(S))); // One less bit is used if the function input is signed, due to how ac_normalize works.
    // Slope and intercept LUT values.
    const ac_fixed<-3 + n_frac_bits, -3, true> m_lut[n_segments_lut] = {
      -.03076171875, -.0298309326171875, -.0289306640625, -.0280914306640625, -.0272674560546875, -.0265045166015625, -.0257568359375, -.0250244140625,
        -.02435302734375, -.0236968994140625, -.0230560302734375, -.0224609375, -.0218658447265625, -.0213165283203125, -.0207672119140625, -.0202484130859375,
        -.019744873046875, -.0192718505859375, -.018798828125, -.0183563232421875, -.0179290771484375, -.0175018310546875, -.0171051025390625, -.0167083740234375,
        -.0163421630859375, -.0159759521484375, -.015625, -.015289306640625, -.01495361328125, -.0146484375, -.0143280029296875, -.0140380859375,
        -.0137481689453125, -.013458251953125, -.0131988525390625, -.0129241943359375, -.0126800537109375, -.012420654296875, -.012176513671875, -.0119476318359375,
        -.01171875, -.011505126953125, -.01129150390625, -.011077880859375, -.0108642578125, -.01068115234375, -.0104827880859375, -.010284423828125,
        -.0101165771484375, -.0099334716796875, -.009765625, -.0095977783203125, -.009429931640625, -.00927734375, -.0091094970703125, -.0089569091796875,
        -.0088043212890625, -.0086669921875, -.0085296630859375, -.0083770751953125, -.0082550048828125, -.00811767578125, -.00799560546875, -.00787353515625
      };
    const ac_fixed<2 + n_frac_bits, 2, false> c_lut[n_segments_lut] = {
      1.99993896484375, 1.96917724609375, 1.9393310546875, 1.910400390625, 1.882293701171875, 1.8550262451171875, 1.828521728515625, 1.802764892578125,
      1.777740478515625, 1.753387451171875, 1.7296905517578125, 1.706634521484375, 1.684173583984375, 1.6623077392578125, 1.6409912109375, 1.6202239990234375,
      1.5999755859375, 1.5802154541015625, 1.560943603515625, 1.542144775390625, 1.5237884521484375, 1.505859375, 1.48834228515625, 1.4712371826171875,
      1.45452880859375, 1.4381866455078125, 1.4221954345703125, 1.4065704345703125, 1.3912811279296875, 1.3763275146484375, 1.3616790771484375, 1.34735107421875,
      1.33331298828125, 1.3195648193359375, 1.3061065673828125, 1.29290771484375, 1.2799835205078125, 1.267303466796875, 1.2548828125, 1.242706298828125,
      1.2307586669921875, 1.2190399169921875, 1.2075347900390625, 1.1962432861328125, 1.1851654052734375, 1.1743011474609375, 1.1636199951171875, 1.15313720703125,
      1.142852783203125, 1.1327362060546875, 1.122802734375, 1.113037109375, 1.1034393310546875, 1.0940093994140625, 1.0847320556640625, 1.07562255859375,
      1.066650390625, 1.0578460693359375, 1.049163818359375, 1.0406341552734375, 1.032257080078125, 1.02398681640625, 1.015869140625, 1.00787353515625
    };
    const ac_fixed<0 + n_frac_bits, 0, false> x_min_lut = .5; // Minimum limit of PWL domain
    const ac_fixed<8 + n_frac_bits, 8, false> sc_constant_lut = 128.0; // Scaling constant

    // End of code outputted by ac_reciprocal_pwl_vha_lutgen.cpp

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

    int expret;
    // If call_normalize is set to true, the design does not assume that the input is already normalized and performs normalization by calling ac_normalize().
    // If call_normalize is set to false, the design assumes that the input is already normalized and does not call ac_normalize, thereby saving on hardware.
    #pragma hls_waive CNS
    if (call_normalize) {
      expret = ac_math::ac_normalize(input_abs_value, normalized_fixed);
    } else {
      if (S && input_abs_value[W - 1] == 0) {
        // If input is a positive signed number, you start slicing from the bit adjacent to the LSB and not the LSB, due to the signed bit being rendered irrelevant
        // after the absolute value operation.
        expret = I - int(S);
        normalized_fixed = 0.0;
        normalized_fixed.set_slc(int(S), input_abs_value.template slc<W - int(S)>(0));
      } else {
        expret = I;
        normalized_fixed.set_slc(0, input_abs_value.template slc<W>(0));
      }
    }

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

    ac_fixed<outW, outI, outS, outQ, outO> output_temp;

    if (input != 0) { // If input is non-zero, De-normalize output by shifting right by expret
      // If input and output are signed, change sign of output_pwl based on whether input is positive or negative.
      if (S && outS) {output_pwl = (input < 0) ? (output_pwl_type)(-output_pwl) : output_pwl;}
      // ac_shift_right function used for denormalization so as to ensure saturation and rounding.
      ac_math::ac_shift_right(output_pwl, expret, output_temp);
    } else {
      // If zero input is encountered, set output to the max possible value.
      output_temp.template set_val<AC_VAL_MAX>();
    }

    output = output_temp;

    #if !defined(__SYNTHESIS__) && defined(AC_RECIPROCAL_PWL_VHA_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "input                   = " << input << std::endl;
    std::cout << "input_abs_value         = " << input_abs_value << std::endl;
    std::cout << "normalized input        = " << normalized_fixed << std::endl;
    std::cout << "expret                  = " << expret << std::endl;
    std::cout << "x_in_sc                 = " << x_in_sc << std::endl;
    std::cout << "x_in_sc_frac            = " << x_in_sc_frac << std::endl;
    std::cout << "output_pwl              = " << output_pwl << std::endl;
    std::cout << "output_temp             = " << output_temp << std::endl;
    std::cout << "output up-scaled by exp = " << output << std::endl;
    #endif
  }

//=========================================================================
// Function: ac_reciprocal_pwl_vha (for ac_float)
//
// Description:
//    Calculation of reciprocal of real inputs, passed as ac_float
//    variables.
//
//    The mantissa of the ac_float number is passed as an ac_fixed variable
//    to the ac_reciprocal_pwl_vha function for ac_fixed numbers. The mantissa of
//    the output is then set to the resultant reciprocal. The exponent
//    of the output ac_float variable is set to the negative of the
//    exponent of the input.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_reciprocal_pwl_vha.h>
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
//      ac_reciprocal_pwl_vha(input, output);
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
  void ac_reciprocal_pwl_vha(
    const ac_float<W, I, E, Q> &input,
    ac_float<outW, outI, outE, outQ> &output
  )
  {
    // Start of code to be changed if the ac_fixed PWL implementation changes.

    const bool is_sc_constant_po2 = true; // is sc_constant_lut used in ac_fixed implementation a power of two?
    const int log2_sc_constant = 7; // log2(sc_constant_lut), if sc_constant_lut is a power-of-two. If not, this value is ignored.
    const int n_frac_bits = 16; // Number of fractional bits used in ac_fixed PWL implementation.

    // End of code to be changed if the ac_fixed PWL implementation changes.

    const int sc_input_frac_bits = is_sc_constant_po2 ? AC_MAX(1, AC_MIN(n_frac_bits, W - log2_sc_constant - 1)) : n_frac_bits;

    const int E_max = 1 << (E - 1);
    const int E_temp_min = ac::nbits<AC_MAX(W - 1 + E_max, 1)>::val;
    const int E_temp = E_temp_min + 1;

    const bool normalize_input = true;
    // input_temp has enough exponent bits to account for worst-case normalization in the input
    // and have a fully normalized mantissa itself.
    ac_float<W, I, E_temp> input_temp(input.mantissa(), input.exp(), normalize_input);

    const int W1 = sc_input_frac_bits + n_frac_bits + 2;
    const int I1 = -I + 3;

    // Find the reciprocal of the mantissa using the ac_fixed implementation.
    ac_fixed<W1, I1, true, outQ> recip_mantissa;

    // Call ac_fixed version to find reciprocal of normalized mantissa.
    // ac_float mantissa is already normalized -> ac_fixed version doesn't have to call ac_normalize.
    const bool call_normalize = false;
    ac_reciprocal_pwl_vha<pwl_Q>(input_temp.mantissa(), recip_mantissa, call_normalize);

    // Find the additive inverse of the input's exponent.
    // Pass it and recip_mantissa to an ac_float constructor that takes care of
    // normalization.
    ac_float<outW, outI, outE, outQ> output_temp(recip_mantissa, -input_temp.exp(), true);

    // If the input is zero, set the temp output to the max. possible value.
    if (input.mantissa() == 0) { output_temp.template set_val<AC_VAL_MAX>(); }

    output = output_temp;

    #if !defined(__SYNTHESIS__) && defined(AC_RECIPROCAL_PWL_VHA_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "input                   = " << input << std::endl;
    std::cout << "input.mantissa()        = " << input.mantissa() << std::endl;
    std::cout << "input.exp()             = " << input.exp() << std::endl;
    std::cout << "input_temp              = " << input_temp << std::endl;
    std::cout << "input_temp.type_name()  = " << input_temp.type_name() << std::endl;
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
// Function: ac_reciprocal_pwl_vha (for ac_complex<ac_fixed>)
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
//    ac_reciprocal_pwl_vha function for ac_fixed numbers.
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
//    #include <ac_math/ac_reciprocal_pwl_vha.h>
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
//      ac_reciprocal_pwl_vha(input, output);
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
  void ac_reciprocal_pwl_vha(
    const ac_complex<ac_fixed<W, I, S, Q, O> > &input,
    ac_complex<ac_fixed<outW, outI, outS, outQ, outO> > &output
  )
  {
    const int W1 = S ? 2*W - 1 : 2*W + 1;
    const int I1 = S ? 2*I - 1 : 2*I + 1;

    // Start of code to be changed if the ac_fixed PWL implementation changes.

    const bool is_sc_constant_po2 = true; // is sc_constant_lut used in ac_fixed implementation a power of two?
    const int log2_sc_constant = 7; // log2(sc_constant_lut), if sc_constant_lut is a power-of-two. If not, this value is ignored.
    const int n_frac_bits = 16; // Number of fractional bits used in ac_fixed PWL implementation.

    // End of code to be changed if the ac_fixed PWL implementation changes.

    const int sc_input_frac_bits = is_sc_constant_po2 ? AC_MAX(1, AC_MIN(n_frac_bits, W1 - log2_sc_constant - int(S))) : n_frac_bits;

    // The derived type for the reciprocal of the mag_sqr() of the input has its
    // bitwidths calculated to ensure a lossless return type for said reciprocal.
    const int W2 = W1 + sc_input_frac_bits + n_frac_bits;
    const int I2 = W1 - I1;

    ac_fixed<W2, I2, false, outQ, outO> recip_mag_sqr;
    ac_reciprocal_pwl_vha<pwl_Q>(input.mag_sqr(), recip_mag_sqr);

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

    #if !defined(__SYNTHESIS__) && defined(AC_RECIPROCAL_PWL_VHA_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "input.mag_sqr() = " << input.mag_sqr() << std::endl;
    std::cout << "recip_mag_sqr   = " << recip_mag_sqr << std::endl;
    std::cout << "output_temp     = " << output_temp << std::endl;
    std::cout << "output          = " << output << std::endl;
    #endif
  }

//=========================================================================
// Function: ac_reciprocal_pwl_vha (for ac_complex<ac_float>)
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
//    ac_reciprocal_pwl_vha function for ac_float numbers.
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
//    #include <ac_math/ac_reciprocal_pwl_vha.h>
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
//      ac_reciprocal_pwl_vha(input, output);
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
  void ac_reciprocal_pwl_vha(
    const ac_complex<ac_float<W, I, E, Q> > &input,
    ac_complex<ac_float<outW, outI, outE, outQ> > &output
  )
  {
    // Start of code to be changed if the ac_fixed PWL implementation changes.

    const bool is_sc_constant_po2 = true; // is sc_constant_lut used in ac_fixed implementation a power of two?
    const int log2_sc_constant = 7; // log2(sc_constant_lut), if sc_constant_lut is a power-of-two. If not, this value is ignored.
    const int n_frac_bits = 16; // Number of fractional bits used in ac_fixed PWL implementation.

    // End of code to be changed if the ac_fixed PWL implementation changes.

    const int sc_input_frac_bits = is_sc_constant_po2 ? AC_MAX(1, AC_MIN(n_frac_bits, 2*W + 1 - log2_sc_constant)) : n_frac_bits;

    // Calculate real^2 + imag^2
    ac_float<W, I, E, Q> input_real = input.r();
    ac_float<W, I, E, Q> input_imag = input.i();

    ac_float<2*W + 1, 2*I + 1, E + 1, Q> input_mag_sqr;

    // Calculate bitwidths for the return type of reciprocal of the input magnitude.
    const int W_r_m_s = sc_input_frac_bits + n_frac_bits + 2;
    const int I_r_m_s = -(2*I + 1) + 3;
    const int E_r_m_s = E + 3;

    ac_float<W_r_m_s, I_r_m_s, E_r_m_s, outQ> recip_mag_sqr;

    // Store value of input_mag_sqr to the variable.
    input_mag_sqr.add(input_real*input_real, input_imag*input_imag);
    ac_reciprocal_pwl_vha<pwl_Q>(input_mag_sqr, recip_mag_sqr);

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

    #if !defined(__SYNTHESIS__) && defined(AC_RECIPROCAL_PWL_VHA_H_DEBUG)
    std::cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << std::endl;
    std::cout << "input                     = " << input << std::endl;
    std::cout << "input_mag_sqr             = " << input_mag_sqr << std::endl;
    std::cout << "input_mag_sqr.to_double() = " << input_mag_sqr.to_double() << std::endl;
    std::cout << "input_mag_sqr.type_name() = " << input_mag_sqr.type_name() << std::endl;
    std::cout << "recip_mag_sqr             = " << recip_mag_sqr << std::endl;
    std::cout << "recip_mag_sqr.to_double() = " << recip_mag_sqr.to_double() << std::endl;
    std::cout << "recip_mag_sqr.type_name() = " << recip_mag_sqr.type_name() << std::endl;
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
// Function: ac_reciprocal_pwl_vha (for ac_std_float)
//
// Description:
//    Calculation of reciprocal of real inputs, passed as ac_std_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_reciprocal_pwl_vha.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_reciprocal_pwl_vha.h>
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
//      ac_reciprocal_pwl_vha(input, output);
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
  void ac_reciprocal_pwl_vha(
    const ac_std_float<W, E> &input,
    ac_std_float<outW, outE> &output
  )
  {
    // This library was designed to have a relative error tolerance of 0.005%. Through empirical testing,
    // it was found that ac_std_float formats with less than 16 bits allocated for the significand could
    // result in relative error percentages greater than this threshold. Warn the user if they're using a
    // format like that.
    #ifndef __SYNTHESIS__
    static bool print_once = true;
    if (outW - outE < 16 && print_once) {
      std::cout << "Warning: Relative error may exceed 0.005%. Please use an output format that assigns at least 19 bits for the significand." << std::endl;
      print_once = false;
    }
    #endif

    ac_float<outW - outE + 1, 2, outE, AC_RND_CONV> output_ac_fl; // Equivalent ac_float representation for output.
    ac_reciprocal_pwl_vha<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_std_float<outW, outE> output_temp(output_ac_fl); // Convert output ac_float to ac_std_float.
    output = output_temp;
  }

//=========================================================================
// Function: ac_reciprocal_pwl_vha (for ac_ieee_float)
//
// Description:
//    Calculation of reciprocal of real inputs, passed as ac_ieee_float
//    variables.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    // IMPORTANT: ac_std_float.h header file must be included in testbench,
//    // before including ac_reciprocal_pwl_vha.h.
//    #include <ac_std_float.h>
//    #include <ac_math/ac_reciprocal_pwl_vha.h>
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
//      ac_reciprocal_pwl_vha(input, output);
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
  void ac_reciprocal_pwl_vha(
    const ac_ieee_float<Format> &input,
    ac_ieee_float<outFormat> &output
  )
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int outW = T_out::width;
    const int outE = T_out::e_width;
    ac_float<outW - outE + 1, 2, outE, AC_RND_CONV> output_ac_fl; // Equivalent ac_float representation for output.
    ac_reciprocal_pwl_vha<pwl_Q>(input.to_ac_float(), output_ac_fl); // Call ac_float version.
    ac_ieee_float<outFormat> output_temp(output_ac_fl); // Convert output ac_float to ac_ieee_float.
    output = output_temp;
  }
  #endif

//=========================================================================
// Version that allows returning of values.
  template<class T_out, ac_q_mode pwl_Q = AC_TRN, class T_in>
  T_out ac_reciprocal_pwl_vha(
    const T_in &input
  )
  {
    // Create an intermediate variable for output and use the pass-by-reference version
    // to evaluate it. This intermediate variable is returned as the output.
    T_out output;
    ac_reciprocal_pwl_vha<pwl_Q>(input, output);
    return output;
  }
}

#endif
