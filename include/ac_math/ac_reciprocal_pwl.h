/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 1.0                                                 *
 *                                                                        *
 *  Release Date    : Thu Mar  1 16:35:45 PST 2018                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 1.0.0                                               *
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
//*****************************************************************************************
// File: ac_reciprocal_pwl.h
//
// Description: Provides piece-wise linear implementations of the
//    reciprocal function for the AC (tm) Datatypes: ac_fixed, ac_float,
//    ac_complex<ac_fixed> and ac_complex<ac_float>.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_fixed.h>
//    #include <ac_math/ac_reciprocal_pwl.h>
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
//      ac_reciprocal_pwl(input, output);
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
// Notes:
//    This file uses C++ function overloading to target implementations
//    specific to each type of data. Attempting to call the function
//    with a type that is not implemented will result in a compile error.
//
//    This file uses the normalization function from ac_normalize.h.
//
// Revision History:
//    Niramay Sanghvi : Aug 10 2017 : Added default parameters for better configurability.
//    Niramay Sanghvi : Aug 09 2017 : Added zero-input handling
//    Niramay Sanghvi : Aug 07 2017 : Used right-shift function from mgc_ac_math
//    Niramay Sanghvi : Jul 27 2017 : Added structs for checking input and output types.
//    Niramay Sanghvi : Jul 06 2017 : Updated header style.
//    Niramay Sanghvi : Jul 05 2017 : Passed output by reference.
//    Niramay Sanghvi : Jun 21 2017 : Made LUT precision configurable.
//    Niramay Sanghvi : Jun 20 2017 : Added header style format.
//    Niramay Sanghvi : Jun 19 2017 : Added support for ac_complex.
//
//*****************************************************************************************

#ifndef _INCLUDED_AC_RECIPROCAL_PWL_H_
#define _INCLUDED_AC_RECIPROCAL_PWL_H_

// The functions use default template parameters, which are only supported by C++11 or later
// compiler standards. Hence, the user should be informed if they are not using those standards.

#if !(__cplusplus >= 201103L)
#error Please use C++11 or a later standard for compilation.
#endif

// Include headers for data types supported by these implementations
#include <ac_int.h>
#include <ac_float.h>
#include <ac_fixed.h>
#include <ac_complex.h>

// Include headers for required functions
#include <ac_math/ac_normalize.h>
#include <ac_math/ac_shift.h>

#if !defined(__SYNTHESIS__) && defined(AC_RECIPROCAL_PWL_DEBUG)
#include <iostream>
using namespace std;
#endif

//=========================================================================
// Function: ac_reciprocal_pwl (for ac_fixed)
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
//    on accuracy. For better accuracy, please use the LUT_Generator.cpp
//    file to get new values for PWL LUTs which utilize more than 7
//    elements.
//
//-------------------------------------------------------------------------

namespace ac_math
{
  template<ac_q_mode pwl_Q = AC_RND,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO>
  void ac_reciprocal_pwl(
    const ac_fixed<W, I, S, Q, O> &input,
    ac_fixed<outW, outI, outS, outQ, outO> &output
  )
  {
    // Use a macro to activate or de-activate the AC_ASSERT
    // If AC_ASSERT is activated: the program will stop running as soon as a zero input
    // is encountered.
    // If AC_ASSERT is not activated: the output will saturate when a zero input is encountered.
    // The functionality behind this is taken care of by other sections of the code.
#ifdef AC_RECIPROCAL_PWL_ASSERT
    AC_ASSERT(!!input, "Reciprocal of zero not supported.");
#endif

    ac_fixed<outW, outI, outS, outQ, outO> output_temp;

    // Initialization for PWL LUT
    static const unsigned n_segments_lut = 6;
    static const ac_fixed<9, 0, true> m_lut[n_segments_lut] = {-.28515625, -.212890625, -.166015625, -.130859375, -.109375, -.08984375};
    static const ac_fixed<10, 1, false> c_lut[n_segments_lut] = {1.994140625, 1.708984375, 1.49609375, 1.330078125, 1.19921875, 1.08984375};
    // Domain of PWL
    static const ac_fixed<1, 0, false> x_min_lut = 0.5;
    static const ac_fixed<1, 1, false> x_max_lut = 1;
    // Scaling constant used later to scale the normalized input from 0 to n_segments_lut
    static const ac_fixed<4, 4, false> sc_constant_lut = n_segments_lut/(x_max_lut - x_min_lut);

    // The absolute value of the input is taken and passed to the normalization function. Initialize variables for the same.
    ac_fixed<W, I, false> input_abs_value;
    ac_fixed<W, 0, false> normalized_fixed;

    // If input is signed, take absolute value and assign to temporary variable.
    if (S) {input_abs_value = ((input >= 0) ? (ac_fixed <W, I, false>)input : (ac_fixed <W, I, false>)(-input));}
    // If input is unsigned, assign value of input to temp. variable.
    else {input_abs_value = input;}

    // Normalize the absolute value. expret stores the value of the returned base 2 exponential.
    int expret_temp = ac_math::ac_normalize(input_abs_value, normalized_fixed);

    // Compute reciprocal using pwl.

    // Scale the normalized input from 0 to n_segments_lut
    ac_fixed<12, 3, true> x_in_sc = (normalized_fixed - x_min_lut)*sc_constant_lut;
    // Take out the fractional bits of the scaled input
    ac_fixed<12 - 3, 0, false> x_in_sc_frac;
    x_in_sc_frac.set_slc(0, x_in_sc.template slc<12 - 3>(0));
    ac_int<3, false> index;
    // The integer part of the input is the index of the LUT table
    index = x_in_sc.to_int();
    // The output of the PWL approximation should have the same signedness as the output of the function.
    // Define the type in such a way that the signedness is the same without affecting precision.
    typedef ac_fixed<20 + int(outS), 2 + int(outS), outS, pwl_Q> output_pwl_type;
    output_pwl_type output_pwl = m_lut[index]*x_in_sc_frac + c_lut[index];

    if (input != 0) { // If input is non-zero, De-normalize output by shifting right by expret_temp
      // If input and output are signed, change sign of output_pwl based on whether input is positive or negative.
      if (S && outS) {output_pwl = (input < 0) ? (output_pwl_type)(-output_pwl) : output_pwl;}
      // ac_shift_right function used for denormalization so as to ensure saturation and rounding.
      ac_math::ac_shift_right(output_pwl, expret_temp, output_temp);
    } else {
      // If zero input is encountered, set output to the max possible value.
      output_temp.template set_val<AC_VAL_MAX>();
    }

    output = output_temp;

#if !defined(__SYNTHESIS__) && defined(AC_RECIPROCAL_PWL_DEBUG)
    cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl;
    cout << "input                     = " << input << endl;
    cout << "input_abs_value           = " << input_abs_value << endl;
    cout << "normalized input          = " << normalized_fixed << endl;
    cout << "expret_temp               = " << expret_temp << endl;
    cout << "x_in_sc                   = " << x_in_sc << endl;
    cout << "x_in_sc_frac              = " << x_in_sc_frac << endl;
    cout << "output_temp               = " << output_temp << endl;
    cout << "output up-scaled by exp   = " << output << endl;
#endif
  }

//=========================================================================
// Function: ac_reciprocal_pwl (for ac_float)
//
// Description:
//    Calculation of reciprocal of real inputs, passed as ac_float
//    variables.
//
//    The mantissa of the ac_float number is passed as an ac_fixed variable
//    to the ac_reciprocal_pwl function for ac_fixed numbers. The mantissa of
//    the output is then set to the resultant reciprocal. The exponent
//    of the output ac_float variable is set to the negative of the
//    exponent of the input.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_float.h>
//    #include <ac_math/ac_reciprocal_pwl.h>
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
//      ac_reciprocal_pwl(input, output);
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

  template<ac_q_mode pwl_Q = AC_RND,
           int W, int I, int E, ac_q_mode Q,
           int outW, int outI, int outE, ac_q_mode outQ>
  void ac_reciprocal_pwl(
    const ac_float<W, I, E, Q> &input,
    ac_float<outW, outI, outE, outQ> &output
  )
  {
    ac_float<outW, outI, outE, outQ> output_temp;

    // Find the reciprocal of the mantissa using the ac_fixed implementation.
    ac_fixed<outW, outI, true, outQ> recip_mantissa;
    ac_reciprocal_pwl<pwl_Q>(input.mantissa(), recip_mantissa);
    output_temp.m = recip_mantissa;

    // If the input is non-zero, assign the additive inverse of the input's
    // exponent to the output's exponent. If the input is zero, set the
    // output exponent to the max. possible value.

    // The output mantissa will have already been set to the max possible value
    // by the ac_fixed version.
    if (input.m != 0) {output_temp.e = -input.exp();}
    else {output_temp.e.template set_val<AC_VAL_MAX>();}

    output = output_temp;

#if !defined(__SYNTHESIS__) && defined(AC_RECIPROCAL_PWL_DEBUG)
    cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl;
    cout << "input             = " << input << endl;
    cout << "input.mantissa()  = " << input.mantissa() << endl;
    cout << "input.exp()       = " << input.exp() << endl;
    cout << "recip_mantissa    = " << recip_mantissa << endl;
    cout << "output_temp       = " << output_temp << endl;
    cout << "output            = " << output << endl;
    cout << "output.mantissa() = " << output.mantissa() << endl;
    cout << "output.exp()      = " << output.exp() << endl;
#endif
  }

//=========================================================================
// Function: ac_reciprocal_pwl (for ac_complex<ac_fixed>)
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
//    ac_reciprocal_pwl function for ac_fixed numbers.
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
//    #include <ac_fixed.h>
//    #include <ac_complex.h>
//    #include <ac_math/ac_reciprocal_pwl.h>
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
//      ac_reciprocal_pwl(input, output);
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

  template<ac_q_mode pwl_Q = AC_RND,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO>
  void ac_reciprocal_pwl(
    const ac_complex<ac_fixed<W, I, S, Q, O> > &input,
    ac_complex<ac_fixed<outW, outI, outS, outQ, outO> > &output
  )
  {
    ac_complex<ac_fixed<outW, outI, outS, outQ, outO> > output_temp;

    ac_fixed<outW, outI, false, outQ, outO> recip_mag_sqr;
    ac_reciprocal_pwl<pwl_Q>(input.mag_sqr(), recip_mag_sqr);

    if (input.mag_sqr() != 0) {
      // Use the formula "1/(a+bi) = (a-bi)/(a^2+b^2)" to assign values to the output.
      output_temp.r() =  input.r() * recip_mag_sqr;
      output_temp.i() = -input.i() * recip_mag_sqr;
    } else {
      // If zero input is passed, then recip_mag_sqr is already set to the max
      // possible value by the ac_fixed implementation. Assign its value to both the
      // real and imaginary parts of the output.
      output_temp.r() = recip_mag_sqr;
      output_temp.i() = recip_mag_sqr;
    }

    output = output_temp;

#if !defined(__SYNTHESIS__) && defined(AC_RECIPROCAL_PWL_DEBUG)
    cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl;
    cout << "input.mag_sqr() = " << input.mag_sqr() << endl;
    cout << "recip_mag_sqr   = " << recip_mag_sqr << endl;
    cout << "output_temp     = " << output_temp << endl;
    cout << "output          = " << output << endl;
#endif
  }

//=========================================================================
// Function: ac_reciprocal_pwl (for ac_complex<ac_float>)
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
//    ac_reciprocal_pwl function for ac_float numbers.
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
//    #include <ac_float.h>
//    #include <ac_complex.h>
//    #include <ac_math/ac_reciprocal_pwl.h>
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
//      ac_reciprocal_pwl(input, output);
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

  template<ac_q_mode pwl_Q = AC_RND,
           int W, int I, int E, ac_q_mode Q,
           int outW, int outI, int outE, ac_q_mode outQ>
  void ac_reciprocal_pwl(
    const ac_complex<ac_float<W, I, E, Q> > &input,
    ac_complex<ac_float<outW, outI, outE, outQ> > &output
  )
  {
    ac_complex<ac_float<outW, outI, outE, outQ> > output_temp;

    // Calculate real^2 + imag^2
    ac_float<W, I, E, Q> input_real = input.r();
    ac_float<W, I, E, Q> input_imag = input.i();

    // Define type for input_mag_sqr.
    // E + 2 is chosen for the exponent width of input_mag_sqr
    // Because it works well with all input exponential widths,
    // as well as corner cases with very low exponential widths (as low as 1 or 2)
    typedef ac_float<2*W + 1, 2*I + 1, E + 2, Q> i_m_s_type;
    i_m_s_type input_mag_sqr;

    ac_float<outW, outI, E + 3, outQ> recip_mag_sqr;

    // Store value of input_mag_sqr to the variable.
    input_mag_sqr = (i_m_s_type)(input_real*input_real);
    input_mag_sqr += (i_m_s_type)(input_imag*input_imag);
    ac_reciprocal_pwl<pwl_Q>(input_mag_sqr, recip_mag_sqr);

    if (input_mag_sqr.mantissa() != 0) {
      // 1/(a+bi) = (a-bi)/(a^2+b^2)
      output_temp.r() =  input.r() * recip_mag_sqr;
      output_temp.i() = -input.i() * recip_mag_sqr;
    } else {
      // If zero input is passed, then recip_mag_sqr is already set to the max
      // possible value by the ac_float implementation.
      output_temp.r() = recip_mag_sqr;
      output_temp.i() = recip_mag_sqr;
    }

    output = output_temp;

#if !defined(__SYNTHESIS__) && defined(AC_RECIPROCAL_PWL_DEBUG)
    cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl;
    cout << "input         = " << input << endl;
    cout << "input_mag_sqr = " << input_mag_sqr << endl;
    cout << "recip_mag_sqr = " << recip_mag_sqr << endl;
    cout << "output_temp   = " << output_temp << endl;
    cout << "output        = " << output << endl;
#endif
  }

//=========================================================================
// Version that allows returning of values.
  template<class T_out, ac_q_mode pwl_Q = AC_RND, class T_in>
  T_out ac_reciprocal_pwl(
    const T_in &input
  )
  {
    // Create a temporary variable for output and use the pass-by-reference version
    // to evaluate it. This temporary variable is returned as the output.
    T_out output;
    ac_reciprocal_pwl<pwl_Q>(input, output);
    return output;
  }
}

#endif

