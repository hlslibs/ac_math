/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 2.0                                                 *
 *                                                                        *
 *  Release Date    : Tue May  1 13:47:52 PDT 2018                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 2.0.2                                               *
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
// **************************************************************************
// File : ac_inverse_sqrt_pwl.h
// 
// Created on: September 17, 2017
// 
// Author: Sachchidanand Deo
// 
// Description: Provides piece-wise linear implementations of the inverse
// square root function for the AC (tm) Datatypes: ac_fixed, ac_float and
// ac_complex<ac_fixed>
// 
// Usage:
//    A sample testbench and its implementation looks like this:
// 
//    #include <ac_math/ac_inverse_sqrt_pwl.h>
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
//      ac_inverse_sqrt_pwl(input,output);
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
//    This file uses C++ function overloading to target implementations specific
//    to each type of data. Attempting to call the function with a type that is
//    not implemented will result in a compile error.
// 
//    This library uses the following functions from other files: ac_normalize()
//    , ac_shift_right() and ac_sqrt_pwl().
// 
// *****************************************************************************

#ifndef _INCLUDED_AC_INVERSE_SQRT_PWL_H_
#define _INCLUDED_AC_INVERSE_SQRT_PWL_H_

#include <ac_int.h>
// Include headers for data types supported by these implementations
#include <ac_fixed.h>
#include <ac_float.h>
#include <ac_complex.h>

// Include headers for required functions
#include <ac_math/ac_normalize.h>
#include <ac_math/ac_shift.h>
#include <ac_math/ac_sqrt_pwl.h>

#if !defined(__SYNTHESIS__) && defined(AC_INVERSE_SQRT_PWL_H_DEBUG)
#include <iostream>
using namespace std;
#endif

// ==============================================================================
// Function: ac_inverse_sqrt_pwl (for ac_fixed)
// 
// Description:
//    Calculation of square root of positive real inputs, passed as ac_fixed
//    variables.
// 
//    Passes inputs for normalization to the function in ac_normalize.h,
//    which gives the exponent and output normalized between 0.5 and 1.
//    The normalized value is then subject to the piecewise linear
//    implementation to calculate the reciprocal of square root.
// 
//    For simplification square root of exponent is simply computed by
//    dividing the exponent value by 2 and square root of normalized value is
//    computed using piecewise linear mechanism.
// 
// ------------------------------------------------------------------------------

namespace ac_math
{
  template<ac_q_mode q_mode_temp = AC_RND, int W1, int I1, ac_q_mode q1, ac_o_mode o1, int W2, int I2, ac_q_mode q2, ac_o_mode o2>
  void ac_inverse_sqrt_pwl(
    const ac_fixed <W1, I1, false, q1, o1> &input,
    ac_fixed <W2, I2, false, q2, o2> &output
  )
  {
    // Use a macro to activate the AC_ASSERT
    // If AC_ASSERT is activated: the program will stop running as soon as a zero input
    // is encountered.
    // If AC_ASSERT is not activated: the output will saturate when a zero input is encountered.
    // The functionality behind this is taken care of by other sections of the code.
#ifdef ASSERT_ON_INVALID_INPUT
    AC_ASSERT(!!input, "Inverse square root of zero not supported.");
#endif

    // npoints are number of points in the PWL(including the extreme points)
    static const unsigned npoints = 9;
    // nsegments are the number of pwl segments (8 segments are used in this implementation)
    static const unsigned nsegments = npoints - 1;

    static const ac_fixed <12, 0, false> inverseroot2 = 0.70703125;
    // normalized_input is basically output of the normalized function that is to be given to the PWL.
    ac_fixed <W1, 0, false, q1, o1> normalized_input;
    int normalized_exp, normalized_exp_temp2; // Temporary variables to store exponents
    // Output variables to store the PWL outputs, q_mode_temp is rounding mode, supplied by user as a template parameter
    ac_fixed <25, 1, false, q_mode_temp> normalized_output;
    ac_fixed <25, 1, false, q_mode_temp> normalized_output_temp;

    normalized_exp = ac_math::ac_normalize (input, normalized_input);
    ac_int <3, false> index;
    // input_sc is scaled value of input, it lies in the range [0, 8)
    ac_fixed <15, 3, false> input_sc;

    static const ac_fixed <1, 0, false> x_min = 0.5;
    static const ac_fixed <1, 1, false> x_max = 1.0;
    ac_fixed <W2, I2, false, q2, o2> m1;
    // proportionality constant that is used for scaling the input
    static const ac_fixed <5, 5, false> prop_constant = (ac_fixed< 4, 4, false>) nsegments/ (x_max - x_min);

    // slope and intercept value array
    static const ac_fixed <12, 0, true> m[nsegments] = {-.080810546875, -.068115234375, -.05859375, -.05126953125, -.045166015625, -.040283203125, -.0361328125, -.03271484375};
    static const ac_fixed <13, 1, false> c[nsegments] = {1.413330078125, 1.33251953125, 1.26416015625, 1.20556640625, 1.154296875, 1.109130859375, 1.06884765625, 1.032470703125};

    // Scaled input is computed from the normalized output value
    input_sc = (normalized_input - x_min) * prop_constant;
    // Take out the fractional bits of the scaled input
    ac_fixed<15 - 3, 0, false> input_sc_frac;
    input_sc_frac.set_slc(0, input_sc.template slc<15 - 3>(0));
    // index is taken as integer part of scaled value and used for selection of m and c values
    index = input_sc.to_int();
    // normalized output provides square root of normalized value
    normalized_output = m[index]*input_sc_frac + c[index];
    // store the initial exponent value in temporary variable
    normalized_exp_temp2 = normalized_exp;
    // Handling of odd exponents
    normalized_output_temp = normalized_output * inverseroot2;
    // Right shift the exponent by 1 to divide by 2
    normalized_exp = normalized_exp >> 1;
    m1 = (normalized_exp_temp2 % 2 == 0) ? normalized_output : 	normalized_output_temp;
    ac_fixed <W2, I2, false, q2, o2> output_temp;
    // "De-normalize" the output by performing a right-shift and cancel out the effects of the previous normalization.
    ac_math::ac_shift_right (m1, normalized_exp, output_temp);

    // If a zero input is encountered, the output must saturate regardless of whether the assert has been activated or not.
    // Assign a variable that stores the saturated output value.
    ac_fixed<W2, I2, false, q2, o2> output_temp_max;
    output_temp_max.template set_val<AC_VAL_MAX>();
    // Use a ternary operator to decide whether the output should store the PWL-calculated value or the saturated value, based
    // on whether a zero was passed or not.
    output = input != 0 ? output_temp : output_temp_max;

#if !defined(__SYNTHESIS__) && defined(AC_INVERSE_SQRT_PWL_H_DEBUG)
    cout << "input                  = " << input << endl;
    cout << "normalized_input       = " << normalized_input << endl;
    cout << "normalized_exp         = " << normalized_exp << endl;
    cout << "prop_constant          = " << prop_constant << endl;
    cout << "input_sc               = " << input_sc << endl;
    cout << "index                  = " << index << endl;
    cout << "normalized_output      = " << normalized_output << endl;
    cout << "normalized_output_temp = " << normalized_output_temp << endl;
    cout << "normalized_exp_temp2   = " << normalized_exp_temp2 << endl;
    cout << "m1                     = " << m1 << endl;
    cout << "normalized_exp         = " << normalized_exp << endl;
    cout << "output                 = " << output << endl;
#endif
  }

// ===========================================================================
// Function: ac_inverse_sqrt_pwl (for ac_complex <ac_fixed>)
// 
// Description:
//    Calculation of square root of fixed point complex inputs,
//    passed as ac_complex <ac_fixed> variables.
// 
//    Uses following formula to compute inverse square root of
//    complex number (given by a+bi)
// 
//    Formula : 1/ sqrt (a+bi) = inverse_sqrt(a^2+b^2) * sqrt (a-bi)
// 
//    Note that, function accepts true as input sign and hence can accept
//    negative/positive real and imaginary signed complex numbers.
// 
//    This requires usage of fixed point implementation of ac_inverse_sqrt_pwl
//    and complex fixed point implementation of ac_sqrt_pwl function.
// 
// ---------------------------------------------------------------------------

  template <ac_q_mode q_mode_temp = AC_RND, int W1, int I1, ac_q_mode q1, ac_o_mode o1, int W2, int I2, ac_q_mode q2, ac_o_mode o2>
  void ac_inverse_sqrt_pwl (const ac_complex <ac_fixed <W1,I1,true, q1, o1> > &input, ac_complex <ac_fixed <W2, I2, true, q2, o2> > &output)
  {
    const unsigned I = (I1 > 0) ? I1 : 1; // I handles the condition for I1<0, when I1<0, then the answer needs max to max 1 integer bit and all other bits can be allocated to fractional part.
    const unsigned W = (W1 > I1) ? W1 : W1 + I1; // W handles condition for W>I and I>0(handled by above statement).
    ac_complex <ac_fixed <W1 + 1, I1 + 1, true, q1, o1> > input_conj = input.conj();
    // Declare variable to store conjugate of input complex number
    ac_complex <ac_fixed <8*W+1, 2*I+1, true, q1, o1> > sqrt_conj;
    // Calculate square root of conjugate of input
    ac_sqrt_pwl<q_mode_temp> (input_conj, sqrt_conj);
    ac_fixed <W2, I2, false, q2, o2> inverse_sqrt;
    // compute ac_inverse_sqrt_pwl (a^2+b^2)
    ac_inverse_sqrt_pwl<q_mode_temp> (input.mag_sqr(), inverse_sqrt);
    ac_complex <ac_fixed <W2, I2, true, q2, o2> > output_temp;
    // compute final result
    output_temp.i() = sqrt_conj.i()*inverse_sqrt;
    output_temp.r() = sqrt_conj.r()*inverse_sqrt;
    // One corner case isn't covered by the above formula. This case happens when the real part of the input is negative and the imaginary part is zero.
    // In such a case, the output imaginary part won't have the correct sign. The line below corrects this.
    output_temp.i() = (input.r() < 0 && input.i() == 0) ? (ac_fixed <W2, I2, true, q2, o2>)-output_temp.i() : output_temp.i();
    // If a zero input is encountered, the real part of the output must saturate regardless of whether the assert has been activated or not.
    // Assign a variable that stores the saturated output value.
    ac_fixed<W2, I2, true, q2, o2> output_temp_max;
    output_temp_max.template set_val<AC_VAL_MAX>();
    bool non_zero_input = input.mag_sqr() != 0;
    // Use a ternary operator to decide whether the output real part should store the PWL-calculated value or the saturated value, based
    // on whether a zero was passed or not at the input.
    output.r() = non_zero_input ? output_temp.r() : output_temp_max; 
    // Use a ternary operator to decide whether the output imaginary part should store the PWL-calculated value or a zero value, based
    // on whether a zero was passed or not at the input.
    output.i() = non_zero_input ? output_temp.i() : 0;

#if !defined(__SYNTHESIS__) && defined(AC_INVERSE_SQRT_PWL_H_DEBUG)
    cout << "W1              = " << W1 << endl;
    cout << "I1              = " << I1 << endl;
    cout << "input           = " << input << endl;
    cout << "input_conj      = " << input_conj << endl;
    cout << "input.mag_sqr   = " << input.mag_sqr() << endl;
    cout << "sqrt_conj       = " << sqrt_conj << endl;
    cout << "inverse_sqrt    = " << inverse_sqrt << endl;
    cout << "output_temp_max = " << output_temp_max << endl;
    cout << "output_temp     = " << output_temp << endl;
    cout << "output          = " << output << endl;
#endif
  }


// =============================================================================
// Function: ac_inverse_sqrt_pwl (for ac_float)
// 
// Description:
//    Calculation of square root of floating point inputs,
//    passed as ac_float variables.
// 
//    Separates mantissa and exponent of floating point number.
//    Gives mantissa to the fixed point implementation of
//    ac_inverse_sqrt_pwl, and halves the exponent. Based on if input exponent
//    is even or odd, final result is multiplied by inverse of root (2) or
//    not.
// 
//    The function accepts only positive real floating point numbers.
// 
//    This function uses the fixed point implementation of ac_inverse_sqrt_pwl.
// 
// ----------------------------------------------------------------------------

  template <ac_q_mode q_mode_temp = AC_RND, int W1, int I1, int E1, ac_q_mode q1, int W2, int I2, int E2, ac_q_mode q2>
  void ac_inverse_sqrt_pwl (const ac_float <W1, I1, E1, q1> &input, ac_float <W2, I2, E2, q2> &output)
  {
    static const ac_fixed <12, 0, false> inverseroot2 = 0.70703125;
    ac_fixed <W2, I2, false, q2> output2;
    ac_fixed <W1, I1, false, q1> m1 = input.m;
    int e1 = input.e;
    int e2 = - (e1 >> 1);
    ac_inverse_sqrt_pwl<q_mode_temp> (m1, output2);
    ac_fixed <W2, I2, false, q2> temp = output2*inverseroot2;
    output.m = (input.e % 2 == 0) ? output2 : temp;
    output.e = e2;

#if !defined(__SYNTHESIS__) && defined(AC_INVERSE_SQRT_PWL_H_DEBUG)
    cout << "m1      = " << m1 << endl;
    cout << "e1      = " << e1 << endl;
    cout << "e2      = " << e2 << endl;
    cout << "output2 = " << output2 << endl;
    cout << "temp    = " << temp << endl;
    cout << "output  = " << output << endl;
#endif
  }

// Function definition to enable return by value.

// =========================================================================
// Version that allows returning of values
  template<class T_out, ac_q_mode q_mode_temp = AC_RND, class T_in>
  T_out ac_inverse_sqrt_pwl(
    const T_in &input
  )
  {
    // Initializing the final output value that is to be returned
    T_out output;
    // Call the function by referencing the output variable. This is call to one of above implementations
    ac_inverse_sqrt_pwl<q_mode_temp>(input, output);
    // Return the final computed output
    return output;
  }

}

#endif // _INCLUDED_AC_INVERSE_SQRT_PWL_H_

