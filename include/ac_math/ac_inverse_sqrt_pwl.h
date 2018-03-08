/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 1.0                                                 *
 *                                                                        *
 *  Release Date    : Thu Mar  8 11:17:22 PST 2018                        *
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
//  *************************************************************************
//  File : ac_inverse_sqrt_pwl.h

//  Created on: September 17, 2017

//  Author: Sachchidanand Deo

//  Description: Provides piece-wise linear implementations of the
//  square root function for the AC (tm) Datatypes: ac_fixed, ac_float.

//  Usage:
//  A sample testbench and its implementation looks like this:
//
//    #include <ac_fixed.h>
//    #include <ac_complex.h>
//    #include <ac_math/ac_inverse_sqrt_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_fixed<16, 8, false, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<16, 8, false, AC_RND, AC_SAT> output_type;
//
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_inverse_sqrt_pwl(input,output);   //can be replaced by output = ac_sqrt_pwl <output_type> (input);
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
//    This file uses C++ function overloading and templates for various datatype
//    implementations. Note that, if the unsupported datatype is provided,
//    then the static assert will be thrown that indicates that the given
//    datatypes aren't supported.
//
//    This file uses following additional functionality that are implemented as separate header files:
//      1. <ac_math/ac_normalize.h> : This is a normalization function and is listed under includes/ac_math.
//
//  *************************************************************************

#ifndef _INCLUDED_AC_INVERSE_SQRT_PWL_H_
#define _INCLUDED_AC_INVERSE_SQRT_PWL_H_

// Include headers for data types supported by these implementations
#include <ac_fixed.h>
#include <ac_complex.h>

// Include headers for required functions
#include <ac_math/ac_sqrt_pwl.h>

// Non-synthesizable debug header files
#ifndef __SYNTHESIS__
#include <iostream>
using namespace std;
// Uncomment this line to enable debug messages
//#define INVERSE_SQRT_DEBUG
#endif

//Additional header file (ac_normalize as mentioned in notes.)
# include <ac_math/ac_normalize.h>

//=========================================================================
// Function: ac_inverse_sqrt_pwl (for ac_fixed)
//
// Description:
//    Calculation of square root of positive real inputs, passed as ac_fixed
//    variables.
//
//    Passes inputs for normalization to the function in xx_normalize.h,
//    which gives the exponent and output normalized between 0.5 and 1.
//    The normalized value is then subject to the piecewise linear
//    implementation to calculate the reciprocal of square root.
//
//    For simplification square root of exponent is simply computed by
//    dividing the exponent value by 2 and square root of normalized value is
//    computed using piecewise linear mechanism. Inverse is taken care by
//    doing right shift instead of left shift as in case of square root function.
//
// Usage:
//    Please check code snippet from above for usage.
//
// Notes:
//    The PWL implementation utilizes 3 elements, which has a small impact
//    on accuracy. The error is found to be less than 1%.
//
//-------------------------------------------------------------------------

namespace ac_math
{
  template <int W1, int I1, ac_q_mode q1, ac_o_mode o1, int W2, int I2, ac_q_mode q2, ac_o_mode o2>
  void ac_inverse_sqrt_pwl (const ac_fixed <W1,I1, false, q1, o1> &input, ac_fixed <W2, I2, false, q2, o2> &output)
  {
    // npoints are number of points in the PWL(including the extreme points)
    static const unsigned npoints = 9;
    // nsegments are the number of pwl segments (8 segments are used in this implementation)
    static const unsigned nsegments = npoints-1;

    // Declaring square root of 2 as constant, with standard 12 bits of precision.
    // 12 bits of precision is computed for ROM, using max(error)
    static const ac_fixed <13, 1, false, AC_RND> inverseroot2 = 0.707106781373095;
    // normalized_input is basically output of the normalized function that is to be given to the PWL.
    ac_fixed <W1, 0, false, q1, o1> normalized_input;
    int normalized_exp, normalized_exp_temp3; // Temporary variables to store exponents
    // Output variables to store the PWL outputs, q_mode_temp is rounding mode, supplied by user as a template parameter
    ac_fixed <13, 1, false, q1> normalized_output;
    ac_fixed <13, 1, false, q1> normalized_output_temp;

    normalized_exp = ac_math::ac_normalize (input, normalized_input);
    // Since 8 segments are used, maximum bits required for index = log2(8) = 3
    ac_int <3, false> index;
    // input_sc is scaled value of input, it varies from 0-7.9999 and hence, 3 integer bits and 12 bits precision = 15 total bits
    ac_fixed <15, 3, false> input_sc;

    static const ac_fixed <1, 0, false> x_min = 0.5;
    static const ac_fixed <1, 1, false> x_max = 1.0;
    ac_fixed <W2, I2, false, q2, o2> m1;
    // prop_constant is going to be constant at 8*2 = 16 for nsemgnets = 8, hence bitwidths are chosen as 4,4 vary as per requirement
    static const ac_fixed <5, 5, false> prop_constant = (ac_fixed< 4, 4, false>) nsegments/ (x_max - x_min);
    // m and c values, max(absolute_error) of PWL is 5.818*e-04, whose ceil(log2) is computed which returns 11,
    // one bit is added to take into account addition of mx+c, giving 12 bit precision

    // Results obtained from the LUT generator -.05859375 -.058837890625 -.080810546875 -.0806884765625
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
    normalized_exp_temp3 = normalized_exp;
    // Handling of odd exponent case

    normalized_output_temp = normalized_output * inverseroot2;
    // Right shift the exponent by 1 to do divide by 1 operation
    normalized_exp = normalized_exp >> 1;
    m1 = (normalized_exp_temp3 % 2 == 0) ? normalized_output : 	normalized_output_temp;
    ac_fixed <W2, I2, false, q2, o2> output_temp;
    ac_math::ac_shift_right (m1, normalized_exp, output_temp);
    output = output_temp;

#if !defined(__SYNTHESIS__) && defined(INVERSE_SQRT_DEBUG)
    cout << "input to normalization function = " << input << endl;
    cout << "output (fractional of normalization function = " << normalized_input << endl;
    cout << "normalized exp = " << normalized_exp_temp3 << endl;
    cout << "prop_constant = " << prop_constant << endl;
    cout << "input_sc = " << input_sc << endl;
    cout << "index of element chosen from ROM = " << index << endl;
    cout << "normalized_output = " << normalized_output << endl;
    cout << "normalized_output_temp = " << normalized_output_temp << endl;
    cout << "normalized_exp_temp3 = " << normalized_exp_temp3 << endl;
    cout << "m1 = " << m1 << endl;
    cout << "normalized_exp = " << normalized_exp << endl;
    cout << "final output" << output << endl;
#endif
  }

//=========================================================================
// Function: ac_inverse_sqrt_pwl (for ac_complex <ac_fixed>)
//
// Description:
//    Calculation of square root of fixed point complex inputs,
//    passed as ac_complex <ac_fixed> variables.
//
//    Uses following formula to compute inverse square root of
//    complex number (given by a+bi)
//
//    Formula : 1/ sqrt (a+bi) = ac_inverse_sqrt_pwl (a^2+b^2) * sqrt (a-bi)
//
//    Note that, function accepts true as input sign and hence can accept
//    negative/positive real and imaginary signed complex numbers.
//
//    This requires usage of fixed point implementation of ac_inverse_sqrt_pwl
//    and complex fixed point implementation of ac_sqrt_pwl function.
//
// Notes:
//    The PWL implementation of fixed point implementations of
//    ac_inverse_sqrt_pwl and ac_sqrt_pwl utilizes 3 elements, which has a
//    small impact on accuracy. The error is found to be less than 1%.
//
//-------------------------------------------------------------------------

  template <int W1, int I1, ac_q_mode q1, ac_o_mode o1, int W2, int I2, ac_q_mode q2, ac_o_mode o2>
  void ac_inverse_sqrt_pwl (const ac_complex <ac_fixed <W1,I1,true, q1, o1> > &input, ac_complex <ac_fixed <W2, I2, true, q2, o2> > &output)
  {

    const unsigned I = (I1 > 0) ? I1 : 1; // I handles the condition for I1<0, when I1<0, then the answer needs max to max 1 integer bit and all other bits can be allocated to fractional part.
    const unsigned W = (W1 > I1) ? W1 : W1 + I1; // W handles condition for W>I and I>0(handled by above statement).
    ac_complex <ac_fixed <W1 + 1, I1 + 1, true, q1, o1> > input_conj = input.conj();
    // Declare variable to store conjugate of input complex number
    ac_complex <ac_fixed <8*W+1, 2*I+1, true, q1, o1> > sqrt_conj;
    // Calculate square root of conjugate of input
    ac_sqrt_pwl (input_conj, sqrt_conj);
    ac_fixed <W2, I2, false, q2, o2> inverse_sqrt;
    // compute ac_inverse_sqrt_pwl (a^2+b^2)
    ac_inverse_sqrt_pwl (input.mag_sqr(), inverse_sqrt);
    ac_complex <ac_fixed <W2, I2, true, q2, o2> > output_temp;
    // compute final result of above formula
    output_temp.i() = sqrt_conj.i()*inverse_sqrt;
    output_temp.r() = sqrt_conj.r()*inverse_sqrt;
    // One corner case isn't covered by the above formula. This case happens when the real part of the input is negative and the imaginary part is zero.
    // In such a case, the output imaginary part won't have the correct sign. The line below corrects this.
    output_temp.i() = (input.r() < 0 && input.i() == 0) ? (ac_fixed <W2, I2, true, q2, o2>)-output_temp.i() : output_temp.i();
    output = output_temp;

#if !defined(__SYNTHESIS__) && defined(INVERSE_SQRT_DEBUG)
    cout << "W1 = " << W1 << endl;
    cout << "I1 = " << I1 << endl;
    cout << "input = " << input << endl;
    cout << "input_conj = " << input_conj << endl;
    cout << "input.mag_sqr = " << input.mag_sqr() << endl;
    cout << "Square root of conjugate = " << sqrt_conj << endl;
    cout << "Output of inverse square root for ac_fixed = " << inverse_sqrt << endl;
    cout << "Temporary output = " << output_temp << endl;
#endif
  }


//=========================================================================
// Function: ac_inverse_sqrt_pwl (for ac_float)
//
// Description:
//    Calculation of square root of floating point inputs,
//    passed as ac_float variables.
//
//    Separates mentissa and exponent of floating point number.
//    Gives mentissa to the fixed point implementation of
//    ac_inverse_sqrt_pwl, exponent is halved. Based on if input exponent
//    is even or odd, final result is multiplied by inverse of root (2) or
//    not.
//
//    Note that, function accepts only positive real floating point numbers.
//
//    This function uses fixed point implementation of ac_inverse_sqrt_pwl.
//
// Notes:
//    The PWL implementation of fixed point implementations of
//    ac_inverse_sqrt_pwl utilizes 3 elements, which has a
//    small impact on accuracy. The error is found to be less than 1%.
//
//-------------------------------------------------------------------------

  template <int W1, int I1, int E1, ac_q_mode q1, int W2, int I2, int E2, ac_q_mode q2>
  void ac_inverse_sqrt_pwl (const ac_float <W1, I1, E1, q1> &input, ac_float <W2, I2, E2, q2> &output)
  {
    static const ac_fixed <12, 0, false, AC_RND, AC_SAT> inverseroot2 = 0.707106781373095;
    ac_fixed <W2, I2, false, q2> output2;
    ac_fixed <W1, I1, false, q1> m1 = input.m;
    int e1 = input.e;
    int e2 = - (e1 >> 1);
    ac_inverse_sqrt_pwl (m1, output2);
    ac_fixed <W2, I2, false, AC_RND> temp = output2*inverseroot2;
    output.m = (input.e % 2 ==0) ? output2 : temp;
    output.e = e2;

#if !defined(__SYNTHESIS__) && defined(INVERSE_SQRT_DEBUG)
    cout << "Mentissa of input = " << m1 << endl;
    cout << "Exponent of input = " << e1 << endl;
    cout << "Negative Half of input exponent =" << e2 << endl;
    cout << "Output of inverse sqrt function call = " << output2 << endl;
    cout << "Temporary output (result for odd exponent =" << temp << endl;
    cout << "Final output =" << output << endl;
#endif
  }
}

#endif

