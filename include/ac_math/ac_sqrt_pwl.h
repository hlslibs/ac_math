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
//  *************************************************************************
//  File : ac_sqrt_pwl.h

//  Created on: Jun 14, 2017

//  Author: Sachchidanand Deo

//  Description: Provides piece-wise linear implementations of the
//  square root function for the AC (tm) Datatypes: ac_fixed, ac_float,
//  ac_complex<ac_fixed> and ac_complex<ac_float>.

//  Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_fixed.h>
//    #include <ac_math/ac_sqrt_pwl.h>
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
//      ac_sqrt_pwl(input,output);   //can be replaced by output = ac_sqrt_pwl <output_type> (input);
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

#ifndef _INCLUDED_AC_SQRT_PWL_H_
#define _INCLUDED_AC_SQRT_PWL_H_

// Enable this debug macro to print the debug statements
//#define SQRT_DEBUG

// Include headers for data types supported by these implementations
#include <ac_fixed.h>
#include <ac_complex.h>

// Include headers for required functions
#include <ac_math/ac_shift.h>
#include <ac_math/ac_normalize.h>

#if !defined(__SYNTHESIS__)
#include <iostream>
using namespace std;
#endif

//=========================================================================
// Function: ac_sqrt_pwl (for ac_fixed)
//
// Description:
//    Calculation of square root of positive real inputs, passed as ac_fixed
//    variables.
//
//    Passes inputs for normalization to the function in xx_normalize.h,
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
// Notes:
//    The PWL implementation utilizes 3 elements, which has a small impact
//    on accuracy. The error is found to be less than 1%.
//
//-------------------------------------------------------------------------

namespace ac_math
{
  template <ac_q_mode q_mode_temp = AC_RND, int input_width, int input_int, ac_q_mode q_mode, ac_o_mode o_mode, int output_width, int output_int, ac_q_mode q_mode_out, ac_o_mode o_mode_out>
  void ac_sqrt_pwl (const ac_fixed < input_width, input_int, false, q_mode, o_mode> input, ac_fixed <output_width, output_int, false, q_mode_out, o_mode_out> &output)
  {

    // npoints are number of points in the PWL(including the extreme points)
    static const unsigned npoints = 5;
    // nsegments are the number of pwl segments (4 segments are used in this implementation)
    static const unsigned nsegments = npoints-1;

    // Declaring square root of 2 as constant, with standard 12 bits of precision.
    // 12 bits of precision is computed for ROM, using max(error)
    static const ac_fixed <13, 1, false, AC_RND> root2 = 1.414213562;

    // normalized_input is basically output of the normalized function that is to be given to the PWL.
    ac_fixed <input_width, 0, false, q_mode, o_mode> normalized_input;
    // Temporary variables to store exponents
    int normalized_exp,normalized_exp_temp1, normalized_exp_temp2, normalized_exp_temp3;
    // Output variables to store the PWL outputs, q_mode_temp is rounding mode, supplied by user as a template parameter
    ac_fixed <13, 1, false, q_mode_temp> normalized_output;
    ac_fixed <13, 1, false, q_mode_temp> normalized_output_temp;
    // Call to normalization function
    normalized_exp = ac_math::ac_normalize (input, normalized_input);
    // Since 4 segments are used, maximum bits required for index = log2(4) = 2

    ac_int <2, false> index;
    // input_sc is scaled value of input, it varies from 0-3.9999 and hence, 2 integer bits and 12 bits precision = 14 total bitss
    ac_fixed <14, 2, false> input_sc;
    // Piece-wise linear implemenation
    static const ac_fixed <1, 0, false> x_min = 0.5;
    static const ac_fixed <1, 1, false> x_max = 1.0;
    ac_fixed <13, 1, false, q_mode_temp> m1;
    // prop_constant is going to be constant at 4*2 = 8 for nsemgnets = 4, hence bitwidths are chosen as 4,4 vary as per requirement
    static const ac_fixed <4, 4, false> prop_constant = (ac_fixed< 4, 4, false>) nsegments/ (x_max - x_min);
    // m and c values, max(absolute_error) of PWL is 5.818*e-04, whose ceil(log2) is computed which returns 11, one bit is added to take into account addition of mx+c, giving 12 bit precision

    static const ac_fixed <12, 0, false> m[nsegments] = {.08349609375, .0751953125, .0693359375, .064453125};
    static const ac_fixed <12, 0, false> c[nsegments] = {.707763671875, .791259765625, .866455078125, .935791015625};

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
    normalized_output_temp = normalized_output * root2;
    // Right shift the exponent by 1 to do divide by 1 operation
    normalized_exp = normalized_exp >> 1;
    m1 = (normalized_exp_temp3 % 2 == 0) ? normalized_output : 	normalized_output_temp;

    // exponent and normalized output are combined to get the final ac_fixed value, which is written at memory location of output
    ac_math::ac_shift_left (m1, normalized_exp, output);
    output = (input == 0) ? 0 : output;

#if !defined(__SYNTHESIS__) && defined(SQRT_DEBUG)
    cout << "input_width = " << input_width << endl;
    cout << "input_int = " << input_int << endl;
    cout << "output_width = " << output_width << endl;
    cout << "output_int = " << output_int << endl;
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
// Function: ac_sqrt_pwl (for ac_float)
//
// Description:
//    Calculation of square root of positive real inputs, passed as ac_float
//    variables.
//
//    This function uses the piecewise linear implementation by separating the
//    mentissa and exponent. Exponent is simply divided by two, if it is even or
//    is made even and then multiplied by square root of 2, which is stored as
//    as a constant, where as mentissa undergoes piecewise linear implementation
//    using helper function defined above.
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_float.h>
//    #include <ac_math/ac_sqrt_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_float<16, 10, 8, AC_RND> input_type;
//    typedef ac_fioat<18, 20, 9, AC_RND> output_type;
//
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_sqrt_pwl(input,output);   //can be replaced by output = ac_sqrt_pwl <output_type> (input);
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
//    The PWL implementation utilizes 3 elements, which has a small impact
//    on accuracy. The error is found to be less than 1%, which makes the function
//    pretty robust and due to PWL, pretty fast.
//
//-------------------------------------------------------------------------

  template <ac_q_mode q_mode_temp= AC_RND, int input_width, int input_int, int input_exp, ac_q_mode q_mode, int output_width, int output_int, int output_exp, ac_q_mode q_mode_out>
  void ac_sqrt_pwl (const ac_float < input_width, input_int, input_exp, q_mode> input, ac_float <output_width, output_int, output_exp, q_mode_out> &output)
  {
    // root2 is assigned as a constant, precision of it is fixed at 16 bits, not significant area change observed with reduction in bitwidth, plus good accuracy at this bitwidth
    static const ac_fixed <16, 1, false, AC_RND> root2 = 1.414213562;

    // mantissa and exponent are separated
    ac_fixed <input_width,input_int, false, q_mode> mantissa = input.m;
    ac_int <input_exp, true> e = input.e;

    // creating enough space in the output of square root of mantissa, for W<I and I<0 cases
    static const int output_width1 = 2*input_width;
    static const int output_int1 = (input_int > 0) ? input_int : 1;

    // declaring variable to store square root of mentissa
    ac_fixed <output_width1,output_int1,false, q_mode_out>	m2;

    // call to ac_fixed implementation to get square root of mantissa
    ac_sqrt_pwl<q_mode_temp>(mantissa, m2);

    // Multiplication by root 2 for odd exponent
    ac_fixed <output_width1, output_int1,false, q_mode_out>	m3 = m2 * root2;

    ac_fixed<output_width1, output_int1, false, q_mode_out> temp;
    // assign temp variable based on even or odd exponent
    temp = (input.e % 2 == 0) ? m2:m3;

    // final exponent computation based on even or odd input exponent
    output.e = input.e >> 1;
    // assign the resultant mantissa to output mantissa
    output.m = temp;

#if !defined(__SYNTHESIS__) && defined(SQRT_DEBUG)
    cout << "final input = " << input << endl;
    cout << "input_width = " <<input_width << endl;
    cout << "input_int = " << input_int << endl;
    cout << "output of call to ac_fixed version of sqrt_pwl = " << m2 << endl;
    cout << "m3 = " << m3.to_string(AC_BIN).c_str() << endl;
    cout << "temp =" << temp.to_string(AC_BIN).c_str() << endl;
    cout << "output.m =" << output.m.to_string(AC_BIN).c_str() << endl;
    cout << "output.m.nonbits =" << output.m << endl;
    cout << "final output = " << output.to_string(AC_BIN).c_str() << endl;
#endif
  }

//=========================================================================
// Function: ac_sqrt_pwl (for ac_complex <ac_fixed>)
//
// Description:
//    Calculation of square root of positive complex inputs, passed as ac_complex
//    data with real and imaginary part as ac_fixed type of data.
//
//    This function uses following mathematical formula for computation of square root:
//
//    output_real_part = square root of ((input_real_part + square root of mod)/2.0)
//    output_imaginary_part = square root of ((-input_real_part + square root of mod)/2.0)
//    where mod of complex number is given by, input_real_part*input_real_part+input_imaginary_part * input_imaginary_part
//    Then square root is calculated by using piecewise linearly implimentation which uses previously implmeneted function for
//    ac_fixed.
//    After testing the function with testbench of about 400 inputs, it was found out that error is always less than 1%.
//    Two structures are overloaded for making sure that only ac_complex <ac_fixed> type are accepted,
//    both the input and output types are same and to compute internal bitwidths to determine internal datatypes.
//
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_fixed.h>
//    #include <ac_complex.h>
//    #include <ac_math/ac_sqrt_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_float<16, 10, 8, AC_RND> input_type;
//    typedef ac_fioat<18, 20, 9, AC_RND> output_type;
//
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_sqrt_pwl(input,output);   //can be replaced by output = ac_sqrt_pwl <output_type> (input);
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
//    The PWL implementation utilizes 3 elements, which has a small impact
//    on accuracy.
//
//-------------------------------------------------------------------------

  template <ac_q_mode q_mode_temp= AC_RND, int input_width, int input_int, ac_q_mode q_mode_in, ac_o_mode o_mode_in, int output_width, int output_int, ac_q_mode q_mode_out, ac_o_mode o_mode_out>
  void ac_sqrt_pwl (const ac_complex <ac_fixed <input_width, input_int, true, q_mode_in, o_mode_in> > input,  ac_complex <ac_fixed <output_width, output_int, true, q_mode_out, o_mode_out> > &output)
  {
    const unsigned I = (input_int > 0) ? input_int : 1; // temp_I handles the condition for I<0, when I<0, then the answer needs max to max 1 integer bit and all other bits can be allocated to fractional part.
    const unsigned W = (input_width > input_int) ? input_width : input_width + input_int; // temp_W handles condition for W>I and I>0(handled by above statement).

    ac_fixed <4*W+1, 2*I+1, false, q_mode_in, o_mode_in> sqrt_mod;
    ac_sqrt_pwl <q_mode_temp> (input.mag_sqr(), sqrt_mod);     // computation of square root of mag_sqr
    ac_fixed <4*W+1, 2*I+1, false, q_mode_in, AC_SAT> temp_real = sqrt_mod + input.r();
    ac_fixed <4*W+1, 2*I+1, false, q_mode_in, AC_SAT> temp_imag = sqrt_mod - input.r();
    ac_fixed <4*W, 2*I, false, q_mode_in, o_mode_in> sqr_real = temp_real >> 1; // calculating square of output real part
    ac_fixed <4*W, 2*I, false, q_mode_in, o_mode_in> sqr_imag = temp_imag >> 1; // calculating square of output imaginary part
    ac_fixed <8*W, 2*I, false, q_mode_in, o_mode_in> x;
    ac_fixed <8*W, 2*I, false, q_mode_in, o_mode_in> y;
    ac_sqrt_pwl (sqr_real, x); // calculating output real part
    ac_sqrt_pwl (sqr_imag, y); // calculating output imaginary part
    // sign adjustment based on quadrant (for non-positive real and/or imaginary part condition)
    ac_fixed <8*W+1, 2*I+1, true, q_mode_in, o_mode_in> y_or = y;
    ac_fixed <8*W+1, 2*I+1, true, q_mode_in, o_mode_in> y_neg = -y;
    output.r() = x;
    output.i() = (input.i() < 0) ? y_neg : y_or; // if imaginary part is less than zero, assign output value as negative otherwise positive

#if !defined(__SYNTHESIS__) && defined(SQRT_DEBUG)
    cout << "initial input = " << input << endl;
    cout << "input_width = " <<input_width << endl;
    cout << "input_int = " << input_int << endl;
    cout << "output_width = " << output_width << endl;
    cout << "output_int = " << output_int << endl;
    cout << "intermediate width = " << W << endl;
    cout << "intermediate integer width = " << I << endl;
    cout << "Value of square root of mag_sqr = " << sqrt_mod << endl;
    cout << "Result of addition = " << temp_real << endl;
    cout << "Result of subtraction = " << temp_imag << endl;
    cout << "Result of square of output real part = " << sqr_real << endl;
    cout << "Result of square of output imaginary part = " << sqr_imag << endl;
    cout << "Absolute value of real part = " << x << endl;
    cout << "Absolute value of imaginary part = " << y << endl;
    cout << "Final value of square root of complex number = " << output << endl;
#endif
  }

//=========================================================================
// Version that allows returning of values
  template<class T_out, ac_q_mode q_mode_temp = AC_RND, class T_in>
  T_out ac_sqrt_pwl(
    const T_in &input
  )
  {
    // Initializing the final output value that is to be returned
    T_out output;
    // Call the function by referencing the output variable. This is call to one of above implementations
    ac_sqrt_pwl<q_mode_temp>(input, output);
    // Return the final computed output
    return output;
  }
}

#endif

