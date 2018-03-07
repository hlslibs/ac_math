/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 1.0                                                 *
 *                                                                        *
 *  Release Date    : Wed Mar  7 13:09:26 PST 2018                        *
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
//*************************************************************************
// File: ac_log_pwl.h
//
// Description: Provides piece-wise linear implementations of the
//    log function for the AC (tm) Datatypes: ac_fixed.
//    Two different functions compute values with bases as 2 and e.
//    'ifdef' directive can be used to choose between these two functions.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_fixed.h>
//    #include <ac_complex.h>
//    #include <ac_math/ac_log_pwl.h>
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
//      #ifdef TEST_LOG2
//      ac_log2_pwl(input,output);
//      #endif
//      #ifdef TEST_LOGE
//      ac_log_pwl(input, output)
//      #endif
//    }
//
//    #ifndef __SYNTHESIS__
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 1.2;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
//    For a more detailed example, please refer to the "test_ac_log_tb.cpp", "test_ac_log.cpp" and
//    "test_ac_log.h" files
//
// Notes:
//    This file contains two functions:
//     1. ac_log2_pwl : It takes ac_fixed type value as input and writes log2 value of it
//     2. ac_log_pwl : It takes ac_fixed type value as input and writes ln value of it
//    This file uses the normalization function from xx_normalize.h.
//
//
//*************************************************************************

#ifndef _INCLUDED_AC_LOG_PWL_H_
#define _INCLUDED_AC_LOG_PWL_H_

//#define AC_LOG_PWL_DEBUG //uncomment this to get access the Debug statements

// Include headers for data types supported by these implementations
#include <ac_fixed.h>
#include <ac_int.h>

// Include headers for required functions
#include <ac_math/ac_normalize.h>

#if !defined(__SYNTHESIS__) && defined(AC_LOG_PWL_DEBUG)
#include <iostream>
using namespace std;
#endif

//=========================================================================
// Function: ac_log2_pwl (for ac_fixed, returns log base 2 value of input provided)
//
// Description:
//    Calculation of log base 2 of real inputs, passed as ac_fixed
//    variables.
//
//    Passes inputs for normalization to the function in xx_normalize.h,
//    which gives the exponent and output normalized between 0.5 and 1.
//    The normalized output is then subjected to piecewise linear implementation of
//    log2 value and returns log2 of normalized value.
//    Since log2 of exponent, is simply the exponent itself. The value is added to
//    the pwl output and is then passed to the calling function, via pass by reference.
//
//
// Usage:
//    See above example code for usage.
//
// Notes:
//    The PWL implementation utilizes 3 elements, which has a small impact
//    on accuracy.
//
//-------------------------------------------------------------------------

namespace ac_math
{
  template <ac_q_mode q_mode_temp = AC_RND, int W1, int I1, ac_q_mode q_mode_in, ac_o_mode o_mode_in, int W2, int I2, ac_q_mode q_mode_out, ac_o_mode o_mode_out>
  void ac_log2_pwl (const ac_fixed <W1, I1, false, q_mode_in, o_mode_in> input, ac_fixed <W2, I2, true, q_mode_out, o_mode_out> &result)
  {
    // input_normalized is used to store the normalized output and is between 0.5-1
    ac_fixed <W1, 0, false, q_mode_in, o_mode_in> input_normalized;
    // exp is used to store the final exponent value
    int exp;
    // call to the ac_normalize function which gives normalized output(input_normalized) alongside the normalized_exponent value
    exp = ac_math::ac_normalize (input, input_normalized);
    // Piecewise linear implementation
    // Define lower and upper limits
    static const ac_fixed <1,0, false> x_min = 0.5;
    static const ac_fixed <1,1, false> x_max = 1.0;
    // input_sc is used to store scaled input, which is always between 0 to nsegments
    ac_fixed <14,3, false> input_sc;
    // PWL constants. Note that these are variables based on the number of segments. User can change this based on his/her precision requirement
    // nsgements = 9 gave really good performance
    static const unsigned npoints = 9;
    static const unsigned nsegments = npoints-1;
    // This variable is used to store the index which is used to extract value out of the ROM. 3 bits is basically, log2(nsegments), replace based on your value of nsegments
    ac_int <3, false> index;
    // prop_constant is basically used to store the integer constant, given by nsegments*2 (for x_min= 0.5, x_min=1.0), hence used bits given by log2(2*nsegments).
    // In this case, log2(ceil(2*8)) = log2(ceil(16)) = 5
    static const ac_fixed <5,5, false> prop_constant = (ac_fixed <5,5, false, AC_RND, AC_SAT>(nsegments))/ (x_max - x_min);
    // m and c ROM values which are configurable: see LUT generator file for more details
    // 11 bit precision logic:
    // After running pwl in MATLAB, it was observed that maximum absolute error of PWL in range 0.5-1 is 0.00125113.
    // Log2(max_error) is -9.64259, ceil of absolute of this value is 10.
    // Taking into consideration the additional 1 bit due to mx+c addition, we get total number of fractional bit requirement = 11, which is used in following m and c precision declaration.
    // c requires one integer bit due to the signed nature of c
    static const ac_fixed <11,0, false> m[nsegments] = {0.169921875, 0.15185546875, 0.13720703125, 0.12548828125, 0.115234375, 0.10693359375, 0.099609375, 0.09326171875};
    static const ac_fixed <12,1, true> c[nsegments] = {-0.99853515625, -0.98046875, -0.95166015625, -0.916015625, -0.8759765625, -0.8330078125, -0.7890625, -0.744140625};
    // compute the scaled input of normalized input
    input_sc = (input_normalized - x_min) * prop_constant;
    // Integer part of scaled input is index
    index = input_sc.to_int();
    // If 0 is supplied as the function input, maximum negative value is returned at the output, given by ac_val_min
    result.template set_val<AC_VAL_MIN>();
    // computation of the pwl output
    ac_fixed <W2, I2, true, q_mode_temp> t = m[index]*input_sc + c[index];
    // Add the exponent to get the final function output
    ac_fixed <W2, I2, true, q_mode_out, o_mode_out> t2 = t+exp;
    // assignment to the final output
    result = (input == 0) ? result : t2;

#if !defined(__SYNTHESIS__) && defined(AC_LOG_PWL_DEBUG)
    cout << __FILE__ << __LINE__ << endl;
    cout << "input_width" << input_width << endl;
    cout << "input_int" << input_int << endl;
    cout << "input = " << input << endl;
    cout << "input to normalization function" << input << endl;
    cout << "output (fractional of normalization function" << input_normalized << endl;
    cout << "normalized exp" << exp << endl;
    cout << "index of element chosen from ROM" << index << endl;
    cout << "final output = " << result << endl;
#endif
  }

//=========================================================================
// Function: ac_log_pwl (for ac_fixed, returns log base e value of input provided)
//
// Description:
//    Calculation of ln (natural log) of real inputs, passed as ac_fixed
//    variables. e = 2.71828
//
//
//    This implementation uses change of base method to compute log base e of
//    input. The value of ln of input is given by,
//    ln(x) = log2(x)/log2(e)
//    1/log2(e) is declared as static constant whose value is declared internally.
//    which makes, ln(x) = log2(x)*constant
//    log2(x) is computed using above piecewise linear implementation.
//    Note that, 1/log2(e) = loge(2), which is declared as constant
//
//
// Usage:
//    See above example code for usage.
//
// Notes:
//    The PWL implementation of log2 utilizes 3 elements, which has a small impact
//    on accuracy. Separate PWL implementation is also possible loge function.
//
//
//-------------------------------------------------------------------------
  template <ac_q_mode q_mode_temp = AC_RND, int W1, int I1, ac_q_mode q_mode_in, ac_o_mode o_mode_in, int W2, int I2, ac_q_mode q_mode_out, ac_o_mode o_mode_out>
  void ac_log_pwl (const ac_fixed <W1, I1, false, q_mode_in, o_mode_in> input, ac_fixed <W2, I2, true, q_mode_out, o_mode_out> &result)
  {
    // 12 bits precision is given by 11 bits precision of pwl + 1 bit precision of addition in mx + c
    static const ac_fixed <12, 0, false, AC_RND> log_constant =  0.69314718056;
    // call to the log base 2 pwl function
    ac_log2_pwl <q_mode_temp> (input, result);
    // multiply the output of log2(input) by 1/log2(e) to get the final result, which is, loge(input): change of base formula.
    // Although, 1/log2(e) = loge(2) which is defined as above log constant
    result = result * log_constant;

#if !defined(__SYNTHESIS__) && defined (AC_LOG_PWL_DEBUG)
    cout << "Input to the log base e function = " << input << endl;
    cout << "constant_width = " << constant_width << endl;
    cout << "Output of the log base 2 call = " << result1 << endl;
    cout << "Final output =" << result << endl;
#endif
  }

//=========================================================================
// Version that allows returning of values for log2. For usage pass output_type as template parameter
  template<class T_out, ac_q_mode q_mode_temp = AC_RND, class T_in>
  T_out ac_log2_pwl(
    const T_in &input
  )
  {
    // create a variable that is to be returned
    T_out output;
    // call above implementation of log base 2
    ac_log2_pwl<q_mode_temp>(input, output);
    // return the final output
    return output;
  }

//=========================================================================
// Version that allows returning of values for ln. For usage pass output_type as template parameter
  template<class T_out, ac_q_mode q_mode_temp = AC_RND, class T_in>
  T_out ac_log_pwl(
    const T_in &input
  )
  {
    // create a variable that is to be returned
    T_out output;
    // call above implementation of log base e/ln
    ac_log_pwl<q_mode_temp>(input, output);
    // return the final output
    return output;
  }
}
#endif

