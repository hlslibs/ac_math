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
//******************************************************************************************
// Function: ac_tan_pwl (for ac_fixed)
//
// Description:
//    Calculation of tangent of real inputs, passed as ac_fixed
//    variables.
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
//    This file uses the ac_reciprocal_pwl() function from ac_reciprocal_pwl.h
//
//******************************************************************************************

#ifndef _INCLUDED_AC_TAN_PWL_H_
#define _INCLUDED_AC_TAN_PWL_H_

// The functions use default template parameters, which are only supported by C++11 or later
// compiler standards. Hence, the user should be informed if they are not using those standards.

#if !(__cplusplus >= 201103L)
#error Please use C++11 or a later standard for compilation.
#endif

#include <ac_int.h>
// Include headers for data types supported by these implementations
#include <ac_fixed.h>

#if !defined(__SYNTHESIS__) && defined(AC_TAN_PWL_H_DEBUG)
#include <iostream>
using namespace std;
#endif

#include <ac_math/ac_reciprocal_pwl.h>

namespace ac_math
{
  template<int W, int I, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  void ac_tan_pwl(
    const ac_fixed<W, I, false, Q, O> &input,
    ac_fixed<outW, outI, false, outQ, outO> &output
  )
  {
    // Store the approximate value of 89.46 degrees in radian
    static const ac_fixed<11, 1, false> sat_limit = 1.5615234375;
    ac_fixed<outW, outI, false, outQ, outO> output_temp;

    // Initialization for PWL LUT
    static const unsigned n_segments_lut = 8;
    static const ac_fixed<10, 0, false> m_lut[n_segments_lut] = {.0986328125, .099609375, .1044921875, .1103515625, .1201171875, .1328125, .15234375, .1787109375};
    static const ac_fixed<10, 0, false> c_lut[n_segments_lut] = {.0, .0986328125, .1982421875, .3037109375, .4140625, .5341796875, .6669921875, .8193359375};
    // Domain of PWL
    static const ac_fixed<10, 0, false> x_max_lut = .78515625;
    // Scaling constant used later to scale the normalized input from 0 to n_segments_lut
    static const ac_fixed<14, 4, false> sc_constant_lut = ((ac_fixed<14, 4, false>)n_segments_lut)/((ac_fixed<14, 4, false>)(x_max_lut));

    // If the input equals or exceeds pi/4, we halve it and use the formula tan(2*x) = 2*tan(x) / (1 - tan(x)^2) to get the tan value.
    ac_fixed<W + 1, I, false, Q, O> input_int;
    bool input_exceeds_pi_by_4;
    // Keep in mind that the input will only exceed pi/4 if the number of integer bits is greater than or equal to zero.
    if(I >= 0) { input_exceeds_pi_by_4 = (input >= x_max_lut) ? true : false; }

    if((I >= 0) && input_exceeds_pi_by_4) { input_int = (ac_fixed<W + 1, I, false, Q, O>)input >> 1; }
    else { input_int = input; }

    // Compute tan using pwl.

    // Scale the input from 0 to n_segments_lut
    ac_fixed<14, 4, false> x_in_sc = (input_int)*sc_constant_lut;
    // Take out the fractional bits of the scaled input
    ac_fixed<14 - 4, 0, false> x_in_sc_frac;
    x_in_sc_frac.set_slc(0, x_in_sc.template slc<14 - 4>(0));
    ac_int<4, false> index;
    // The integer part of the input is the index of the LUT table
    index = x_in_sc.to_int();
    typedef ac_fixed<20, 0, false> output_pwl_type;
    output_pwl_type output_pwl = m_lut[index]*x_in_sc_frac + c_lut[index];

    // As mentioned earlier, if the input equals or exceeds pi/4, we use the formula tan(2*x) = 2*tan(x) / (1 - tan(x)^2) to get the tan value.
    if((I >= 0) && input_exceeds_pi_by_4) {
      ac_fixed<22, 6, false> recip_value;
      ac_fixed<20, 0, false> one_minus_tan_theta_by_2_sqr = 1 - output_pwl * output_pwl;
      // Use the reciprocal_pwl function to calculate 1 / (1 - tan(x)^2)
      ac_math::ac_reciprocal_pwl<AC_TRN>(one_minus_tan_theta_by_2_sqr, recip_value);
      output_temp = output_pwl * recip_value;
      output_temp <<= 1;
    } else {
      output_temp = output_pwl;
    }

    // If input crosses or equals 89.46 degrees (roughly), set the output to saturate
    if(input >= sat_limit) { output_temp.template set_val<AC_VAL_MAX>(); }

    output = output_temp;

#if !defined(__SYNTHESIS__) && defined(AC_TAN_PWL_H_DEBUG)
    cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl;
    cout << "input           = " << input << endl;
    cout << "sc_constant_lut = " << sc_constant_lut << endl;
    cout << "input_int       = " << input_int << endl;
    cout << "x_in_sc         = " << x_in_sc << endl;
    cout << "output_pwl      = " << output_pwl << endl;
    cout << "output_temp     = " << output_temp << endl;
    cout << "output          = " << output << endl;
#endif
  }

  // The following version enables a return-by-value.
  template<class T_out, class T_in>
  T_out ac_tan_pwl(const T_in &input)
  {
    T_out output;
    ac_tan_pwl(input, output);
    return output;
  }

}

#endif // _INCLUDED_AC_TAN_PWL_H_