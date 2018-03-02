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
// File: ac_pow_pwl.h
//
// Description:
//    Provides piece-wise linear implementations of the
//    base 2 and base e exponential functions for the AC (tm) Datatype:
//    ac_fixed.
//
// Usage:
//    A sample testbench and its implementation look like this:
//
//    #include <ac_fixed.h>
//    #include <ac_math/ac_pow_pwl.h>
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
//      ac_pow2_pwl(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 2.5;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN(0);
//    }
//    #endif
//
// Notes:
//    Attempting to call the function with a type that is not implemented will result
//    in a compile error.
//
// Revision History:
//    Niramay Sanghvi : Aug 15 2017 : Default template parameters for configurability added
//    Niramay Sanghvi : Aug 07 2017 : Used ac_shift_left function for RND and SAT support.
//    Niramay Sanghvi : Jul 31 2017 : Added structs for type-checking.
//    Niramay Sanghvi : Jul 12 2017 : Added header style format.
//    Niramay Sanghvi : Jul 11 2017 : Added support for all possible values of W & I.
//    Niramay Sanghvi : Jul 05 2017 : Passed output by reference.
//    Niramay Sanghvi : Jun 29 2017 : Renamed header files and functions.
//
//*****************************************************************************************

#ifndef _INCLUDED_AC_POW_PWL_H_
#define _INCLUDED_AC_POW_PWL_H_

// The functions use default template parameters, which are only supported by C++11 or later
// compiler standards. Hence, the user should be informed if they are not using those standards.

#if __cplusplus < 201103L
#error Please use C++11 or a later standard for compilation.
#endif

// Include headers for data types supported by these implementations
#include <ac_int.h>
#include <ac_float.h>
#include <ac_fixed.h>
#include <ac_complex.h>

// Include headers for required functions
#include <ac_math/ac_shift.h>

#if !defined(__SYNTHESIS__)
#include <iostream>
using namespace std;
#endif

namespace pow_pwl
{
// This struct computes precision of input_inter variable ("pii") for base e exponent. It also makes sure
// that there are a set minimum no. of fractional bits to represent the multiplication of x with log2(e)
// (this is decided by the n_f_b variable).
  template <int W, int I, bool S, ac_q_mode Q, ac_o_mode O, int n_f_b>
  struct comp_pii_exp {
    enum {
      pit_i       = I + 1,
      pit_w_inter = W + 1,
      pit_w       = (W - I) > n_f_b ? pit_w_inter : pit_i + n_f_b
    };
    typedef ac_fixed<pit_w, pit_i, S, Q, O> pit_t;
  };
};

//=========================================================================
// Function: ac_pow2_pwl (for ac_fixed)
//
// Description:
//    Calculation of base 2 exponential of real inputs, passed as ac_fixed
//    variables.
//
//    Separates input into integer and fractional part, the fractional part
//    is passed to the PWL approximation. The output is then left-shifted
//    by the value of the integer part, in order to de-normalize.
//
// Usage:
//    See above example code for usage.
//
// Notes:
//    The PWL implementation utilizes 3 elements, which has a small impact
//    on accuracy. For better accuracy, please use the LUT_Generator.cpp
//    file to get new values for PWL LUT values which utilize more than 3
//    elements. The relevant arrays and variables with information for the
//    new LUT will be copied to a text file in C++ syntax.
//
//-------------------------------------------------------------------------

namespace ac_math
{
  template<ac_q_mode pwl_Q = AC_RND,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  void ac_pow2_pwl(
    const ac_fixed<W, I, S, Q, O> &input,
    ac_fixed<outW, outI, false, outQ, outO> &output
  )
  {
    // Stores the fractional part of the input. By default it is set to 0
    ac_fixed<AC_MAX(W - I, 1), 0, false> input_frac_part = 0;

    // Take out the fractional part of the input
    // This serves as a sort of normalization, with the fractional part being
    // the normalized data (can only vary from 0 to 0.9999...)

    // Only carry out slicing if the input has a fractional component.
    // If the input doesn't have a fractional part, the default value of input_frac_part, i.e. 0,
    // is suitable to be used in later calculations.
    if (W > I) {input_frac_part.set_slc(0, input.template slc<AC_MAX(W - I, 1)>(0));}

    // Initialization for PWL LUT
    const unsigned n_segments_lut = 3;
    // Initializing the LUT arrays
    static const ac_fixed<9, 0, false> m_lut[n_segments_lut] = {.2587890625, .326171875, .412109375};
    static const ac_fixed<10, 1, false> c_lut[n_segments_lut] = {.998046875, 1.255859375, 1.58203125};
    // Domain of the PWL implementation
    static const ac_fixed<1, 0, false> x_min_lut = 0;
    static const ac_fixed<1, 1, false> x_max_lut = 1;
    // Scaling constant used later to scale the normalized input from 0 to n_segments_lut
    static const ac_fixed<2, 2, false> sc_constant_lut = n_segments_lut/(x_max_lut - x_min_lut);

    // Compute power of two using pwl
    // Scale the normalized input from 0 to n_segments_lut
    ac_fixed<11, 2, false> x_in_sc = (input_frac_part - x_min_lut)*sc_constant_lut;
    ac_fixed<11 - 2, 0, false> x_in_sc_frac;
    // Slice out the fractional part from the scaled input, store it in another variable.
    x_in_sc_frac.set_slc(0, x_in_sc.template slc<11 - 2>(0));
    // The integer part of the scaled input is the index of the LUT table
    ac_int<2, false> index = x_in_sc.to_int();
    typedef ac_fixed<20, 2, false, pwl_Q> output_pwl_type;
    output_pwl_type output_pwl = m_lut[index] * x_in_sc_frac + c_lut[index];

    // Shift left by the integer part of the input to multiply by (2^input_integer_part)
    ac_math::ac_shift_left(output_pwl, input.to_int(), output);

#if !defined(__SYNTHESIS__) && defined(AC_POW_PWL_DEBUG)
    cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl;
    cout << "Actual input              = " << input << endl;
    cout << "normalized input          = " << input_frac_part << endl;
    cout << "output up-scaled by exp   = " << output << endl;
    cout << "index                     = " << index  << endl;
#endif
  }

//=============================================================================
// Version that allows the return of values.
  template<class T_out, ac_q_mode pwl_Q = AC_RND, class T_in>
  T_out ac_pow2_pwl(
    const T_in &input
  )
  {
    // Create a temporary variable for output and use the pass-by-reference version
    // to evaluate it. This temporary variable is returned as the output.
    T_out output;
    ac_pow2_pwl<pwl_Q>(input, output);
    return output;
  }

//=============================================================================
// Function: ac_exp_pwl (for ac_fixed)
//
// Description:
//    Calculation of base e exponential of real inputs, passed as ac_fixed
//    variables.
//
//    Separates input into integer and fractional part, the fractional part
//    is passed to the PWL approximation. The output is then left-shifted
//    by the value of the integer part, in order to de-normalize.
//
// Usage:
//    A sample testbench and its implementation look like this:
//
//    #include <ac_fixed.h>
//    #include <ac_math/ac_pow_pwl.h>
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
//      ac_exp_pwl(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 2.5;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN(0);
//    }
//    #endif
//
// Notes:
//    This function relies on the ac_pow2_pwl function for its computation. It
//    does this by multiplying the input with log2(e), then passing it to
//    the ac_pow2_pwl function. In doing so, we also make sure that the
//    multiplied value has enough precision to store the result of
//    input*log2(e)
//
//-----------------------------------------------------------------------------

//n_f_b = minimum no of fractional bits used in storing the result of multiplication by log2(e)
  template<int n_f_b = 9, ac_q_mode pwl_Q = AC_RND,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  void ac_exp_pwl(
    const ac_fixed<W, I, S, Q, O> &input,
    ac_fixed<outW, outI, false, outQ, outO> &output
  )
  {
    static const ac_fixed<17, 3, true> log2e = 1.44269504089;
    // Find type of intermediate variable used to store output of x*log2(e)
    typedef class pow_pwl::comp_pii_exp<W, I, S, Q, O, n_f_b>::pit_t input_inter_type;
    input_inter_type input_inter;
    // e^x = 2^(x*log2(e))
    input_inter = input*log2e;
    ac_pow2_pwl<pwl_Q>(input_inter, output);

#if !defined(__SYNTHESIS__) && defined(AC_POW_PWL_DEBUG)
    cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl;
    cout << "input_inter.width       = " << input_inter.width << endl;
    cout << "input_inter.i_width     = " << input_inter.i_width << endl;
    cout << "input (power_exp)       = " << input << endl;
    cout << "log2e (power_exp)       = " << log2e << endl;
    cout << "input_inter (power_exp) = " << input_inter << endl;
    cout << "output (power_exp)      = " << output << endl;
#endif
  }

//=============================================================================
// Version that allows the return of values.
  template<class T_out, int n_f_b = 9, ac_q_mode pwl_Q = AC_RND, class T_in>
  T_out ac_exp_pwl(
    const T_in &input
  )
  {
    // Create a temporary variable for output and use the pass-by-reference version
    // to evaluate it. This temporary variable is returned as the output.
    T_out output;
    ac_exp_pwl<n_f_b, pwl_Q>(input, output);
    return output;
  }
}

#endif

