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
 *  Copyright 2018 Siemens                                                *
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
//******************************************************************************************
// Function: ac_softplus_pwl (for ac_fixed)
//
// Description:
//    Provides piece-wise linear approximation of the softplus function
//    for the ac_fixed datatype
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_softplus_pwl.h>
//    using namespace ac_math;
//
//    typedef ac_fixed<10, 4, false, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<20, 2, false, AC_RND, AC_SAT> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_softplus_pwl(input,output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 3.5;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    The softplus activation function is equal to log(1+e^x). This function relies on 
//    the pow_pwl and log_pwl functions for its computation. e^x is computed using 
//    e^x = 2^(x*log2(e)). So here the input is multiplied with log2(e) and the passed to 
//    the ac_pow2_pwl function. And then finally 1 + e^x is passed to the ac_log_pwl function. 
//    As the below function relies on the ac_pow_pwl function, make sure to use the same value 
//    of constants n_segments_lut and n_frac_bits used below as defined in the ac_fixed 
//    implementation for ac_pow_pwl.h.    
//*************************************************************************************************

#ifndef _INCLUDED_AC_SOFTPLUS_PWL_H_
#define _INCLUDED_AC_SOFTPLUS_PWL_H_

// The below functions use default template parameters, which are only supported by C++11 or later
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
#include <ac_math/ac_pow_pwl.h>
#include <ac_math/ac_log_pwl.h>

#if !defined(__SYNTHESIS__)
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
#endif

//=========================================================================
// Function: ac_softplus_pwl (for ac_fixed)
//
// Description:
//    Softplus function for real inputs, passed as ac_fixed
//    variables.
//
// Usage:
//    See above example for usage.
//
//-------------------------------------------------------------------------

namespace ac_math
{
  template<ac_q_mode pwl_Q = AC_TRN,
           int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, ac_q_mode outQ, ac_o_mode outO>
  void ac_softplus_pwl(
    const ac_fixed<W, I, S, Q, O> &input,
    ac_fixed<outW, outI, false, outQ, outO> &output
  )
  {
    // n_segments_lut and n_frac_bits should be the same as the corresponding values in the ac_fixed
    // implementation for ac_exp_pwl.
    // The bitwidth calculations which use these constants are designed with a PWL domain of [0, 1) in
    // mind. They will still work for other PWL domains, but at a possibly lower-than-optimal accuracy.
    //
    // NOTE: Change these constants if you change the either of the corresponding values in the
    // ac_fixed implementation for the ac_exp_pwl function.    
    const unsigned n_segments_lut = 4;
    const int n_frac_bits = 10;

    const bool is_n_seg_po2 = !bool(n_segments_lut & (n_segments_lut - 1));
    const int extra_f_bits = is_n_seg_po2 ? ac::nbits<n_segments_lut - 1>::val : 0;

    const int I_width = (I<0) ? -I:I;

    // Find type of intermediate variable used to store output of x*log2(e)
    typedef class comp_pii_exp<W, I_width, S, n_frac_bits + extra_f_bits>::pit_t input_inter_type;

    const int out_pow_width  = input_inter_type::width;
    const int out_pow_iwidth = input_inter_type::i_width;  

    // Since scaling constant is a positive power-of-two, multiplication with it is the same as left-shifting by 2.
    // Accordingly, the scaled normalized input will have 2 less fractional bits than the normalized input, provided that this
    // number of fractional bits is lesser than n_frac_bits. If not, the number of fractional bits in the scaled input is set to n_frac_bits.
    const int sc_input_frac_bits = AC_MAX(1, AC_MIN(n_frac_bits, out_pow_width - out_pow_iwidth - 2));

    typedef ac_fixed<((out_pow_iwidth + 1)*out_pow_iwidth) + sc_input_frac_bits + n_frac_bits, ((out_pow_iwidth + 1)*out_pow_iwidth), false, AC_TRN, AC_WRAP> out_pow_type;

    out_pow_type out_pow;

    ac_exp_pwl(input, out_pow);

    typedef ac_fixed<out_pow_type::width + 1, out_pow_type::i_width + 1, false, AC_TRN, AC_WRAP> inp_log_type;

    inp_log_type inp_log = out_pow + 1;

    ac_log_pwl<pwl_Q>(inp_log, output);
  }

  // The following version enables a return-by-value.
  template<class T_out,
           ac_q_mode pwl_Q = AC_TRN,
           class T_in>
  T_out ac_softplus_pwl(const T_in &input)
  {
    T_out output;
    ac_softplus_pwl<pwl_Q>(input, output);
    return output;
  }  
}
#endif    
