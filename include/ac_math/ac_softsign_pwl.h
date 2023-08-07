/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 3.5                                                 *
 *                                                                        *
 *  Release Date    : Sun Jul 23 16:34:46 PDT 2023                        *
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
// Function: ac_softsign_pwl (for ac_fixed)
//
// Description:
//    Provides piece-wise linear approximation of the softsign function
//    for the ac_fixed datatype
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_softsign_pwl.h>
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
//      ac_softsign_pwl(input,output);
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
//*************************************************************************************************

#ifndef _INCLUDED_AC_SOFTSIGN_PWL_H_
#define _INCLUDED_AC_SOFTSIGN_PWL_H_

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
#include <ac_math/ac_reciprocal_pwl.h>

#if !defined(__SYNTHESIS__)
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
#endif

//=========================================================================
// Function: ac_softsign_pwl (for ac_fixed)
//
// Description:
//    Softsign function for real inputs, passed as ac_fixed
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
           int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO>
  void ac_softsign_pwl(
    const ac_fixed<W, I, S, Q, O> &input,
    ac_fixed<outW, outI, outS, outQ, outO> &output
  )
  {
    // n_frac_bits should be the same as the corresponding values in the ac_fixed implementation for ac_reciprocal_pwl.
    // The bitwidth calculations which use these constants are designed with a PWL domain of [0.5, 1) in
    // mind. They will still work for other PWL domains, but at a possibly lower-than-optimal accuracy.
    //
    // NOTE: Change these constants if you change the either of the corresponding values in the
    // ac_fixed implementation for the ac_reciprocal_pwl function.    
    const int n_frac_bits = 10;

    // Since scaling constant is a positive power-of-two, multiplication with it is the same as left-shifting by 4.
    // Accordingly, the scaled normalized input will have 4 less fractional bits than the normalized input, provided that this
    // number of fractional bits is lesser than n_frac_bits. If not, the number of fractional bits in the scaled input is set to n_frac_bits.
    const int sc_input_frac_bits = AC_MAX(1, AC_MIN(n_frac_bits, W - 4));

    ac_fixed<W, I, false> input_abs_value;

    ac_fixed<W, 0, false> normalized_fixed;
    
    #pragma hls_waive CNS
    if (S) {
      input_abs_value = ((input >= 0) ? (ac_fixed <W, I, false>)input : (ac_fixed <W, I, false>)(-input));
    }
    // If input is unsigned, assign value of input to intermediate variable.
    else {
      input_abs_value = input;
    }

    const int I_width = (I<0) ? -I:I;

    const int W_width = (W<I_width) ? W+I_width:W;

    typedef ac_fixed<W_width + 1, I_width + 1, false, AC_TRN, AC_WRAP> rec_inp_type;

    rec_inp_type rec_inp = input_abs_value + 1;

    typedef ac_fixed<sc_input_frac_bits + n_frac_bits + 1, 1, outS, pwl_Q> rec_out_type;

    rec_out_type out_rec;

    ac_reciprocal_pwl(rec_inp, out_rec);

    output = input * out_rec;
  }

  // The following version enables a return-by-value.
  template<class T_out,
           ac_q_mode pwl_Q = AC_TRN,
           class T_in>
  T_out ac_softsign_pwl(const T_in &input)
  {
    T_out output;
    ac_softsign_pwl<pwl_Q>(input, output);
    return output;
  }  
}
#endif    
