/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 3.8                                                 *
 *                                                                        *
 *  Release Date    : Tue May 13 15:34:32 PDT 2025                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.8.1                                               *
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
//******************************************************************************************
// Function: ac_gelu_pwl (for ac_fixed)
//
// Description:
//    Provides piece-wise linear approximation of the gelu function
//    for the ac_fixed datatype
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_gelu_pwl.h>
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
//      ac_gelu_pwl(input,output);
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

#ifndef _INCLUDED_AC_GELU_PWL_H_
#define _INCLUDED_AC_GELU_PWL_H_

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
#include <ac_math/ac_tanh_pwl.h>

#if !defined(__SYNTHESIS__)
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
#endif

//=========================================================================
// Function: ac_gelu_pwl (for ac_fixed)
//
// Description:
//    Gelu function for real inputs, passed as ac_fixed
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
  void ac_gelu_pwl(
    const ac_fixed<W, I, S, Q, O> &input,
    ac_fixed<outW, outI, outS, outQ, outO> &output
  )
  {
    const int n_frac_bits = 10;

    typedef ac_fixed<outW, outI, outS, outQ, outO> T_out;

    const int tanh_intwidth = ((I<W) && (I==1))?1:((I<W) && (I==0))?0:(I<0)?(I):(3*I - 3);

    typedef ac_fixed<W + 2*I + 10, tanh_intwidth, outS, outQ, outO> input_tanh_type;

    typedef ac_fixed<n_frac_bits*2 + 1, 1, true, pwl_Q> output_tanh_type;

    typedef ac_fixed<n_frac_bits*2 + 1 + I, I + 1, true, pwl_Q> output_gelu_type;

    ac_fixed<20, 1, false> aconst = 0.7978845608;

    ac_fixed<20, 1, false> bconst = 0.044715;

    ac_fixed<3, 1, false> cconst = 0.5;

    input_tanh_type tanh_input; 

    output_tanh_type tanh_output; 

    output_gelu_type gelu_output;

    tanh_input = aconst * (input + (bconst * input * input * input));

    tanh_output = ac_tanh_pwl<output_tanh_type, pwl_Q>(tanh_input);   

    output = cconst * input * (1 + tanh_output);

  }

  // The following version enables a return-by-value.
  template<class T_out,
           ac_q_mode pwl_Q = AC_TRN,
           class T_in>
  T_out ac_gelu_pwl(const T_in &input)
  {
    T_out output;
    ac_gelu_pwl<pwl_Q>(input, output);
    return output;
  }  
}
#endif    
