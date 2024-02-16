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
// Function: ac_leakyrelu (for ac_fixed)
//
// Description:
//    Provides an implementation for the leakyrelu function
//    for the ac_fixed datatype
//
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    #include <ac_math/ac_leakyrelu.h>
//    using namespace ac_math;
//
//    typedef ac_fixed<10, 4, false, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<11, 4, false, AC_RND, AC_SAT> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      alpha_type &alpha
//    )
//    {
//      ac_leakyrelu(input,output,alpha);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 3.5;
//      output_type output;
//      CCS_DESIGN(project)(input, output, alpha);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:    
//*************************************************************************************************

#ifndef _INCLUDED_AC_LEAKYRELU_H_
#define _INCLUDED_AC_LEAKYRELU_H_

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

#if !defined(__SYNTHESIS__)
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
#endif

//=========================================================================
// Function: ac_leakyrelu (for ac_fixed)
//
// Description:
//    Selu function for real inputs, passed as ac_fixed
//    variables.
//
// Usage:
//    See above example for usage.
//
//-------------------------------------------------------------------------

namespace ac_math
{
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O,
           int outW, int outI, bool outS, ac_q_mode outQ, ac_o_mode outO,
           int alphaW, int alphaI, bool alphaS, ac_q_mode alphaQ, ac_o_mode alphaO>
  void ac_leakyrelu(const ac_fixed<W, I, S, Q, O> &input,
    ac_fixed<outW, outI, outS, outQ, outO> &output,
    const ac_fixed<alphaW, alphaI, alphaS, alphaQ, alphaO> &alpha
  )
  { 
    if (input>0)
      output = input;
    else 
      output = alpha*input;
  }

  // The following version enables a return-by-value.
  template<class T_out,
           class T_in,
           class T_alpha>
  T_out ac_leakyrelu(const T_in &input, const T_alpha &alpha)
  {
    T_out output;
    ac_leakyrelu(input, output, alpha);
    return output;
  }  
}
#endif    
