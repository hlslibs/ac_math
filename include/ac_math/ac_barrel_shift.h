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
// File: ac_tanh_pwl.h
//
// Description:
//     Provides Barrel Shift Right implimentation for ac_int datatype, Barrel Left Shift can
//     be conversly implemented with same hardware.
// Usage:
//    A sample testbench and its implementation looks like this:
//
//    void design( ac_int<NB,false> &data_out,
//             ac_int<NB,false> data_inp,
//             ac_int< ac::nbits< NB - 1 >::val,false> s_in  ){
//    Calling of ac_barrel_shift templetized funtion
//    data_out = ac_barrel_shift(data_inp, s_in);
//    }
//
// Notes:
//    Attempting to call this file with a datatype that is not implemented will
//    result in a compile-time error.
//
// Revision History:
// 3.2.4  -  Barrel Shift Added
//*************************************************************************************************

#ifndef __INCLUDED_AC_BARREL_SHIFT__
#define __INCLUDED_AC_BARREL_SHIFT__

#include<ac_int.h>
//template< int NUM_BITS > ac_int<NUM_BITS, false> barrel_shift( ac_int<NUM_BITS,false> , ac_int< ac::nbits< NUM_BITS >::val,false> );

// Barrel Shift block is high throughput block to provide optimum latency for Dynamic rotate right shift
// ac_barrel_shift class can be instantiated as a block or Combo-CCore  based on requirements
// snippet of funtion instantiation is below
// Shifted_Output = ac_barrel_shift< Number_of_Bits >(Input, Control_Bits);
// IO-port of barrel shifter will be as follows
//  Port:
//    Data_Input : Number_of_Bits" wide port for input data
//    Shift_IN   : log2("Number_of_Bits") Bits wide port for right shift number
//    Data_Output: Number_of_Bits" wide port for output data
namespace ac_math
{
//CCORE is uniquified and constant propagates the index "i"
//RTL characterization is enabled along with delay contraint to crunch the CCORE
  #pragma hls_design ccore
  #pragma hls_ccore_type combinational
  template< int N_BITS >
  void shift_layer(ac_int<N_BITS,false> din_in,  ac_int<N_BITS,false> &din_tmp_out,  ac_int< ac::nbits< N_BITS >::val,false> s, ac_int< ac::nbits< N_BITS >::val,false> i)
  {
    ac_int<N_BITS*2,false> din_tmp;
    din_tmp.set_slc(0,din_in);
    din_tmp.set_slc(N_BITS,din_in);
    if (s[i] == 1)
    { din_tmp = (din_tmp >> (1<<i)); }
    din_tmp_out = din_tmp. template slc <N_BITS>(0);
  }

  //All layers are CCOREd and clocke overhead is set to 0
  #pragma hls_design ccore
  template< int NUM_BITS >
  ac_int<NUM_BITS,false> ac_barrel_shift(ac_int<NUM_BITS,false> din,
                                         ac_int< ac::nbits< NUM_BITS - 1 >::val,false> s_in)
  {
    ac_int< ac::nbits< NUM_BITS -1 >::val,false> s = s_in;
    ac_int<NUM_BITS,false> din_tmp = 0;
    din_tmp = din;
#pragma hls_unroll yes
    for (int i=0; i< ac::nbits< NUM_BITS -1 >::val; i++) {
      shift_layer(din_tmp,din_tmp,s, i);
    }
    return din_tmp;
  }
} //namespace ac_math
#endif


