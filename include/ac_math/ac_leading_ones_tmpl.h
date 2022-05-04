/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 3.4                                                 *
 *                                                                        *
 *  Release Date    : Wed May  4 10:47:29 PDT 2022                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.4.3                                               *
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
//*****************************************************************************************
// File: ac_leading_ones_tmpl.h
//
// Description: An example of using C++ template recursion to implement a leading ones
//    detection block. Note that this results in single-bit logic, so as it get bigger
//    timing will get worse. Placing this in a C-CORE will help as long as it can fit
//    in one clock cycle.
//
// Usage:
//
// Notes:
//
// Revision History:
//    3.4.0  - Copied from HLS Bluebook examples
//
//*****************************************************************************************

#ifndef _INCLUDED_AC_LEADING_ONES_TMPL_H_
#define _INCLUDED_AC_LEADING_ONES_TMPL_H_

#include <ac_int.h>

template<int N_BITS>
bool ac_leading_ones(ac_int<N_BITS,false> &din,
                     ac_int<ac::log2ceil<N_BITS>::val,false> &dout)
{
  enum {
    P2 = ac::nbits<(N_BITS+1)/2>::val
  };
  ac_int<N_BITS-P2,false> upper;
  ac_int<P2,false> lower;
  ac_int<ac::log2ceil<N_BITS>::val,0> idx=0;
  ac_int<ac::log2ceil<N_BITS-P2>::val,0> idxu=0;
  ac_int<ac::log2ceil<P2>::val,0> idxl=0;
  static bool flag = false;

  upper.set_slc(0, din.template slc<N_BITS-P2>(P2));
  lower.set_slc(0, din.template slc<P2>(0));
  flag = (din!=0) ? 1 : 0;
  if (upper) {
    ac_leading_ones<N_BITS-P2>(upper,idxu);
    idx = idxu | P2;
  } else {
    ac_leading_ones<P2>(lower,idxl);
    idx = idxl;
  }
  dout = idx;
  return flag;
}

template<>
bool ac_leading_ones<1>(ac_int<1,false> &din,
                        ac_int<1,false> &dout)
{
  dout = 0;
  return din[0];
}

#endif

