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
// Description:
// Usage:
// Notes:
// Revision History:
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//******************************************************************************************

#ifndef _INCLUDED_AC_LEADING_H_
#define _INCLUDED_AC_LEADING_H_

#if (defined(__GNUC__) && (__cplusplus < 201103L))
#error Please use C++11 or a later standard for compilation.
#endif
#if (defined(_MSC_VER) && (_MSC_VER < 1920) && !defined(__EDG__))
#error Please use Microsoft VS 2019 or a later standard for compilation.
#endif

#include <ac_int.h>
#include <ac_fixed.h>

namespace ac_math
{


  /*
  #################################################################
  For Unisgned ac_int implementation
  #################################################################
  */

//ac_int Unsigned number implementation Function to find the position of the leading 1's bit position from the LSB side.
  template<int W>
  ac_int<ac::nbits<W-1>::val,false> leading1(ac_int<W,false> din, bool &flag)
  {
    return ~(din.leading_sign(flag));
  }

//ac_int  Unsigned number implementation Function to find the position of the trailing 1's bit position from the LSB side.
  template<int W>
  ac_int<ac::nbits<W-1>::val,false> trailing1(ac_int<W,false> din, bool &flag)
  {
    ac_int<W,false> rev;
#pragma hls_unroll yes
    for (int i=0; i<W; i++)
    { rev[i] = din[W-1-i]; }
    return rev.leading_sign(flag);
  }

//ac_int Unsigned number implementation Function to find the bit position of the leading 0 from the LSB side.
  template<int W>
  ac_int<ac::nbits<W-1>::val,false> leading0(ac_int<W,false> din, bool &flag)
  {
    ac_int<W,false> cmp = ~din;
    return ~(cmp.leading_sign(flag));
  }

//ac_int Unsigned number implementation Function to find the bit position of the trailing 0 from the LSB side.
  template<int W>
  ac_int<ac::nbits<W-1>::val,false> trailing0(ac_int<W,false> din, bool &flag)
  {
    ac_int<W,false> rev;
#pragma hls_unroll yes
    for (int i=0; i<W; i++)
    { rev[i] = !din[W-1-i]; }
    return rev.leading_sign(flag);
  }

  /*
  #################################################################
  For Signed ac_int implementation
  #################################################################
  */

//ac_int  Signed number implementation Function to find the position of the leading 1's bit position from the LSB side.
  template<int W>
  ac_int<ac::nbits<W-1>::val,false> leading1(ac_int<W,true> din, bool &flag)
  {
    din[W-1]=0;
    ac_int<W,false> reg = din;
    return ~(reg.leading_sign(flag));
  }

//ac_int  Signed number implementation Function to find the position of the trailing 1's bit position from the LSB side.
  template<int W>
  ac_int<ac::nbits<W-1>::val,false> trailing1(ac_int<W,true> din, bool &flag)
  {
    ac_int<W-1,false> trail_reg;
#pragma hls_unroll yes
    for (int i=0; i<W-1; i++)
    { trail_reg[i] = din[W-1-1-i]; }
    return trail_reg.leading_sign(flag);
  }

//ac_int Signed number implementation Function to find the bit position of the leading 0 from the LSB side.
  template<int W>
  ac_int<ac::nbits<W-1>::val,false> leading0(ac_int<W,true> din, bool &flag)
  {
    ac_int<W,false> reg = reg.set_slc(0,~(din.template slc<31>(0)));
    reg[W-1]=0;
    return ~(reg.leading_sign(flag));
  }

//ac_int Signed number implementation Function to find the bit position of the trailing 0 from the LSB side.
  template<int W>
  ac_int<ac::nbits<W-1>::val,false> trailing0(ac_int<W,true> din, bool &flag)
  {
    ac_int<W-1,false> trail_reg;
#pragma hls_unroll yes
    for (int i=0; i<W-1; i++)
    { trail_reg[i] = !(din[W-2-i]); }
    return trail_reg.leading_sign(flag);
  }

  /*
  #################################################################
  For Unisgned ac_fixed implementation
  #################################################################
  */

// Unsigned Fixed number implementation Function to find the position of the leading 1's bit position from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W-1>::val,false> leading1(ac_fixed<W, I, false> din, bool &flag)
  {
    return ~(din.leading_sign(flag));
  }

  // Unsigned Fixed number implementation Function to find the position of the trailing 1's bit position from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W-1>::val,false> trailing1(ac_fixed<W, I, false> din, bool &flag)
  {
    ac_fixed<W, I, false> rev;
#pragma hls_unroll yes
    for (int i=0; i<W; i++)
    { rev[i] = din[W-1-i]; }
    return rev.leading_sign(flag);
  }

//Fixed Unsigned number implementation Function to find the bit position of the leading 0 from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W-1>::val,false> leading0(ac_fixed<W, I, false> din, bool &flag)
  {
    ac_fixed<W, I, false> cmp = ~din;
    return ~(cmp.leading_sign(flag));
  }

//Fixed Unsigned number implementation Function to find the bit position of the trailing 0 from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W-1>::val,false> trailing0(ac_fixed<W, I, false> din, bool &flag)
  {
    ac_fixed<W, I, false> rev;
#pragma hls_unroll yes
    for (int i=0; i<W; i++)
    { rev[i] = !din[W-1-i]; }
    return rev.leading_sign(flag);
  }

  /*
  #################################################################
  For Signed ac_fixed implementation
  #################################################################
  */
// Fixed Signed number implementation Function to find the position of the leading 1's bit position from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W-1>::val,false> leading1(ac_fixed<W, I, true> din, bool &flag)
  {
    din[W-1]=0;
    ac_fixed<W, I, false> reg = din;
    return ~(reg.leading_sign(flag));
  }

//Fixed Signed number implementation Function to find the position of the trailing 1's bit position from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W-1>::val,false> trailing1(ac_fixed<W, I, true> din, bool &flag)
  {
    ac_fixed<W-1, I, false> trail_reg;
#pragma hls_unroll yes
    for (int i=0; i<=W-2; i++)
    { trail_reg[i] = din[W-2-i]; }
    return trail_reg.leading_sign(flag);
  }

//Fixed Signed number implementation Function to find the bit position of the leading 0 from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W-1>::val,false> leading0(ac_fixed<W, I, true> din, bool &flag)
  {
    ac_fixed<W, I, false> reg = reg.set_slc(0,~(din.template slc<31>(0)));
    reg[W-1]=0;
    return ~(reg.leading_sign(flag));
  }

// Fixed Signed number implementation Function to find the bit position of the trailing 0 from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W-1>::val,false> trailing0(ac_fixed<W, I, true> din, bool &flag)
  {
    ac_fixed<W-1, I, false> trail_reg;
#pragma hls_unroll yes
    for (int i=0; i<=W-2; i++)
    { trail_reg[i] = !(din[W-2-i]); }
    return trail_reg.leading_sign(flag);
  }

};

#endif
