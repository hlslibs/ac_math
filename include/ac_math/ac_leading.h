/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 3.6                                                 *
 *                                                                        *
 *  Release Date    : Tue Nov 12 23:14:00 PST 2024                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.6.0                                               *
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

#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>

#include <ac_int.h>
#include <ac_fixed.h>
#include <mc_scverify.h>
using namespace std;

namespace ac_math
{
  /*
  #################################################################
  For Unisgned ac_int implementation
  #################################################################
  */

//ac_int Unsigned number implementation Function to find the position of the leading 1's bit position from the LSB side.
  template<int W>
  ac_int<ac::nbits<W>::val,false> leading1(ac_int<W,false> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<ac::nbits<W>::val,false> lead_sign=din.leading_sign(lead_sign_flag);
    // std::cout << "lead_sign: "<< lead_sign << std::endl;
    flag=1;
    if (W==lead_sign) {flag=0;}
    if (lead_sign_flag) { return 0; }
    else { return (W-1-lead_sign); }
  }

//ac_int  Unsigned number implementation Function to find the position of the trailing 1's bit position from the LSB side.
  template<int W>
  ac_int<ac::nbits<W-1>::val,false> trailing1(ac_int<W,false> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<W,false> rev;
    #pragma hls_unroll yes
    for (int i=0; i<=W-1; i++) {
      rev[i] = din[W-1-i];
    }
    ac_int<ac::nbits<W>::val,false> lead_sign=rev.leading_sign(lead_sign_flag);
    // std::cout << "lead_sign: "<< lead_sign << std::endl;
    flag=1;
    if (W==lead_sign) {flag=0;}
    if (lead_sign_flag) { return 0; }
    else { return (lead_sign); }
  }

//ac_int Unsigned number implementation Function to find the bit position of the leading 0 from the LSB side.
  template<int W>
  ac_int<ac::nbits<W-1>::val,false> leading0(ac_int<W,false> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<W,false> cmp = ~din;
    ac_int<ac::nbits<W>::val,false> lead_sign=cmp.leading_sign(lead_sign_flag);
    // std::cout << "lead_sign: "<< lead_sign << std::endl;
    flag=1;
    if (W==lead_sign) {flag=0;}
    if (lead_sign_flag) { return 0; }
    else { return (W-1-lead_sign); }
  }

//ac_int Unsigned number implementation Function to find the bit position of the trailing 0 from the LSB side.
  template<int W>
  ac_int<ac::nbits<W-1>::val,false> trailing0(ac_int<W,false> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<W,false> rev;
    #pragma hls_unroll yes
    for (int i=0; i<=W-1; i++) {
      rev[i] = !din[W-1-i];
    }
    ac_int<ac::nbits<W>::val,false> lead_sign=rev.leading_sign(lead_sign_flag);
    // std::cout << "lead_sign: "<< lead_sign << std::endl;
    flag=1;
    if (W==lead_sign) {flag=0;}
    if (lead_sign_flag) { return 0; }
    else { return (lead_sign); }
  }

  /*
  #################################################################
  For Signed ac_int implementation
  #################################################################
  */

//ac_int Signed number implementation Function to find the position of the leading 1's bit position from the LSB side.
  template<int W>
  ac_int<ac::nbits<W>::val,false> leading1(ac_int<W,true> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<W-1,false> truncated;
    truncated.set_slc(0,(din.template slc<W-1>(0)));
    ac_int<ac::nbits<W>::val,false> lead_sign=truncated.leading_sign(lead_sign_flag);
    // std::cout << "lead_sign: "<< lead_sign << std::endl;
    flag=1;
    if (W-1==lead_sign) {flag=0;}
    if (lead_sign_flag) { return 0; }
    else { return (W-2-lead_sign); }
  }

//ac_int Signed number implementation Function to find the position of the trailing 1's bit position from the LSB side.
  template<int W>
  ac_int<ac::nbits<W-1>::val,false> trailing1(ac_int<W,true> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<W-1,false> rev,truncated;
    truncated.set_slc(0,(din.template slc<W-1>(0)));

    #pragma hls_unroll yes
    for (int i=0; i<=W-2; i++) {
      rev[i] = truncated[W-2-i];
    }
    ac_int<ac::nbits<W>::val,false> lead_sign=rev.leading_sign(lead_sign_flag);
    // std::cout << "lead_sign: "<< lead_sign << std::endl;
    flag=1;
    if (W-1==lead_sign) {flag=0;}
    if (lead_sign_flag) { return 0; }
    else { return (lead_sign); }
  }

// //ac_int Signed number implementation Function to find the bit position of the leading 0 from the LSB side.
  template<int W>
  ac_int<ac::nbits<W-1>::val,false> leading0(ac_int<W,true> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<W-1,false> cmp;
    cmp.set_slc(0,~(din.template slc<W-1>(0)));
    // Find the leading one in the one's complement, which is equivalent to the leading zero in truncated_din
    ac_int<ac::nbits<W-1>::val,false> lead_sign = cmp.leading_sign(lead_sign_flag);
    // Output adjustment based on leading_sign output
    flag=1;
    if (W-1==lead_sign) {flag=0;}
    if (lead_sign_flag) {
      return 0; // If flag is true, all bits in cmp are 1, implying all were 0 in truncated_din
    } else {
      return W-2 - lead_sign; // Adjusted for the fact that we're working with W-1 bits now
    }
  }

//ac_int Signed number implementation Function to find the bit position of the trailing 0 from the LSB side.
  template<int W>
  ac_int<ac::nbits<W-1>::val,false> trailing0(ac_int<W,true> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<W-1,false> cmp,rev;
    cmp.set_slc(0,~(din.template slc<W-1>(0)));
    #pragma hls_unroll yes
    for (int i=0; i<=W-2; i++) {
      rev[i] = cmp[W-2-i];
    }
    ac_int<ac::nbits<W>::val,false> lead_sign=rev.leading_sign(lead_sign_flag);
    // std::cout << "lead_sign: "<< lead_sign << std::endl;
    flag=1;
    if (W-1==lead_sign) {flag=0;}
    if (lead_sign_flag) { return 0; }
    else { return (lead_sign); }
  }

  /*
  #################################################################
  For Unisgned ac_fixed implementation
  #################################################################
  */

//ac_fixed Unsigned number implementation Function to find the position of the leading 1's bit position from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W>::val,false> leading1(ac_fixed<W,I,false> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<ac::nbits<W>::val,false> lead_sign=din.leading_sign(lead_sign_flag);
    // std::cout << "lead_sign: "<< lead_sign << std::endl;
    flag=1;
    if (W==lead_sign) {flag=0;}
    if (lead_sign_flag) { return 0; }
    else { return (W-1-lead_sign); }
  }

//ac_fixed  Unsigned number implementation Function to find the position of the trailing 1's bit position from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W-1>::val,false> trailing1(ac_fixed<W,I,false> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<W,false> rev;
    #pragma hls_unroll yes
    for (int i=0; i<=W-1; i++) {
      rev[i] = din[W-1-i];
    }
    ac_int<ac::nbits<W>::val,false> lead_sign=rev.leading_sign(lead_sign_flag);
    // std::cout << "lead_sign: "<< lead_sign << std::endl;
    flag=1;
    if (W==lead_sign) {flag=0;}
    if (lead_sign_flag) { return 0; }
    else { return (lead_sign); }
  }

//ac_fixed Unsigned number implementation Function to find the bit position of the leading 0 from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W-1>::val,false> leading0(ac_fixed<W,I,false> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_fixed<W,I,false> cmp=~din;
    ac_int<ac::nbits<W>::val,false> lead_sign=cmp.leading_sign(lead_sign_flag);
    // std::cout << "lead_sign: "<< lead_sign << std::endl;
    flag=1;
    if (W==lead_sign) {flag=0;}
    if (lead_sign_flag) { return 0; }
    else { return (W-1-lead_sign); }
  }

//ac_fixed Unsigned number implementation Function to find the bit position of the trailing 0 from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W-1>::val,false> trailing0(ac_fixed<W,I,false> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<W,false> rev;
    #pragma hls_unroll yes
    for (int i=0; i<=W-1; i++) {
      rev[i] = !din[W-1-i];
    }
    ac_int<ac::nbits<W>::val,false> lead_sign=rev.leading_sign(lead_sign_flag);
    // std::cout << "lead_sign: "<< lead_sign << std::endl;
    flag=1;
    if (W==lead_sign) {flag=0;}
    if (lead_sign_flag) { return 0; }
    else { return (lead_sign); }
  }

  /*
  #################################################################
  For Signed ac_fixed implementation
  #################################################################
  */

//ac_fixed Signed number implementation Function to find the position of the leading 1's bit position from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W>::val,false> leading1(ac_fixed<W,I,true> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<W-1,false> truncated;
    truncated.set_slc(0,(din.template slc<W-1>(0)));
    ac_int<ac::nbits<W>::val,false> lead_sign=truncated.leading_sign(lead_sign_flag);
    // std::cout << "lead_sign: "<< lead_sign << std::endl;
    flag=1;
    if (W-1==lead_sign) {flag=0;}
    if (lead_sign_flag) { return 0; }
    else { return (W-2-lead_sign); }
  }

//ac_fixed Signed number implementation Function to find the position of the trailing 1's bit position from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W-1>::val,false> trailing1(ac_fixed<W,I,true> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<W-1,false> rev,truncated;
    truncated.set_slc(0,(din.template slc<W-1>(0)));

    #pragma hls_unroll yes
    for (int i=0; i<=W-2; i++) {
      rev[i] = truncated[W-2-i];
    }
    ac_int<ac::nbits<W>::val,false> lead_sign=rev.leading_sign(lead_sign_flag);
    // std::cout << "lead_sign: "<< lead_sign << std::endl;
    flag=1;
    if (W-1==lead_sign) {flag=0;}
    if (lead_sign_flag) { return 0; }
    else { return (lead_sign); }
  }

//ac_fixed Signed number implementation Function to find the bit position of the leading 0 from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W-1>::val,false> leading0(ac_fixed<W,I,true> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<W-1,false> cmp;
    cmp.set_slc(0,~(din.template slc<W-1>(0)));
    // Find the leading one in the one's complement, which is equivalent to the leading zero in truncated_din
    ac_int<ac::nbits<W-1>::val,false> lead_sign = cmp.leading_sign(lead_sign_flag);
    // Output adjustment based on leading_sign output
    flag=1;
    if (W-1==lead_sign) {flag=0;}
    if (lead_sign_flag) {
      return 0; // If flag is true, all bits in cmp are 1, implying all were 0 in truncated_din
    } else {
      return W-2 - lead_sign; // Adjusted for the fact that we're working with W-1 bits now
    }
  }

//ac_fixed Signed number implementation Function to find the bit position of the trailing 0 from the LSB side.
  template<int W, int I>
  ac_int<ac::nbits<W-1>::val,false> trailing0(ac_fixed<W,I,true> din, bool &flag)
  {
    bool lead_sign_flag;
    ac_int<W-1,false> cmp,rev;
    cmp.set_slc(0,~(din.template slc<W-1>(0)));
    #pragma hls_unroll yes
    for (int i=0; i<=W-2; i++) {
      rev[i] = cmp[W-2-i];
    }
    ac_int<ac::nbits<W>::val,false> lead_sign=rev.leading_sign(lead_sign_flag);
    // std::cout << "lead_sign: "<< lead_sign << std::endl;
    flag=1;
    if (W-1==lead_sign) {flag=0;}
    if (lead_sign_flag) { return 0; }
    else { return (lead_sign); }
  }



};

#endif
