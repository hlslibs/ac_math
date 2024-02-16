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
//*****************************************************************************************
// File: ac_sqrt.h
//
// Description: Provides square root functions for AC datatypes
//    + input: argument
//    + output: square root
//
//    * Integer Unsigned
//      void ac_sqrt(ac_int<XW,false> x, ac_int<OW,OS> &sqrt)
//
//    * Fixed Point Unsigned
//      void ac_sqrt(ac_fixed<XW,XI,false,XQ,XO> x, ac_fixed<OW,OI,false,OQ,OO> &sqrt)
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_sqrt.h>
//    using namespace ac_math;
//
//    typedef ac_int<20, false> input_type;
//    typedef ac_int<24, false> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_sqrt(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 4;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    This file uses C++ function overloading to target implementations
//    specific to each type of data. Attempting to call the function
//    with a type that is not implemented will result in a compile error.
//
//    This header file also uses the ac_shift_right() and ac_shift_left() functions from
//    the ac_shift header file.
//
// Revision History:
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//    3.4.1  - Added special input handling for ac_std_float and ac_ieee_float implementations.
//    25.1.0 - [CAT-28829] Added ac_float, ac_std_float and ac_ieee_float support.
//    3.3.0  - [CAT-25798] Added CDesignChecker fixes/waivers for code check and Synthesis-simulation mismatch/violations in ac_math PWL and Linear Algebra IPs.
//    3.1.0  - Fixed lines that could cause OVL violations.
//    2.0.10 - Official open-source release as part of the ac_math library.
//
//*****************************************************************************************

#ifndef _INCLUDED_AC_SQRT_H_
#define _INCLUDED_AC_SQRT_H_

// Include headers for data types supported by these implementations
#include <ac_int.h>
#include <ac_fixed.h>
#include <ac_float.h>
#include <ac_std_float.h>

// Include header for required functions
#include <ac_math/ac_shift.h>

// The floating point functions use default template parameters, which are only supported by C++11 or
// later compiler standards. Hence, the user should be informed if they are not using those standards.
#if (defined(__GNUC__) && (__cplusplus < 201103L))
#error Please use C++11 or a later standard for compilation.
#endif
#if (defined(_MSC_VER) && (_MSC_VER < 1920) && !defined(__EDG__))
#error Please use Microsoft VS 2019 or a later standard for compilation.
#endif

namespace ac_math
{

//=========================================================================
// Function: ac_sqrt (for ac_int)
//
// Description:
//    Calculation of square root of real inputs, passed as ac_int
//    variables.
//    + input: unsigned ac_int argument
//    + output: ac_int square root
//
//-------------------------------------------------------------------------

  template< int XW, int OW, bool OS >
  void ac_sqrt(
    ac_int<XW,false> x,
    ac_int<OW,OS> &sqrt
  )
  {
    const int RW = (XW+1)/2;
    // masks used only to hint synthesis on precision
    ac_int<RW+2, false> mask_d = 0;

    ac_int<RW+2, false> d = 0;
    ac_int<RW, false> r = 0;
    ac_int<2*RW, false> z = x;

    // needs to pick 2 bits of z for each iteration starting from
    // the 2 MSB bits. Inside loop, z will be shifted left by 2 each
    // iteration. The bits of interest are always on the same
    // position (z_shift+1 downto z_shift)
    unsigned int z_shift = (RW-1)*2;

    for (int i = RW-1; i >= 0; i--) {
      r <<= 1;

      mask_d = (mask_d << 2) | 0x3;
      d = (mask_d & (d << 2)) | ((z >> z_shift) & 0x3 );

      ac_int<RW+2, false> t = (ac_int<RW+2, false>)(d - (( ((ac_int<RW+1, false>)r) << 1) | 0x1));
      if ( !t[RW+1] ) { // since t is unsigned, look at MSB
        r |= 0x1;
        d = mask_d & t;
      }
      z <<= 2;
    }
    sqrt = r;
  }

//=========================================================================
// Function: ac_sqrt (for ac_fixed)
//
// Description:
//    Calculation of square root of real inputs, passed as ac_fixed
//    variables.
//    + input: unsigned ac_fixed argument
//    + output: ac_fixed square root
//
//-------------------------------------------------------------------------

  template< int XW, int XI, ac_q_mode XQ, ac_o_mode XO,
            int OW, int OI, ac_q_mode OQ, ac_o_mode OO >
  void ac_sqrt(
    ac_fixed<XW,XI,false,XQ,XO> x,
    ac_fixed<OW,OI,false,OQ,OO> &sqrt
  )
  {
    const int RBIT = (OQ == AC_TRN || OQ == AC_TRN_ZERO) ? 0 : 1;
    const int OF = (OW-OI) + RBIT;

    const int RI = (XI+1)/2;

    #pragma hls_waive CNS
    if (RI-1 < -OF) {
      // MSB of result is smaller than LSB of requested output
      sqrt = 0.0;
      return;
    }

    // max is used to avoid compilation problems with non pos bitwidth
    const int RF = AC_MAX(OF, -RI+1);  // OF may be negative
    const int RW = RI + RF;

    // store relevant bits of x in z
    const int ZF = 2*AC_MIN((XW-XI+1)/2, RF);
    const int ZW = 2*RI+ZF;
    ac_fixed<ZW,ZW,false> z_fx;
    ac_math::ac_shift_left(x, ZF, z_fx);
    ac_int<ZW,false> z = z_fx.template slc<ZW>(0);

    // masks used only to hint synthesis on precision
    ac_int<RW+2,false> mask_d = 0;

    ac_int<RW+2,false> d = 0;
    ac_int<RW,false> r = 0;

    // needs to pick 2 bits of z for each iteration starting from
    // the 2 MSB bits. Inside loop, z will be shifted left by 2 each
    // iteration. The bits of interest are always on the same
    // position (z_shift+1 downto z_shift)
    unsigned int z_shift = ZW-2;

    for (int i = RW-1; i >= 0; i--) {
      r <<= 1;

      mask_d = (mask_d << 2) | 0x3;
      d = (mask_d & (d << 2)) | ((z >> z_shift) & 0x3 );

      ac_int<RW+2,false> t = (ac_int<RW+2, false>)(d - (( ((ac_int<RW+1,false>)r) << 1) | 0x1));
      if ( !t[RW+1] ) { // since t is unsigned, look at MSB
        r |= 0x1;
        d = mask_d & t;
      }
      z <<= 2;
    }

    ac_fixed<RW+1,RW,false> r2 = (ac_fixed<RW+1,RW,false>) r;
    #pragma hls_waive CNS

    if (OQ == AC_RND_ZERO || OQ == AC_RND_MIN_INF || OQ == AC_RND_CONV ||  OQ == AC_RND_CONV_ODD) {
      bool rem = (d != 0) || ((z >> 2*RW) != 0);
      if (ZF < (XW-XI)) {
        // max is to used to avoid compilation problems with non pos bitwidth
        const int rbits = AC_MAX((XW-XI)-ZF,1);
        ac_fixed<rbits,-ZF,false> zr = x;
        rem |= !! zr;
      }
      r2[0] =  rem ? 1 : 0;
    }
    ac_math::ac_shift_right(r2, RF, sqrt);
  }

//=========================================================================
// Function: ac_sqrt (for ac_float)
//
// Description:
//    Calculation of square root of real inputs, passed as ac_float
//    variables.
//    + input: ac_float argument
//    + output: ac_float square root
//
//-------------------------------------------------------------------------

  template<bool OR_TF = false, // Override default fractional bitwidth for temp variables?
           int TF_ = 32, // Template argument for fractional bitwidth to override with.
           ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
           int XW, int XI, int XE, ac_q_mode XQ,
           int OW, int OI, int OE, ac_q_mode OQ>
  void ac_sqrt(
    ac_float<XW,XI,XE,XQ> x,
    ac_float<OW,OI,OE,OQ> &sqrt
  )
  {
    #ifdef ASSERT_ON_INVALID_INPUT
    AC_ASSERT(x >= 0, "Negative input not supported.");
    #endif

    // While we don't consider the sign bit for the input mantissa, we still need to add an extra
    // MSB because we might need to left-shift by 1 later.
    ac_fixed<XW, XI, false, XQ> x_mant = x.mantissa();
    // Left shift by 1 if exponent is odd, to aid in denormalization later.
    x_mant <<= x.exp()[0];

    // If OR_TF is set to true, fractional bitwidth used in temp ac_fixed variables (sqrt_mant and sqrt_mant2) is set to TF_.
    // If OR_TF is set to false, fractional bitwidth of the temp ac_fixed variables is set to that of the output mantissa.
    const int TF = OR_TF ? TF_ : OW - OI;
    // Number of integer bits required to store square root of input mantissa = ceil(integer_bits_in_x_mant/2) = ceil(XI/2)
    const int TI = XI%2 == 0 ? XI/2 : (XI + 1)/2;
    ac_fixed<TF + TI, TI, false, TQ> sqrt_mant;

    ac_sqrt(x_mant, sqrt_mant);
    // sqrt(mant*2^exp) = sqrt(mant)*2^(exp/2)
    // Since we've already taken into account odd exponentials by left-shifting the mantissa if needed,
    // we only need to right shift the exponent by 1 to take the "2^(exp/2)" factor into account, even
    // for odd exponential values.
    ac_float<OW, OI, OE, OQ> sqrt_temp(sqrt_mant, x.exp() >> 1, true);

    sqrt = sqrt_temp;
  }

//=========================================================================
// Function: ac_sqrt (for ac_std_float)
//
// Description:
//    Calculation of square root of real inputs, passed as ac_std_float
//    variables.
//    + input: ac_std_float argument
//    + output: ac_std_float square root
//
//    This implementation is a wrapper around the ac_float version, with
//    temporary variables added to enable compatibility between it and the
//    ac_float implementation.
//
//    The template arguments with default values (i.e. OR_TF, TF_ and TQ)
//    apply to the associated ac_float implementation.
//
//-------------------------------------------------------------------------

  template<bool OR_TF = false, // Override default fractional bitwidth for temp variables?
           int TF_ = 32, // Template argument for fractional bitwidth to override with.
           ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
           int XW, int XE, int OW, int OE>
  void ac_sqrt(
    ac_std_float<XW, XE> x,
    ac_std_float<OW, OE> &sqrt
  )
  {
    ac_float<OW - OE + 1, 2, OE> sqrt_ac_fl; // Equivalent ac_float representation for output.
    ac_sqrt<OR_TF, TF_, TQ>(x.to_ac_float(), sqrt_ac_fl); // Call ac_float version.
    ac_std_float<OW, OE> sqrt_temp(sqrt_ac_fl); // Convert output ac_float to ac_std_float.

    // Start of special input handling (modeled after the sqrt() function in math.h)

    // If input is -0.0, output is also set to -0.0.
    if (x == x.zero() && x.signbit()) {
      sqrt_temp.set_signbit(true);
    }
    #ifdef AC_SQRT_NAN_SUPPORTED
    // If input is nan or a negative, non-zero number, the output is set to nan. The sign bit
    // of the output is the same as that of the input.
    else if (x.isnan() || x.signbit()) {
      sqrt_temp = sqrt_temp.nan();
      sqrt_temp.set_signbit(x.signbit());
    }
    #endif

    // End of special input handling.

    sqrt = sqrt_temp;
  }

//=========================================================================
// Function: ac_sqrt (for ac_ieee_float)
//
// Description:
//    Calculation of square root of real inputs, passed as ac_ieee_float
//    variables.
//    + input: ac_ieee_float argument
//    + output: ac_ieee_float square root
//
//    This implementation is a wrapper around the ac_float version, with
//    temporary variables added to enable compatibility between it and the
//    ac_float implementation.
//
//    The template arguments with default values (i.e. OR_TF, TF_ and TQ)
//    apply to the associated ac_float implementation.
//
//-------------------------------------------------------------------------

  template<bool OR_TF = false, // Override default fractional bitwidth for temp variables?
           int TF_ = 32, // Template argument for fractional bitwidth to override with.
           ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
           ac_ieee_float_format Format,
           ac_ieee_float_format outFormat>
  void ac_sqrt(
    ac_ieee_float<Format> x,
    ac_ieee_float<outFormat> &sqrt
  )
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int OW = T_out::width;
    const int OE = T_out::e_width;
    ac_float<OW - OE + 1, 2, OE> sqrt_ac_fl; // Equivalent ac_float representation for output.
    ac_sqrt<OR_TF, TF_, TQ>(x.to_ac_float(), sqrt_ac_fl); // Call ac_float version.
    ac_ieee_float<outFormat> sqrt_temp(sqrt_ac_fl); // Convert output ac_float to ac_ieee_float.

    // Start of special input handling (modeled after the sqrt() function in math.h)

    // If input is -0.0, output is also set to -0.0.
    if (x == x.zero() && x.signbit()) {
      sqrt_temp.set_signbit(true);
    }
    #ifdef AC_SQRT_NAN_SUPPORTED
    // If input is nan or a negative, non-zero number, the output is set to nan. The sign bit
    // of the output is the same as that of the input.
    else if (x.isnan() || x.signbit()) {
      sqrt_temp = sqrt_temp.nan();
      sqrt_temp.set_signbit(x.signbit());
    }
    #endif

    // End of special input handling.

    sqrt = sqrt_temp;
  }
}

#endif
