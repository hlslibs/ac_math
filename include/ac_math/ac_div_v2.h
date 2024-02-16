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
 *  Copyright  Siemens                                                *
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
/******************************************************************************************
 File: ac_div_v2.h for Unsigned and Signed integers
******************************************************************************************/

/******************************************************************************************

 Port Description:
      + inputs: dividend, divisor
      + outputs: quotient and remainder


 Design Parameters
 		+ NW: Numerator Width
		+ DW: Denominator Width
		+ QW: Quotient Width
		+ RW: Remainder Width
		+ PW: Pad Width (value to be withing 1-> DW), this determines range of the divisor values the minimum value of the divisor supported.

		Only divisor values above 2^(DW-PW) are supported
		if PW=1 , area shall be least, if PW = DW, area shall be maximum.

    * Integer Unsigned
      void ac_div(ac_int<NW,false> dividend, ac_int<DW,false> divisor, ac_int<QW,false> &quotient, ac_int<RW,false> &remainder)
	
 Usage:
    A sample testbench for unsigned inputs and output look like this:

    #include <ac_math/ac_div_v2.h>
    using namespace ac_math;

	enum {
	  MAX_NW =  5,
	  MAX_DW =  3,
	  QUO_W  =  5,
	  REM_W	 =  5,
	  PAD_W  = 1, // Value to be padded to the numerator width, that varies the range of value divisor can accept.
	};

    typedef ac_int<MAX_NW,false> num_type;
    typedef ac_int<MAX_DW,false> den_type;
    typedef ac_int<QUO_W,false> quo_type;
    typedef ac_int<REM_W,false> rem_type

    #pragma hls_design top
    void project(
      const num_type &num,
      const den_type &den,
      quo_type &quo,
	rem_type &rem
    )
    {

      //to replicate: quo = num / den and rem= num%den	  
      ac_div<NW,DW,QW,RW>(num, den, quo, rem);
    }

    #ifndef __SYNTHESIS__
    #include <mc_scverify.h>

    CCS_MAIN(int arg, char **argc)
    {
      num_type num = 32;
      den_type den = 6;
      quo_type quo;
	rem_type rem;
      CCS_DESIGN(project)(num, den, quo,rem);
      CCS_RETURN (0);
    }
    #endif

********************************************************************************************/

#ifndef _INCLUDED_AC_DIV_V2_H_
#define _INCLUDED_AC_DIV_V2_H_

#ifndef __cplusplus
#error C++ is required to include this header file
#endif

// Include headers for data types supported by these implementations
#include <ac_int.h>

// Include header file for ac_shift functions
#include <ac_math/ac_shift.h>

template< int NW, int DW, int QW, int RW, int PW>
void ac_div_v2(
  ac_int<NW,false> dividend,
  ac_int<DW,false> divisor,
  ac_int<QW,false> &quotient,
  ac_int<RW,false> &remainder
)
{
  const unsigned int int_nw=NW+PW;
  const unsigned int int_qw=int_nw-DW;

  ac_int<DW,false> min_div_value=0;
  min_div_value[DW-PW]=1;
  ac_int<int_nw+1,false> temp_numer;
  ac_int<DW+1,true > rem;
  ac_int<int_qw,  false> quot;
  ac_int<int_nw+1, false> const_mask;

  const_mask.template set_val<AC_VAL_MAX>();

  AC_ASSERT(!!divisor, "Division by zero.");
  AC_ASSERT(!!(divisor>= min_div_value), "Divisor is Below the range supported by the synthesized design, check for Padding width");

  temp_numer = (ac_int<int_nw+1,false>) dividend;
  quot = 0;
  #pragma hls_unroll yes
  /*UNROLL it to save on latency, Roll it save on Area*/
  AC_DIV_CORE_LOOP :for (ac_int<int_qw,false> cnt=0; cnt<int_qw; cnt++) {
    rem = (ac_int<int_nw+1,true>) ((temp_numer<<1)>>int_qw) - (ac_int<int_nw+1,true>) divisor;
    quot = (rem < 0) ? ((quot<<1)|0) : ((quot<<1)|1);
    temp_numer = (rem < 0) ? (temp_numer<<1) : (((ac_int<int_nw+1,false>) rem)<<int_qw)|((temp_numer<<1)&(const_mask>>(DW+1)));
  }
  remainder= (ac_int<RW,false>)(rem>0)? (ac_int<int_nw+1,false>)rem : (ac_int<int_nw+1,false>)(temp_numer>>int_qw);
  quotient = (ac_int<QW,false>)quot;
}

/******************************************************************************************

 Port Description:
      + inputs: dividend, divisor
      + outputs: quotient and remainder


 Design Parameters
 		+ NW: Numerator Width
		+ DW: Denominator Width
		+ QW: Quotient Width
		+ RW: Remainder Width
		+ PW: Pad Width (value to be withing 2 to DW), this determines range of the divisor values the minimum value of the divisor supported.

		Range of Divisor Values supported for signed input/output is (-2^(DW-1)/(2^(PW-1) to (2^(DW-1) - 1)/(2^(PW-1)
		if PW=2 , area shall be least, if PW = DW, area shall be maximum.

	* Integer Signed  
      void ac_div(ac_int<NW,true> dividend, ac_int<DW,true> divisor, ac_int<QW,true> &quotient, ac_int<RW,true> &remainder)


 Usage:
    A sample testbench for unsigned inputs and output look like this:

    #include <ac_math/ac_div_v2.h>
    using namespace ac_math;

	enum {
	  MAX_NW =  5,
	  MAX_DW =  3,
	  QUO_W  =  5,
	  REM_W	 =  5,
	  PAD_W  = 1, // Value to be padded to the numerator width, that varies the range of value divisor can accept.
	};

    typedef ac_int<MAX_NW,true> num_type;
    typedef ac_int<MAX_DW,true> den_type;
    typedef ac_int<QUO_W,true> quo_type;
    typedef ac_int<REM_W,true> rem_type

    #pragma hls_design top
    void project(
      const num_type &num,
      const den_type &den,
      quo_type &quo,
	rem_type &rem
    )
    {
      //to replicate: quo = num / den and rem= num%den	  
      ac_div<NW,DW,QW,RW>(num, den, quo, rem);
    }

    #ifndef __SYNTHESIS__
    #include <mc_scverify.h>

    CCS_MAIN(int arg, char **argc)
    {
      num_type num = 32;
      den_type den = 6;
      quo_type quo;
	  rem_type rem;
      CCS_DESIGN(project)(num, den, quo,rem);
      CCS_RETURN (0);
    }
    #endif

********************************************************************************************/

template< int NW, int DW, int QW, int RW, int PW >
void ac_div_v2(
  ac_int<NW,true> dividend,
  ac_int<DW,true> divisor,
  ac_int<QW,true> &quotient,
  ac_int<RW,true> &remainder
)
{
  AC_ASSERT(!!divisor, "Division by zero.");
  bool neg_dividend = dividend < 0;
  ac_int<NW,false> uN = neg_dividend ? (ac_int<NW,false>) -dividend : (ac_int<NW,false>) dividend;
  bool neg_divisor = divisor < 0;
  ac_int<DW,false> uD = neg_divisor ? (ac_int<DW,false>) -divisor : (ac_int<DW,false>) divisor;
  ac_int<QW,false> uQ;
  ac_int<RW,false> uR;
  ac_div_v2<NW,DW,QW,RW,PW>(uN, uD, uQ, uR);
  // std::cout << "Resulting Quotient is: " <<uQ << " [ " << std::bitset<QW>(uQ) << " ] "<< std::endl;

  ac_int<QW,true> quotient_temp = neg_dividend == neg_divisor ? (ac_int<QW,true>) uQ : (ac_int<QW,true>) -uQ;

  quotient = quotient_temp;

  // std::cout << divisor << std::endl;

  ac_int<RW,true> rem = neg_dividend ? (ac_int<DW,true>) -uR : (ac_int<DW,true>) uR;
  remainder = rem;
}

#endif

