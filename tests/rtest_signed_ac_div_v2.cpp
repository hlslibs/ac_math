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
// =========================TESTBENCH=======================================
// This testbench file contains a stand-alone testbench that exercises the
// ac_div() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_signed_ac_div_v2.cpp -o design
//   ./design

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <bitset>

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_div_v2.h>

using namespace std;

const int numWint = 10;
const int denWint = 5;
const int padWint= 2;

// Pad width minimum is 2 for Signed values and maximum is divisor width.
// For Unsigned minium is 1 and max is divisor width.

const int quoWint = numWint-denWint+1+padWint;
const int remWint = numWint;

// Define types for ac_int<> signed usage
typedef ac_int<numWint, true> s_real_num_int_type;
typedef ac_int<denWint, true> s_real_den_int_type;
typedef ac_int<quoWint, true> s_real_quo_int_type;
typedef ac_int<remWint, true> s_real_rem_int_type;

#pragma hls_design top
void rtest_ac_div(
  const s_real_num_int_type &in1,
  const s_real_den_int_type &in2,
  s_real_quo_int_type &out1,
  s_real_rem_int_type &out2
)
{
  ac_div_v2<numWint,denWint,quoWint,remWint,padWint>(in1,in2,out1,out2);
}

int main(int argc, char *argv[])
{

//================== 	Variable Declared =================================
  bool ERROR=false;
  s_real_num_int_type	Dividend,max_num,min_num;
  s_real_den_int_type	Divisor,max_den,min_den,mask=0;
  s_real_quo_int_type	DUT_Quotient;
  s_real_rem_int_type	DUT_Remainder;
  int32 REF_Quotient,   REF_Remainder;
  mask[padWint-1]=1;
//================================================
  max_den.template set_val<AC_VAL_MAX>();
  max_num.template set_val<AC_VAL_MAX>();
  min_num.template set_val<AC_VAL_MIN>();
  min_den.template set_val<AC_VAL_MIN>();

#ifdef DEBUG

  cout << "Max value of the numerator   : " << max_num << " [ " << bitset<numWint>(max_num) << " ] "<< endl;
  cout << "Max value of the denominator : " << max_den << " [ " << bitset<denWint>(max_den) << " ] "<< endl;
  cout << "Min value of the numerator   : " << min_num << " [ " << bitset<numWint>(min_num) << " ] "<< endl;
  cout << "Min value of the denominator : " << min_den << " [ " << bitset<denWint>(min_den) << " ] "<< endl;;
  cout << "Range of Divisor values not supported by synthesized design are from : " << (min_den/mask)-1 <<	" to " << (max_den/mask)+1 << endl;

#endif

  for (double i=min_den; i<=(double)max_den; i++) {
    for (double j =min_num; j<=(double)max_num; j++) {
      // cout << "----------------------------------------" << endl;

      Dividend=j;
      Divisor=i;

      if (Divisor==0) {
        // cout << "Skipping the Zero Value input " << endl;
        continue;
      }

      else if (Divisor>(min_den/mask)-1 && Divisor < (max_den/mask)+1) {
        // cout << "Skipping the Divisor values out of range " << endl;
        continue;
      }

      /* Feed the random values into the design */
      rtest_ac_div(Dividend,Divisor,DUT_Quotient,DUT_Remainder);

      /* Golden Reference of the division model */
      REF_Quotient=Dividend/Divisor;												 // Reference Version
      REF_Remainder=Dividend%Divisor;

#ifdef DEBUG
      /* Debug System Outs */
      cout << "Dividend : "  << Dividend       << " [ " << bitset<numWint>(Dividend) << " ] "<< endl;
      cout << "Divisor  : "  << Divisor        << " [ " << bitset<denWint>(Divisor) << " ] "<< endl;
      cout << "Quotient  : " << "New Design: " <<  DUT_Quotient  << " || ReF Value:  " << REF_Quotient  << endl;
      cout << "Remainder : " << "New Design: " <<  DUT_Remainder << " || ReF Value:  " << REF_Remainder << endl;
#endif 

      /* Compare the DUT values against Reference value */
      if (DUT_Quotient.to_int() != REF_Quotient) {
        cout << "ERROR !!!!!!!!!!! in Quotient " << endl;
        ERROR= true;
		return 0;
      }
      if (DUT_Remainder.to_int() != REF_Remainder) {
        cout << "ERROR !!!!!!!!!!! in Remainder " << endl;
        ERROR = true;
		return 0;
      }
    }
  }
  if (!ERROR) {
    cout << " No Error , All Clear " << endl;
  }

  return 0;
}
