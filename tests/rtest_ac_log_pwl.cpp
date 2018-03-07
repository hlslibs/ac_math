/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 1.0                                                 *
 *                                                                        *
 *  Release Date    : Wed Mar  7 13:09:26 PST 2018                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 1.0.0                                               *
 *                                                                        *
 *  Copyright , Mentor Graphics Corporation,                     *
 *                                                                        *
 *  All Rights Reserved.                                                  *
 *  
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
// ac_reciprocal_pwl() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_log_pwl.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_log_pwl.h>
using namespace ac_math;

//#include <string>
//#include <stdlib.h>
//#include <fstream>

//==============================================================================

#include <math.h>
#include <iostream>
using namespace std;

//------------------------------------------------------------------------------
// Helper functions for AC_FIXED
// Function to convert test input data (double) into specific type
template<int input_width, int input_int, bool input_S, ac_q_mode input_Q, ac_o_mode input_O>
void double_to_type(
  const double                                                   double_value,
  ac_fixed<input_width, input_int, input_S, input_Q, input_O>   &type_value)
{
  type_value = double_value;
}

// Overloaded function to convert type specific test input data to double
template<int input_width, int input_int, bool input_S, ac_q_mode input_Q, ac_o_mode input_O>
double type_to_double(
  ac_fixed<input_width, input_int, input_S, input_Q, input_O>   &type_value)
{
  return type_value.to_double();
}

//==============================================================================
// Test Design
//   This simple function allows executing the ac_log_pwl() function
//   using multiple data types at the same time. Template parameters are
//   used to configure the bit-widths of the types.

template <int W_IN, int I_IN, bool S_IN, int W_OUT, int I_OUT, bool S_OUT>
void test_ac_log_pwl(
  const ac_fixed<W_IN, I_IN, S_IN>  &in,
  ac_fixed<W_OUT,I_OUT,S_OUT> &loge_out
)
{
  ac_log_pwl(in, loge_out);
}

template <int W_IN, int I_IN, bool S_IN, int W_OUT, int I_OUT, bool S_OUT>
void test_driver()
{
  ac_fixed<W_IN+2,I_IN+1,S_IN> i; // make loop variable slightly larger
  ac_fixed<W_IN,I_IN,S_IN> input;
  ac_fixed<W_IN,I_IN,S_IN> last;
  ac_fixed<W_OUT,I_OUT,S_OUT> loge_out;

  ac_fixed<W_IN,I_IN,S_IN> lower_limit = input.template set_val<AC_VAL_MIN>().to_double();
  ac_fixed<W_IN,I_IN,S_IN> upper_limit = input.template set_val<AC_VAL_MAX>().to_double();
  ac_fixed<W_IN,I_IN,S_IN> step        = input.template set_val<AC_VAL_QUANTUM>().to_double();
  if (step < 0.01) { step = 0.01; }

  printf("TEST: ac_log_pwl() INPUT: ac_fixed<%2d,%2d,%5s,%7s,%7s> OUTPUT: ac_fixed<%2d,%2d,%5s,%7s,%7s>  RESULT: ",
         W_IN,I_IN,(S_IN?"true":"false"),"AC_TRN","AC_WRAP",W_OUT,I_OUT,(S_OUT?"true":"false"),"AC_TRN","AC_WRAP");
//#ifdef DEBUG
//  cout << "    lower_limit  = " << lower_limit << endl;
//  cout << "    upper_limit  = " << upper_limit << endl;
//  cout << "    step         = " << step << endl;
//#endif

  bool passed = true;
  double allowed_error = 0.00160;
  double max_loge_error = 0.0;
  for (i = lower_limit; i < upper_limit; i += step) {
    input = i;
    if (type_to_double(input) == 0) { continue; }

    // call reference log() with fixed-pt value converted back to double
    double expected_loge = log(type_to_double(input));

    // call DUT with fixed-pt value
    test_ac_log_pwl(input,loge_out);

    double actual_loge = type_to_double(loge_out);

    double diff_loge = abs(expected_loge-actual_loge);
    if (diff_loge > allowed_error) {
      cout << "    diff_loge = " << diff_loge << endl;
      passed = false;
      assert(0);
    }
    if (diff_loge > max_loge_error) { max_loge_error = diff_loge; }
  }
  if (passed) { printf("PASSED , max err (%f)\n", max_loge_error); }
  else        { printf("FAILED , max err (%f)\n", max_loge_error); }
}

int main(int argc, char *argv[])
{
  cout << "Testbench start" << endl;
  //          W_IN   I_IN   S_IN   W_OUT  I_OUT  S_OUT         lower   upper    step
  test_driver<  32,    12, false,     32,    16,  true >(); //
  test_driver<  31,    12, false,     32,    16,  true >(); //
  test_driver<  30,    12, false,     32,    16,  true >(); //
  test_driver<  29,    12, false,     32,    16,  true >(); //
  test_driver<  28,    12, false,     32,    16,  true >(); //
  test_driver<  27,    12, false,     32,    16,  true >(); //
  test_driver<  26,    12, false,     32,    16,  true >(); //
  test_driver<  32,    -3, false,     32,    16,  true >(); //
  test_driver<  32,    -2, false,     32,    16,  true >(); //
  test_driver<  32,    -1, false,     32,    16,  true >(); //
  test_driver<  32,     0, false,     32,    16,  true >(); //
  test_driver<  32,     1, false,     32,    16,  true >(); //
  test_driver<  32,     2, false,     32,    16,  true >(); //
  test_driver<  32,     3, false,     32,    16,  true >(); //
  test_driver<  14,    15, false,     32,    16,  true >(); //
  test_driver<  14,    16, false,     32,    16,  true >(); //
  test_driver<  14,    17, false,     32,    16,  true >(); //
  test_driver<  14,    18, false,     32,    16,  true >(); //
  test_driver<  14,    19, false,     32,    16,  true >(); //
  test_driver<  14,    20, false,     32,    16,  true >(); //

#if 0
  cout << "Testbench finished" << endl;
  cout << "max_loge_error      = " << max_loge_error << endl;
  if (max_loge_error > allowed_error) {
    cout << "Error tolerance of " << allowed_error << " percent error exceeded - FAIL" << endl;
    cout << "===================================================" << endl;
    return (-1);
  }
#endif
  cout << "===================================================" << endl;
  return (0);
}
