/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 1.0                                                 *
 *                                                                        *
 *  Release Date    : Fri Mar  2 14:27:58 PST 2018                        *
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
// ac_sincos() function using a variety of bit-widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_sincos_lut.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_sincos.h>
using namespace ac_math;

//==============================================================================
// Test Design
//   This simple function allows executing the ac_sincos() function.
//   Template parameters are used to configure the bit-widths of the types.

template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi, bool outSfi>
void test_ac_sincos(
  const                ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP>    &in,
  ac_complex<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> >   &out
)
{
  ac_sincos(in, out);
}

//==============================================================================

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

//==============================================================================
// Function: test_driver()
// Description: A templatized function that can be configured for certain bit-
//   widths of the fixed point AC datatype. It uses the type information to
//   iterate through a range of valid values on that type in order to compare
//   the precision of the lookup table sinecosine model with the sinecosine using
//   the standard C math library. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <int Wfi, int Ifi, int Sfi, int outWfi, int outIfi, bool outSfi>
int test_driver(
  double &cummulative_max_error_sine,
  double &cummulative_max_error_cosine,
  const double allowed_error,
  bool details = false
)
{
  bool passed = true;
  double max_error_sine   = 0.0; // reset for this run
  double max_error_cosine = 0.0; // reset for this run

  ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP>                input_angle;
  ac_complex<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> >  output;

  double lower_limit, upper_limit, step;

  // set ranges and step size for the testbench
  lower_limit   = input_angle.template set_val<AC_VAL_MIN>().to_double();
  upper_limit   = input_angle.template set_val<AC_VAL_MAX>().to_double();
  step          = input_angle.template set_val<AC_VAL_QUANTUM>().to_double();

  printf("TEST: ac_sincos_lut() INPUT: ac_fixed<%2d,%2d,%5s,%7s,%7s> OUTPUT: ac_fixed<%2d,%2d,%5s,%7s,%7s>  RESULT: ",
         Wfi,Ifi,(Sfi?"true":"false"),"AC_TRN","AC_WRAP",outWfi,outIfi,(outSfi?"true":"false"),"AC_TRN","AC_WRAP");

  // Dump the test details
  if (details) {
    cout << endl;
    cout << "  Ranges for input types:" << endl;
    cout << "    lower_limit    = " << lower_limit << endl;
    cout << "    upper_limit    = " << upper_limit << endl;
    cout << "    step           = " << step << endl;
  }

  for (double i = lower_limit; i < upper_limit; i += step) {
    // Set values for input.
    input_angle = i;
    test_ac_sincos(input_angle, output);

    double expected_sine_value   = sin(input_angle.to_double()*2*M_PI);
    double actual_sine_value     = output.i().to_double();
    double expected_cosine_value = cos(input_angle.to_double()*2*M_PI);
    double actual_cosine_value   = output.r().to_double();
    double this_error_sine, this_error_cosine;

    //When the C++ math library is supposed to return 0 for the cases when sine and cosine
    //values are 0, the C++ math library returns a small non-zero value. Hence for these cases
    //the values are rounded off to zero.
    if (fmod(abs(input_angle.to_double()),0.25) == 0 && fmod(abs(input_angle.to_double()),0.5) != 0)
    { expected_cosine_value = 0.0; }

    if (fmod(abs(input_angle.to_double()),0.5) == 0)
    { expected_sine_value = 0.0; }

    if (expected_sine_value != 0)
    { this_error_sine   = abs( ( expected_sine_value - actual_sine_value) / expected_sine_value ) * 100.0; }
    else
    { this_error_sine   = abs(expected_sine_value - actual_sine_value) * 100.0; }

    if (expected_cosine_value != 0)
    { this_error_cosine = abs( ( expected_cosine_value - actual_cosine_value) / expected_cosine_value ) * 100.0; }
    else
    { this_error_cosine = abs(expected_cosine_value - actual_cosine_value) * 100.0; }

    if (this_error_sine > max_error_sine) {max_error_sine = this_error_sine;}

    if (this_error_cosine > max_error_cosine) {max_error_cosine = this_error_cosine;}

  }

  if (passed) { printf("PASSED , max err (%f sin) (%f cos)\n", max_error_sine, max_error_cosine); }
  else        { printf("FAILED , max err (%f sin) (%f cos)\n", max_error_sine, max_error_cosine); }

  if (max_error_sine>cummulative_max_error_sine) { cummulative_max_error_sine = max_error_sine; }
  if (max_error_cosine>cummulative_max_error_cosine) { cummulative_max_error_cosine = max_error_cosine; }

  return 0;
}


int main(int argc, char *argv[])
{
  double max_error_sine = 0, max_error_cosine = 0;
  double allowed_error = 0.1;
  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_sincos() - Allowed error " << allowed_error << endl;

  // template <int Wfi, int Ifi, int Sfi>
  test_driver< 12,  1, false, 24, 2, true>(max_error_sine, max_error_cosine, allowed_error);
  test_driver<  2,  0, false, 27, 2, true>(max_error_sine, max_error_cosine, allowed_error);
  test_driver< 12,  0,  true, 24, 3, true>(max_error_sine, max_error_cosine, allowed_error);
  test_driver<  4,  1,  true, 25, 3, true>(max_error_sine, max_error_cosine, allowed_error);
  test_driver<  5, -2,  true, 30, 4, true>(max_error_sine, max_error_cosine, allowed_error);
  test_driver<  3, -2,  true, 28, 2, true>(max_error_sine, max_error_cosine, allowed_error);
  test_driver<  2, -2,  true, 32, 2, true>(max_error_sine, max_error_cosine, allowed_error);
  test_driver< 11,  1,  true, 34, 2, true>(max_error_sine, max_error_cosine, allowed_error);
  test_driver<  8, -3, false, 25, 2, true>(max_error_sine, max_error_cosine, allowed_error);
  test_driver<  9, -3, false, 26, 2, true>(max_error_sine, max_error_cosine, allowed_error);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum errors observed across all bit-width variations:" << endl;
  cout << "    max_error_sine       = " << max_error_sine   << endl;
  cout << "    max_error_cosine     = " << max_error_cosine << endl;

  // If error limits on any tested datatype have been crossed, the test has failed
  bool test_fail = (max_error_sine > allowed_error) || (max_error_cosine > allowed_error);

  // Notify the user that the test was a failure.
  if (test_fail) {
    cout << "  ac_sincos - FAILED - Error tolerance(s) exceeded" << endl;
    cout << "=============================================================================" << endl;
    return -1;
  } else {
    cout << "  ac_sincos - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }
  return 0;
}


