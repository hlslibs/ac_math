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
// ac_pow2_pwl() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_pow2_pwl.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_pow_pwl.h>
using namespace ac_math;

//==============================================================================
// Test Design
//   This simple function allows executing the ac_pow2_pwl() function.
//   Template parameters are used to configure the bit-widths of the
//   ac_fixed inputs.

template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi, bool outSfi>
void test_ac_pow2_pwl(
  const ac_fixed<   Wfi,    Ifi,    Sfi, AC_TRN, AC_WRAP> &in,
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> &out_pow2
)
{
  ac_pow2_pwl(in, out_pow2);
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
//   widths of AC datatypes. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   piece-wise linear power model with the computed power using a
//   standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi, bool outSfi>
int test_driver(
  double &cumulative_max_error_pow2,
  const double allowed_error,
  const double threshold,
  bool details = false
)
{
  bool passed = true;
  bool check_monotonic = true;
  double max_error_pow2 = 0.0; // reset for this run

  ac_fixed<   Wfi,    Ifi,    Sfi, AC_TRN, AC_WRAP> input;
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> output_pow2;

  typedef ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> T_out;

  double lower_limit, upper_limit, step;

  // set ranges and step size for fixed point testbench
  lower_limit = input.template set_val<AC_VAL_MIN>().to_double();
  upper_limit = input.template set_val<AC_VAL_MAX>().to_double();
  step        = input.template set_val<AC_VAL_QUANTUM>().to_double();

  printf("TEST: ac_pow2_pwl() INPUT: ac_fixed<%2d,%2d,%5s,%7s,%7s> OUTPUT: ac_fixed<%2d,%2d,%5s,%7s,%7s>  RESULT: ",
         Wfi,Ifi,(Sfi?"true":"false"),"AC_TRN","AC_WRAP",outWfi,outIfi,(outSfi?"true":"false"),"AC_TRN","AC_WRAP");

  // Dump the test details
  if (details) {
    cout << endl;
    cout << "  Ranges for input types:" << endl;
    cout << "    lower_limit          = " << lower_limit << endl;
    cout << "    upper_limit          = " << upper_limit << endl;
    cout << "    step                 = " << step << endl;
  }

  double old_output_pow2;
  bool compare_pow2 = false;

  for (double i = lower_limit; i < upper_limit; i += step) {
    //Set values for input.
    input = i;
    test_ac_pow2_pwl(input, output_pow2);

    double expected_value_pow2 = ((T_out)pow(2, input.to_double())).to_double();
    double actual_value_pow2   = output_pow2.to_double();

    double this_error_pow2;

    //If expected value of either output falls below the threshold, calculate absolute error instead of relative
    if (expected_value_pow2 > threshold) {this_error_pow2 = abs( (expected_value_pow2 - actual_value_pow2) / expected_value_pow2 ) * 100.0;}
    else {this_error_pow2 = abs(expected_value_pow2 - actual_value_pow2) * 100.0;}

    if (check_monotonic) {
      // MONOTONIC: Make sure that function is monotonic. Compare old value (value of previous iteration) with current value. Since the exponential function we
      //are testing is an increasing function, and our testbench value keeps incrementing or remains the same (in case of saturation), we expect the
      //old value to be lesser than or equal to the current one.

      //This comparison is only carried out once there is an old value to compare with, for the base 2 exponential.
      if (compare_pow2) {
        //if by any chance the function output has dropped in value, print out at what point the problem has occured and throw a runtime assertion.
        if (old_output_pow2 > actual_value_pow2) {
          cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl;
          cout << "pow2 output not monotonic at :" << endl;
          cout << "x = " << input << endl;
          cout << "y = " << output_pow2 << endl;
          cout << "old_output_pow2 = " << old_output_pow2 << endl;
          //assert(false);   //Uncomment if you want the program to stop whenever monotonicity is violated.
        }
      }
      //Update the old value
      old_output_pow2 = actual_value_pow2;
      //Once an old value has been stored, i.e. towards the end of the first iteration, this value is set to true.
      compare_pow2 = true;
    }

    if (this_error_pow2 > max_error_pow2) {max_error_pow2 = this_error_pow2;}
  }
  if (passed) { printf("PASSED , max err (%f pow2)\n", max_error_pow2); }
  else        { printf("FAILED , max err (%f pow2)\n", max_error_pow2); }

  if (max_error_pow2>cumulative_max_error_pow2) { cumulative_max_error_pow2 = max_error_pow2; }

  return 0;
}


int main(int argc, char *argv[])
{
  double max_error_pow2 = 0;
  double allowed_error = 0.5;
  double threshold = 0.005;
  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_pow2_pwl() - Allowed error " << allowed_error << endl;

  // template <int Wfi, int Ifi, int Sfi>
  test_driver< 12,  3,  true, 64, 32, false>(max_error_pow2, allowed_error, threshold);
  test_driver<  4,  2,  true, 64, 32, false>(max_error_pow2, allowed_error, threshold);
  test_driver<  3,  5,  true, 64, 32, false>(max_error_pow2, allowed_error, threshold);
  test_driver<  4, -2,  true, 64, 32, false>(max_error_pow2, allowed_error, threshold);
  test_driver<  3,  5,  true, 64, 32, false>(max_error_pow2, allowed_error, threshold);
  test_driver<  2,  5,  true, 64, 32, false>(max_error_pow2, allowed_error, threshold);
  test_driver< 16,  5,  true, 64, 32, false>(max_error_pow2, allowed_error, threshold);
  test_driver< 12,  4, false, 64, 32, false>(max_error_pow2, allowed_error, threshold);
  test_driver<  4,  2, false, 64, 32, false>(max_error_pow2, allowed_error, threshold);
  test_driver<  4, -2, false, 60, 30, false>(max_error_pow2, allowed_error, threshold);
  test_driver<  3,  4, false, 64, 32, false>(max_error_pow2, allowed_error, threshold);
  test_driver<  1,  5, false, 61, 33, false>(max_error_pow2, allowed_error, threshold);
  test_driver<  2,  5, false, 64, 32, false>(max_error_pow2, allowed_error, threshold);
  test_driver< 16,  4, false, 64, 32, false>(max_error_pow2, allowed_error, threshold);
  test_driver< 16,  0, false, 64, 32, false>(max_error_pow2, allowed_error, threshold);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum errors observed across all bit-width variations:" << endl;
  cout << "    max_error_pow2 = " << max_error_pow2 << endl;

  // If error limits on any tested datatype have been crossed, the test has failed
  bool test_fail = (max_error_pow2 > allowed_error);

  // Notify the user that the test was a failure.
  if (test_fail) {
    cout << "  ac_pow2_pwl - FAILED - Error tolerance(s) exceeded" << endl;
    cout << "=============================================================================" << endl;
    return -1;
  } else {
    cout << "  ac_pow2_pwl - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }
  return 0;
}




































