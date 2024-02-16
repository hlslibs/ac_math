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
 *  Copyright 2018 Siemens                                                *
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
// ac_softplus_pwl() function using a variety of bit-widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_softplus_pwl.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_softplus_pwl.h>
#include <ac_math/ac_log_pwl.h>
using namespace ac_math;

// ==============================================================================
// Test Designs
//   These simple function allow executing the ac_softplus_pwl() function.
//   Template parameters are used to configure the bit-widths of the
//   inputs and outputs.

// Test design for real fixed point values.
template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi>
void test_ac_softplus_pwl_fixed(
  const ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP>  &in,
  ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> &out
)
{
  ac_softplus_pwl(in, out);
}

// ==============================================================================

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

// ==============================================================================
// Function: test_driver_fixed()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_fixed inputs. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   piece-wise linear softplus model with the computed softplus using a
//   standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi>
int test_driver_fixed(
  double &cumulative_max_error_fixed,
  const double allowed_error,
  const double threshold_fixed,
  bool details = false
)
{
  bool check_monotonic = true;
  double max_error_fixed = 0.0; // reset for this run

  ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP> input_fixed;
  typedef ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> T_out;
  T_out output;

  double lower_limit, upper_limit, step;

  // set ranges and step size for fixed point testbench
  lower_limit = input_fixed.template set_val<AC_VAL_MIN>().to_double();
  upper_limit = input_fixed.template set_val<AC_VAL_MAX>().to_double();
  step        = input_fixed.template set_val<AC_VAL_QUANTUM>().to_double();

  cout << "TEST: ac_softplus_pwl() INPUT: ";
  cout.width(38);
  cout << left << input_fixed.type_name();
  cout << "OUTPUT: ";
  cout.width(38);
  cout << left << output.type_name();
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl; // LCOV_EXCL_LINE
    cout << "  Ranges for input types:" << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit = " << lower_limit << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit = " << upper_limit << endl; // LCOV_EXCL_LINE
    cout << "    step        = " << step << endl; // LCOV_EXCL_LINE
  }

  double old_output;
  bool compare = false;

  for (double i = lower_limit; i <= upper_limit; i += step) {
    // Set values for input.
    input_fixed = i;
    test_ac_softplus_pwl_fixed(input_fixed, output);

    double expected_value = ((T_out)log(1 + exp(input_fixed.to_double()))).to_double();
    double actual_value   = output.to_double();

    double this_error;

    // If expected value of output falls below the threshold, calculate absolute error instead of relative
    if (expected_value > threshold_fixed) {
      this_error = abs( (expected_value - actual_value) / expected_value ) * 100.0;
    } else {
      this_error = abs(expected_value - actual_value) * 100.0;
    }


    if (check_monotonic) {
      // MONOTONIC: Make sure that function is monotonic. Compare old value (value of previous iteration) with current value. Since the softplus function we
      // are testing is an increasing function, and our testbench value keeps incrementing, we expect the
      // old value to be lesser than or equal to the current one.

      // This comparison is only carried out once there is an old value to compare with, for the base 2 exponential.
      if (compare) {
        // if by any chance the function output has dropped in value, print out at what point the problem has occured and throw a runtime assertion.
        if (old_output > actual_value) {
          cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl; 
          cout << "softplus output not monotonic at :" << endl; 
          cout << "x = " << input_fixed << endl; 
          cout << "y = " << output << endl; 
          cout << "old_output = " << old_output << endl; 
          assert(false); 
        }
      }
      // Update the old value
      old_output = actual_value;
      // Once an old value has been stored, i.e. towards the end of the first iteration, this value is set to true.
      compare = true;
    }


#ifdef DEBUG
    if (this_error > allowed_error) {
      cout << endl;
      cout << "Error exceeds tolerance" << endl;
      cout << "input_fixed         = " << input_fixed << endl;
      cout << "expected_value = " << expected_value << endl;
      cout << "actual_value   = " << actual_value << endl;
      cout << "this_error     = " << this_error << endl;
      cout << "threshold_fixed     = " << threshold_fixed << endl;
      assert(false);
    }
#endif

    if (this_error > max_error_fixed) {max_error_fixed = this_error;}
  }

  bool passed = (max_error_fixed < allowed_error);

  if (passed) { printf("PASSED , max err (%f)\n", max_error_fixed); }
  else        { printf("FAILED , max err (%f)\n", max_error_fixed); } 

  if (max_error_fixed>cumulative_max_error_fixed) { cumulative_max_error_fixed = max_error_fixed; }

  return 0;
}

int main(int argc, char *argv[])
{
  double max_error_fixed = 0.0;
  double allowed_error = 2.0;
  double threshold = 0.005;

  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_softplus_pwl() - Allowed error " << allowed_error << endl;


  // template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi>
  test_driver_fixed< 10,  3, false, 64, 32>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<  9,  4, false, 64, 32>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<  4,  2, false, 64, 32>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<  4, -2, false, 60, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<  3,  4, false, 64, 32>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<  1,  3, false, 61, 33>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<  2,  3, false, 64, 32>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 18,  4, false, 64, 32>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 20,  0, false, 64, 32>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 18,  3,  true, 64, 32>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 19,  2,  true, 64, 32>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 21, -4,  true, 64, 32>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 20,  3,  true, 64, 32>(max_error_fixed, allowed_error, threshold);


  // If error limits on any tested datatype have been crossed, the test has failed
  bool test_fail = (max_error_fixed > allowed_error);

  // Notify the user that the test was a failure.
  if (test_fail) {
    cout << "  ac_softplus_pwl - FAILED - Error tolerance(s) exceeded" << endl; 
    cout << "=============================================================================" << endl;
    return -1; 
  } else {
    cout << "  ac_softplus_pwl - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }  
  return (0);
}
