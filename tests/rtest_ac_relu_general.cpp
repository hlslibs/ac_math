/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 2026.3                                              *
 *                                                                        *
 *  Release Date    : Wed Sep  2 19:51:33 PDT 2026                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 2026.3.0                                            *
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
// =========================TESTBENCH=======================================
// This testbench file contains a stand-alone testbench that exercises the
// ac_relu_general() function using a variety of bit-widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_relu_general.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_relu_general.h>
using namespace ac_math;

// ==============================================================================
// Test Designs
//   This simple function allow executing the ac_relu_general() function.
//   Template parameters are used to configure the bit-widths of the
//   inputs and outputs.

// Test design for real fixed point values.

template <int Wfi, int Ifi, bool Sfi, 
          int outWfi, int outIfi, bool outSfi, 
          int maxWfi, int maxIfi, bool maxSfi,
          int negWfi, int negIfi, bool negSfi,
          int thrWfi, int thrIfi, bool thrSfi>
void test_ac_relu_general_fixed_max(
  const ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP>  &in,
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> &out,
  const ac_fixed<maxWfi, maxIfi, maxSfi, AC_TRN, AC_WRAP> &max_value,
  const ac_fixed<negWfi, negIfi, negSfi, AC_TRN, AC_WRAP> &neg_slope,
  const ac_fixed<thrWfi, thrIfi, thrSfi, AC_TRN, AC_WRAP> &threshold
)
{
  ac_relu_general(in, out, max_value, neg_slope, threshold);
}

template <int Wfi, int Ifi, bool Sfi, 
          int outWfi, int outIfi, bool outSfi, 
          int negWfi, int negIfi, bool negSfi,
          int thrWfi, int thrIfi, bool thrSfi>
void test_ac_relu_general_fixed_no_max(
  const ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP>  &in,
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> &out,
  const ac_fixed<negWfi, negIfi, negSfi, AC_TRN, AC_WRAP> &neg_slope,
  const ac_fixed<thrWfi, thrIfi, thrSfi, AC_TRN, AC_WRAP> &threshold
)
{
  ac_relu_general(in, out, neg_slope, threshold);
}

// ==============================================================================

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

inline double custom_abs(const double x) {
  return (x < 0) ? -x : x;
}

// ==============================================================================
// Function: test_driver_fixed()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_fixed inputs. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   piece-wise linear relu model with the computed relu using a
//   standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi, bool outSfi>
int test_driver_fixed(
  double &cumulative_max_error_fixed,
  const double allowed_error,
  const double threshold_fixed,
  int max_value_type,
  bool details = false
)
{
  bool check_monotonic = true;
  double max_error_fixed = 0.0; // reset for this run

  typedef ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP> T_in;
  typedef ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> T_out;
  T_in input_fixed;
  T_out output;
  T_in max_value;
  T_in neg_slope;
  T_in threshold;

  double lower_limit, upper_limit, step;

  // set ranges and step size for fixed point testbench
  lower_limit = input_fixed.template set_val<AC_VAL_MIN>().to_double();
  upper_limit = input_fixed.template set_val<AC_VAL_MAX>().to_double();
  step        = input_fixed.template set_val<AC_VAL_QUANTUM>().to_double();

  cout << "TEST: ac_relu_general() INPUT: ";
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

  max_value = 6;
  neg_slope = 0.2;
  threshold = 0.1;
  for (double i = lower_limit; i <= upper_limit; i += step) {
    // Set values for input.
    input_fixed = i;
    double expected_value;
    if (max_value_type==1) {
      test_ac_relu_general_fixed_max(input_fixed, output, max_value, neg_slope, threshold);
      expected_value = (input_fixed >= max_value)?max_value.to_double():
                       (input_fixed < threshold)?(neg_slope.to_double()) *(input_fixed-threshold).to_double():
                        input_fixed.to_double();
    }
    else {
      test_ac_relu_general_fixed_no_max(input_fixed, output, neg_slope, threshold);
      expected_value = (input_fixed < threshold)?(neg_slope.to_double()) *(input_fixed-threshold).to_double():
                        input_fixed.to_double();
    }
    
    double actual_value   = output.to_double();

    double this_error;

    // If expected value of output falls below the threshold, calculate absolute error instead of relative
    if (expected_value > threshold_fixed) {
      this_error = custom_abs( (expected_value - actual_value) / expected_value ) * 100.0;
    } else {
      this_error = custom_abs(expected_value - actual_value) * 100.0;
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
          cout << "relu_general output not monotonic at :" << endl; 
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
  cout << "Testing function: ac_relu_general() - Allowed error " << allowed_error << endl;

  // template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi>

  for (int max_value_type = 0; max_value_type < 2; max_value_type++) {
    test_driver_fixed< 10,  3, false, 64, 32, true>(max_error_fixed, allowed_error, threshold, max_value_type);
    test_driver_fixed<  9,  4, false, 64, 32, true>(max_error_fixed, allowed_error, threshold, max_value_type);
    test_driver_fixed<  4,  2, false, 64, 32, true>(max_error_fixed, allowed_error, threshold, max_value_type);
    test_driver_fixed<  3,  4, false, 64, 32, true>(max_error_fixed, allowed_error, threshold, max_value_type);
    test_driver_fixed<  1,  3, false, 61, 33, true>(max_error_fixed, allowed_error, threshold, max_value_type);
    test_driver_fixed<  2,  3, false, 64, 32, true>(max_error_fixed, allowed_error, threshold, max_value_type);
    test_driver_fixed< 18,  4, false, 64, 32, true>(max_error_fixed, allowed_error, threshold, max_value_type);
    test_driver_fixed< 20,  0, false, 64, 32, true>(max_error_fixed, allowed_error, threshold, max_value_type);
    test_driver_fixed< 18,  4,  true, 64, 32,  true>(max_error_fixed, allowed_error, threshold, max_value_type);
    test_driver_fixed< 19,  2,  true, 64, 32,  true>(max_error_fixed, allowed_error, threshold, max_value_type);
    test_driver_fixed< 21, -5,  true, 64, 32,  true>(max_error_fixed, allowed_error, threshold, max_value_type);
    test_driver_fixed< 20,  4,  true, 64, 32,  true>(max_error_fixed, allowed_error, threshold, max_value_type);

    // If error limits on any tested datatype have been crossed, the test has failed
    bool test_fail = (max_error_fixed > allowed_error);

    // Notify the user that the test was a failure.
    if (test_fail) {
        if (max_value_type == 1) 
            cout << "  ac_relu_general - FAILED - Error tolerance(s) exceeded for max_value present" << endl;
        else 
            cout << "  ac_relu_general - FAILED - Error tolerance(s) exceeded for max_value None" << endl;
        cout << "=============================================================================" << endl;
        return -1; 
    } else {
        if (max_value_type == 1) 
            cout << "  ac_relu_general - PASSED for max_value present" << endl;
        else
            cout << "  ac_relu_general - PASSED for max_value None" << endl;
        cout << "=============================================================================" << endl;
    } 
  } 
  return (0);
}
