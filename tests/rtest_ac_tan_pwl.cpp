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
// =========================TESTBENCH=======================================
// This testbench file contains a stand-alone testbench that exercises the
// ac_tan_pwl() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_tan_pwl.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_tan_pwl.h>
using namespace ac_math;

// ==============================================================================
// Test Designs
//   These simple function allow executing the ac_tan_pwl() function.
//   Template parameters are used to configure the bit-widths of the
//   inputs and outputs.

// Test design for real fixed point values.
template <int Wfi, int Ifi, int outWfi, int outIfi>
void test_ac_tan_pwl_fixed(
  const ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP> &in,
  ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> &tan_out
)
{
  tan_out = ac_tan_pwl<ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> >(in);
}

// Test design for real floating point values.
template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
void test_ac_tan_pwl_float(
  const ac_float<Wfl, Ifl, Efl, AC_TRN> &in,
  ac_float<outWfl, outIfl, outEfl, AC_TRN> &tan_out
)
{
  tan_out = ac_tan_pwl<ac_float<outWfl, outIfl, outEfl, AC_TRN> >(in);
}

// ------------------------------------------------------------------------------
// Helper function for absolute value calculation. This can avoid any naming conflicts
// with other absolute value functions.

double abs_double(double x)
{
  return x >= 0 ? x : -x;
}

// ==============================================================================

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

// ==============================================================================
// Functions: test_driver functions
// Description: Templatized functions that can be configured for certain bit-
//   widths of AC datatypes. They use the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   piece-wise linear tangent model with the computed tangent using a C++
//   math library with the standard C++ double type. The maximum error for each
//   type is accumulated in variables defined in the calling function.

// ===============================================================================
// Function: test_driver_fixed()
// Description: test_driver function for ac_fixed inputs and outputs.

template <int Wfi, int Ifi, int outWfi, int outIfi>
int test_driver_fixed(
  double &cumulative_max_error_fixed,
  const double allowed_error,
  const double threshold,
  const double upper_limit_degrees,
  bool details = false
)
{
  ac_fixed<Wfi + 2, Ifi + 1, false, AC_TRN, AC_WRAP> i; // make loop variable slightly larger
  ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP> input;
  ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP> last;
  typedef ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> T_out;
  T_out tan_out;

  // set ranges and step size for testbench
  ac_fixed<Wfi, Ifi> lower_limit = input.template set_val<AC_VAL_MIN>();
  double upper_limit_type = input.template set_val<AC_VAL_MAX>().to_double();
  double upper_limit_rad  = upper_limit_degrees * M_PI / 180;
  ac_fixed<Wfi, Ifi> upper_limit = upper_limit_type < upper_limit_rad ? upper_limit_type : upper_limit_rad;
  ac_fixed<Wfi, Ifi> step        = input.template set_val<AC_VAL_QUANTUM>();

  cout << "TEST: ac_tan_pwl() INPUT: ";
  cout.width(38);
  cout << left << input.type_name();
  cout << "upper_limit_degrees = " << upper_limit_degrees << " ";
  cout << "OUTPUT: ";
  cout.width(38);
  cout << left << tan_out.type_name();
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl; // LCOV_EXCL_LINE
    cout << "  Ranges for input types:" << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit = " << lower_limit << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit = " << upper_limit << endl; // LCOV_EXCL_LINE
    cout << "    step        = " << step << endl; // LCOV_EXCL_LINE
  }

  bool passed = true;
  double max_tan_error = 0.0;

  bool check_monotonic = true;
  bool compare_tan = false;
  double old_output_tan;

  for (i = lower_limit; i <= upper_limit; i += step) {
    // Set values for input.
    input = i;

    // call reference tan() with fixed-pt value converted back to double
    // an additional step of typecasting is required in order to perform
    // quantization on the expected output.
    double expected_tan = ((T_out)tan(input.to_double())).to_double();

    // call DUT with fixed-pt value
    test_ac_tan_pwl_fixed(input, tan_out);

    double actual_tan = tan_out.to_double();
    double this_error_tan;

    // If expected value of either output falls below the threshold, calculate absolute error instead of relative
    if (expected_tan > threshold) {this_error_tan = abs_double( (expected_tan - actual_tan) / expected_tan ) * 100.0;}
    else {this_error_tan = abs_double(expected_tan - actual_tan) * 100.0;}

    if (check_monotonic) {
      // MONOTONIC: Make sure that function is monotonic. Compare old value (value of previous iteration) with current value. Since the tangent function we
      // are testing is an increasing function, and our testbench value keeps incrementing or remains the same (in case of saturation), we expect the
      // old value to be lesser than or equal to the current one.

      // This comparison is only carried out once there is an old value to compare with
      if (compare_tan) {
        // Figuring out what the normalized value was for the input is a good way to figure out where the discontinuity occured w.r.t. the PWL segments.
        ac_fixed<Wfi, 0, false, AC_TRN, AC_WRAP> norm_input;
        ac_normalize(input, norm_input);
        // if by any chance the function output has dropped in value, print out at what point the problem has occured and throw a runtime assertion.
        if (old_output_tan > actual_tan) {
          cout << endl; // LCOV_EXCL_LINE
          cout << "  tan output not monotonic at :" << endl; // LCOV_EXCL_LINE
          cout << "  x = " << input << endl; // LCOV_EXCL_LINE
          cout << "  y = " << tan_out << endl; // LCOV_EXCL_LINE
          cout << "  old_output_tan = " << old_output_tan << endl; // LCOV_EXCL_LINE
          assert(false); // LCOV_EXCL_LINE
        }
      }
      // Update the old value
      old_output_tan = actual_tan;
      // Once an old value has been stored, i.e. towards the end of the first iteration, this value is set to true.
      compare_tan = true;
    }

#ifdef DEBUG
    double input_degrees = input.to_double() * 180 / M_PI;
    if (this_error_tan > allowed_error) {
      cout << endl;
      cout << "  Error exceeds tolerance" << endl;
      cout << "  input          = " << input << endl;
      cout << "  input_degrees  = " << input_degrees << endl;
      cout << "  expected_tan   = " << expected_tan << endl;
      cout << "  actual_tan     = " << actual_tan << endl;
      cout << "  this_error_tan = " << this_error_tan << endl;
      cout << "  threshold      = " << threshold << endl;
      assert(false);
    }
#endif

    if (this_error_tan > max_tan_error) { max_tan_error = this_error_tan; }
  }

  if (max_tan_error > cumulative_max_error_fixed) { cumulative_max_error_fixed = max_tan_error; }

  passed = (max_tan_error < allowed_error);

  if (passed) { printf("PASSED , max err (%f)\n", max_tan_error); }
  else        { printf("FAILED , max err (%f)\n", max_tan_error); } // LCOV_EXCL_LINE

  return 0;
}

// ===============================================================================
// Function: test_driver_float()
// Description: test_driver function for ac_float inputs and outputs.

template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
int test_driver_float(
  double &cumulative_max_error_float,
  const double allowed_error_float,
  const double threshold,
  const double upper_limit_degrees,
  bool details = false
)
{
  double max_error_float = 0.0; // Reset with every function call.

  typedef ac_float<Wfl, Ifl, Efl, AC_TRN> T_in;
  typedef ac_float<outWfl, outIfl, outEfl, AC_TRN> T_out;

  // Since ac_float values are normalized, the bit adjacent to the sign bit in the mantissa
  // will always be set to 1. We will hence cycle through all the bit patterns that correspond to the last (Wfl - 2)
  // bits in the mantissa.
  ac_int<Wfl - 2, false> sample_mantissa_slc;
  // Set the lower limit, upper limit and step size of the test iterations.
  ac_int<Wfl - 2, false> lower_limit_it = 0;
  ac_int<Wfl - 2, false> upper_limit_it = sample_mantissa_slc.template set_val<AC_VAL_MAX>().to_double();
  ac_int<Wfl - 2, false> step_it = 1; // Since sample_mantissa_slc is an integer.

  cout << "TEST: ac_tan_pwl() INPUT: ";
  cout.width(38);
  cout << left << T_in::type_name();
  cout << "OUTPUT: ";
  cout.width(38);
  cout << left << T_out::type_name();
  cout << "RESULT: ";

  // sample_exponent_array stores all values of exponent to be tested.
  const int exp_arr_size = 2*(Efl - 1) + 3;
  ac_int<Efl, true> sample_exponent_array[exp_arr_size];
  ac_int<Efl, true> sample_exponent_value;

  // The first element of the array is the minimum exponent value, the middle element is a zero exponent, and
  // the last element is the maximum possible value.
  sample_exponent_array[0]                = sample_exponent_value.template set_val<AC_VAL_MIN>();
  sample_exponent_array[Efl]              = 0;
  sample_exponent_array[exp_arr_size - 1] = sample_exponent_value.template set_val<AC_VAL_MAX>();

  // All the other elements are set to values that correspond to a one-hot encoding scheme, in which only one
  // bit of the absolute value of the exponent is set to one. Both negative and positive values are encoded this way.
  for (int i = (Efl - 2); i >= 0; i--) {
    sample_exponent_value = 0;
    sample_exponent_value[i] = 1;
    sample_exponent_array[Efl + i + 1] = sample_exponent_value;
    sample_exponent_array[Efl - i - 1] = -sample_exponent_value;
  }

  // Dump the test details
  if (details) {
    cout << endl;
    cout << "  Ranges for testing iterations:" << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit_it       = " << lower_limit_it << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit_it       = " << upper_limit_it << endl; // LCOV_EXCL_LINE
    cout << "    step_it              = " << step_it << endl; // LCOV_EXCL_LINE
    cout << "    allowed_error_float  = " << allowed_error_float << endl; // LCOV_EXCL_LINE
  }

  const double upper_limit_rad = upper_limit_degrees* M_PI/180;

  for (int i = 0; i < exp_arr_size; i++) {
    // For a particular exponent value, go through every possible value that can be represented by the mantissa.
    // The iteration variable has a bitwidth that is 1 higher (bitwidth = Wfl - 1) than the slice of the mantissa
    // we'll be changing from one iteration to the other (bitwidth of slice = Wfl - 2), to ensure that the iteration
    // variable does not overflow.
    for (ac_int<Wfl - 1, false> mant_i = lower_limit_it; mant_i <= upper_limit_it; mant_i += step_it) {
      ac_fixed<Wfl, Ifl, true> input_mant = 0; // Initializing this variable avoids possible compiler warnings.
      // Set the sign bit to zero to ensure a positive value.
      input_mant[Wfl - 1] = 0;
      // Set the bit adjacent to the sign bit to 1 to ensure a normalized mantissa
      input_mant[Wfl - 2] = 1;
      // Set the remaining bits to the bit pattern stored in the last (Wfl - 2) bits in mant_i.
      input_mant.set_slc(0, mant_i.template slc<Wfl - 2>(0));
      // Use a parameterized ac_float constructor to set the mantissa and exponent of the temporary floating point input.
      T_in input_float(input_mant, sample_exponent_array[i]);
      // Make sure that input_mant was normalized and that the mantissa and exponent values haven't changed after calling the constructor.
      if (input_float.mantissa() != input_mant || input_float.exp() != sample_exponent_array[i]) {
        cout << "input_mant was not normalized correctly." << endl;
        assert(false);
      }
      if (input_float > upper_limit_rad) {
        break;
      }
      T_out output_float;
      test_ac_tan_pwl_float(input_float, output_float);
      double actual_float = output_float.to_double();
      double expected_float = T_out(tan(input_float.to_double())).to_double();
      double this_error_float;

      // If expected value of either output falls below the threshold, calculate absolute error instead of relative
      if (expected_float > threshold) { this_error_float = abs_double((expected_float - actual_float)/expected_float)* 100.0; }
      else { this_error_float = abs_double(expected_float - actual_float)* 100.0; }

      if (this_error_float > max_error_float) { max_error_float = this_error_float; }

#ifdef DEBUG
      double input_degrees = input_float.to_double() * 180 / M_PI;
      if (this_error_float > allowed_error_float) {
        cout << endl;
        cout << "  Error exceeds tolerance" << endl;
        cout << "  input_float      = " << input_float << endl;
        cout << "  input_degrees    = " << input_degrees << endl;
        cout << "  expected_float   = " << expected_float << endl;
        cout << "  actual_float     = " << actual_float << endl;
        cout << "  this_error_float = " << this_error_float << endl;
        cout << "  threshold        = " << threshold << endl;
      }
#endif
    }
  }

  bool passed = (max_error_float < allowed_error_float);

  if (passed) { printf("PASSED , max err (%f) \n", max_error_float); }
  else { printf("FAILED , max err (%f) \n", max_error_float); } // LCOV_EXCL_LINE

  if (max_error_float > cumulative_max_error_float) { cumulative_max_error_float = max_error_float; }

  return 0;
}

int main(int argc, char *argv[])
{
  double max_error_fixed = 0.0, max_error_float = 0.0;
  double allowed_error = 4.0;
  double threshold = 0.1;

  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_tan_pwl() - Allowed error " << allowed_error << endl;

  // template <int Wfi, int Ifi, int outWfi, int outIfi>
  test_driver_fixed< 8,  2, 64, 20>(max_error_fixed, allowed_error, threshold, 45);
  test_driver_fixed< 9,  1, 64, 20>(max_error_fixed, allowed_error, threshold, 50);
  test_driver_fixed<10,  3, 64, 20>(max_error_fixed, allowed_error, threshold, 55);
  test_driver_fixed< 8,  2, 64, 20>(max_error_fixed, allowed_error, threshold, 75);
  test_driver_fixed<23,  5, 64, 20>(max_error_fixed, allowed_error, threshold, 45);
  test_driver_fixed<23,  5, 64, 20>(max_error_fixed, allowed_error, threshold, 75);
  test_driver_fixed<25,  2, 64, 20>(max_error_fixed, allowed_error, threshold, 45);
  test_driver_fixed<25,  2, 64, 20>(max_error_fixed, allowed_error, threshold, 60);
  test_driver_fixed<25,  2, 64, 20>(max_error_fixed, allowed_error, threshold, 80);
  test_driver_fixed<25,  2, 64, 20>(max_error_fixed, allowed_error, threshold, 89.18);
  test_driver_fixed<22, -2, 64, 20>(max_error_fixed, allowed_error, threshold, 45);
  test_driver_fixed<22, -3, 64, 20>(max_error_fixed, allowed_error, threshold, 45);

  // template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
  test_driver_float<20, 2, 6, 20, 2, 9>(max_error_float, allowed_error, threshold, 45);
  test_driver_float<20, 2, 6, 20, 2, 9>(max_error_float, allowed_error, threshold, 70);
  test_driver_float<20, 2, 6, 20, 2, 9>(max_error_float, allowed_error, threshold, 89.18);
  test_driver_float<17, 2, 8, 25, 2, 9>(max_error_float, allowed_error, threshold, 45);
  test_driver_float<15, 2, 9, 25, 2, 9>(max_error_float, allowed_error, threshold, 70);
  test_driver_float< 9, 2, 5, 25, 2, 9>(max_error_float, allowed_error, threshold, 89.18);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum errors observed across all bit-width variations:" << endl;
  cout << "    max_error_fixed = " << max_error_fixed << endl;
  cout << "    max_error_float = " << max_error_float << endl;

  // If error limits on any test value have been crossed, the test has failed
  // Notify the user that the test was a failure if that is the case.
  if (max_error_fixed > allowed_error || max_error_float > allowed_error) {
    cout << "  ac_tan_pwl - FAILED - Error tolerance(s) exceeded" << endl; // LCOV_EXCL_LINE
    cout << "=============================================================================" << endl; // LCOV_EXCL_LINE
    return (-1); // LCOV_EXCL_LINE
  }

  cout << "  ac_tan_pwl - PASSED" << endl;
  cout << "=============================================================================" << endl;
  return (0);
}
