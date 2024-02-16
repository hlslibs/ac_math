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
// ac_atan_pwl_vha() function using a variety of data types and bit-
// widths.

// To compile satandalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_atan_pwl_vha.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_atan_pwl_vha.h>
using namespace ac_math;

// ==============================================================================
// Test Designs
//   This simple function allow executing the ac_atan_pwl_vha() function.
//   Template parameters are used to configure the bit-widths of the
//   inputs and outputs.

// Test design for real fixed point values.
template <int Wfi, int Ifi, int outWfi, int outIfi>
void test_ac_atan_pwl_vha(
  const ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP> &in,
  ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> &atan_out
)
{
  atan_out = ac_atan_pwl_vha<ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> >(in);
}

// Test Design for real floating point values.
template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
void test_ac_atan_pwl_vha_float(
  const ac_float<Wfl, Ifl, Efl, AC_TRN>    &in1,
  ac_float<outWfl, outIfl, outEfl, AC_TRN> &out1
)
{
  out1 = ac_atan_pwl_vha<ac_float<outWfl, outIfl, outEfl, AC_TRN> >(in1);
}

// ------------------------------------------------------------------------------
// Helper function for absolute value calculation. This can avoid any naming conflicts
// with other absolute value functions.

double abs_double(double x)
{
  return x >= 0 ? x : -x;
}

// ==============================================================================

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

// ===============================================================================
// Function: test_driver_fixed()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_fixed inputs. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   piecewise linear atan model with the computed arctangent using a
//   standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <int Wfi, int Ifi, int outWfi, int outIfi>
int test_driver_fixed(
  double &cumulative_max_error_fixed,
  const double allowed_error,
  bool details = false
)
{
  double i; // make loop variable slightly larger
  ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP> input;
  typedef ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> T_out;
  T_out atan_out;

  // set ranges and step size for testbench
  double lower_limit = input.template set_val<AC_VAL_MIN>().to_double();
  double upper_limit = input.template set_val<AC_VAL_MAX>().to_double();
  double step        = input.template set_val<AC_VAL_QUANTUM>().to_double();

  cout << "TEST: ac_atan_pwl_vha() INPUT: ";
  cout.width(38);
  cout << left << input.type_name();
  cout << "OUTPUT: ";
  cout.width(38);
  cout << left << atan_out.type_name();
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
  double max_atan_error = 0.0;

  bool check_monotonic = true;
  bool compare_atan = false;
  double old_output_atan;

  for (i = lower_limit; i <= upper_limit; i += step) {
    // Set values for input.
    input = i;

    // call reference atan() with fixed-pt value converted back to double
    // an additional step of typecasting is required in order to perform
    // quantization on the expected output.
    double expected_atan = ((T_out)atan(input.to_double())).to_double();

    // call DUT with fixed-pt value
    test_ac_atan_pwl_vha(input, atan_out);

    double actual_atan = atan_out.to_double();
    double this_error_fixed;

    // Calculate absolute error.
    this_error_fixed = abs_double(expected_atan - actual_atan) * 100.0;

    if (check_monotonic) {
      // MONOTONIC: Make sure that function is monotonic. Compare old value (value of previous iteration) with current value. Since the arctangent function we
      // are testing is an increasing function, and our testbench value keeps incrementing or remains the same (in case of saturation), we expect the
      // old value to be lesser than or equal to the current one.

      // This comparison is only carried out once there is an old value to compare with
      if (compare_atan) {
        // if by any chance the function output has dropped in value, print out at what point the problem has occured and throw a runtime assertion.
        if (old_output_atan > actual_atan) {
          cout << endl; // LCOV_EXCL_LINE
          cout << "  atan output not monotonic at :" << endl; // LCOV_EXCL_LINE
          cout << "  x = " << input << endl; // LCOV_EXCL_LINE
          cout << "  y = " << atan_out << endl; // LCOV_EXCL_LINE
          cout << "  old_output_atan = " << old_output_atan << endl; // LCOV_EXCL_LINE
          assert(false); // LCOV_EXCL_LINE
        }
      }
      // Update the old value
      old_output_atan = actual_atan;
      // Once an old value has been stored, i.e. towards the end of the first iteration, this value is set to true.
      compare_atan = true;
    }

#ifdef DEBUG
    double output_degrees = actual_atan * 180 / M_PI;
    if (this_error_fixed > allowed_error) {
      cout << endl;
      cout << "  Error exceeds tolerance" << endl;
      cout << "  input           = " << input << endl;
      cout << "  expected_atan   = " << expected_atan << endl;
      cout << "  actual_atan     = " << actual_atan << endl;
      cout << "  output_degrees  = " << output_degrees << endl;
      cout << "  this_error_fixed = " << this_error_fixed << endl;
      assert(false);
    }
#endif

    if (this_error_fixed > max_atan_error) { max_atan_error = this_error_fixed; }
  }

  if (max_atan_error > cumulative_max_error_fixed) { cumulative_max_error_fixed = max_atan_error; }

  passed = (max_atan_error < allowed_error);

  if (passed) { printf("PASSED , max err (%f)\n", max_atan_error); }
  else        { printf("FAILED , max err (%f)\n", max_atan_error); } // LCOV_EXCL_LINE

  return 0;
}

// =================================================================================
// Function: test_driver_float()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_float inputs. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   piecewise linear atan model with the computed arctangent using a
//   standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
int test_driver_float(
  double &cumulative_max_error_float,
  const double allowed_error_float,
  bool details = false
)
{
  bool passed = true;
  double max_error_float = 0.0; // reset for this run

  typedef ac_float<   Wfl,    Ifl,    Efl> T_in;
  typedef ac_float<outWfl, outIfl, outEfl> T_out;

  // Since ac_float values are normalized, the bit adjacent to the sign bit in the mantissa
  // will always be set to 1. We will hence cycle through all the bit patterns that correspond to the last (Wfl - 2)
  // bits in the mantissa.
  ac_int<Wfl - 2, false> sample_mantissa_slc;
  // Set the lower limit, upper limit and step size of the test iterations.
  ac_int<Wfl - 2, false> lower_limit_it = 0;
  ac_int<Wfl - 2, false> upper_limit_it = sample_mantissa_slc.template set_val<AC_VAL_MAX>().to_double();
  ac_int<Wfl - 2, false> step_it = 1; // Since sample_mantissa_slc is an integer.

  cout << "TEST: ac_atan_pwl_vha() INPUT: ";
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
    cout << endl << "  Ranges for testing iterations:" << endl; // LCOV_EXCL_LINE
    cout         << "    lower_limit_it       = " << lower_limit_it << endl; // LCOV_EXCL_LINE
    cout         << "    upper_limit_it       = " << upper_limit_it << endl; // LCOV_EXCL_LINE
    cout         << "    step_it              = " << step_it << endl; // LCOV_EXCL_LINE
    cout         << "    allowed_error_float  = " << allowed_error_float << endl; // LCOV_EXCL_LINE
  }

  for (int i = 0; i < exp_arr_size; i++) {
    // For a particular exponent value, go through every possible value that can be represented by the mantissa.
    // The iteration variable has a bitwidth that is 1 higher (bitwidth = Wfl - 1) than the slice of the mantissa
    // we'll be changing from one iteration to the other (bitwidth of slice = Wfl - 2), to ensure that the iteration
    // variable does not overflow.
    for (ac_int<Wfl - 1, false> mant_i = lower_limit_it; mant_i <= upper_limit_it; mant_i += step_it) {
      ac_fixed<Wfl, Ifl, true> input_mant;
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
      T_out output_float;
      test_ac_atan_pwl_vha_float(input_float, output_float);
      double actual_float = output_float.to_double();
      double expected_float = T_out(atan(input_float.to_double())).to_double();
      double this_error_float = abs_double(expected_float - actual_float)*100;
      if (this_error_float > max_error_float) { max_error_float = this_error_float; }

#ifdef DEBUG
      double output_degrees = actual_float * 180 / M_PI;
      if (this_error_float > allowed_error_float) {
        cout << endl;
        cout << "  Error exceeds tolerance" << endl;
        cout << "  input_float      = " << input_float << endl;
        cout << "  expected_float   = " << expected_float << endl;
        cout << "  actual_float     = " << output_float << endl;
        cout << "  output_degrees   = " << output_degrees << endl;
        cout << "  this_error_float = " << this_error_float << endl;
        assert(false);
      }
#endif
    }
  }

  passed = (max_error_float < allowed_error_float);

  if (passed) { printf("PASSED , max err (%f) \n", max_error_float); }
  else        { printf("FAILED , max err (%f) \n", max_error_float); } // LCOV_EXCL_LINE

  if (max_error_float > cumulative_max_error_float) { cumulative_max_error_float = max_error_float; }

  return 0;
}

int main(int argc, char *argv[])
{
  double max_error_fixed = 0.0;
  double allowed_error = 0.005;

  double max_error_float = 0.0;

  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_atan_pwl_vha() - Allowed error " << allowed_error << endl;

  // template <int Wfi, int Ifi, int outWfi, int outIfi>
  test_driver_fixed< 5,  2, 33,  1>(max_error_fixed, allowed_error);
  test_driver_fixed< 5, -2, 33,  1>(max_error_fixed, allowed_error);
  test_driver_fixed< 2, -2, 33,  1>(max_error_fixed, allowed_error);
  test_driver_fixed< 7,  3, 33,  1>(max_error_fixed, allowed_error);
  test_driver_fixed<10,  5, 33,  1>(max_error_fixed, allowed_error);
  test_driver_fixed<20,  6, 33,  1>(max_error_fixed, allowed_error);
  test_driver_fixed<20, -6, 33,  1>(max_error_fixed, allowed_error);
  test_driver_fixed<22,  8, 33,  1>(max_error_fixed, allowed_error);
  test_driver_fixed<20,  0, 33,  1>(max_error_fixed, allowed_error);
  test_driver_fixed< 8, 12, 33,  1>(max_error_fixed, allowed_error);
  test_driver_fixed<20,  6, 42, 10>(max_error_fixed, allowed_error);
  test_driver_fixed<20, -6, 42, 10>(max_error_fixed, allowed_error);
  test_driver_fixed<22,  8, 42, 10>(max_error_fixed, allowed_error);
  test_driver_fixed<20,  0, 42, 10>(max_error_fixed, allowed_error);
  test_driver_fixed< 8, 12, 42, 10>(max_error_fixed, allowed_error);

  test_driver_float<20, 12, 10, 32,  1, 15>(max_error_float, allowed_error);
  test_driver_float<20,  4, 10, 32,  1, 15>(max_error_float, allowed_error);
  test_driver_float<20,  2, 10, 32,  1, 15>(max_error_float, allowed_error);
  test_driver_float<20,  1, 10, 32,  1, 15>(max_error_float, allowed_error);
  test_driver_float<20,  0, 10, 32,  1, 15>(max_error_float, allowed_error);
  test_driver_float<20, -1, 10, 32,  1, 15>(max_error_float, allowed_error);
  test_driver_float<20, -2, 10, 32,  1, 15>(max_error_float, allowed_error);
  test_driver_float<20, 12, 10, 32,  6, 15>(max_error_float, allowed_error);
  test_driver_float<20,  4, 10, 32,  6, 15>(max_error_float, allowed_error);
  test_driver_float<20,  2, 10, 32,  6, 15>(max_error_float, allowed_error);
  test_driver_float<20,  1, 10, 32,  6, 15>(max_error_float, allowed_error);
  test_driver_float<20,  0, 10, 32,  6, 15>(max_error_float, allowed_error);
  test_driver_float<20, -1, 10, 32,  6, 15>(max_error_float, allowed_error);
  test_driver_float<20, -2, 10, 32,  6, 15>(max_error_float, allowed_error);
  test_driver_float<20, 12, 10, 32, -4, 15>(max_error_float, allowed_error);
  test_driver_float<20,  4, 10, 32, -4, 15>(max_error_float, allowed_error);
  test_driver_float<20,  2, 10, 32, -4, 15>(max_error_float, allowed_error);
  test_driver_float<20,  1, 10, 32, -4, 15>(max_error_float, allowed_error);
  test_driver_float<20,  0, 10, 32, -4, 15>(max_error_float, allowed_error);
  test_driver_float<20, -1, 10, 32, -4, 15>(max_error_float, allowed_error);
  test_driver_float<20, -2, 10, 32, -4, 15>(max_error_float, allowed_error);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum errors observed across all data type / bit-width variations:" << endl;
  cout << "    max_error_fixed = " << max_error_fixed << endl;
  cout << "    max_error_float = " << max_error_float << endl;

  bool test_fail = (max_error_fixed > allowed_error) || (max_error_float > allowed_error);

  // If error limits on any test value have been crossed, the test has failed
  // Notify the user that the test was a failure if that is the case.
  if (test_fail) {
    cout << "  ac_atan_pwl_vha - FAILED - Error tolerance(s) exceeded" << endl; // LCOV_EXCL_LINE
    cout << "=============================================================================" << endl; // LCOV_EXCL_LINE
    return (-1); // LCOV_EXCL_LINE
  }

  cout << "  ac_atan_pwl_vha - PASSED" << endl;
  cout << "=============================================================================" << endl;
  return (0);
}
