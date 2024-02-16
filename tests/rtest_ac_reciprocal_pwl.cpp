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
// ac_reciprocal_pwl() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_reciprocal_pwl.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_reciprocal_pwl.h>
using namespace ac_math;

// ==============================================================================
// Test Designs
//   These simple functions allow executing the ac_reciprocal_pwl() function
//   using multiple data types at the same time. Template parameters are
//   used to configure the bit-widths of the types.

// Test Design for real and complex fixed point values.
template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi, bool outSfi>
void test_ac_reciprocal_pwl_fixed(
  const ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP>   &in1,
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP>   &out1,
  const ac_complex<ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP> > &in2,
  ac_complex<ac_fixed<outWfi, outIfi, true, AC_TRN, AC_WRAP> > &out2
)
{
  out1 = ac_reciprocal_pwl<ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> >(in1);
  out2 = ac_reciprocal_pwl<ac_complex<ac_fixed<outWfi, outIfi, true, AC_TRN, AC_WRAP> > >(in2);
}

// Test Design for real fixed point values.
template <int Wfi, int Ifi, int outWfi, int outIfi>
void test_ac_reciprocal_pwl_real_fixed(
  const  ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP>   &in1,
  ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP>   &out1
)
{
  out1 = ac_reciprocal_pwl<ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> >(in1);
}

// Test Design for real and complex floating point values.
template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
void test_ac_reciprocal_pwl_float(
  const ac_float<Wfl, Ifl, Efl, AC_TRN>   &in1,
  ac_float<outWfl, outIfl, outEfl, AC_TRN>   &out1,
  const ac_complex<ac_float<Wfl, Ifl, Efl, AC_TRN> > &in2,
  ac_complex<ac_float<outWfl, outIfl, outEfl, AC_TRN> > &out2
)
{
  out1 = ac_reciprocal_pwl<ac_float<outWfl, outIfl, outEfl, AC_TRN> >(in1);
  out2 = ac_reciprocal_pwl<ac_complex<ac_float<outWfl, outIfl, outEfl, AC_TRN> > >(in2);
}

// Test Design for real and complex floating point values.
template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
void test_ac_reciprocal_pwl_real_float(
  const ac_float<Wfl, Ifl, Efl, AC_TRN>   &in1,
  ac_float<outWfl, outIfl, outEfl, AC_TRN>   &out1
)
{
  out1 = ac_reciprocal_pwl<ac_float<outWfl, outIfl, outEfl, AC_TRN> >(in1);
}

// ==============================================================================

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

// ------------------------------------------------------------------------------
// Helper functions for error calculation and monotonicity checks.

// Calculating error for real datatype.
template <class T_in, class T_out>
double err_calc(
  const T_in input,
  double &actual_value,
  const T_out output,
  const double allowed_error,
  const double threshold
)
{
  double expected_value;

  if (input.to_double() != 0) {
    // The typecasting is done in order to provide quantization on the expected output.
    expected_value = ((T_out)(1.0 / input.to_double())).to_double();
  } else {
    // If input is zero, saturate the expected output according to the type of the real/imaginary part.
    T_out output_max;
    expected_value = output_max.template set_val<AC_VAL_MAX>().to_double();
  }

  actual_value = output.to_double();
  double this_error;

  // If expected value is greater than a particular threshold, calculate relative error, else, calculate absolute error.
  if (abs(expected_value) > threshold) {
    this_error = abs( (expected_value - actual_value) / expected_value ) * 100.0;
  } else {
    this_error = abs(expected_value - actual_value) * 100.0;
  }

#ifdef DEBUG
  if (this_error > allowed_error) {
    cout << endl;
    cout << "  Error tolerance exceeded for following values: " << endl;
    cout << "  input          = " << input << endl;
    cout << "  output         = " << output << endl;
    cout << "  expected_value = " << expected_value << endl;
    cout << "  this_error     = " << this_error << endl;
    assert(false);
  }
#endif

  return this_error;
}

// Calculating error for complex datatype.
template <class T_in, class T_out>
double cmplx_err_calc(
  const ac_complex<T_in> input,
  const ac_complex<T_out> output,
  const double allowed_error,
  const double threshold
)
{
  double in_r = input.r().to_double();
  double in_i = input.i().to_double();

  // Declare variables to store the expected output value, store the difference between
  // expected and actual output values and store the actual output value (converted to
  // double)
  ac_complex<double> exp_op, diff_op, act_op;

  // Convert actual output to double and store it in a separate complex variable.
  act_op.r() = output.r().to_double();
  act_op.i() = output.i().to_double();

  // Calculate expected value of real and imaginary parts.
  exp_op.r() =  in_r / (in_r * in_r + in_i * in_i);
  exp_op.i() = -in_i / (in_r * in_r + in_i * in_i);

  if (input.r() != 0 || input.i() != 0) {
    // The typecasting is done in order to provide quantization on the expected output.
    exp_op.r() = ((T_out)exp_op.r()).to_double();
    exp_op.i() = ((T_out)exp_op.i()).to_double();
  } else {
    // If input is zero, saturate the expected output according to the type of the real/imaginary part.
    T_out output_max;
    exp_op.r() = output_max.template set_val<AC_VAL_MAX>().to_double();
    exp_op.i() = output_max.template set_val<AC_VAL_MAX>().to_double();
  }

  diff_op = exp_op - act_op;
  double this_error;

  // If magnitude of expected value is greater than a particular threshold, calculate relative error, else, calculate absolute error.
  if (sqrt(exp_op.mag_sqr()) > threshold) {this_error = sqrt((diff_op / exp_op).mag_sqr()) * 100;}
  else {this_error = sqrt(diff_op.mag_sqr()) * 100;}

#ifdef DEBUG
  if (this_error > allowed_error) {
    cout << endl;
    cout << "  Error tolerance exceeded for following values: " << endl;
    cout << "  input          = " << input << endl;
    cout << "  output         = " << output << endl;
    cout << "  exp_op         = " << exp_op << endl;
    cout << "  this_error     = " << this_error << endl;
    assert(false);
  }
#endif

  return this_error;
}

// Function for monotonicity checking in ac_fixed inputs.
template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi, bool outSfi>
void monotonicity_check(
  double &old_real_output,
  const double actual_value_fixed,
  bool &compare,
  ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP> input_fixed,
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> output_fixed
)
{
  // MONOTONIC: Make sure that function is monotonic. Compare old value (value of previous iteration) with current value. Since the reciprocal function we
  // are testing is a decreasing function, and our testbench value keeps incrementing or remains the same (in case of saturation), we expect the
  // old value to be greater than or equal to the current one.
  // Also, since the reciprocal function has a large discontinuity at x = 0; we make sure we don't compare values when we cross this point.
  // We do this by checking the signage of the old output vs that of the new output. Since we aren't checking for zero inputs, crossing x = 0
  // will mean that the old output is negative and the new one is positive, in case of an increasing testbench.
  bool sign_same = (old_real_output > 0 && actual_value_fixed > 0) || (old_real_output < 0 && actual_value_fixed < 0);
  if (compare && sign_same) {
    // Figuring out what the normalized value was for the input is a good way to figure out where the discontinuity occured w.r.t. the PWL segments.
    ac_fixed<Wfi, int(Sfi), Sfi, AC_TRN, AC_WRAP> norm_input_fixed;
    ac_normalize(input_fixed, norm_input_fixed);
    if (old_real_output < actual_value_fixed) {
      cout << endl; // LCOV_EXCL_LINE
      cout << "  Real, fixed point output not monotonic at :" << endl; // LCOV_EXCL_LINE
      cout << "  input_fixed = " << input_fixed << endl; // LCOV_EXCL_LINE
      cout << "  output_fixed = " << output_fixed << endl; // LCOV_EXCL_LINE
      cout << "  old_real_output = " << old_real_output << endl; // LCOV_EXCL_LINE
      cout << "  normalized x    = " << norm_input_fixed << endl; // LCOV_EXCL_LINE
      assert(false); // LCOV_EXCL_LINE
    }
  }
  // Update the variable for old_real_output.
  old_real_output = actual_value_fixed;
  // By setting compare to true, we make sure that once there is an old value stored, we can start comparing for monotonicity.
  compare = true;
}

// ==============================================================================
// Functions: test_driver functions
// Description: Templatized functions that can be configured for certain bit-
//   widths of AC datatypes. They use the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   piece-wise linear reciprocal model with the computed reciprocal using a
//   standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

// ==============================================================================
// Function: test_driver_fixed()
// Description: test_driver function for ac_fixed and ac_complex<ac_fixed> inputs
//   and outputs.
//   Number of iterations per run: 2^(2*Wfi)

template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi, bool outSfi>
int test_driver_fixed(
  double &cumulative_max_error_fixed,
  double &cumulative_max_error_cmplx_fixed,
  const double allowed_error_fixed,
  const double threshold,
  bool details = false
)
{
  bool passed = true;
  bool check_monotonic = true;
  double max_error_fixed = 0.0; // reset for this run
  double max_error_cmplx_fixed = 0.0; // reset for this run

  ac_fixed<   Wfi,    Ifi,    Sfi, AC_TRN, AC_WRAP>   input_fixed;
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP>   output_fixed;
  ac_complex<ac_fixed<   Wfi,    Ifi,  Sfi, AC_TRN, AC_WRAP> > cmplx_input_fixed;
  ac_complex<ac_fixed<outWfi, outIfi, true, AC_TRN, AC_WRAP> >  cmplx_output_fixed;

  double lower_limit_fixed, upper_limit_fixed, step_fixed;

  // set ranges and step size for fixed point testbench
  lower_limit_fixed   = input_fixed.template set_val<AC_VAL_MIN>().to_double();
  upper_limit_fixed   = input_fixed.template set_val<AC_VAL_MAX>().to_double();
  step_fixed          = input_fixed.template set_val<AC_VAL_QUANTUM>().to_double();

  cout << "TEST: ac_reciprocal_pwl() INPUT: ";
  cout.width(38);
  cout << left << input_fixed.type_name();
  cout << "        OUTPUTS: ";
  cout.width(38);
  cout << left << output_fixed.type_name();
  cout.width(50);
  cout << left << cmplx_output_fixed.type_name();
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << "  Ranges for input types:" << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit_fixed    = " << lower_limit_fixed << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit_fixed    = " << upper_limit_fixed << endl; // LCOV_EXCL_LINE
    cout << "    step_fixed           = " << step_fixed << endl; // LCOV_EXCL_LINE
  }

  double old_real_output;
  double actual_value_fixed;

  // test fixed-point real and complex.
  for (double i = lower_limit_fixed; i <= upper_limit_fixed; i += step_fixed) {
    bool compare = false;
    for (double j = lower_limit_fixed; j <= upper_limit_fixed; j += step_fixed) {
      cmplx_input_fixed.r() = i;
      cmplx_input_fixed.i() = j;
      input_fixed = j;
      test_ac_reciprocal_pwl_fixed(input_fixed, output_fixed, cmplx_input_fixed, cmplx_output_fixed);

      double this_error_fixed = err_calc(input_fixed, actual_value_fixed, output_fixed, allowed_error_fixed, threshold);
      double this_error_complex = cmplx_err_calc(cmplx_input_fixed, cmplx_output_fixed, allowed_error_fixed, threshold);

      if (check_monotonic) { monotonicity_check(old_real_output, actual_value_fixed, compare, input_fixed, output_fixed); }
      if (this_error_fixed > max_error_fixed) {max_error_fixed = this_error_fixed;}
      if (this_error_complex > max_error_cmplx_fixed) {max_error_cmplx_fixed = this_error_complex;}
    }
  }

  passed = (max_error_fixed < allowed_error_fixed) && (max_error_cmplx_fixed < allowed_error_fixed);

  if (passed) { printf("PASSED , max err (%f) (%f complex)\n", max_error_fixed, max_error_cmplx_fixed); }
  else        { printf("FAILED , max err (%f) (%f complex)\n", max_error_fixed, max_error_cmplx_fixed); } // LCOV_EXCL_LINE

  if (max_error_fixed>cumulative_max_error_fixed) { cumulative_max_error_fixed = max_error_fixed; }
  if (max_error_cmplx_fixed>cumulative_max_error_cmplx_fixed) { cumulative_max_error_cmplx_fixed = max_error_cmplx_fixed; }

  return 0;
}

// =================================================================================
// Function: test_driver_real_fixed()
// Description: A specialized case of the above function which only tests real
//   ac_fixed values, reducing the number of iterations and allowing the user to
//   test larger bitwidths than test_driver_fixed. 
//   Number of iterations per run = 2^Wfi.

template <int Wfi, int Ifi, int outWfi, int outIfi>
int test_driver_real_fixed(
  double &cumulative_max_error_fixed,
  const double allowed_error_fixed,
  const double threshold,
  bool details = false
)
{
  bool passed = true;
  bool check_monotonic = true;
  double max_error_fixed = 0.0; // reset for this run

  ac_fixed<   Wfi,        Ifi, false, AC_TRN, AC_WRAP>   input_fixed;
  ac_fixed<outWfi,     outIfi, false, AC_TRN, AC_WRAP>   output_fixed;

  double lower_limit_fixed, upper_limit_fixed, step_fixed;

  // set ranges and step size for fixed point testbench
  step_fixed        = input_fixed.template set_val<AC_VAL_QUANTUM>().to_double();
  lower_limit_fixed = input_fixed.template set_val<AC_VAL_MIN>().to_double();
  upper_limit_fixed = input_fixed.template set_val<AC_VAL_MAX>().to_double();

  cout << "TEST: ac_reciprocal_pwl() INPUTS: ";
  cout.width(38);
  cout << left << input_fixed.type_name();
  cout << "OUTPUTS: ";
  cout.width(38);
  cout << left << output_fixed.type_name();
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl; // LCOV_EXCL_LINE
    cout << "  Ranges for input types:" << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit_fixed = " << lower_limit_fixed << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit_fixed = " << upper_limit_fixed << endl; // LCOV_EXCL_LINE
    cout << "    step_fixed        = " << step_fixed << endl; // LCOV_EXCL_LINE
  }

  double old_real_output;
  double actual_value_fixed;

  bool compare = false;
  // Fix the real part of the input at 0, and iterate through all possible values of imaginary part based on its type.
  for (double i = lower_limit_fixed; i <= upper_limit_fixed; i += step_fixed) {
    input_fixed = i;

    // Pass all inputs at one go
    test_ac_reciprocal_pwl_real_fixed(input_fixed, output_fixed);

    double this_error_fixed = err_calc(input_fixed, actual_value_fixed, output_fixed, allowed_error_fixed, threshold);

    if (check_monotonic) { monotonicity_check(old_real_output, actual_value_fixed, compare, input_fixed, output_fixed); }
    if (this_error_fixed > max_error_fixed) {max_error_fixed = this_error_fixed;}
  }

  passed = (max_error_fixed < allowed_error_fixed);

  if (passed) { printf("PASSED , max error (%f)\n", max_error_fixed); }
  else        { printf("FAILED , max error (%f)\n", max_error_fixed); } // LCOV_EXCL_LINE

  if (max_error_fixed>cumulative_max_error_fixed) { cumulative_max_error_fixed = max_error_fixed; }

  return 0;
}

// ==============================================================================
// Function: test_driver_float()
// Description: test_driver function for ac_float and ac_complex<ac_float> inputs
//   and outputs.
//   Number of iterations per run: 2^(2*Wfl)

template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
int test_driver_float(
  double &cumulative_max_error_float,
  double &cumulative_max_error_cmplx_float,
  const double allowed_error_float,
  const double threshold,
  bool details = false
)
{
  bool passed = true;
  double max_error_float = 0.0; // reset for this run
  double max_error_cmplx_float = 0.0; // reset for this run

  typedef ac_float<   Wfl,    Ifl,    Efl, AC_TRN> T_in;
  T_in input_float;
  ac_float<outWfl, outIfl, outEfl, AC_TRN>   output_float;
  ac_complex<ac_float<   Wfl,    Ifl,    Efl, AC_TRN> > cmplx_input_float;
  ac_complex<ac_float<outWfl, outIfl, outEfl, AC_TRN> > cmplx_output_float;

  double lower_limit_mantissa, upper_limit_mantissa, step_mantissa;
  double actual_value_float;

  // Declare an ac_fixed variable of same type as mantissa
  typedef ac_fixed<Wfl, Ifl, true> T_mant;
  T_mant sample_mantissa;

  lower_limit_mantissa = sample_mantissa.template set_val<AC_VAL_MIN>().to_double();
  upper_limit_mantissa = sample_mantissa.template set_val<AC_VAL_MAX>().to_double();
  step_mantissa        = sample_mantissa.template set_val<AC_VAL_QUANTUM>().to_double();

  string empty_str = "";

  cout << "TEST: ac_reciprocal_pwl() AC_FLOAT INPUT: ";
  cout.width(38);
  cout << left << input_float.type_name();
  cout << "AC_FLOAT OUTPUT: ";
  cout.width(38);
  cout << left << output_float.type_name();
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl << "  Ranges for input types:" << endl; // LCOV_EXCL_LINE
    cout         << "    lower_limit_mantissa = " << lower_limit_mantissa << endl; // LCOV_EXCL_LINE
    cout         << "    upper_limit_mantissa = " << upper_limit_mantissa << endl; // LCOV_EXCL_LINE
    cout         << "    step_mantissa        = " << step_mantissa << endl; // LCOV_EXCL_LINE
  }

  // Declare arrays to store all values of exponent to be tested.
  const int exp_arr_size = 2*(Efl - 1) + 3;
  ac_int<Efl, true> sample_exponent;
  ac_int<Efl, true> sample_exponent_array[exp_arr_size];

  // The first element of the array is the minimum exponent value, the middle element is a zero exponent, and
  // the last element is the maximum possible value.
  sample_exponent_array[0].template set_val<AC_VAL_MIN>();
  sample_exponent_array[Efl] = 0;
  sample_exponent_array[exp_arr_size - 1].template set_val<AC_VAL_MAX>();

  // All the other elements are set to values that correspond to a one-hot encoding scheme, in which only one
  // bit of the absolute value of the exponent is set to one. Both negative and positive values are encoded this way.
  for (int i = (Efl - 2); i >= 0; i--) {
    sample_exponent = 0;
    sample_exponent[i] = 1;
    sample_exponent_array[Efl + i + 1] = sample_exponent;
    sample_exponent_array[Efl - i - 1] = -sample_exponent;
  }

  for (int i = 0; i < exp_arr_size; i++) {
    for (double j = lower_limit_mantissa; j <= upper_limit_mantissa; j += step_mantissa) {
      for (double k = lower_limit_mantissa; k <= upper_limit_mantissa; k += step_mantissa) {
        // Normalize real and imaginary mantissas before passing it to testing function.
        T_in input_float_real((T_mant)j, sample_exponent_array[i]);
        T_in input_float_imag((T_mant)k, sample_exponent_array[i]);
        input_float = input_float_real;
        cmplx_input_float.r() = input_float_real;
        cmplx_input_float.i() = input_float_imag;
        test_ac_reciprocal_pwl_float(input_float, output_float, cmplx_input_float, cmplx_output_float);

        double this_error_float = err_calc(input_float, actual_value_float, output_float, allowed_error_float, threshold);
        double this_error_complex = cmplx_err_calc(cmplx_input_float, cmplx_output_float, allowed_error_float, threshold);

        if (this_error_float > max_error_float) {max_error_float = this_error_float;}
        if (this_error_complex > max_error_cmplx_float) {max_error_cmplx_float = this_error_complex;}
      }
    }
  }

  passed = (max_error_float < allowed_error_float) && (max_error_cmplx_float < allowed_error_float);

  if (passed) { printf("PASSED , max err (%f) (%f complex)\n", max_error_float, max_error_cmplx_float); }
  else        { printf("FAILED , max err (%f) (%f complex)\n", max_error_float, max_error_cmplx_float); } // LCOV_EXCL_LINE

  if (max_error_float>cumulative_max_error_float) { cumulative_max_error_float = max_error_float; }
  if (max_error_cmplx_float>cumulative_max_error_cmplx_float) { cumulative_max_error_cmplx_float = max_error_cmplx_float; }

  return 0;
}

// =================================================================================
// Function: test_driver_real_float()
// Description: A specialized case of the above function which only tests real
//   ac_float values, reducing the number of iterations allowing the user to test 
//   larger bitwidths than test_driver_float().
//   Number of iterations per run = 2^Wfl.

template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
int test_driver_real_float(
  double &cumulative_max_error_float,
  const double allowed_error_float,
  const double threshold,
  bool details = false
)
{
  bool passed = true;
  double max_error_float = 0.0; // reset for this run

  typedef ac_float<Wfl, Ifl, Efl, AC_TRN> T_in;
  typedef ac_float<outWfl, outIfl, outEfl, AC_TRN> T_out;

  // Since ac_float values are normalized, the bit adjacent to the sign bit in the mantissa
  // will always be set to 1 except for certain negative values. We will hence cycle through all the
  // bit patterns that correspond to the last (Wfl - 2) bits in the mantissa.
  ac_int<Wfl - 2, false> sample_mantissa_slc;
  // Set the lower limit, upper limit and step size of the test iterations.
  ac_int<Wfl - 2, false> lower_limit_it = 0;
  ac_int<Wfl - 2, false> upper_limit_it = sample_mantissa_slc.template set_val<AC_VAL_MAX>().to_double();
  ac_int<Wfl - 2, false> step_it = 1; // Since sample_mantissa_slc is an integer.

  // Dump the test details
  if (details) {
    cout << endl; // LCOV_EXCL_LINE
    cout << "  Ranges for testing iterations:" << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit_it = " << lower_limit_it << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit_it = " << upper_limit_it << endl; // LCOV_EXCL_LINE
    cout << "    step_it        = " << step_it << endl; // LCOV_EXCL_LINE
  }

  // Declare arrays to store all values of exponent to be tested.
  const int exp_arr_size = 2*(Efl - 1) + 3;
  ac_int<Efl, true> sample_exponent;
  ac_int<Efl, true> sample_exponent_array[exp_arr_size];

  // The first element of the array is the minimum exponent value, the middle element is a zero exponent, and
  // the last element is the maximum possible value.
  sample_exponent_array[0].template set_val<AC_VAL_MIN>();
  sample_exponent_array[Efl] = 0;
  sample_exponent_array[exp_arr_size - 1].template set_val<AC_VAL_MAX>();

  // All the other elements are set to values that correspond to a one-hot encoding scheme, in which only one
  // bit of the absolute value of the exponent is set to one. Both negative and positive values are encoded this way.
  for (int i = (Efl - 2); i >= 0; i--) {
    sample_exponent = 0;
    sample_exponent[i] = 1;
    sample_exponent_array[Efl + i + 1] = sample_exponent;
    sample_exponent_array[Efl - i - 1] = -sample_exponent;
  }

  string empty_str = "";

  cout << "TEST: ac_reciprocal_pwl() AC_FLOAT INPUT: ";
  cout.width(38);
  cout << left << T_in::type_name();
  cout << "AC_FLOAT OUTPUT: ";
  cout.width(38);
  cout << left << T_out::type_name();
  cout << "RESULT: ";

  for (int i = 0; i < exp_arr_size; i++) {
    // For a particular exponent value, go through every possible value that can be represented by the mantissa.
    // The iteration variable has a bitwidth that is 1 higher (bitwidth = Wfl - 1) than the slice of the mantissa
    // we'll be changing from one iteration to the other (bitwidth of slice = Wfl - 2), to ensure that the loop variable
    // doesn't overflow.
    for (ac_int<Wfl - 1, false> mant_i = lower_limit_it; mant_i <= upper_limit_it; mant_i += step_it) {
      ac_fixed<Wfl, Ifl, true> input_mant;
      // Set the sign bit to zero to ensure a positive value.
      input_mant[Wfl - 1] = 0;
      // Set the bit adjacent to the sign bit to 1 to ensure a normalized mantissa
      input_mant[Wfl - 2] = 1;
      // Set the remaining bits to the bit pattern stored in the last (Wfl - 2) bits in mant_i.
      input_mant.set_slc(0, mant_i.template slc<Wfl - 2>(0));
      // Use a parameterized ac_float constructor to set the mantissa and exponent of the floating point input.
      T_in input_float(input_mant, sample_exponent_array[i]);
      // Make sure that input_mant was normalized and that the mantissa and exponent values haven't changed after calling the constructor.
      if (input_float.mantissa() != input_mant || input_float.exp() != sample_exponent_array[i]) {
        cout << "input_mant was not normalized correctly." << endl;
        assert(false);
      }
      T_out output_float;
      test_ac_reciprocal_pwl_real_float(input_float, output_float);
      double actual_value_float;
      double this_error_float = err_calc(input_float, actual_value_float, output_float, allowed_error_float, threshold);
      if (this_error_float > max_error_float) { max_error_float = this_error_float; }
      // Pass the negative value of the input and test the output too, in a similar manner.
      ac_fixed<Wfl, Ifl, true> input_mant_neg = -input_mant;
      T_in input_float_neg(input_mant_neg, sample_exponent_array[i]);
      if ((input_float_neg.mantissa() != input_mant_neg || input_float_neg.exp() != sample_exponent_array[i]) && mant_i != 0) {
        cout << "input_mant_neg was not normalized correctly." << endl;
        assert(false);
      }
      test_ac_reciprocal_pwl_real_float(input_float_neg, output_float);
      double actual_value_float_neg;
      double this_error_float_neg = err_calc(input_float_neg, actual_value_float_neg, output_float, allowed_error_float, threshold);
      if (this_error_float_neg > max_error_float) { max_error_float = this_error_float_neg; }
    }
  }

  passed = (max_error_float < allowed_error_float);

  if (passed) { printf("PASSED , max err (%f)\n", max_error_float); }
  else        { printf("FAILED , max err (%f)\n", max_error_float); } // LCOV_EXCL_LINE

  if (max_error_float > cumulative_max_error_float) { cumulative_max_error_float = max_error_float; }

  return 0;
}


int main(int argc, char *argv[])
{
  double max_error_fixed = 0, cmplx_max_error_fixed = 0, max_error_float = 0, cmplx_max_error_float = 0;
  double allowed_error_fixed = 0.5;
  double allowed_error_float = 0.5;
  const double threshold = 0.005;
  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_reciprocal_pwl() - Allowed error " << allowed_error_fixed << " (fixed pt), " << allowed_error_float << " (float pt)" << endl;

  // template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi, bool outSfi>
  test_driver_fixed< 10,  3,  true, 64, 32,  true>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, threshold);
  test_driver_fixed< 10,  1,  true, 64, 32,  true>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, threshold);
  test_driver_fixed< 10,  0, false, 64, 32, false>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, threshold);
  test_driver_fixed< 10,  2, false, 64, 32,  true>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, threshold);
  test_driver_fixed<  4,  9,  true, 64, 32,  true>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, threshold);
  test_driver_fixed<  4, -2,  true, 64, 32,  true>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, threshold);
  test_driver_fixed<  5,  8, false, 64, 32, false>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, threshold);
  test_driver_fixed<  4, -2, false, 60, 30, false>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, threshold);
  test_driver_fixed< 10,  4,  true, 64, 32,  true>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, threshold);
  test_driver_fixed< 10,  3, false, 64, 32, false>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, threshold);
  test_driver_fixed<  9,  4,  true, 60, 30,  true>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, threshold);
  test_driver_fixed<  9,  2, false, 64, 32, false>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, threshold);

  // template <int Wfi, int Ifi, int outWfi, int outIfi>
  test_driver_real_fixed< 21,  8, 64, 32>(max_error_fixed, allowed_error_fixed, threshold);
  test_driver_real_fixed< 21,  7, 64, 32>(max_error_fixed, allowed_error_fixed, threshold);
  test_driver_real_fixed< 21,  2, 60, 30>(max_error_fixed, allowed_error_fixed, threshold);
  test_driver_real_fixed< 21,  1, 64, 32>(max_error_fixed, allowed_error_fixed, threshold);
  test_driver_real_fixed< 21,  0, 64, 32>(max_error_fixed, allowed_error_fixed, threshold);
  test_driver_real_fixed< 21, -1, 60, 30>(max_error_fixed, allowed_error_fixed, threshold);
  test_driver_real_fixed< 21, -2, 63, 33>(max_error_fixed, allowed_error_fixed, threshold);
  test_driver_real_fixed< 11, 20, 64, 32>(max_error_fixed, allowed_error_fixed, threshold);
  test_driver_real_fixed< 11, 21, 64, 32>(max_error_fixed, allowed_error_fixed, threshold);

  // template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
  test_driver_float<5,  3, 3, 64, 32, 10>(max_error_float, cmplx_max_error_float, allowed_error_float, threshold);
  test_driver_float<5,  1, 8, 64, 32, 10>(max_error_float, cmplx_max_error_float, allowed_error_float, threshold);
  test_driver_float<5,  0, 3, 64, 32, 10>(max_error_float, cmplx_max_error_float, allowed_error_float, threshold);
  test_driver_float<5, -2, 3, 60, 30, 11>(max_error_float, cmplx_max_error_float, allowed_error_float, threshold);
  test_driver_float<5,  9, 6, 64, 32, 10>(max_error_float, cmplx_max_error_float, allowed_error_float, threshold);
  test_driver_float<5,  5, 3, 64, 32, 10>(max_error_float, cmplx_max_error_float, allowed_error_float, threshold);
  test_driver_float<5,  5, 1, 60, 30, 11>(max_error_float, cmplx_max_error_float, allowed_error_float, threshold);
  test_driver_float<9,  5, 4, 64, 32, 10>(max_error_float, cmplx_max_error_float, allowed_error_float, threshold);
  test_driver_float<5,  3, 5, 64, 32, 10>(max_error_float, cmplx_max_error_float, allowed_error_float, threshold);

  // template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
  test_driver_real_float< 21,  9, 8, 21,  9,  8>(max_error_float, allowed_error_float, threshold);
  test_driver_real_float< 21,  8, 8, 21,  8,  8>(max_error_float, allowed_error_float, threshold);
  test_driver_real_float< 21,  0, 8, 21,  0,  8>(max_error_float, allowed_error_float, threshold);
  test_driver_real_float< 21,  0, 8, 21, -8,  8>(max_error_float, allowed_error_float, threshold);
  test_driver_real_float< 21,  1, 8, 21,  1,  8>(max_error_float, allowed_error_float, threshold);
  test_driver_real_float< 21,  2, 8, 21,  8,  8>(max_error_float, allowed_error_float, threshold);
  test_driver_real_float< 21, -1, 8, 21, -1,  8>(max_error_float, allowed_error_float, threshold);
  test_driver_real_float< 21, -2, 8, 21,  8,  8>(max_error_float, allowed_error_float, threshold);
  test_driver_real_float< 21, 11, 8, 21, -1,  8>(max_error_float, allowed_error_float, threshold);
  test_driver_real_float< 21, 12, 8, 21,  8,  8>(max_error_float, allowed_error_float, threshold);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum errors observed across all data type / bit-width variations:" << endl;
  cout << "    max_error_fixed       = " << max_error_fixed << endl;
  cout << "    cmplx_max_error_fixed = " << cmplx_max_error_fixed << endl;
  cout << "    max_error_float       = " << max_error_float << endl;
  cout << "    cmplx_max_error_float = " << cmplx_max_error_float << endl;

  // If error limits on any tested datatype have been crossed, the test has failed
  bool test_fail = (max_error_fixed > allowed_error_fixed) || (cmplx_max_error_fixed > allowed_error_fixed) || (max_error_float > allowed_error_float) || (cmplx_max_error_float > allowed_error_float);

  // Notify the user whether or not the test was a failure.
  if (test_fail) {
    cout << "  ac_reciprocal_pwl - FAILED - Error tolerance(s) exceeded" << endl; // LCOV_EXCL_LINE
    cout << "=============================================================================" << endl; // LCOV_EXCL_LINE
    return -1; // LCOV_EXCL_LINE
  } else {
    cout << "  ac_reciprocal_pwl - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }
  return 0;
}


