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
// ac_abs() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_abs.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench

#include <ac_math/ac_abs.h>
using namespace ac_math;

// ==============================================================================
// Test Designs
//   These simple functions allow executing the ac_abs() function
//   using multiple data types at the same time. Template parameters are used
//   to configure the bit-widths of the types.

// Test for signed and unsigned ac_int outputs.
template <int Wint, int outWint>
void test_ac_abs_int(
  const ac_int<Wint, true> &in1,
  ac_int<outWint, false> &out1,
  ac_int<outWint + 1, true> &out2
)
{
  ac_abs(in1, out1);
  ac_abs(in1, out2);
}

//Test for signed and unsigned ac_fixed outputs.
template <int Wfi, int Ifi, int outWfi, int outIfi>
void test_ac_abs_fixed(
  const ac_fixed<Wfi, Ifi, true, AC_TRN, AC_WRAP> &in2,
  ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> &out3,
  ac_fixed<outWfi + 1, outIfi + 1, true, AC_TRN, AC_WRAP> &out4
)
{
  ac_abs(in2, out3);
  ac_abs(in2, out4);
}

// Test for ac_float input and output.
template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
void test_ac_abs_float(
  const ac_float<Wfl, Ifl, Efl, AC_TRN> &in3,
  ac_float<outWfl, outIfl, outEfl, AC_TRN> &out5
)
{
  ac_abs(in3, out5);
}

//Test for ac_complex input and ac_fixed output.
template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi>
void test_ac_abs_complex(
  const ac_complex<ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP> > &in,
  ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> &out
)
{
  ac_abs(in, out);
}

// ==============================================================================

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

// ------------------------------------------------------------------------------
// Helper function to check if output is correct

template <class T_in, class T_out>
bool output_check(
  const T_in input,
  const T_out output
)
{
  bool correct = abs(input.to_double()) == output.to_double();

#ifdef DEBUG
  if (!correct) {
    cout << endl;
    cout << "  Output not correct" << endl;
    cout << "  input  = " << input << endl;
    cout << "  output = " << output << endl;
    assert(false);
  }
#endif

  return correct;
}

// Calculating error for complex datatype.
template <class T_in, class T_out>
double cmplx_err_calc(
  const ac_complex<T_in> input_complex,
  const T_out output_complex,
  const double allowed_error_cmplx,
  const double threshold
)
{
  double this_error;
  double a = pow((input_complex.r()).to_double(), 2) + pow((input_complex.i()).to_double(), 2);
  double actual_value = sqrt(a);
  double expected_value = output_complex.to_double();

  // If magnitude of expected value is greater than a particular threshold, calculate relative error, else, calculate absolute error.
  if (abs(expected_value) > threshold) {
    this_error = ((abs(expected_value - actual_value)) / expected_value ) * 100;
  } else {
    this_error = abs(expected_value - actual_value);
  }

#ifdef DEBUG
  if (this_error > allowed_error_cmplx) {
    cout << endl;
    cout << "  Error tolerance exceeded for following values: " << endl;
    cout << "  input          = " << input_complex << endl;
    cout << "  output         = " << output_complex << endl;
    cout << "  exp_op         = " << expected_value << endl;
    cout << "  this_error     = " << this_error << endl;
    assert(false);
  }
#endif

  return this_error;
}
// ==============================================================================
// Functions: test_driver functions
// Description: Templatized functions that can be configured for certain bit-
//   widths of AC datatypes. They use the type information to iterate through a
//   range of valid values on that type and make sure that the input and
//   output of the ac_abs function match each other in terms of absolute value.

// ==============================================================================
// Function: test_driver_int()
// Description: test_driver function for ac_int inputs and outputs.

template <int Wint, int outWint>
int test_driver_int(
  bool &all_tests_pass,
  bool details = false
)
{
  ac_int<       Wint,  true> input_s_int;
  ac_int<    outWint, false> output_us_int;
  ac_int<outWint + 1,  true> output_s_int;

  double lower_limit, upper_limit, step;

  // set ranges and step size for integer testbench
  lower_limit = input_s_int.template set_val<AC_VAL_MIN>().to_double();
  upper_limit = input_s_int.template set_val<AC_VAL_MAX>().to_double();
  step        = input_s_int.template set_val<AC_VAL_QUANTUM>().to_double();

  cout << "TEST: ac_abs() INPUT: ";
  cout.width(50);
  cout << left << input_s_int.type_name();
  cout << "OUTPUTS: ";
  cout.width(38);
  cout << left << output_us_int.type_name();
  cout.width(38);
  cout << left << output_s_int.type_name();
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl; // LCOV_EXCL_LINE
    cout << "  Ranges for input types:" << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit = " << lower_limit << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit = " << upper_limit << endl; // LCOV_EXCL_LINE
    cout << "    step        = " << step << endl; // LCOV_EXCL_LINE
  }

  bool correct = true;

  // test integer values.

  for (double i = lower_limit; i <= upper_limit; i += step) {
    input_s_int = i;
    test_ac_abs_int(input_s_int, output_us_int, output_s_int);
    bool correct_iteration = output_check(input_s_int, output_us_int) && output_check(input_s_int, output_s_int);
    // If any iteration does not produce the correct value for the output, then the "correct" variable will be set to false.
    correct = correct && correct_iteration;
  }

  if (correct) { printf("PASSED\n"); }
  else         { printf("FAILED\n"); } // LCOV_EXCL_LINE

  all_tests_pass = all_tests_pass && correct;

  return 0;
}

// ==============================================================================
// Function: test_driver_fixed()
// Description: test_driver function for ac_fixed inputs and outputs.

template <int Wfi, int Ifi, int outWfi, int outIfi>
int test_driver_fixed(
  bool &all_tests_pass,
  bool details = false
)
{
  ac_fixed<       Wfi,        Ifi,  true, AC_TRN, AC_WRAP> input_s_fixed;
  ac_fixed<    outWfi,     outIfi, false, AC_TRN, AC_WRAP> output_us_fixed;
  ac_fixed<outWfi + 1, outIfi + 1,  true, AC_TRN, AC_WRAP> output_s_fixed;

  double lower_limit, upper_limit, step;

  // set ranges and step size for fixed point testbench
  lower_limit = input_s_fixed.template set_val<AC_VAL_MIN>().to_double();
  upper_limit = input_s_fixed.template set_val<AC_VAL_MAX>().to_double();
  step        = input_s_fixed.template set_val<AC_VAL_QUANTUM>().to_double();

  cout << "TEST: ac_abs() INPUT: ";
  cout.width(50);
  cout << left << input_s_fixed.type_name();
  cout << "OUTPUTS: ";
  cout.width(38);
  cout << left << output_us_fixed.type_name();
  cout.width(38);
  cout << left << output_s_fixed.type_name();
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl; // LCOV_EXCL_LINE
    cout << "  Ranges for input types:" << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit = " << lower_limit << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit = " << upper_limit << endl; // LCOV_EXCL_LINE
    cout << "    step        = " << step << endl; // LCOV_EXCL_LINE
  }

  bool correct = true;

  // test fixed-point values.

  for (double i = lower_limit; i <= upper_limit; i += step) {
    input_s_fixed = i;
    test_ac_abs_fixed(input_s_fixed, output_us_fixed, output_s_fixed);
    bool correct_iteration = output_check(input_s_fixed, output_us_fixed) && output_check(input_s_fixed, output_s_fixed);
    // If any iteration does not produce the correct value for the output, then the "correct" variable will be set to false.
    correct = correct && correct_iteration;
  }

  if (correct) { printf("PASSED\n"); }
  else         { printf("FAILED\n"); } // LCOV_EXCL_LINE

  all_tests_pass = all_tests_pass && correct;

  return 0;
}

// ==============================================================================
// Function: test_driver_float()
// Description: test_driver function for ac_float inputs and outputs.

template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
int test_driver_float(
  bool &all_tests_pass,
  bool details = false
)
{
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
  // bit of the absolute value of the exponent is set to one. Both negative and positive values are encoded this way,
  // and all the bits are covered
  for (int i = (Efl - 2); i >= 0; i--) {
    sample_exponent = 0;
    sample_exponent[i] = 1;
    sample_exponent_array[Efl + i + 1] = sample_exponent;
    sample_exponent_array[Efl - i - 1] = -sample_exponent;
  }

  string empty_str = "";

  cout << "TEST: ac_abs() INPUT: ";
  cout.width(50);
  cout << left << T_in::type_name();
  cout << "OUTPUT:  ";
  cout.width(38);
  cout << left << T_out::type_name();
  cout.width(38);
  cout << left << empty_str;
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl; // LCOV_EXCL_LINE
    cout << "  Ranges for testing iterations:" << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit_it = " << lower_limit_it << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit_it = " << upper_limit_it << endl; // LCOV_EXCL_LINE
    cout << "    step_it        = " << step_it << endl; // LCOV_EXCL_LINE
  }

  bool correct = true;

  // test floating-point values.

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
      test_ac_abs_float(input_float, output_float);
      // If any iteration does not produce the correct value for the output, then the "correct" variable will be set to false.
      bool correct_iteration = output_check(input_float, output_float);
      correct = correct && correct_iteration;
      // Pass the negative value of the input and test the output too, in a similar manner.
      ac_fixed<Wfl, Ifl, true> input_mant_neg = -input_mant;
      T_in input_float_neg(input_mant_neg, sample_exponent_array[i]);
      if ((input_float_neg.mantissa() != input_mant_neg || input_float_neg.exp() != sample_exponent_array[i]) && mant_i != 0) {
        cout << "input_mant_neg was not normalized correctly." << endl;
        assert(false);
      }
      test_ac_abs_float(input_float_neg, output_float);
      bool correct_iteration_neg = output_check(input_float_neg, output_float);
      correct = correct && correct_iteration_neg;
    }
  }

  if (correct) { printf("PASSED\n"); }
  else         { printf("FAILED\n"); } // LCOV_EXCL_LINE

  all_tests_pass = all_tests_pass && correct;

  return 0;
}

// ==============================================================================
// Function: test_driver_complex()
// Description: test_driver function for ac_fixed inputs and outputs.

template <int Wfi, int Ifi, bool S, int outWfi, int outIfi>
int test_driver_complex(
  double max_error_cmplx,
  const double allowed_error_cmplx,
  const double threshold,
  bool &all_tests_pass,
  bool details = false
)
{
  ac_fixed< Wfi, Ifi, S, AC_TRN, AC_WRAP> input_fixed;
  ac_complex<ac_fixed< Wfi, Ifi, S, AC_TRN, AC_WRAP> > input_complex;
  ac_fixed< outWfi, outIfi, false, AC_TRN, AC_WRAP> output_complex;

  double lower_limit, upper_limit, step, this_error;
  bool correct_iteration;

  // set ranges and step size for fixed point testbench
  lower_limit = input_fixed.template set_val<AC_VAL_MIN>().to_double();
  upper_limit = input_fixed.template set_val<AC_VAL_MAX>().to_double();
  step        = input_fixed.template set_val<AC_VAL_QUANTUM>().to_double();

  cout << "TEST: ac_abs() INPUT: ";
  cout.width(50);
  cout << left << input_complex.type_name();
  cout << "OUTPUT:  ";
  cout.width(38);
  cout << left << output_complex.type_name();
  cout.width(38);
  cout << left << " ";
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl; // LCOV_EXCL_LINE
    cout << "  Ranges for input types:" << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit = " << lower_limit << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit = " << upper_limit << endl; // LCOV_EXCL_LINE
    cout << "    step        = " << step << endl; // LCOV_EXCL_LINE
  }

  bool correct = true;

  // test fixed-point complex values.

  for (double i = lower_limit; i <= upper_limit; i += step) {
    for (double j = lower_limit; j <= upper_limit; j += step) {
      input_complex.r() = i;
      input_complex.i() = j;
      test_ac_abs_complex(input_complex, output_complex);
      double this_error_complex = cmplx_err_calc(input_complex, output_complex, allowed_error_cmplx, threshold);
      if (this_error_complex > max_error_cmplx) {max_error_cmplx = this_error_complex;}
      if (max_error_cmplx < allowed_error_cmplx)
      { correct_iteration = true; }
      else
      { correct_iteration = false; }
      correct = correct && correct_iteration;
    }
  }
  if (correct) { printf("PASSED\n"); }
  else         { printf("FAILED\n"); } // LCOV_EXCL_LINE

  all_tests_pass = all_tests_pass && correct;

  return 0;
}

int main(int argc, char *argv[])
{
  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_abs()" << endl;

  bool all_tests_pass = true;

  // If any of the tests fail, the all_tests_pass variable will be set to false

  double max_error_cmplx = 0;
  double allowed_error_cmplx = 0.5;

  const double threshold = 0.005;

  // template <int Wint, int outWint>
  test_driver_int<7, 7>(all_tests_pass);
  test_driver_int<7, 9>(all_tests_pass);
  test_driver_int<5, 5>(all_tests_pass);
  test_driver_int<5, 8>(all_tests_pass);
  test_driver_int<16, 16>(all_tests_pass);
  test_driver_int<16, 19>(all_tests_pass);

  // template <int Wfi, int Ifi, int outWfi, int outIfi>
  test_driver_fixed<16, 16, 16, 16>(all_tests_pass);
  test_driver_fixed<10,  6, 10,  6>(all_tests_pass);
  test_driver_fixed<11,  8, 14,  9>(all_tests_pass);
  test_driver_fixed<16,  7, 32,  9>(all_tests_pass);
  test_driver_fixed<12, -5, 12, -5>(all_tests_pass);
  test_driver_fixed<12, -5, 15, -2>(all_tests_pass);
  test_driver_fixed< 4,  9,  4,  9>(all_tests_pass);
  test_driver_fixed< 7,  9,  7,  9>(all_tests_pass);
  test_driver_fixed< 7,  9,  8,  9>(all_tests_pass);
  test_driver_fixed< 7,  9, 14, 10>(all_tests_pass);

  // template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
  test_driver_float<12,  5, 10, 12,  5, 11>(all_tests_pass);
  test_driver_float<14,  4,  5, 16,  5,  7>(all_tests_pass);
  test_driver_float<12, -5, 10, 12, -5, 11>(all_tests_pass);
  test_driver_float<12, -5, 10, 13, -4, 12>(all_tests_pass);
  test_driver_float< 4,  9,  5,  4,  9,  6>(all_tests_pass);
  test_driver_float< 4,  9,  5,  6,  9,  6>(all_tests_pass);
  test_driver_float< 4,  9,  5, 12,  9,  6>(all_tests_pass);

  // template <int Wfi, int Ifi, int Sfi, int outWfi, int outIfi>
  test_driver_complex<6, 4, true,  32, 6>(max_error_cmplx, allowed_error_cmplx, threshold, all_tests_pass);
  test_driver_complex<6, 4, false, 32, 7>(max_error_cmplx, allowed_error_cmplx, threshold, all_tests_pass);
  test_driver_complex<12, 5, true,  32, 7>(max_error_cmplx, allowed_error_cmplx, threshold, all_tests_pass);
  test_driver_complex<8, 5, false, 32, 8>(max_error_cmplx, allowed_error_cmplx, threshold, all_tests_pass);
  test_driver_complex<9, 3, true,  32, 5>(max_error_cmplx, allowed_error_cmplx, threshold, all_tests_pass);
  test_driver_complex<9, 3, false, 32, 6>(max_error_cmplx, allowed_error_cmplx, threshold, all_tests_pass);
  test_driver_complex<10, 4, true,  32, 6>(max_error_cmplx, allowed_error_cmplx, threshold, all_tests_pass);
  test_driver_complex<13, 5, false, 32, 7>(max_error_cmplx, allowed_error_cmplx, threshold, all_tests_pass);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished." << endl;

  // Notify the user whether or not the test was a failure.
  if (!all_tests_pass) {
    cout << "  ac_abs - FAILED - Output not correct for all test values" << endl; // LCOV_EXCL_LINE
    cout << "=============================================================================" << endl; // LCOV_EXCL_LINE
    return -1; // LCOV_EXCL_LINE
  } else {
    cout << "  ac_abs - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }

  return 0;

}

