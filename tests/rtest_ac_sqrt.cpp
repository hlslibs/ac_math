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
// ac_sqrt() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_sqrt.cpp -o design
//   ./design

// This testbench tests NaN inputs/outputs as well. Make sure NaN support in ac_sqrt.h
// is enabled by defining the AC_SQRT_NAN_SUPPORTED macro below, and make sure it's
// defined BEFORE including <ac_math/ac_sqrt.h>.
// If you do not define the macro below, NaN testing will be disabled.
#define AC_SQRT_NAN_SUPPORTED

// Include the AC Math function that is exercised with this testbench

#include <ac_math/ac_sqrt.h>
using namespace ac_math;

// ==============================================================================
// Test Designs
//   These simple functions allow executing the ac_sqrt() function
//   using multiple data types at the same time. Template parameters are used to
//   configure the bit-widths of the types.

// Test for ac_int inputs and outputs.
template <int Wint, int outWint>
void test_ac_sqrt_int(
  const ac_int<Wint, false> &in1,
  ac_int<outWint, false> &out1
)
{
  ac_sqrt(in1, out1);
}

// Test for ac_fixed inputs and outputs.
template <int Wfi, int Ifi, int outWfi, int outIfi>
void test_ac_sqrt_fixed(
  const ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP> &in2,
  ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> &out2
)
{
  ac_sqrt(in2, out2);
}

// Test for ac_float inputs and outputs.
template <bool OR_TF, int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
void test_ac_sqrt_float(
  const ac_float<Wfl, Ifl, Efl, AC_TRN> &in3,
  ac_float<outWfl, outIfl, outEfl, AC_TRN> &out3
)
{
  ac_sqrt<OR_TF>(in3, out3);
}

// Test for ac_std_float inputs and outputs.
template <bool OR_TF, int Wstfl, int Estfl, int outWstfl, int outEstfl>
void test_ac_sqrt_stfloat(
  const ac_std_float<Wstfl, Estfl> &in4,
  ac_std_float<outWstfl, outEstfl> &out4
)
{
  ac_sqrt<OR_TF>(in4, out4);
}

// Test for ac_ieee_float inputs and outputs.
template <bool OR_TF, ac_ieee_float_format in_format, ac_ieee_float_format out_format>
void test_ac_sqrt_ifloat(
  const ac_ieee_float<in_format> &in5,
  ac_ieee_float<out_format> &out5
)
{
  ac_sqrt<OR_TF>(in5, out5);
}

// ==============================================================================

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

// ------------------------------------------------------------------------------
// Helper structs for printing out the type name for ac_std_floats

// Generic struct, enables template specialization
template <typename T>
struct type_string_st { };

// Specialized struct, handles ac_std_floats
template <int W, int E>
struct type_string_st<ac_std_float<W, E> > {
  static string type_string() {
    string format_string = "ac_std_float<";
    format_string += ac_int<32,true>(W).to_string(AC_DEC);
    format_string += ",";
    format_string += ac_int<32,true>(E).to_string(AC_DEC);
    format_string += ">";

    return format_string;
  }
};

// Specialized struct, handles ac_ieee_floats
template <ac_ieee_float_format Format>
struct type_string_st<ac_ieee_float<Format> > {
  static string type_string() {
    string format_string = "ac_ieee_float<";
    if (Format == binary16)  { format_string += "binary16"; }
    if (Format == binary32)  { format_string += "binary32"; }
    if (Format == binary64)  { format_string += "binary64"; }
    if (Format == binary128) { format_string += "binary128"; }
    if (Format == binary256) { format_string += "binary256"; }
    format_string += ">";

    return format_string;
  }
};

// -------------------------------------------------------
// Helper functions for output-matching/ error calculation

// See if output is correct for real ac_int inputs.

template <int Wint, int outWint>
bool output_check_int(
  const ac_int<Wint, false> input,
  ac_int<outWint, false> output
)
{
  bool correct = output == (ac_int<outWint, false>)(sqrt(input.to_double()));

  #ifdef DEBUG
  if (!correct) {
    cout << "Output not correct" << endl;
    cout << "input = " << input << endl;
    cout << "output = " << output << endl;
    assert(false);
  }
  #endif

  return correct;
}

// Calculating error for real, ac_fixed datatype.

template <int Wfi, int Ifi, int outWfi, int outIfi>
double err_calc(
  const ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP> input,
  const ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> output,
  const double allowed_error,
  const double threshold
)
{
  // The typecasting is done in order to provide quantization on the expected output.
  double expected_value = ((ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP>)(sqrt(input.to_double()))).to_double();
  double actual_value = output.to_double();
  double this_error;

  // If expected value is greater than a particular threshold, calculate relative error, else, calculate absolute error.
  if (abs(expected_value) > threshold) {
    this_error = abs( (expected_value - actual_value) / expected_value ) * 100.0;
  } else {
    this_error = abs(expected_value - actual_value) * 100.0;
  }

  return this_error;
}

// Calculating error for real, ac_float datatype.

template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
double err_calc(
  const ac_float<Wfl, Ifl, Efl, AC_TRN> input,
  ac_float<outWfl, outIfl, outEfl, AC_TRN> output,
  const double allowed_error,
  const double threshold
)
{
  // The typecasting is done in order to provide quantization on the expected output.
  double expected_value = ((ac_float<outWfl, outIfl, outEfl, AC_TRN>)(sqrt(input.to_double()))).to_double();
  double actual_value = output.to_double();
  double this_error;

  // If expected value is greater than a particular threshold, calculate relative error, else, calculate absolute error.
  if (abs(expected_value) > threshold) {
    this_error = abs( (expected_value - actual_value) / expected_value ) * 100.0;
  } else {
    this_error = abs(expected_value - actual_value) * 100.0;
  }

  return this_error;
}

template<class T_in_acfl, int Wstfl, int Estfl, int outWstfl, int outEstfl>
void err_calc_stfloat(
  const ac_std_float<Wstfl, Estfl> input_stfloat,
  const ac_std_float<outWstfl, outEstfl> output_stfloat,
  const T_in_acfl input_acfl_norm,
  const double allowed_error_stfloat,
  const double threshold,
  double &max_error_stfloat
)
{
  typedef ac_std_float<outWstfl, outEstfl> T_out;

  // The typecasting to ac_std_float<64, 11> is done in order to provide quantization on the expected output.
  double actual_stfloat = ac_std_float<64, 11>(output_stfloat).to_double();
  double expected_stfloat = ac_std_float<64, 11>(T_out(sqrt(input_acfl_norm.to_double()))).to_double();
  double this_error_stfloat;
  // If expected value is greater than a particular threshold, calculate relative error, else, calculate absolute error.
  if (expected_stfloat > threshold) {
    this_error_stfloat = abs((expected_stfloat - actual_stfloat)/expected_stfloat)*100.0;
  } else {
    this_error_stfloat = abs(expected_stfloat - actual_stfloat)*100.0;
  }

  if (this_error_stfloat > max_error_stfloat) {
    max_error_stfloat = this_error_stfloat;
  }

  #ifdef DEBUG
  if (this_error_stfloat > allowed_error_stfloat) {
    cout << endl;
    cout << "  Error exceeds tolerance" << endl;
    cout << "  input_acfl_norm   = " << input_acfl_norm << endl;
    cout << "  input_stfloat      = " << input_stfloat.to_ac_float().to_double() << endl;
    cout << "  expected_stfloat   = " << expected_stfloat << endl;
    cout << "  actual_stfloat     = " << actual_stfloat << endl;
    cout << "  this_error_stfloat = " << this_error_stfloat << endl;
    assert(false);
  }
  #endif
}

template<class T_in_acfl, ac_ieee_float_format in_format, ac_ieee_float_format out_format>
void err_calc_ifloat(
  const ac_ieee_float<in_format> input_ifloat,
  const ac_ieee_float<out_format> output_ifloat,
  const T_in_acfl input_acfl_norm,
  const double allowed_error_ifloat,
  const double threshold,
  double &max_error_ifloat
)
{
  typedef ac_ieee_float<out_format> T_out;

  // The typecasting to ac_ieee_float<binary64> is done in order to provide quantization on the expected output.
  double actual_ifloat = ac_ieee_float<binary64>(output_ifloat).to_double();
  double expected_ifloat = ac_ieee_float<binary64>(T_out(sqrt(input_acfl_norm.to_double()))).to_double();
  double this_error_ifloat;
  // If expected value is greater than a particular threshold, calculate relative error, else, calculate absolute error.
  if (expected_ifloat > threshold) {
    this_error_ifloat = abs((expected_ifloat - actual_ifloat)/expected_ifloat)*100.0;
  } else {
    this_error_ifloat = abs(expected_ifloat - actual_ifloat)*100.0;
  }

  if (this_error_ifloat > max_error_ifloat) {
    max_error_ifloat = this_error_ifloat;
  }

  #ifdef DEBUG
  if (this_error_ifloat > allowed_error_ifloat) {
    cout << endl;
    cout << "  Error exceeds tolerance" << endl;
    cout << "  input_acfl_norm   = " << input_acfl_norm << endl;
    cout << "  input_ifloat      = " << input_ifloat.to_ac_float().to_double() << endl;
    cout << "  expected_ifloat   = " << expected_ifloat << endl;
    cout << "  actual_ifloat     = " << actual_ifloat << endl;
    cout << "  this_error_ifloat = " << this_error_ifloat << endl;
    assert(false);
  }
  #endif
}

// ==============================================================================
// Function: test_driver_int()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_int values. It uses the type information to iterate through a
//   range of valid values on that type in order to make sure that the
//   output of the ac_sqrt function is as expected.

template <int Wint, int outWint>
int test_driver_int(
  bool &all_tests_pass
)
{
  ac_int<Wint, false> input_int;
  ac_int<outWint, false> output_int;

  cout << "TEST: ac_sqrt() INPUT: ";
  cout.width(38);
  cout << left << input_int.type_name();
  cout << "OUTPUT: ";
  cout.width(38);
  cout << left << output_int.type_name();
  cout << "RESULT: ";

  bool correct = true;

  // set ranges and step size for integer testbench
  double lower_limit_int   = input_int.template set_val<AC_VAL_MIN>().to_double();
  double upper_limit_int   = input_int.template set_val<AC_VAL_MAX>().to_double();
  double step_int          = input_int.template set_val<AC_VAL_QUANTUM>().to_double();

  // Test integer inputs
  for (double i = lower_limit_int; i <= upper_limit_int; i += step_int) {
    // Set values for input.
    input_int = i;
    test_ac_sqrt_int(input_int, output_int);
    bool correct_iteration = output_check_int(input_int, output_int);
    correct = correct_iteration && correct;
  }

  if (correct) { printf("PASSED\n"); }
  else         { printf("FAILED\n"); } // LCOV_EXCL_LINE

  all_tests_pass = all_tests_pass && correct;

  return 0;
}

// ==============================================================================
// Function: test_driver_fixed()
// Description: A templatized function that can be configured for certain bit-
//   widths of unsigned ac_fixed values. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   output of the ac_sqrt() function with the computed square root output
//   using a standard C double type. The maximum error for each type is
//   accumulated in variables defined in the calling function.

template <int Wfi, int Ifi, int outWfi, int outIfi>
int test_driver_fixed(
  double &cumulative_max_error_fixed,
  const double allowed_error_fixed,
  const double threshold,
  bool details = false
)
{
  bool passed;
  double max_error_fixed = 0.0; // reset for this run

  ac_fixed<   Wfi,    Ifi, false, AC_TRN, AC_WRAP> input_fixed;
  ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> output_fixed;

  double lower_limit_fixed, upper_limit_fixed, step_fixed;

  // set ranges and step size for fixed point testbench
  lower_limit_fixed   = input_fixed.template set_val<AC_VAL_MIN>().to_double();
  upper_limit_fixed   = input_fixed.template set_val<AC_VAL_MAX>().to_double();
  step_fixed          = input_fixed.template set_val<AC_VAL_QUANTUM>().to_double();

  cout << "TEST: ac_sqrt() INPUT: ";
  cout.width(38);
  cout << left << input_fixed.type_name();
  cout << "OUTPUT: ";
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

  // Test unsigned fixed point inputs
  for (double i = lower_limit_fixed; i <= upper_limit_fixed; i += step_fixed) {
    // Set values for input.
    input_fixed = i;
    test_ac_sqrt_fixed(input_fixed, output_fixed);
    double this_error_fixed = err_calc(input_fixed, output_fixed, allowed_error_fixed, threshold);
    if (this_error_fixed > max_error_fixed) {max_error_fixed = this_error_fixed;}
  }

  passed = max_error_fixed < allowed_error_fixed;

  if (passed) { printf("PASSED, max err (%f)\n", max_error_fixed); }
  else        { printf("FAILED, max err (%f)\n", max_error_fixed); } // LCOV_EXCL_LINE

  if (max_error_fixed>cumulative_max_error_fixed) { cumulative_max_error_fixed = max_error_fixed; }

  return 0;
}

// ==============================================================================
// Function: test_driver_float()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_float values. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   output of the ac_sqrt() function with the computed square root output
//   using a standard C double type. The maximum error for each type is
//   accumulated in variables defined in the calling function.

template <bool OR_TF, int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
int test_driver_float(
  double &cumulative_max_error_float,
  const double allowed_error_float,
  const double threshold,
  bool details = false
)
{
  double max_error_float = 0.0; // reset for this run

  typedef ac_float<   Wfl,    Ifl,    Efl> T_in;
  typedef ac_float<outWfl, outIfl, outEfl> T_out;

  // Since ac_float values are normalized, the bit adjacent to the sign bit in the mantissa
  // will always be set to 1. We will hence cycle through all the bit patterns that correspond to the last (Wfl - 2)
  // bits in the mantissa.
  ac_int<Wfl - 2, false> sample_mantissa_slc;
  // Set the lower limit, upper limit and step size of the test iterations.
  ac_int<Wfl - 2, false> lower_limit_it = 0; // Set to zero because only positive inputs are supported.
  ac_int<Wfl - 2, false> upper_limit_it = sample_mantissa_slc.template set_val<AC_VAL_MAX>().to_double();
  ac_int<Wfl - 2, false> step_it = 1; // Since sample_mantissa_slc is an integer.

  cout << "TEST: ac_sqrt() INPUT: ";
  cout.width(32);
  cout << left << T_in::type_name();
  cout << "OUTPUT: ";
  cout.width(32);
  cout << left << T_out::type_name();
  cout << "OR_TF = ";
  cout.width(7);
  cout << left << (OR_TF ? "true" : "false");
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

  // Alternate between odd and even values for the other elements. The even values are powers of two,
  // while the odd values are one less than the nearest power of two. Exponent values = +/-1 are left
  // as they are.
  bool odd_elem = true;
  for (int i = (Efl - 2); i >= 0; i--) {
    sample_exponent_value = 0;
    sample_exponent_value[i] = 1;
    sample_exponent_array[Efl + i + 1] = sample_exponent_value - int(odd_elem && i != 0);
    sample_exponent_array[Efl - i - 1] = -(sample_exponent_value - int(odd_elem && i != 0));
    odd_elem = !odd_elem;
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
      ac_fixed<Wfl, Ifl, true> input_mant = 0;
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
      test_ac_sqrt_float<OR_TF>(input_float, output_float);
      double this_error_float = err_calc(input_float, output_float, allowed_error_float, threshold);
      if (this_error_float > max_error_float) { max_error_float = this_error_float; }
    }
  }

  bool passed = (max_error_float < allowed_error_float);

  if (passed) { printf("PASSED , max err (%f) \n", max_error_float); }
  else        { printf("FAILED , max err (%f) \n", max_error_float); } // LCOV_EXCL_LINE

  if (max_error_float > cumulative_max_error_float) { cumulative_max_error_float = max_error_float; }

  return 0;
}

// ==============================================================================
// Function: test_driver_stfloat()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_std_float values. It uses the type information to iterate
//   through a range of valid values on that type in order to compare the
//   precision of the output of the ac_sqrt() function with the computed square
//   root output using a standard C double type. The maximum error for each type
//   is accumulated in variables defined in the calling function.

template <bool OR_TF, int Wstfl, int Estfl, int outWstfl, int outEstfl>
int test_driver_stfloat(
  double &cumulative_max_error_stfloat,
  const double allowed_error_stfloat,
  const double threshold,
  bool &all_tests_pass,
  bool details = false
)
{
  // In order to insure correct quantization in err_calc_stfloat, the input bitwidth
  // is limited to 64 or less, while the input exponent width is limited to 11 or less.
  static_assert(Wstfl <= 64, "Input ac_std_float bitwidth must not be greater than 64.");
  static_assert(Estfl <= 11, "Input ac_std_float exponent width must not be greater than 11.");

  double max_error_stfloat = 0.0; // reset for this run

  typedef ac_std_float<Wstfl, Estfl> T_in;
  typedef ac_std_float<outWstfl, outEstfl> T_out;

  const int W2acfl = Wstfl - Estfl + 1;
  const int I2acfl = 2;
  typedef ac_float<W2acfl, I2acfl, Estfl> T_in_acfl;

  // Set the lower limit, upper limit and step size of the test iterations.
  ac_int<W2acfl - 2, false> sample_mantissa_slc_stfloat;
  ac_int<W2acfl - 2, false> lower_limit_stfloat = 0; // Set to zero because only +ve values are supported.
  ac_int<W2acfl - 2, false> upper_limit_stfloat = sample_mantissa_slc_stfloat.template set_val<AC_VAL_MAX>().to_double();
  ac_int<W2acfl - 2, false> step_stfloat = 0;
  // For larger ac_std_float bitwidths, the step size is increased by increasing the number of non-zero
  // LSBs, so as to reduce the time taken.
  for (int i = 0; i < AC_MAX(W2acfl - 2 - 16, 1); i++) { step_stfloat[i] = 1; }

  cout << "TEST: ac_sqrt() INPUT: ";
  cout.width(32);
  cout << left << type_string_st<T_in>::type_string();
  cout << "OUTPUT: ";
  cout.width(32);
  cout << left << type_string_st<T_out>::type_string();
  cout << "OR_TF = ";
  cout.width(7);
  cout << left << (OR_TF ? "true" : "false");
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl << "  Ranges for testing iterations:" << endl; // LCOV_EXCL_LINE
    cout         << "    lower_limit_stfloat   = " << lower_limit_stfloat << endl; // LCOV_EXCL_LINE
    cout         << "    upper_limit_stfloat   = " << upper_limit_stfloat << endl; // LCOV_EXCL_LINE
    cout         << "    step_stfloat          = " << step_stfloat << endl; // LCOV_EXCL_LINE
    cout         << "    allowed_error_stfloat = " << allowed_error_stfloat << endl; // LCOV_EXCL_LINE
  }

  // sample_exponent_array_stfloat stores all values of ac_std_float exponent to be tested.
  const int exp_arr_size_stfloat = 2*(Estfl - 1) + 3;
  ac_int<Estfl, true> sample_exponent_array_stfloat[exp_arr_size_stfloat];
  ac_int<Estfl, true> sample_exponent_value_stfloat;

  const int bias = (1 << (Estfl - 1)) - 1;
  // The first element of the array is the minimum exponent value, the middle element is a zero exponent, and
  // the last element is the maximum possible value.
  sample_exponent_array_stfloat[0]                       = -(bias - 1);
  sample_exponent_array_stfloat[Estfl]                    = 0;
  sample_exponent_array_stfloat[exp_arr_size_stfloat - 1] = bias;

  // Alternate between odd and even values for the other elements. The even values are powers of two,
  // while the odd values are one less than the nearest power of two. Exponent values = +/-1 are left
  // as they are.
  bool odd_elem = true;

  for (int i = (Estfl - 2); i >= 0; i--) {
    sample_exponent_value_stfloat = 0;
    sample_exponent_value_stfloat[i] = 1;
    sample_exponent_array_stfloat[Estfl + i + 1] = sample_exponent_value_stfloat - int(odd_elem&&i!=0);
    sample_exponent_array_stfloat[Estfl - i - 1] = -(sample_exponent_value_stfloat - int(odd_elem&&i!=0));
    odd_elem = !odd_elem;
  }

  // For large ac_std_float bitwidths, an increment value > 1 might lead to the testbench not testing the largest possible
  // input mantissa value. In such a case, we test the max. mantissa value separately.
  bool test_max_mant_val_sep = true;

  for (int i = 0; i < exp_arr_size_stfloat; i++) {
    for (ac_int<W2acfl - 1, false> mant_i = lower_limit_stfloat; mant_i <= upper_limit_stfloat; mant_i += step_stfloat) {
      test_max_mant_val_sep = test_max_mant_val_sep && mant_i != upper_limit_stfloat;
      ac_fixed<W2acfl, I2acfl, true> input_mant = 0;
      input_mant[W2acfl - 1] = 0;
      input_mant[W2acfl - 2] = 1; // Bit next to sign bit is set to 1 to ensure normalization.
      input_mant.set_slc(0, mant_i.template slc<W2acfl - 2>(0));
      T_in_acfl input_acfl_norm(input_mant, sample_exponent_array_stfloat[i]);
      // If mantissa and exponent values have changed after normalization, that means that we didn't normalize the
      // input value correctly to begin with.
      if (input_acfl_norm.mantissa() != input_mant || input_acfl_norm.exp() != sample_exponent_array_stfloat[i]) {
        cout << "Input value was not normalized correctly." << endl;
        assert(false);
      }
      // Assign normalized float value to input ac_std_float.
      T_in input_stfloat(input_acfl_norm);
      T_out output_stfloat;
      test_ac_sqrt_stfloat<OR_TF>(input_stfloat, output_stfloat);
      // Update max. error if necessary.
      err_calc_stfloat(input_stfloat, output_stfloat, input_acfl_norm, allowed_error_stfloat, threshold, max_error_stfloat);
    }
    if (test_max_mant_val_sep) {
      // Test maximum mantissa value separately.
      ac_fixed<W2acfl, I2acfl, true> input_mant;
      input_mant.template set_val<AC_VAL_MAX>();
      T_in_acfl input_acfl_norm(input_mant, sample_exponent_array_stfloat[i]);
      // Assign normalized float value to input ac_std_float.
      T_in input_stfloat(input_acfl_norm);
      T_out output_stfloat;
      test_ac_sqrt_stfloat<OR_TF>(input_stfloat, output_stfloat);
      // Update max. error if necessary.
      err_calc_stfloat(input_stfloat, output_stfloat, input_acfl_norm, allowed_error_stfloat, threshold, max_error_stfloat);
    }
  }

  if (max_error_stfloat > cumulative_max_error_stfloat) { cumulative_max_error_stfloat = max_error_stfloat; }

  bool passed = (max_error_stfloat < allowed_error_stfloat);

  T_in input_neg_zero = -T_in::zero();
  T_out output_neg_zero;
  test_ac_sqrt_stfloat<OR_TF>(input_neg_zero, output_neg_zero);

  if (output_neg_zero != T_out::zero() || !output_neg_zero.signbit()) {
    printf("Negative Zero testing FAILED\n");
    all_tests_pass = false;
    return -1;
  }

  #ifdef AC_SQRT_NAN_SUPPORTED
  // Test if negative input (-1.0) gives a -nan output.
  T_in input_neg = -T_in::one();
  T_out output_nan;
  test_ac_sqrt_stfloat<OR_TF>(input_neg, output_nan);
  bool nan_testing_passed = output_nan.isnan() && output_nan.signbit();

  // +nan input should give +nan output.
  T_in input_nan = T_in::nan();
  test_ac_sqrt_stfloat<OR_TF>(input_nan, output_nan);
  nan_testing_passed = nan_testing_passed && output_nan.isnan() && !output_nan.signbit();

  // -nan input should give -nan output.
  input_nan = -T_in::nan();
  test_ac_sqrt_stfloat<OR_TF>(input_nan, output_nan);
  nan_testing_passed = nan_testing_passed && output_nan.isnan() && output_nan.signbit();

  if (!nan_testing_passed) {
    printf("NaN testing FAILED\n");
    all_tests_pass = false;
    return -1;
  }
  #endif

  all_tests_pass = all_tests_pass && passed;

  if (passed) { printf("PASSED , max err (%f) \n", max_error_stfloat); }
  else        { printf("FAILED , max err (%f) \n", max_error_stfloat); } // LCOV_EXCL_LINE

  return 0;
}

// ==============================================================================
// Function: test_driver_ifloat()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_ieee_float values. It uses the type information to iterate
//   through a range of valid values on that type in order to compare the
//   precision of the output of the ac_sqrt() function with the computed square
//   root output using a standard C double type. The maximum error for each type
//   is accumulated in variables defined in the calling function.

template <bool OR_TF, ac_ieee_float_format in_format, ac_ieee_float_format out_format>
int test_driver_ifloat(
  double &cumulative_max_error_ifloat,
  const double allowed_error_ifloat,
  const double threshold,
  bool &all_tests_pass,
  bool details = false
)
{
  // In order to insure correct quantization in err_calc_ifloat, we do not accept input
  // formats of binary128 or binary256.
  static_assert(in_format != binary128 && in_format != binary256, "binary128 and binary256 input ac_ieee_float formats not supported.");

  double max_error_ifloat = 0.0; // reset for this run

  typedef ac_ieee_float<in_format> T_in;
  typedef ac_ieee_float<out_format> T_out;

  const int Wifl = T_in::width;
  const int Eifl = T_in::e_width;
  const int W2acfl = Wifl - Eifl + 1;
  const int I2acfl = 2;
  typedef ac_float<W2acfl, I2acfl, Eifl> T_in_acfl;

  // Set the lower limit, upper limit and step size of the test iterations.
  ac_int<W2acfl - 2, false> sample_mantissa_slc_ifloat;
  ac_int<W2acfl - 2, false> lower_limit_ifloat = 0; // Set to zero because only +ve values are supported.
  ac_int<W2acfl - 2, false> upper_limit_ifloat = sample_mantissa_slc_ifloat.template set_val<AC_VAL_MAX>().to_double();
  ac_int<W2acfl - 2, false> step_ifloat = 0;
  // For larger ac_ieee_float bitwidths, the step size is increased by increasing the number of non-zero
  // LSBs, so as to reduce the time taken.
  for (int i = 0; i < AC_MAX(W2acfl - 2 - 16, 1); i++) { step_ifloat[i] = 1; }

  cout << "TEST: ac_sqrt() INPUT: ";
  cout.width(32);
  cout << left << type_string_st<T_in>::type_string();
  cout << "OUTPUT: ";
  cout.width(32);
  cout << left << type_string_st<T_out>::type_string();
  cout << "OR_TF = ";
  cout.width(7);
  cout << left << (OR_TF ? "true" : "false");
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl << "  Ranges for testing iterations:" << endl; // LCOV_EXCL_LINE
    cout         << "    lower_limit_ifloat   = " << lower_limit_ifloat << endl; // LCOV_EXCL_LINE
    cout         << "    upper_limit_ifloat   = " << upper_limit_ifloat << endl; // LCOV_EXCL_LINE
    cout         << "    step_ifloat          = " << step_ifloat << endl; // LCOV_EXCL_LINE
    cout         << "    allowed_error_ifloat = " << allowed_error_ifloat << endl; // LCOV_EXCL_LINE
  }

  // sample_exponent_array_ifloat stores all values of ac_ieee_float exponent to be tested.
  const int exp_arr_size_ifloat = 2*(Eifl - 1) + 3;
  ac_int<Eifl, true> sample_exponent_array_ifloat[exp_arr_size_ifloat];
  ac_int<Eifl, true> sample_exponent_value_ifloat;

  const int bias = (1 << (Eifl - 1)) - 1;
  // The first element of the array is the minimum exponent value, the middle element is a zero exponent, and
  // the last element is the maximum possible value.
  sample_exponent_array_ifloat[0]                       = -(bias - 1);
  sample_exponent_array_ifloat[Eifl]                    = 0;
  sample_exponent_array_ifloat[exp_arr_size_ifloat - 1] = bias;

  // Alternate between odd and even values for the other elements. The even values are powers of two,
  // while the odd values are one less than the nearest power of two. Exponent values = +/-1 are left
  // as they are.
  bool odd_elem = true;

  // All the other elements are set to values that correspond to a one-hot encoding scheme, in which only one
  // bit of the absolute value of the exponent is set to one. Both negative and positive values are encoded this way.
  for (int i = (Eifl - 2); i >= 0; i--) {
    sample_exponent_value_ifloat = 0;
    sample_exponent_value_ifloat[i] = 1;
    sample_exponent_array_ifloat[Eifl + i + 1] = sample_exponent_value_ifloat - int(odd_elem&&i!=0);
    sample_exponent_array_ifloat[Eifl - i - 1] = -(sample_exponent_value_ifloat - int(odd_elem&&i!=0));
    odd_elem = !odd_elem;
  }

  // For large ac_ieee_float bitwidths, an increment value > 1 might lead to the testbench not testing the largest possible
  // input mantissa value. In such a case, we test the max. mantissa value separately.
  bool test_max_mant_val_sep = true;

  for (int i = 0; i < exp_arr_size_ifloat; i++) {
    for (ac_int<W2acfl - 1, false> mant_i = lower_limit_ifloat; mant_i <= upper_limit_ifloat; mant_i += step_ifloat) {
      test_max_mant_val_sep = test_max_mant_val_sep && mant_i != upper_limit_ifloat;
      ac_fixed<W2acfl, I2acfl, true> input_mant = 0;
      input_mant[W2acfl - 1] = 0;
      input_mant[W2acfl - 2] = 1; // Bit next to sign bit is set to 1 to ensure normalization.
      input_mant.set_slc(0, mant_i.template slc<W2acfl - 2>(0));
      T_in_acfl input_acfl_norm(input_mant, sample_exponent_array_ifloat[i]);
      // If mantissa and exponent values have changed after normalization, that means that we didn't normalize the
      // input value correctly to begin with.
      if (input_acfl_norm.mantissa() != input_mant || input_acfl_norm.exp() != sample_exponent_array_ifloat[i]) {
        cout << "Input value was not normalized correctly." << endl;
        assert(false);
      }
      // Assign normalized float value to input ac_ieee_float.
      T_in input_ifloat(input_acfl_norm);
      T_out output_ifloat;
      test_ac_sqrt_ifloat<OR_TF>(input_ifloat, output_ifloat);
      // Update max. error if necessary.
      err_calc_ifloat(input_ifloat, output_ifloat, input_acfl_norm, allowed_error_ifloat, threshold, max_error_ifloat);
    }
    if (test_max_mant_val_sep) {
      // Test maximum mantissa value separately.
      ac_fixed<W2acfl, I2acfl, true> input_mant;
      input_mant.template set_val<AC_VAL_MAX>();
      T_in_acfl input_acfl_norm(input_mant, sample_exponent_array_ifloat[i]);
      // Assign normalized float value to input ac_ieee_float.
      T_in input_ifloat(input_acfl_norm);
      T_out output_ifloat;
      test_ac_sqrt_ifloat<OR_TF>(input_ifloat, output_ifloat);
      // Update max. error if necessary.
      err_calc_ifloat(input_ifloat, output_ifloat, input_acfl_norm, allowed_error_ifloat, threshold, max_error_ifloat);
    }
  }

  if (max_error_ifloat > cumulative_max_error_ifloat) { cumulative_max_error_ifloat = max_error_ifloat; }

  bool passed = (max_error_ifloat < allowed_error_ifloat);

  // Start special input testing.

  // -0.0 input should give -0.0 output.
  T_in input_neg_zero = -T_in::zero();
  T_out output_neg_zero;
  test_ac_sqrt_ifloat<OR_TF>(input_neg_zero, output_neg_zero);

  if (output_neg_zero != T_out::zero() || !output_neg_zero.signbit()) {
    printf("Negative Zero testing FAILED\n");
    all_tests_pass = false;
    return -1;
  }

  #ifdef AC_SQRT_NAN_SUPPORTED
  // Test if negative input (-1.0) gives a -nan output.
  T_in input_neg = -T_in::one();
  T_out output_nan;
  test_ac_sqrt_ifloat<OR_TF>(input_neg, output_nan);
  bool nan_testing_passed = output_nan.isnan() && output_nan.signbit();

  // +nan input should give +nan output.
  T_in input_nan = T_in::nan();
  test_ac_sqrt_ifloat<OR_TF>(input_nan, output_nan);
  nan_testing_passed = nan_testing_passed && output_nan.isnan() && !output_nan.signbit();

  // -nan input should give -nan output.
  input_nan = -T_in::nan();
  test_ac_sqrt_ifloat<OR_TF>(input_nan, output_nan);
  nan_testing_passed = nan_testing_passed && output_nan.isnan() && output_nan.signbit();

  if (!nan_testing_passed) {
    printf("NaN testing FAILED\n");
    all_tests_pass = false;
    return -1;
  }
  #endif

  // End special input testing.

  all_tests_pass = all_tests_pass && passed;

  if (passed) { printf("PASSED , max err (%f) \n", max_error_ifloat); }
  else        { printf("FAILED , max err (%f) \n", max_error_ifloat); } // LCOV_EXCL_LINE

  return 0;
}

int main(int argc, char *argv[])
{
  double max_error_fixed = 0.0, max_error_float = 0.0, max_error_stfloat = 0.0, max_error_ifloat = 0.0;

  // Set tolerance
  double allowed_error_fixed = 0.005, allowed_error_float = 0.005, allowed_error_stfloat = 0.005, allowed_error_ifloat = 0.005;

  // threshold below which we calculate absolute error instead of relative for fixed point
  double threshold_fixed = 0.00005, threshold_float = 0.00005, threshold_stfloat = 0.00005, threshold_ifloat = 0.00005;

  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_sqrt() - Allowed error " << allowed_error_fixed << " (fixed pt) " << allowed_error_float << " (float pt) " << allowed_error_float << " (std float pt)" << endl;

  bool all_tests_pass = true;

  // If any of the tests fail for ac_int, the all_tests_pass variable will be set to false
  // template <int Wint, int outWint>
  test_driver_int<13, 24>(all_tests_pass);
  test_driver_int<10, 24>(all_tests_pass);
  test_driver_int<16, 32>(all_tests_pass);
  test_driver_int< 8, 16>(all_tests_pass);
  test_driver_int<11, 16>(all_tests_pass);

  // template <int Wfi, int Ifi, int outWfi, int outIfi>
  test_driver_fixed< 12,  0, 64, 32>(max_error_fixed, allowed_error_fixed, threshold_fixed);
  test_driver_fixed<  8, -2, 64, 32>(max_error_fixed, allowed_error_fixed, threshold_fixed);
  test_driver_fixed< 10, -4, 64, 32>(max_error_fixed, allowed_error_fixed, threshold_fixed);
  test_driver_fixed< 12,  8, 64, 32>(max_error_fixed, allowed_error_fixed, threshold_fixed);
  test_driver_fixed< 16,  8, 64, 32>(max_error_fixed, allowed_error_fixed, threshold_fixed);
  test_driver_fixed<  9,  4, 60, 30>(max_error_fixed, allowed_error_fixed, threshold_fixed);
  test_driver_fixed<  4,  9, 60, 30>(max_error_fixed, allowed_error_fixed, threshold_fixed);
  test_driver_fixed< 12, 12, 64, 32>(max_error_fixed, allowed_error_fixed, threshold_fixed);
  test_driver_fixed< 12,  5, 64, 32>(max_error_fixed, allowed_error_fixed, threshold_fixed);

  // template <bool OR_TF, int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
  test_driver_float<false, 18, 12, 10, 32,  1, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<false, 18, 11, 10, 32,  1, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<false, 18,  2, 10, 32,  1, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<false, 18,  1, 10, 32,  1, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<false, 18,  0, 10, 32,  1, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<true,  18, -4, 10, 32,  1, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<true,  18, -5, 10, 32,  1, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<true,  18, 12, 10, 32,  6, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<true,  18, 11, 10, 32,  6, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<true,  18,  2, 10, 32,  6, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<true,  18,  1, 10, 32,  6, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<true,  18,  0, 10, 32,  6, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<true,  18, -4, 10, 32,  6, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<true,  18, -5, 10, 32,  6, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<false, 18, 12, 10, 32, -4, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<false, 18, 11, 10, 32, -4, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<false, 18,  2, 10, 32, -4, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<false, 18,  1, 10, 32, -4, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<false, 18,  0, 10, 32, -4, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<false, 18, -4, 10, 32, -4, 15>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<false, 18, -5, 10, 32, -4, 15>(max_error_float, allowed_error_float, threshold_float);

  // template <bool OR_TF, int Wstfl, int Estfl, int outWstfl, int outEstfl>
  test_driver_stfloat<true, 16,  5, 48,  9>(max_error_stfloat, allowed_error_stfloat, threshold_stfloat, all_tests_pass);
  test_driver_stfloat<true, 24,  7, 48,  9>(max_error_stfloat, allowed_error_stfloat, threshold_stfloat, all_tests_pass);
  test_driver_stfloat<true, 32,  8, 48,  9>(max_error_stfloat, allowed_error_stfloat, threshold_stfloat, all_tests_pass);
  test_driver_stfloat<true, 48,  9, 48,  9>(max_error_stfloat, allowed_error_stfloat, threshold_stfloat, all_tests_pass);
  test_driver_stfloat<true, 48, 11, 64, 11>(max_error_stfloat, allowed_error_stfloat, threshold_stfloat, all_tests_pass);

  // template <bool OR_TF, ac_ieee_float_format in_format, ac_ieee_float_format out_format>
  test_driver_ifloat<true, binary16, binary64>(max_error_ifloat, allowed_error_ifloat, threshold_ifloat, all_tests_pass);
  test_driver_ifloat<true, binary32, binary64>(max_error_ifloat, allowed_error_ifloat, threshold_ifloat, all_tests_pass);
  test_driver_ifloat<true, binary64, binary64>(max_error_ifloat, allowed_error_ifloat, threshold_ifloat, all_tests_pass);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum error observed across all data type/bit-width variations:" << endl;
  cout << "    max_error_fixed   = " << max_error_fixed << endl;
  cout << "    max_error_float   = " << max_error_float << endl;
  cout << "    max_error_stfloat = " << max_error_stfloat << endl;
  cout << "    max_error_ifloat  = " << max_error_ifloat << endl;

  // If error limits on any ac_fixed/ac_float/ac_std_float/ac_ieee_float test value have been crossed, or the output
  // for ac_int datatypes is not correct, the all_test_pass flag is set to false.
  all_tests_pass = (max_error_fixed > allowed_error_fixed) || (max_error_float > allowed_error_float) || (!all_tests_pass);

  // Notify the user that the test was a failure.
  if (all_tests_pass) {
    cout << "  ac_sqrt - FAILED - Error tolerance(s) exceeded and/or special input Testing Failed." << endl; // LCOV_EXCL_LINE
    cout << "=============================================================================" << endl; // LCOV_EXCL_LINE
    return -1; // LCOV_EXCL_LINE
  } else {
    cout << "  ac_sqrt - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }

  return 0;
}
