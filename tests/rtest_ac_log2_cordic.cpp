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
// ac_log2_cordic() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_log2_cordic.cpp -o design
//   ./design

// This testbench tests NaN inputs/outputs as well. Make sure NaN support in ac_hcordic.h
// is enabled by defining the AC_HCORDIC_NAN_SUPPORTED macro below, and make sure it's
// defined BEFORE including <ac_math/ac_hcordic.h>.
// If you do not define the macro below, NaN testing will be disabled.
#define AC_HCORDIC_NAN_SUPPORTED

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_hcordic.h>
using namespace ac_math;

// ==============================================================================
// Test Design
//   This simple function allows executing the ac_log2_cordic() function.
//   Template parameters are used to configure the bit-widths of the
//   ac_fixed inputs.

template <int Wfi, int Ifi, int outWfi, int outIfi, bool outSfi>
void test_ac_log2_cordic_fixed(
  const ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP>  &in,
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> &out
)
{
  ac_log2_cordic(in, out);
}

// Test design for real floating point values.
template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
void test_ac_log2_cordic_float(
  const ac_float<Wfl, Ifl, Efl, AC_TRN> &in,
  ac_float<outWfl, outIfl, outEfl, AC_TRN> &out
)
{
  ac_log2_cordic(in, out);
}

// Test for ac_std_float input and output.
template <bool OR_TF, int Wstfl, int Estfl, int outWstfl, int outEstfl>
void test_ac_log2_cordic_stfloat(
  const ac_std_float<Wstfl, Estfl> &in,
  ac_std_float<outWstfl, outEstfl> &out
)
{
  ac_log2_cordic<OR_TF>(in, out);
}

// Test for ac_ieee_float input and output.
template <bool OR_TF, ac_ieee_float_format in_format, ac_ieee_float_format out_format>
void test_ac_log2_cordic_ifloat(
  const ac_ieee_float<in_format> &in,
  ac_ieee_float<out_format> &out
)
{
  ac_log2_cordic<OR_TF>(in, out);
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

// ==============================================================================
// Functions: test_driver functions
// Description: Templatized functions that can be configured for certain bit-
//   widths of AC datatypes. They use the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   cordic table base 2 logarithm model with the computed base 2 logarithm
//   using a C++ math library with the standard C++ double type. The maximum
//   error for each type is accumulated in variables defined in the calling function.

// ===============================================================================
// Function: test_driver_fixed()
// Description: test_driver function for ac_fixed inputs and outputs.

template <int Wfi, int Ifi, int outWfi, int outIfi, bool outSfi>
int test_driver_fixed(
  double &cumulative_max_error_fixed,
  const double allowed_error,
  bool details = false
)
{
  ac_fixed<Wfi + 2, Ifi + 1, false, AC_TRN, AC_WRAP> i; // make loop variable slightly larger
  ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP> input;
  ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP> last;
  typedef ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> T_out;
  T_out log2_out;

  // set ranges and step size for testbench
  ac_fixed<Wfi, Ifi, false> lower_limit;
  // If output is unsigned, make sure that the input is always greater than or equal to 1
  // by setting the lower_limit of the testbench to 1 or AC_VAL_QUANTUM, whichever is
  // greater. This is because input values lesser than 1 correspond to a negative output,
  // which can't be stored in unsigned variables.
  if(outSfi) { lower_limit = input.template set_val<AC_VAL_MIN>(); }
  else {
    input.template set_val<AC_VAL_QUANTUM>();
    lower_limit = AC_MAX(1.0, input);
  }
  ac_fixed<Wfi, Ifi, false> upper_limit = input.template set_val<AC_VAL_MAX>();
  ac_fixed<Wfi, Ifi, false> step        = input.template set_val<AC_VAL_QUANTUM>();

  cout << "TEST: ac_log2_cordic() INPUT: ";
  cout.width(38);
  cout << left << input.type_name();
  cout << "OUTPUT: ";
  cout.width(38);
  cout << left << log2_out.type_name();
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl; // LCOV_EXCL_LINE
    cout << "  Ranges for input types:" << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit          = " << lower_limit << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit          = " << upper_limit << endl; // LCOV_EXCL_LINE
    cout << "    step                 = " << step << endl; // LCOV_EXCL_LINE
  }

  bool passed = true;
  double max_error_fixed = 0.0;

  for (i = lower_limit; i <= upper_limit; i += step) {
    // Set values for input.
    input = i;
    if (input.to_double() == 0) { continue; }

    // call reference log2() with fixed-pt value converted back to double
    // an additional step of typecasting is required in order to perform
    // quantization on the expected output.
    double expected_fixed = ((T_out)log2(input.to_double())).to_double();

    // call DUT with fixed-pt value
    test_ac_log2_cordic_fixed(input,log2_out);

    double actual_fixed = log2_out.to_double();

    // Since denormalization for ac_log2_cordic is always done via addition and not shifting, absolute
    // error is the best metric to use.
    double this_error_fixed = abs_double(expected_fixed - actual_fixed) * 100.0;

#ifdef DEBUG
    if (this_error_fixed > allowed_error) {
      cout << endl;
      cout << "  Error exceeds tolerance" << endl;
      cout << "  input            = " << input << endl;
      cout << "  expected_fixed   = " << expected_fixed << endl;
      cout << "  actual_fixed     = " << actual_fixed << endl;
      cout << "  this_error_fixed = " << this_error_fixed << endl;
      assert(false);
    }
#endif

    if (this_error_fixed > max_error_fixed) { max_error_fixed = this_error_fixed; }
  }

  if (max_error_fixed > cumulative_max_error_fixed) { cumulative_max_error_fixed = max_error_fixed; }

  passed = (max_error_fixed < allowed_error);

  if (passed) { printf("PASSED , max err (%f)\n", max_error_fixed); }
  else        { printf("FAILED , max err (%f)\n", max_error_fixed); } // LCOV_EXCL_LINE

  return 0;
}

// ===============================================================================
// Function: test_driver_float()
// Description: test_driver function for ac_float inputs and outputs.

template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
int test_driver_float(
  double &cumulative_max_error_float,
  const double allowed_error_float,
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

  cout << "TEST: ac_log2_cordic() INPUT: ";
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
      T_out output_float;
      test_ac_log2_cordic_float(input_float, output_float);
      double actual_float = output_float.to_double();
      double expected_float = T_out(log2(input_float.to_double())).to_double();
      double this_error_float = abs_double(expected_float - actual_float)* 100.0;

      if (this_error_float > max_error_float) { max_error_float = this_error_float; }

#ifdef DEBUG
      if (this_error_float > allowed_error_float) {
        cout << endl;
        cout << "  Error exceeds tolerance" << endl;
        cout << "  input_float      = " << input_float << endl;
        cout << "  expected_float   = " << expected_float << endl;
        cout << "  actual_float     = " << actual_float << endl;
        cout << "  this_error_float = " << this_error_float << endl;
      }
#endif
    }
  }

  bool passed = (max_error_float < allowed_error_float);

  if (passed) {
    printf("PASSED , max err (%f) \n", max_error_float);
  } else {
    printf("FAILED , max err (%f) \n", max_error_float); // LCOV_EXCL_LINE
  }

  if (max_error_float > cumulative_max_error_float) { cumulative_max_error_float = max_error_float; }

  return 0;
}

// ==============================================================================
// Function: test_driver_stfloat()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_std_float inputs. It uses the type information to iterate
//   through a range of valid values on that type in order to compare the
//   precision of the cordic log2 model with the computed log2 using
//   a standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <bool OR_TF, int Wstfl, int Estfl, int outWstfl, int outEstfl>
int test_driver_stfloat(
  double &cumulative_max_error_stfloat,
  const double allowed_error,
  bool details = false
) {
  // In order to insure correct quantization while calculating the error, the output bitwidth
  // is limited to 64 or less, while the output exponent width is limited to 11 or less.
  static_assert(outWstfl <= 64, "Output ac_std_float bitwidth must not be greater than 64.");
  static_assert(outEstfl <= 11, "Output ac_std_float exponent width must not be greater than 11.");

  double max_error_stfloat = 0.0; // reset for this run

  typedef ac_std_float<Wstfl, Estfl> T_in;
  typedef ac_std_float<outWstfl, outEstfl> T_out;

  const int W2acfl = Wstfl - Estfl + 1;
  const int I2acfl = 2;
  // ac_float equivalent of input ac_std_float type.
  typedef ac_float<W2acfl, I2acfl, Estfl> T_in_acfl;

  // Set the lower limit, upper limit and step size of the test iterations.
  ac_int<W2acfl - 2, false> sample_mantissa_slc_stfloat;
  ac_int<W2acfl - 2, false> lower_limit_stfloat = 0; // Set to zero because only +ve values are supported.
  ac_int<W2acfl - 2, false> upper_limit_stfloat = sample_mantissa_slc_stfloat.template set_val<AC_VAL_MAX>().to_double();
  ac_int<W2acfl - 2, false> step_stfloat = 0;
  // For larger ac_std_float bitwidths, the step size is increased by increasing the number of non-zero
  // LSBs, so as to reduce the time taken.
  for (int i = 0; i < AC_MAX(W2acfl - 2 - 18, 1); i++) { step_stfloat[i] = 1; }

  cout << "TEST: ac_log2_cordic() INPUT: ";
  cout.width(32);
  cout << left << T_in::type_name();
  cout << "OUTPUT: ";
  cout.width(32);
  cout << left << T_out::type_name();
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl << "  Ranges for testing iterations:" << endl; // LCOV_EXCL_LINE
    cout         << "    lower_limit_stfloat = " << lower_limit_stfloat << endl; // LCOV_EXCL_LINE
    cout         << "    upper_limit_stfloat = " << upper_limit_stfloat << endl; // LCOV_EXCL_LINE
    cout         << "    step_stfloat        = " << step_stfloat << endl; // LCOV_EXCL_LINE
    cout         << "    allowed_error       = " << allowed_error << endl; // LCOV_EXCL_LINE
  }

  // sample_exponent_array_stfloat stores all values of ac_std_float exponent to be tested.
  const int exp_arr_size_stfloat = 2*(Estfl - 1) + 3;
  ac_int<Estfl, true> sample_exponent_array_stfloat[exp_arr_size_stfloat];
  ac_int<Estfl, true> sample_exponent_value_stfloat;

  const int bias = T_in::exp_bias;
  // The first element of the array is the minimum exponent value, the middle element is a zero exponent,
  // and the last element is the maximum possible value.
  sample_exponent_array_stfloat[0]                        = -(bias - 1);
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
  
  for (int i = 0; i < exp_arr_size_stfloat; i++) {
    for (ac_int<W2acfl - 1, false> mant_i = lower_limit_stfloat; mant_i <= upper_limit_stfloat; mant_i += step_stfloat) {
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
      
      test_ac_log2_cordic_stfloat<OR_TF>(input_stfloat, output_stfloat);
      
      double expected_value_stfloat = ((T_out)log2(input_acfl_norm.to_double())).to_double();
      double actual_value_stfloat   = output_stfloat.to_double();
      double this_error_stfloat = abs(expected_value_stfloat - actual_value_stfloat)*100.0;
        
      if (this_error_stfloat > max_error_stfloat) {
        max_error_stfloat = this_error_stfloat;
      }
    }
  }

  if (max_error_stfloat > cumulative_max_error_stfloat) {
    cumulative_max_error_stfloat = max_error_stfloat;
  }

  bool passed = (max_error_stfloat < allowed_error);
  
  #ifdef AC_HCORDIC_NAN_SUPPORTED
  // +nan input should give +nan output.
  T_in input_nan = T_in::nan();
  T_out output_nan;
  test_ac_log2_cordic_stfloat<OR_TF>(input_nan, output_nan);
  bool nan_testing_passed = output_nan.isnan() && !output_nan.signbit();

  // -nan input should give -nan output.
  input_nan = -T_in::nan();
  test_ac_log2_cordic_stfloat<OR_TF>(input_nan, output_nan);
  nan_testing_passed = nan_testing_passed && output_nan.isnan() && output_nan.signbit();
  
  AC_ASSERT(nan_testing_passed, "NaN testing FAILED\n");
  #endif

  if (passed) { printf("PASSED , max err (%f) \n", max_error_stfloat); }
  else        { printf("FAILED , max err (%f) \n", max_error_stfloat); }
  
  return 0;
}

// ==============================================================================
// Function: test_driver_ifloat()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_ieee_float inputs. It uses the type information to iterate
//   through a range of valid values on that type in order to compare the
//   precision of the cordic log2 model with the computed log2 using
//   a standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <bool OR_TF, ac_ieee_float_format in_format, ac_ieee_float_format out_format>
int test_driver_ifloat(
  double &cumulative_max_error_ifloat,
  const double allowed_error,
  bool details = false
) {
  // In order to insure correct quantization while calculating the error, the output format must not
  // be binary128 or binary256.
  static_assert(out_format != binary128 && out_format != binary256, "binary128 and binary256 output formats not supported.");

  double max_error_ifloat = 0.0; // reset for this run

  typedef ac_ieee_float<in_format> T_in;
  typedef ac_ieee_float<out_format> T_out;

  const int Wifl = T_in::width;
  const int Eifl = T_in::e_width;
  const int W2acfl = Wifl - Eifl + 1;
  const int I2acfl = 2;
  // ac_float equivalent of input ac_ieee_float type.
  typedef ac_float<W2acfl, I2acfl, Eifl> T_in_acfl;

  // Set the lower limit, upper limit and step size of the test iterations.
  ac_int<W2acfl - 2, false> sample_mantissa_slc_ifloat;
  ac_int<W2acfl - 2, false> lower_limit_ifloat = 0; // Set to zero because only +ve values are supported.
  ac_int<W2acfl - 2, false> upper_limit_ifloat = sample_mantissa_slc_ifloat.template set_val<AC_VAL_MAX>().to_double();
  ac_int<W2acfl - 2, false> step_ifloat = 0;
  // For larger ac_ieee_float bitwidths, the step size is increased by increasing the number of non-zero
  // LSBs, so as to reduce the time taken.
  for (int i = 0; i < AC_MAX(W2acfl - 2 - 18, 1); i++) { step_ifloat[i] = 1; }

  cout << "TEST: ac_log2_cordic() INPUT: ";
  cout.width(32);
  cout << left << T_in::type_name();
  cout << "OUTPUT: ";
  cout.width(32);
  cout << left << T_out::type_name();
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl << "  Ranges for testing iterations:" << endl; // LCOV_EXCL_LINE
    cout         << "    lower_limit_ifloat = " << lower_limit_ifloat << endl; // LCOV_EXCL_LINE
    cout         << "    upper_limit_ifloat = " << upper_limit_ifloat << endl; // LCOV_EXCL_LINE
    cout         << "    step_ifloat        = " << step_ifloat << endl; // LCOV_EXCL_LINE
    cout         << "    allowed_error      = " << allowed_error << endl; // LCOV_EXCL_LINE
  }

  // sample_exponent_array_ifloat stores all values of ac_ieee_float exponent to be tested.
  const int exp_arr_size_ifloat = 2*(Eifl - 1) + 3;
  ac_int<Eifl, true> sample_exponent_array_ifloat[exp_arr_size_ifloat];
  ac_int<Eifl, true> sample_exponent_value_ifloat;

  const int bias = T_in::exp_bias;
  // The first element of the array is the minimum exponent value, the middle element is a zero exponent,
  // and the last element is the maximum possible value.
  sample_exponent_array_ifloat[0]                       = -(bias - 1);
  sample_exponent_array_ifloat[Eifl]                    = 0;
  sample_exponent_array_ifloat[exp_arr_size_ifloat - 1] = bias;

  // Alternate between odd and even values for the other elements. The even values are powers of two,
  // while the odd values are one less than the nearest power of two. Exponent values = +/-1 are left
  // as they are.
  bool odd_elem = true;

  for (int i = (Eifl - 2); i >= 0; i--) {
    sample_exponent_value_ifloat = 0;
    sample_exponent_value_ifloat[i] = 1;
    sample_exponent_array_ifloat[Eifl + i + 1] = sample_exponent_value_ifloat - int(odd_elem&&i!=0);
    sample_exponent_array_ifloat[Eifl - i - 1] = -(sample_exponent_value_ifloat - int(odd_elem&&i!=0));
    odd_elem = !odd_elem;
  }
  
  for (int i = 0; i < exp_arr_size_ifloat; i++) {
    for (ac_int<W2acfl - 1, false> mant_i = lower_limit_ifloat; mant_i <= upper_limit_ifloat; mant_i += step_ifloat) {
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
      
      test_ac_log2_cordic_ifloat<OR_TF>(input_ifloat, output_ifloat);
      
      double expected_value_ifloat = ((T_out)log2(input_acfl_norm.to_double())).to_ac_float().to_double();
      double actual_value_ifloat   = output_ifloat.to_ac_float().to_double();
      double this_error_ifloat = abs(expected_value_ifloat - actual_value_ifloat)*100.0;
        
      if (this_error_ifloat > max_error_ifloat) {
        max_error_ifloat = this_error_ifloat;
      }
    }
  }

  if (max_error_ifloat > cumulative_max_error_ifloat) {
    cumulative_max_error_ifloat = max_error_ifloat;
  }

  bool passed = (max_error_ifloat < allowed_error);
  
  #ifdef AC_HCORDIC_NAN_SUPPORTED
  // +nan input should give +nan output.
  T_in input_nan = T_in::nan();
  T_out output_nan;
  test_ac_log2_cordic_ifloat<OR_TF>(input_nan, output_nan);
  bool nan_testing_passed = output_nan.isnan() && !output_nan.signbit();

  // -nan input should give -nan output.
  input_nan = -T_in::nan();
  test_ac_log2_cordic_ifloat<OR_TF>(input_nan, output_nan);
  nan_testing_passed = nan_testing_passed && output_nan.isnan() && output_nan.signbit();
  
  AC_ASSERT(nan_testing_passed, "NaN testing FAILED\n");
  #endif

  if (passed) { printf("PASSED , max err (%f) \n", max_error_ifloat); }
  else        { printf("FAILED , max err (%f) \n", max_error_ifloat); }
  
  return 0;
}

int main(int argc, char *argv[])
{
  double max_error_fixed = 0.0, max_error_float = 0.0, max_error_stfloat = 0.0, max_error_ifloat = 0.0;
  double allowed_error = 0.005;

  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_log2_cordic() - Allowed error " << allowed_error << endl;

  //template <int Wfi, int Ifi, int outWfi, int outIfi, bool outSfi>
  test_driver_fixed<16,  1, 40, 12, false>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  2, 40, 12, false>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  3, 40, 12, false>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  4, 40, 12, false>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  5, 40, 12, false>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  6, 40, 12, false>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  7, 40, 12, false>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  8, 40, 12, false>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  9, 40, 12, false>(max_error_fixed, allowed_error);
  test_driver_fixed< 8, 12, 40, 12, false>(max_error_fixed, allowed_error);

  test_driver_fixed<16, -5, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed<16, -4, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed<16, -3, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed<16, -2, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed<16, -1, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  0, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  1, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  2, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  3, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  4, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  5, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  6, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  7, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  8, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed<16,  9, 40, 12,  true>(max_error_fixed, allowed_error);
  test_driver_fixed< 8, 12, 40, 12,  true>(max_error_fixed, allowed_error);

  //template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
  test_driver_float<14,  0, 7, 32, 2, 10>(max_error_float, allowed_error);
  test_driver_float<15, -3, 6, 32, 2, 10>(max_error_float, allowed_error);
  test_driver_float<15,  4, 6, 32, 2, 10>(max_error_float, allowed_error);
  test_driver_float< 8, 10, 6, 32, 2, 10>(max_error_float, allowed_error);
  test_driver_float<20,  0, 7, 32, 2, 10>(max_error_float, allowed_error);
  test_driver_float<20, -3, 6, 32, 2, 10>(max_error_float, allowed_error);
  test_driver_float<20,  4, 6, 32, 2, 10>(max_error_float, allowed_error);

  // template <bool OR_TF, int Wstfl, int Estfl, int outWstfl, int outEstfl>
  test_driver_stfloat<false, 16,  5, 32, 8>(max_error_stfloat, allowed_error);
  test_driver_stfloat<false, 32,  8, 32, 8>(max_error_stfloat, allowed_error);
  test_driver_stfloat<true,  16,  5, 64, 11>(max_error_stfloat, allowed_error);
  test_driver_stfloat<true,  32,  8, 64, 11>(max_error_stfloat, allowed_error);
  test_driver_stfloat<true,  64, 11, 64, 11>(max_error_stfloat, allowed_error);

  // template <bool OR_TF, ac_ieee_float_format in_format, ac_ieee_float_format out_format>
  test_driver_ifloat<false, binary16, binary32>(max_error_ifloat, allowed_error);
  test_driver_ifloat<false, binary32, binary32>(max_error_ifloat, allowed_error);
  test_driver_ifloat<true,  binary16, binary64>(max_error_ifloat, allowed_error);
  test_driver_ifloat<true,  binary32, binary64>(max_error_ifloat, allowed_error);
  test_driver_ifloat<true,  binary64, binary64>(max_error_ifloat, allowed_error);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum errors observed across all bit-width variations:" << endl;
  cout << "    max_error_fixed   = " << max_error_fixed << endl;
  cout << "    max_error_float   = " << max_error_float << endl;
  cout << "    max_error_stfloat = " << max_error_stfloat << endl;
  cout << "    max_error_ifloat  = " << max_error_ifloat << endl;

  // If error limits on any test value have been crossed, the test has failed
  // Notify the user that the test was a failure if that is the case.
  if ((max_error_fixed > allowed_error) || (max_error_float > allowed_error) || (max_error_stfloat > allowed_error) || (max_error_ifloat > allowed_error)) {
    cout << "  ac_log2_cordic - FAILED - Error tolerance(s) exceeded" << endl; // LCOV_EXCL_LINE
    cout << "=============================================================================" << endl; // LCOV_EXCL_LINE
    return (-1); // LCOV_EXCL_LINE
  }

  cout << "  ac_log2_cordic - PASSED" << endl;
  cout << "=============================================================================" << endl;
  return (0);
}
