/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 3.6                                                 *
 *                                                                        *
 *  Release Date    : Tue Nov 12 23:14:00 PST 2024                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.6.0                                               *
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
// ac_exp2_cordic() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_exp2_cordic.cpp -o design
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
//   This simple function allows executing the ac_exp2_cordic() function.
//   Template parameters are used to configure the bit-widths of the
//   ac_fixed inputs.

template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi>
void test_ac_exp2_cordic_fixed(
  const ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP> &in,
  ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> &out
)
{
  ac_exp2_cordic(in, out);
}

// Test for ac_float input and output.
template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
void test_ac_exp2_cordic_float(
  const ac_float<Wfl, Ifl, Efl, AC_TRN> &in,
  ac_float<outWfl, outIfl, outEfl, AC_TRN> &out
)
{
  ac_exp2_cordic(in, out);
}

// Test for ac_std_float input and output.
template <bool OR_TF, int Wstfl, int Estfl, int outWstfl, int outEstfl>
void test_ac_exp2_cordic_stfloat(
  const ac_std_float<Wstfl, Estfl> &in,
  ac_std_float<outWstfl, outEstfl> &out
)
{
  ac_exp2_cordic<OR_TF>(in, out);
}

// Test for ac_ieee_float input and output.
template <bool OR_TF, ac_ieee_float_format in_format, ac_ieee_float_format out_format>
void test_ac_exp2_cordic_ifloat(
  const ac_ieee_float<in_format> &in,
  ac_ieee_float<out_format> &out
)
{
  ac_exp2_cordic<OR_TF>(in, out);
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
//   cordic table exp2 model with the computed power of two using a
//   standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi>
int test_driver_fixed(
  double &cumulative_max_error_fixed,
  const double allowed_error,
  const double threshold,
  bool details = false
)
{
  bool passed = true;
  double max_error_fixed = 0.0; // reset for this run

  ac_fixed<   Wfi,    Ifi,   Sfi, AC_TRN, AC_WRAP> input;
  typedef ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> T_out;
  T_out output_fixed;

  double lower_limit, upper_limit, step;

  // set ranges and step size for testbench
  lower_limit = input.template set_val<AC_VAL_MIN>().to_double();
  upper_limit = input.template set_val<AC_VAL_MAX>().to_double();
  step        = input.template set_val<AC_VAL_QUANTUM>().to_double();

  cout << "TEST: ac_exp2_cordic() INPUT: ";
  cout.width(38);
  cout << left << input.type_name();
  cout << "OUTPUT: ";
  cout.width(38);
  cout << left << output_fixed.type_name();
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl; // LCOV_EXCL_LINE
    cout << "  Ranges for input types:" << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit          = " << lower_limit << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit          = " << upper_limit << endl; // LCOV_EXCL_LINE
    cout << "    step                 = " << step << endl; // LCOV_EXCL_LINE
  }

  for (double i = lower_limit; i <= upper_limit; i += step) {
    // Set values for input.
    input = i;
    test_ac_exp2_cordic_fixed(input, output_fixed);

    // call reference exp2() with fixed-pt value converted back to double
    // an additional step of typecasting is required in order to perform
    // quantization on the expected output.
    double expected_value_fixed = ((T_out)exp2(input.to_double())).to_double();
    double actual_value_fixed   = output_fixed.to_double();

    double this_error_fixed;
    // If expected value of output falls below the threshold, calculate absolute error instead of relative
    if (expected_value_fixed > threshold) {this_error_fixed = abs_double( (expected_value_fixed - actual_value_fixed) / expected_value_fixed ) * 100.0;}
    else {this_error_fixed = abs_double(expected_value_fixed - actual_value_fixed) * 100.0;}

#ifdef DEBUG
    if (this_error_fixed > allowed_error) {
      cout << endl;
      cout << "  Error exceeds tolerance" << endl;
      cout << "  input                = " << input << endl;
      cout << "  expected_value_fixed = " << expected_value_fixed << endl;
      cout << "  actual_value_fixed   = " << actual_value_fixed << endl;
      cout << "  this_error_fixed     = " << this_error_fixed << endl;
      cout << "  threshold            = " << threshold << endl;
      assert(false);
    }
#endif

    if (this_error_fixed > max_error_fixed) {max_error_fixed = this_error_fixed;}
  }

  passed = (max_error_fixed < allowed_error);

  if (passed) { printf("PASSED , max err %f\n", max_error_fixed); }
  else        { printf("FAILED , max err %f\n", max_error_fixed); } // LCOV_EXCL_LINE

  if (max_error_fixed>cumulative_max_error_fixed) { cumulative_max_error_fixed = max_error_fixed; }

  return 0;
}

// ==============================================================================
// Function: test_driver_float()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_float inputs. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   cordic exp2 model with the computed power using a standard C double type.
//   The maximum error for each type is accumulated in variables defined in the
//   calling function.

template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
int test_driver_float(
  double &cumulative_max_error_float,
  const double allowed_error,
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

  cout << "TEST: ac_exp2_cordic() INPUT: ";
  cout.width(32);
  cout << left << T_in::type_name();
  cout << "OUTPUT: ";
  cout.width(32);
  cout << left << T_out::type_name();
  cout << "RESULT: ";

  // sample_exponent_array stores all values of exponent to be tested.
  const int exp_arr_size = 2*(Efl - 1) + 3;
  ac_int<Efl, true> sample_exponent_array[exp_arr_size];
  ac_int<Efl, true> sample_exponent_value;

  // The first element of the array is the minimum exponent value, the middle element is a zero exponent,
  // and the last element is the maximum possible value.
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
    cout         << "    lower_limit_it = " << lower_limit_it << endl; // LCOV_EXCL_LINE
    cout         << "    upper_limit_it = " << upper_limit_it << endl; // LCOV_EXCL_LINE
    cout         << "    step_it        = " << step_it << endl; // LCOV_EXCL_LINE
    cout         << "    allowed_error  = " << allowed_error << endl; // LCOV_EXCL_LINE
  }
  
  ac_float<outWfl, outIfl, outEfl> output_max;
  output_max.template set_val<AC_VAL_MAX>();
  // input_max is the max input that can be passed to the function without causing the output to saturate.
  const double input_max = log2(output_max.to_double());
  
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
      // Use a parameterized ac_float constructor to set the mantissa and exponent of the floating point input.
      T_in input_float(input_mant, sample_exponent_array[i]);
      // Make sure that input_mant was normalized and that the mantissa and exponent values haven't changed after calling the constructor.
      if (input_float.mantissa() != input_mant || input_float.exp() != sample_exponent_array[i]) {
        cout << "input_mant was not normalized correctly." << endl;
        assert(false);
      }
      T_out output_float;
      double expected_value_float, actual_value_float, this_error_float;
      
      if (input_float < input_max) {
        test_ac_exp2_cordic_float(input_float, output_float);

        expected_value_float = ((T_out)exp2(input_float.to_double())).to_double();
        actual_value_float   = output_float.to_double();

        // If expected output falls below the threshold, calculate absolute error instead of relative
        if (expected_value_float > threshold) {
          this_error_float = abs((expected_value_float - actual_value_float)/expected_value_float)*100.0;
        } else {
          this_error_float = abs(expected_value_float - actual_value_float) * 100.0;
        }
        if (this_error_float > max_error_float) { max_error_float = this_error_float; }
      }
      
      // Repeat the same, but this time with a negative input.
      ac_fixed<Wfl, Ifl, true> input_mant_neg = -input_mant;
      T_in input_float_neg(input_mant_neg, sample_exponent_array[i]);
      if ((input_float_neg.mantissa() != input_mant_neg || input_float_neg.exp() != sample_exponent_array[i]) && mant_i != 0) {
        cout << "input_mant_neg was not normalized correctly." << endl;
        assert(false);
      }
      test_ac_exp2_cordic_float(input_float_neg, output_float);

      expected_value_float = ((T_out)exp2(input_float_neg.to_double())).to_double();
      actual_value_float   = output_float.to_double();

      if (expected_value_float > threshold) {
        this_error_float = abs((expected_value_float - actual_value_float)/expected_value_float)*100.0;
      } else {
        this_error_float = abs(expected_value_float - actual_value_float) * 100.0;
      }
      if (this_error_float > max_error_float) { max_error_float = this_error_float; }
    }
  }

  bool passed = (max_error_float < allowed_error);

  if (passed) { printf("PASSED , max err (%f)\n", max_error_float); }
  else        { printf("FAILED , max err (%f)\n", max_error_float); } // LCOV_EXCL_LINE

  if (max_error_float > cumulative_max_error_float) {
    cumulative_max_error_float = max_error_float;
  }
  
  return 0;
}

// ==============================================================================
// Function: test_driver_stfloat()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_std_float inputs. It uses the type information to iterate
//   through a range of valid values on that type in order to compare the
//   precision of the cordic exp2 model with the computed exp2 using
//   a standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <bool OR_TF, int Wstfl, int Estfl, int outWstfl, int outEstfl>
int test_driver_stfloat(
  double &cumulative_max_error_stfloat,
  const double allowed_error,
  const double threshold,
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

  cout << "TEST: ac_exp2_cordic() INPUT: ";
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
  
  T_out output_max;
  output_max = T_out::max();
  // input_max is the max input that can be passed to the function without causing the output to saturate.
  const T_in input_max(log2(output_max.to_double()));
  
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
      double expected_value_stfloat, actual_value_stfloat, this_error_stfloat;
      
      if (input_stfloat < input_max) {
        test_ac_exp2_cordic_stfloat<OR_TF>(input_stfloat, output_stfloat);
      
        expected_value_stfloat = ((T_out)exp2(input_acfl_norm.to_double())).to_double();
        actual_value_stfloat   = output_stfloat.to_double();

        // If expected output falls below the threshold, calculate absolute error instead of relative
        if (expected_value_stfloat > threshold) {
          this_error_stfloat = abs((expected_value_stfloat - actual_value_stfloat)/expected_value_stfloat)*100.0;
        } else {
          this_error_stfloat = abs(expected_value_stfloat - actual_value_stfloat)*100.0;
        }
        if (this_error_stfloat > max_error_stfloat) {
          max_error_stfloat = this_error_stfloat;
        }
      }
      
      // Repeat the same, but this time with a negative input.
      ac_fixed<W2acfl, I2acfl, true> input_mant_neg = -input_mant;
      T_in_acfl input_acfl_norm_neg(input_mant_neg, sample_exponent_array_stfloat[i]);
      if ((input_acfl_norm_neg.mantissa() != input_mant_neg || input_acfl_norm_neg.exp() != sample_exponent_array_stfloat[i]) && mant_i != 0) {
        cout << "Input value was not normalized correctly." << endl;
        assert(false);
      }
      T_in input_stfloat_neg(input_acfl_norm_neg);
      test_ac_exp2_cordic_stfloat<OR_TF>(input_stfloat_neg, output_stfloat);
      
      expected_value_stfloat = ((T_out)exp2(input_acfl_norm_neg.to_double())).to_double();
      actual_value_stfloat   = output_stfloat.to_double();

      if (expected_value_stfloat > threshold) {
        this_error_stfloat = abs((expected_value_stfloat - actual_value_stfloat)/expected_value_stfloat)*100.0;
      } else {
        this_error_stfloat = abs(expected_value_stfloat - actual_value_stfloat)*100.0;
      }
      if (this_error_stfloat > max_error_stfloat) {
        max_error_stfloat = this_error_stfloat;
      }
    }
  }

  if (max_error_stfloat > cumulative_max_error_stfloat) { cumulative_max_error_stfloat = max_error_stfloat; }

  bool passed = (max_error_stfloat < allowed_error);
  
  #ifdef AC_HCORDIC_NAN_SUPPORTED
  // +nan input should give +nan output.
  T_in input_nan = T_in::nan();
  T_out output_nan;
  test_ac_exp2_cordic_stfloat<OR_TF>(input_nan, output_nan);
  bool nan_testing_passed = output_nan.isnan() && !output_nan.signbit();

  // -nan input should give -nan output.
  input_nan = -T_in::nan();
  test_ac_exp2_cordic_stfloat<OR_TF>(input_nan, output_nan);
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
//   precision of the cordic exp2 model with the computed exp2 using
//   a standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <bool OR_TF, ac_ieee_float_format in_format, ac_ieee_float_format out_format>
int test_driver_ifloat(
  double &cumulative_max_error_ifloat,
  const double allowed_error,
  const double threshold,
  bool details = false
) {
  typedef ac_ieee_float<in_format> T_in;
  typedef ac_ieee_float<out_format> T_out;
  
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

  cout << "TEST: ac_exp2_cordic() INPUT: ";
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
  
  T_out output_max;
  output_max = T_out::max();
  // input_max is the max input that can be passed to the function without causing the output to saturate.
  const T_in input_max(log2(output_max.to_ac_float().to_double()));
  
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
      double expected_value_ifloat, actual_value_ifloat, this_error_ifloat;
      
      if (input_ifloat < input_max) {
        T_out output_ifloat;
        test_ac_exp2_cordic_ifloat<OR_TF>(input_ifloat, output_ifloat);
      
        expected_value_ifloat = ((T_out)exp2(input_acfl_norm.to_double())).to_ac_float().to_double();
        actual_value_ifloat   = output_ifloat.to_ac_float().to_double();

        // If expected output falls below the threshold, calculate absolute error instead of relative
        if (expected_value_ifloat > threshold) {
          this_error_ifloat = abs((expected_value_ifloat - actual_value_ifloat)/expected_value_ifloat)*100.0;
        } else {
          this_error_ifloat = abs(expected_value_ifloat - actual_value_ifloat)*100.0;
        }
        if (this_error_ifloat > max_error_ifloat) {
          max_error_ifloat = this_error_ifloat;
        }
      }
      
      // Repeat the same, but this time with a negative input.
      ac_fixed<W2acfl, I2acfl, true> input_mant_neg = -input_mant;
      T_in_acfl input_acfl_norm_neg(input_mant_neg, sample_exponent_array_ifloat[i]);
      if ((input_acfl_norm_neg.mantissa() != input_mant_neg || input_acfl_norm_neg.exp() != sample_exponent_array_ifloat[i]) && mant_i != 0) {
        cout << "Input value was not normalized correctly." << endl;
        assert(false);
      }
      T_in input_ifloat_neg(input_acfl_norm_neg);
      test_ac_exp2_cordic_ifloat<OR_TF>(input_ifloat_neg, output_ifloat);
      
      expected_value_ifloat = ((T_out)exp2(input_acfl_norm_neg.to_double())).to_ac_float().to_double();
      actual_value_ifloat   = output_ifloat.to_ac_float().to_double();

      if (expected_value_ifloat > threshold) {
        this_error_ifloat = abs((expected_value_ifloat - actual_value_ifloat)/expected_value_ifloat)*100.0;
      } else {
        this_error_ifloat = abs(expected_value_ifloat - actual_value_ifloat)*100.0;
      }
      if (this_error_ifloat > max_error_ifloat) {
        max_error_ifloat = this_error_ifloat;
      }
    }
  }

  if (max_error_ifloat > cumulative_max_error_ifloat) { cumulative_max_error_ifloat = max_error_ifloat; }

  bool passed = (max_error_ifloat < allowed_error);
  
  #ifdef AC_HCORDIC_NAN_SUPPORTED
  // +nan input should give +nan output.
  T_in input_nan = T_in::nan();
  T_out output_nan;
  test_ac_exp2_cordic_ifloat<OR_TF>(input_nan, output_nan);
  bool nan_testing_passed = output_nan.isnan() && !output_nan.signbit();

  // -nan input should give -nan output.
  input_nan = -T_in::nan();
  test_ac_exp2_cordic_ifloat<OR_TF>(input_nan, output_nan);
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
  double allowed_error = 0.5;
  double threshold = 0.005;

  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_exp2_cordic() - Allowed error " << allowed_error << endl;

  // template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi>
  test_driver_fixed< 16,  4, false, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16,  3, false, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16,  2, false, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16,  1, false, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16,  0, false, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16, -1, false, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16, -2, false, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16, -3, false, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16, -4, false, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<  2,  4, false, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<  3,  4, false, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<  1,  3, false, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<  5,  4, false, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16,  5,  true, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16,  4,  true, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16,  3,  true, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16,  2,  true, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16,  1,  true, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16,  0,  true, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16, -1,  true, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16, -2,  true, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16, -3,  true, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 16, -4,  true, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<  2,  4,  true, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<  3,  4,  true, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<  3,  3,  true, 50, 30>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<  5,  4,  true, 50, 30>(max_error_fixed, allowed_error, threshold);

  // template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
  test_driver_float<16,  4, 6, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<16,  5, 6, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<16, -3, 6, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<16, -2, 6, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<16,  4, 8, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<16,  5, 8, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<16, -3, 8, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<16, -2, 8, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<20,  4, 6, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<20,  5, 6, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<20, -3, 6, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<20, -2, 6, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<20,  4, 8, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<20,  5, 8, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<20, -3, 8, 32, 5, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<20, -2, 8, 32, 5, 10>(max_error_float, allowed_error, threshold);

  // template <bool OR_TF, int Wstfl, int Estfl, int outWstfl, int outEstfl>
  test_driver_stfloat<false, 16,  5, 32, 8>(max_error_stfloat, allowed_error, threshold);
  test_driver_stfloat<false, 32,  8, 32, 8>(max_error_stfloat, allowed_error, threshold);
  test_driver_stfloat<true,  16,  5, 64, 11>(max_error_stfloat, allowed_error, threshold);
  test_driver_stfloat<true,  32,  8, 64, 11>(max_error_stfloat, allowed_error, threshold);
  test_driver_stfloat<true,  64, 11, 64, 11>(max_error_stfloat, allowed_error, threshold);

  // template <bool OR_TF, ac_ieee_float_format in_format, ac_ieee_float_format out_format>
  test_driver_ifloat<false, binary16, binary32>(max_error_ifloat, allowed_error, threshold);
  test_driver_ifloat<false, binary32, binary32>(max_error_ifloat, allowed_error, threshold);
  test_driver_ifloat<true,  binary16, binary64>(max_error_ifloat, allowed_error, threshold);
  test_driver_ifloat<true,  binary32, binary64>(max_error_ifloat, allowed_error, threshold);
  test_driver_ifloat<true,  binary64, binary64>(max_error_ifloat, allowed_error, threshold);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum error observed across all bit-width variations:" << endl;
  cout << "    max_error_fixed   = " << max_error_fixed << endl;
  cout << "    max_error_float   = " << max_error_float << endl;
  cout << "    max_error_stfloat = " << max_error_stfloat << endl;
  cout << "    max_error_ifloat  = " << max_error_ifloat << endl;

  bool test_fail = (max_error_fixed > allowed_error) || (max_error_float > allowed_error) || (max_error_stfloat > allowed_error) || (max_error_ifloat > allowed_error);

  // If error limits on any test value have been crossed, the test has failed
  // Notify the user that the test was a failure if that is the case.
  if (test_fail) {
    cout << "  ac_exp2_cordic - FAILED - Error tolerance(s) exceeded" << endl; // LCOV_EXCL_LINE
    cout << "=============================================================================" << endl; // LCOV_EXCL_LINE
    return -1; // LCOV_EXCL_LINE
  } else {
    cout << "  ac_exp2_cordic - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }

  return 0;
}
