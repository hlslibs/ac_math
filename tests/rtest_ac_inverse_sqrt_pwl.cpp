/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 1.0                                                 *
 *                                                                        *
 *  Release Date    : Thu Mar  8 11:17:22 PST 2018                        *
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
// ac_inverse_sqrt_pwl() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_inverse_sqrt_pwl.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench

#include <ac_math/ac_inverse_sqrt_pwl.h>
using namespace ac_math;

// ==============================================================================
// Test Design
//   This simple function allows executing the ac_inverse_sqrt_pwl() function
//   using multiple data types at the same time (in this case, ac_fixed and
//   ac_complex<ac_fixed>). Template parameters are used to configure the
//   bit-widths of the types.

template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi, bool outSfi>
void test_ac_inverse_sqrt_pwl_fixed(
  const ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP> &in1,
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> &out1,
  const ac_complex<ac_fixed<Wfi, Ifi, true, AC_TRN, AC_WRAP> > &in2,
  ac_complex<ac_fixed<outWfi, outIfi, true, AC_TRN, AC_WRAP> > &out2
)
{
  ac_inverse_sqrt_pwl(in1, out1);
  ac_inverse_sqrt_pwl(in2, out2);
}

// ==============================================================================
// Test Design
//   This simple function allows executing the ac_inverse_sqrt_pwl() function
//   using the ac_float datatype.
//   Template parameters are used to configure the bit-widths of the types.

template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
void test_ac_inverse_sqrt_pwl_float(
  const ac_float<   Wfl,    Ifl,    Efl, AC_RND>   &in3,
  ac_float<outWfl, outIfl, outEfl, AC_RND>   &out3
)
{
  ac_inverse_sqrt_pwl(in3, out3);
}

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
  // The typecasting is done in order to provide quantization on the expected output.
  double expected_value = ((T_out)(1.0 / sqrt(input.to_double()))).to_double();
  actual_value   = output.to_double();
  double this_error;

  // If expected value is greater than a particular threshold, calculate relative error, else, calculate absolute error.
  if (abs(expected_value) > threshold) {
    this_error = abs( (expected_value - actual_value) / expected_value ) * 100.0;
  } else {
    this_error = abs(expected_value - actual_value) * 100.0;
  }

  return this_error;
}

// Calculate error for complex outputs

template <class T_in, class T_out>
double cmplx_err_calc(
  const ac_complex<T_in> input,
  const ac_complex<T_out> output,
  const double allowed_error,
  const double threshold
)
{
  ac_complex<double> input_double(input.r().to_double(), input.i().to_double());
  ac_complex<double> output_sqrt, exp_op, diff_op, act_op;

  // Calculate the accurate value of the square root using double precision.
  double mod = sqrt(input_double.mag_sqr());
  double int_sqr_x = (mod + input_double.r()) / 2;
  double int_sqr_y = (mod - input_double.r()) / 2;
  output_sqrt.r() = sqrt(int_sqr_x);
  output_sqrt.i() = sqrt(int_sqr_y);
  output_sqrt.i() = input_double.i() < 0 ? -output_sqrt.i() : output_sqrt.i();
  exp_op = 1 / output_sqrt;

#ifdef DEBUG
  // Store the value before quantization. This can come in handy when debugging.
  ac_complex<double> exp_op_no_quant = exp_op;
#endif

  // Perform quantization on the expected output by converting it to the output of the expected
  // value, and then converting the quantized output back to a double.
  exp_op.r() = ((T_out)exp_op.r()).to_double();
  exp_op.i() = ((T_out)exp_op.i()).to_double();

  act_op.r() = output.r().to_double();
  act_op.i() = output.i().to_double();
  // Store the difference between the expected, accurate value, vs the actual, approximate output.
  diff_op = exp_op - act_op;

  double error;

  if (sqrt(exp_op.mag_sqr()) > threshold) {error = sqrt((diff_op / exp_op).mag_sqr()) * 100;}
  else {error = sqrt(diff_op.mag_sqr()) * 100;}

#ifdef DEBUG
  if (error > allowed_error) {
    cout << endl;
    cout << "FAILED, complex error exceeded" << endl;
    cout << "error           = " << error << endl;
    cout << "input           = " << input << endl;
    cout << "input_double    = " << input_double << endl;
    cout << "mod             = " << mod << endl;
    cout << "int_sqr_x       = " << int_sqr_x << endl;
    cout << "int_sqr_y       = " << int_sqr_y << endl;
    cout << "output_sqrt     = " << output_sqrt << endl;
    cout << "exp_op_no_quant = " << exp_op_no_quant << endl;
    cout << "exp_op          = " << exp_op << endl;
    cout << "act_op          = " << act_op << endl;
    cout << "output          = " << output << endl;
    assert(false);
  }

  // If in case the calculations for the expected value are wrong, then the double output without
  // quantization is quite likely to have "nan" as the value for the real and/or imaginary part.
  // The assert below takes care of that.
  if (isnan(float(exp_op_no_quant.r())) || isnan(float(exp_op_no_quant.i()))) {
    cout << "Real and/or imaginary parts of the calculated expected output were set to nan. Please check your calculations." << endl;
    assert(false);
  }
#endif

  // Return the error in the real part or the one in the imaginary part, whichever is greater.
  return error;

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
  // MONOTONIC: Make sure that function is monotonic. Compare old value (value of previous iteration) with current value. Since the inverse square root function we
  // are testing is a decreasing function, and our testbench value keeps incrementing or remains the same (in case of saturation), we expect the
  // old value to be greater than or equal to the current one.
  // Also, since the inverse square root function has a large discontinuity at x = 0; we make sure we don't compare values when we cross this point.
  // We do this by checking the signage of the old output vs that of the new output. Since we aren't checking for zero inputs, crossing x = 0
  // will mean that the old output is negative and the new one is positive, in case of an increasing testbench.
  bool sign_same = (old_real_output > 0 && actual_value_fixed > 0) || (old_real_output < 0 && actual_value_fixed < 0);
  if (compare && sign_same) {
    // Figuring out what the normalized value was for the input is a good way to figure out where the discontinuity occured w.r.t. the PWL segments.
    ac_fixed<Wfi, int(Sfi), Sfi, AC_TRN, AC_WRAP> norm_input_fixed;
    ac_normalize(input_fixed, norm_input_fixed);
    if (old_real_output < actual_value_fixed) {
      cout << endl;
      cout << "  Real, fixed point output not monotonic at :" << endl;
      cout << "  input_fixed = " << input_fixed << endl;
      cout << "  output_fixed = " << output_fixed << endl;
      cout << "  old_real_output = " << old_real_output << endl;
      cout << "  normalized x    = " << norm_input_fixed << endl;
      assert(false);
    }
  }
  // Update the variable for old_real_output.
  old_real_output = actual_value_fixed;
  // By setting compare to true, we make sure that once there is an old value stored, we can start comparing for monotonicity.
  compare = true;
}

// ==============================================================================
// Function: test_driver_fixed()
// Description: A templatized function that can be configured for certain bit-
//   widths of AC datatypes. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   piece-wise linear inverse_sqrt model with the computed inverse square root
//   using a standard C double type. The maximum error for each type is
//   accumulated in variables defined in the calling function.

template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi, bool outSfi>
int test_driver_fixed(
  double &cumulative_max_error_fixed,
  double &cumulative_max_error_cmplx_fixed,
  const double allowed_error_fixed,
  const double allowed_error_complex,
  const double threshold,
  bool details = false
)
{
  bool passed;
  bool check_monotonic = true;
  double max_error_fixed = 0.0; // reset for this run
  double max_error_cmplx_fixed = 0.0; // reset for this run
  double old_max_error_cmplx_fixed = 0.0; // used later to help in finding max error

  ac_fixed<   Wfi,    Ifi,    Sfi, AC_TRN, AC_WRAP>   input_fixed;
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP>   output_fixed;
  ac_complex<ac_fixed<   Wfi,    Ifi, true, AC_TRN, AC_WRAP> > cmplx_input_fixed;
  ac_complex<ac_fixed<outWfi, outIfi, true, AC_TRN, AC_WRAP> > cmplx_output_fixed;

  double lower_limit_fixed, upper_limit_fixed, step_fixed;

  // set ranges and step size for fixed point testbench
  lower_limit_fixed   = input_fixed.template set_val<AC_VAL_MIN>().to_double();
  upper_limit_fixed   = input_fixed.template set_val<AC_VAL_MAX>().to_double();
  step_fixed          = input_fixed.template set_val<AC_VAL_QUANTUM>().to_double();

  printf("TEST: ac_inverse_sqrt_pwl() INPUT: ac_fixed<%2d,%2d,%5s,%7s,%7s> OUTPUT: ac_fixed<%2d,%2d,%5s,%7s,%7s>  RESULT: ",
         Wfi,Ifi,(Sfi?"true":"false"),"AC_TRN","AC_WRAP",outWfi,outIfi,(outSfi?"true":"false"),"AC_TRN","AC_WRAP");

  // Dump the test details
  if (details) {
    cout << endl;
    cout << "  Ranges for input types:" << endl;
    cout << "    lower_limit_fixed    = " << lower_limit_fixed << endl;
    cout << "    upper_limit_fixed    = " << upper_limit_fixed << endl;
    cout << "    step_fixed           = " << step_fixed << endl;
  }

  double old_real_output;
  bool compare = false;
  double actual_value_fixed;

  // test fixed-point real and complex.

  // Fix the real part of the input at 0, and iterate through all possible values of imaginary part based on its type.
  for(double i = lower_limit_fixed; i <= upper_limit_fixed; i += step_fixed) {
    cmplx_input_fixed.r() = 0;
    cmplx_input_fixed.i() = i;
    input_fixed = i;

    if(input_fixed != 0) {
      // Pass all inputs at one go
      test_ac_inverse_sqrt_pwl_fixed(input_fixed, output_fixed, cmplx_input_fixed, cmplx_output_fixed);

      double this_error_fixed = err_calc(input_fixed, actual_value_fixed, output_fixed, allowed_error_fixed, threshold);
      double this_error_complex = cmplx_err_calc(cmplx_input_fixed, cmplx_output_fixed, allowed_error_complex, threshold);

      if(check_monotonic) { monotonicity_check(old_real_output, actual_value_fixed, compare, input_fixed, output_fixed); }
      if (this_error_fixed > max_error_fixed) {max_error_fixed = this_error_fixed;}
      if (this_error_complex > max_error_cmplx_fixed) {max_error_cmplx_fixed = this_error_complex;}
    }
  }

  // Store max_error variable in another variable for later comparisons.
  old_max_error_cmplx_fixed = max_error_cmplx_fixed;

  // Do the same thing as above, but while keeping the imaginary part fixed at 0
  for(double i = lower_limit_fixed; i <= upper_limit_fixed; i += step_fixed) {
    cmplx_input_fixed.r() = i;
    cmplx_input_fixed.i() = 0;
    input_fixed = i;

    if(input_fixed != 0) {
      // Pass all inputs at one go
      test_ac_inverse_sqrt_pwl_fixed(input_fixed, output_fixed, cmplx_input_fixed, cmplx_output_fixed);
      // All possible real fixed point values have already been tested in the first set of iterations.
      double this_error_complex = cmplx_err_calc(cmplx_input_fixed, cmplx_output_fixed, allowed_error_complex, threshold);
      // Note: there's no need to check monotonicity, because that was already taken care of in the set of iterations before this.
      if (this_error_complex > max_error_cmplx_fixed) {max_error_cmplx_fixed = this_error_complex;}
    }
  }

  // If the old value for max_error is smaller than the current one, store the current value in the variable for the old value.
  if(max_error_cmplx_fixed > old_max_error_cmplx_fixed) { old_max_error_cmplx_fixed = max_error_cmplx_fixed; }

  // Now, keep the real part at the maximum value while iterating through all the possible values for the imaginary part.
  for(double i = lower_limit_fixed; i <= upper_limit_fixed; i += step_fixed) {
    cmplx_input_fixed.r() = upper_limit_fixed;
    cmplx_input_fixed.i() = i;
    input_fixed = upper_limit_fixed;

    if(input_fixed != 0) {
      // Pass all inputs at one go
      test_ac_inverse_sqrt_pwl_fixed(input_fixed, output_fixed, cmplx_input_fixed, cmplx_output_fixed);
      double this_error_complex = cmplx_err_calc(cmplx_input_fixed, cmplx_output_fixed, allowed_error_complex, threshold);
      if (this_error_complex > max_error_cmplx_fixed) {max_error_cmplx_fixed = this_error_complex;}
    }
  }

  // If the old values for max_error were smaller than the current ones, store the current values in the variables for the old value.
  if(max_error_cmplx_fixed > old_max_error_cmplx_fixed) { old_max_error_cmplx_fixed = max_error_cmplx_fixed; }

  // If the real/imaginary part is unsigned, minimum limit is zero. Hence, we need not go through this stage of testing if the real part
  // is unsigned, because it has always been covered in the first two sets of iterations, in which we fix the real/imaginary value at zero.
  if(Sfi) {
    // Now, keep the real part at the minimum value while iterating through all the possible values for the imaginary part.
    for(double i = lower_limit_fixed; i <= upper_limit_fixed; i += step_fixed) {
      cmplx_input_fixed.r() = lower_limit_fixed;
      cmplx_input_fixed.i() = i;
      input_fixed = lower_limit_fixed;

      if(input_fixed != 0) {
        // Pass all inputs at one go
        test_ac_inverse_sqrt_pwl_fixed(input_fixed, output_fixed, cmplx_input_fixed, cmplx_output_fixed);
        double this_error_complex = cmplx_err_calc(cmplx_input_fixed, cmplx_output_fixed, allowed_error_complex, threshold);
        if (this_error_complex > max_error_cmplx_fixed) {max_error_cmplx_fixed = this_error_complex;}
      }
    }

    // If the old values for max_error were smaller than the current ones, store the current values in the variables for the old value.
    if(max_error_cmplx_fixed > old_max_error_cmplx_fixed) { old_max_error_cmplx_fixed = max_error_cmplx_fixed; }

    // Now, keep the real part at the minimum value while iterating through all the possible values for the imaginary part.
    for(double i = lower_limit_fixed; i <= upper_limit_fixed; i += step_fixed) {
      cmplx_input_fixed.r() = i;
      cmplx_input_fixed.i() = lower_limit_fixed;
      input_fixed = lower_limit_fixed;

      if(input_fixed != 0) {
      // Pass all inputs at one go
      test_ac_inverse_sqrt_pwl_fixed(input_fixed, output_fixed, cmplx_input_fixed, cmplx_output_fixed);
      double this_error_complex = cmplx_err_calc(cmplx_input_fixed, cmplx_output_fixed, allowed_error_complex, threshold);
      if (this_error_complex > max_error_cmplx_fixed) {max_error_cmplx_fixed = this_error_complex;}
      }
    }

    // If the old values for max_error were smaller than the current ones, store the current values in the variables for the old value.
    if(max_error_cmplx_fixed > old_max_error_cmplx_fixed) { old_max_error_cmplx_fixed = max_error_cmplx_fixed; }
  }

  // Now, iterate through all the possible combinations of one-hot encodings of the real/imaginary parts

  if(!Sfi) {
    // Give the input a non-zero dummy value
    input_fixed[0] = 1;

    for(int i = 0; i < Wfi; i++) {
      cmplx_input_fixed.r() = 0;
      cmplx_input_fixed.r()[i] = 1;
      for(int j = 0; j < Wfi; j++) {
        cmplx_input_fixed.r() = 0;
        cmplx_input_fixed.r()[i] = 1;

        // The inputs can never be of zero magnitude in this case. Hence, skip zero checking for inputs.
        // Pass all inputs at one go
        test_ac_inverse_sqrt_pwl_fixed(input_fixed, output_fixed, cmplx_input_fixed, cmplx_output_fixed);
        double this_error_complex = cmplx_err_calc(cmplx_input_fixed, cmplx_output_fixed, allowed_error_complex, threshold);
        if (this_error_complex > max_error_cmplx_fixed) {max_error_cmplx_fixed = this_error_complex;}
      }
    }

    // If the old values for max_error were smaller than the current ones, store the current values in the variables for the old value.
    if(max_error_cmplx_fixed > old_max_error_cmplx_fixed) { old_max_error_cmplx_fixed = max_error_cmplx_fixed; }
  } else {
    // If input is signed, two separate sets of iterations are required: one to go through the negative
    // inputs and the other to go through the positive inputs, all of which are positive/negative powers
    // of two due to one-hot encoding.
    for(int i = 0; i < Wfi - 1; i++) {
      cmplx_input_fixed.r() = 0;
      cmplx_input_fixed.r()[i] = 1;
      for(int j = 0; j < Wfi - 1; j++) {
        cmplx_input_fixed.i() = 0;
        cmplx_input_fixed.i()[j] = 1;

        // The inputs can never be of zero magnitude in this case. Hence, skip zero checking for inputs.
        // Pass all inputs at one go
        test_ac_inverse_sqrt_pwl_fixed(input_fixed, output_fixed, cmplx_input_fixed, cmplx_output_fixed);
        double this_error_complex = cmplx_err_calc(cmplx_input_fixed, cmplx_output_fixed, allowed_error_complex, threshold);
        if (this_error_complex > max_error_cmplx_fixed) {max_error_cmplx_fixed = this_error_complex;}
      }
    }

    // If the old values for max_error were smaller than the current ones, store the current values in the variables for the old value.
    if(max_error_cmplx_fixed > old_max_error_cmplx_fixed) { old_max_error_cmplx_fixed = max_error_cmplx_fixed; }

    for(int i = 0; i < Wfi - 1; i++) {
      cmplx_input_fixed.r() = 0;
      cmplx_input_fixed.r()[i] = 1;
      cmplx_input_fixed.r() = -cmplx_input_fixed.r();
      for(int j = 0; j < Wfi - 1; j++) {
        cmplx_input_fixed.i() = 0;
        cmplx_input_fixed.i()[j] = 1;
        cmplx_input_fixed.i() = -cmplx_input_fixed.i();

        // The inputs can never be of zero magnitude in this case. Hence, skip zero checking for inputs.
        // Pass all inputs at one go
        test_ac_inverse_sqrt_pwl_fixed(input_fixed, output_fixed, cmplx_input_fixed, cmplx_output_fixed);
        double this_error_complex = cmplx_err_calc(cmplx_input_fixed, cmplx_output_fixed, allowed_error_complex, threshold);
        if (this_error_complex > max_error_cmplx_fixed) {max_error_cmplx_fixed = this_error_complex;}
      }
    }

    // If the old values for max_error were smaller than the current ones, store the current values in the variables for the old value.
    if(max_error_cmplx_fixed > old_max_error_cmplx_fixed) { old_max_error_cmplx_fixed = max_error_cmplx_fixed; }
  }

  max_error_cmplx_fixed = old_max_error_cmplx_fixed;

  passed = (max_error_fixed < allowed_error_fixed) && (max_error_cmplx_fixed < allowed_error_complex);

  if (passed) { printf("PASSED , max err (%f) (%f complex)\n", max_error_fixed, max_error_cmplx_fixed); }
  else        { printf("FAILED , max err (%f) (%f complex)\n", max_error_fixed, max_error_cmplx_fixed); }

  if (max_error_fixed>cumulative_max_error_fixed) { cumulative_max_error_fixed = max_error_fixed; }
  if (max_error_cmplx_fixed>cumulative_max_error_cmplx_fixed) { cumulative_max_error_cmplx_fixed = max_error_cmplx_fixed; }

  return 0;
}

// ==============================================================================
// Function: test_driver_float()
// Description: A templatized function that can be configured for certain bit-
//   widths of AC datatypes. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   piece-wise linear inverse_sqrt model with the computed inverse square root
//   using a standard C double type. The maximum error for each type is
//   accumulated in variables defined in the calling function.

template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
int test_driver_float(
  double &cumulative_max_error_float,
  const double allowed_error_float,
  const double threshold,
  bool details = false
)
{
  bool passed = true;
  bool check_monotonic = true;
  double max_error_float = 0.0; // reset for this run

  ac_float<   Wfl,    Ifl,    Efl, AC_RND> input_float;
  ac_float<outWfl, outIfl, outEfl, AC_RND> output_float;

  // Declare an ac_fixed variable of same type as mantissa
  ac_fixed<Wfl, Ifl, true> sample_mantissa;
  double lower_limit_mantissa, upper_limit_mantissa, step_mantissa;

  lower_limit_mantissa = sample_mantissa.template set_val<AC_VAL_MIN>().to_double();
  upper_limit_mantissa = sample_mantissa.template set_val<AC_VAL_MAX>().to_double();
  step_mantissa        = sample_mantissa.template set_val<AC_VAL_QUANTUM>().to_double();

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

  printf("TEST: ac_inverse_sqrt_pwl()  INPUT: ac_float<%2d,%2d,%2d,%7s>  OUTPUT: ac_float<%2d,%2d,%2d,%7s>    RESULT: ",
         Wfl,Ifl,Efl,"AC_RND",outWfl,outIfl,outEfl,"AC_RND");

  // Dump the test details
  if (details) {
    cout << endl << "  Ranges for input types:" << endl;
    cout         << "    lower_limit_mantissa = " << lower_limit_mantissa << endl;
    cout         << "    upper_limit_mantissa = " << upper_limit_mantissa << endl;
    cout         << "    step_mantissa        = " << step_mantissa << endl;
    cout         << "    allowed_error_float  = " << allowed_error_float << endl;
  }

  double old_real_output;

  for (int i = 0; i < exp_arr_size; i++) {
    // Extract a value to be tested for the exponent part.
    input_float.e = sample_exponent_array[i];
    bool compare = false;

    // For that particular exponent value, go through every possible value that can be represented by the mantissa.
    for (double mant_i = 0; mant_i <= upper_limit_mantissa; mant_i += step_mantissa) {
      input_float.m = mant_i;
      if (input_float != 0) {
        test_ac_inverse_sqrt_pwl_float(input_float, output_float);
        double expected_value_float   = 1.0/ sqrt(input_float.to_double());
        double actual_value_float     = output_float.to_double();

        double this_error_float;

        // If expected value is greater than a particular threshold, calculate relative error, else, calculate absolute error.
        if (abs(expected_value_float) > threshold) {this_error_float = abs( (expected_value_float - actual_value_float) / expected_value_float ) * 100.0;}
        else {this_error_float = abs(expected_value_float - actual_value_float) * 100.0;}

        if (check_monotonic) {
          // This function has the same basic working as the fixed point version, the only difference being is that now, ac_float values are considered
          bool sign_same = (old_real_output > 0 && actual_value_float > 0) || (old_real_output < 0 && actual_value_float < 0);
          if (compare && sign_same) {
            // Figure out what the input mantissa was normalized to
            ac_fixed<input_float.m.width, int(input_float.m.sign), input_float.m.sign, input_float.m.q_mode> norm_input_mant_fixed;
            ac_normalize(input_float.m, norm_input_mant_fixed);
            if (old_real_output < actual_value_float) {
              cout << endl;
              cout << "  Real, floating point output not monotonic at :" << endl;
              cout << "x = " << input_float << endl;
              cout << "y = " << output_float << endl;
              cout << "old_real_output = " << old_real_output << endl;
              cout << "normalized mantissa of x = " << norm_input_mant_fixed << endl;
              assert(false);
            }
          }
          // Update the variable for old_real_output.
          old_real_output = actual_value_float;
          // By setting compare to true, we make sure that once there is an old value stored, we can start comparing for monotonicity.
          compare = true;
        }

        if (this_error_float > max_error_float) {max_error_float = this_error_float;}
      }
    }
  }

  passed = (max_error_float < allowed_error_float);

  if (passed) { printf("PASSED , max error (%f)\n", max_error_float); }
  else        { printf("FAILED , max error (%f)\n", max_error_float); }

  if (max_error_float > cumulative_max_error_float) { cumulative_max_error_float = max_error_float; }

  return 0;

}

int main(int argc, char *argv[])
{
  double max_error_fixed = 0, cmplx_max_error_fixed = 0, max_error_float = 0;

  // Set tolerance
  double allowed_error_fixed = 0.5;
  // ac_complex values for input and output often give a higher error than ac_fixed values,
  // owing to certain intermediate computations in the inverse_sqrt function.
  // Hence, it is best to define a separate error tolerance for them.
  double allowed_error_complex = 3;
  // threshold below which we calculate absolute error instead of relative for fixed point
  double threshold_fixed = 0.005;

  // Set tolerance
  double allowed_error_float = 0.5;
  // threshold below which we calculate absolute error instead of relative for floating point
  double threshold_float = 0.005;

  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_inverse_sqrt_pwl() - Allowed error " << allowed_error_fixed << " (fixed pt), " << allowed_error_complex << " (complex fixed pt), " << allowed_error_float << " (float pt)" << endl;

  // template <int Wfi, int Ifi, int Sfi, int outWfi, int outIfi, intoutSfi>
  test_driver_fixed< 12,  0, false, 64, 32, false>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, allowed_error_complex, threshold_fixed);
  test_driver_fixed<  8, -2, false, 64, 32, false>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, allowed_error_complex, threshold_fixed);
  test_driver_fixed< 10, -4, false, 64, 32, false>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, allowed_error_complex, threshold_fixed);
  test_driver_fixed< 12,  8, false, 64, 32, false>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, allowed_error_complex, threshold_fixed);
  test_driver_fixed< 16,  8, false, 64, 32, false>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, allowed_error_complex, threshold_fixed);
  test_driver_fixed<  9,  4, false, 60, 30, false>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, allowed_error_complex, threshold_fixed);
  test_driver_fixed<  4,  9, false, 60, 30, false>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, allowed_error_complex, threshold_fixed);
  test_driver_fixed< 12, 12, false, 64, 32, false>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, allowed_error_complex, threshold_fixed);
  test_driver_fixed< 12,  5, false, 64, 32, false>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed, allowed_error_complex, threshold_fixed);

  // template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
  test_driver_float<  5,  3, 3, 64, 32, 10>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<  5,  1, 3, 64, 32, 10>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<  5,  0, 3, 64, 32, 10>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<  5, -2, 3, 64, 32, 10>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<  5,  9, 3, 60, 30, 11>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<  5,  5, 3, 64, 32, 10>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<  5,  5, 1, 64, 32, 10>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float< 10,  5, 4, 61, 33, 11>(max_error_float, allowed_error_float, threshold_float);
  test_driver_float<  5,  3, 5, 64, 32, 10>(max_error_float, allowed_error_float, threshold_float);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum errors observed across all data type / bit-width variations:" << endl;
  cout << "    max_error_fixed       = " << max_error_fixed << endl;
  cout << "    cmplx_max_error_fixed = " << cmplx_max_error_fixed << endl;
  cout << "    max_error_float       = " << max_error_float << endl;

  // If error limits on any tested datatype have been crossed, the test has failed
  bool test_fail = (max_error_fixed > allowed_error_fixed) || (cmplx_max_error_fixed > allowed_error_complex) || (max_error_float > allowed_error_float);

  // Notify the user that the test was a failure.
  if (test_fail) {
    cout << "  ac_inverse_sqrt_pwl - FAILED - Error tolerance(s) exceeded" << endl;
    cout << "=============================================================================" << endl;
    return -1;
  } else {
    cout << "  ac_inverse_sqrt_pwl - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }
  return 0;
}























