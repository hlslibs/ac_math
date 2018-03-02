/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 1.0                                                 *
 *                                                                        *
 *  Release Date    : Thu Mar  1 16:35:45 PST 2018                        *
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
// ac_sqrt_pwl() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_sqrt_pwl.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_sqrt_pwl.h>
using namespace ac_math;

//==============================================================================
// Test Design
//   This simple function allows executing the ac_sqrt_pwl() function
//   using multiple data types at the same time. Template parameters are
//   used to configure the bit-widths of the types.

template <int Wfi, int Ifi, int outWfi, int outIfi>
void test_ac_sqrt_pwl_fixed(
  const  ac_fixed<Wfi, Ifi, false, AC_TRN, AC_WRAP>   &in1,
  ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP>   &out1,
  const  ac_complex<ac_fixed<Wfi, Ifi, true, AC_TRN, AC_WRAP> > &in2,
  ac_complex<ac_fixed<outWfi, outIfi, true, AC_TRN, AC_WRAP> > &out2
)
{
  ac_sqrt_pwl(in1, out1);
  ac_sqrt_pwl(in2, out2);
}

//==============================================================================
// Test Design
//   This simple function allows executing the ac_sqrt_pwl() function
//   using the float datatype. Template parameters are
//   used to configure the bit-widths of the types.

template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
void test_ac_sqrt_pwl_float(
  const ac_float<Wfl, Ifl, Efl, AC_TRN>     &in3,
  ac_float<outWfl, outIfl, outEfl, AC_TRN>  &out3
)
{
  ac_sqrt_pwl(in3, out3);
}

//==============================================================================

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

//------------------------------------------------------------------------------
// Helper functions for working with AC_COMPLEX
// Overloaded function to convert test input data (double) into specific type
template<class T>
void double_to_complex(
  const double double_value,
  ac_complex<T> &type_value)
{
  type_value.r() = double_value;
  type_value.i() = double_value;
}

template <class T_in, class T_out>
double cmplx_err_calc(
  const ac_complex<T_in> input,
  const ac_complex<T_out> output,
  const double allowed_error
)
{
  double in_r = input.r().to_double();
  double in_i = input.i().to_double();
  double mod = in_r*in_r + in_i*in_i;

  //Declare variables to store the expected output value, store the difference between
  //expected and actual output values and store the actual output value (converted to
  //double)
  ac_complex<double> exp_op, diff_op, act_op;

  //Convert actual output to double and store it in a separate complex variable.
  act_op.r() = output.r().to_double();
  act_op.i() = output.i().to_double();

  //Calculate the expected value using the C++ math library
  exp_op.r() = sqrt ((sqrt(mod) + in_r)/2);
  exp_op.i() = sqrt ((abs (-sqrt(mod) + in_r))/2);
  double y_neg = -exp_op.i();
  //handling sign of square root of complex number
  exp_op.i() = (in_i) < 0? y_neg : exp_op.i();

  //Perform quantization on the expected output by converting it to the output of the expected
  //value, and then converting the quantized output back to a double.
  exp_op.r() = ((T_out)exp_op.r()).to_double();
  exp_op.i() = ((T_out)exp_op.i()).to_double();

  diff_op = exp_op - act_op;
  double error;

  if (sqrt(exp_op.mag_sqr()) != 0) {error = sqrt((diff_op / exp_op).mag_sqr()) * 100;}
  else {error = sqrt(diff_op.mag_sqr()) * 100;}

  return error;
}

//==============================================================================
// Function: test_driver_fixed()
// Description: A templatized function that can be configured for certain bit-
//   widths of AC datatypes. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   piece-wise linear square root model with the computed square root using a
//   standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <int Wfi, int Ifi, int outWfi, int outIfi>
int test_driver_fixed(
  double &cumulative_max_error_fixed,
  double &cumulative_max_error_cmplx_fixed,
  const double allowed_error_fixed,
  bool details = false
)
{
  bool passed = true;
  bool check_monotonic = true;
  double max_error_fixed = 0.0; // reset for this run
  double max_error_cmplx_fixed = 0.0; // reset for this run

  ac_fixed<   Wfi,    Ifi,  false, AC_TRN, AC_WRAP>   input_fixed;
  ac_fixed<outWfi, outIfi,  false, AC_TRN, AC_WRAP>   output_fixed;
  ac_complex<ac_fixed<   Wfi,    Ifi,    true, AC_TRN, AC_WRAP> >    cmplx_input_fixed;
  ac_complex<ac_fixed<outWfi, outIfi, true, AC_TRN, AC_WRAP> >      cmplx_output_fixed;

  typedef ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP>   T_out_fixed;

  double lower_limit_fixed, upper_limit_fixed, step_fixed;

  // set ranges and step size for fixed point testbench
  lower_limit_fixed   = input_fixed.template set_val<AC_VAL_MIN>().to_double();
  upper_limit_fixed   = input_fixed.template set_val<AC_VAL_MAX>().to_double();
  step_fixed          = input_fixed.template set_val<AC_VAL_QUANTUM>().to_double();

  printf("TEST: ac_sqrt_pwl()  INPUT: ac_fixed<%2d,%2d,%5s,%7s,%7s>  OUTPUT: ac_fixed<%2d,%2d,%5s,%7s,%7s>   RESULT: ",
         Wfi,Ifi,"false","AC_TRN","AC_WRAP",outWfi,outIfi,"false","AC_TRN","AC_WRAP");

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

  // test fixed-point
  for (double i = lower_limit_fixed; i < upper_limit_fixed; i += step_fixed) {
    // Set values for real and complex fixed point inputs.
    input_fixed = i;
    double_to_complex(i, cmplx_input_fixed);

    test_ac_sqrt_pwl_fixed(input_fixed, output_fixed, cmplx_input_fixed, cmplx_output_fixed);

    double expected_value_fixed   = ((T_out_fixed)(sqrt(input_fixed.to_double()))).to_double();
    double actual_value_fixed     = output_fixed.to_double();
    double this_error_fixed;

    // If expected value is not zero, calculate relative error, else, calculate absolute error.
    if (abs(expected_value_fixed) != 0) {
      this_error_fixed = abs( (expected_value_fixed - actual_value_fixed) / expected_value_fixed ) * 100.0;
    } else {
      this_error_fixed = abs(expected_value_fixed - actual_value_fixed) * 100.0;
    }
    double this_error_complex = cmplx_err_calc(cmplx_input_fixed, cmplx_output_fixed, allowed_error_fixed);

    if (check_monotonic) {
      // MONOTONIC: Make sure that function is monotonic. Compare old value (value of previous iteration) with current value. Since the sqrt function we
      // are testing is a increasing function, and our testbench value keeps incrementing or remains the same (in case of saturation), we expect the
      // old value to be lesser than or equal to the current one.

      if (compare == true && abs(old_real_output) > abs(actual_value_fixed)) {
        cout << "  Real, fixed point output not monotonic at :" << endl;
        assert(false);
      }

      // Update the variable for old_real_output.
      old_real_output = actual_value_fixed;
      // By setting compare to true, we make sure that once there is an old value stored, we can start comparing for monotonicity.
      compare = true;
    }

    if (this_error_fixed > max_error_fixed) {max_error_fixed = this_error_fixed;}
    if (this_error_complex > max_error_cmplx_fixed) {max_error_cmplx_fixed = this_error_complex;}

  }
  if (passed) { printf("PASSED , max error (%f) (%f complex)\n", max_error_fixed, max_error_cmplx_fixed); }
  else        { printf("FAILED , max error (%f) (%f complex)\n", max_error_fixed, max_error_cmplx_fixed); }

  if (max_error_fixed>cumulative_max_error_fixed) { cumulative_max_error_fixed = max_error_fixed; }
  if (max_error_cmplx_fixed>cumulative_max_error_cmplx_fixed) { cumulative_max_error_cmplx_fixed = max_error_cmplx_fixed; }

  return 0;
}

//==============================================================================
// Function: test_driver_float()
// Description: A templatized function that can be configured for certain bit-
//   widths of AC datatypes. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   piece-wise linear reciprocal model with the computed reciprocal using a
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
  bool check_monotonic = true;
  double max_error_float = 0.0; // reset for this run

  ac_float<   Wfl,    Ifl,    Efl, AC_TRN>   input_float;
  ac_float<outWfl, outIfl, outEfl, AC_TRN>   output_float;
  typedef ac_float<outWfl, outIfl, outEfl, AC_TRN> T_out;

  double lower_limit_mantissa, upper_limit_mantissa, step_mantissa;
  double old_real_output;
  bool compare = false;

  // Declare an ac_fixed variable of same type as mantissa
  ac_fixed<Wfl, Ifl, true> sample_mantissa;

  lower_limit_mantissa = 0;
  upper_limit_mantissa = sample_mantissa.template set_val<AC_VAL_MAX>().to_double();
  step_mantissa        = sample_mantissa.template set_val<AC_VAL_QUANTUM>().to_double();

  printf("TEST: ac_sqrt_pwl()  INPUT: ac_float<%2d,%2d,%2d,%7s>  OUTPUT: ac_float<%2d,%2d,%2d,%7s>   RESULT: ",
         Wfl,Ifl,Efl,"AC_RND",outWfl,outIfl,outEfl,"AC_RND");

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

  // Dump the test details
  if (details) {
    cout << endl << "  Ranges for input types:" << endl;
    cout         << "    lower_limit_mantissa = " << lower_limit_mantissa << endl;
    cout         << "    upper_limit_mantissa = " << upper_limit_mantissa << endl;
    cout         << "    step_mantissa        = " << step_mantissa << endl;
    cout         << "    allowed_error_float  = " << allowed_error_float << endl;
  }

  for (int i = 0; i < exp_arr_size; i++) {
    // Extract a value to be tested for the exponent part.
    input_float.e = sample_exponent_array[i];
    compare = false;

    // For that particular exponent value, go through every possible value that can be represented by the mantissa.
    for (double mant_i = lower_limit_mantissa; mant_i <= upper_limit_mantissa; mant_i += step_mantissa) {
      input_float.m = mant_i;
      if (input_float != 0) {
        test_ac_sqrt_pwl_float(input_float, output_float);
        double expected_value_float   = ((T_out)(sqrt(input_float.to_double()))).to_double();
        double actual_value_float     = output_float.to_double();

        double this_error_float;
        // If expected value is not equal to zero, calculate relative error, else, calculate absolute error.
        if (abs(expected_value_float) != 0) {this_error_float = abs( (expected_value_float - actual_value_float) / expected_value_float ) * 100.0;}
        else {this_error_float = abs(expected_value_float - actual_value_float) * 100.0;}

        if (check_monotonic) {
          // This function has the same basic working as the fixed point version, the only difference being is that now, ac_float values are considered
          if (compare == true && old_real_output > actual_value_float) {
            cout << "  Real, floating point output not monotonic at :" << endl;
            cout << "  x = " << input_float << endl;
            cout << "  y = " << output_float << endl;
            cout << "  old_real_output = " << old_real_output << endl;
            assert(false);
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
  if (passed) { printf("PASSED , max error (%f) \n", max_error_float); }
  else        { printf("FAILED , max error (%f) \n", max_error_float); }

  if (max_error_float>cumulative_max_error_float) { cumulative_max_error_float = max_error_float; }

  return 0;
}


int main(int argc, char *argv[])
{
  double max_error_fixed = 0, cmplx_max_error_fixed = 0, max_error_float = 0;
  double allowed_error_fixed = 0.5;
  double allowed_error_float = 0.5;
  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_sqrt_pwl() - Allowed error " << allowed_error_fixed << " (fixed pt), " << allowed_error_float << " (float pt)" << endl;

  // template <int Wfi, int Ifi, int Sfi>
  test_driver_fixed< 12,  3, 64, 32>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed);
  test_driver_fixed<  4,  9, 64, 32>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed);
  test_driver_fixed<  4, -2, 64, 32>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed);
  test_driver_fixed<  5,  8, 64, 32>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed);
  test_driver_fixed<  4, -2, 60, 30>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed);
  test_driver_fixed< 13,  4, 64, 32>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed);
  test_driver_fixed< 14,  3, 64, 32>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed);
  test_driver_fixed<  9,  4, 60, 30>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed);
  test_driver_fixed<  9,  2, 64, 32>(max_error_fixed, cmplx_max_error_fixed, allowed_error_fixed);

  // template <int Wfl, int Ifl, int Efl>
  test_driver_float< 10,  5, 3, 32, 16, 3>(max_error_float, allowed_error_float);
  test_driver_float< 11,  4, 2, 32, 16, 4>(max_error_float, allowed_error_float);
  test_driver_float<  9,  2, 3, 32, 15, 3>(max_error_float, allowed_error_float);
  test_driver_float<  8,  2, 4, 60, 30, 11>(max_error_float, allowed_error_float);
  test_driver_float<  7,  2, 5, 64, 32, 10>(max_error_float, allowed_error_float);
  test_driver_float<  8,  3, 4, 64, 32, 10>(max_error_float, allowed_error_float);
  test_driver_float<  6,  3, 5, 60, 30, 11>(max_error_float, allowed_error_float);
  test_driver_float< 10,  5, 8, 64, 32, 10>(max_error_float, allowed_error_float);
  test_driver_float<  9,  3, 1, 32, 16, 4>(max_error_float, allowed_error_float);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum errors observed across all data type / bit-width variations:" << endl;
  cout << "    max_error_fixed       = " << max_error_fixed << endl;
  cout << "    cmplx_max_error_fixed = " << cmplx_max_error_fixed << endl;
  cout << "    max_error_float       = " << max_error_float << endl;

  // If error limits on any tested datatype have been crossed, the test has failed
  bool test_fail = (max_error_fixed > allowed_error_fixed) || (cmplx_max_error_fixed > allowed_error_fixed) || (max_error_float > allowed_error_float);

  // Notify the user that the test was a failure.
  if (test_fail) {
    cout << "  ac_sqrt_pwl - FAILED - Error tolerance(s) exceeded" << endl;
    cout << "=============================================================================" << endl;
    return -1;
  } else {
    cout << "  ac_sqrt_pwl - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }
  return 0;
}


