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
// ac_arccos_cordic() function using a variety of bit-widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_arccos_cordic.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_arccos_cordic.h>
using namespace ac_math;

//==============================================================================
// Test Design
//   This simple function allows executing the ac_arccos_cordic() function.
//   Template parameters are used to configure the bit-widths of the types.

template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi, bool outSfi>
void test_ac_arccosine_cordic_fixed(
  const    ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP>       &in,
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP>  &out_cos
)
{
  ac_arccos_cordic(in, out_cos);
}

// Test Design for real floating point values.
template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
void test_ac_arccosine_cordic_float(
  const ac_float<Wfl, Ifl, Efl, AC_TRN>    &in,
  ac_float<outWfl, outIfl, outEfl, AC_TRN> &out_sin
)
{
  ac_arccos_cordic(in, out_sin);
}

//==============================================================================

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

//==============================================================================
// Function: test_driver_fixed()
// Description: A templatized function that can be configured for certain bit-
//   widths of the fixed point AC datatype. It uses the type information to
//   iterate through a range of valid values on that type in order to compare
//   the precision of the cordic based arccosine model with arccosine using
//   the standard C math library. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <int Wfi, int Ifi, int Sfi, int outWfi, int outIfi, bool outSfi>
int test_driver_fixed(
  double &cummulative_max_error_arccosine,
  const double allowed_error_fixed,
  bool details = false
)
{
  bool passed = true;
  double max_error_arccosine   = 0.0; // reset for this run

  ac_fixed<Wfi, Ifi, Sfi, AC_TRN, AC_WRAP>                 input;
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP>      out_cos;

  typedef ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> T_out;

  double lower_limit, upper_limit, step;

  // set ranges and step size for the testbench
  lower_limit   = input.template set_val<AC_VAL_MIN>().to_double();
  upper_limit   = input.template set_val<AC_VAL_MAX>().to_double();
  step          = input.template set_val<AC_VAL_QUANTUM>().to_double();

  cout << "TEST: ac_arccos_cordic() INPUT: ";
  cout.width(38);
  cout << left << input.type_name();
  cout << "        OUTPUTS: ";
  cout.width(38);
  cout << left << out_cos.type_name();
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl; // LCOV_EXCL_LINE
    cout << "  Ranges for input types:" << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit    = " << lower_limit << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit    = " << upper_limit << endl; // LCOV_EXCL_LINE
    cout << "    step           = " << step << endl; // LCOV_EXCL_LINE
  }

  for (double i = -1; i < 1; i += step) {
    // Set values for input.
    input = i;
    test_ac_arccosine_cordic_fixed(input, out_cos);

    double expected_arccosine_value   = ((T_out)(acos(input.to_double())/M_PI)).to_double();
    double actual_arccosine_value       = out_cos.to_double();
    double this_error_arccosine;

    this_error_arccosine   = fabs(expected_arccosine_value - actual_arccosine_value) * 100.0;

    if (this_error_arccosine > max_error_arccosine) {max_error_arccosine = this_error_arccosine;}

  }

  passed = (max_error_arccosine < allowed_error_fixed);

  if (passed) { printf("PASSED , max err (%f arccos) \n", max_error_arccosine); }
  else        { printf("FAILED , max err (%f arccos) \n", max_error_arccosine); } // LCOV_EXCL_LINE

  if (max_error_arccosine>cummulative_max_error_arccosine) { cummulative_max_error_arccosine = max_error_arccosine; }

  return 0;
}

// ===============================================================================
// Function: test_driver_float()
// Description: test_driver function for ac_float inputs and outputs.

template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
int test_driver_float(
  double &cummulative_max_error_arccosine_float,
  const double allowed_error_float,
  bool details = false
)
{
  bool passed = true;
  double max_error_arccosine   = 0.0; // reset for this run
  double this_error_arccosine;

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

  cout << "TEST: ac_arccos_cordic() INPUT: ";
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

        T_out arccos; 
        double expected_arccosine_value, actual_arccosine_value;
        if (input_float <= 1) 
        { 
          test_ac_arccosine_cordic_float(input_float, arccos);
          expected_arccosine_value   = ((T_out)(acos(input_float.to_double())/M_PI)).to_double();
          actual_arccosine_value     = arccos.to_double();

          this_error_arccosine   = fabs(expected_arccosine_value - actual_arccosine_value) * 100.0;

          if (this_error_arccosine > max_error_arccosine) {max_error_arccosine = this_error_arccosine;}
        }

        // Repeat the same, but this time with a negative input.
        ac_fixed<Wfl, Ifl, true> input_mant_neg = -input_mant;
        T_in input_float_neg(input_mant_neg, sample_exponent_array[i]);
        if ((input_float_neg.mantissa() != input_mant_neg || input_float_neg.exp() != sample_exponent_array[i]) && mant_i != 0) {
          cout << "input_mant_neg was not normalized correctly." << endl;
          assert(false);
        }

        if (input_float_neg >= -1)
        {
          test_ac_arccosine_cordic_float(input_float, arccos);
          expected_arccosine_value   = ((T_out)(acos(input_float.to_double())/M_PI)).to_double();
          actual_arccosine_value     = arccos.to_double();

          this_error_arccosine   = fabs(expected_arccosine_value - actual_arccosine_value) * 100.0;

          if (this_error_arccosine > max_error_arccosine) {max_error_arccosine = this_error_arccosine;}
        }


  #ifdef DEBUG
        if (this_error_arccosine > allowed_error_float) {
          cout << endl;
          cout << "  Error exceeds tolerance" << endl;
          cout << "  input_float      = " << input_float << endl;
          cout << "  expected_arccosine_value   = " << expected_arccosine_value << endl;
          cout << "  actual_arccosine_value     = " << actual_arccosine_value << endl;
          cout << "  this_error_arccosine = " << this_error_arccosine << endl;
        }
  #endif
      }
  }

      passed = (max_error_arccosine < allowed_error_float);

      if (passed) { printf("PASSED , max err (%f arccos) \n", max_error_arccosine); }
      else        { printf("FAILED , max err (%f arccos) \n", max_error_arccosine); } // LCOV_EXCL_LINE

      if (max_error_arccosine>cummulative_max_error_arccosine_float) { cummulative_max_error_arccosine_float = max_error_arccosine; }

      return 0;
}

int main(int argc, char *argv[])
{
  double max_error_arccosine_fixed = 0, max_error_arccosine_float = 0;
  double allowed_error = 0.1;
  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_arccos_cordic() - Allowed error " << allowed_error << endl;

  // template <int Wfi, int Ifi, int Sfi>
  test_driver_fixed< 12,  1,  true, 24, 2, false>(max_error_arccosine_fixed, allowed_error);
  test_driver_fixed<  2,  0,  true, 27, 2, false>(max_error_arccosine_fixed, allowed_error);
  test_driver_fixed< 12,  0,  true, 24, 3, false>(max_error_arccosine_fixed, allowed_error);
  test_driver_fixed<  4,  1,  true, 25, 3, false>(max_error_arccosine_fixed, allowed_error);
  test_driver_fixed<  5, -2,  true, 30, 4, false>(max_error_arccosine_fixed, allowed_error);
  test_driver_fixed<  3, -2,  true, 28, 2, false>(max_error_arccosine_fixed, allowed_error);
  test_driver_fixed<  2, -2,  true, 32, 2, false>(max_error_arccosine_fixed, allowed_error);
  test_driver_fixed< 11,  1,  true, 34, 2, false>(max_error_arccosine_fixed, allowed_error);
  test_driver_fixed<  8, -3,  true, 25, 2, false>(max_error_arccosine_fixed, allowed_error);
  test_driver_fixed<  9, -3,  true, 26, 2, false>(max_error_arccosine_fixed, allowed_error);
  test_driver_fixed< 14,  2,  true, 24, 2, false>(max_error_arccosine_fixed, allowed_error);
  test_driver_fixed< 12,  2,  true, 24, 2, false>(max_error_arccosine_fixed, allowed_error);

  // template <int Wfl, int Ifl, int Efl, int outWfl, int outIfl, int outEfl>
  test_driver_float<10,-2, 1, 20, 3, 9>(max_error_arccosine_float, allowed_error);
  test_driver_float< 9,-3, 1, 20, 3, 9>(max_error_arccosine_float, allowed_error);
  test_driver_float<12, 3, 1, 24, 3, 9>(max_error_arccosine_float, allowed_error);
  test_driver_float< 5,-1, 1, 25, 3, 9>(max_error_arccosine_float, allowed_error);
  test_driver_float< 3,-1, 1, 25, 3, 9>(max_error_arccosine_float, allowed_error);
  test_driver_float< 8,-2, 1, 25, 3, 9>(max_error_arccosine_float, allowed_error);


  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum errors observed across all bit-width variations:" << endl;
  cout << "    max_error_arccosine_fixed       = " << max_error_arccosine_fixed   << endl;
  cout << "    max_error_arccosine_float       = " << max_error_arccosine_float   << endl;

  // If error limits on any tested datatype have been crossed, the test has failed
  bool test_fail = max_error_arccosine_fixed > allowed_error || max_error_arccosine_float > allowed_error;

  // Notify the user that the test was a failure.
  if (test_fail) {
    cout << "  ac_arccos_cordic - FAILED - Error tolerance(s) exceeded" << endl; // LCOV_EXCL_LINE
    cout << "=============================================================================" << endl; // LCOV_EXCL_LINE
    return -1; // LCOV_EXCL_LINE
  } else {
    cout << "  ac_arccos_cordic - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }
  return 0;
}


