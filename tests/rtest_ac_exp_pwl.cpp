/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 1.0                                                 *
 *                                                                        *
 *  Release Date    : Fri Mar  2 16:26:42 PST 2018                        *
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
// ac_exp_pwl() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_exp_pwl.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_pow_pwl.h>
using namespace ac_math;

//==============================================================================
// Test Design
//   This simple function allows executing the ac_exp_pwl() function. Template
//   parameters are used to configure the bit-widths of the ac_fixed inputs.

template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi, bool outSfi>
void test_ac_exp_pwl(
  const ac_fixed<   Wfi,    Ifi,    Sfi, AC_TRN, AC_WRAP> &in,
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> &out_exp
)
{
  ac_exp_pwl(in, out_exp);
}

//==============================================================================

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

//==============================================================================
// Function: test_driver()
// Description: A templatized function that can be configured for certain bit-
//   widths of AC datatypes. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   piece-wise linear power model with the computed power using a
//   standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <int Wfi, int Ifi, bool Sfi, int outWfi, int outIfi, bool outSfi>
int test_driver(
  double &cumulative_max_error_exp,
  const double allowed_error,
  const double threshold,
  bool details = false
)
{
  bool passed = true;
  bool check_monotonic = true;
  double max_error_exp  = 0.0; // reset for this run

  ac_fixed<   Wfi,    Ifi,    Sfi, AC_TRN, AC_WRAP> input;
  ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> output_exp;

  typedef ac_fixed<outWfi, outIfi, outSfi, AC_TRN, AC_WRAP> T_out;

  double lower_limit, upper_limit, step;

  // set ranges and step size for fixed point testbench
  lower_limit = input.template set_val<AC_VAL_MIN>().to_double();
  upper_limit = input.template set_val<AC_VAL_MAX>().to_double();
  step        = input.template set_val<AC_VAL_QUANTUM>().to_double();

  printf("TEST: ac_exp_pwl() INPUT: ac_fixed<%2d,%2d,%5s,%7s,%7s> OUTPUT: ac_fixed<%2d,%2d,%5s,%7s,%7s>  RESULT: ",
         Wfi,Ifi,(Sfi?"true":"false"),"AC_TRN","AC_WRAP",outWfi,outIfi,(outSfi?"true":"false"),"AC_TRN","AC_WRAP");

  // Dump the test details
  if (details) {
    cout << endl;
    cout << "  Ranges for input types:" << endl;
    cout << "    lower_limit          = " << lower_limit << endl;
    cout << "    upper_limit          = " << upper_limit << endl;
    cout << "    step                 = " << step << endl;
  }

  double old_output_exp;
  bool compare_exp = false;

  for (double i = lower_limit; i < upper_limit; i += step) {
    //Set values for input.
    input = i;
    test_ac_exp_pwl(input, output_exp);

    double expected_value_exp  = ((T_out)exp(input.to_double())).to_double();
    double actual_value_exp    = output_exp.to_double();

    double this_error_exp;

    //If expected value of either output falls below the threshold, calculate absolute error instead of relative

    if (expected_value_exp > threshold) {this_error_exp = abs( (expected_value_exp - actual_value_exp) / expected_value_exp ) * 100.0;}
    else {this_error_exp = abs(expected_value_exp - actual_value_exp) * 100.0;}

    if (check_monotonic) {
      // MONOTONIC: Make sure that function is monotonic. Compare old value (value of previous iteration) with current value. Since the exponential function we
      //are testing is an increasing function, and our testbench value keeps incrementing or remains the same (in case of saturation), we expect the
      //old value to be lesser than or equal to the current one.

      //Update the old value
      old_output_exp = actual_value_exp;
      //Once an old value has been stored, i.e. towards the end of the first iteration, this value is set to true.
      compare_exp = true;

      //same thing as above, but for the natural exponential.
      if (compare_exp) {
        if (old_output_exp > actual_value_exp) {
          cout << "FILE : " << __FILE__ << ", LINE : " << __LINE__ << endl;
          cout << "exp output not monotonic at :" << endl;
          cout << "x = " << input << endl;
          cout << "old_output_exp = " << old_output_exp << endl;
          //assert(false);   //Uncomment if you want the program to stop whenever monotonicity is violated.
        }
      }
      old_output_exp = actual_value_exp;
      compare_exp = true;
    }

    if (this_error_exp > max_error_exp) {max_error_exp = this_error_exp;}
  }
  if (passed) { printf("PASSED , max err (%f exp)\n", max_error_exp); }
  else        { printf("FAILED , max err (%f exp)\n", max_error_exp); }

  if (max_error_exp>cumulative_max_error_exp) { cumulative_max_error_exp = max_error_exp; }

  return 0;
}


int main(int argc, char *argv[])
{
  double max_error_exp = 0;
  double allowed_error = 0.5;
  double threshold = 0.005;
  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_exp_pwl() - Allowed error " << allowed_error << endl;

  // template <int Wfi, int Ifi, int Sfi>
  test_driver< 12,  3,  true, 64, 32, false>(max_error_exp, allowed_error, threshold);
  test_driver<  4,  2,  true, 64, 32, false>(max_error_exp, allowed_error, threshold);
  test_driver<  3,  5,  true, 64, 32, false>(max_error_exp, allowed_error, threshold);
  test_driver<  4, -2,  true, 64, 32, false>(max_error_exp, allowed_error, threshold);
  test_driver<  3,  5,  true, 64, 32, false>(max_error_exp, allowed_error, threshold);
  test_driver<  2,  5,  true, 64, 32, false>(max_error_exp, allowed_error, threshold);
  test_driver< 16,  5,  true, 64, 32, false>(max_error_exp, allowed_error, threshold);
  test_driver< 12,  4, false, 64, 32, false>(max_error_exp, allowed_error, threshold);
  test_driver<  4,  2, false, 64, 32, false>(max_error_exp, allowed_error, threshold);
  test_driver<  4, -2, false, 60, 30, false>(max_error_exp, allowed_error, threshold);
  test_driver<  3,  4, false, 64, 32, false>(max_error_exp, allowed_error, threshold);
  test_driver<  1,  5, false, 61, 33, false>(max_error_exp, allowed_error, threshold);
  test_driver<  2,  5, false, 64, 32, false>(max_error_exp, allowed_error, threshold);
  test_driver< 16,  4, false, 64, 32, false>(max_error_exp, allowed_error, threshold);
  test_driver< 16,  0, false, 64, 32, false>(max_error_exp, allowed_error, threshold);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum errors observed across all bit-width variations:" << endl;
  cout << "    max_error_exp  = " << max_error_exp  << endl;

  // If error limits on any tested datatype have been crossed, the test has failed
  bool test_fail = (max_error_exp > allowed_error);

  // Notify the user that the test was a failure.
  if (test_fail) {
    cout << "  ac_exp_pwl - FAILED - Error tolerance(s) exceeded" << endl;
    cout << "=============================================================================" << endl;
    return -1;
  } else {
    cout << "  ac_exp_pwl - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }
  return 0;
}




































