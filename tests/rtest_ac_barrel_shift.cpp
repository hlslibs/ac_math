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
// ac_barrel_shift() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_barrel_shift.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench

#include <ac_math/ac_barrel_shift.h>
using namespace ac_math;

// ==============================================================================
// Test Designs
//   These simple functions allow executing the ac_barrel_shift function
//   using multiple data types at the same time. Template parameters are used to
//   configure the bit-widths of the types.

// Test for barrel shifts on unsigned ac_int inputs and outputs.
template <int W >
void test_ac_barrel_shift_us(
  const ac_int<W, false> &us_int_input,
  ac_int<W, false> &us_int_output,
  const ac_int< ac::nbits< W - 1 >::val,false> &us_shift
)
{
 us_int_output =ac_barrel_shift<W>(us_int_input, (ac_int< ac::nbits< W - 1 >::val,false>)us_shift);
}

// ==============================================================================

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

// ------------------------------------------------------------
// Helper functions for error detection and output checking.

#ifdef DEBUG

// Helper functions to pin-point error.

template <int W, int I>
void print_origin_error(
  const ac_int<W, false> input,
  const int n
)
{
  cout << "error coming from barrel shift function" << endl;
}
#endif

// Check output and make sure it is correct for real values and shifts

template <int W >
bool output_check_barrel_shift(
  const ac_int<W,false> input,
  const  ac_int< ac::nbits< W - 1 >::val,false> n,
  const ac_int<W,false> output
)
{
  // Since the output of ac_barrel_shift operations have an arithmetic
  // basis and take into account rounding and saturation, test
  // the output against the input multiplied by 2^(shift_count)

  ac_int<W,false> expected_output = (input >> (n%W)) | (input << (W-(n%W)));
  bool correct = (output == expected_output);

#ifdef DEBUG
  if (!correct) {
    cout << endl;
    cout << "barrel shift" << endl;
    print_origin_error(input, n);
    cout << "The output is not as expected." << endl;
    cout << "input           = " << input << endl;
    cout << "expected_output = " << expected_output << endl;
    cout << "output          = " << output << endl;
    cout << "n               = " << n << endl;
    assert(false);
  }
#endif

  return correct;
}

// ==============================================================================
// Functions: test_driver functions
// Description: Templatized functions that can be configured for certain bit-
//   widths of AC datatypes. They use the type information to iterate through a
//   range of valid values on that type and make sure that the output of the
//   ac_barrel_shift function is correct

// ==============================================================================
// Function: test_driver_us()
// Description: test_driver function for unsigned ac_int and ac_int inputs
//   and outputs.

template <int W>
int test_driver_us(
  bool &all_tests_pass,
  bool details = false
)
{
  ac_int<W, false>   us_int_input;
  ac_int<W, false>   us_int_output;

  double lower_limit, upper_limit, step;

  // set ranges and step size for int point testbench
  lower_limit = us_int_input.template set_val<AC_VAL_MIN>().to_double();
  upper_limit = us_int_input.template set_val<AC_VAL_MAX>().to_double();
  step        = us_int_input.template set_val<AC_VAL_QUANTUM>().to_double()*upper_limit/50;
  step        = step==0?1:step;

  cout << "TEST: ac_barrel_shift() INPUT: ";
  cout.width(38);
  cout << left << us_int_input.type_name();
  cout << "OUTPUT: ";
  cout.width(50);
  cout << left << us_int_output.type_name();
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

  for (int n = 0; n <= W; n++) {
    for (double i = lower_limit; i <= upper_limit; i += step) {
        us_int_input = i;
        test_ac_barrel_shift_us<W>(us_int_input, us_int_output, n);
        bool correct_iteration_barrel_shift = output_check_barrel_shift( us_int_input, n, us_int_output ); 

        correct = correct && correct_iteration_barrel_shift;
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
  cout << "Testing function: ac_barrel_shift(), unsigned ac_int types." << endl;

  bool all_tests_pass = true;

  // If any of the tests fail, the all_tests_pass variable will be set to false

  test_driver_us<10>(all_tests_pass);
  test_driver_us<15>(all_tests_pass);
  test_driver_us<13>(all_tests_pass);
  test_driver_us<18>(all_tests_pass);
  test_driver_us<20>(all_tests_pass);
 
  cout << "=============================================================================" << endl;
  cout << "  Testbench finished." << endl;

  // Notify the user if the test was a failure.
  if (!all_tests_pass) {
    cout << "  ac_barrel_shift - FAILED - output not correct for all test values" << endl; // LCOV_EXCL_LINE
    cout << "=============================================================================" << endl; // LCOV_EXCL_LINE
    return -1; // LCOV_EXCL_LINE
  } else {
    cout << "  ac_barrel_shift - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }

  return 0;
}




























