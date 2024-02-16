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
// ac_softmax_pwl() function using a variety of bitwidths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_softmax_pwl.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_softmax_pwl.h>
using namespace ac_math;

// ==============================================================================

#include <ac_math/ac_random.h>
#include <math.h>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;

// ------------------------------------------------------------------------------
// Helper function for absolute value calculation. This can avoid any naming conflicts
// with other absolute value functions.

double abs_double(const double &input)
{
  return input < 0 ? -input : input;
}

// ------------------------------------------------------------------------------
// Helper function for softmax calculation
template <int W, int I, bool S, ac_q_mode Q, ac_o_mode O, unsigned K_tb>
void softmax_tb(const ac_fixed<W, I, S, Q, O> (&input)[K_tb], double (&output)[K_tb])
{
  double exp_in[K_tb];
  double sum_exp = 0;
  for (unsigned i = 0; i < K_tb; i++) {
    exp_in[i] = exp(input[i].to_double());
    sum_exp += exp_in[i];
  }
  for (unsigned i = 0; i < K_tb; i++) { output[i] = exp_in[i]/sum_exp; }
}

#ifdef DEBUG
// ------------------------------------------------------------------------------
// Helper function to print ac_fixed arrays.
template<unsigned K_tb, int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
void print_arr(const ac_fixed<W, I, S, Q, O> (&input)[K_tb])
{
  for (unsigned i = 0; i < K_tb; i++) { cout << input[i].to_double() << ", "; }
  cout << endl;
}

// ------------------------------------------------------------------------------
// Helper function to print double arrays.
template <unsigned K_tb>
void print_arr(const double input[K_tb])
{
  for (unsigned i = 0; i < K_tb; i++) { cout << input[i] << ", "; }
  cout << endl;
}
#endif

// ==============================================================================
// Test Design
//   This simple function allows executing the ac_softmax_pwl() function.
//   Template parameters are used to configure the bit-widths of the
//   ac_fixed inputs/outputs.

template <int Wfi, int Ifi, int outWfi, int outIfi, unsigned K_tb>
void test_ac_softmax_pwl(
  const ac_fixed<Wfi, Ifi, true, AC_TRN, AC_WRAP> (&input)[K_tb],
  ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> (&output)[K_tb]
)
{
  ac_softmax_pwl(input, output);
}


// ===============================================================================
// Function: test_driver()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_fixed inputs. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   ac_softmax_pwl() model with the computed softmax using a standard C double type.
//   The maximum error for each type is accumulated in variables defined in the
//   calling function.

template <int Wfi, int Ifi, int outWfi, int outIfi, unsigned K_tb, unsigned n_vectors_tb>
int test_driver (
  double &cumulative_max_error_softmax,
  const double allowed_error,
  const double score_threshold
)
{
  ac_fixed<Wfi, Ifi, true, AC_TRN, AC_WRAP> input[K_tb];
  typedef ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> output_type;
  output_type output[K_tb];
  double output_expected[K_tb];

  cout << "TEST: ac_softmax_pwl() INPUT: ";
  cout.width(38);
  cout << left << input[0].type_name();
  cout << "OUTPUT: ";
  cout.width(38);
  cout << left << output[0].type_name();
  cout << "RESULT: ";

  double max_error_softmax = 0.0;

  for (unsigned i = 0; i < n_vectors_tb; i++) {
    // Write randomized testing values into the input array.
    for (unsigned j = 0; j < K_tb; j++) { ac_random(input[j]); }
    // Call ac_softmax_pwl through the test design.
    test_ac_softmax_pwl(input, output);
    // Call testbench on softmax input
    softmax_tb(input, output_expected);
    // Iterate through all logit outputs to find the expected and actual maximum output.
    for (unsigned j = 0; j < K_tb; j++) {
      // Calculate absolute error.
      double error_it = abs_double(output[j].to_double() - ((output_type)output_expected[j]).to_double());
      // If error during current iteration exceeds max_error_softmax value, store the current
      // iteration's error in the max_error_softmax variable.
      max_error_softmax = AC_MAX(error_it, max_error_softmax);
#ifdef DEBUG
      // Test should fail if absolute error exceeds tolerance.
      if (error_it > allowed_error) {
        cout << "i = " << i << endl;
        cout << "input : " << endl;
        print_arr(input);
        cout << "output : " << endl;
        print_arr(output);
        cout << "output_expected : " << endl;
        print_arr(output_expected);
        AC_ASSERT(false, "Error value exceeds threshold.");
      }
#endif
    } // for (unsigned j = 0; j < K_tb; j++)
  } // for (unsigned i = 0; i < n_vectors_tb; i++) {

  bool passed = !(max_error_softmax > score_threshold);

  if (passed) { printf("PASSED , max err (%f)\n", max_error_softmax); }
  else        { printf("FAILED , max err (%f)\n", max_error_softmax); } // LCOV_EXCL_LINE

  cumulative_max_error_softmax = AC_MAX(cumulative_max_error_softmax, max_error_softmax);

  return 0;
}

int main(int argc, char *argv[])
{
  // Define the variables used to store the max error of the softmax function, the error tolerance and the threshold which
  // the output must cross in order for the bounding box to be deemed valid.
  // Note that error metrics are all absolute value metrics.
  double max_error_softmax = 0.0;
  const double allowed_error = 0.005;
  const double score_threshold = 0.3;

  cout << "================================================================================" << endl;
  cout << "Testing function: ac_softmax_pwl() - Allowed error (absolute) " << allowed_error << endl;

  test_driver<32, 6, 64, 32, 20, 845>(max_error_softmax, allowed_error, score_threshold);
  test_driver<32, 6, 33,  1, 25, 500>(max_error_softmax, allowed_error, score_threshold);
  test_driver<22, 6, 64, 32, 15, 900>(max_error_softmax, allowed_error, score_threshold);
  test_driver<22, 6, 33,  1, 30, 300>(max_error_softmax, allowed_error, score_threshold);
  test_driver<32, 6, 32, 16, 20, 700>(max_error_softmax, allowed_error, score_threshold);
  test_driver<32, 6, 17,  1, 30, 550>(max_error_softmax, allowed_error, score_threshold);
  test_driver<22, 6, 32, 16, 25, 800>(max_error_softmax, allowed_error, score_threshold);
  test_driver<22, 6, 17,  1, 25, 400>(max_error_softmax, allowed_error, score_threshold);

  cout << "================================================================================" << endl;
  cout << "Testbench finished. Maximum error observed across all bit-width variations:" << endl;
  cout << "max_error_softmax = " << max_error_softmax << endl;

  if (max_error_softmax > score_threshold) {
    cout << "ac_softmax_pwl - FAILED - Max error exceeds threshold. Use -DDEBUG for more information." << endl; // LCOV_EXCL_LINE
    cout << "================================================================================" << endl; // LCOV_EXCL_LINE
    return (-1); // LCOV_EXCL_LINE
  }

  cout << "ac_softmax_pwl - PASSED" << endl;
  cout << "================================================================================" << endl;

  return 0;
}
