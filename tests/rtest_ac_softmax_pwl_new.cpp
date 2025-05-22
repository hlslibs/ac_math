/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 3.8                                                 *
 *                                                                        *
 *  Release Date    : Tue May 13 15:34:32 PDT 2025                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.8.1                                               *
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
// ac_softmax_pwl_new() function using a variety of bitwidths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_softmax_pwl_new.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_softmax_pwl_new.h>
using namespace ac_math;

// ==============================================================================

#include <limits>
#include <random>
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

#ifdef DEBUG
// ------------------------------------------------------------------------------
// Helper function to print ac_fixed arrays.
template<unsigned K_tb, int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
void print_arr(const ac_fixed<W, I, S, Q, O> (&input)[K_tb])
{
  for (const auto& input_elem : input) { cout << input_elem.to_double() << ", "; }
  cout << endl;
}

// ------------------------------------------------------------------------------
// Helper function to print double arrays.
template <unsigned K_tb>
void print_arr(const double (&input)[K_tb])
{
  for (const auto& input_elem : input) { cout << input_elem << ", "; }
  cout << endl;
}
#endif

// ------------------------------------------------------------------------------
// Helper function for softmax calculation
template <int W, int I, bool S, ac_q_mode Q, ac_o_mode O, unsigned K_tb>
#ifdef DEBUG
// If debugging, we may need to print the exp_in and sum_exp arrays externally.
void softmax_tb(const ac_fixed<W, I, S, Q, O> (&input)[K_tb], double (&exp_in)[K_tb], double &sum_exp,  double (&output)[K_tb])
#else
void softmax_tb(const ac_fixed<W, I, S, Q, O> (&input)[K_tb], double (&output)[K_tb])
#endif
{
  #ifndef DEBUG
  double exp_in[K_tb];
  double sum_exp;
  #endif
  sum_exp = 0.0;
  for (unsigned i = 0; i < K_tb; i++) {
    exp_in[i] = exp(input[i].to_double());
    sum_exp += exp_in[i];
  }
  for (unsigned i = 0; i < K_tb; i++) { output[i] = exp_in[i]/sum_exp; }
}

// ==============================================================================
// Test Design
//   This simple function allows executing the ac_softmax_pwl_new() function.
//   Template parameters are used to configure the bit-widths of the
//   ac_fixed inputs/outputs.

template <int Wfi, int Ifi, int outWfi, int outIfi, unsigned K_tb>
void test_ac_softmax_pwl_new(
  const ac_fixed<Wfi, Ifi, true, AC_TRN, AC_WRAP> (&input)[K_tb],
  ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> (&output)[K_tb]
)
{
  ac_softmax_pwl_new(input, output);
}

// ==============================================================================
// Function: assign_in_val
//   Assign values to input vector.

template <int W, int I, ac_q_mode Q, ac_o_mode O, unsigned K_tb>
void assign_in_val(const int i, ac_fixed<W, I, true, Q, O> (&in)[K_tb]) {
  typedef ac_fixed<W, I, true, Q, O> input_type;

  static default_random_engine generator;

  double type_supported_max = value<AC_VAL_MAX>(input_type{}).to_double();
  double dbl_supported_max = log(std::numeric_limits<double>::max()/double(K_tb));
  double type_supported_min = value<AC_VAL_MIN>(input_type{}).to_double();
  double dbl_supported_min = log(std::numeric_limits<double>::denorm_min());

  double max_in_val_1 = AC_MIN(type_supported_max, dbl_supported_max);
  double min_in_val_1 = max_in_val_1/4.0;
  double min_in_val_2 = AC_MAX(type_supported_min, dbl_supported_min);
  double max_in_val_2 = min_in_val_2/4.0;
  double min_in_val_3 = -2.0;
  double max_in_val_3 = 2.0;

  static uniform_real_distribution<double> distribution_1(min_in_val_1, max_in_val_1);
  static uniform_real_distribution<double> distribution_2(min_in_val_2, max_in_val_2);
  static uniform_real_distribution<double> distribution_3(min_in_val_3, max_in_val_3);

  int mod_den = (type_supported_max < dbl_supported_max && type_supported_min > dbl_supported_min) ? 4 : 3;

  for (input_type& in_elem : in) {
    switch (i%mod_den) {
      case 0:
        in_elem = distribution_1(generator);
        break;
      case 1:
        in_elem = distribution_2(generator);
        break;
      case 2:
        in_elem = distribution_3(generator);
        break;
      default: // 3
        ac_random(in_elem);
    }
  }
}

// ===============================================================================
// Function: test_driver()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_fixed inputs. It uses the type information to iterate through a
//   range of valid values on that type in order to compare the precision of the
//   ac_softmax_pwl_new() model with the computed softmax using a standard C double type.
//   The maximum error for each type is accumulated in variables defined in the
//   calling function.

template <int Wfi, int Ifi, int outWfi, int outIfi, unsigned K_tb, unsigned n_vectors_tb>
int test_driver (
  double &cumulative_max_error_softmax,
  const double allowed_error,
  const double score_threshold
)
{
  typedef ac_fixed<Wfi, Ifi, true, AC_TRN, AC_WRAP> input_type;
  typedef ac_fixed<outWfi, outIfi, false, AC_TRN, AC_WRAP> output_type;
  input_type input[K_tb];
  output_type output[K_tb];
  double output_expected[K_tb];

  cout << "TEST: ac_softmax_pwl_new() INPUT: ";
  cout.width(38);
  cout << left << input[0].type_name();
  cout << "OUTPUT: ";
  cout.width(38);
  cout << left << output[0].type_name();
  cout << "K_tb = " << K_tb << "  ";
  cout << "RESULT: ";

  double max_error_softmax = 0.0;

  for (unsigned i = 0; i < n_vectors_tb; i++) {
    assign_in_val(i, input);
    // Call ac_softmax_pwl_new through the test design.
    test_ac_softmax_pwl_new(input, output);
    #ifdef DEBUG
    double exp_in_debug[K_tb];
    double sum_exp_debug;
    #endif
    // Call testbench reference function with softmax input
    #ifdef DEBUG
    softmax_tb(input, exp_in_debug, sum_exp_debug, output_expected);
    #else
    softmax_tb(input, output_expected);
    #endif
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
        cout << "j = " << j << endl;
        cout << "output[j] = " << output[j] << endl;
        cout << "output_expected[j] = " << output_expected[j] << endl;
        cout << "error_it = " << error_it << endl;
        cout << "allowed_error = " << allowed_error << endl;
        cout << "input : " << endl;
        print_arr(input);
        cout << "output : " << endl;
        print_arr(output);
        cout << "exp_in_debug : " << endl;
        print_arr(exp_in_debug);
        cout << "sum_exp_debug = " << sum_exp_debug << endl;
        cout << "output_expected : " << endl;
        print_arr(output_expected);
        AC_ASSERT(false, "Error value exceeds threshold.");
      }
      #endif
    }
  }

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
  const double allowed_error = 0.002;
  const double score_threshold = 0.3;

  cout << "================================================================================" << endl;
  cout << "Testing function: ac_softmax_pwl_new() - Allowed error (absolute) " << allowed_error << endl;
  
  // template <int Wfi, int Ifi, int outWfi, int outIfi, unsigned K_tb, unsigned n_vectors_tb>
  test_driver<32, 6, 64, 32, 15, 845>(max_error_softmax, allowed_error, score_threshold);
  test_driver<32, 6, 64, 32, 20, 845>(max_error_softmax, allowed_error, score_threshold);
  test_driver<32, 6, 33,  1, 25, 500>(max_error_softmax, allowed_error, score_threshold);
  test_driver<22, 6, 64, 32, 15, 900>(max_error_softmax, allowed_error, score_threshold);
  test_driver<22, 6, 33,  1, 30, 300>(max_error_softmax, allowed_error, score_threshold);
  test_driver<32, 6, 32, 16, 20, 700>(max_error_softmax, allowed_error, score_threshold);
  test_driver<32, 6, 17,  1, 30, 550>(max_error_softmax, allowed_error, score_threshold);
  test_driver<22, 6, 32, 16, 25, 800>(max_error_softmax, allowed_error, score_threshold);
  test_driver<22, 6, 17,  1, 25, 400>(max_error_softmax, allowed_error, score_threshold);

  test_driver<32, 16, 64, 32, 15, 845>(max_error_softmax, allowed_error, score_threshold);
  test_driver<32, 16, 64, 32, 20, 845>(max_error_softmax, allowed_error, score_threshold);
  test_driver<32, 16, 33,  1, 25, 500>(max_error_softmax, allowed_error, score_threshold);
  test_driver<22, 16, 64, 32, 15, 900>(max_error_softmax, allowed_error, score_threshold);
  test_driver<22, 16, 33,  1, 30, 300>(max_error_softmax, allowed_error, score_threshold);
  test_driver<32, 16, 32, 16, 20, 700>(max_error_softmax, allowed_error, score_threshold);
  test_driver<32, 16, 17,  1, 30, 550>(max_error_softmax, allowed_error, score_threshold);
  test_driver<22, 16, 32, 16, 25, 800>(max_error_softmax, allowed_error, score_threshold);
  test_driver<22, 16, 17,  1, 25, 400>(max_error_softmax, allowed_error, score_threshold);

  cout << "================================================================================" << endl;
  cout << "Testbench finished. Maximum error observed across all bit-width variations:" << endl;
  cout << "max_error_softmax = " << max_error_softmax << endl;

  if (max_error_softmax > score_threshold) {
    cout << "ac_softmax_pwl_new - FAILED - Max error exceeds threshold. Use -DDEBUG for more information." << endl; // LCOV_EXCL_LINE
    cout << "================================================================================" << endl; // LCOV_EXCL_LINE
    return (-1); // LCOV_EXCL_LINE
  }

  cout << "ac_softmax_pwl_new - PASSED" << endl;
  cout << "================================================================================" << endl;

  return 0;
}
