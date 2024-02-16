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
// ac_array() class.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_array.cpp -o design
//   ./design

#include <ac_int.h>

// Include the AC Math function that is exercised with this testbench
#include <ac_array.h>

// ==============================================================================
// Test Designs
//   These simple functions allow executing the ac_abs() function
//   using multiple data types at the same time. Template parameters are used
//   to configure the bit-widths of the types.

typedef ac_int<8,false> dtype_t;

template <unsigned int D1, unsigned int D2, unsigned int D3>
void test_ac_array(
  ac_array<dtype_t,D1>       array1D,
  ac_array<dtype_t,D1,D2>    array2D,
  ac_array<dtype_t,D1,D2,D3> array3D,
  dtype_t&                   sum_out)
{
  dtype_t sum = 0;
  for (int i=0; i<D1; i++ ) {
    for (int j=0; j<D2; j++) {
      for (int k=0; k<D3; k++) {
        sum += array3D[i][j][k];
      }
      sum += array2D[i][j];
    }
    sum += array1D[i];
  }
  sum_out = sum;
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
//   range of valid values on that type and make sure that the input and
//   output of the ac_abs function match each other in terms of absolute value.

// ==============================================================================
// Function: test_driver()
// Description: test_driver function for ac_int inputs and outputs.

template <unsigned int D1, unsigned int D2, unsigned int D3>
int test_driver(
  bool &all_tests_pass,
  bool details = false
)
{
  ac_array<dtype_t,D1>       a_1d;
  ac_array<dtype_t,D1,D2>    a_2d;
  ac_array<dtype_t,D1,D2,D3> a_3d;
  for (int i=0;i<D1;i++) { a_1d[i] = i; }
  for (int i=0;i<D1;i++) for (int j=0;j<D2;j++) { a_2d[i][j] = i*D2+j; }
  for (int i=0;i<D1;i++) for (int j=0;j<D2;j++) for (int k=0;k<D3;k++) { a_3d[i][j][k] = (i*D2+j)*D3+k; }

  cout << "TEST: ac_array() INPUT: " << std::endl;;
  cout << "  Print 1D: a_1d: " << a_1d << std::endl;
  cout << "  Print 2D: a_2d: " << a_2d << std::endl;
  cout << "  Print 3D: a_3d: " << a_3d << std::endl;

  bool correct = true;
  ac_array<dtype_t,D1>        b_1d;
  ac_array<dtype_t,D1,D2>     b_2d;
  ac_array<dtype_t,D1,D2,D3>  b_3d;
  b_1d = a_1d;
  b_2d = a_2d;
  b_3d = a_3d;

  // operator==
  if (a_1d == b_1d) {
    std::cout << "  1D operator==() correct" << std::endl;
  } else {
    correct = correct && false;
    std::cerr << "  1D operator==() error" << std::endl;
  }
  if (a_2d == b_2d) {
    std::cout << "  2D operator==() correct" << std::endl;
  } else {
    correct = correct && false;
    std::cerr << "  2D operator==() error" << std::endl;
  }
  if (a_3d == b_3d) {
    std::cout << "  3D operator==() correct" << std::endl;
  } else {
    correct = false;
    std::cerr << "  3D operator==() error" << std::endl;
  }

  // operator!=
  b_1d[2] = 255;
  b_2d[2][2] = 255;
  b_3d[2][2][1] = 255;
  if (a_1d != b_1d) {
    std::cout << "  1D operator!=() correct" << std::endl;
  } else {
    correct = false;
    std::cerr << "  1D operator!=() error" << std::endl;
  }
  if (a_2d != b_2d) {
    std::cout << "  2D operator!=() correct" << std::endl;
  } else {
    correct = false;
    std::cerr << "  2D operator!=() error" << std::endl;
  }
  if (a_3d != b_3d) {
    std::cout << "  3D operator!=() correct" << std::endl;
  } else {
    correct = false;
    std::cerr << "  3D operator!=() error" << std::endl;
  }
  dtype_t a_sum;
  test_ac_array(a_1d,a_2d,a_3d,a_sum);
  if (a_sum == 92) {
    std::cout << "  test_ac_array<" << D1 << "," << D2 << "," << D3 << "> returned 92 - correct" << std::endl;
  } else {
    correct = false;
    std::cout << "  test_ac_array<" << D1 << "," << D2 << "," << D3 << "> returned " << a_sum << " - error" << std::endl;
  }

  // operator=
  a_1d = 7; a_2d = 9; a_3d = 11;
  bool set_1d, set_2d, set_3d;
  set_1d = set_2d = set_3d = true;
  for (int i=0;i<D1;i++) { if (a_1d[i] != 7) set_1d = false; }
  for (int i=0;i<D1;i++) for (int j=0;j<D2;j++) { if (a_2d[i][j] != 9) set_2d = false; }
  for (int i=0;i<D1;i++) for (int j=0;j<D2;j++) for (int k=0;k<D3;k++) { if (a_3d[i][j][k] != 11) set_3d = false; }
  if (set_1d) {
    cout << "  1D operator=() correct" << endl;
  } else {
    correct = false;
    cerr << "  1D operator=() incorrect" << endl;
  }
  if (set_2d) {
    cout << "  2D operator=() correct" << endl;
  } else {
    correct = false;
    cerr << "  2D operator=() incorrect" << endl;
  }
  if (set_3d) {
    cout << "  3D operator=() correct" << endl;
  } else {
    correct = false;
    cerr << "  3D operator=() incorrect" << endl;
  }

  // clearAll()
  a_1d.clearAll(); a_2d.clearAll(); a_3d.clearAll();
  bool clear_1d, clear_2d, clear_3d;
  clear_1d = clear_2d = clear_3d = true;
  for (int i=0;i<D1;i++) { if (a_1d[i] != 0) clear_1d = false; }
  for (int i=0;i<D1;i++) for (int j=0;j<D2;j++) { if (a_2d[i][j] != 0) clear_2d = false; }
  for (int i=0;i<D1;i++) for (int j=0;j<D2;j++) for (int k=0;k<D3;k++) { if (a_3d[i][j][k] != 0) clear_3d = false; }
  if (clear_1d) {
    cout << "  1D clearAll() correct" << endl;
  } else {
    correct = false;
    cerr << "  1D clearAll() incorrect"  << endl;
  }
  if (clear_2d) {
    cout << "  2D clearAll() correct" << endl;
  } else {
    correct = false;
    cerr << "  2D clearAll() incorrect" << endl;
  }
  if (clear_3d) {
    cout << "  3D clearAll() correct" << endl;
  } else {
    correct = false;
    cerr << "  3D clearAll() incorrect" << endl;
  }

  if (correct) { printf("PASSED\n"); }
  else         { printf("FAILED\n"); } // LCOV_EXCL_LINE

  all_tests_pass = all_tests_pass && correct;

  return 0;
}

int main(int argc, char *argv[])
{
  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_array()" << endl;

  // If any of the tests fail, the all_tests_pass variable will be set to false
  bool all_tests_pass = true;

  test_driver<4,3,2>(all_tests_pass);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished." << endl;

  // Notify the user whether or not the test was a failure.
  if (!all_tests_pass) {
    cout << "  ac_array - FAILED - Output not correct for all test values" << endl; // LCOV_EXCL_LINE
    cout << "=============================================================================" << endl; // LCOV_EXCL_LINE
    return -1; // LCOV_EXCL_LINE
  } else {
    cout << "  ac_array - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }

  return 0;
}

