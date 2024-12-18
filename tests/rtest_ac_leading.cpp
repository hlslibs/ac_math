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

#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>

#include <mc_scverify.h>
#include <ac_math/ac_leading.h>

#include <ac_math/ac_random.h>

using namespace std;

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Display the binary version of the test number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
template<class IN_TYPE, unsigned W>
void bin(IN_TYPE input_number)
{
  for (unsigned int i = W; i>0; --i)
  { cout << std::setw(2) << i-1 << " " ; }
  cout << endl;
  for (unsigned int i = W; i>0; --i)
  { cout << std::setw(2) << input_number[i-1] << " " ; }
  cout << endl;
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Generate Random number to test the design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
template<class IN_TYPE>
static void throw_dice(IN_TYPE &a)
{
  ac_random(a);
  cout << "Test Number " << a << endl;
}

//********************************************************************INTEGER******************************************************************//
//--------------------------------------------------------------------UNSIGNED-----------------------------------------------------------------//
/*##########################################
Reference design
###########################################*/
template<int W>
unsigned int leading1_uint_ref(ac_int<W,false> input)
{
  for (int i=W-1; i>=0; i--) {
    if (input[i]==1) { return i; }
  }
  return 0;
}

template<int W>
unsigned int leading0_uint_ref(ac_int<W,false> input)
{
  for (int i=W-1; i>=0; i--) {
    if (input[i]==0) { return i; }
  }
  return 0;
}

template<int W>
unsigned int trailing1_uint_ref(ac_int<W,false> input)
{
  for (int i=0; i<=W-1; i++) {
    if (input[i]==1) { return i; }
  }
  return 0;
}

template<int W>
unsigned int trailing0_uint_ref(ac_int<W,false> input)
{
  for (int i=0; i<=W-1; i++) {
    if (input[i]==0) { return i; }
  }
  return 0;
}

/*######### SIGNED INT ##############*/
template<int W>
unsigned int leading1_sint_ref(ac_int<W,true> input)
{
  for (int i=W-2; i>=0; i--) {
    if (input[i]==1) { return i; }
  }
  return 0;
}

template<int W>
unsigned int leading0_sint_ref(ac_int<W,false> input)
{
    // Start the loop from W-2 to ignore the MSB which is at W-1
    for (int i = W-2; i >= 0; i--) {
        if (input[i] == 0) {
            return i;  // Return the index of the leading zero from LSB, excluding the MSB
        }
    }
}

template<int W>
unsigned int trailing1_sint_ref(ac_int<W,false> input)
{
  for (int i=0; i<=W-2; i++) {
    if (input[i]==1) { return i; }
  }
  return 0;
}

template<int W>
unsigned int trailing0_sint_ref(ac_int<W,false> input)
{
  for (int i=0; i<=W-2; i++) {
    if (input[i]==0) { return i; }
  }
  return 0;
}

//########### FIXED UNSIGNED #############
template<int W, int I>
unsigned int leading1_ufix_ref(ac_fixed<W,I,false> input)
{
  for (int i=W-1; i>=0; i--) {
    if (input[i]==1) { return i; }
  }
  return 0;
}

template<int W, int I>
unsigned int leading0_ufix_ref(ac_fixed<W,I,false> input)
{
  for (int i=W-1; i>=0; i--) {
    if (input[i]==0) { return i; }
  }
  return 0;
}

template<int W, int I>
unsigned int trailing1_ufix_ref(ac_fixed<W,I,false> input)
{
  for (int i=0; i<=W-1; i++) {
    if (input[i]==1) { return i; }
  }
  return 0;
}

template<int W, int I>
unsigned int trailing0_ufix_ref(ac_fixed<W,I,false> input)
{
  for (int i=0; i<=W-1; i++) {
    if (input[i]==0) { return i; }
  }
  return 0;
}

/*######### SIGNED Fixed ##############*/
template<int W, int I>
unsigned int leading1_sfix_ref(ac_fixed<W,I,true> input)
{
  for (int i=W-2; i>=0; i--) {
    if (input[i]==1) { return i; }
  }
  return 0;
}

template<int W, int I>
unsigned int leading0_sfix_ref(ac_fixed<W,I,true> input)
{
    // Start the loop from W-2 to ignore the MSB which is at W-1
    for (int i = W-2; i >= 0; i--) {
        if (input[i] == 0) {
            return i;  // Return the index of the leading zero from LSB, excluding the MSB
        }
    }
}

template<int W, int I>
unsigned int trailing1_sfix_ref(ac_fixed<W,I,true> input)
{
  for (int i=0; i<=W-2; i++) {
    if (input[i]==1) { return i; }
  }
  return 0;
}

template<int W, int I>
unsigned int trailing0_sfix_ref(ac_fixed<W,I,true> input)
{
  for (int i=0; i<=W-2; i++) {
    if (input[i]==0) { return i; }
  }
  return 0;
}


template<int W>
void test_uint()
{
  bool test_pass = true;
  ac_int<W,false> input, MAX_VAL;
  MAX_VAL.template set_val<AC_VAL_MAX>();
  bool flag;
  // std::cout << "MAX_VAL = " << MAX_VAL << std::endl;
  for (unsigned int i = 0; i <= MAX_VAL; i++) {
    input =i;
    // std::cout << "Input = " << input << std::endl;
    int result_leading_1 = ac_math::leading1<W>(input,flag);
    int result_leading_0 = ac_math::leading0<W>(input,flag);
    int result_trailing_1 = ac_math::trailing1<W>(input,flag);
    int result_trailing_0 = ac_math::trailing0<W>(input,flag);

    int expected_leading_1 = leading1_uint_ref<W>(input);
    int expected_leading_0 = leading0_uint_ref<W>(input);
    int expected_trailing_1 = trailing1_uint_ref<W>(input);
    int expected_trailing_0 = trailing0_uint_ref<W>(input);
	
	// std::cout << " ---------------- " << std::endl;
    // bin<ac_int<W,false>,W>(input);
    // // Printing results and expected values
    // std::cout << "Input: " << input << std::endl;
    // std::cout << "Actual vs Expected Leading 1 index: " << result_leading_1 << " vs " << expected_leading_1 << std::endl;
    // std::cout << "Actual vs Expected Leading 0 index: " << result_leading_0 << " vs " << expected_leading_0 << std::endl;
    // std::cout << "Actual vs Expected Trailing 1 index: " << result_trailing_1 << " vs " << expected_trailing_1 << std::endl;
    // std::cout << "Actual vs Expected Trailing 0 index: " << result_trailing_0 << " vs " << expected_trailing_0 << std::endl;
	
    if (result_leading_1 != expected_leading_1) {
      test_pass = false;
      std::cout << "Test failed for leading1 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_leading_1 << ", got " << result_leading_1 << std::endl;
	  break;				
    }
    if (result_trailing_1 != expected_trailing_1) {
      test_pass = false;
      std::cout << "Test failed for trailing1 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_trailing_1 << ", got " << result_trailing_1 << std::endl;
	  break;				

    }
    if (result_leading_0 != expected_leading_0) {
      test_pass = false;
      std::cout << "Test failed for leading0 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_leading_0 << ", got " << result_leading_0 << std::endl;
	  break;				

    }
    if (result_trailing_0 != expected_trailing_0) {
      test_pass = false;
      std::cout << "Test failed for trailing0 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_trailing_0 << ", got " << result_trailing_0 << std::endl;
	  break;				
    }
  }
  if (test_pass) {
    std::cout << "All tests for unsigned integer values passed for W = " << W << std::endl;
  }
}

template<int W>
void test_sigint()
{
  bool test_pass = true;
  ac_int<W,true> input, MAX_VAL;
  MAX_VAL.template set_val<AC_VAL_MAX>();
  bool flag;
  // std::cout << "MAX_VAL = " << MAX_VAL << std::endl;
  for (int i = -MAX_VAL; i <= MAX_VAL; i++) {
    input =i;
    // std::cout << "Input = " << input << std::endl;
    int result_leading_1 = ac_math::leading1<W>(input,flag);
    int result_leading_0 = ac_math::leading0<W>(input,flag);
    int result_trailing_1 = ac_math::trailing1<W>(input,flag);
    int result_trailing_0 = ac_math::trailing0<W>(input,flag);

    int expected_leading_1 = leading1_sint_ref<W>(input);
    int expected_leading_0 = leading0_sint_ref<W>(input);
    int expected_trailing_1 = trailing1_sint_ref<W>(input);
    int expected_trailing_0 = trailing0_sint_ref<W>(input);
	
	// std::cout << " ---------------- " << std::endl;
    // bin<ac_int<W,true>,W>(input);
    // std::cout << "Input: " << input << std::endl;
    // std::cout << "Actual vs Expected Leading 1 index: " << result_leading_1 << " vs " << expected_leading_1 << std::endl;
    // std::cout << "Actual vs Expected Leading 0 index: " << result_leading_0 << " vs " << expected_leading_0 << std::endl;
    // std::cout << "Actual vs Expected Trailing 1 index: " << result_trailing_1 << " vs " << expected_trailing_1 << std::endl;
    // std::cout << "Actual vs Expected Trailing 0 index: " << result_trailing_0 << " vs " << expected_trailing_0 << std::endl;

    if (result_leading_1 != expected_leading_1) {
      test_pass = false;
      std::cout << "Test failed for leading1 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_leading_1 << ", got " << result_leading_1 << std::endl;
	  break;				
    }
    if (result_trailing_1 != expected_trailing_1) {
      test_pass = false;
      std::cout << "Test failed for trailing1 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_trailing_1 << ", got " << result_trailing_1 << std::endl;
	  break;				

    }
    if (result_leading_0 != expected_leading_0) {
      test_pass = false;
      std::cout << "Test failed for leading0 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_leading_0 << ", got " << result_leading_0 << std::endl;
	  break;				

    }
    if (result_trailing_0 != expected_trailing_0) {
      test_pass = false;
      std::cout << "Test failed for trailing0 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_trailing_0 << ", got " << result_trailing_0 << std::endl;
	  break;				
    }
  }
  if (test_pass) {
    std::cout << "All tests for signed integer values passed for W = " << W << std::endl;
  }
}

template<int W, int I>
void test_sigfix()
{
  bool test_pass = true;
  ac_fixed<W,I,true> input, MAX_VAL;
  MAX_VAL.template set_val<AC_VAL_MAX>();
  bool flag;
  // std::cout << "MAX_VAL = " << MAX_VAL << std::endl;
  for (ac_fixed<32,16,true> i = -MAX_VAL; i <= MAX_VAL; i++) {
    input =i;
    // std::cout << "Input = " << input << std::endl;
    int result_leading_1 = ac_math::leading1<W,I>(input,flag);
    int result_leading_0 = ac_math::leading0<W,I>(input,flag);
    int result_trailing_1 = ac_math::trailing1<W,I>(input,flag);
    int result_trailing_0 = ac_math::trailing0<W,I>(input,flag);

    int expected_leading_1 = leading1_sfix_ref<W,I>(input);
    int expected_leading_0 = leading0_sfix_ref<W,I>(input);
    int expected_trailing_1 = trailing1_sfix_ref<W,I>(input);
    int expected_trailing_0 = trailing0_sfix_ref<W,I>(input);
	
	// std::cout << " ---------------- " << std::endl;
    // bin<ac_fixed<W,I,true>,W>(input);
    // std::cout << "Input: " << input << std::endl;
    // std::cout << "Actual vs Expected Leading 1 index: " << result_leading_1 << " vs " << expected_leading_1 << std::endl;
    // std::cout << "Actual vs Expected Leading 0 index: " << result_leading_0 << " vs " << expected_leading_0 << std::endl;
    // std::cout << "Actual vs Expected Trailing 1 index: " << result_trailing_1 << " vs " << expected_trailing_1 << std::endl;
    // std::cout << "Actual vs Expected Trailing 0 index: " << result_trailing_0 << " vs " << expected_trailing_0 << std::endl;
    if (result_leading_1 != expected_leading_1) {
      test_pass = false;
      std::cout << "Test failed for leading1 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_leading_1 << ", got " << result_leading_1 << std::endl;
	  break;				
    }
    if (result_trailing_1 != expected_trailing_1) {
      test_pass = false;
      std::cout << "Test failed for trailing1 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_trailing_1 << ", got " << result_trailing_1 << std::endl;
	  break;				

    }
    if (result_leading_0 != expected_leading_0) {
      test_pass = false;
      std::cout << "Test failed for leading0 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_leading_0 << ", got " << result_leading_0 << std::endl;
	  break;				

    }
    if (result_trailing_0 != expected_trailing_0) {
      test_pass = false;
      std::cout << "Test failed for trailing0 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_trailing_0 << ", got " << result_trailing_0 << std::endl;
	  break;				
    }
  }
  if (test_pass) {
    std::cout << "All tests for signed fixed values passed for W = " << W << std::endl;
  }
}

template<int W, int I>
void test_ufix()
{
  bool test_pass = true;
  ac_fixed<W,I,false> input, MAX_VAL;
  MAX_VAL.template set_val<AC_VAL_MAX>();
  bool flag;
  // std::cout << "MAX_VAL = " << MAX_VAL << std::endl;
  for (ac_fixed<32,16,true> i = 0; i <= MAX_VAL; i++) {
    input =i;
    // std::cout << "Input = " << input << std::endl;
    int result_leading_1 = ac_math::leading1<W,I>(input,flag);
    int result_leading_0 = ac_math::leading0<W,I>(input,flag);
    int result_trailing_1 = ac_math::trailing1<W,I>(input,flag);
    int result_trailing_0 = ac_math::trailing0<W,I>(input,flag);

    int expected_leading_1 = leading1_ufix_ref<W,I>(input);
    int expected_leading_0 = leading0_ufix_ref<W,I>(input);
    int expected_trailing_1 = trailing1_ufix_ref<W,I>(input);
    int expected_trailing_0 = trailing0_ufix_ref<W,I>(input);
	
	// std::cout << " ---------------- " << std::endl;
    // bin<ac_fixed<W,I,true>,W>(input);
    // std::cout << "Input: " << input << std::endl;
    // std::cout << "Actual vs Expected Leading 1 index: " << result_leading_1 << " vs " << expected_leading_1 << std::endl;
    // std::cout << "Actual vs Expected Leading 0 index: " << result_leading_0 << " vs " << expected_leading_0 << std::endl;
    // std::cout << "Actual vs Expected Trailing 1 index: " << result_trailing_1 << " vs " << expected_trailing_1 << std::endl;
    // std::cout << "Actual vs Expected Trailing 0 index: " << result_trailing_0 << " vs " << expected_trailing_0 << std::endl;

    if (result_leading_1 != expected_leading_1) {
      test_pass = false;
      std::cout << "Test failed for leading1 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_leading_1 << ", got " << result_leading_1 << std::endl;
	  break;				
    }
    if (result_trailing_1 != expected_trailing_1) {
      test_pass = false;
      std::cout << "Test failed for trailing1 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_trailing_1 << ", got " << result_trailing_1 << std::endl;
	  break;				

    }
    if (result_leading_0 != expected_leading_0) {
      test_pass = false;
      std::cout << "Test failed for leading0 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_leading_0 << ", got " << result_leading_0 << std::endl;
	  break;				

    }
    if (result_trailing_0 != expected_trailing_0) {
      test_pass = false;
      std::cout << "Test failed for trailing0 with W = " << W << ",Signed input = " << input
                << ". Expected " << expected_trailing_0 << ", got " << result_trailing_0 << std::endl;
	  break;				
    }
  }
  if (test_pass) {
    std::cout << "All tests for unsigned fix values passed for W = " << W << std::endl;
  }
}



/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Main block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

int main(int argc, char *argv[])
{
  test_uint<8>();
  test_uint<9>();
  test_sigint<8>();
  test_sigint<9>();
  test_sigfix<8,3>();
  test_sigfix<9,4>();
  test_ufix<8,3>();
  test_ufix<9,4>();
  
  return 0;
}


