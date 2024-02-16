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
  Reference design to find the position of the leading 1 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int leading1_ref(ac_int<W,false> input) {
    for (int i=W-1; i>=0; i--) {
      if (input[i]==1)
      { return i; }
    }
  }

  /*##########################################
  Reference design to find the position of the trailing 1 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int trailing1_ref(ac_int<W,false> input) {
    for (int i=0; i<=W-1; i++) {
      if (input[i]==1)
      { return i; }
    }
  }

  /*##########################################
  Reference design to find the position of the leading 0 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int leading0_ref(ac_int<W,false> input) {
    for (int i=W-1; i>=0; i--) {
      if (input[i]==0)
      { return i; }
    }
  }

  /*##########################################
  Reference design to find the position of the trailing 0 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int trailing0_ref(ac_int<W,false> input) {
    for (int i=0; i<=W-1; i++) {
      if (input[i]==0)
      { return i; }
    }
  }

//--------------------------------------------------------------------SIGNED-----------------------------------------------------------------//
  /*##########################################
  Reference design to find the position of the leading 1 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int leading1_ref(ac_int<W,true> input) {
    for (int i=W-2; i>=0; i--) {
      if (input[i]==1)
      { return i; }
    }
  }

  /*##########################################
  Reference design to find the position of the trailing 1 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int trailing1_ref(ac_int<W,true> input) {
    for (int i=0; i<=W-2; i++) {
      if (input[i]==1)
      { return i; }
    }
  }

  /*##########################################
  Reference design to find the position of the leading 0 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int leading0_ref(ac_int<W,true> input) {
    for (int i=W-2; i>=0; i--) {
      if (input[i]==0)
      { return i; }
    }
  }

  /*##########################################
  Reference design to find the position of the trailing 0 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int trailing0_ref(ac_int<W,true> input) {
    for (int i=0; i<=W-2; i++) {
      if (input[i]==0)
      { return i; }
    }
  }

//********************************************************************FIXED******************************************************************//
//--------------------------------------------------------------------UNSIGNED-----------------------------------------------------------------//
  /*##########################################
  Reference design to find the position of the leading 1 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int leading1_ref(ac_fixed<W, I, false> input) {
    for (int i=W-1; i>=0; i--) {
      if (input[i]==1)
      { return i; }
    }
  }

  /*##########################################
  Reference design to find the position of the trailing 1 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int trailing1_ref(ac_fixed<W, I, false> input) {
    for (int i=0; i<=W-1; i++) {
      if (input[i]==1)
      { return i; }
    }
  }

  /*##########################################
  Reference design to find the position of the leading 0 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int leading0_ref(ac_fixed<W, I, false> input) {
    for (int i=W-1; i>=0; i--) {
      if (input[i]==0)
      { return i; }
    }
  }

  /*##########################################
  Reference design to find the position of the trailing 0 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int trailing0_ref(ac_fixed<W, I, false> input) {
    for (int i=0; i<=W-1; i++) {
      if (input[i]==0)
      { return i; }
    }
  }

//--------------------------------------------------------------------SIGNED-----------------------------------------------------------------//
  /*##########################################
  Reference design to find the position of the leading 1 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int leading1_ref(ac_fixed<W, I, true> input) {
    for (int i=W-2; i>=0; i--) {
      if (input[i]==1)
      { return i; }
    }
  }

  /*##########################################
  Reference design to find the position of the trailing 1 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int trailing1_ref(ac_fixed<W, I, true> input) {
    for (int i=0; i<=W-2; i++) {
      if (input[i]==1)
      { return i; }
    }
  }

  /*##########################################
  Reference design to find the position of the leading 0 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int leading0_ref(ac_fixed<W, I, true> input) {
    for (int i=W-2; i>=0; i--) {
      if (input[i]==0)
      { return i; }
    }
  }

  /*##########################################
  Reference design to find the position of the trailing 0 bit from the LSB
  ###########################################*/
  template<int W, int I>
  unsigned int trailing0_ref(ac_fixed<W, I, true> input) {
    for (int i=0; i<=W-2; i++) {
      if (input[i]==0)
      { return i; }
    }
  }
  

template <class input_type, class output_type, unsigned width, unsigned i_width>
 int driver_function(input_type input) {
	
	bool test_fail=0;
	bool flag;
	/*##############################
	LEADING_1
	###############################*/
    if (leading1_ref<width,i_width>(input)==ac_math::leading1<width>(input, flag)) { // leading one
      cout << "leading one  Results Match" << endl;
    } else {
      cout << "leading one  Results DONOT Match" << endl;
	  test_fail=1;
    }
    cout << "Design    Value: Leading 1 bit position from the LSB side " << ac_math::leading1<width>(input, flag) << endl;
    cout << "Reference Value: Leading 1 bit position from the LSB side " << leading1_ref<width,i_width>(input) << endl;

	/*##############################
	TRAILING_1
	###############################*/
    if (trailing1_ref<width,i_width>(input)==ac_math::trailing1<width>(input, flag)) { // leading one
      cout << "Trailing one  Results Match" << endl;
    } else {
      cout << "Trailing one  Results DONOT Match" << endl;
  	  test_fail=1;
    }
    cout << "Design    Value: Trailing 1 bit position from the LSB side " << ac_math::trailing1<width>(input, flag) << endl;
    cout << "Reference Value: Trailing 1 bit position from the LSB side " << trailing1_ref<width,i_width>(input) << endl;
	
	
	/*##############################
	LEADING_0
	###############################*/
    if (leading0_ref<width,i_width>(input)==ac_math::leading0<width>(input, flag)) { // leading one
      cout << "leading one  Results Match" << endl;
    } else {
      cout << "leading one  Results DONOT Match" << endl;
  	  test_fail=1;
    }
    cout << "Design    Value: Leading 0 bit position from the LSB side " << ac_math::leading0<width>(input, flag) << endl;
    cout << "Reference Value: Leading 0 bit position from the LSB side " << leading0_ref<width,i_width>(input) << endl;
	
	
	/*##############################
	TRAILING_0
	###############################*/
    if (trailing0_ref<width,i_width>(input)==ac_math::trailing0<width>(input, flag)) { // leading one
      cout << "Trailing 0  Results Match" << endl;
    } else {
      cout << "Trailing 0  Results DONOT Match" << endl;
	  test_fail=1;
	}
    cout << "Design    Value: Trailing 0 bit position from the LSB side " << ac_math::trailing0<width>(input, flag) << endl;
    cout << "Reference Value: Trailing 0 bit position from the LSB side " << trailing0_ref<width,i_width>(input) << endl;
   
  if (test_fail) {
    cout << "  ac_leading - FAILED - Mismatch found" << endl; 
    cout << "=============================================================================" << endl; 
    return -1; 
  } else {
    cout << "  ac_leading - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }
  
  return 0;
  }
  
  
  
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Main block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

int main(int argc, char *argv[])
{
  enum {
    DATA_WIDTH = 32,
    INTEGER_WIDTH = 8,
	TEST_NUMBERS = 10,
  };
  
  typedef ac_int<ac::nbits<DATA_WIDTH-1>::val,false> OUT_TYP;
  OUT_TYP out;
  bool status;
  
  typedef ac_int<DATA_WIDTH,false> IN_TYP1;
  IN_TYP1 in1;
  cout << "===================INT_ UNSIGNED================================" << endl;
  for (int i=0; i<=TEST_NUMBERS; i++)
  {
  throw_dice<IN_TYP1>(in1);
  bin<IN_TYP1,DATA_WIDTH>(in1);
  driver_function<IN_TYP1, OUT_TYP, DATA_WIDTH, INTEGER_WIDTH>(in1);
  cout << "---------------------------" << endl;

  }

  cout << "===================INT_ NSIGNED================================" << endl;
  typedef ac_int<DATA_WIDTH,true> IN_TYP2;
  IN_TYP2 in2;
  for (int i=0; i<=TEST_NUMBERS; i++)
  {
  throw_dice<IN_TYP2>(in2);
  bin<IN_TYP2,DATA_WIDTH>(in2);
  driver_function<IN_TYP2, OUT_TYP, DATA_WIDTH, INTEGER_WIDTH>(in1);
  cout << "---------------------------" << endl;
  }

  cout << "===================FIX_UNSIGNED================================" << endl;

  typedef ac_fixed<DATA_WIDTH,INTEGER_WIDTH,false> IN_TYP3;
  IN_TYP3 in3;
  for (int i=0; i<=TEST_NUMBERS; i++)
  {
  throw_dice<IN_TYP3>(in3);
  bin<IN_TYP3,DATA_WIDTH>(in3);
  driver_function<IN_TYP3, OUT_TYP, DATA_WIDTH, INTEGER_WIDTH>(in1);
  cout << "---------------------------" << endl;
  }

  cout << "===================FIX_SIGNED================================" << endl;

  typedef ac_fixed<DATA_WIDTH,INTEGER_WIDTH,true> IN_TYP4;
  IN_TYP4 in4;
  for (int i=0; i<=TEST_NUMBERS; i++)
  {
  throw_dice<IN_TYP4>(in4);
  bin<IN_TYP4,DATA_WIDTH>(in4);
  driver_function<IN_TYP4, OUT_TYP, DATA_WIDTH, INTEGER_WIDTH>(in1);
  cout << "---------------------------" << endl;
  }

  return 0;
}


