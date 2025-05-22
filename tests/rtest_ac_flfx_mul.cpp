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
// ac_flfx_mul() function using a variety of bitwidths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include -O3 rtest_ac_flfx_mul.cpp -o design
//   ./design
//
// Using the -O3 flag isn't necessary but the tests will take a lot longer
// without it.

#include <ac_math/ac_flfx_mul.h>

template <bool subn_support, class out_type, class fx_in_type, class fl_in_type>
out_type test_ac_flfx_mul (
  const fx_in_type &fx_in,
  const fl_in_type &fl_in
)
{
  out_type out = 0.0;
  ac_math::ac_flfx_mul<subn_support>(fx_in, fl_in, out);
  return out;
}

template <class out_type, class fx_in_type, class fl_in_type>
out_type calc_ref(const fx_in_type &fx_in, const fl_in_type &fl_in)
{
  auto out_ref_full = fx_in*fl_in; // Reference output with full precision.
  // Construct reference output with the mantissa and exponent of the full precision output. Make
  // sure it is normalized, as the * operator for ac_fixed*ac_float multiplication doesn't give a
  // normalized output by default.
  out_type out_ref(out_ref_full.mantissa(), out_ref_full.exp(), true);
  return out_ref;
}

template <class fl_in_type>
bool out_zero_checking(const fl_in_type &out)
{
  // Both the output mantissa and exponent should be zero. If they are not, this
  // function returns true.
  // Zero checking for floating point values generally only means checking
  // the mantissa but this test is a little more thorough.
  return (out.mantissa() != 0 || out.exp() != 0);
}

template <bool subn_support, class fl_in_type, class out_type>
bool check_output(const fl_in_type &fl_in, const out_type &out, const out_type &out_ref)
{
  if (out_ref == 0) {
    if (out_zero_checking(out)) {
      return true;
    }
  } else if (fl_in.is_subn() || out_ref.is_subn()) {
    if (subn_support) {
      if (out != out_ref || out.is_norm() != out_ref.is_norm()) {
        return true;
      }
    } else if (out_zero_checking(out)) {
      return true;
    }
  } else if (out != out_ref || !out.is_norm()) {
    return true;
  }

  return false;
}

template <class fx_in_type, class fl_in_type, class out_type>
void debug_statements(
  const bool &test_fail,
  const fx_in_type &fx_in, const fl_in_type &fl_in,
  const out_type &out, const out_type &out_ref
)
{
  #ifdef DEBUG
  if (test_fail) {
    std::cout << "fx_in   = " << fx_in << std::endl;
    std::cout << "fl_in   = " << fl_in << std::endl;
    std::cout << "out     = " << out << std::endl;
    std::cout << "out_ref = " << out_ref << std::endl;
    AC_ASSERT(false, "Test FAILED!");
  }
  #endif
}

template <bool subn_support, class out_type, class fx_in_type, class fl_in_type>
bool calc_and_check_outputs(const fx_in_type &fx_in, const fl_in_type &fl_in)
{
  out_type out = test_ac_flfx_mul<subn_support, out_type>(fx_in, fl_in);
  out_type out_ref = calc_ref<out_type>(fx_in, fl_in);
  
  bool output_mismatch = check_output<subn_support>(fl_in, out, out_ref);
  debug_statements(output_mismatch, fx_in, fl_in, out, out_ref); // Print debug statements if needed.
  return output_mismatch;
}

#include <vector>

template <
  bool subn_support,
  int Wfx, int Ifx, bool Sfx, int Wfl, int Ifl, int Efl,
  int outW, int outI, int outE>
bool test_driver(bool details = false)
{
  typedef ac_fixed<Wfx, Ifx, Sfx> fx_in_type;
  typedef ac_float<Wfl, Ifl, Efl> fl_in_type;
  typedef ac_float<outW, outI, outE> out_type;

  std::cout << "TEST: ac_flfx_mul() ";
  std::cout << "subn_support = " << subn_support;
  std::cout << " INPUTS: ";
  std::cout.width(39);
  std::cout << std::left << fx_in_type::type_name();
  std::cout.width(31);
  std::cout << std::left << fl_in_type::type_name();
  std::cout << "OUTPUT: ";
  std::cout.width(32);
  std::cout << std::left << out_type::type_name();
  std::cout << "RESULT: ";

  ac_fixed<Wfx, Ifx, Sfx> fx_min = 0.0, fx_max = 0.0, fx_step = 0.0;
  // Number of LSBs to set to 1 for the step. If only the last LSB is set to 1, the testing is exhaustive,
  // but if Wfx > 10, (Wfx - 9) LSBs will be set to 1 and the testing is not exhaustive, to save on time.
  constexpr int num_fx_1_lsbs = AC_MAX(Wfx - 10 + 1, 1);
  fx_step.set_slc(0, ac_int<num_fx_1_lsbs, false>(-1));
  fx_min.template set_val<AC_VAL_MIN>();
  fx_max.template set_val<AC_VAL_MAX>();

  ac_fixed<Wfl, Ifl, true> mant_min = 0.0, mant_max = 0.0, mant_step = 0.0;
  // The number of LSBs set to 1 for the mantissa testing is calculated in a manner similar to the fixed
  // point testing above.
  constexpr int num_mant_1_lsbs = AC_MAX(Wfl - 14 + 1, 1);
  mant_step.set_slc(0, ac_int<num_mant_1_lsbs, false>(-1));
  mant_min.template set_val<AC_VAL_MIN>();
  mant_max.template set_val<AC_VAL_MAX>();

  if (details) {
    std::cout << "Testing limits:" << std::endl;
    std::cout << "fx_min    = " << fx_min << std::endl;
    std::cout << "fx_max    = " << fx_max << std::endl;
    std::cout << "fx_step   = " << fx_step << std::endl;
    std::cout << "mant_min  = " << mant_min << std::endl;
    std::cout << "mant_max  = " << mant_max << std::endl;
    std::cout << "mant_step = " << mant_step << std::endl;
  }

  constexpr int IMIN_EXP = fl_in_type::MIN_EXP;
  constexpr int IMAX_EXP = fl_in_type::MAX_EXP;
  // Vector of floating point inputs to test with.
  std::vector<ac_int<Efl, true> > infl_exp_vec = {IMIN_EXP, IMIN_EXP/2, 0, IMAX_EXP/2, IMAX_EXP};

  bool all_outputs_correct = true;
  
  for (ac_int<Efl, true> infl_exp : infl_exp_vec) {
    for (ac_fixed<Wfl + 1, Ifl + 1, true> mant_it = mant_min; mant_it <= mant_max; mant_it += mant_step) {
      ac_fixed<Wfl, Ifl, true> infl_mant = mant_it;
      // Construct input float. Bear in mind that the mantissa value isn't already normalized. As a
      // result, both mantissa and exponent values may change after normalization.
      fl_in_type fl_in(mant_it, infl_exp, true);
      
      for (ac_fixed<Wfx + 1, Ifx + 1, Sfx> fx_it = fx_min; fx_it <= fx_max; fx_it += fx_step) {
        fx_in_type fx_in = fx_it;
        bool output_mismatch = calc_and_check_outputs<subn_support, out_type>(fx_in, fl_in);
        if (output_mismatch) { all_outputs_correct = false; }
      }
    }
  }

  if (all_outputs_correct) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
  }

  return all_outputs_correct;
}

int main(int, char **)
{
  std::cout << "=============================================================================" << std::endl;
  std::cout << "Testing function: ac_flfx_mul()" << std::endl;
  
  bool all_tests_pass = true;

  all_tests_pass = test_driver<true, 10,  2, true, 14,  2, 6, 14,  2,  6>() && all_tests_pass;
  all_tests_pass = test_driver<true, 10,  2, true, 14,  2, 6, 24,  4,  6>() && all_tests_pass;
  all_tests_pass = test_driver<true, 10,  2, true, 14,  2, 6, 10, 12,  6>() && all_tests_pass;
  all_tests_pass = test_driver<true, 10,  2, true, 14,  2, 6, 10, -4,  6>() && all_tests_pass;
  all_tests_pass = test_driver<true, 10,  2, true, 14,  2, 6, 10, 12, 10>() && all_tests_pass;
  all_tests_pass = test_driver<true, 10,  2, true, 14,  2, 6, 10, -4, 10>() && all_tests_pass;

  std::cout << std::endl;
  
  all_tests_pass = test_driver<true, 10,  2, false, 14,  2, 6, 24,  4,  6>() && all_tests_pass;
  all_tests_pass = test_driver<true, 10,  2, false, 14,  2, 6, 10, 12,  6>() && all_tests_pass;
  all_tests_pass = test_driver<true, 10,  2, false, 14,  2, 6, 10, -4,  6>() && all_tests_pass;

  std::cout << std::endl;

  all_tests_pass = test_driver<true, 50,  2, true, 54,  2, 8, 54,  2, 10>() && all_tests_pass;
  all_tests_pass = test_driver<true, 10, -2, true, 12, -2, 8, 54,  2, 10>() && all_tests_pass;
  all_tests_pass = test_driver<true, 10, 15, true, 12, 17, 8, 54,  2, 10>() && all_tests_pass;

  std::cout << std::endl;

  all_tests_pass = test_driver<false, 10,  2, true, 14,  2, 6, 14,  2,  6>() && all_tests_pass;
  all_tests_pass = test_driver<false, 10,  2, true, 14,  2, 6, 24,  4,  6>() && all_tests_pass;
  all_tests_pass = test_driver<false, 10,  2, true, 14,  2, 6, 10, 12,  6>() && all_tests_pass;
  all_tests_pass = test_driver<false, 10,  2, true, 14,  2, 6, 10, -4,  6>() && all_tests_pass;
  all_tests_pass = test_driver<false, 10,  2, true, 14,  2, 6, 10, 12, 10>() && all_tests_pass;
  all_tests_pass = test_driver<false, 10,  2, true, 14,  2, 6, 10, -4, 10>() && all_tests_pass;

  std::cout << std::endl;
  
  all_tests_pass = test_driver<false, 10,  2, false, 14,  2, 6, 24,  4,  6>() && all_tests_pass;
  all_tests_pass = test_driver<false, 10,  2, false, 14,  2, 6, 10, 12,  6>() && all_tests_pass;
  all_tests_pass = test_driver<false, 10,  2, false, 14,  2, 6, 10, -4,  6>() && all_tests_pass;

  std::cout << std::endl;

  all_tests_pass = test_driver<false, 50,  2, true, 54,  2, 8, 54,  2, 10>() && all_tests_pass;
  all_tests_pass = test_driver<false, 10, -2, true, 12, -2, 8, 54,  2, 10>() && all_tests_pass;
  all_tests_pass = test_driver<false, 10, 15, true, 12, 17, 8, 54,  2, 10>() && all_tests_pass;
  
  std::cout << "=============================================================================" << std::endl;
  std::cout << "  Testbench finished." << std::endl;
  
  // Notify the user whether or not the test was a failure.
  if (!all_tests_pass) {
    std::cout << "  ac_flfx_mul - FAILED - Output not correct for all test values" << std::endl;
    std::cout << "=============================================================================" << std::endl;
    return -1;
  } else {
    std::cout << "  ac_flfx_mul - PASSED" << std::endl;
    std::cout << "=============================================================================" << std::endl;
  }

  return 0;
}
