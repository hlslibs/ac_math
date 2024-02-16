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
 *  Copyright  Siemens                                                *
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
#include <ac_float_add_tree.h>
#include <random>
#include <ac_math/ac_random.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

using namespace ac_math;

template<int W, int I, int E, ac_q_mode Q>
void random_flt_gen(ac_float<W,I,E,Q> &x) {
  ac_fixed<W, I, true> x_mant = 0.0;
  ac_random(x_mant);
  ac_int<E, true> x_exp = 0;
  ac_random(x_exp);
  x.set_mantissa(x_mant);
  x.set_exp(x_exp);
  x.normalize();
}

#ifdef DEBUG
template <int N, class T>
void print_arr_to_file(T a[N]) {
  std::ofstream my_file;
  my_file.open("error_causing_input.txt");

  for (int i = 0; i < N; i++) {
    my_file << a[i] << std::endl;
  }
  
  my_file.close();
}
#endif

enum out_val_format {
  AC_FX_OUTPUT,
  AC_FL_OUTPUT,
  AC_OTHER_OUTPUT,
};

// Specifies what format the input vectors follow. Description of each format is as follows:
// ALL_MAX_VALS : All values are maximum input values.
// ALL_MIN_VALS : All values are minimum input values.
// ALL_QUANT_VALS : All values correspond to quantum input values.
// ALL_QUANT_FX_VALS : All values correspond to quantum output value, if using fixed point outputs.
// RANDOM_ZERO_VALS :
//   If using floating point inputs/outputs, the input vector is randomly set to either a zero value
//   or the quantum input value, with a 50% chance of either.
//   If using floating point inputs+fixed point outputs, the input vector is randomly set to either a
//   zero value or the quantum output value, with a 50% chance of either.
// UNCONSTRAINED_RANDOM_VALS : Random unconstrained values.
// CONSTRAINED_RANDOM_VALS : values constrained to not trigger saturation in output, if using fixed
//   point outputs.

enum test_vector_format {
  ALL_MAX_VALS,
  ALL_MIN_VALS,
  ALL_QUANT_VALS,
  ALL_QUANT_FX_VALS,
  RANDOM_ZERO_VALS,
  UNCONSTRAINED_RANDOM_VALS,
  CONSTRAINED_RANDOM_VALS,
};

template <int N, class T>
struct helper_st {
  static const out_val_format outvf = AC_OTHER_OUTPUT;
};

template <int N, int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
struct helper_st<N, ac_fixed<W, I, S, Q, O> > {
  static const out_val_format outvf = AC_FX_OUTPUT;
  
  static void test_for_invld_formats(const test_vector_format &tvf) {
    if (O == AC_WRAP) {
      AC_ASSERT(tvf != ALL_MAX_VALS && tvf != ALL_MIN_VALS && tvf != UNCONSTRAINED_RANDOM_VALS, "tvf set to incorrect value for non-saturating output type.");
    }
  }
  
  static void gen_test_vals(
    ac_fixed<W, I, S, Q, O> &quant_out_val,
    double &min_rand_limit, double &max_rand_limit
  ) {
    typedef ac_fixed<W, I, S, Q, O> T_fx;
    // Truncation mode is chosen to make sure the constrained random values are
    // constrained correctly.
    typedef ac_fixed<T_fx::width, T_fx::i_width, T_fx::sign, AC_TRN> T_fx_trn_wrap;
    T_fx_trn_wrap max_fx_val = 0.0, min_fx_val = 0.0;
    min_fx_val.template set_val<AC_VAL_MIN>();
    max_fx_val.template set_val<AC_VAL_MAX>();
    
    quant_out_val.template set_val<AC_VAL_QUANTUM>();
    min_rand_limit = (min_fx_val/N).to_double();
    max_rand_limit = (max_fx_val/N).to_double();
  }
};

template <int N, int W, int I, int E, ac_q_mode Q>
struct helper_st<N, ac_float<W, I, E, Q> > {
  static const out_val_format outvf = AC_FL_OUTPUT;
  
  static void test_for_invld_formats(const test_vector_format &tvf) {
    AC_ASSERT(tvf != ALL_QUANT_FX_VALS && tvf != CONSTRAINED_RANDOM_VALS, "tvf set to incorrect value.");
  }
  
  static void gen_test_vals(
    ac_float<W, I, E, Q> &quant_out_val,
    double &min_rand_limit, double &max_rand_limit
  ) {
    // All of these values will be ignored for floating point testing.
    quant_out_val = 0.0;
    min_rand_limit = 0.0;
    max_rand_limit = 0.0;
  }
};

template<typename T_out, int N, int W, int I, int E, ac_q_mode Q>
void gen_test_vectors(ac_float<W,I,E,Q> x[N], test_vector_format tvf)
{
  typedef helper_st<N, T_out> hst_type;
  
  static_assert(hst_type::outvf != AC_OTHER_OUTPUT, "Invalid output type.");
  hst_type::test_for_invld_formats(tvf);
  
  // quant_fx_val, min_rand_limit and max_rand_limit will all be set to 0 and ignored if the output
  // type is ac_float.
  T_out quant_fx_val;
  double min_rand_limit, max_rand_limit;
  hst_type::gen_test_vals(quant_fx_val, min_rand_limit, max_rand_limit);
  
  static std::default_random_engine generator;
  static std::uniform_real_distribution<double> distribution(min_rand_limit, max_rand_limit);
  
  ac_float<W,I,E,Q> quant_in_val;
  quant_in_val.template set_val<AC_VAL_QUANTUM>();
  double quant_val_zero_testing = (hst_type::outvf == AC_FX_OUTPUT ? quant_fx_val.to_double() : quant_in_val.to_double());

  for (int i = 0; i < N; i++) {
    switch (tvf) {
      case ALL_MAX_VALS:
        x[i].template set_val<AC_VAL_MAX>();
        break;
      case ALL_MIN_VALS:
        x[i].template set_val<AC_VAL_MIN>();
        break;
      case ALL_QUANT_VALS:
        x[i].template set_val<AC_VAL_QUANTUM>();
        break;
      case ALL_QUANT_FX_VALS:
        x[i] = quant_fx_val;
        break;
      case RANDOM_ZERO_VALS:
        x[i] = (rand()%2 == 0) ? 0.0 : quant_val_zero_testing;
        break;
      case CONSTRAINED_RANDOM_VALS:
        x[i] = distribution(generator);
        break;
      default: // UNCONSTRAINED_RANDOM_VALS
        random_flt_gen(x[i]);
        break;
    }
  }
}

template <typename T_out, int N, typename T>
double gen_ref_val(T a[N]) {
  double acc_ref = 0.0;
  for (int i=0; i<N; i++) {
    acc_ref += a[i].to_double();
  }
  acc_ref = (T_out(acc_ref)).to_double();
  
  return acc_ref;
}


template <int N, typename T_out, typename T_in>
void err_calc(
  const T_in (&a)[N],
  const T_out &acc_at, const T_out &acc_bat, const double &acc_ref,
  const double &allowed_error_at, const double &allowed_error_bat,
  double &max_error_at, double &max_error_bat
)
{
  double this_error_at = 0.0, this_error_bat = 0.0;

  if (acc_ref != 0) {
    this_error_at = (acc_at.to_double() - acc_ref)/acc_ref;
    this_error_bat = (acc_bat.to_double() - acc_ref)/acc_ref;
  } else {
    this_error_at = acc_at.to_double() - acc_ref;
    this_error_bat = acc_bat.to_double() - acc_ref;
  }

  this_error_at *= 100.0;
  if (this_error_at < 0) {
    this_error_at = -this_error_at;
  }
  if (this_error_at > max_error_at) {
    max_error_at = this_error_at;
  }

  this_error_bat *= 100.0;
  if (this_error_bat < 0) {
    this_error_bat = -this_error_bat;
  }
  if (this_error_bat > max_error_bat) {
    max_error_bat = this_error_bat;
  }

  #ifdef DEBUG
  bool this_tolerance_exceeded = (this_error_at > allowed_error_at) || (this_error_bat > allowed_error_bat);

  if (this_tolerance_exceeded) {
    print_arr_to_file<N>(a);
    std::cout << "acc_at            = " << acc_at << std::endl;
    std::cout << "acc_bat           = " << acc_bat << std::endl;
    std::cout << "acc_ref           = " << acc_ref << std::endl;
    std::cout << "this_error_at     = " << this_error_at << std::endl;
    std::cout << "this_error_bat    = " << this_error_bat << std::endl;
    std::cout << "allowed_error_at  = " << allowed_error_at << std::endl;
    std::cout << "allowed_error_bat = " << allowed_error_bat << std::endl;
    AC_ASSERT(false, "Test failed!");
  }
  #endif
}

template<int N, int AccExtraLSBs, typename elem_flt, typename out_fx>
void test_fx(double &cumul_max_error_at, double &cumul_max_error_bat, const double &allowed_error_at, const double &allowed_error_bat) {
  std::cout << "test_fx.  N = ";
  std::cout.width(5);
  std::cout << std::left << N;
  
  std::cout << "ELSBs = ";
  std::cout.width(5);
  std::cout << std::left << AccExtraLSBs;
  
  std::cout << "elem_flt : ";
  std::cout.width(30);
  std::cout << std::left << elem_flt::type_name();
  
  std::cout << "out_fx : ";
  std::cout.width(37);
  std::cout << std::left << out_fx::type_name();
  
  // Vector of special test vector formats.
  std::vector<test_vector_format> sp_tvfs;
  
  sp_tvfs.push_back(ALL_QUANT_VALS);
  sp_tvfs.push_back(ALL_QUANT_FX_VALS);
  
  if (out_fx::o_mode != AC_WRAP) {
    sp_tvfs.push_back(ALL_MAX_VALS);
    sp_tvfs.push_back(ALL_MIN_VALS);
  }
  
  double max_error_at = 0.0, max_error_bat = 0.0;

  for (test_vector_format tvf : sp_tvfs) {
    elem_flt a[N];
    gen_test_vectors<out_fx, N>(a, tvf);
    
    out_fx acc_at = 0.0, acc_bat = 0.0;
    // add_tree output
    add_tree_ptr<N, AccExtraLSBs>(a, acc_at);
    // block_add_tree output
    block_add_tree_ptr<N, AccExtraLSBs>(a, acc_bat);

    double acc_ref = gen_ref_val<out_fx, N>(a);
    
    err_calc(a, acc_at, acc_bat, acc_ref, allowed_error_at, allowed_error_bat, max_error_at, max_error_bat);
  }

  // Vector of random test vector formats.
  std::vector<test_vector_format> rn_tvfs;
  
  rn_tvfs.push_back(RANDOM_ZERO_VALS);
  rn_tvfs.push_back(CONSTRAINED_RANDOM_VALS);
  
  if (out_fx::o_mode != AC_WRAP) {
    rn_tvfs.push_back(UNCONSTRAINED_RANDOM_VALS);
  }
  
  constexpr int num_rn_test_vectors = 500;
  
  for (test_vector_format tvf : rn_tvfs) {
    for (int k = 0; k < num_rn_test_vectors; k++) {
      elem_flt a[N];
      gen_test_vectors<out_fx, N>(a, tvf);
    
      out_fx acc_at = 0.0, acc_bat = 0.0;
      // add_tree output
      add_tree_ptr<N, AccExtraLSBs>(a, acc_at);
      // block_add_tree output
      block_add_tree_ptr<N, AccExtraLSBs>(a, acc_bat);

      double acc_ref = gen_ref_val<out_fx, N>(a);
      
      err_calc(a, acc_at, acc_bat, acc_ref, allowed_error_at, allowed_error_bat, max_error_at, max_error_bat);
    }
  }
  
  if (max_error_at > cumul_max_error_at) {
    cumul_max_error_at = max_error_at;
  }
  
  if (max_error_bat > cumul_max_error_bat) {
    cumul_max_error_bat = max_error_bat;
  }
  
  std::cout << " add_tree error = ";
  std::cout.width(12);
  std::cout << max_error_at;
  std::cout << " block_add_tree error = ";
  std::cout.width(12);
  std::cout << max_error_bat;
  
  if (max_error_at > allowed_error_at || max_error_bat > allowed_error_bat) {
    std::cout << " FAILED." << std::endl;
  } else {
    std::cout << " PASSED." << std::endl;
  }
}

template<int N, int AccExtraLSBs, typename elem_flt>
void test_fl(double &cumul_max_error_at, double &cumul_max_error_bat, const double &allowed_error_at, const double &allowed_error_bat) {
  std::cout << "test_fl.  N = ";
  std::cout.width(5);
  std::cout << std::left << N;
  
  std::cout << "ELSBs = ";
  std::cout.width(5);
  std::cout << std::left << AccExtraLSBs;
  
  std::cout << "elem_flt : ";
  std::cout.width(30);
  std::cout << std::left << elem_flt::type_name();
  
  constexpr int acc_W = elem_flt::width;
  constexpr int acc_I = elem_flt::i_width;
  constexpr ac_q_mode acc_Q = elem_flt::q_mode;
  
  constexpr int elem_E = elem_flt::e_width;
  // acc_E_calc_val = intermediate value necessary to calculate acc_E.
  constexpr int acc_E_calc_val = (1 << (elem_E - 1)) + ac::log2_ceil<N>::val;
  // acc_E should have enough bits to store the maximum possible exponent value after accumulation.
  constexpr int acc_E = ac::log2_ceil<acc_E_calc_val>::val + 1;
  
  typedef ac_float<acc_W, acc_I, acc_E, acc_Q> acc_flt;
  
  std::cout << "acc_flt : ";
  std::cout.width(30);
  std::cout << std::left << acc_flt::type_name();
  
  // Vector of special test vector formats.
  std::vector<test_vector_format> sp_tvfs;
  
  sp_tvfs.push_back(ALL_QUANT_VALS);
  sp_tvfs.push_back(ALL_MAX_VALS);
  sp_tvfs.push_back(ALL_MIN_VALS);
  
  double max_error_at = 0.0, max_error_bat = 0.0;

  for (test_vector_format tvf : sp_tvfs) {
    elem_flt a[N];
    gen_test_vectors<acc_flt, N>(a, tvf);
    
    acc_flt acc_at = 0.0, acc_bat = 0.0;
    // add_tree output
    add_tree_ptr<N, AccExtraLSBs>(a, acc_at);
    // block_add_tree output
    block_add_tree_ptr<N, AccExtraLSBs>(a, acc_bat);

    double acc_ref = gen_ref_val<acc_flt, N>(a);
    
    err_calc(a, acc_at, acc_bat, acc_ref, allowed_error_at, allowed_error_bat, max_error_at, max_error_bat);
  }

  // Vector of random test vector formats.
  std::vector<test_vector_format> rn_tvfs;
  
  rn_tvfs.push_back(RANDOM_ZERO_VALS);
  rn_tvfs.push_back(UNCONSTRAINED_RANDOM_VALS);
  
  constexpr int num_rn_test_vectors = 500;
  
  for (test_vector_format tvf : rn_tvfs) {
    for (int k = 0; k < num_rn_test_vectors; k++) {
      elem_flt a[N];
      gen_test_vectors<acc_flt, N>(a, tvf);
    
      acc_flt acc_at = 0.0, acc_bat = 0.0;
      // add_tree output
      add_tree_ptr<N, AccExtraLSBs>(a, acc_at);
      // block_add_tree output
      block_add_tree_ptr<N, AccExtraLSBs>(a, acc_bat);

      double acc_ref = gen_ref_val<acc_flt, N>(a);
      
      err_calc(a, acc_at, acc_bat, acc_ref, allowed_error_at, allowed_error_bat, max_error_at, max_error_bat);
    }
  }
  
  if (max_error_at > cumul_max_error_at) {
    cumul_max_error_at = max_error_at;
  }
  
  if (max_error_bat > cumul_max_error_bat) {
    cumul_max_error_bat = max_error_bat;
  }
  
  std::cout << "add_tree error = ";
  std::cout.width(12);
  std::cout << max_error_at;
  std::cout << " block_add_tree error = ";
  std::cout.width(12);
  std::cout << max_error_bat;
  
  if (max_error_at > allowed_error_at || max_error_bat > allowed_error_bat) {
    std::cout << " FAILED." << std::endl;
  } else {
    std::cout << " PASSED." << std::endl;
  }
}

int main(int argc, char *argv[]) {
  std::cout << "===========================================================================" << std::endl;
  std::cout << "------------------- Running rtest_ac_float_add_tree.cpp -------------------" << std::endl;
  std::cout << "===========================================================================" << std::endl;
  
  double max_error_at = 0, max_error_bat = 0, allowed_error_at = 0.05, allowed_error_bat = 0.5;
  
  test_fx< 16,   0, ac_float<25, 2, 8>,  ac_fixed<64, 32, true, AC_TRN,  AC_SAT> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fx< 16,   8, ac_float<25, 2, 8>,  ac_fixed<64, 32, true, AC_TRN,  AC_SAT> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fx< 16,  32, ac_float<25, 2, 8>,  ac_fixed<64, 32, true, AC_TRN,  AC_SAT> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  std::cout << std::endl;
  
  test_fx<128,   0, ac_float<25, 2, 8>,  ac_fixed<64, 32, true, AC_TRN,  AC_SAT> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fx<128,   8, ac_float<25, 2, 8>,  ac_fixed<64, 32, true, AC_TRN,  AC_SAT> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fx<128,  32, ac_float<25, 2, 8>,  ac_fixed<64, 32, true, AC_TRN,  AC_SAT> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  std::cout << std::endl;
  
  test_fx< 16,   0, ac_float<25, 2, 8>,  ac_fixed<64, 32, true, AC_TRN, AC_WRAP> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fx< 16,   0, ac_float<25, 2, 8>, ac_fixed<64, 32, false, AC_TRN, AC_WRAP> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fx< 16,   0, ac_float<25, 2, 8>, ac_fixed<64, 32, false, AC_TRN,  AC_SAT> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  std::cout << std::endl;
  
  test_fl< 16,   0, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl< 16,   8, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl< 16,  16, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl< 16,  32, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  
  std::cout << std::endl;
  
  test_fl< 32,   0, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl< 32,   8, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl< 32,  16, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl< 32,  32, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  
  std::cout << std::endl;
  
  test_fl< 64,   0, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl< 64,   8, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl< 64,  16, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl< 64,  32, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl< 64,  64, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  
  std::cout << std::endl;
  
  test_fl<128,   0, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl<128,   8, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl<128,  16, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl<128,  32, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl<128,  64, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  test_fl<128, 128, ac_float<25, 2, 8> >(max_error_at, max_error_bat, allowed_error_at, allowed_error_bat);
  
  std::cout << std::endl;
  
  std::cout << "max_error_at      = " << max_error_at << std::endl;
  std::cout << "allowed_error_at  = " << allowed_error_at << std::endl;
  std::cout << "max_error_bat     = " << max_error_bat << std::endl;
  std::cout << "allowed_error_bat = " << allowed_error_bat << std::endl;
  
  if (max_error_at > allowed_error_at || max_error_bat > allowed_error_bat) {
    std::cout << "Test FAILED." << std::endl;
    return -1;
  }
  
  std::cout << "Test PASSED." << std::endl;
  
  return 0;
}
