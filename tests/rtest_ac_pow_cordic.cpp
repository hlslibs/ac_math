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
// ac_pow_cordic() function using a variety of data types and bit-
// widths.

// To compile standalone and run:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_pow_cordic.cpp -o design
//   ./design

// Include the AC Math function that is exercised with this testbench
#include <ac_math/ac_hcordic.h>
using namespace ac_math;

// ==============================================================================
// Test Design
//   This simple function allows executing the ac_pow_cordic() function.
//   Template parameters are used to configure the bit-widths of the
//   input and output.

// Test for ac_fixed input and output.
template <int baseW, int baseI, int expW, int expI, bool expS, int outW, int outI>
void test_ac_pow_cordic_fixed(
  const ac_fixed<baseW, baseI, false, AC_TRN, AC_WRAP> &base_fixed,
  const ac_fixed<expW, expI, expS, AC_TRN, AC_WRAP> &expon_fixed,
  ac_fixed<outW, outI, false, AC_TRN, AC_WRAP> &output_fixed
)
{
  ac_pow_cordic(base_fixed, expon_fixed, output_fixed);
}

// Test for ac_float input and output.
template <int baseW, int baseI, int baseE, int expW, int expI, int expE, int outW, int outI, int outE>
void test_ac_pow_cordic_float(
  const ac_float<baseW, baseI, baseE> &base_float,
  const ac_float<expW, expI, expE> &expon_float,
  ac_float<outW, outI, outE> &output_float
)
{
  ac_pow_cordic(base_float, expon_float, output_float);
}

// Test for ac_std_float input and output.
template <bool OR_TF, int baseW, int baseE, int expW, int expE, int outW, int outE>
void test_ac_pow_cordic_stfloat(
  const ac_std_float<baseW, baseE> &base_stfloat,
  const ac_std_float<expW, expE> &expon_stfloat,
  ac_std_float<outW, outE> &output_stfloat
)
{
  ac_pow_cordic<OR_TF>(base_stfloat, expon_stfloat, output_stfloat);
}

// Test for ac_ieee_float input and output.
template <bool OR_TF, ac_ieee_float_format base_format, ac_ieee_float_format expon_format, ac_ieee_float_format out_format>
void test_ac_pow_cordic_ifloat(
  const ac_ieee_float<base_format> &base_ifloat,
  const ac_ieee_float<expon_format> &expon_ifloat,
  ac_ieee_float<out_format> &output_ifloat
)
{
  ac_pow_cordic<OR_TF>(base_ifloat, expon_ifloat, output_ifloat);
}

// ------------------------------------------------------------------------------
// Helper function for absolute value calculation. This can avoid any naming conflicts
// with other absolute value functions.

double abs_double(double x)
{
  return x >= 0 ? x : -x;
}

// ==============================================================================

#include <ac_math/ac_random.h>
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

// ==============================================================================
// Function: test_driver_fixed()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_fixed datatypes. It uses the type information to iterate
//   through a range of valid values on that type in order to compare the
//   precision of the cordic table pow model with the computed power output using
//   a standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

template <int baseW, int baseI, int expW, int expI, bool expS, int outW, int outI>
int test_driver_fixed(
  double &cumulative_max_error_fixed,
  const double allowed_error,
  const double threshold,
  bool details = false
)
{
  double max_error_fixed = 0.0; // reset for this run

  ac_fixed< baseW,  baseI, false, AC_TRN, AC_WRAP> base_fixed;
  ac_fixed<  expW,   expI,  expS, AC_TRN, AC_WRAP> expon_fixed;

  typedef ac_fixed<outW, outI, false, AC_TRN, AC_WRAP> T_out;
  T_out output_fixed;

  double lower_limit_base, upper_limit_base, step_base, lower_limit_expon, upper_limit_expon, step_expon;

  // set ranges and step size for testbench
  lower_limit_base = base_fixed.template set_val<AC_VAL_MIN>().to_double();
  upper_limit_base = base_fixed.template set_val<AC_VAL_MAX>().to_double();
  step_base        = base_fixed.template set_val<AC_VAL_QUANTUM>().to_double();

  lower_limit_expon = expon_fixed.template set_val<AC_VAL_MIN>().to_double();
  upper_limit_expon = expon_fixed.template set_val<AC_VAL_MAX>().to_double();
  step_expon        = expon_fixed.template set_val<AC_VAL_QUANTUM>().to_double();

  cout << "TEST: ac_pow_cordic() INPUTS: ";
  cout.width(38);
  cout << left << base_fixed.type_name();
  cout.width(38);
  cout << left << expon_fixed.type_name();
  cout << "OUTPUT: ";
  cout.width(38);
  cout << left << output_fixed.type_name();
  cout << "RESULT: ";

  // Dump the test details
  if (details) {
    cout << endl; // LCOV_EXCL_LINE
    cout << "  Ranges for input types:" << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit_base     = " << lower_limit_base << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit_base     = " << upper_limit_base << endl; // LCOV_EXCL_LINE
    cout << "    step_base            = " << step_base << endl; // LCOV_EXCL_LINE
    cout << "    lower_limit_expon     = " << lower_limit_expon << endl; // LCOV_EXCL_LINE
    cout << "    upper_limit_expon     = " << upper_limit_expon << endl; // LCOV_EXCL_LINE
    cout << "    step_expon            = " << step_expon << endl; // LCOV_EXCL_LINE
  }

  for (double i = lower_limit_base; i <= upper_limit_base; i += step_base) {
    for (double j = lower_limit_expon; j < upper_limit_expon; j += step_expon) {
      // Set values for inputs
      base_fixed = i;
      expon_fixed = j;
      if (base_fixed == 0 || expon_fixed == 0) { continue; }
      test_ac_pow_cordic_fixed(base_fixed, expon_fixed, output_fixed);
      double expected_value_fixed = ((T_out)pow(base_fixed.to_double(), expon_fixed.to_double())).to_double();
      double actual_value_fixed   = output_fixed.to_double();
      double this_error_fixed;

      // If expected value of output falls below the threshold, calculate absolute error instead of relative
      if (expected_value_fixed > threshold) {
        this_error_fixed = abs_double((expected_value_fixed - actual_value_fixed)/ expected_value_fixed)*100.0;
      } else {
        this_error_fixed = abs_double(expected_value_fixed - actual_value_fixed) * 100.0;
      }

      #ifdef DEBUG
      if (this_error_fixed > allowed_error) {
        cout << endl;
        cout << "  Error exceeds tolerance" << endl;
        cout << "  base_fixed           = " << base_fixed << endl;
        cout << "  expon_fixed          = " << expon_fixed << endl;
        cout << "  expected_value_fixed = " << expected_value_fixed << endl;
        cout << "  actual_value_fixed   = " << actual_value_fixed << endl;
        cout << "  this_error_fixed     = " << this_error_fixed << endl;
        cout << "  threshold            = " << threshold << endl;
        assert(false);
      }
      #endif

      if (this_error_fixed > max_error_fixed) { max_error_fixed = this_error_fixed; }
    }
  }

  bool passed = (max_error_fixed < allowed_error);

  if (passed) { printf("PASSED , max err %f\n", max_error_fixed); }
  else        { printf("FAILED , max err %f\n", max_error_fixed); } // LCOV_EXCL_LINE

  if (max_error_fixed > cumulative_max_error_fixed) { cumulative_max_error_fixed = max_error_fixed; }

  return 0;
}

// ==============================================================================
// Function: test_driver_float()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_float datatypes. It uses the type information to iterate
//   through a range of valid values on that type in order to compare the
//   precision of the cordic table pow model with the computed power output using
//   a standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

// num_it_per_loop: Number of iterations per test loop.
template <int num_it_per_loop, int baseW, int baseI, int baseE, int expW, int expI, int expE, int outW, int outI, int outE>
int test_driver_float(
  double &cumulative_max_error_float,
  const double allowed_error,
  const double threshold
)
{
  double max_error_float = 0.0; // reset for this run

  // typedefs for inputs and output.
  typedef ac_float<baseW, baseI, baseE> T_base;
  typedef ac_float<expW, expI, expE> T_expon;
  typedef ac_float<outW, outI, outE> T_out;

  cout << "TEST: ac_pow_cordic() INPUTS: ";
  cout.width(38);
  cout << left << T_base::type_name();
  cout.width(38);
  cout << left << T_expon::type_name();
  cout << "OUTPUT: ";
  cout.width(38);
  cout << left << T_out::type_name();
  cout << "RESULT: ";

  T_base sp_vals_base_quant;
  sp_vals_base_quant.template set_val<AC_VAL_QUANTUM>();
  T_base sp_vals_base_max;
  sp_vals_base_max.template set_val<AC_VAL_MAX>();

  // Always test the special base values: zero, quantum base value, unity and max base value.
  T_base sp_vals_base[] = {0.0, sp_vals_base_quant, 1.0, sp_vals_base_max};
  const int n_sp_vals_base = sizeof(sp_vals_base)/sizeof(T_base); // Number of special base values.

  for (int i = 0; i < num_it_per_loop; i++) {
    T_base base_float;
    if (i < n_sp_vals_base) {
      base_float = sp_vals_base[i];
    } else {
      // Every bit except for the first two MSBs are randomized.
      ac_int<baseW - 2, false> base_m_slc;
      ac_random(base_m_slc);
      ac_fixed<baseW, baseI, true> base_m = 0;
      base_m[baseW - 2] = 1; // Normalize base mantissa.
      base_m.template set_slc(0, base_m_slc.template slc<baseW - 2>(0));
      ac_int<baseE, true> base_e;
      ac_random(base_e);
      T_base base_float_temp(base_m, base_e, true);
      if (base_float_temp.mantissa() != base_m || base_float_temp.exp() != base_e) {
        AC_ASSERT(false, "base_float_temp not normalized correctly.");
      }
      base_float = base_float_temp;
    }
    T_expon sp_vals_expon_min;
    sp_vals_expon_min.template set_val<AC_VAL_MIN>();
    T_expon sp_vals_expon_quant;
    sp_vals_expon_quant.template set_val<AC_VAL_QUANTUM>();
    T_expon sp_vals_expon_max;
    sp_vals_expon_max.template set_val<AC_VAL_MAX>();
    T_expon sp_vals_expon[] = {
      sp_vals_expon_min, -sp_vals_expon_quant, 0.0, sp_vals_expon_quant, sp_vals_expon_max
    };
    const int n_sp_vals_expon = sizeof(sp_vals_expon)/sizeof(T_expon);
    
    for (int j = 0; j < num_it_per_loop; j++) {
      T_expon expon_float;
      if (j < n_sp_vals_expon) {
        expon_float = sp_vals_expon[j];
      } else {
        ac_int<expW - 2, false> expon_m_slc;
        ac_random(expon_m_slc);
        ac_fixed<expW, expI, true> expon_m = 0;
        expon_m[expW - 2] = 1; // Normalize exponent mantissa.
        expon_m.template set_slc(0, expon_m_slc.template slc<expW - 2>(0));
        // Randomly negate exponent mantissa.
        bool negate_expon_m = (rand()%2 == 0);
        expon_m = negate_expon_m ? ac_fixed<expW, expI, true>(-expon_m) : expon_m;
        ac_int<expE, true> expon_e;
        ac_random(expon_e);
        T_expon expon_float_tmp(expon_m, expon_e, true);
        if ((expon_float_tmp.mantissa() != expon_m || expon_float_tmp.exp() != expon_e) && expon_m_slc != 0) {
          AC_ASSERT(false, "expon_float_tmp not normalized correctly.");
        }
        expon_float = expon_float_tmp;
      }
      if (base_float == 0.0 && expon_float <= 0.0) {
        continue;
      }
      T_out output_max;
      output_max.template set_val<AC_VAL_MAX>();
      double expon_val_max = log(output_max.to_double())/log(base_float.to_double());
      if (isfinite(expon_val_max)) {
        if ((expon_val_max > 0.0 && expon_float > expon_val_max) || (expon_val_max < 0.0 && expon_float < expon_val_max)) {
          continue;
        }
      }
      T_out output_float;
      test_ac_pow_cordic_float(base_float, expon_float, output_float);
      double math_fn_pow = pow(base_float.to_double(), expon_float.to_double());
      bool is_pow_fin = isfinite(math_fn_pow);
      AC_ASSERT(is_pow_fin, "Infinite pow result.");
      double expected_value_float = ((T_out)math_fn_pow).to_double();
      double actual_value_float   = output_float.to_double();
      double this_error_float;

      // If expected value of output falls below the threshold, calculate absolute error instead of relative
      if (expected_value_float > threshold) {
        this_error_float = abs_double((expected_value_float - actual_value_float)/ expected_value_float)*100.0;
      } else {
        this_error_float = abs_double(expected_value_float - actual_value_float) * 100.0;
      }

      #ifdef DEBUG
      if (this_error_float > allowed_error) {
        cout << endl;
        cout << "  Error exceeds tolerance" << endl;
        cout << "  base_float           = " << base_float << endl;
        cout << "  expon_float          = " << expon_float << endl;
        cout << "  expected_value_float = " << expected_value_float << endl;
        cout << "  actual_value_float   = " << actual_value_float << endl;
        cout << "  this_error_float     = " << this_error_float << endl;
        cout << "  threshold            = " << threshold << endl;
        assert(false);
      }
      #endif

      if (this_error_float > max_error_float) {
        max_error_float = this_error_float;
      }
    }
  }

  bool passed = (max_error_float < allowed_error);

  if (passed) { printf("PASSED , max err %f\n", max_error_float); }
  else        { printf("FAILED , max err %f\n", max_error_float); } // LCOV_EXCL_LINE

  if (max_error_float > cumulative_max_error_float) { cumulative_max_error_float = max_error_float; }

  return 0;
}

// ==============================================================================
// Function: test_driver_stfloat()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_std_float datatypes. It uses the type information to iterate
//   through a range of valid values on that type in order to compare the
//   precision of the cordic table pow model with the computed power output using
//   a standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

// OR_TF: Value of OR_TF template parameter for ac_pow_cordic function.
// num_it_per_loop: Number of iterations per test loop.
template <bool OR_TF, int num_it_per_loop, int baseW, int baseE, int expW, int expE, int outW, int outE>
int test_driver_stfloat(
  double &cumulative_max_error_stfloat,
  const double allowed_error,
  const double threshold
)
{
  double max_error_stfloat = 0.0; // reset for this run

  typedef ac_std_float<baseW, baseE> T_base;
  typedef ac_std_float<expW, expE> T_expon;
  typedef ac_std_float<outW, outE> T_out;

  cout << "TEST: ac_pow_cordic() INPUTS: ";
  cout.width(38);
  cout << left << T_base::type_name();
  cout.width(38);
  cout << left << T_expon::type_name();
  cout << "OUTPUT: ";
  cout.width(38);
  cout << left << T_out::type_name();
  cout << "RESULT: ";
  
  // Always test the special base values: zero, min. denorm value, unity and max base value.
  T_base sp_vals_base[] = {
    T_base::zero(), T_base::denorm_min(), T_base::one(), T_base::min(), T_base::max()
  };
  const int n_sp_vals_base = sizeof(sp_vals_base)/sizeof(T_base); // Number of special base values.

  for (int i = 0; i < num_it_per_loop; i++) {
    T_base base_stfloat;
    if (i < n_sp_vals_base) {
      base_stfloat = sp_vals_base[i];
    } else {
      ac_int<baseW - baseE - 1, false> base_m;
      // Randomize base_m.
      ac_random(base_m);
      // Range of base_e values: 0 to bias*2.
      ac_int<baseE, false> base_e = rand()%(2*T_base::exp_bias + 1);
      ac_int<baseW, true> base_data = 0;
      base_data.set_slc(0, base_m);
      base_data.set_slc(baseW - baseE - 1, base_e);
      base_stfloat.set_data(base_data);
    }
    T_expon sp_vals_expon[] = {
      -T_expon::max(), -T_expon::min(), -T_expon::denorm_min(), T_expon::zero(), T_expon::denorm_min(),
      -T_expon::min(), T_expon::max()
    };
    const int n_sp_vals_expon = sizeof(sp_vals_expon)/sizeof(T_expon);
    
    for (int j = 0; j < num_it_per_loop; j++) {
      T_expon expon_stfloat;
      if (j < n_sp_vals_expon) {
        expon_stfloat = sp_vals_expon[j];
      } else {
        ac_int<expW - expE - 1, false> expon_m;
        // Randomize expon_m.
        ac_random(expon_m);
        // Range of expon_e values: 0 to bias*2.
        ac_int<expE, false> expon_e = rand()%(2*T_expon::exp_bias + 1);
        ac_int<expW, true> expon_data;
        expon_data[expW - 1] = rand()%2; // Randomly negate exponent value.
        expon_data.set_slc(0, expon_m);
        expon_data.set_slc(expW - expE - 1, expon_e);
        expon_stfloat.set_data(expon_data);
      }
      if (base_stfloat == T_base::zero() && expon_stfloat <= T_expon::zero()) {
        continue;
      }
      double expon_val_max = log((T_out::max()).to_double())/log(base_stfloat.to_double());
      if (isfinite(expon_val_max)) {
        bool cond_1 = (expon_val_max > 0.0 && expon_stfloat >= T_expon(expon_val_max));
        bool cond_2 = (expon_val_max < 0.0 && expon_stfloat <= T_expon(expon_val_max));
        if (cond_1 || cond_2) {
          continue;
        }
      }
      T_out output_stfloat;
      test_ac_pow_cordic_stfloat<OR_TF>(base_stfloat, expon_stfloat, output_stfloat);
      double math_fn_pow = pow(base_stfloat.to_double(), expon_stfloat.to_double());
      bool is_pow_fin = isfinite(math_fn_pow);
      #if 0
      if (!is_pow_fin) {
        bool cond_1 = (expon_val_max > 0.0 && expon_stfloat >= T_expon(expon_val_max));
        bool cond_2 = (expon_val_max < 0.0 && expon_stfloat <= T_expon(expon_val_max));
        bool cond_3 = (expon_stfloat < T_expon(expon_val_max));
        T_expon expon_val_conv = T_expon(expon_val_max);
        int x23123 = 2;
      }
      #endif
      AC_ASSERT(is_pow_fin, "Infinite pow result.");
      double expected_value_stfloat = ((T_out)math_fn_pow).to_double();
      double actual_value_stfloat   = output_stfloat.to_double();
      double this_error_stfloat;

      // If expected value of output falls below the threshold, calculate absolute error instead of relative
      if (expected_value_stfloat > threshold) {
        this_error_stfloat = abs_double((expected_value_stfloat - actual_value_stfloat)/ expected_value_stfloat)*100.0;
      } else {
        this_error_stfloat = abs_double(expected_value_stfloat - actual_value_stfloat) * 100.0;
      }

      #ifdef DEBUG
      if (this_error_stfloat > allowed_error) {
        cout << endl;
        cout << "  Error exceeds tolerance" << endl;
        cout << "  base_stfloat           = " << base_stfloat << endl;
        cout << "  expon_stfloat          = " << expon_stfloat << endl;
        cout << "  expected_value_stfloat = " << expected_value_stfloat << endl;
        cout << "  actual_value_stfloat   = " << actual_value_stfloat << endl;
        cout << "  this_error_stfloat     = " << this_error_stfloat << endl;
        cout << "  threshold              = " << threshold << endl;
        assert(false);
      }
      #endif

      if (this_error_stfloat > max_error_stfloat) {
        max_error_stfloat = this_error_stfloat;
      }
    }
  }

  bool passed = (max_error_stfloat < allowed_error);

  if (passed) { printf("PASSED , max err %f\n", max_error_stfloat); }
  else        { printf("FAILED , max err %f\n", max_error_stfloat); } // LCOV_EXCL_LINE

  if (max_error_stfloat > cumulative_max_error_stfloat) {
    cumulative_max_error_stfloat = max_error_stfloat;
  }

  return 0;
}

// ==============================================================================
// Function: test_driver_ifloat()
// Description: A templatized function that can be configured for certain bit-
//   widths of ac_ieee_float datatypes. It uses the type information to iterate
//   through a range of valid values on that type in order to compare the
//   precision of the cordic table pow model with the computed power output using
//   a standard C double type. The maximum error for each type is accumulated
//   in variables defined in the calling function.

// OR_TF: Value of OR_TF template parameter for ac_pow_cordic function.
// num_it_per_loop: Number of iterations per test loop.
template <bool OR_TF, int num_it_per_loop, ac_ieee_float_format base_format, ac_ieee_float_format expon_format, ac_ieee_float_format out_format>
int test_driver_ifloat(
  double &cumulative_max_error_ifloat,
  const double allowed_error,
  const double threshold
) {
  double max_error_ifloat = 0.0; // reset for this run

  typedef ac_ieee_float<base_format> T_base;
  typedef ac_ieee_float<expon_format> T_expon;
  typedef ac_ieee_float<out_format> T_out;
  
  const int baseW = T_base::width;
  const int baseE = T_base::e_width;
  const int expW = T_expon::width;
  const int expE = T_expon::e_width;

  cout << "TEST: ac_pow_cordic() INPUTS: ";
  cout.width(38);
  cout << left << T_base::type_name();
  cout.width(38);
  cout << left << T_expon::type_name();
  cout << "OUTPUT: ";
  cout.width(38);
  cout << left << T_out::type_name();
  cout << "RESULT: ";
  
  // Always test the special base values: zero, min. denorm value, unity and max base value.
  T_base sp_vals_base[] = {
    T_base::zero(), T_base::denorm_min(), T_base::one(), T_base::min(), T_base::max()
  };
  const int n_sp_vals_base = sizeof(sp_vals_base)/sizeof(T_base); // Number of special base values.

  for (int i = 0; i < num_it_per_loop; i++) {
    T_base base_ifloat;
    if (i < n_sp_vals_base) {
      base_ifloat = sp_vals_base[i];
    } else {
      ac_int<baseW - baseE - 1, false> base_m;
      // Randomize base_m.
      ac_random(base_m);
      // Range of base_e values: 0 to bias*2.
      ac_int<baseE, false> base_e = rand()%(2*T_base::exp_bias + 1);
      ac_int<baseW, true> base_data = 0;
      base_data.set_slc(0, base_m);
      base_data.set_slc(baseW - baseE - 1, base_e);
      base_ifloat.set_data(base_data);
    }
    T_expon sp_vals_expon[] = {
      -T_expon::max(), -T_expon::min(), -T_expon::denorm_min(), T_expon::zero(), T_expon::denorm_min(),
      -T_expon::min(), T_expon::max()
    };
    const int n_sp_vals_expon = sizeof(sp_vals_expon)/sizeof(T_expon);
    for (int j = 0; j < num_it_per_loop; j++) {
      T_expon expon_ifloat;
      if (j < n_sp_vals_expon) {
        expon_ifloat = sp_vals_expon[j];
      } else {
        ac_int<expW - expE - 1, false> expon_m;
        // Randomize expon_m.
        ac_random(expon_m);
        // Range of expon_e values: 0 to bias*2.
        ac_int<expE, false> expon_e = rand()%(2*T_expon::exp_bias + 1);
        ac_int<expW, true> expon_data;
        expon_data[expW - 1] = rand()%2; // Randomly negate exponent value.
        expon_data.set_slc(0, expon_m);
        expon_data.set_slc(expW - expE - 1, expon_e);
        expon_ifloat.set_data(expon_data);
      }
      if (base_ifloat == T_base::zero() && expon_ifloat <= T_expon::zero()) {
        continue;
      }
      double expon_val_max = log((T_out::max()).to_ac_float().to_double())/log(base_ifloat.to_ac_float().to_double());
      if (isfinite(expon_val_max)) {
        bool cond_1 = (expon_val_max > 0.0 && expon_ifloat >= T_expon(expon_val_max));
        bool cond_2 = (expon_val_max < 0.0 && expon_ifloat <= T_expon(expon_val_max));
        if (cond_1 || cond_2) {
          continue;
        }
      }
      T_out output_ifloat;
      test_ac_pow_cordic_ifloat<OR_TF>(base_ifloat, expon_ifloat, output_ifloat);
      double math_fn_pow = pow(base_ifloat.to_ac_float().to_double(), expon_ifloat.to_ac_float().to_double());
      bool is_pow_fin = isfinite(math_fn_pow);
      AC_ASSERT(is_pow_fin, "Infinite pow result.");
      double expected_value_ifloat = ((T_out)math_fn_pow).to_ac_float().to_double();
      double actual_value_ifloat   = output_ifloat.to_ac_float().to_double();
      double this_error_ifloat;

      // If expected value of output falls below the threshold, calculate absolute error instead of relative
      if (expected_value_ifloat > threshold) {
        this_error_ifloat = abs_double((expected_value_ifloat - actual_value_ifloat)/ expected_value_ifloat)*100.0;
      } else {
        this_error_ifloat = abs_double(expected_value_ifloat - actual_value_ifloat) * 100.0;
      }

      #ifdef DEBUG
      if (this_error_ifloat > allowed_error) {
        cout << endl;
        cout << "  Error exceeds tolerance" << endl;
        cout << "  base_ifloat           = " << base_ifloat << endl;
        cout << "  expon_ifloat          = " << expon_ifloat << endl;
        cout << "  expected_value_ifloat = " << expected_value_ifloat << endl;
        cout << "  actual_value_ifloat   = " << actual_value_ifloat << endl;
        cout << "  this_error_ifloat     = " << this_error_ifloat << endl;
        cout << "  threshold              = " << threshold << endl;
        assert(false);
      }
      #endif

      if (this_error_ifloat > max_error_ifloat) {
        max_error_ifloat = this_error_ifloat;
      }
    }
  }

  bool passed = (max_error_ifloat < allowed_error);

  if (passed) { printf("PASSED , max err %f\n", max_error_ifloat); }
  else        { printf("FAILED , max err %f\n", max_error_ifloat); } // LCOV_EXCL_LINE

  if (max_error_ifloat > cumulative_max_error_ifloat) {
    cumulative_max_error_ifloat = max_error_ifloat;
  }
  
  return 0;
}

int main(int argc, char *argv[])
{
  double max_error_fixed = 0.0, max_error_float = 0.0, max_error_stfloat = 0.0, max_error_ifloat = 0.0;
  double allowed_error = 0.5;
  double threshold = 0.005;

  cout << "=============================================================================" << endl;
  cout << "Testing function: ac_pow_cordic() - Allowed error " << allowed_error << endl;

  #if 0
  // template <int baseW, int baseI, int expW, int expI, bool expS, int outW, int outI>
  test_driver_fixed< 8, 3, 10, 2,  true, 40, 20>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 8, 2, 10, 2,  true, 40, 20>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 8, 2, 10, 2, false, 40, 20>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed< 8, 3, 10, 2, false, 40, 20>(max_error_fixed, allowed_error, threshold);
  test_driver_fixed<10, 4, 10, 2, false, 40, 20>(max_error_fixed, allowed_error, threshold);

  // template <int num_it_per_loop, int baseW, int baseI, int baseE, int expW, int expI, int expE, int outW, int outI, int outE>
  test_driver_float<1000, 16,  5, 6, 15,  4, 6, 32, 2, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<1000, 16, -3, 6, 15, -2, 6, 32, 2, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<1000, 20,  5, 6, 20,  4, 6, 32, 2, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<1000, 20, -3, 6, 20, -2, 6, 32, 2, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<1000, 16,  5, 8, 15,  4, 8, 32, 2, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<1000, 16, -3, 8, 15, -2, 8, 32, 2, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<1000, 20,  5, 8, 20,  4, 8, 32, 2, 10>(max_error_float, allowed_error, threshold);
  test_driver_float<1000, 20, -3, 8, 20, -2, 8, 32, 2, 10>(max_error_float, allowed_error, threshold);

  // template <bool OR_TF, int num_it_per_loop, int baseW, int baseE, int expW, int expE, int outW, int outE>
  test_driver_stfloat<false, 1000, 16, 5, 16, 5, 32,  8>(max_error_stfloat, allowed_error, threshold);
  test_driver_stfloat<false, 1000, 32, 8, 32, 8, 32,  8>(max_error_stfloat, allowed_error, threshold);
  #endif
  
  test_driver_stfloat<true,  1000, 16, 5, 16, 5, 64, 11>(max_error_stfloat, allowed_error, threshold);
  test_driver_stfloat<true,  1000, 32, 8, 16, 5, 64, 11>(max_error_stfloat, allowed_error, threshold);
  test_driver_stfloat<true,  1000, 16, 5, 32, 8, 64, 11>(max_error_stfloat, allowed_error, threshold);
  test_driver_stfloat<true,  1000, 32, 8, 32, 8, 64, 11>(max_error_stfloat, allowed_error, threshold);

  // template <bool OR_TF, int num_it_per_loop, ac_ieee_float_format base_format, ac_ieee_float_format expon_format, ac_ieee_float_format out_format>
  test_driver_ifloat<false, 1000, binary16, binary16, binary32>(max_error_ifloat, allowed_error, threshold);
  test_driver_ifloat<false, 1000, binary32, binary32, binary32>(max_error_ifloat, allowed_error, threshold);
  test_driver_ifloat<true,  1000, binary16, binary16, binary64>(max_error_ifloat, allowed_error, threshold);
  test_driver_ifloat<true,  1000, binary32, binary16, binary64>(max_error_ifloat, allowed_error, threshold);
  test_driver_ifloat<true,  1000, binary16, binary32, binary64>(max_error_ifloat, allowed_error, threshold);
  test_driver_ifloat<true,  1000, binary32, binary32, binary64>(max_error_ifloat, allowed_error, threshold);

  cout << "=============================================================================" << endl;
  cout << "  Testbench finished. Maximum errors observed across all bit-width variations:" << endl;
  cout << "    max_error_fixed   = " << max_error_fixed << endl;
  cout << "    max_error_float   = " << max_error_float << endl;
  cout << "    max_error_stfloat = " << max_error_stfloat << endl;
  cout << "    max_error_ifloat  = " << max_error_ifloat << endl;

  bool test_fail = (max_error_fixed > allowed_error) || (max_error_float > allowed_error) || (max_error_stfloat > allowed_error) || (max_error_ifloat > allowed_error);

  // If error limits on any test value have been crossed, the test has failed
  // Notify the user that the test was a failure if that is the case.
  if (test_fail) {
    cout << "  ac_pow_cordic - FAILED - Error tolerance(s) exceeded" << endl; // LCOV_EXCL_LINE
    cout << "=============================================================================" << endl; // LCOV_EXCL_LINE
    return -1; // LCOV_EXCL_LINE
  } else {
    cout << "  ac_pow_cordic - PASSED" << endl;
    cout << "=============================================================================" << endl;
  }

  return 0;
}
