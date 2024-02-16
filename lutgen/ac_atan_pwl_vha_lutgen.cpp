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
//******************************************************************************************
// Description:
//  This is another version of the ac_atan_pwl_lutgen file that is accuracy-targeted.
//  The user must supply a target error value, and the lutgen file will iterate over different
//  PWL implementations--each having a set of fixed PWL end points while also having a
//  different number of segments--until it finds one that meets the accuracy requirement.
//
// Usage:
//   $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include -O3 ac_atan_pwl_vha_lutgen.cpp -o ac_atan_pwl_vha_lutgen
//   ./ac_atan_pwl_vha_lutgen
//  results in a text file ac_atan_pwl_vha_lut_values.txt which can be pasted into
//  a locally modified version of ac_atan_pwl_vha.h.
// Notes:
// Revision History:
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//******************************************************************************************

// ac_reciprocal_pwl_vha.h is a header file that is included by default if the C++ compiler
// is C++11 or a later version. (earlier C++ versions will result in compilation errors).
// If the user does not wish to include/use the ac_reciprocal_pwl_vha.h header file for later C++
// versions, they must define the DO_NOT_USE_AC_RECIPROCAL_PWL_VHA macro, preferably as a
// compilation argument.
// (The ac_reciprocal_pwl_vha function is used to calculate input reciprocals and extend the domain
// of the ac_atan_pwl function from [0, 1) to [1, inf) through the following formula:
// atan(x) = pi/2 - atan(1/x) )
#if (defined(__GNUC__) && (__cplusplus < 201103L)) && !defined(DO_NOT_USE_AC_RECIPROCAL_PWL_VHA)
#define DO_NOT_USE_AC_RECIPROCAL_PWL_VHA
#endif
#if (defined(_MSC_VER) && (_MSC_VER < 1920) && !defined(__EDG__)) && !defined(DO_NOT_USE_AC_RECIPROCAL_PWL_VHA)
#define DO_NOT_USE_AC_RECIPROCAL_PWL_VHA
#endif

#ifndef DO_NOT_USE_AC_RECIPROCAL_PWL_VHA
#include <ac_math/ac_reciprocal_pwl_vha.h>
using namespace ac_math;
#endif

#include <ac_fixed.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#define _USE_MATH_DEFINES // Defined to enable usage of the "M_PI_2" macro.
#include <math.h>
#include <climits>
#include <cstddef>
#include <iostream>
using namespace std;

#include "helper_functions.h"

// This struct allows you to choose between a different number of max. iterations, based
// on whether the focus is on finding power-of-2 LUTs or not.
template <bool po2_only>
struct max_it_st {
  enum {
    // Since using power of two LUTs can cover a large range of LUT sizes pretty quickly,
    // We use a max. of only 10 iterations. (max. size = 1024 if the lower testing limit is 2)
    max_it = 10,
    upper_it = 1 << max_it,
  };
};

template <>
struct max_it_st<false> {
  enum {
    max_it = 127,
    upper_it = 1 + max_it,
  };
};

int main()
{
  // Define your target absolute error (this is NOT a percent value)
  const double t_err = 5e-5;

  enum {
    f_b_n_i = 15, // Number of fractional bits used for normalized inputs.
    w_pi_2 = 16, // Bitwidth used for variable that stores pi/2 value.
    po2_only = true, // Set to true -> program only looks at power of 2 (po2) LUTs.
    // max_it: Max. number of iterations the program goes through before it decides
    // that the target accuracy is unattainable.
    max_it = max_it_st<po2_only>::max_it,
    lower_it = 2, // Lower LUT size limit for testing.
    upper_it = max_it_st<po2_only>::upper_it, // Upper LUT size limit for testing
    nfrac_bits_u = sizeof(int)*CHAR_BIT - 2, // Max number of fractional bits used for fixed point LUTs.
    // The max number of fractional bits is constrained by the size of the int datatype, to allow for
    // the correct functioning of the o_ac_f function in the helper_functions header.
  };

  // Quick aside: Power-of-two LUTs often greatly reduce area, thanks not only to the
  // power-of-two size but also due to the elimination of multipliers needed for scaling.
  // (these multipliers should get replaced by a bit slicer which occupies less area)
  // The user is encouraged to give power-of-two LUTs a shot.

  // Arctangent domain limits.
  const double x_min = 0.0;
  const double x_max = 1.0;

  FILE *fp;

  // The output ROM values will be printed to a file in c++ syntax.
  // Define the filename below. The new file will be in the same folder
  // as this .cpp file
  const char filename[] = "ac_atan_pwl_vha_lut_values.txt";

  double x[upper_it + 1], x_sc[upper_it + 1], y[upper_it + 1], m[upper_it], c[upper_it];
  ac_fixed<nfrac_bits_u + 1, 1, false> m_fixed[upper_it], c_fixed[upper_it];

  // The program sets this flag to true if/when the target accuracy is achieved.
  bool target_achieved = false;
  unsigned nsegments_final, nfrac_bits_final;
  double sc_constant;
  double abs_error_max = 0, input_error_max;
#ifndef DO_NOT_USE_AC_RECIPROCAL_PWL_VHA
  typedef ac_fixed<w_pi_2, 1, false> pi_by_2_type;
  const pi_by_2_type pi_by_2 = M_PI_2;
#endif

  // LUT size = nsegments
  // lower_it and upper_it: Lower and upper iteration limits, respectively.
  for (unsigned nsegments = lower_it; nsegments <= upper_it;) {
    cout << "nsegments = " << nsegments << endl;
    const unsigned npoints = nsegments + 1;
    sc_constant = nsegments/(x_max - x_min);
    double x_val_inc = x_min;

    // Find slope and intercept for each segment.
    for (unsigned i = 0; i < npoints; i++) {
      y[i] = atan(x_val_inc);
      x[i] = x_val_inc;
      x_sc[i] = (x_val_inc - x_min)*sc_constant;
      x_val_inc += (1/sc_constant);
    }

    // Shift segments against the direction of function concavity, to improve accuracy.
    for (unsigned i = 0; i < nsegments; i++) {
      m[i] = y[i + 1] - y[i];
      c[i] = y[i];
      double x_mid = 0.5*(x[i + 1] + x[i]);
      double max_diff = pwl_new(x_mid, m, c, sc_constant, x_min, nsegments) - atan(x_mid);
      c[i] = c[i] - 0.5*(max_diff);
    }

    double y1_new, y2_new;

    // Correct slopes and intercepts in order for monotonicity to be maintained.
    for (unsigned i = 1; i < nsegments; i++) {
      y1_new = m[i - 1] + c[i - 1];
      y2_new = m[i] + c[i];
      m[i] = y2_new - y1_new;
      c[i] = y1_new;
    }

    //Use a negative power of two as an increment.
    const double increment = pow(2, -16.0);
    abs_error_max = 0;

    // Find the maximum PWL error over [x_min, x_max).
    for (double input_tb = x_min; input_tb < x_max; input_tb += increment) {
      const double expected = atan(input_tb);
      const double actual = pwl_new(input_tb, m, c, sc_constant, x_min, nsegments);
      const double abs_error = abs(expected - actual);
      if (abs_error > abs_error_max) {
        abs_error_max = abs_error;
        // input_error_max: input for which output absolute error is maximum.
        input_error_max = input_tb;
      }
    }
    
    cout << "abs_error_max (double version) = " << abs_error_max << endl;

    if (abs_error_max > t_err) {
      // If the absolute error you get using double PWL datatypes is higher than your target tolerance,
      // then it's almost guaranteed that fixed point PWL datatypes will fail to meet your tolerance
      // requirements as well. Instead of trying out fixed point implementations, update the loop iterator
      // and move on to the next loop iteration.
      nsegments = po2_only ? nsegments*2 : nsegments + 1;
      continue;
    }

    // The higher the accuracy of the double implementation, the more fractional bits you will need to
    // encode the PWL fixed point LUT values with to ensure an accuracy close to the double implementation.
    // Modify nfrac_bits_l accordingly to make sure you start at an appropriate lower limit for iterations and
    // avoid unncessary iterations that might waste runtime.
    const unsigned nfrac_bits_l = abs(floor(log2(abs_error_max))) + 1;
    double old_abs_error_max;

    for (unsigned nfrac_bits = nfrac_bits_l; nfrac_bits <= nfrac_bits_u; nfrac_bits++) {
      cout << "nfrac_bits = " << nfrac_bits << endl;
      // Quantize the double LUT values according the the number of fractional bits desired.
      for (unsigned i = 0; i < nsegments; i++) {
        m_fixed[i] = o_ac_f(m[i], nfrac_bits);
        c_fixed[i] = o_ac_f(c[i], nfrac_bits);
      }
      // Do the same for the domain and scaling constant variables.
      double x_min_q = o_ac_f(x_min, nfrac_bits).to_double();
      double x_max_q = o_ac_f(x_max, nfrac_bits).to_double();
      double sc_constant_q = o_ac_f(sc_constant, nfrac_bits).to_double();
      const ac_fixed<nfrac_bits_u, 0, false> quant_val = pow(2.0, -double(nfrac_bits)); // quantum value

#ifndef DO_NOT_USE_AC_RECIPROCAL_PWL_VHA
      if (x_max_q == 1) {
        // Find right- and left-hand limits at x = 1.
        ac_fixed<32, 1, false> rhl_1 = pi_by_2 - (m_fixed[nsegments - 1] + c_fixed[nsegments - 1]);
        ac_fixed<32, 1, false> lhl_1 = m_fixed[nsegments - 1] + c_fixed[nsegments - 1];
        // If the right-hand limit is lesser than the left hand limit, that indicates that the function is
        // non-monotonic at x = 1.
        if (rhl_1 < lhl_1) {
          // Adjust the slope value of last segment. The net effect is to lower the right end point.
          m_fixed[nsegments - 1] = (pi_by_2 >> 1) - c_fixed[nsegments - 1];
          m_fixed[nsegments - 1] = o_ac_f_trn(m_fixed[nsegments - 1].to_double(), nfrac_bits);
          rhl_1 = pi_by_2 - (m_fixed[nsegments - 1] + c_fixed[nsegments - 1]);
          lhl_1 = m_fixed[nsegments - 1] + c_fixed[nsegments - 1];
          if (rhl_1 < lhl_1) {
            // If the function is still non-monotonic at x = 1, that indicates faulty quantization. Lower
            // the slope value by one quantum, to fix this issue.
            m_fixed[nsegments - 1] = m_fixed[nsegments - 1] - quant_val;
          }
        }
      }
#endif

      // Check for monotonicity at every boundary between segments, reducing the slope if necessary.
      for (int i = 0; i < nsegments - 1; i++) {
        ac_fixed<32, 1, false> lhl = m_fixed[i] + c_fixed[i];
        ac_fixed<32, 1, false> rhl = c_fixed[i + 1];
        if (rhl < lhl) {
          m_fixed[i] = m_fixed[i] - quant_val;
        }
      }

      // Store the final fixed point values in a double array, to enable us to pass them to the pwl_new
      // function (the function does not accept fixed point inputs)
      double m_quant[upper_it], c_quant[upper_it];

      for (int i = 0; i < nsegments; i++) {
        m_quant[i] = m_fixed[i].to_double();
        c_quant[i] = c_fixed[i].to_double();
      }

      abs_error_max = 0;

#ifndef DO_NOT_USE_AC_RECIPROCAL_PWL_VHA
      // If your upper PWL domain limit is 1 and if you wish to use the ac_reciprocal_pwl_vha.h file, this variable
      // will enable testing for inputs > 1, all the way upto 100, using the reciprocal function for normalization into atan domain.
      const double input_tb_u = x_max_q == 1 ? 100 : x_max_q;
#else
      // If the reciprocal function isn't being used, we will only carry out testing for the domain of arctangent inputs,
      // i.e. [0, x_max_q)
      const double input_tb_u = x_max_q;
#endif
        //cout << "input_tb_u = " << input_tb_u << endl;

      // The increment value can be small enough to not update input_tb if we use double inputs. Instead, we use long double inputs
      // to get around this problem.
      for (long double input_tb = x_min_q; input_tb < input_tb_u; input_tb += increment) {
        const double input_tb_d = input_tb; // Convert iterator variable back into double, to enable compatibility with functions.
        const double expected = atan(input_tb_d);
#ifndef DO_NOT_USE_AC_RECIPROCAL_PWL_VHA
        ac_fixed<33, 1, false> recip_norm = ac_reciprocal_pwl_vha<ac_fixed<33, 1, false> >(ac_fixed<26, 10, false>(input_tb_d));
#endif
        double actual;
        if (input_tb_d < 1) {
          actual = pwl_new(input_tb_d, m_quant, c_quant, sc_constant_q, x_min_q, nsegments);
        } else {
#ifndef DO_NOT_USE_AC_RECIPROCAL_PWL_VHA
          // atan(x) = pi/2 - atan(1/x)
          actual = pi_by_2.to_double() - pwl_new(recip_norm.to_double(), m_quant, c_quant, sc_constant_q, x_min_q, nsegments);
#endif
        }
        const double abs_error = abs(expected - actual);
        if (abs_error > abs_error_max) {
          abs_error_max = abs_error;
          input_error_max = input_tb_d;
        }
      }
      
      cout << "abs_error_max = " << abs_error_max << endl;

      if (abs_error_max < t_err) {
        target_achieved = true; // We have achieved our target accuracy.
      }

      nfrac_bits_final = nfrac_bits; // nfrac_bits_final: Available outside for loop.
      // If the absolute error increases as the number of fractional bits increases, it indicates that the absolute error is
      // unlikely to get any lower than it already is. Break out of the current nfrac_bits iteration and try with a new nsegments value
      // to reduce runtime.
      if (nfrac_bits != nfrac_bits_l && old_abs_error_max < abs_error_max) { break; }
      old_abs_error_max = abs_error_max;
      // Break out of inner for loop if you've achieved your target accuracy.
      if (target_achieved) { break; }
    }

    nsegments_final = nsegments; // nsegments_final: Available outside for loop.
    if (target_achieved) { break; } // break out of the outer for loop as well, if you've achieved your target accuracy.
    // Update nsegments value.
    nsegments = po2_only ? nsegments*2 : nsegments + 1;
  }

  if (!target_achieved) {
    cout << "Target accuracy not achieved." << endl;
    return -1;
  }

  ostringstream mstrstream;
  ostringstream cstrstream;
  string mstr="";
  string cstr="";

  // Add elements to Objects that contain declaration of LUT arrays
  // in C++ syntax.
  for (int i = 0; i < nsegments_final; i++) {
    if (i == 0) {
      mstrstream << "{" << m_fixed[i] << ", ";
      cstrstream << "{" << c_fixed[i] << ", ";
    } else if (i == nsegments_final - 1) {
      mstrstream << m_fixed[i] << "}";
      cstrstream << c_fixed[i] << "}";
    } else {
      mstrstream << m_fixed[i] << ", ";
      cstrstream << c_fixed[i] << ", ";
    }
  }
  mstr = mstrstream.str();
  cstr = cstrstream.str();

  // Find max value in array, and see if it has any negative values. This helps figure out
  // the number of integer bits to use to store array values.
  double m_max_val, c_max_val;

  bool is_neg_m = is_neg_max_array(m, nsegments_final, m_max_val);
  bool is_neg_c = is_neg_max_array(c, nsegments_final, c_max_val);

  string is_neg_m_s = is_neg_m ? "true" : "false";
  string is_neg_c_s = is_neg_c ? "true" : "false";

  // Calculate number of int bits required to store various PWL variables.
  int m_int_bits = int_bits_calc(m_max_val, is_neg_m);
  int c_int_bits = int_bits_calc(c_max_val, is_neg_c);
  int x_min_int_bits = int_bits_calc(x_min, x_min < 0);
  int x_max_int_bits = int_bits_calc(x_max, x_max < 0);
  int sc_int_bits = int_bits_calc(sc_constant, false);

  string is_neg_x_min_s = (x_min < 0) ? "true" : "false";
  string is_neg_x_max_s = (x_max < 0) ? "true" : "false";

  // Find out whether the scaling constant is a power of two, which is only true its log2 value is an integer.
  double log2_pc = log2(sc_constant);
  bool pc_is_po2 = (trunc(log2_pc) == log2_pc);

  std::ofstream outfile(filename);
  outfile << "const int f_b_n_i = " << f_b_n_i << "; // Number of fractional bits to be assigned for the normalized input." << endl;
  outfile << "const unsigned n_segments_lut = " << nsegments_final << "; // Number of PWL segments." <<endl;
  outfile << "const int n_frac_bits = " << nfrac_bits_final << "; // Number of fractional bits used for storage of PWL values." << endl;
  if (pc_is_po2) {
    outfile << "// Since scaling constant is a positive power-of-two, multiplication with it is the same as left-shifting by " << log2_pc << "." << endl;
    outfile << "// Accordingly, the scaled normalized input will have " << log2_pc << " less fractional bits than the normalized input." << endl;
    outfile << "const int sc_input_frac_bits = f_b_n_i - " << log2_pc << ";" << endl;
  } else {
    outfile << "const int sc_input_frac_bits = f_b_n_i; // Number of fractional bits in scaled input" << endl;
  }
#ifndef DO_NOT_USE_AC_RECIPROCAL_PWL_VHA
  if (x_max == 1) {
    outfile << "const ac_fixed<" << pi_by_2.width << ", " << pi_by_2.i_width << ", false> pi_by_2 = " << pi_by_2 << ";" << endl;
  }
#endif
  outfile << "// Slope and intercept LUTs." << endl;
  outfile << "const ac_fixed<" << m_int_bits << " + n_frac_bits, " << m_int_bits << ", " << is_neg_m_s << "> m_lut[n_segments_lut] = " << mstr << ";" << endl;
  outfile << "const ac_fixed<" << c_int_bits << " + n_frac_bits, " << c_int_bits << ", " << is_neg_c_s << "> c_lut[n_segments_lut] = " << cstr << ";" << endl;
  outfile << "const ac_fixed<" << x_min_int_bits << " + n_frac_bits, " << x_min_int_bits << ", " << is_neg_x_min_s << "> x_min_lut = " << o_ac_f(x_min, nfrac_bits_final) << "; // Minimum limit of PWL domain." << endl;
  outfile << "const ac_fixed<" << sc_int_bits << " + n_frac_bits, " << sc_int_bits << ", false> sc_constant_lut = " << o_ac_f(sc_constant, nfrac_bits_final) << "; // Scaling constant." << endl;
  outfile.close();

  cout << endl;
  cout << __FILE__ << ", " << __LINE__ << ": Values are written, check \""
       << filename << "\" for " << "the required ROM values" << endl ;
  cout << "abs_error_max = " << abs_error_max << endl;
  cout << "input_error_max = " << input_error_max << endl;
  cout << "nfrac_bits_final = " << nfrac_bits_final << endl;

  return 0;
}

