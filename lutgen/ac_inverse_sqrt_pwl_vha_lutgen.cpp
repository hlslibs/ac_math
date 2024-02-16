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
// Usage:
//   g++ -I$MGC_HOME/shared/include ac_inverse_sqrt_pwl_vha_lutgen.cpp -o ac_inverse_sqrt_pwl_vha_lutgen
//   ./ac_inverse_sqrt_pwl_vha_lutgen
// results in a text file ac_inverse_sqrt_pwl_vha_lut_values.txt which can be pasted into
// a locally modified version of ac_inverse_sqrt_pwl_vha.h.

#include<ac_fixed.h>
#include<stdio.h>
#include<string>
#include<sstream>
#include<fstream>
#include<cmath>
#include<math.h>
#include <climits>
#include <cstddef>
#include<iostream>
using namespace std;

#include "helper_functions.h"

template <bool po2_only>
struct max_it_st {
  enum {
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
  // Define your target absolute error (this is NOT a percentage value)
  const double t_err = 0.00007071067;

  enum {
    po2_only = true, // Set to true -> program only looks at power of 2 (po2) LUTs.
    // max_it: Max. number of iterations the program goes through before it decides
    // that the target accuracy is unattainable.
    max_it = max_it_st<po2_only>::max_it,
    lower_it = 2,
    upper_it = max_it_st<po2_only>::upper_it,
    nfrac_bits_u = sizeof(int)*CHAR_BIT - 2,
  };

  // Domain limits.
  const double x_min = 0.5;
  const double x_max = 1.0;

  FILE *fp;

  // The output LUT values will be printed to a file in c++ syntax.
  // Define the filename below. The new file will be in the same folder
  // as this .cpp file
  const char filename[] = "ac_inverse_sqrt_pwl_vha_lut_values.txt";

  double x[upper_it + 1], x_sc[upper_it + 1], y[upper_it + 1], m[upper_it], c[upper_it];
  ac_fixed<nfrac_bits_u + 3, 3, true> m_fixed[upper_it], c_fixed[upper_it];

  // The program sets this flag to true if/when the target accuracy is achieved.
  bool target_achieved = false;
  unsigned nsegments_final, nfrac_bits_final;
  double prop_constant;
  double abs_error_max = 0, input_error_max;

  // LUT size = nsegments
  // lower_it and upper_it: Lower and upper iteration limits, respectively.
  for (unsigned nsegments = lower_it; nsegments <= upper_it;) {
    cout << "nsegments = " << nsegments << endl;
    const unsigned npoints = nsegments + 1;
    prop_constant = nsegments/(x_max - x_min);
    double x_val_inc = x_min;

    // Find slope and intercept for each segment.
    for (unsigned i = 0; i < npoints; i++) {
      y[i] = 1/sqrt(x_val_inc);
      x[i] = x_val_inc;
      x_sc[i] = (x_val_inc - x_min)*prop_constant;
      x_val_inc += (1/prop_constant);
    }

    // Shift segments against the direction of function concavity, to improve accuracy.
    for (unsigned i = 0; i < nsegments; i++) {
      m[i] = y[i + 1] - y[i];
      c[i] = y[i];
      double x_mid = 0.5*(x[i + 1] + x[i]);
      double max_diff = pwl_new(x_mid, m, c, prop_constant, x_min, nsegments) - 1/sqrt(x_mid);
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
      const double expected = 1/sqrt(input_tb);
      const double actual = pwl_new(input_tb, m, c, prop_constant, x_min, nsegments);
      const double abs_error = abs(expected - actual);
      if (abs_error > abs_error_max) {
        abs_error_max = abs_error;
        // input_error_max: input for which output absolute error is maximum.
        input_error_max = input_tb;
      }
    }

    cout << "abs_error_max = " << abs_error_max << endl;

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
      double prop_constant_q = o_ac_f(prop_constant, nfrac_bits).to_double();
      ac_fixed<nfrac_bits_u, 0, false> quant_val = 0;
      quant_val[nfrac_bits_u - nfrac_bits] = 1;

      // Note: The montonicity checks in the if branch below depend upon the PWL end points being 0.5 and 1.
      // For any other end points, you may have to manually change slope-intercept values to maintain monotonicity.
      if (x_min == 0.5 && x_max == 1) {
        // The precision of this value, and the value itself, depend on the precision/value used for the corresponding
        // inverseroot2 value in the ac_inverse_sqrt_pwl header.
        // If you change the precision or the value in the header file, make sure you change it here too.
        const ac_fixed <32, 0, false> inverseroot2 = 0.707106781186547524400844362;
        // Find the left- and right-hand limits as x tends to 1.
        double lhl_1 = (m_fixed[nsegments - 1] + c_fixed[nsegments - 1]).to_double();
        // The precision/rounding for the value of the right-hand limit as x tends to 1 depends on the precision/rounding
        // used in the fixed point header.
        // (ac_fixed<2*n_frac_bits, 0, false, q_mode_temp>, where q_mode_temp is set to AC_TRN by default)
        // If the user changes the precision/rounding from the default in the pwl header, they must change
        // it here too.
        const double rhl_1 = o_ac_f_trn((c_fixed[0]*inverseroot2).to_double(), 2*nfrac_bits).to_double();

        while (rhl_1 > lhl_1) {
          // Keep raising the end point of the last segment by adding a quantum value to the slope, until
          // the left hand limit is greater than/equal to the right hand limit.
          m_fixed[nsegments - 1] = m_fixed[nsegments - 1] + quant_val;
          lhl_1 = (m_fixed[nsegments - 1] + c_fixed[nsegments - 1]).to_double();
        }
      }

      // Check for monotonicity at every boundary between segments, increasing the slope if necessary.
      for (int i = 0; i < nsegments - 1; i++) {
        ac_fixed<32, 1, false> lhl = m_fixed[i] + c_fixed[i];
        ac_fixed<32, 1, false> rhl = c_fixed[i + 1];
        if (rhl > lhl) {
          m_fixed[i] = m_fixed[i] + quant_val;
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

      for (double input_tb = x_min_q; input_tb < x_max_q; input_tb += increment) {
        const double expected = 1/sqrt(input_tb);
        double actual = pwl_new(input_tb, m_quant, c_quant, prop_constant_q, x_min_q, nsegments);
        const double abs_error = abs(expected - actual);
        if (abs_error > abs_error_max) {
          abs_error_max = abs_error;
          input_error_max = input_tb;
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

  // Find out whether the proportionality constant is a power of two, which is only true its log2 value is an integer.
  double log2_pc = log2(prop_constant);
  bool pc_is_po2 = (trunc(log2_pc) == log2_pc);

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
  int sc_int_bits = int_bits_calc(prop_constant, false);

  string is_neg_x_min_s = (x_min < 0) ? "true" : "false";
  string is_neg_x_max_s = (x_max < 0) ? "true" : "false";

  std::ofstream outfile(filename);
  outfile << "const unsigned n_segments_lut = " << nsegments_final << "; // Number of PWL segments." <<endl;
  outfile << "const int n_frac_bits = " << nfrac_bits_final << "; // Number of fractional bits" << endl;
  if (pc_is_po2) {
    outfile << "// Since scaling constant is a positive power-of-two, multiplication with it is the same as left-shifting by " << log2_pc << "." << endl;
    outfile << "// Accordingly, the scaled normalized input will have " << log2_pc << " less fractional bits than the normalized input, provided that this" << endl;
    outfile << "// number of fractional bits is lesser than n_frac_bits. If not, the number of fractional bits in the scaled input is set to n_frac_bits." << endl;
    outfile << "const int sc_input_frac_bits = AC_MAX(1, AC_MIN(n_frac_bits, W1 - " << log2_pc << ")); // One less bit is used if the function input is signed, due to how ac_normalize works." << endl;
  } else {
    outfile << "const int sc_input_frac_bits = n_frac_bits; // Number of fractional bits in scaled input" << endl;
  }
  outfile << "// Slope and intercept LUT values." << endl;
  outfile << "const ac_fixed<" << m_int_bits << " + n_frac_bits, " << m_int_bits << ", " << is_neg_m_s << "> m_lut[n_segments_lut] = " << mstr << ";" << endl;
  outfile << "const ac_fixed<" << c_int_bits << " + n_frac_bits, " << c_int_bits << ", " << is_neg_c_s << "> c_lut[n_segments_lut] = " << cstr << ";" << endl;
  outfile << "const ac_fixed<" << x_min_int_bits << " + n_frac_bits, " << x_min_int_bits << ", " << is_neg_x_min_s << "> x_min_lut = " << o_ac_f(x_min, nfrac_bits_final) << "; // Minimum limit of PWL domain" << endl;
  outfile << "const ac_fixed<" << sc_int_bits << " + n_frac_bits, " << sc_int_bits << ", false> sc_constant_lut = " << o_ac_f(prop_constant, nfrac_bits_final) << "; // Scaling constant" << endl;
  outfile.close();

  cout << endl;
  cout << __FILE__ << ", " << __LINE__ << ": Values are written, check \""
       << filename << "\" for " << "the required ROM values" << endl ;
  cout << "abs_error_max = " << abs_error_max << endl;
  cout << "input_error_max = " << input_error_max << endl;
  cout << "nfrac_bits_final = " << nfrac_bits_final << endl;

  return 0;
}
