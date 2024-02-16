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
//   g++ -I$MGC_HOME/shared/include ac_reciprocal_pwl_lutgen.cpp -o ac_reciprocal_pwl_lutgen
//   ./ac_reciprocal_pwl_lutgen
// results in a text file ac_reciprocal_pwl_lut_values.txt which can be pasted into
// a locally modified version of ac_reciprocal_pwl.h.

#include<ac_fixed.h>
#include<stdio.h>
#include<string>
#include<sstream>
#include<fstream>
#include<cmath>
#include<math.h>
#include<iostream>
using namespace std;

#include "helper_functions.h"

int main()
{
  FILE *fp;

  //Define the number of points in your LUT.
  const unsigned npoints = 9;
  const unsigned nsegments = npoints - 1;
  const double x_min = 0.5;
  const double x_max = 1.0;
  const double prop_constant = nsegments/(x_max - x_min);
  double m[nsegments];
  double c[nsegments];

  //The output ROM values will be printed to a file in c++ syntax.
  //Define the filename below. The new file will be in the same folder
  //as this .cpp file
  const char filename[] = "ac_reciprocal_pwl_lut_values.txt";

  //Define your domain for calculation. In case you want to use the
  //value of pi, use M_PI, a math.h macro that has the value of pi.
  //(e.g. "x_max = M_PI/4")
  double x[npoints];
  double x_sc[npoints];
  double y[npoints];
  double x_val_inc = x_min;
  ostringstream mstrstream;
  ostringstream cstrstream;
  string mstr="";
  string cstr="";

  //Find slope and intercept for each segment.
  for (int i = 0; i < npoints; i++) {
    y[i] = 1.0 / x_val_inc;
    x[i] = x_val_inc;
    x_sc[i] = (x_val_inc - x_min)*prop_constant;
    x_val_inc += (1/prop_constant);
  }

  //Shift segments downward
  for (int i = 0; i < nsegments; i++) {
    m[i] = y[i + 1] - y[i];
    c[i] = y[i];
    double x_mid = 0.5*(x[i + 1] + x[i]);
    double max_diff = (pwl_new(x_mid, m, c, prop_constant, x_min, nsegments) - (1.0 / x_mid));
    c[i] = c[i] - 0.5*(max_diff);
    //cout << "max_diff = " << max_diff << endl;
  }

  double y1_new, y2_new;

  //Correct slopes and intercepts in order for monotonicity to
  //be maintained.
  for (int i = 1; i < nsegments; i++) {
    y1_new = m[i - 1] + c[i - 1];
    y2_new = m[i] + c[i];
    m[i] = y2_new - y1_new;
    c[i] = y1_new;
  }

  //Use a negative power of two as an increment.
  double increment = 1.0 / 65536.0;
  double abs_error, abs_error_max = 0, rel_error_max = 0, input_error_max;

  for (double input_tb = x_min; input_tb < x_max; input_tb += increment) {
    double expected = (1.0 / input_tb);
    double actual = pwl_new(input_tb, m, c, prop_constant, x_min, nsegments);
    double rel_error = abs( (expected - actual) / expected) * 100;
    if (rel_error > rel_error_max) { rel_error_max = rel_error; }
    abs_error = abs(expected - actual);
    if (abs_error > abs_error_max) {
      abs_error_max = abs_error;
      input_error_max = input_tb;
    }
  }

  int nfrac_bits = abs(floor(log2(abs_error_max))) + 1;

  ac_fixed<128, 64, true> m_fixed[nsegments], c_fixed[nsegments];

  // Quantize double values according to fixed point precision.
  for (int i = 0; i < nsegments; i++) {
    m_fixed[i] = o_ac_f(m[i], nfrac_bits);
    c_fixed[i] = o_ac_f(c[i], nfrac_bits);
  }

  // Find the quantum value of an fixed point variable with nfrac_bits number of fractional bits.
  const ac_fixed<128, 64, true> quant_val = pow(2, double(-nfrac_bits));

  // Note: The montonicity checks in the if branch below depend upon the PWL end points being 0.5 and 1.
  // For any other end points, you may have to manually change slope-intercept values to maintain monotonicity.
  if (x_min == 0.5 && x_max == 1) {
    // Find right- and left-hand limits at x = 1
    const double rhl_1 = 0.5*c_fixed[0].to_double();
    double lhl_1 = (m_fixed[nsegments - 1] + c_fixed[nsegments - 1]).to_double();
    while (rhl_1 > lhl_1) {
      // Keep raising the end point of the last segment by adding a quantum value to the slope, until
      // the right hand limit is lesser than/equal to the left hand limit.
      m_fixed[nsegments - 1] = m_fixed[nsegments - 1] + quant_val;
      lhl_1 = (m_fixed[nsegments - 1] + c_fixed[nsegments - 1]).to_double();
    }
  }

  // Check left- and right-hand limits of the PWL function for each segment boundary.
  for (int i = 0; i < nsegments - 1; i++) {
    double lhl = (m_fixed[i] + c_fixed[i]).to_double();
    double rhl  = (c_fixed[i + 1]).to_double();
    if (lhl < rhl) {
      // If the PWL output is not decreasing at the segment boundary, increase the slope of the preceding
      // segment so as to raise the left hand limit and make the PWL output decrease across the boundary.
      m_fixed[i] = m_fixed[i] + quant_val;
    }
  }

  //Add elements to Objects that contain declaration of LUT arrays
  //in C++ syntax.
  for (int i = 0; i < nsegments; i++) {
    if (i == 0) {
      mstrstream << "{" << m_fixed[i] << ", ";
      cstrstream << "{" << c_fixed[i] << ", ";
    } else if (i == nsegments - 1) {
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
  bool is_neg_m, is_neg_c;

  is_neg_m = is_neg_max_array(m, nsegments, m_max_val);
  is_neg_c = is_neg_max_array(c, nsegments, c_max_val);

  string is_neg_m_s = is_neg_m ? "true" : "false";
  string is_neg_c_s = is_neg_c ? "true" : "false";

  int m_int_bits = int_bits_calc(m_max_val, is_neg_m);
  int c_int_bits = int_bits_calc(c_max_val, is_neg_c);
  int x_min_int_bits = int_bits_calc(x_min, x_min < 0);
  int p_c_int_bits = int_bits_calc(prop_constant, false);

  string is_neg_x_min_s = (x_min < 0) ? "true" : "false";

  std::ofstream outfile(filename);
  outfile << "const unsigned n_segments_lut = " << nsegments << "; // Number of PWL segments." <<endl;
  outfile << "const int n_frac_bits = " << nfrac_bits << "; // Number of fractional bits" << endl;
  if (pc_is_po2) {
    outfile << "// Since scaling constant is a positive power-of-two, multiplication with it is the same as left-shifting by " << log2_pc << "." << endl;
    outfile << "// Accordingly, the scaled normalized input will have " << log2_pc << " less fractional bits than the normalized input, provided that this" << endl;
    outfile << "// number of fractional bits is lesser than n_frac_bits. If not, the number of fractional bits in the scaled input is set to n_frac_bits." << endl;
    outfile << "const int sc_input_frac_bits = AC_MAX(1, AC_MIN(n_frac_bits, W - " << log2_pc << " - int(S))); // One less bit is used if the function input is signed, due to how ac_normalize works." << endl;
  } else {
    outfile << "const int sc_input_frac_bits = n_frac_bits; // Number of fractional bits in scaled input" << endl;
  }
  outfile << "// Slope and intercept LUT values." << endl;
  outfile << "const ac_fixed<" << m_int_bits << " + n_frac_bits, " << m_int_bits << ", " << is_neg_m_s << "> m_lut[n_segments_lut] = " << mstr << ";" << endl;
  outfile << "const ac_fixed<" << c_int_bits << " + n_frac_bits, " << c_int_bits << ", " << is_neg_c_s << "> c_lut[n_segments_lut] = " << cstr << ";" << endl;
  outfile << "const ac_fixed<" << x_min_int_bits << " + n_frac_bits, " << x_min_int_bits << ", " << is_neg_x_min_s << "> x_min_lut = " << o_ac_f(x_min, nfrac_bits) << "; // Minimum limit of PWL domain" << endl;
  outfile << "const ac_fixed<" << p_c_int_bits << " + n_frac_bits, " << p_c_int_bits << ", false> sc_constant_lut = " << o_ac_f(prop_constant, nfrac_bits) << "; // Scaling constant" << endl;
  outfile.close();

  cout << endl;
  cout << __FILE__ << ", " << __LINE__ << ": Values are written, check \""
       << filename << "\" for " << "the required ROM values" << endl ;
  cout << "abs_error_max = " << abs_error_max << endl;
  cout << "input_error_max = " << input_error_max << endl;
  cout << "nfrac_bits = " << nfrac_bits << endl;
  cout << "rel_error_max = " << rel_error_max << endl;

  return 0;
}

