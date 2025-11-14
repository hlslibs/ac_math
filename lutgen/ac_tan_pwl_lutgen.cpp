/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 2025.4                                              *
 *                                                                        *
 *  Release Date    : Tue Nov 11 17:44:22 PST 2025                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 2025.4.0                                            *
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
//   g++ -I$MGC_HOME/shared/include -Wall -fsanitize=address ac_tan_pwl_lutgen.cpp -o ac_tan_pwl_lutgen
//   ./ac_tan_pwl_lutgen
// results in a text file ac_tan_pwl_lut_values.txt which can be pasted into
// a locally modified version of ac_tan_pwl.h.

#include<ac_fixed.h>
#include<stdio.h>
#include<string>
#include<sstream>
#include<fstream>
#include<cmath>
#include<math.h>
#include<vector>
#include<iostream>

#include "helper_functions.h"

#define PRINT_EXPR(expr) std::cout << #expr << " = " << (expr) << std::endl

template <class T_stream>
void show_help(T_stream &stream, char* argv[]) {
  stream << "Usage : " << argv[0] << " [command-line options] <optional filename>" << std::endl;
  stream << "Available command-line options:" << std::endl;
  stream << "  -h : Display this message." << std::endl;
  stream << "  -m : Change LUT values slightly to ensure a monotonically non-decreasing output." << std::endl;
}

int parse_opts(const int argc, char* argv[], bool &mono_guarantee, std::string &filename) {
  mono_guarantee = false;
  filename = "";

  for (int i = 1; i < argc; ++i) {
    std::string opt_str = argv[i];
    
    if (opt_str == std::string("-h")) {
      show_help(std::cout, argv);
      return 1; // Exit code of 1 -> Quit C++ execution with a zero exit code.
    } else if (opt_str == std::string("-m")) {
      std::cout << "PWL approximation is guaranteed to give a monotonic output for inputs in the domain [0, pi/4)." << std::endl;
      std::cout << "You may need to set the extra_monotonicity_check template argument for ac_tan_pwl(...)" << std::endl;
      std::cout << "to true in order to ensure monotonicity for inputs in the domain [0, pi/2)." << std::endl;
      std::cout << "This setting may not give as high an accuracy as manually adjusting LUT values for monotonicity." << std::endl;
      mono_guarantee = true;
    } else {
      if (!filename.empty()) {
        std::cerr << "Attempted to assign filename twice." << std::endl;
        std::cerr << "Previous filename : " << filename << std::endl;
        std::cerr << "Current filename : " << opt_str << std::endl;
        return 2;
      } // Exit code of 2 -> Quit C++ execution with a non-zero exit code.
      filename = opt_str;
    }
  }

  return 0;
}

template <int size>
void calc_err_and_frac_bits(
  const double x_min, const double x_max,
  const double (&m)[size], const double (&c)[size],
  const double prop_constant,
  double &abs_error_max, double &rel_error_max, double &input_error_max
) {
  //Use an inverse power of two (2^-16) as an increment.
  double increment = 1.0 / 65536.0;
  abs_error_max = 0.0;
  rel_error_max = 0.0;

  for (double input_tb = x_min; input_tb < x_max; input_tb += increment) {
    double expected = tan(input_tb);
    double actual = pwl_new(input_tb, m, c, prop_constant, x_min, size);
    double rel_error = abs( (expected - actual) / expected) * 100;
    if (rel_error > rel_error_max) { rel_error_max = rel_error; }
    double abs_error = abs(expected - actual);
    if (abs_error > abs_error_max) {
      abs_error_max = abs_error;
      input_error_max = input_tb;
    }
  }
}

int main(int argc, char* argv[])
{
  //Define the number of points in your LUT.
  //Note that this is declared on a separate line to aid in internal testing.
  const unsigned npoints = 9;
  const unsigned nsegments = npoints - 1;
  double x_min = 0;
  double x_max = M_PI_4;
  double prop_constant = nsegments/(x_max - x_min);
  double m[nsegments];
  double c[nsegments];

  bool mono_guarantee;
  std::string filename;
  int ret_code = parse_opts(argc, argv, mono_guarantee, filename);
  if (ret_code == 1) { return  0; }
  else if (ret_code != 0) { return -1; }

  // If the output filename wasn't already specified via the command line, assign a default value.
  if (filename.empty()) {
    filename = "ac_tan_pwl_lut_values.txt";
  }

  if (!mono_guarantee) {
    std::cout << "PWL approximation is not guaranteed to be monotonic for inputs in the domain [0, pi/4)." << std::endl;
    std::cout << "You may need to manually adjust the slopes for some of the segments to ensure monotonicity." << std::endl;
    std::cout << "If you want a PWL implementation which gives a non-monotonic output, consider using the" << std::endl;
    std::cout << "command-line option -m" << std::endl;
  }

  double x[npoints];
  double y[npoints];
  double x_val_inc = x_min;
  std::ostringstream mstrstream;
  std::ostringstream cstrstream;
  std::string mstr="";
  std::string cstr="";

  //Find slope and intercept for each segment.
  for (int i = 0; i < int(npoints); i++) {
    y[i] = tan(x_val_inc);
    x[i] = x_val_inc;
    x_val_inc += (1/prop_constant);
  }

  //Shift segments downward
  for (int i = 0; i < int(nsegments); i++) {
    m[i] = y[i + 1] - y[i];
    c[i] = y[i];
    double x_mid = 0.5*(x[i + 1] + x[i]);
    double max_diff = (pwl_new(x_mid, m, c, prop_constant, x_min, nsegments) - tan(x_mid));
    c[i] = c[i] - 0.5*(max_diff);
  }

  //Correct slopes and intercepts in order for monotonicity to
  //be maintained.
  mono_adjust(m, c);
  
  double abs_error_max, rel_error_max, input_error_max;
  calc_err_and_frac_bits(x_min, x_max, m, c, prop_constant, abs_error_max, rel_error_max, input_error_max);

  int nfrac_bits = abs(floor(log2(abs_error_max))) + 1;

  // Convert LUT values to fixed point precision.
  conv_arr_to_fxpt(m, nfrac_bits);
  conv_arr_to_fxpt(c, nfrac_bits);
  std::vector<double> x_min_vals;
  std::vector<double> x_max_vals;
  std::vector<double> prop_constant_vals;

  // Store LUT parameter (x_min, x_max and prop_constant) values with and without rounding.
  // The set of parameter values which leads to the lowest maximum absolute error over
  // the domain of [x_min, x_max) is the final set that's chosen.
  double x_min_trn = o_ac_f_trn(x_min, nfrac_bits).to_double();
  double x_min_rnd = o_ac_f(x_min, nfrac_bits).to_double();
  double x_max_trn = o_ac_f_trn(x_max, nfrac_bits).to_double();
  double x_max_rnd = o_ac_f(x_max, nfrac_bits).to_double();
  double prop_constant_trn = o_ac_f_trn(prop_constant, nfrac_bits).to_double();
  double prop_constant_rnd = o_ac_f(prop_constant, nfrac_bits).to_double();

  x_min_vals.push_back(x_min_trn);
  if (x_min_trn != x_min_rnd) { x_min_vals.push_back(x_min_rnd); }
  x_max_vals.push_back(x_max_trn);
  if (x_max_trn != x_max_rnd) { x_max_vals.push_back(x_max_rnd); }
  prop_constant_vals.push_back(prop_constant_trn);
  if (prop_constant_trn != prop_constant_rnd) { prop_constant_vals.push_back(prop_constant_rnd); }

  for (const double x_min_ : x_min_vals) {
    for (const double x_max_ : x_max_vals) {
      for (const double prop_constant_ : prop_constant_vals) {
        //ac_fixed<n_frac_bits + int_bits, int_bits, false> x_in_sc = (input_int - x_min_lut) * sc_constant_lut;
        double x_max_sc = (x_max_ - x_min_)*prop_constant_;
        AC_ASSERT(x_max_sc >= 0, "Maximum scaled value must not be negative.");
        // If the maximum scaled input isn't lesser than nsegments, we can get an out-of-bounds access on the LUT arrays.
        if (x_max_sc >= double(nsegments)) { continue; }
        double abs_error_max_, rel_error_max_, input_error_max_;
        calc_err_and_frac_bits(x_min_, x_max_, m, c, prop_constant_, abs_error_max_, rel_error_max_, input_error_max_);
        if (abs_error_max_ < abs_error_max) {
          abs_error_max = abs_error_max_;
          rel_error_max = rel_error_max_;
          input_error_max = input_error_max_;
          x_min = x_min_;
          x_max = x_max_;
          prop_constant = prop_constant_;
        }
      }
    }
  }

  if (mono_guarantee) {
    // Adjust for monotonicity again.
    mono_adjust(m, c);
  }

  dump_array_to_stream(m, mstrstream);
  dump_array_to_stream(c, cstrstream);

  mstr = mstrstream.str();
  cstr = cstrstream.str();

  // Find max value in array, and see if it has any negative values. This helps figure out
  // the number of integer bits to use to store array values.
  double m_max_val, c_max_val;
  bool is_neg_m, is_neg_c;

  is_neg_m = is_neg_max_array(m, nsegments, m_max_val);
  is_neg_c = is_neg_max_array(c, nsegments, c_max_val);

  std::string is_neg_m_s = is_neg_m ? "true" : "false";
  std::string is_neg_c_s = is_neg_c ? "true" : "false";

  int m_int_bits = int_bits_calc(m_max_val, is_neg_m);
  int c_int_bits = int_bits_calc(c_max_val, is_neg_c);
  int x_min_int_bits = int_bits_calc(x_min, x_min < 0);
  int x_max_int_bits = int_bits_calc(x_max, x_max < 0);
  int p_c_int_bits = int_bits_calc(prop_constant, false);

  std::string is_neg_x_min_s = (x_min < 0) ? "true" : "false";
  std::string is_neg_x_max_s = (x_max < 0) ? "true" : "false";

  std::ofstream outfile(filename.c_str());
  outfile << "const unsigned n_segments_lut = " << nsegments << ";" << std::endl;
  outfile << "const int n_frac_bits = " << nfrac_bits << ";" << std::endl;
  outfile << "const ac_fixed<" << m_int_bits << " + n_frac_bits, " << m_int_bits << ", " << is_neg_m_s << "> m_lut[n_segments_lut] = " << mstr << ";" << std::endl;
  outfile << "const ac_fixed<" << c_int_bits << " + n_frac_bits, " << c_int_bits << ", " << is_neg_c_s << "> c_lut[n_segments_lut] = " << cstr << ";" << std::endl;
  outfile << "const ac_fixed<" << x_min_int_bits << " + n_frac_bits, " << x_min_int_bits << ", " << is_neg_x_min_s << "> x_min_lut = " << ha_fxpt_type(x_min) << ";" << std::endl;
  outfile << "const ac_fixed<" << x_max_int_bits << " + n_frac_bits, " << x_max_int_bits << ", " << is_neg_x_max_s << "> x_max_lut = " << ha_fxpt_type(x_max) << ";" << std::endl;
  outfile << "const ac_fixed<" << p_c_int_bits << " + n_frac_bits, " << p_c_int_bits << ", false> sc_constant_lut = " << ha_fxpt_type(prop_constant) << ";" << std::endl;
  outfile.close();

  std::cout << std::endl;
  std::cout << __FILE__ << ", " << __LINE__ << ": Values are written, check \"" << filename << "\" for " << "the required ROM values" << std::endl;
  PRINT_EXPR(abs_error_max);
  PRINT_EXPR(input_error_max);
  PRINT_EXPR(nfrac_bits);
  PRINT_EXPR(rel_error_max);

  return 0;
}
