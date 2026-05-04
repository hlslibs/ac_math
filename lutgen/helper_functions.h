/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 2026.1                                              *
 *                                                                        *
 *  Release Date    : Wed Mar 11 20:39:39 PDT 2026                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 2026.1.1                                            *
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
#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

template <int size>
double pwl_new(const double x_in, const double (&m)[size], const double (&c)[size], const double prop_constant, const double x_min, const unsigned nsegments)
{
  unsigned index;
  double x_in_sc = (x_in - x_min) * prop_constant;
  if (x_in_sc < nsegments) {index = floor(x_in_sc);}
  else {index = nsegments - 1;}
  double y_out = m[index] * (x_in_sc - index) + c[index];
  return y_out;
}

// Function that returns max value and presence of negative elements in an array.
template <int size>
bool is_neg_max_array(const double (&array)[size], const unsigned nsegments, double &max_val)
{
  // This variable is set to true if even a single element is negative.
  bool is_neg = (array[0] < 0);
  max_val = abs(array[0]);
  for(unsigned i = 1; i < nsegments; i++) {
    if(array[i] < 0) { is_neg = true; }
    if(abs(array[i]) > max_val) { max_val = abs(array[i]); }
  }
  return is_neg;
}

// Make a number non-zero, useful for log calculations.
double make_non_zero(double input)
{
  if(input == 0) {input = 1;}
  return input;
}

int int_bits_calc(double val, bool S)
{
  return ceil(log2(make_non_zero(abs(val)))) + int(S) + 1;
}

// High-accuracy fixed point type.
typedef ac_fixed<128, 64, true> ha_fxpt_type;

// This function takes a double variable and performs an operation that mimics the quantization of the same double variable into an ac_fixed variable with "nfrac_bits" number of fractional
// bits and rounding turned on (AC_RND).
ha_fxpt_type o_ac_f(double input, int nfrac_bits)
{
  return (ha_fxpt_type)((double)rint(input * (1 << nfrac_bits)) * pow(2, (double)(-nfrac_bits)));
}

// This function takes a double variable and performs an operation that mimics the quantization of the same double variable into an ac_fixed variable with "nfrac_bits" number of fractional
// bits and rounding turned off (AC_TRN).
ha_fxpt_type o_ac_f_trn(double input, int nfrac_bits)
{
  return (ha_fxpt_type)((double)floor(input * (1 << nfrac_bits)) * pow(2, (double)(-nfrac_bits)));
}

// Find out whether a double number is an integer or not.
bool is_int(double &input) {
  return trunc(input) == input;
}

// Convert values in arrays to a fixed point precision, where nfrac_bits is the number of fractional bits.
template <int size>
void conv_arr_to_fxpt(double (&array)[size], int nfrac_bits)
{
  for (int i = 0; i < size; ++i) { array[i] = o_ac_f(array[i], nfrac_bits).to_double(); }
}

// Adjust m and c arrays to ensure monotonicity.
template <int size>
void mono_adjust(double (&m)[size], double (&c)[size])
{
  for (int i = 1; i < size; i++) {
    double y1_new = m[i - 1] + c[i - 1];
    double y2_new = m[i] + c[i];
    m[i] = y2_new - y1_new;
    c[i] = y1_new;
  }
}

template <int size, typename stream_type>
void dump_array_to_stream(double (&array)[size], stream_type &stream, int max_elems_per_line = 8) {
  bool separate_across_lines = (size > max_elems_per_line);

  for (int i = 0; i < size; i++) {
    if (i == 0) {
      stream << "{";
      if (separate_across_lines) { stream << std::endl << "  "; }
      stream << ha_fxpt_type(array[i]) << ", ";
    } else if (i == size - 1) {
      stream << ha_fxpt_type(array[i]);
      if (separate_across_lines) { stream << std::endl; }
      stream << "}";
    } else {
      bool start_new_line = ((i%max_elems_per_line) == (max_elems_per_line - 1));
      stream << ha_fxpt_type(array[i]) << ",";
      if (separate_across_lines && start_new_line) { stream << std::endl << "  "; }
      else { stream << " "; }
    }
  }
}

#endif // HELPER_FUNCTIONS_H
