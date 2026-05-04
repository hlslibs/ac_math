/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 2026.1                                              *
 *                                                                        *
 *  Release Date    : Wed Feb 11 11:04:12 PST 2026                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 2026.1.0                                            *
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
//***************************************************************************
// File: ac_random.h
//
// Description: This file provides random number generators, with and
//  without min/max constraints. These functions are used in the
//  automatically generated testbench.
//
// Revision History:
//    4.0.0  - CAT-41533 - Added constrained randomization.
//    3.6.0  - Added ac_float support.
//
//***************************************************************************

#ifndef _INCLUDED_AC_RANDOM_H_
#define _INCLUDED_AC_RANDOM_H_

#include <stdlib.h>
#include <limits.h>

// Check for macro definitions that will conflict with template parameter names in this file
#if defined(T)
#define T 0
#error The macro name 'T' conflicts with a Template Parameter name.
#error The compiler should have produced a warning about redefinition of 'T' giving the location of the previous definition.
#endif
#if defined(size)
#define size 0
#error The macro name 'size' conflicts with a Template Parameter name.
#error The compiler should have produced a warning about redefinition of 'size' giving the location of the previous definition.
#endif
#if defined(SIZE)
#define SIZE 0
#error The macro name 'SIZE' conflicts with a Template Parameter name.
#error The compiler should have produced a warning about redefinition of 'SIZE' giving the location of the previous definition.
#endif
#if defined(SMALL)
#define SMALL 0
#error The macro name 'SMALL' conflicts with a Template Parameter name.
#error The compiler should have produced a warning about redefinition of 'SMALL' giving the location of the previous definition.
#endif
#if defined(W)
#define W 0
#error The macro name 'W' conflicts with a Template Parameter name.
#error The compiler should have produced a warning about redefinition of 'W' giving the location of the previous definition.
#endif
#if defined(S)
#define S 0
#error The macro name 'S' conflicts with a Template Parameter name.
#error The compiler should have produced a warning about redefinition of 'S' giving the location of the previous definition.
#endif
#if defined(I)
#define I 0
#error The macro name 'I' conflicts with a Template Parameter name.
#error The compiler should have produced a warning about redefinition of 'I' giving the location of the previous definition.
#endif
#if defined(Q)
#define Q 0
#error The macro name 'Q' conflicts with a Template Parameter name.
#error The compiler should have produced a warning about redefinition of 'Q' giving the location of the previous definition.
#endif
#if defined(O)
#define O 0
#error The macro name 'O' conflicts with a Template Parameter name.
#error The compiler should have produced a warning about redefinition of 'O' giving the location of the previous definition.
#endif

template <typename T, int size>
struct ac_random_c_builtin_s {
  void operator()(T &v);
};

template <typename T, int size>
void ac_random_c_builtin_s<T, size>::operator()(T &v)
{
  const int long mask = 0xffff;
  v = (rand() & mask);
  for (int i = 1; i < size; ++i) {
    v <<= 16;
    v |= (rand() & mask);
  }
}

template <typename T>
struct ac_random_c_builtin_s<T,0> {
  void operator()(T &v) { v = rand(); }
};

template <typename T>
struct ac_random_c_builtin_s<T,1> {
  void operator()(T &v) { v = rand(); }
};

template <typename T>
inline void ac_random_c_builtin(T &v)
{
  enum { size = (sizeof(T) * CHAR_BIT) / 16 };
  ac_random_c_builtin_s<T,size>()(v);
}

// ======================================================================
// Native C Type Support
//

inline void ac_random(long long int &v) { ac_random_c_builtin(v); }
inline void ac_random(long long unsigned &v) { ac_random_c_builtin(v); }
inline void ac_random(long int &v) { ac_random_c_builtin(v); }
inline void ac_random(long unsigned &v) { ac_random_c_builtin(v); }
inline void ac_random(int &v) { ac_random_c_builtin(v); }
inline void ac_random(unsigned &v) { ac_random_c_builtin(v); }
inline void ac_random(short int &v) { ac_random_c_builtin(v); }
inline void ac_random(short unsigned &v) { ac_random_c_builtin(v); }
inline void ac_random(char &v) { ac_random_c_builtin(v); }
inline void ac_random(signed char &v) { ac_random_c_builtin(v); }
inline void ac_random(unsigned char &v) { ac_random_c_builtin(v); }
inline void ac_random(bool &v) { v = rand() & 1; }

#if (defined(SYSTEMC_INCLUDED) || defined(SYSTEMC_H)) && !defined(MC_SYSTEMC)
#define AC_SYSTEMC
#endif

// ======================================================================
// SystemC Type Support
//

#if defined(AC_SYSTEMC) && !defined(AC_RANDOM_H_INC_SYSTEMC_INCLUDED)
#define AC_RANDOM_H_INC_SYSTEMC_INCLUDED
template <typename T>
void ac_systemc_builtin_it(T &v, const int size)
{
  const int long mask = 0xffff;
  int lo = 0;
  int hi = 15;
  for (; hi < size; lo += 16, hi += 16) {
    v.range(hi, lo) = (rand() & mask);
  }
  if (lo < size) {
    v.range(size - 1, lo) = (rand() & mask);
  }
}

inline void ac_random(sc_int_base &v) { ac_systemc_builtin_it(v, v.length()); }
inline void ac_random(sc_uint_base &v) { ac_systemc_builtin_it(v, v.length()); }
inline void ac_random(sc_signed &v) { ac_systemc_builtin_it(v, v.length()); }
inline void ac_random(sc_unsigned &v) { ac_systemc_builtin_it(v, v.length()); }
template <int Tlen>
inline void ac_random(sc_bv<Tlen> &v) { ac_systemc_builtin_it(v, Tlen); }
template <int Tlen>
inline void ac_random(sc_lv<Tlen> &v) { ac_systemc_builtin_it(v, Tlen); }
#endif // AC_RANDOM_H_INC_SYSTEMC_INCLUDED

#if defined(SC_INCLUDE_FX) && defined(AC_SYSTEMC) && !defined(AC_RANDOM_H_INC_SC_INCLUDE_FX)
#define AC_RANDOM_H_INC_SC_INCLUDE_FX

void ac_random(sc_fxnum &v) {  ac_systemc_builtin_it(v, v.wl()); }

#endif // AC_RANDOM_H_INC_SC_INCLUDE_FX

// ======================================================================
// AC Datatype "ac_int" Support
//

#if defined(__AC_INT_H) && !defined(AC_RANDOM_H_INC__AC_INT_H)
#define AC_RANDOM_H_INC__AC_INT_H

template <class T, int SIZE, bool SMALL>
struct ac_random_ac_s {
  void operator()(T &v);
};

template <class T, int SIZE, bool SMALL>
void ac_random_ac_s<T,SIZE,SMALL>::operator()(T &v)
{
  const ac_int<6> ac_16 = 16;
  ac_int<16> rv = rand();
  v = 0;
  v.set_slc(0, rv);
  for (int i = 16; i < SIZE; i += 16) {
    v <<= ac_16;
    rv = rand();
    v.set_slc(0, rv);
  }
}

template <class T, int SIZE>
struct ac_random_ac_s<T,SIZE,true> {
  void operator()(T &v) { v.set_slc(0, ac_int<SIZE>(rand())); }
};

template <int W, bool S>
inline void ac_random(ac_int<W,S> &v)
{
  // This variable is static so that if we're using VRA, it's only constructed
  // once, thereby reducing the stacktracing overhead that's incurred every time
  // an ac_int variable is constructed while using VRA. This is particularly
  // beneficial if, say, ac_random is called inside a loop with many iterations.
  static ac_int<W, S> v_temp;
  ac_random_ac_s<ac_int<W,S>, W, W<=16>()(v_temp);
  v = v_temp;
}

// Called if a min/max random value is specified for constrained randomization.
template <int W, bool S>
inline void ac_random(
  ac_int<W,S> &v, // Range of output : [min, max]
  const ac_int<W,S> &min,
  const ac_int<W,S> &max
) {
  AC_ASSERT(min <= max, "Min. limit must not exceed max. limit.");
  typedef ac_int<W, false> v_temp_type;
  // The random number generation operates with static variables to reduce VRA runtime,
  // similar to what's done during regular ac_random operation.
  static v_temp_type v_temp;
  ac_random_ac_s<v_temp_type, v_temp.width, v_temp.width<=16>()(v_temp);
  // Constrain the final randomized output by performing
  // a mod operation first and then adding an offset.
  ac_int<W + 1, false> mod_den = max - min + 1;
  v = (v_temp%mod_den) + min;
}

#endif // AC_RANDOM_H_INC__AC_INT_H

// ======================================================================
// AC Datatype "ac_fixed" Support
//

#if defined(__AC_FIXED_H) && !defined(AC_RANDOM_H_INC__AC_FIXED_H)
#define AC_RANDOM_H_INC__AC_FIXED_H

template <int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
inline void ac_random(ac_fixed<W,I,S,Q,O> &v)
{
  // This variable is static so that if we're using VRA, it's only constructed
  // once, thereby reducing the stacktracing overhead that's incurred every time
  // an ac_int variable is constructed while using VRA. This is particularly
  // beneficial if, say, ac_random is called inside a loop with many iterations.
  static ac_fixed<W, I, S, Q, O> v_temp;
  ac_random_ac_s<ac_fixed<W,I,S,Q,O>, W, W<=16>()(v_temp);
  v = v_temp;
}

// Called if a min/max random value is specified for constrained randomization.
template <int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
inline void ac_random(
  ac_fixed<W,I,S,Q,O> &v, // Range of output : [min, max]
  const ac_fixed<W,I,S,Q,O> &min,
  const ac_fixed<W,I,S,Q,O> &max
) {
  AC_ASSERT(min <= max, "Min. limit must not exceed max. limit.");
  // The random number generation operates with static variables to reduce VRA runtime,
  // similar to what's done during regular ac_random operation.
  static ac_int<W, S> v_temp;
  static ac_int<W, S> min_int, max_int;
  min_int.set_slc(0, min.template slc<W>(0));
  max_int.set_slc(0, max.template slc<W>(0));
  ac_random(v_temp, min_int, max_int); // Call constrained ac_int RNG.
  v.set_slc(0, v_temp.template slc<W>(0));
  #ifdef AC_FIXED_VRA
  v = v; // This otherwise redundant assignment is needed for VRA tracking.
  #endif
}

#endif // AC_RANDOM_H_INC__AC_FIXED_H

// ======================================================================
// AC Datatype "ac_float" Support
//

#if defined(__AC_FLOAT_H) && !defined(AC_RANDOM_H_INC__AC_FLOAT_H)
#define AC_RANDOM_H_INC__AC_FLOAT_H

// Called if a min/max random value is specified for constrained randomization.
template <int W, int I, int E, ac_q_mode Q>
inline void ac_random(ac_float<W, I, E, Q> &v)
{
  // As explained earlier, the internal ac_fixed and ac_int variable
  // are static to reduce VRA runtime penalties.
  static ac_fixed<W, I, true> v_mant;
  static ac_int<E, true> v_exp;
  // Use ac_fixed and ac_int ac_random versions for mantissa and
  // exponent, respectively.
  ac_random(v_mant);
  ac_random(v_exp);
  ac_float<W, I, E, Q> v_temp(v_mant, v_exp, true);
  
  v = v_temp;
}

#if (defined(__GNUC__) && (__cplusplus >= 201103L)) || (defined(_MSC_VER) && (_MSC_VER >= 1920)) || defined(__EDG__)
#include <limits>
#include <random>

// Called if a min/max random value is specified for constrained randomization.
//
// Since constrained ac_float randomization uses the std::uniform_real_distribution
// library which was made available since C++11, the user must make sure they use
// C++11 or a later standard.
//
// Unlike constrained randomization for ac_int and ac_fixed, the output range for this
// random function does not included the max. specified.
template <int W, int I, int E, ac_q_mode Q>
inline void ac_random(
  ac_float<W, I, E, Q> &v, // Range of output : [min, max)
  const ac_float<W, I, E, Q> &min,
  const ac_float<W, I, E, Q> &max
) {
  AC_ASSERT(min < max, "Min. limit must be lesser than max. limit.");
  ac_float<54, 2, 11> max_lim_d = std::numeric_limits<double>::max();
  ac_float<54, 2, 11> min_lim_d = std::numeric_limits<double>::lowest();
  ac_float<W + 1, I, E + 1> mm_diff = max;
  mm_diff -= min;
  double min_d = min.to_double();
  double max_d = max.to_double();

  bool mmd_within_d_lims = (mm_diff <= max_lim_d);
  bool min_max_within_d_lims = (min >= min_lim_d) && (max <= max_lim_d);
  // Are we using the denormalization inside the C++ library function?
  bool use_clib_denorm = mmd_within_d_lims && min_max_within_d_lims && (max_d != min_d);

  ac_float<W, I, E> v_temp;

  if (use_clib_denorm) {
    // If the min./max. are small enough that the standard library can handle the
    // denormalization without overflowing, we let the standard library do it.
    static std::uniform_real_distribution<double> distribution(min_d, max_d);
    static std::default_random_engine generator;
    v_temp = distribution(generator);
  } else {
    // If the min./max. too large for the internal library to denormalize
    // without overflowing, we obtain the normalized random output {range: [0, 1)}
    // and do the denormalization with ac_float types which are large enough to
    // not overflow.
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::default_random_engine generator;
    ac_float<54, 2, 11> norm_random = distribution(generator);
    constexpr int temp_F = AC_MAX(W - I + 1, 52);
    constexpr int temp_I = AC_MAX(I, 2);
    constexpr int temp_W = temp_F + temp_I;
    constexpr int temp_E = AC_MAX(E + 1, 11);
    // Denormalize random output, following this formula: norm_random*(max - min) + min
    ac_float<temp_W, temp_I, temp_E> denorm_random = norm_random;
    denorm_random *= mm_diff;
    denorm_random += min;
    ac_float<W, I, E> v_temp_2(denorm_random.mantissa(), denorm_random.exp(), true);
    v_temp = v_temp_2;
  }

  // There are buggy cases where the output of std::uniform_real_distribution is equal to the max
  // limit specified. Out of an abundance of caution, the condition below sets v_temp to a value
  // slightly lesser than max if it isn't lesser than max already.
  if (v_temp >= max) {
    ac_fixed<W, I, true> mant_q = value<AC_VAL_QUANTUM>(ac_fixed<W, I, true>{});
    ac_fixed<W, I, true> max_allowed_mant = max.mantissa() - mant_q;
    ac_float<W, I, E> max_allowed(max_allowed_mant, max.exp(), true);
    v_temp = max_allowed;
  }

  v = v_temp;
}
#endif

#endif // AC_RANDOM_H_INC__AC_FLOAT_H

// ======================================================================
// AC Datatype "ac_channel" Support
//

#if defined(__AC_CHANNEL_H) && !defined(AC_RANDOM_H_INC__AC_CHANNEL_H)
#define AC_RANDOM_H_INC__AC_CHANNEL_H

// Used to extract the sub type of the channel
template <class Tchantype> struct ac_extract_chan_subtype;
template <class Tsubtype>
struct ac_extract_chan_subtype<ac_channel<Tsubtype> > {
  typedef Tsubtype chan_subtype;
};

#endif

#endif

