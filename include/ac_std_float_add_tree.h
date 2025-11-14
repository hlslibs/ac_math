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
//*****************************************************************************************
// File: ac_std_float_add_tree.h
//
// Description: Fused adder tree implementations for ac_std_float, ac_ieee_float and 
//    ac::bfloat16 datatypes. All functions accept an array of inputs following the 
//    referenced datatypes and produce an output result in the same datatype as the inputs.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_std_float_add_tree.h>
//    using namespace ac_math;
//    
//    typedef ac_std_float<32, 8> input_type;
//    typedef ac_std_float<32, 8> output_type;
//    constexpr int N_ELEMS = 4;
//    
//    #pragma hls_design top
//    void project(
//      const input_type (&input)[N_ELEMS],
//      output_type &output
//    )
//    {
//      fadd_tree(input, output);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int, char **)
//    {
//      input_type input[N_ELEMS];
//      input[0] =  ac_std_float<32, 8>(1.5);
//      input[1] =  ac_std_float<32, 8>(2.5);
//      input[2] =  ac_std_float<32, 8>(0.5);
//      input[3] =  ac_std_float<32, 8>(-1.0);
//      output_type output;
//      
//      CCS_DESIGN(project)(input, output);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    This file uses C++ function overloading to target implementations
//    specific to each type of data. Attempting to call the function
//    with a type that is not implemented will result in a compile error.
//
//    The pass-by-pointer functions are differentiated from their
//    pass-by-reference counterparts with the "_ptr" suffix in their name.
//    Overloading is not used here in order to minimize ambiguity.
//
//*****************************************************************************************
#ifndef _INCLUDED_AC_STD_FLOAT_ADD_TREE_H
#define _INCLUDED_AC_STD_FLOAT_ADD_TREE_H

#include <ac_std_float.h>

namespace ac_math {

  template<int N>
  struct max_s {
    template<typename T>
    static T max(T *a) {
      T m0 = max_s<N/2>::max(a);
      T m1 = max_s<N-N/2>::max(a+N/2);

      return m0 > m1 ? m0 : m1;
    }
  };

  template<> 
  struct max_s<1> {
    template<typename T>
    static T max(T *a) {
      return a[0];
    }
  };

  template<int N, typename T>
  T max(T *a) {
    return max_s<N>::max(a);
  };


  // This function extracts the special values and the sign, exponent and mantissa fiels,
  // from the ac_int representation of an ac_std_float datatype.
  template<int W, int E, typename mu_t, typename e_t>
  inline void extract(const ac_int<W,true> &d, mu_t &m, e_t &e, bool &sign, bool &normal, bool &zero, bool &inf, bool &nan, bool biased_exp=false, bool no_subnormals=false) {

    static const int mant_bits = W - E - 1;
    static const int exp_bias = (1 << (E-1)) - 1;

    e = d.template slc<E>(mant_bits);
    bool exception = e == -1;
    normal = !!e | no_subnormals;
    m = d;
    bool m_zero = !m.template slc<mant_bits>(0);
    zero = (!e) & (no_subnormals | m_zero);
    m[mant_bits] = !!e;
    if(!biased_exp) {
      e -= exp_bias;
      e += !normal;
    }
    sign = d[W-1];
    inf = exception & m_zero;
    nan = exception & !m_zero;
  };

  template<int N, typename T>
  T perf_add(const T a[N]) {
    T r = T(0);

    _Pragma ("hls_unroll yes")
    AddOps: for (int i=0; i < N; i++) {
      r+= a[i];
    }
    return r;
  };

  template<
    ac_q_mode QR=AC_RND_CONV, 
    bool No_SubNormals=false, 
    int AccExtraLSBs=3,
    int W, int E,
    int N_ELEMS
  >
  void fadd_tree_core(
    const ac_std_float<W,E> (&x)[N_ELEMS],
    ac_std_float<W,E> &acc
  ) {
    ac_private::check_supported<QR>();

    static const int mant_bits = W - E - 1;
    static const int mu_bits = mant_bits + 1;
    static const int exp_bias = (1 << (E-1)) - 1;
    static const int max_exp = exp_bias;

    typedef ac_int<E,true> e_t;
    typedef ac_int<E,false> eu_t;
    typedef ac_int<mu_bits,false> mu_t;

    e_t op_e[N_ELEMS];
    eu_t op_e_b[N_ELEMS];
    mu_t op_mu[N_ELEMS];
    bool op_normal[N_ELEMS], op_sign[N_ELEMS], op_zero[N_ELEMS];
    bool op_inf[N_ELEMS], op_nan[N_ELEMS];

    // For each input extract the exponent, mantissa, sign and special values
    _Pragma ("hls_unroll yes")
    Extract: for (int i=0; i < N_ELEMS; i++) {
      extract<W,E,mu_t,e_t>(x[i].data(), op_mu[i], op_e[i], op_sign[i], op_normal[i], op_zero[i], op_inf[i], op_nan[i], true, No_SubNormals);
      op_e_b[i] = eu_t(op_e[i]);
    }
    
    // Find the maximum of the exponents
    eu_t op_e_b_max = max<N_ELEMS>(op_e_b);

    // Fix the MSB in cases of Subnormal or zero
    // and align mantisssa in case of subnormals
    _Pragma ("hls_unroll yes")
    SubNormMan: for (int i=0; i < N_ELEMS; i++) {
      if (!op_normal[i]) {
        op_mu[i] <<= 1;
        if (No_SubNormals)
          op_mu[i] = mu_t(0);
      } else if (No_SubNormals & op_zero[i]) {
        op_mu[i] = 0;
      }
    }

    bool neg_inf = false;
    bool pos_inf = false;
    bool is_inf = false;
    bool is_nan = false;
    _Pragma ("hls_unroll yes")
    InfNan: for (int i=0; i < N_ELEMS; i++) {
      is_inf |= op_inf[i];
      is_nan |= op_nan[i];

      // Check if there exist both a +inf and a -inf input. 
      // This would result in a NaN result
      neg_inf |= op_inf[i] & op_sign[i];
      pos_inf |= op_inf[i] & !op_sign[i];

    }
    is_nan |= (neg_inf & pos_inf);


    // saturate e_dif (if bigger than mu_bits, result will be the same as mu_bits+1 because shift will be bigger than number of bits)
    typedef ac_int<ac::nbits<mu_bits+1>::val,false> e_dif_sat_t;
    e_dif_sat_t e_dif_sat[N_ELEMS];

    // For each exponent compute the difference to the maximum. 
    // The value will be used to align the mantissas.
    eu_t e_dif[N_ELEMS];
    _Pragma ("hls_unroll yes")
    ExpDiff: for (int i=0; i < N_ELEMS; i++) {
      e_dif[i] = op_e_b_max - op_e_b[i];
      e_dif_sat[i] = ac_fixed<ac::nbits<mu_bits+1>::val,ac::nbits<mu_bits+1>::val,false,AC_TRN,AC_SAT>(e_dif[i]).to_ac_int();
    }

    const int logN_MSBs = ac::nbits<N_ELEMS-1>::val; //MSBs required for safe addition without overflow 
    const int op_bits = mu_bits + logN_MSBs; 
    typedef ac_fixed<   op_bits + AccExtraLSBs,    op_bits, false, QR> uadd_GRS_t;
    typedef ac_fixed<1+ op_bits + AccExtraLSBs, 1+ op_bits, true,  QR> add_GRS_t;

    uadd_GRS_t add_op[N_ELEMS];
    add_GRS_t to_add[N_ELEMS];

    _Pragma ("hls_unroll yes")
    Allign: for (int i=0; i < N_ELEMS; i++) {
      add_op[i] = uadd_GRS_t(op_mu[i]);

      // Compute the sticky bit for each operand by ORing the bits that will be shifted out
      // The +(AccExtraLSBs-1) bits on mask work as a guard in the case of e_dif_sat > mu_bits
      // in order to keep the mask correctly aligned due to (>>=(AccExtraLSBs-1))
      ac_int<mu_bits+1+((AccExtraLSBs == 0) ? 0 : AccExtraLSBs-1),false> sticky_bit_mask = -1; 
      sticky_bit_mask <<= e_dif_sat[i];
      sticky_bit_mask = ~sticky_bit_mask;
      sticky_bit_mask >>= (AccExtraLSBs == 0) ? 0 : AccExtraLSBs-1;
      bool sticky_bit = !!(op_mu[i] & sticky_bit_mask);

      // Shift the smaller mantissas according to the exponent difference 
      add_op[i] >>= e_dif_sat[i];
      add_op[i][0] = sticky_bit;

      // Convert to two's complement representation
      to_add[i] = (op_sign[i]) ? add_GRS_t(-add_op[i]) : add_GRS_t(add_op[i]);
    }

    // Add the mantissas
    add_GRS_t res_mant_signed = perf_add<N_ELEMS, add_GRS_t>(to_add);

    // The sign of the result is taken from the result's MSB
    bool r_sign = res_mant_signed[op_bits + AccExtraLSBs]; 

    // Invert the result if it has a negative value
    uadd_GRS_t res_mant = (r_sign) ? uadd_GRS_t(-res_mant_signed) : uadd_GRS_t(res_mant_signed); 

    // Find the leading zeros to normalize
    bool all_sign;
    ac_int<ac::nbits<op_bits+AccExtraLSBs>::val,false> ls = (res_mant.template slc<op_bits+AccExtraLSBs>(0)).leading_sign(all_sign); 
    bool r_zero = all_sign; // we know res_mant >= 0, so all_sign means equal zero

    bool r_normal = true;
    bool r_inf = false;
    bool r_nan = false;
  
    // compute left shift and res exponent
    int resexp = op_e_b_max - ls + logN_MSBs-1;
    bool incr = 0;
    if (resexp < 0) {
      ls = op_e_b_max;
      op_e_b_max = -1;
      r_normal = false;
    } else {
      op_e_b_max = resexp;
      r_normal = true;
    }
    res_mant <<= ls;
    res_mant >>= 1;

    // Round the result
    typedef ac_fixed<mu_bits+1,op_bits,false,QR> add_rounded_t;
    add_rounded_t res_rounded = res_mant;
    
    // Check for overflow in rounding. In the case of truncation 
    // the rounding can't result in ovwerflow 
    if (res_rounded[mu_bits] & (QR != AC_TRN_ZERO)) { 
      res_rounded >>= 1;
      incr = 1;
      r_normal |= resexp == -1;
    }

    op_e_b_max += incr + 1;
    r_inf = op_e_b_max > max_exp + exp_bias + 1;
    bool exp_max = eu_t(op_e_b_max) == (max_exp + exp_bias+1);
    r_zero |= No_SubNormals & !r_normal;

    ac_int<mant_bits,false> m_r;
    m_r = res_rounded.template slc<mant_bits>(0);

    // special case when AC_TRN_ZERO : infinity is replaced by max value
    if ((r_inf|exp_max) & QR==AC_TRN_ZERO) {
      op_e_b_max = max_exp + exp_bias; // saturate res
      r_inf = false;
      m_r |= ac_int<1,true>(-1);  // saturate (set all bits to 1)
    }

    r_inf |= is_inf;
    r_nan |= is_nan;

    // In case of infinity, select the correct sign
    r_sign = (neg_inf | pos_inf) ? (neg_inf & !pos_inf) : r_sign;
    
    // compute flags and assign result
    bool exception = r_nan | r_inf;
    ac_int<E,true> e_r = op_e_b_max;
    if(exception | r_zero)
      e_r = ac_int<E,true>(-1)*exception;
    exception |= (exp_max & (QR!=AC_TRN_ZERO));
    if(exception | r_zero)
      m_r = ac_int<mant_bits,false>(-1)*r_nan;
    ac_int<W,true> d_r = m_r;
    d_r.set_slc(mant_bits, e_r);
    d_r[W-1] = r_sign;
    ac_std_float<W,E> r;
    r.set_data(d_r);
    acc = r;
  }

  //=========================================================================
  // Function: fadd_tree (ac_std_float inputs, ac_std_float output)
  //
  // Description:
  //    "Fused" implementation of adder tree, where the input mantissas are 
  //    aligned all at once with respect to the maximum exponent value and
  //    then added via an unrolled loop. The addition of the mantissas is
  //    performed with two's complement fixed point values. The result of the
  //    addition is then normalized and rounded. Finally the output is packed 
  //    to ac_std_float by considering any possible special values.
  //
  // Usage:
  //    A sample testbench and its implementation look like
  //    this:
  //
  //    #include <ac_std_float_add_tree.h>
  //    using namespace ac_math;
  //    
  //    typedef ac_std_float<32, 8> input_type;
  //    typedef ac_std_float<32, 8> output_type;
  //    constexpr int N_ELEMS = 4;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const input_type (&input)[N_ELEMS],
  //      output_type &output
  //    )
  //    {
  //      fadd_tree(input, output);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //
  //    CCS_MAIN(int, char **)
  //    {
  //      input_type input[N_ELEMS];
  //      input[0] =  ac_std_float<32, 8>(1.5);
  //      input[1] =  ac_std_float<32, 8>(2.5);
  //      input[2] =  ac_std_float<32, 8>(0.5);
  //      input[3] =  ac_std_float<32, 8>(-1.0);
  //      output_type output;
  //      
  //      CCS_DESIGN(project)(input, output);
  //      
  //      CCS_RETURN (0);
  //    }
  //    #endif
  //
  // Notes:
  //
  //-------------------------------------------------------------------------
  template<
    ac_q_mode QR=AC_RND_CONV, 
    bool No_SubNormals=false, 
    int AccExtraLSBs=3,
    int W, int E,
    int N_ELEMS
  >
  void fadd_tree(
    const ac_std_float<W,E> (&x)[N_ELEMS],
    ac_std_float<W,E> &acc
  ) {

    fadd_tree_core<QR, No_SubNormals, AccExtraLSBs>(x, acc);
  }
  
  //=========================================================================
  // Function: fadd_tree_ptr (ac_ieee_float inputs, ac_ieee_float output)
  //
  // Description:
  //    A variant of the fadd_tree function defined for the ac_ieee_float
  //    datatype.
  //
  // Usage: 
  //    A sample testbench and its implementation look like
  //    this:
  //
  //    #include <ac_std_float_add_tree.h>
  //    using namespace ac_math;
  //    
  //    typedef ac_ieee_float<binary32> input_type;
  //    typedef ac_ieee_float<binary32> output_type;
  //    constexpr int N_ELEMS = 4;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const input_type (&input)[N_ELEMS],
  //      output_type &output
  //    )
  //    {
  //      fadd_tree(input, output);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //
  //    CCS_MAIN(int, char **)
  //    {
  //      input_type input[N_ELEMS];
  //      input[0] =  ac_ieee_float<binary32>(1.5);
  //      input[1] =  ac_ieee_float<binary32>(2.5);
  //      input[2] =  ac_ieee_float<binary32>(0.5);
  //      input[3] =  ac_ieee_float<binary32>(-1.0);
  //      output_type output;
  //      
  //      CCS_DESIGN(project)(input, output);
  //      
  //      CCS_RETURN (0);
  //    }
  //    #endif
  //
  // Notes:
  //    helper_t datatype is the ac_std_float equivalent defined inside the
  //    the ac_ieee_float class.
  //
  //-------------------------------------------------------------------------
  template<
    ac_q_mode QR=AC_RND_CONV, 
    bool No_SubNormals=false, 
    int AccExtraLSBs=3,
    ac_ieee_float_format Format,
    int N_ELEMS
  >
  void fadd_tree(
    const ac_ieee_float<Format> (&x)[N_ELEMS],
    ac_ieee_float<Format> &acc
  ) {
    typename ac_ieee_float<Format>::helper_t cpyx[N_ELEMS];
    typename ac_ieee_float<Format>::helper_t cpyacc;

    _Pragma ("hls_unroll yes")
    CONV_TO_STD: for (int k = 0; k < N_ELEMS; k++) {
      cpyx[k] = x[k].to_ac_std_float();
    }

    fadd_tree_core<QR, No_SubNormals, AccExtraLSBs>(cpyx, cpyacc);

    acc = ac_ieee_float<Format>(cpyacc);
  }

  //=========================================================================
  // Function: fadd_tree_ptr (ac::bfloat16 inputs, ac::bfloat16 output)
  //
  // Description:
  //    A variant of the fadd_tree function defined for the ac::bfloat16
  //    datatype.
  //
  // Usage: 
  //    A sample testbench and its implementation look like
  //    this:
  //
  //    #include <ac_std_float_add_tree.h>
  //    using namespace ac_math;
  //    
  //    typedef ac::bfloat16 input_type;
  //    typedef ac::bfloat16 output_type;
  //    constexpr int N_ELEMS = 4;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const input_type (&input)[N_ELEMS],
  //      output_type &output
  //    )
  //    {
  //      fadd_tree(input, output);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //
  //    CCS_MAIN(int, char **)
  //    {
  //      input_type input[N_ELEMS];
  //      input[0] =  ac::bfloat16(1.5);
  //      input[1] =  ac::bfloat16(2.5);
  //      input[2] =  ac::bfloat16(0.5);
  //      input[3] =  ac::bfloat16(-1.0);
  //      output_type output;
  //      
  //      CCS_DESIGN(project)(input, output);
  //      
  //      CCS_RETURN (0);
  //    }
  //    #endif
  //
  // Notes:
  //    helper_t datatype is the ac_std_float equivalent defined inside the
  //    the ac::bfloat16 class. The default rounding mode on this case is 
  //    AC_TRN_ZERO to comply with the default rounding mode of all ac::bfloat16
  //    operators in the ac_std_float.h library.
  //
  //-------------------------------------------------------------------------
  template<
    ac_q_mode QR=AC_TRN_ZERO, 
    bool No_SubNormals=false, 
    int AccExtraLSBs=3,
    int N_ELEMS
  >
  void fadd_tree(
    const ac::bfloat16 (&x)[N_ELEMS],
    ac::bfloat16 &acc
  ) {
    ac::bfloat16::helper_t cpyx[N_ELEMS];
    ac::bfloat16::helper_t cpyacc;

    _Pragma ("hls_unroll yes")
    CONV_TO_STD: for (int k = 0; k < N_ELEMS; k++) {
      cpyx[k] = x[k].to_ac_std_float();
    }

    fadd_tree_core<QR, No_SubNormals, AccExtraLSBs>(cpyx, cpyacc);

    acc = ac::bfloat16(cpyacc);
  }

  //=========================================================================
  // Function: fadd_tree_ptr (ac_std_float inputs, ac_std_float output)
  //
  // Description:
  //    A variant of the equivalent fadd_tree function defined above,
  //    with the only major difference being that the input array is passed
  //    as a pointer rather than a reference. Since size information for the
  //    pass-by-pointer version is lost, the user must explicitly specify
  //    N_ELEMS in the function call.
  //
  // Usage: 
  //    A sample testbench and its implementation look like
  //    this:
  //
  //    #include <ac_std_float_add_tree.h>
  //    using namespace ac_math;
  //    
  //    typedef ac_std_float<32, 8> input_type;
  //    typedef ac_std_float<32, 8> output_type;
  //    constexpr int N_ELEMS = 4;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const input_type input[N_ELEMS],
  //      output_type &output
  //    )
  //    {
  //      fadd_tree_ptr<N_ELEMS>(input, output);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //
  //    CCS_MAIN(int, char **)
  //    {
  //      input_type input[N_ELEMS];
  //      input[0] =  ac_std_float<32, 8>(1.5);
  //      input[1] =  ac_std_float<32, 8>(2.5);
  //      input[2] =  ac_std_float<32, 8>(0.5);
  //      input[3] =  ac_std_float<32, 8>(-1.0);
  //      output_type output;
  //      
  //      CCS_DESIGN(project)(input, output);
  //      
  //      CCS_RETURN (0);
  //    }
  //    #endif
  //
  // Notes:
  //    The inputs are copied to a temporary array in an unrolled loop, and
  //    this array is then passed to the pass-by-reference function.
  //
  //-------------------------------------------------------------------------
  template<
    int N_ELEMS,
    ac_q_mode QR=AC_RND_CONV, 
    bool No_SubNormals=false, 
    int AccExtraLSBs=3,
    int W, int E
  >
  void fadd_tree_ptr(
    const ac_std_float<W,E> x[N_ELEMS],
    ac_std_float<W,E> &acc
  ) {
    ac_std_float<W,E> cpyx[N_ELEMS];

    _Pragma ("hls_unroll yes")
    CPY_IN_ARR: for (int k = 0; k < N_ELEMS; k++) {
      cpyx[k] = x[k];
    }

    fadd_tree_core<QR, No_SubNormals, AccExtraLSBs>(cpyx, acc);
  }

  //=========================================================================
  // Function: fadd_tree_ptr (ac_ieee_float inputs, ac_ieee_float output)
  //
  // Description:
  //    A variant of the fadd_tree_ptr function defined for the ac_ieee_float
  //    datatype.
  //
  // Usage: 
  //    A sample testbench and its implementation look like
  //    this:
  //
  //    #include <ac_std_float_add_tree.h>
  //    using namespace ac_math;
  //    
  //    typedef ac_ieee_float<binary32> input_type;
  //    typedef ac_ieee_float<binary32> output_type;
  //    constexpr int N_ELEMS = 4;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const input_type input[N_ELEMS],
  //      output_type &output
  //    )
  //    {
  //      fadd_tree_ptr<N_ELEMS>(input, output);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //
  //    CCS_MAIN(int, char **)
  //    {
  //      input_type input[N_ELEMS];
  //      input[0] =  ac_ieee_float<binary32>(1.5);
  //      input[1] =  ac_ieee_float<binary32>(2.5);
  //      input[2] =  ac_ieee_float<binary32>(0.5);
  //      input[3] =  ac_ieee_float<binary32>(-1.0);
  //      output_type output;
  //      
  //      CCS_DESIGN(project)(input, output);
  //      
  //      CCS_RETURN (0);
  //    }
  //    #endif
  //
  // Notes:
  //    helper_t datatype is the ac_std_float equivalent defined inside the
  //    the ac_ieee_float class.
  //
  //-------------------------------------------------------------------------
  template<
    int N_ELEMS,
    ac_q_mode QR=AC_RND_CONV, 
    bool No_SubNormals=false, 
    int AccExtraLSBs=3,
    ac_ieee_float_format Format
  >
  void fadd_tree_ptr(
    const ac_ieee_float<Format> x[N_ELEMS],
    ac_ieee_float<Format> &acc
  ) {
    typename ac_ieee_float<Format>::helper_t cpyx[N_ELEMS];
    typename ac_ieee_float<Format>::helper_t cpyacc;

    _Pragma ("hls_unroll yes")
    CONV_TO_STD: for (int k = 0; k < N_ELEMS; k++) {
      cpyx[k] = x[k].to_ac_std_float();
    }

    fadd_tree_core<QR, No_SubNormals, AccExtraLSBs>(cpyx, cpyacc);

    acc = ac_ieee_float<Format>(cpyacc);
  }

  //=========================================================================
  // Function: fadd_tree_ptr (ac::bfloat16 inputs, ac::bfloat16 output)
  //
  // Description:
  //    A variant of the fadd_tree_ptr function defined for the ac::bfloat16
  //    datatype.
  //
  // Usage: 
  //    A sample testbench and its implementation look like
  //    this:
  //
  //    #include <ac_std_float_add_tree.h>
  //    using namespace ac_math;
  //    
  //    typedef ac::bfloat16 input_type;
  //    typedef ac::bfloat16 output_type;
  //    constexpr int N_ELEMS = 4;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const input_type input[N_ELEMS],
  //      output_type &output
  //    )
  //    {
  //      fadd_tree_ptr<N_ELEMS>(input, output);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //
  //    CCS_MAIN(int, char **)
  //    {
  //      input_type input[N_ELEMS];
  //      input[0] =  ac::bfloat16(1.5);
  //      input[1] =  ac::bfloat16(2.5);
  //      input[2] =  ac::bfloat16(0.5);
  //      input[3] =  ac::bfloat16(-1.0);
  //      output_type output;
  //      
  //      CCS_DESIGN(project)(input, output);
  //      
  //      CCS_RETURN (0);
  //    }
  //    #endif
  //
  // Notes:
  //    helper_t datatype is the ac_std_float equivalent defined inside the
  //    the ac::bfloat16 class. The default rounding mode on this case is 
  //    AC_TRN_ZERO to comply with the default rounding mode of all ac::bfloat16
  //    operators in the ac_std_float.h library.
  //
  //-------------------------------------------------------------------------
  template<
    int N_ELEMS,
    ac_q_mode QR=AC_TRN_ZERO, 
    bool No_SubNormals=false, 
    int AccExtraLSBs=3
  >
  void fadd_tree_ptr(
    const ac::bfloat16 x[N_ELEMS],
    ac::bfloat16 &acc
  ) {
    ac::bfloat16::helper_t cpyx[N_ELEMS];
    ac::bfloat16::helper_t cpyacc;

    _Pragma ("hls_unroll yes")
    CONV_TO_STD: for (int k = 0; k < N_ELEMS; k++) {
      cpyx[k] = x[k].to_ac_std_float();
    }

    fadd_tree_core<QR, No_SubNormals, AccExtraLSBs>(cpyx, cpyacc);

    acc = ac::bfloat16(cpyacc);
  }
 
}

#endif

