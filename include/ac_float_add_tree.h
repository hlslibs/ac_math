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
//*****************************************************************************************
// File: ac_float_add_tree.h
//
// Description: Adder tree implementations for ac_float datatypes. All functions
//    accept an array of ac_float inputs, and can produce ac_float or ac_fixed outputs.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_float_add_tree.h>
//    using namespace ac_math;
//    
//    typedef ac_float<25, 2, 8, AC_TRN> input_type;
//    typedef ac_fixed<64, 32, true, AC_TRN, AC_SAT> output_type;
//    constexpr int N_ELEMS = 4;
//    
//    #pragma hls_design top
//    void project(
//      const input_type (&input)[N_ELEMS],
//      output_type &output
//    )
//    {
//      add_tree(input, output);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int, char **)
//    {
//      input_type input[N_ELEMS];
//      input[0] =  1.5;
//      input[1] =  2.5;
//      input[2] =  0.5;
//      input[3] = -1.0;
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
// Revision History:
//    3.6.0  - Added usage examples.
//           - [CAT-36472] Added input array size to parameter list for *_ptr functions.
//    3.5.0  - [CAT-26464] Added library to ac_math subproject.
//
//*****************************************************************************************

#ifndef _INCLUDED_AC_FLOAT_ADD_TREE_H_
#define _INCLUDED_AC_FLOAT_ADD_TREE_H_

#include <ac_float.h>
#include <ac_math/ac_shift.h>

namespace ac_math {
  // This function sets the exponent value to the minimum value if the mantissa 
  // is zero and if ZeroHandling is set to true.
  template<bool ZeroHandling, int W1, int W2, int I>
  ac_int<W1, true> zero_handle_exp(const ac_int<W1, true> &exp, const ac_fixed<W2, I, true> &mant) {
    ac_int<W1, true> out_exp;
    
    if (ZeroHandling) {
      bool mant_zero = !mant;
      ac_int<W1 - 1, false> lsb_mask = -int(!mant_zero);
      out_exp = (exp.template slc<W1 - 1>(0)) & lsb_mask;
      out_exp[W1 - 1] = mant_zero | exp[W1 - 1];
    } else {
      out_exp = exp;
    }
    
    return out_exp;
  }
  
  // This function converts mantissa and exponent values to ac_fixed values, with rounding and saturation
  // as needed. It's specifically used for the adder tree functions which provide an ac_fixed accumulator
  // output.
  template<bool use_ac_sl, int W, int I, int E, int WR, int IR, bool SR, ac_q_mode QR, ac_o_mode OR>
  void acc_fxpt_conv(
    const ac_fixed<W, I, true> &mant, const ac_int<E, true> &exp,
    ac_fixed<WR, IR, SR, QR, OR> &acc
  ) {  
    if (use_ac_sl) {
      ac_fixed<WR + int(!SR), IR + int(!SR), true, QR, OR> fxpt_temp;
      ac_math::ac_shift_left(mant, exp.to_int(), fxpt_temp);
      acc = fxpt_temp;
    } else {
      typename ac_float<W, I, E>::rt_unary::to_ac_fixed_t fxpt_temp = mant;
      fxpt_temp <<= exp;
      acc = fxpt_temp;
    }
  }

  // This function assumes that IR is larger than I1 and I2 (no saturation when
  //   addition takes place)
  // Bits that are lost due to exponent alignment are truncated 
  template<
    int W1, int I1, ac_q_mode Q1,
    int W2, int I2, ac_q_mode Q2,
    int E,
    int WR, int IR, int ER, ac_q_mode QR
  >
  void sum2(
    const ac_fixed<W1,I1,true,Q1> &op1_m,
    const ac_int<E,true> &op1_e,
    const ac_fixed<W2,I2,true,Q2> &op2_m,
    const ac_int<E,true> &op2_e,
    ac_fixed<WR,IR,true,QR> &acc_mant,
    ac_int<ER,true> &acc_exp
  ) {
    typedef ac_fixed<WR,IR,true> m_t;
    int e_dif = op1_e - op2_e;
    bool e1_lt_e2 = e_dif < 0;
    bool op1_zero = !op1_m;
    bool op2_zero = !op2_m;
    e_dif = (op1_zero | op2_zero) ? 0 : e1_lt_e2 ? -e_dif : e_dif;
    m_t op_rshift = e1_lt_e2 ? (m_t) op1_m : (m_t) op2_m; 
    m_t op_no_shift = e1_lt_e2 ? (m_t) op2_m : (m_t) op1_m;
    // truncate bits (no rounding)
    op_rshift >>= e_dif;
    acc_mant = op_no_shift + op_rshift;
    acc_exp = !op2_zero && (op1_zero || e1_lt_e2) ? op2_e : op1_e;
  }

  template<int N_ELEMS>
  struct add_r {
    template<
      int WA, int IA, ac_q_mode QA,
      int WE, int IE, int EE, ac_q_mode QE
    >
    static void sum(
      const ac_float<WE,IE,EE,QE> *x,
      ac_fixed<WA,IA,true,QA> &acc_mant,
      ac_int<EE,true> &acc_exp
    ) {
      enum {
        N_ELEMS_L = N_ELEMS/2, N_ELEMS_H = N_ELEMS-N_ELEMS_L, 
        Growth_L = ac::log2_ceil<N_ELEMS_L>::val,
        Growth_H = ac::log2_ceil<N_ELEMS_H>::val,
        AccExtraLSBs = WA-IA-WE+IE,
        WA_L = WE+Growth_L+AccExtraLSBs,
        IA_L = IE+Growth_L,
        WA_H = WE+Growth_H+AccExtraLSBs,
        IA_H = IE+Growth_H
      };
      ac_fixed<WA_L,IA_L,true,QA> acc_mant_l;
      ac_int<EE,true> acc_exp_l;
      add_r<N_ELEMS_L>::sum(x, acc_mant_l, acc_exp_l);
    
      ac_fixed<WA_H,IA_H,true,QA> acc_mant_h;
      ac_int<EE,true> acc_exp_h;
      add_r<N_ELEMS_H>::sum(x+N_ELEMS_L, acc_mant_h, acc_exp_h);
    
      sum2(acc_mant_l, acc_exp_l, acc_mant_h, acc_exp_h, acc_mant, acc_exp);
    }
  };
  
  template<> struct add_r<1> {
    template<
      int WA, int IA, ac_q_mode QA,
      int WE, int IE, int EE, ac_q_mode QE
    >
    static void sum(
      const ac_float<WE,IE,EE,QE> *x,
      ac_fixed<WA,IA,true,QA> &acc_mant,
      ac_int<EE,true> &acc_exp
    ) {
      acc_mant = x->mantissa();
      acc_exp = x->exp();
    }
  };
  
  template<
    int N_ELEMS,
    int WA, int IA, ac_q_mode QA,
    int WE, int IE, int EE, ac_q_mode QE
  >
  void add_tree(
    const ac_float<WE,IE,EE,QE> *x,
    ac_fixed<WA,IA,true,QA> &acc_mant,
    ac_int<EE,true> &acc_exp
  ) {
    add_r<N_ELEMS>::sum(&x[0], acc_mant, acc_exp); 
  }
  

  //=========================================================================
  // Function: add_tree (ac_float inputs, ac_fixed output)
  //
  // Description:
  //    Adder tree implementation which uses recursive templates to model
  //    all the adder tree stages. The mantissa alignment and addition is
  //    done pairwise, and zero handling is built-in. The input array is
  //    passed as a pointer. The adder tree output is converted to an
  //    ac_fixed output. The user can adjust the saturation/overflow
  //    handling for this conversion by changing the saturation/overflow
  //    mode of the output type.
  //
  // Usage:
  //    See above example code for usage.
  //
  // Notes:
  //    The output conversion to ac_fixed may or may not use the
  //    ac_shift_left() function from ac_shift.h, and this choice
  //    may decrease or increase area based on the input exponent width.
  //    Refer to the AC Math reference manual for more details.
  //
  //-------------------------------------------------------------------------
  
  template<
    int AccExtraLSBs = 0,
    bool use_ac_sl = true,
    int WA, int IA, bool SA, ac_q_mode QA, ac_o_mode OA,
    int WE, int IE, int EE, ac_q_mode QE,
    int N_ELEMS
  >
  void add_tree(
    const ac_float<WE,IE,EE,QE> (&x)[N_ELEMS],
    ac_fixed<WA,IA,SA,QA,OA> &acc
  ) {
    enum {
      Growth = ac::log2_ceil<N_ELEMS>::val,
      WFX = WE+Growth+AccExtraLSBs,
      IFX = IE+Growth
    };
    ac_fixed<WFX,IFX,true> mant;
    ac_int<EE> exp;
    add_r<N_ELEMS>::sum(&x[0], mant, exp);
    acc_fxpt_conv<use_ac_sl>(mant, exp, acc);
  }

  //=========================================================================
  // Function: add_tree (ac_float inputs, ac_float output)
  //
  // Description:
  //    Adder tree implementation which uses recursive templates to model
  //    all the adder tree stages. The mantissa alignment and addition is
  //    done pairwise, and zero handling is built-in. The input array is
  //    passed as a reference.
  //
  // Usage:
  //    A sample testbench and its implementation look like
  //    this:
  //
  //    #include <ac_float_add_tree.h>
  //    using namespace ac_math;
  //    
  //    typedef ac_float<25, 2, 8, AC_TRN> input_type;
  //    typedef ac_float<25, 2, 9, AC_TRN> output_type;
  //    constexpr int N_ELEMS = 4;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const input_type (&input)[N_ELEMS],
  //      output_type &output
  //    )
  //    {
  //      add_tree(input, output);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //
  //    CCS_MAIN(int, char **)
  //    {
  //      input_type input[N_ELEMS];
  //      input[0] =  1.5;
  //      input[1] =  2.5;
  //      input[2] =  0.5;
  //      input[3] = -1.0;
  //      output_type output;
  //      
  //      CCS_DESIGN(project)(input, output);
  //      
  //      CCS_RETURN (0);
  //    }
  //    #endif
  //
  // Notes:
  //    The output conversion to ac_float has normalization turned on.
  //
  //-------------------------------------------------------------------------
  
  template<
    int AccExtraLSBs = 0,
    int WA, int IA, int EA, ac_q_mode QA,
    int WE, int IE, int EE, ac_q_mode QE,
    int N_ELEMS
  >
  void add_tree(
    const ac_float<WE,IE,EE,QE> (&x)[N_ELEMS],
    ac_float<WA,IA,EA,QA> &acc
  ) {
    enum {
      Growth = ac::log2_ceil<N_ELEMS>::val,
      WFX = WE+Growth+AccExtraLSBs,
      IFX = IE+Growth
    };
    ac_fixed<WFX,IFX,true> mant;
    ac_int<EE> exp;
    add_r<N_ELEMS>::sum(&x[0], mant, exp);
    acc = ac_float<WA,IA,EA,QA>(mant,exp,true); 
  }
  
  //=========================================================================
  // Function: add_tree_ptr (ac_float inputs, ac_fixed output)
  //
  // Description:
  //    A variant of the equivalent add_tree function defined above, with
  //    the only major difference being that the input array is passed as
  //    a pointer rather than a reference. Since size information for the
  //    pass-by-pointer version is lost, the user must explicitly specify
  //    N_ELEMS in the function call.
  //
  // Usage:
  //    A sample testbench and its implementation look like
  //    this:
  //
  //    #include <ac_float_add_tree.h>
  //    using namespace ac_math;
  //    
  //    typedef ac_float<25, 2, 8, AC_TRN> input_type;
  //    typedef ac_fixed<64, 32, true, AC_TRN, AC_SAT> output_type;
  //    constexpr int N_ELEMS = 4;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const input_type input[N_ELEMS],
  //      output_type &output
  //    )
  //    {
  //      add_tree_ptr<N_ELEMS>(input, output);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //
  //    CCS_MAIN(int, char **)
  //    {
  //      input_type input[N_ELEMS];
  //      input[0] =  1.5;
  //      input[1] =  2.5;
  //      input[2] =  0.5;
  //      input[3] = -1.0;
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
    int AccExtraLSBs = 0,
    bool use_ac_sl = true,
    int WA, int IA, bool SA, ac_q_mode QA, ac_o_mode OA,
    int WE, int IE, int EE, ac_q_mode QE
  >
  void add_tree_ptr(
    const ac_float<WE,IE,EE,QE> x[N_ELEMS],
    ac_fixed<WA,IA,SA,QA,OA> &acc
  ) {
    ac_float<WE,IE,EE,QE> x2[N_ELEMS];
    
    _Pragma ("hls_unroll yes")
    CPY_IN_ARR_AT_FXPT: for (int k = 0; k < N_ELEMS; k++) {
      x2[k] = x[k];
    }
    
    add_tree<AccExtraLSBs, use_ac_sl>(x2, acc);
  }
  
  //=========================================================================
  // Function: add_tree_ptr (ac_float inputs, ac_float output)
  //
  // Description:
  //    A variant of the equivalent add_tree function defined above, with
  //    the only major difference being that the input array is passed as
  //    a pointer rather than a reference. Since size information for the
  //    pass-by-pointer version is lost, the user must explicitly specify
  //    N_ELEMS in the function call.
  //
  // Usage:
  //    A sample testbench and its implementation look like
  //    this:
  //
  //    #include <ac_float_add_tree.h>
  //    using namespace ac_math;
  //    
  //    typedef ac_float<25, 2, 8, AC_TRN> input_type;
  //    typedef ac_float<25, 2, 9, AC_TRN> output_type;
  //    constexpr int N_ELEMS = 4;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const input_type input[N_ELEMS],
  //      output_type &output
  //    )
  //    {
  //      add_tree_ptr<N_ELEMS>(input, output);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //
  //    CCS_MAIN(int, char **)
  //    {
  //      input_type input[N_ELEMS];
  //      input[0] =  1.5;
  //      input[1] =  2.5;
  //      input[2] =  0.5;
  //      input[3] = -1.0;
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
    int AccExtraLSBs = 0,
    int WA, int IA, int EA, ac_q_mode QA,
    int WE, int IE, int EE, ac_q_mode QE
  >
  void add_tree_ptr(
    const ac_float<WE,IE,EE,QE> x[N_ELEMS],
    ac_float<WA,IA,EA,QA> &acc
  ) {
    ac_float<WE,IE,EE,QE> x2[N_ELEMS];
    
    _Pragma ("hls_unroll yes")
    CPY_IN_ARR_AT: for (int k = 0; k < N_ELEMS; k++) {
      x2[k] = x[k];
    }
    
    add_tree<AccExtraLSBs>(x2, acc);
  }
  
  // recursive template to compute max from array of ac_int
  template<int N_ELEMS>
  struct max_r {
    template<int W, bool S>
    static ac_int<W,S> max(const ac_int<W,S> *x) {
      enum { N_ELEMS_L = N_ELEMS/2, N_ELEMS_H = N_ELEMS-N_ELEMS_L};
      ac_int<W,S> max_l = max_r<N_ELEMS_L>::max(x); 
      ac_int<W,S> max_h = max_r<N_ELEMS_H>::max(x+N_ELEMS_L); 
      return (max_h > max_l) ? max_h : max_l;
    }
  };

  template<> struct max_r<1> {
    template<int W, bool S>
    static ac_int<W,S> max(const ac_int<W,S> *x) {
      return *x;
    }
  };

  template<
    int AccExtraLSBs = 0,
    bool ZeroHandling = true,
    int WE, int IE, int EE, ac_q_mode QE1,
    ac_q_mode QE2, ac_o_mode OE,
    int N_ELEMS
  >
  void exp_equalizer(
    const ac_float<WE,IE,EE,QE1> (&x)[N_ELEMS],
    ac_fixed<WE+AccExtraLSBs,IE,true,QE2,OE> (&y_mant)[N_ELEMS],
    ac_int<EE,true> &y_exp
  ) {
    ac_int<EE,true> y_exp_ary[N_ELEMS];
    
    _Pragma ("hls_unroll yes")
    CPY_EXP_VALS: for(int k=0; k < N_ELEMS; k++) {
      y_exp_ary[k] = zero_handle_exp<ZeroHandling>(x[k].exp(), x[k].mantissa());
    }
    
    y_exp = max_r<N_ELEMS>::max(&y_exp_ary[0]);
    
    _Pragma ("hls_unroll yes")
    EQUALIZE_EXPS: for(int k=0; k < N_ELEMS; k++) {
      ac_fixed<WE+AccExtraLSBs,IE,true> elem_mant = x[k].mantissa();
      ac_int<EE,false> exp_dif = y_exp - x[k].exp();
      elem_mant >>= exp_dif;
      y_mant[k] = elem_mant;
    }
  }
  
  template<
    int AccExtraLSBs = 0,
    bool ZeroHandling = true,
    int WA, int IA,
    int WE, int IE, int EE, ac_q_mode QE,
    int N_ELEMS
  >
  void block_add_tree_core(
    const ac_float<WE,IE,EE,QE> (&x)[N_ELEMS],
    ac_fixed<WA,IA,true> &acc_fx,
    ac_int<EE,true> &max_exp
  ) {
    enum {
      Growth = ac::log2_ceil<N_ELEMS>::val,
      WFX = WE+Growth+AccExtraLSBs,
      IFX = IE+Growth
    };
    
    static_assert(WA == WFX, "Output bitwidth must be equal to WE+log2_ceil(N_ELEMS)+AccExtraLSBs.");
    static_assert(IA == IFX, "Input bitwidth must be equal to IE+log2_ceil(N_ELEMS).");
    
    ac_fixed<WE+AccExtraLSBs,IE,true,AC_TRN,AC_WRAP> mants_algn[N_ELEMS];
    
    exp_equalizer<AccExtraLSBs,ZeroHandling,WE,IE,EE,QE,AC_TRN,AC_WRAP,N_ELEMS>(x, mants_algn, max_exp);
    
    ac_fixed<WFX,IFX,true> acc_fx_temp = 0.0;

    _Pragma("cluster addtree")
    _Pragma("cluster_type combinational")
    {
      _Pragma ("hls_unroll yes")
      BAT_FX_SUM: for (int k=0; k<N_ELEMS; k++) acc_fx_temp += mants_algn[k];
    }
    
    acc_fx = acc_fx_temp;
  }
  
  //=========================================================================
  // Function: block_add_tree (ac_float inputs, ac_fixed output)
  //
  // Description:
  //    "Block" implementation of adder tree, where the input mantissas are 
  //    aligned all at once with respect to the maximum exponent value and
  //    then added via an unrolled loop. The actual adder tree only operates
  //    with fixed point values, thereby saving on MUXes in its datapath
  //    and reducing the critical path. The adder tree output is combined
  //    with the maximum exponent value and converted to an ac_fixed output.
  //    The user can adjust the saturation/overflow handling for this
  //    conversion by changing the saturation/overflow mode of the output
  //    type.
  //
  // Usage:
  //    A sample testbench and its implementation look like
  //    this:
  //
  //    #include <ac_float_add_tree.h>
  //    using namespace ac_math;
  //    
  //    typedef ac_float<25, 2, 8, AC_TRN> input_type;
  //    typedef ac_fixed<64, 32, true, AC_TRN, AC_SAT> output_type;
  //    constexpr int N_ELEMS = 4;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const input_type (&input)[N_ELEMS],
  //      output_type &output
  //    )
  //    {
  //      block_add_tree(input, output);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //
  //    CCS_MAIN(int, char **)
  //    {
  //      input_type input[N_ELEMS];
  //      input[0] =  1.5;
  //      input[1] =  2.5;
  //      input[2] =  0.5;
  //      input[3] = -1.0;
  //      output_type output;
  //      
  //      CCS_DESIGN(project)(input, output);
  //      
  //      CCS_RETURN (0);
  //    }
  //    #endif
  //
  // Notes:
  //    The output conversion to ac_fixed may or may not use the
  //    ac_shift_left() function from ac_shift.h, and this choice
  //    may decrease or increase area based on the input exponent width.
  //    Refer to the AC Math reference manual for more details.
  //
  //-------------------------------------------------------------------------
  
  template<
    int AccExtraLSBs = 0,
    bool ZeroHandling = true,
    bool use_ac_sl = true,
    int WA, int IA, bool SA, ac_q_mode QA, ac_o_mode OA,
    int WE, int IE, int EE, ac_q_mode QE,
    int N_ELEMS
  >
  void block_add_tree(
    const ac_float<WE,IE,EE,QE> (&x)[N_ELEMS],
    ac_fixed<WA,IA,SA,QA,OA> &acc
  ) {
    enum {
      Growth = ac::log2_ceil<N_ELEMS>::val,
      WFX = WE+Growth+AccExtraLSBs,
      IFX = IE+Growth
    };
    
    ac_fixed<WFX,IFX,true> acc_fx;
    ac_int<EE,true> max_exp;
    block_add_tree_core<AccExtraLSBs, ZeroHandling>(x, acc_fx, max_exp);
    
    acc_fxpt_conv<use_ac_sl>(acc_fx, max_exp, acc);
  }
  
  //=========================================================================
  // Function: block_add_tree (ac_float inputs, ac_float output)
  //
  // Description:
  //    "Block" implementation of adder tree, where the input mantissas are 
  //    aligned all at once with respect to the maximum exponent value and
  //    then added via an unrolled loop. The actual adder tree only operates
  //    with fixed point values, thereby saving on MUXes in its datapath
  //    and reducing the critical path. The adder tree output is combined
  //    with the maximum exponent value and converted to an ac_float output.
  //
  // Usage:
  //    A sample testbench and its implementation look like
  //    this:
  //
  //    #include <ac_float_add_tree.h>
  //    using namespace ac_math;
  //    
  //    typedef ac_float<25, 2, 8, AC_TRN> input_type;
  //    typedef ac_float<25, 2, 9, AC_TRN> output_type;
  //    constexpr int N_ELEMS = 4;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const input_type (&input)[N_ELEMS],
  //      output_type &output
  //    )
  //    {
  //      block_add_tree(input, output);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //
  //    CCS_MAIN(int, char **)
  //    {
  //      input_type input[N_ELEMS];
  //      input[0] =  1.5;
  //      input[1] =  2.5;
  //      input[2] =  0.5;
  //      input[3] = -1.0;
  //      output_type output;
  //      
  //      CCS_DESIGN(project)(input, output);
  //      
  //      CCS_RETURN (0);
  //    }
  //    #endif
  //
  // Notes:
  //    The output conversion to ac_float has normalization turned on if the
  //    norm template parameter is set to true.
  //
  //-------------------------------------------------------------------------
  
  template<
    int AccExtraLSBs = 0,
    bool ZeroHandling = true,
    bool norm = true,
    int WA, int IA, int EA, ac_q_mode QA,
    int WE, int IE, int EE, ac_q_mode QE,
    int N_ELEMS
  >
  void block_add_tree(
    const ac_float<WE,IE,EE,QE> (&x)[N_ELEMS],
    ac_float<WA,IA,EA,QA> &acc
  ) {
    enum {
      Growth = ac::log2_ceil<N_ELEMS>::val,
      WFX = WE+Growth+AccExtraLSBs,
      IFX = IE+Growth
    };
    
    ac_fixed<WFX,IFX,true> acc_fx;
    ac_int<EE,true> max_exp;
    
    block_add_tree_core<AccExtraLSBs, ZeroHandling>(x, acc_fx, max_exp);
    
    acc = ac_float<WA,IA,EA,QA>(acc_fx,max_exp,norm);
  }
  
  //=========================================================================
  // Function: block_add_tree_ptr (ac_float inputs, ac_fixed output)
  //
  // Description:
  //    A variant of the equivalent block_add_tree function defined above,
  //    with the only major difference being that the input array is passed
  //    as a pointer rather than a reference. Since size information for the
  //    pass-by-pointer version is lost, the user must explicitly specify
  //    N_ELEMS in the function call.
  //
  // Usage:
  //    A sample testbench and its implementation look like
  //    this:
  //
  //    #include <ac_float_add_tree.h>
  //    using namespace ac_math;
  //    
  //    typedef ac_float<25, 2, 8, AC_TRN> input_type;
  //    typedef ac_fixed<64, 32, true, AC_TRN, AC_SAT> output_type;
  //    constexpr int N_ELEMS = 4;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const input_type input[N_ELEMS],
  //      output_type &output
  //    )
  //    {
  //      block_add_tree_ptr<N_ELEMS>(input, output);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //
  //    CCS_MAIN(int, char **)
  //    {
  //      input_type input[N_ELEMS];
  //      input[0] =  1.5;
  //      input[1] =  2.5;
  //      input[2] =  0.5;
  //      input[3] = -1.0;
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
    int AccExtraLSBs = 0,
    bool ZeroHandling = true,
    bool use_ac_sl = true,
    int WA, int IA, bool SA, ac_q_mode QA, ac_o_mode OA,
    int WE, int IE, int EE, ac_q_mode QE
  >
  void block_add_tree_ptr(
    const ac_float<WE,IE,EE,QE> x[N_ELEMS],
    ac_fixed<WA,IA,SA,QA,OA> &acc
  ) {
    ac_float<WE,IE,EE,QE> x2[N_ELEMS];
    
    _Pragma ("hls_unroll yes")
    CPY_IN_ARR_BAT_FXPT: for (int k = 0; k < N_ELEMS; k++) {
      x2[k] = x[k];
    }
    
    block_add_tree<AccExtraLSBs, ZeroHandling, use_ac_sl>(x2, acc);
  }
  
  //=========================================================================
  // Function: block_add_tree_ptr (ac_float inputs, ac_float output)
  //
  // Description:
  //    A variant of the equivalent block_add_tree function defined above,
  //    with the only major difference being that the input array is passed
  //    as a pointer rather than a reference. Since size information for the
  //    pass-by-pointer version is lost, the user must explicitly specify
  //    N_ELEMS in the function call.
  //
  // Usage:
  //    A sample testbench and its implementation look like
  //    this:
  //
  //    #include <ac_float_add_tree.h>
  //    using namespace ac_math;
  //    
  //    typedef ac_float<25, 2, 8, AC_TRN> input_type;
  //    typedef ac_float<25, 2, 9, AC_TRN> output_type;
  //    constexpr int N_ELEMS = 4;
  //    
  //    #pragma hls_design top
  //    void project(
  //      const input_type input[N_ELEMS],
  //      output_type &output
  //    )
  //    {
  //      block_add_tree_ptr<N_ELEMS>(input, output);
  //    }
  //    
  //    #ifndef __SYNTHESIS__
  //    #include <mc_scverify.h>
  //
  //    CCS_MAIN(int, char **)
  //    {
  //      input_type input[N_ELEMS];
  //      input[0] =  1.5;
  //      input[1] =  2.5;
  //      input[2] =  0.5;
  //      input[3] = -1.0;
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
    int AccExtraLSBs = 0,
    bool ZeroHandling = true,
    bool norm = true,
    int WA, int IA, int EA, ac_q_mode QA,
    int WE, int IE, int EE, ac_q_mode QE
  >
  void block_add_tree_ptr(
    const ac_float<WE,IE,EE,QE> x[N_ELEMS],
    ac_float<WA,IA,EA,QA> &acc
  ) {
    ac_float<WE,IE,EE,QE> x2[N_ELEMS];
    
    _Pragma ("hls_unroll yes")
    CPY_IN_ARR_BAT: for (int k = 0; k < N_ELEMS; k++) {
      x2[k] = x[k];
    }
    
    block_add_tree<AccExtraLSBs, ZeroHandling, norm>(x2, acc);
  }
}

#endif
