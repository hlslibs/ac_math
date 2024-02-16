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
  
  template<
    int N_ELEMS,
    int AccExtraLSBs = 0,
    bool use_ac_sl = true,
    int WA, int IA, bool SA, ac_q_mode QA, ac_o_mode OA,
    int WE, int IE, int EE, ac_q_mode QE
  >
  void add_tree_ptr(
    const ac_float<WE,IE,EE,QE> *x,
    ac_fixed<WA,IA,SA,QA,OA> &acc
  ) {
    ac_float<WE,IE,EE,QE> x2[N_ELEMS];
    
    _Pragma ("hls_unroll yes")
    CPY_IN_ARR_AT_FXPT: for (int k = 0; k < N_ELEMS; k++) {
      x2[k] = x[k];
    }
    
    add_tree<AccExtraLSBs, use_ac_sl>(x2, acc);
  }
  
  template<
    int N_ELEMS,
    int AccExtraLSBs = 0,
    int WA, int IA, int EA, ac_q_mode QA,
    int WE, int IE, int EE, ac_q_mode QE
  >
  void add_tree_ptr(
    const ac_float<WE,IE,EE,QE> *x,
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
  
  template<
    int N_ELEMS,
    int AccExtraLSBs = 0,
    bool ZeroHandling = true,
    bool use_ac_sl = true,
    int WA, int IA, bool SA, ac_q_mode QA, ac_o_mode OA,
    int WE, int IE, int EE, ac_q_mode QE
  >
  void block_add_tree_ptr(
    const ac_float<WE,IE,EE,QE> *x,
    ac_fixed<WA,IA,SA,QA,OA> &acc
  ) {
    ac_float<WE,IE,EE,QE> x2[N_ELEMS];
    
    _Pragma ("hls_unroll yes")
    CPY_IN_ARR_BAT_FXPT: for (int k = 0; k < N_ELEMS; k++) {
      x2[k] = x[k];
    }
    
    block_add_tree<AccExtraLSBs, ZeroHandling, use_ac_sl>(x2, acc);
  }
  
  template<
    int N_ELEMS,
    int AccExtraLSBs = 0,
    bool ZeroHandling = true,
    bool norm = true,
    int WA, int IA, int EA, ac_q_mode QA,
    int WE, int IE, int EE, ac_q_mode QE
  >
  void block_add_tree_ptr(
    const ac_float<WE,IE,EE,QE> *x,
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
