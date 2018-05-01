/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 2.0                                                 *
 *                                                                        *
 *  Release Date    : Tue May  1 13:47:52 PDT 2018                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 2.0.2                                               *
 *                                                                        *
 *  Copyright , Mentor Graphics Corporation,                     *
 *                                                                        *
 *  All Rights Reserved.                                                  *
 *  
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
// File: ac_hcordic.h
//
// Description: Hyperbolic CORDIC implementations of synthesizable
//              log/exp/log2/exp2/pow functions for AC fixed point
//              datatypes
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//
//    typedef ac_fixed<16,  8, false, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<32, 16,  true, AC_RND, AC_SAT> output_type;
//
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output
//    )
//    {
//      ac_log2_cordic(input, output);
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 2.5;
//      output_type output;
//      CCS_DESIGN(project)(input, output);
//      CCS_RETURN (0);
//    }
//    #endif
//
// Notes:
//    This file uses C++ function overloading to target implementations
//    specific to each type of data. Attempting to call the function
//    with a type that is not implemented will result in a compile error.
//
//*****************************************************************************************

#ifndef _INCLUDED_AC_HCORDIC_H_
#define _INCLUDED_AC_HCORDIC_H_

#ifndef __cplusplus
#error C++ is required to include this header file
#endif

#if defined(AC_HCORDIC_H_DEBUG) && !defined(__SYNTHESIS__)
#include <iostream>
using namespace std;
#endif

#include <ac_int.h>
#include <ac_fixed.h>

namespace ac_math
{

  template <bool b> struct AcHtrigAssert { /* Compile error, no 'test' symbol. */ };
  template <> struct AcHtrigAssert<true> { enum { test }; };

  //  Multi-precision approximation of tanh(2^-j).
  //  53 + 3 + 1 -- double-precision + extra-iters + padding
  // + 6         -- ceil(log2(57))
  //  63
  static const ac_fixed<63,0,false> hcordic_table[] = {
    ac_fixed<63,0,false>(0.54930614433405478) + ac_fixed<63,0,false>(6.5665816287507998e-17),
    ac_fixed<63,0,false>(0.25541281188299530) + ac_fixed<63,0,false>(3.6245743399179330e-17),
    ac_fixed<63,0,false>(0.12565721414045302) + ac_fixed<63,0,false>(1.1446974688604729e-17),
    ac_fixed<63,0,false>(0.062581571477002995) + ac_fixed<63,0,false>(1.1963745773904336e-17),
    ac_fixed<63,0,false>(0.031260178490666992) + ac_fixed<63,0,false>(2.0406997556843409e-18),
    ac_fixed<63,0,false>(0.015626271752052209) + ac_fixed<63,0,false>(2.0718863563586340e-18),
    ac_fixed<63,0,false>(0.0078126589515404194) + ac_fixed<63,0,false>(1.4328684251271785e-18),
    ac_fixed<63,0,false>(0.0039062698683968253) + ac_fixed<63,0,false>(7.1043409533913290e-19),
    ac_fixed<63,0,false>(0.0019531274835325497) + ac_fixed<63,0,false>(2.4678511188227481e-19),
    ac_fixed<63,0,false>(0.00097656281044103572) + ac_fixed<63,0,false>(1.1576923652829265e-19),
    ac_fixed<63,0,false>(0.00048828128880511276) + ac_fixed<63,0,false>(5.7825061215340352e-20),
    ac_fixed<63,0,false>(0.00024414062985063855) + ac_fixed<63,0,false>(2.8912065318488310e-20),
    ac_fixed<63,0,false>(0.00012207031310632979) + ac_fixed<63,0,false>(1.4456029024172932e-20),
    ac_fixed<63,0,false>(0.000061035156325791220) + ac_fixed<63,0,false>(4.6869156419245725e-21),
    ac_fixed<63,0,false>(0.000030517578134473900) + ac_fixed<63,0,false>(2.2640484819353284e-21),
    ac_fixed<63,0,false>(0.000015258789063684236) + ac_fixed<63,0,false>(1.1295426991282718e-21),
    ac_fixed<63,0,false>(7.6293945313980291e-6) + ac_fixed<63,0,false>(5.6469380138169552e-22),
    ac_fixed<63,0,false>(3.8146972656435034e-6) + ac_fixed<63,0,false>(2.8234447731014682e-22),
    ac_fixed<63,0,false>(1.9073486328148128e-6) + ac_fixed<63,0,false>(1.4117216292442651e-22),
    ac_fixed<63,0,false>(9.5367431640653904e-7) + ac_fixed<63,0,false>(7.0586079095630541e-23),
    ac_fixed<63,0,false>(4.7683715820316110e-7) + ac_fixed<63,0,false>(3.5293039473859560e-23),
    ac_fixed<63,0,false>(2.3841857910156699e-7) + ac_fixed<63,0,false>(1.7646519734618664e-23),
    ac_fixed<63,0,false>(1.1920928955078180e-7) + ac_fixed<63,0,false>(8.8232598672371098e-24),
    ac_fixed<63,0,false>(5.9604644775390691e-8) + ac_fixed<63,0,false>(4.4116299336162979e-24),
    ac_fixed<63,0,false>(2.9802322387695319e-8) + ac_fixed<63,0,false>(2.2058149668080784e-24),
    ac_fixed<63,0,false>(1.4901161193847656e-8) + ac_fixed<63,0,false>(1.1029074834040370e-24),
    ac_fixed<63,0,false>(7.4505805969238281e-9) + ac_fixed<63,0,false>(1.3786343542550460e-25),
    ac_fixed<63,0,false>(3.7252902984619140e-9) + ac_fixed<63,0,false>(1.7232929428188075e-26),
    ac_fixed<63,0,false>(1.8626451492309570e-9) + ac_fixed<63,0,false>(2.1541161785235094e-27),
    ac_fixed<63,0,false>(9.3132257461547851e-10) + ac_fixed<63,0,false>(2.6926452231543868e-28),
    ac_fixed<63,0,false>(4.6566128730773925e-10) + ac_fixed<63,0,false>(3.3658065289429835e-29),
    ac_fixed<63,0,false>(2.3283064365386962e-10) + ac_fixed<63,0,false>(4.2072581611787293e-30),
    ac_fixed<63,0,false>(1.1641532182693481e-10) + ac_fixed<63,0,false>(5.2590727014734117e-31),
    ac_fixed<63,0,false>(5.8207660913467407e-11) + ac_fixed<63,0,false>(6.5738408768417646e-32),
    ac_fixed<63,0,false>(2.9103830456733703e-11) + ac_fixed<63,0,false>(8.2173010960522058e-33),
    ac_fixed<63,0,false>(1.4551915228366851e-11) + ac_fixed<63,0,false>(1.0271626370065257e-33),
    ac_fixed<63,0,false>(7.2759576141834259e-12) + ac_fixed<63,0,false>(1.2839532962581571e-34),
    ac_fixed<63,0,false>(3.6379788070917129e-12) + ac_fixed<63,0,false>(1.6049416203226964e-35),
    ac_fixed<63,0,false>(1.8189894035458564e-12) + ac_fixed<63,0,false>(2.0061770254033705e-36),
    ac_fixed<63,0,false>(9.0949470177292823e-13) + ac_fixed<63,0,false>(2.5077212817542132e-37),
    ac_fixed<63,0,false>(4.5474735088646411e-13) + ac_fixed<63,0,false>(3.1346516021927665e-38),
    ac_fixed<63,0,false>(2.2737367544323205e-13) + ac_fixed<63,0,false>(3.9183145027409581e-39),
    ac_fixed<63,0,false>(1.1368683772161602e-13) + ac_fixed<63,0,false>(4.8978931284261976e-40),
    ac_fixed<63,0,false>(5.6843418860808014e-14) + ac_fixed<63,0,false>(6.1223664105327470e-41),
    ac_fixed<63,0,false>(2.8421709430404007e-14) + ac_fixed<63,0,false>(7.6529580131659338e-42),
    ac_fixed<63,0,false>(1.4210854715202003e-14) + ac_fixed<63,0,false>(9.5661975164574173e-43),
    ac_fixed<63,0,false>(7.1054273576010018e-15) + ac_fixed<63,0,false>(1.1957746895571771e-43),
    ac_fixed<63,0,false>(3.5527136788005009e-15) + ac_fixed<63,0,false>(1.4947183619464714e-44),
    ac_fixed<63,0,false>(1.7763568394002504e-15) + ac_fixed<63,0,false>(1.8683979524330893e-45),
    ac_fixed<63,0,false>(8.8817841970012523e-16) + ac_fixed<63,0,false>(2.3354974405413616e-46),
    ac_fixed<63,0,false>(4.4408920985006261e-16) + ac_fixed<63,0,false>(2.9193718006767020e-47),
    ac_fixed<63,0,false>(2.2204460492503130e-16) + ac_fixed<63,0,false>(3.6492147508458775e-48),
    ac_fixed<63,0,false>(1.1102230246251565e-16) + ac_fixed<63,0,false>(4.5615184385573469e-49),
    ac_fixed<63,0,false>(5.5511151231257827e-17) + ac_fixed<63,0,false>(5.7018980481966837e-50),
    ac_fixed<63,0,false>(2.7755575615628913e-17) + ac_fixed<63,0,false>(7.1273725602458546e-51),
    ac_fixed<63,0,false>(1.3877787807814456e-17) + ac_fixed<63,0,false>(8.9092157003073183e-52),
    ac_fixed<63,0,false>(6.9388939039072283e-18) + ac_fixed<63,0,false>(1.1136519625384147e-52),
  };

  // Multi-precision approximation of tanh(2^-j)/ln(2)
  static const ac_fixed<63,0,false> hcordic_table_inv_ln2[] = {
    ac_fixed<63,0,false>(0.79248125036057803) + ac_fixed<63,0,false>(5.2898906200562771e-17),
    ac_fixed<63,0,false>(0.36848279708310305) + ac_fixed<63,0,false>(3.0181969648116828e-17),
    ac_fixed<63,0,false>(0.18128503969235412) + ac_fixed<63,0,false>(3.2751601209121701e-19),
    ac_fixed<63,0,false>(0.090286122820910433) + ac_fixed<63,0,false>(6.1807566423978258e-18),
    ac_fixed<63,0,false>(0.045098904485789112) + ac_fixed<63,0,false>(1.9394826157873123e-18),
    ac_fixed<63,0,false>(0.022543944764269015) + ac_fixed<63,0,false>(3.0409189283261096e-18),
    ac_fixed<63,0,false>(0.011271284325544132) + ac_fixed<63,0,false>(7.4172706383649422e-19),
    ac_fixed<63,0,false>(0.0056355561675100838) + ac_fixed<63,0,false>(8.4327977335471815e-19),
    ac_fixed<63,0,false>(0.0028177673347163502) + ac_fixed<63,0,false>(2.0019188336311420e-19),
    ac_fixed<63,0,false>(0.0014088823237398710) + ac_fixed<63,0,false>(1.8534395143051122e-19),
    ac_fixed<63,0,false>(0.00070444099391800798) + ac_fixed<63,0,false>(1.9981201974068463e-20),
    ac_fixed<63,0,false>(0.00035222047596502426) + ac_fixed<63,0,false>(4.1468140580709961e-20),
    ac_fixed<63,0,false>(0.00017611023535826504) + ac_fixed<63,0,false>(3.9132915644949836e-21),
    ac_fixed<63,0,false>(0.000088055117351101637) + ac_fixed<63,0,false>(7.4640093593241846e-21),
    ac_fixed<63,0,false>(0.000044027558634546956) + ac_fixed<63,0,false>(6.4581814039897067e-21),
    ac_fixed<63,0,false>(0.000022013779312147997) + ac_fixed<63,0,false>(1.8865372223735620e-21),
    ac_fixed<63,0,false>(0.000011006889655433313) + ac_fixed<63,0,false>(1.1993015349962144e-21),
    ac_fixed<63,0,false>(5.5034448276365712e-6) + ac_fixed<63,0,false>(2.0814889794509639e-22),
    ac_fixed<63,0,false>(2.7517224138082749e-6) + ac_fixed<63,0,false>(5.5137043047108037e-23),
    ac_fixed<63,0,false>(1.3758612069028860e-6) + ac_fixed<63,0,false>(1.8027003363582673e-22),
    ac_fixed<63,0,false>(6.8793060345128661e-7) + ac_fixed<63,0,false>(5.6283146948640298e-23),
    ac_fixed<63,0,false>(3.4396530172562377e-7) + ac_fixed<63,0,false>(1.0675199949815282e-23),
    ac_fixed<63,0,false>(1.7198265086280942e-7) + ac_fixed<63,0,false>(2.3006637985929779e-23),
    ac_fixed<63,0,false>(8.5991325431404408e-8) + ac_fixed<63,0,false>(1.0403226294140314e-23),
    ac_fixed<63,0,false>(4.2995662715702171e-8) + ac_fixed<63,0,false>(1.0101788439922535e-25),
    ac_fixed<63,0,false>(2.1497831357851078e-8) + ac_fixed<63,0,false>(1.8944763720248383e-24),
    ac_fixed<63,0,false>(1.0748915678925539e-8) + ac_fixed<63,0,false>(3.5055350218754517e-25),
    ac_fixed<63,0,false>(5.3744578394627697e-9) + ac_fixed<63,0,false>(1.0069116561566334e-25),
    ac_fixed<63,0,false>(2.6872289197313848e-9) + ac_fixed<63,0,false>(4.1022384623068012e-26),
    ac_fixed<63,0,false>(1.3436144598656924e-9) + ac_fixed<63,0,false>(1.9345792538438548e-26),
    ac_fixed<63,0,false>(6.7180722993284621e-10) + ac_fixed<63,0,false>(9.5272212975823422e-27),
    ac_fixed<63,0,false>(3.3590361496642310e-10) + ac_fixed<63,0,false>(4.7454012773365546e-27),
    ac_fixed<63,0,false>(1.6795180748321155e-10) + ac_fixed<63,0,false>(2.3704244672364501e-27),
    ac_fixed<63,0,false>(8.3975903741605777e-11) + ac_fixed<63,0,false>(1.1849277121892467e-27),
    ac_fixed<63,0,false>(4.1987951870802888e-11) + ac_fixed<63,0,false>(5.9242829091600106e-28),
    ac_fixed<63,0,false>(2.0993975935401444e-11) + ac_fixed<63,0,false>(2.9620969981067276e-28),
    ac_fixed<63,0,false>(1.0496987967700722e-11) + ac_fixed<63,0,false>(1.4810429419942040e-28),
    ac_fixed<63,0,false>(5.2484939838503610e-12) + ac_fixed<63,0,false>(7.4052077636470704e-29),
    ac_fixed<63,0,false>(2.6242469919251805e-12) + ac_fixed<63,0,false>(3.7026030135330413e-29),
    ac_fixed<63,0,false>(1.3121234959625902e-12) + ac_fixed<63,0,false>(1.8513013982302090e-29),
    ac_fixed<63,0,false>(6.5606174798129513e-13) + ac_fixed<63,0,false>(9.2565068554806557e-30),
    ac_fixed<63,0,false>(3.2803087399064756e-13) + ac_fixed<63,0,false>(4.6282534107815294e-30),
    ac_fixed<63,0,false>(1.6401543699532378e-13) + ac_fixed<63,0,false>(2.3141267032709147e-30),
    ac_fixed<63,0,false>(8.2007718497661891e-14) + ac_fixed<63,0,false>(1.1570633513704762e-30),
    ac_fixed<63,0,false>(4.1003859248830945e-14) + ac_fixed<63,0,false>(5.7853167565211543e-31),
    ac_fixed<63,0,false>(2.0501929624415472e-14) + ac_fixed<63,0,false>(2.8926583782191736e-31),
    ac_fixed<63,0,false>(1.0250964812207736e-14) + ac_fixed<63,0,false>(1.4463291891044114e-31),
    ac_fixed<63,0,false>(5.1254824061038682e-15) + ac_fixed<63,0,false>(7.2316459455155881e-32),
    ac_fixed<63,0,false>(2.5627412030519341e-15) + ac_fixed<63,0,false>(3.6158229727569855e-32),
    ac_fixed<63,0,false>(1.2813706015259670e-15) + ac_fixed<63,0,false>(1.8079114863783915e-32),
    ac_fixed<63,0,false>(6.4068530076298352e-16) + ac_fixed<63,0,false>(9.0395574318918317e-33),
    ac_fixed<63,0,false>(3.2034265038149176e-16) + ac_fixed<63,0,false>(4.5197787159459001e-33),
    ac_fixed<63,0,false>(1.6017132519074588e-16) + ac_fixed<63,0,false>(2.2598893579729480e-33),
    ac_fixed<63,0,false>(8.0085662595372941e-17) + ac_fixed<63,0,false>(1.1299446789864738e-33),
    ac_fixed<63,0,false>(4.0042831297686470e-17) + ac_fixed<63,0,false>(5.6497233949323683e-34),
    ac_fixed<63,0,false>(2.0021415648843235e-17) + ac_fixed<63,0,false>(2.8248616974661841e-34),
    ac_fixed<63,0,false>(1.0010707824421617e-17) + ac_fixed<63,0,false>(1.4124308487330920e-34),
  };

  static const int shift_dist_table[60] = {
    1, 2, 3, 4, 4, 5, 6, 7, 8, 9,
    10,11,12,13,13,14,15,16,17,18,
    19,20,21,22,23,24,25,26,27,28,
    29,30,31,32,33,34,35,36,37,38,
    39,40,40,41,42,43,44,45,46,47,
    48,49,50,51,52,53,54,55,56,57
  };

  // Statically computed log2(x)
  template <int W>
  struct slog2 {
    enum {
      floor = 1 + slog2< (W >> 1) >::floor,
      pow2 = (1 << floor) == W,
      ceil = pow2 ? floor : floor + 1
    };
  };

  template <>
  struct slog2<1> {
    enum {
      floor = 0,
      pow2 = 1,
      ceil = 0
    };
  };

  template<int J>
  struct XtraIters {
    static const int valid = AcHtrigAssert<(J <= 60)>::test;
    enum {
      // sum-of-cordic-terms for {i | i > j + 4*B2+2*B1+B0} < ulp(iteration(j))
      // j >=14 && j <= 60: B0 + B1 + B2 = 3
      // j >= 4 && j <= 13: B0 + B1 + B2 = 2
      // j >= 1 && j <=  3: B0 + B1 + B2 = 1
      B0 = (1 << 0)*((J >= 1 && J <= 3) || (J >= 14)),
      B1 = (1 << 1)*(J >= 4 && J <= 13),
      B2 = (1 << 2)*(J >= 14)
    };
  };

  static const ac_fixed<159,0,false> _ln2 =
    ac_fixed<159,0,false>(0.69314718055994528) +
    ac_fixed<159,0,false>(2.3190468138462995e-17) +
    ac_fixed<159,0,false>(5.7077084384162112e-34);

  template<int W>
  ac_fixed<W,0,false> ln2_function()
  {
    (void)AcHtrigAssert< (W <= 159) >::test;
    return _ln2;
  }

  static const ac_fixed<159,1,false> _inv_ln2 =
    ac_fixed<159,1,false>(1.4426950408889633) +
    ac_fixed<159,1,false>(2.0355273740931030e-17) +
    ac_fixed<159,1,false>(2.0200219154078514e-33);

  template<int W>
  ac_fixed<W,1,false> inv_ln2_function()
  {
    (void)AcHtrigAssert< (W <= 158) >::test;
    return _inv_ln2;
  }

  static const ac_fixed<159,1,false> _inv_K =
    ac_fixed<159,1,false>(1.2074970677630720) +
    ac_fixed<159,1,false>(4.3877290122698160e-17) +
    ac_fixed<159,1,false>(7.4569481788958734e-34);

  template<int W>
  ac_fixed<W+1,1,false> inv_K()
  {
    (void)AcHtrigAssert< (W <= 158) >::test;
    return _inv_K;
  }

  // Find the positition of the msb (zero-indexed) in an unsigned integer x.
  // r     : msb position.
  // valid : equal to zero if all bits of x are zero, otherwise 1.
  template <int W, int WO>
  struct LOD {
    static void get(ac_int<W,false> x, ac_int<WO,false> &r, bool &valid)
    {
      const int N = slog2<W>::ceil - 1;
      ac_int<W, true> mask(-1);
      mask <<= (1 << N);
      bool sel = (x & mask);
      valid |= sel;
      r |= sel << N;
      LOD<(1 << N), WO>::get(sel ? (x >> (1 << N)) : x, r, valid);
    }
  };

  template <int WO>
  struct LOD<1, WO> {
    static void get(ac_int<1,false> x, ac_int<WO,false> &r, bool &valid)
    {
      valid |= (x & 1);
    }
  };

  template< int ZW >
  static ac_fixed< ZW, 0, false > hcordic_table_function(int i)
  {
    (void)AcHtrigAssert< ZW <= 56 >::test;
    return hcordic_table[i];
  }

  template< int ZW >
  static int shift_dist(int i)
  {
    (void)AcHtrigAssert< ZW <= 56 >::test;
    return shift_dist_table[i];
  }


  template< int ZW >
  static ac_fixed< ZW, 0, false > hcordic_table_inv_ln2_function(int i)
  {
    (void)AcHtrigAssert< ZW <= 56 >::test;
    return hcordic_table_inv_ln2[i];
  }

  struct HCordicConstants {
    enum mode {
      ROT_Y,
      ROT_Z,
    };
    enum scale {
      SCALE_1,
      SCALE_LN2,
    };
  };

  template< int W, int I, bool S, ac_q_mode Q, ac_o_mode V,
            int ZW, int ZI, bool ZS, ac_q_mode ZQ, ac_o_mode ZV,
            enum HCordicConstants::mode mode,
            enum HCordicConstants::scale scale >
  void ac_hcordic(ac_fixed<W,I,S,Q,V> &x,
                  ac_fixed<W,I,S,Q,V> &y,
                  ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
  {
    const int L = ZW - ZI + (XtraIters<ZW-ZI>::B0 +
                             XtraIters<ZW-ZI>::B1 +
                             XtraIters<ZW-ZI>::B2);
    const int LW = L+slog2<L>::ceil;
    typedef ac_fixed<LW,I+1,true> dp_t;
    dp_t xi = x;
    dp_t yi = y;
    dp_t zi = z;
    dp_t yi_d, xi_d;
    dp_t xi_n, yi_n, zi_n;
    for (int j = 0; j < L; j++) {
      int step = shift_dist<LW>(j);
      dp_t table;
      if (scale == HCordicConstants::SCALE_1) {
        table = hcordic_table_function<LW>(step - 1);
      }
      if (scale == HCordicConstants::SCALE_LN2) {
        table = hcordic_table_inv_ln2_function<LW>(step - 1);
      }
      bool dir = false; // true -> +, false -> -
      if (mode == HCordicConstants::ROT_Y) {
        dir = yi < 0;
      }
      if (mode == HCordicConstants::ROT_Z) {
        dir = zi >= 0;
      }
      xi_d = xi >> step;
      yi_d = yi >> step;
      if (dir) {
        xi += yi_d;
        yi += xi_d;
        zi -= table;
      } else {
        xi -= yi_d;
        yi -= xi_d;
        zi += table;
      }
    }
    x = xi;
    y = yi;
    z = zi;
  }

  // Range Reduced to: 0.5 <= x < 1, -.69 < z < 0
  template <int AW,int AI,bool AS,ac_q_mode AQ,ac_o_mode AV,
            int ZW,int ZI,bool ZS,ac_q_mode ZQ,ac_o_mode ZV>
  void ac_ln_rr(const ac_fixed<AW,AI,AS,AQ,AV> &x, ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
  {
    const int EW = (AW-AI) > (ZW-ZI) ? (AW-AI) : (ZW-ZI);
    ac_fixed<EW+2,2,true> xc = x;
    xc += 1;
    ac_fixed<EW+2,2,true> yc = x;
    yc -= 1;
    ac_fixed<EW+3,2,true> zc = 0; // Post-multiply by 2. Compute an extra bit.
    ac_hcordic<EW+2,2,true,AC_TRN,AC_WRAP,
               EW+3,2,true,AC_TRN,AC_WRAP,
               HCordicConstants::ROT_Y,
               HCordicConstants::SCALE_1>(xc,yc,zc);
    z = 2*zc;
  }

  // Range reduced to 0.5 <= x < 1, -1 <= z < 1
  template <int AW, int AI, bool AS, ac_q_mode AQ, ac_o_mode AV,
            int ZW, int ZI, bool ZS, ac_q_mode ZQ, ac_o_mode ZV>
  void ac_log2_rr(const ac_fixed<AW,AI,AS,AQ,AV> &x, ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
  {
    const int EW = (AW-AI) > (ZW-ZI) ? (AW-AI) : (ZW-ZI);
    ac_fixed<EW+2,2,true> xc = x;
    xc += 1;
    ac_fixed<EW+2,2,true> yc = x;
    yc -= 1;
    ac_fixed<EW+3,2,true> zc = 0; // Post-multiply by 2. Compute an extra bit.
    ac_hcordic<EW+2,2,true,AC_TRN,AC_WRAP,
               EW+3,2,true,AC_TRN,AC_WRAP,
               HCordicConstants::ROT_Y,
               HCordicConstants::SCALE_LN2>(xc,yc,zc);
    z = 2*zc;
  }

  // Range Reduced to: |x| < ln(2)  -0.69 < x < 0.69, -2 < z < 2
  template <int AW,int AI,bool AS,ac_q_mode AQ,ac_o_mode AV,
            int ZW,int ZI,bool ZS,ac_q_mode ZQ,ac_o_mode ZV>
  void ac_exp_rr(const ac_fixed<AW,AI,AS,AQ,AV> &x, ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
  {
    const int EW = ZW-ZI; // Evaluation width.
    ac_fixed<EW+3,3,true> xc = 1.0;
    ac_fixed<EW+3,3,true> yc = 1.0;
    ac_fixed<EW+3,3,true> zc = x;
    ac_hcordic<EW+3,3,true,AC_TRN,AC_WRAP,
               EW+3,3,true,AC_TRN,AC_WRAP,
               HCordicConstants::ROT_Z,
               HCordicConstants::SCALE_1>(xc,yc,zc);
    ac_fixed<EW+3,3,true> k = inv_K<EW>();
    yc *= k;
    z = yc;
  }

  // Range Reduced to: 0 <= |x| < 1  0 <= x < 1, 1 <= z < 2
  template <int AW,int AI,bool AS,ac_q_mode AQ,ac_o_mode AV,
            int ZW,int ZI,bool ZS,ac_q_mode ZQ,ac_o_mode ZV>
  void ac_exp2_rr(const ac_fixed<AW,AI,AS,AQ,AV> &x, ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
  {
    const int EW = AW > ZW ? AW : ZW;
    ac_fixed<EW+3,3,true> xc = 1.0;
    ac_fixed<EW+3,3,true> yc = 1.0;
    ac_fixed<EW+3,3,true> ln2 = ln2_function<EW>();
    ac_fixed<EW+3,3,true> zc = x*ln2;
    ac_hcordic<EW+3,3,true,AC_TRN,AC_WRAP,
               EW+3,3,true,AC_TRN,AC_WRAP,
               HCordicConstants::ROT_Z,
               HCordicConstants::SCALE_1>(xc,yc,zc);
    ac_fixed<EW+3,3,true> k = inv_K<EW>();
    yc *= k;
    z = yc;
  }

  struct AcLogRR {
    enum base {
      BASE_E,
      BASE_2,
    };
  };

  template <int AW, int AI, bool AS, ac_q_mode AQ, ac_o_mode AV,
            int ZW, int ZI, bool ZS, ac_q_mode ZQ, ac_o_mode ZV, enum AcLogRR::base BASE >
  void ac_log_(const ac_fixed<AW,AI,AS,AQ,AV> &x, ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
  {
    // BASE_E: RR: ln( m 2^q ) = ln(m) + q*ln(2)
    // BASE_2: RR: ln( m 2^q ) = log2(m) + q*log2(2) = log2(m) + q
    const int LODW = slog2<AW>::ceil+1;
    ac_int<LODW,false> leading_one = 0;
    bool no_bits_set;
    ac_int<AW,false> x_int;
    x_int.set_slc(0, x.template slc<AW>(0));
    LOD<AW,LODW>::get(x_int, leading_one, no_bits_set);
    ac_fixed<AW,0,false,AQ,AV> x_norm;
    int dist;
    bool right_shift;
    leading_one += 1; // From index to bits.

    // The expressions below take into account the fact that the input integer width
    // may be negative, or that the input integer width may be greater than bitwidth.
    const int AW_m = (AI > 0) ? (AI > AW ? AI : AW) : (AW - AI);

    if (leading_one > AW - AI) {
      dist = leading_one - (AW - AI);
      ac_fixed<2*AW_m,AW_m,false> x_shift = x;
      x_norm = x_shift >> dist;
      right_shift = true;
    } else {
      dist = (AW - AI) - leading_one;
      ac_fixed<2*AW_m,AW_m,false> x_shift = x;
      x_norm = x_shift << dist;
      right_shift = false;
    }
    ac_fixed<ZW+2,2,true> zc;
    if (BASE == AcLogRR::BASE_E) {
      ac_ln_rr(x_norm, zc);
    }
    if (BASE == AcLogRR::BASE_2) {
      ac_log2_rr(x_norm, zc);
    }
    // Max shift-distance is S = max(AW-AI,AI).
    //   BASE_E: max-offset = S*ln(2)
    //   BASE_2: max-offset = S

    const int OFW = ac::nbits<(AW-AI>(AI-AS)?AW-AI:(AI-AS))>::val;
    ac_fixed<ZW+1+OFW+1,OFW+1,true> offset;

    if (BASE == AcLogRR::BASE_E) {
      offset = ln2_function<ZW+1>();
      if (right_shift) {
        offset *= dist;
      } else {
        offset *= -dist;
      }
    }
    if (BASE == AcLogRR::BASE_2) {
      if (right_shift)
      { offset = dist; }
      else
      { offset = -dist; }
    }
    z = offset;
    z += zc;

#if defined(AC_HCORDIC_H_DEBUG) && !defined(__SYNTHESIS__)
    cout << "x = " << x << endl;
    cout << "leading_one = " << leading_one << endl;
    cout << "x_norm = " << x_norm << endl;
    cout << "dist = " << dist << endl;
    cout << "zc   = " << zc << endl;
    cout << "offset = " << offset << endl;
    cout << "right_shift = " << right_shift << endl;
    cout << "offset.type_name() = " << offset.type_name() << endl;
#endif

  }

  template <int AW, int AI, bool AS, ac_q_mode AQ, ac_o_mode AV,
            int ZW, int ZI, bool ZS, ac_q_mode ZQ, ac_o_mode ZV>
  void ac_log_cordic(const ac_fixed<AW,AI,AS,AQ,AV> &x, ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
  {
    ac_log_<AW,AI,AS,AQ,AV,ZW,ZI,ZS,ZQ,ZV,AcLogRR::BASE_E>(x, z);
  }

  template <int AW, int AI, bool AS, ac_q_mode AQ, ac_o_mode AV,
            int ZW, int ZI, bool ZS, ac_q_mode ZQ, ac_o_mode ZV>
  void ac_log2_cordic(const ac_fixed<AW,AI,AS,AQ,AV> &x, ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
  {
    ac_log_<AW,AI,AS,AQ,AV,ZW,ZI,ZS,ZQ,ZV,AcLogRR::BASE_2>(x, z);
  }

  // The result is expected to have a range which accomodates all
  // resulting values exp(x) for inputs x.
  template <int AW, int AI, bool AS, ac_q_mode AQ, ac_o_mode AV,
            int ZW, int ZI, bool ZS, ac_q_mode ZQ, ac_o_mode ZV>
  void ac_exp_cordic(const ac_fixed<AW,AI,AS,AQ,AV> &x, ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
  {
    // Range reduction: prescale x by 1/ln(2)
    // Let R = trunc(1/ln(2), 2^-(AW - AI))
    // Q = trunc(R*x)
    // Then x = Q*ln(2) + M, where |M| < ln(2)
    //
    // exp(Q*ln(2) + M) = cosh(Q*ln(2) + M) + sinh(Q*ln(2) + M)
    //                  = [exp(Q*ln(2) + M) + exp(-Q*ln(2) - M)]/2 +
    //                    [exp(Q*ln(2) + M) - exp(-Q*ln(2) - M)]/2
    //                  = [exp(M)exp(Q*ln(2) + exp(-M)*exp(-Q*ln(2))]/2 +
    //                    [exp(M)exp(Q*ln(2) - exp(-M)*exp(-Q*ln(2))]/2
    // Given exp(Q*ln(2)) = (exp(ln(2)))^Q = 2^Q
    //                  = [2^Q*exp(M) + 2^-Q*exp(-M)]/2 + [2^Q*exp(M) - 2^-Q*exp(-M)]/2
    // Cancel terms 2^-Q*exp(M)/2 and -2^-Q*exp(-M)/2.
    // Add 2^Q*exp(-M)/2 - 2^Q*exp(-M)/2 = 0
    //                  = [2^Q*exp(M) + 2^Q*exp(-M)]/2 + [2^Q*exp(M) - 2^Q*exp(-M)]/2
    //                  = 2^Q[cosh(M) + sinh(M)]
    //
    // and exp(x) = 2^Q exp(M), |M| < 0.69

    // Output should be at least wide enough to represent exp(2^(AI-AS)), i.e.,
    //  log2(exp(2^(AI-AS)))
    //  = 2^(AI-AS)*log2(exp(1))
    const int QWE = (int)(1.443*(1 << AC_MAX(AI-AS, 0)))+1+1;
    //const int QW = slog2<QWE>::ceil+AS;
    const int QW = ac::nbits<QWE>::val + AS;
    ac_int<QW+1,true> q_int;
    ac_fixed<ZW-ZI+QWE+QW+1,QW+1,true> m;
    ac_fixed<ZW-ZI,1,false> inv_ln2 = inv_ln2_function<ZW-ZI>();
    ac_fixed<QW+1,QW+1,true> q = x*inv_ln2;
    q_int = 0;
    q_int.set_slc(0, q.template slc<QW+1>(0));
    const int MW = ZW-ZI > AW-AI ? ZW-ZI : AW-AI;
    m = ln2_function<MW+QWE>();
    m *= -q;
    m += x;
    // May need up to QWE extra bits forpost-scaling by 2^Q.
    ac_fixed<MW+QWE+3,3,true> zc;
    ac_exp_rr(m, zc);
    // May shift left or right by up to QWE positions.
    ac_fixed<MW+QWE+3+QWE,3+QWE,true> zs = zc;
    if (q_int < 0) {
      q_int = -q_int;
      zs = zs >> q_int;
    } else {
      zs = zs << q_int;
    }
    z = zs;
  }

  // Example range-reduction algorithm which invokes the ac_exp_rr
  // routine.
  template <int AW, int AI, bool AS, ac_q_mode AQ, ac_o_mode AV,
            int ZW, int ZI, bool ZS, ac_q_mode ZQ, ac_o_mode ZV>
  void ac_exp2_cordic(const ac_fixed<AW,AI,AS,AQ,AV> &x, ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
  {
    // If exp(t) = 2^x, then t = x ln(2), i.e., we can compute 2^x
    // using exp(x ln(2)).
    //
    // Recalling our original range reduction form, x = Qln(2) + M.
    // We see that t = x ln(2) = (xI + xF) ln(2), where xI is the
    // integer part of x and xF is the fractional part.
    //
    // Therefore t = xI ln(2) + xF ln(2), i.e., our Q value is
    // directly obtained from xI. Similarly our M value is obtained
    // from multiplying xF by ln(2), i.e., M = xF*ln(2).
    //const int QW = AI-AS;
    const int QW = AC_MAX(AI-AS, 0);
    const int QWE = (1<<QW)+1;
    //ac_int<AI+1,true> q_int;
    ac_int<AC_MAX(AI+1, 2),true> q_int;
    ac_fixed<AW-AI+QWE+1,1,true> m;
    //ac_fixed<AI+1,AI+1,true> q = x;
    ac_fixed<AC_MAX(AI+1, 2), AC_MAX(AI+1, 2),true> q = x;
    q_int.set_slc(0, q.template slc<AC_MAX(AI+1, 2)>(0));
    m = 0;
    m.set_slc(QWE, x.template slc<AW-AI>(0));
    const int MW = ZW-ZI > AW-AI ? ZW-ZI : AW-AI;
    ac_fixed<MW+QWE+1,1,true> mw = m; // Widen m if required.
    mw *= ln2_function<MW+QWE>();
    ac_fixed<ZW-ZI+QWE+3,3,true> zc;
    ac_exp_rr(mw, zc);
    ac_fixed<ZW-ZI+QWE+3+QWE,3+QWE,true> zs = zc;
    if (q_int < 0) {
      q_int = -q_int;
      zs = zs >> q_int;
    } else {
      zs = zs << q_int;
    }
    z = zs;

  }

  // This implementation tries to be as general as possible without
  // taking into account the actual range of arguments 'a', and 'b'. A
  // specific implementation should take full advantage of argument
  // ranges to size the intermediate results of the computation.
  //
  // The following function implements pow(a,b)=a^b using the following
  // decomposition.
  //
  // a^b = exp2(log2(a^b)) = exp2(b*log2(a))
  // t = log(a)
  // q = qI + qF = b*t
  // z = exp2(qI + qF) = exp2(qI)*exp2(qF)
  //
  // This bounds the argument qF of exp2(qF) between 0 and 1.
  // Factor exp2(qI) is accounted for with a final shift.
  template <int AW, int AI, bool AS, ac_q_mode AQ, ac_o_mode AV,
            int BW, int BI, bool BS, ac_q_mode BQ, ac_o_mode BV,
            int ZW, int ZI, bool ZS, ac_q_mode ZQ, ac_o_mode ZV>
  void ac_pow_cordic(const ac_fixed<AW,AI,AS,AQ,AV> &a,
                     const ac_fixed<BW,BI,BS,BQ,BV> &b,
                     ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
  {
    // log2(min(a)) = log(2^-(AW-AI))
    const int TI_MIN = slog2<AW-AI>::ceil;
    // log2(max(a)) = log(2^AI)
    const int TI_MAX = slog2<AI>::ceil;
    const int TI = (TI_MIN < TI_MAX ? TI_MAX : TI_MIN)+1;
    const int TW = ZW-ZI+TI;
    ac_fixed<TW,TI,true> t;
    ac_log2_cordic(a, t);
    const int QI = TI+(BI-BS);
    // A left shift results in loss of precision if the bits to be shifted-in
    // aren't computed. This is possible when 'a' has a large number of fractional
    // positions and 'b' can be negative.
    //
    // * Case b < 0: If TI is obtained from TI_MIN, i.e.,
    //   log2(ulp(a)) < 0, then log2(a)*b can produce a 2^(QI-1) shift distance.
    //
    // * Case b >= 0: If TI is obtained from TI_MAX, i.e.,
    //   log2(max(a)), then b*log2(a) <=  b*max(log2(a)) <= b*2^TI_MAX can produce
    //   a 2^(BI+TI_MAX) shift distance.
    //
    // When argument 'b' is an unsigned value, then the maximum shift distance is
    // dependent on the maximum value of 'a', and not on its number of fractional
    // bits.
    const int SHIFT_W = BS?(1<<(QI-1)):(1<<(BI+TI_MAX));
    const int QW = ZW-ZI+SHIFT_W+QI;
    ac_fixed<QW,QI,true> q = b*t;
    ac_fixed<QI,QI,true> q_i = 0;
    q_i.set_slc(0, q.template slc<QI>(QW-QI));
    ac_fixed<QW-QI,0,false> q_f = 0;
    q_f.set_slc(0, q.template slc<QW-QI>(0));
    ac_int<QI,true> q_int = 0;
    q_int.set_slc(0, q.template slc<QI>(QW-QI));
    // 0 <= q_f < 1
    const int ZCI = 2;
    const int ZCW = QW-QI+2;
    ac_fixed<ZCW,ZCI,false> zc;
    ac_exp2_cordic(q_f, zc);
    ac_fixed<ZCW+SHIFT_W,SHIFT_W+2,false> zc_shift = zc;
    if (q_int < 0) {
      q_int = -q_int;
      zc_shift >>= q_int;
    } else {
      zc_shift <<= q_int;
    }
    z = zc_shift;
  }

}

#endif // _INCLUDED_AC_HCORDIC_H_
