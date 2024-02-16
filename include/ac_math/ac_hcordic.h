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
//*****************************************************************************************
// File: ac_hcordic.h
//
// Description: Hyperbolic CORDIC implementations of synthesizable
//              log/exp/log2/exp2/pow functions for ac_fixed, ac_float,
//              ac_std_float and ac_ieee_float datatypes
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_fixed<16,  8, false, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<32, 10,  true, AC_RND, AC_SAT> output_type;
//    
//    typedef ac_fixed<50,  8, false, AC_RND, AC_SAT> input_2_type;
//    typedef ac_fixed<55, 10,  true, AC_RND, AC_SAT> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_log2_cordic(input, output);
//      // The bitwidth of input_2 and fractional width of output_2 are so large
//      // that they would trigger static_asserts in the code if used as-is, on
//      // account of the default type widths for the temp. variables in
//      // ac_log_. To avoid this, we override the default type widths
//      // by setting the "OR_TF" template parameter to "true" and the "TF_"
//      // parameter to 30. Refer to the documentation for more details.
//      ac_log2_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 2.5;
//      output_type output;
//      input_2_type input_2 = 3.5;
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
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
// Revision History:
//    3.4.4  - [CAT-30998] Added comment examples, improved compiler warnings and input range.
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//    3.4.2  - [CAT-28828] Added floating point support.
//    2.0.10 - Official open-source release as part of the ac_math library.
//
//*****************************************************************************************

#ifndef _INCLUDED_AC_HCORDIC_H_
#define _INCLUDED_AC_HCORDIC_H_

#ifndef __cplusplus
#error C++ is required to include this header file
#endif

// The functions use default template parameters and static_asserts, which are only supported by C++11 or
// later compiler standards. Hence, the user should be informed if they are not using those standards.
#if (defined(__GNUC__) && (__cplusplus < 201103L))
#error Please use C++11 or a later standard for compilation.
#endif
#if (defined(_MSC_VER) && (_MSC_VER < 1920) && !defined(__EDG__))
#error Please use Microsoft VS 2019 or a later standard for compilation.
#endif

#if defined(AC_HCORDIC_H_DEBUG) && !defined(__SYNTHESIS__)
#include <string.h>
#include <iostream>
#endif

#include <ac_int.h>
#include <ac_fixed.h>
#include <ac_float.h>
#include <ac_std_float.h>
#include <ac_math/ac_normalize.h>
#include <ac_math/ac_shift.h>

namespace ac_math
{

  template <bool b>
  struct AcHtrigAssert { 
    //Compile error, no 'test' symbol.
  };
  
  template <>
  struct AcHtrigAssert<true> {
    enum { test };
  };

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

  template<int J>
  struct XtraIters {
    static const int valid = AcHtrigAssert<(J <= 60)>::test;
    enum {
      // sum-of-cordic-terms for {i | i > j + 4*B2+2*B1+B0} < ulp(iteration(j))
      // j >=14 && j <= 60: B0 + B1 + B2 = 5
      // j >= 4 && j <= 13: B0 + B1 + B2 = 2
      // j >= 1 && j <=  3: B0 + B1 + B2 = 1
      B0 = (1 << 0)*int((J >= 1 && J <= 3) || (J >= 14)),
      B1 = (1 << 1)*int(J >= 4 && J <= 13),
      B2 = (1 << 2)*int(J >= 14)
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
    
    const int LW = L + ac::nbits<L>::val;
    
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
  template <int AW, int AI, ac_q_mode AQ, ac_o_mode AV,
            int ZW, int ZI, bool ZS, ac_q_mode ZQ, ac_o_mode ZV>
  void ac_ln_rr(const ac_fixed<AW,AI,false,AQ,AV> &x, ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
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
  template <int AW, int AI, ac_q_mode AQ, ac_o_mode AV,
            int ZW, int ZI, bool ZS, ac_q_mode ZQ, ac_o_mode ZV>
  void ac_log2_rr(const ac_fixed<AW,AI,false,AQ,AV> &x, ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
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
            int ZW,int ZI,ac_q_mode ZQ,ac_o_mode ZV>
  void ac_exp2_rr(const ac_fixed<AW,AI,AS,AQ,AV> &x, ac_fixed<ZW,ZI,false,ZQ,ZV> &z)
  {
    const int EW = AC_MAX(AW, ZW);
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

  template <enum AcLogRR::base BASE,
            bool OR_TF, // Override default fractional bitwidth for temp variables?
            int TF_, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ, // Rounding mode allocated for x_norm_2 variable.
            int AW, int AI, ac_q_mode AQ, ac_o_mode AV,
            int ZW, int ZI, bool ZS, ac_q_mode ZQ, ac_o_mode ZV>
  void ac_log_(const ac_fixed<AW,AI,false,AQ,AV> &x, ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
  {
    // BASE_E: RR: ln(m * 2^q) = ln(m) + q*ln(2)
    // BASE_2: RR: log2(m * 2^q) = log2(m) + q*log2(2) = log2(m) + q

    ac_fixed<AW,0,false,AQ,AV> x_norm;
    int expret = ac_normalize(x, x_norm);
    
    // Both the x_norm variable and the zc variable eventually have to be passed to ac_ln_rr/ac_log2_rr,
    // and we might need to reduce the fractional variables in both to reduce the chances of the
    // compiler throwing an error in the shift_dist function. To do so, use the "OR_TF" and "TF_"
    // template parameters to override the default number of fractional bits assigned to both.
    const int x_norm_2_f_bits = OR_TF ? TF_ : AW;
    
    static_assert(x_norm_2_f_bits <= 44, "Number of fractional bits in x_norm_2 must be no more than 44. Consider overriding the default fractional bits allocated with the OR_TF and TF_ template parameters, so as to satisfy this condition. Refer to the documentation for more information.");
    
    ac_fixed<x_norm_2_f_bits, 0, false, TQ> x_norm_2 = x_norm;
    
    const int zc_f_bits = OR_TF ? TF_ : ZW - ZI;
    
    static_assert(zc_f_bits <= 44, "Number of fractional bits in zc must be no more than 44. Consider overriding the default fractional bits allocated with the OR_TF and TF_ template parameters, so as to satisfy this condition. Refer to the documentation for more information.");

    // Max shift-distance is S = max(AW-AI,AI).
    //   BASE_E: max-offset = S*ln(2)
    //   BASE_2: max-offset = S

    const int OFW = ac::nbits<AC_MAX(AW - AI, AI)>::val;
    
    // The intermediate variable in this case has to account for the fact that
    // for inputs of unity, the magnitude of "zc" can slightly exceed that for
    // "offset", with "zc" being negative. Hence, the addition of "zc" and "offset"
    // results in a slightly negative value, which is enough to cause a highly
    // undesirable output if the output type is unsigned.
    // To take that into account, the intermediate variable storing the addition result
    // is of a saturating type if the output type is unsigned, in order to saturate to 0
    // when negative values are assigned to it.
    const ac_o_mode z_inter_o_mode = ZS ? AC_WRAP : AC_SAT;
    ac_fixed<ZW,ZI,ZS,ZQ,z_inter_o_mode> z_inter;

    if (BASE == AcLogRR::BASE_E) {
      ac_fixed<zc_f_bits+1,1,true> zc;
      // Range Reduced to: 0.5 <= x < 1, -.69 < z < 0
      ac_ln_rr(x_norm_2, zc);
      ac_fixed<zc_f_bits+1+OFW+1,OFW+1,true> offset = ln2_function<zc_f_bits+1>();
      offset *= expret;
      z_inter = offset + zc;
      
      #if defined(AC_HCORDIC_H_DEBUG) && !defined(__SYNTHESIS__)
      std::string base_string = (BASE == AcLogRR::BASE_E) ? "Base e logarithm" : "Base 2 logarithm";
      std::cout << base_string << std::endl;
      std::cout << "x = " << x << std::endl;
      std::cout << "x_norm = " << x_norm << std::endl;
      std::cout << "x_norm_2 = " << x_norm_2 << std::endl;
      std::cout << "zc = " << zc << std::endl;
      std::cout << "offset = " << offset << std::endl;
      std::cout << "z_inter = " << z_inter << std::endl;
      #endif
    }
    if (BASE == AcLogRR::BASE_2) {
      ac_fixed<zc_f_bits+2,2,true> zc;
      // Range reduced to 0.5 <= x < 1, -1 <= z < 1
      ac_log2_rr(x_norm_2, zc);
      ac_fixed<OFW+1,OFW+1,true> offset = expret;
      z_inter = offset + zc;
      
      #if defined(AC_HCORDIC_H_DEBUG) && !defined(__SYNTHESIS__)
      std::string base_string = (BASE == AcLogRR::BASE_E) ? "Base e logarithm" : "Base 2 logarithm";
      std::cout << base_string << std::endl;
      std::cout << "x = " << x << std::endl;
      std::cout << "x_norm = " << x_norm << std::endl;
      std::cout << "x_norm_2 = " << x_norm_2 << std::endl;
      std::cout << "zc = " << zc << std::endl;
      std::cout << "offset = " << offset << std::endl;
      std::cout << "z_inter = " << z_inter << std::endl;
      #endif
    }
    
    z = z_inter;

    #if defined(AC_HCORDIC_H_DEBUG) && !defined(__SYNTHESIS__)
    std::cout << "z = " << z << std::endl;
    #endif

  }

//=========================================================================
// Function: ac_log2_cordic (for ac_fixed)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable log2 function for
//    ac_fixed datatypes.
//
//    This serves as a wrapper around the ac_log_ function for ac_fixed
//    types, calling it with the AcLogRR::BASE_2 enum as a template
//    parameter to select log2 computations.
//
// Usage:
//    See example code at the beginning of this header file for usage.
//
//-------------------------------------------------------------------------

  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for x_norm_2 variable.
            int AW, int AI, ac_q_mode AQ, ac_o_mode AV,
            int ZW, int ZI, bool ZS, ac_q_mode ZQ, ac_o_mode ZV>
  void ac_log2_cordic(const ac_fixed<AW,AI,false,AQ,AV> &x, ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
  {
    ac_log_<AcLogRR::BASE_2,OR_TF,TF_,TQ,AW,AI,AQ,AV,ZW,ZI,ZS,ZQ,ZV>(x, z);
  }

//=========================================================================
// Function: ac_log_cordic (for ac_fixed)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable log function for
//    ac_fixed datatypes.
//
//    This serves as a wrapper around the ac_log_ function for ac_fixed
//    types, calling it with the AcLogRR::BASE_E enum as a template
//    parameter to select loge computations.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_fixed<16,  8, false, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<32, 10,  true, AC_RND, AC_SAT> output_type;
//    
//    typedef ac_fixed<50,  8, false, AC_RND, AC_SAT> input_2_type;
//    typedef ac_fixed<55, 10,  true, AC_RND, AC_SAT> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_log_cordic(input, output);
//      // The bitwidth of input_2 and fractional width of output_2 are so large
//      // that they would trigger static_asserts in the code if used as-is, on
//      // account of the default type widths for the temp. variables in
//      // ac_log_. To avoid this, we override the default type widths
//      // by setting the "OR_TF" template parameter to "true" and the "TF_"
//      // parameter to 30. Refer to the documentation for more details.
//      ac_log_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 2.5;
//      output_type output;
//      input_2_type input_2 = 3.5;
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for x_norm_2 variable.
            int AW, int AI, ac_q_mode AQ, ac_o_mode AV,
            int ZW, int ZI, bool ZS, ac_q_mode ZQ, ac_o_mode ZV>
  void ac_log_cordic(const ac_fixed<AW,AI,false,AQ,AV> &x, ac_fixed<ZW,ZI,ZS,ZQ,ZV> &z)
  {
    ac_log_<AcLogRR::BASE_E,OR_TF,TF_,TQ,AW,AI,AQ,AV,ZW,ZI,ZS,ZQ,ZV>(x, z);
  }

  template <enum AcLogRR::base BASE,
            bool OR_TF, // Override default fractional bitwidth for temp variables?
            int TF_, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ, // Rounding mode allocated for x_norm_2 variable.
            int AW, int AI, int AE, ac_q_mode AQ,
            int ZW, int ZI, int ZE, ac_q_mode ZQ>
  void ac_log_(const ac_float<AW,AI,AE,AQ> &x, ac_float<ZW,ZI,ZE,ZQ> &z)
  {
    // BASE_E: RR: ln(m * 2^q * 2^exp) = ln(m) + (q + exp)*ln(2)
    // BASE_2: RR: log2(m * 2^q * 2^exp) = log2(m) + (q + exp)*log2(2) = log2(m) + q + exp
    
    #ifdef ASSERT_ON_INVALID_INPUT
    AC_ASSERT(x >= 0, "Negative input not supported.");
    #endif
    
    ac_fixed<AW - 1, AI - 1, false> x_mant = x.mantissa();
    ac_int<AE, true> x_exp = x.exp();

    ac_fixed<AW - 1, 0, false> x_norm;
    int x_norm_exp = ac_normalize(x_mant, x_norm);
    
    // Both the x_norm variable and the zc variable eventually have to be passed to ac_ln_rr/ac_log2_rr,
    // and we might need to reduce the fractional variables in both to reduce the chances of the
    // compiler throwing an error in the shift_dist function. To do so, use the "OR_TF" and "TF_"
    // template parameters to override the default number of fractional bits assigned to both.
    
    const int x_norm_2_f_bits = OR_TF ? TF_ : AW - 1;
    
    static_assert(x_norm_2_f_bits <= 44, "Number of fractional bits in x_norm_2 must be no more than 44. Consider overriding the default fractional bits allocated with the OR_TF and TF_ template parameters, so as to satisfy this condition. Refer to the documentation for more information.");
    
    ac_fixed<x_norm_2_f_bits, 0, false, TQ> x_norm_2 = x_norm;
    
    const int zc_f_bits = OR_TF ? TF_ : ZW - ZI;
    
    static_assert(zc_f_bits <= 44, "Number of fractional bits in zc must be no more than 44. Consider overriding the default fractional bits allocated with the OR_TF and TF_ template parameters, so as to satisfy this condition. Refer to the documentation for more information.");
    
    const int OFW = ac::nbits<AC_MAX(AI - 1, AW - AI - 1) + (1 << (AE - 1))>::val;
    
    ac_float<ZW, ZI, ZE, ZQ> z_inter;
    
    if (BASE == AcLogRR::BASE_E) {
      ac_fixed<zc_f_bits + 1, 1, true> zc;
      // Range Reduced to: 0.5 <= x < 1, -.69 < z < 0
      ac_ln_rr(x_norm_2, zc);
      ac_fixed<zc_f_bits + 1 + OFW + 1, OFW + 1, true> offset = ln2_function<zc_f_bits + 1>();
      offset *= (x_exp + x_norm_exp);
      z_inter = ac_float<ZW, ZI, ZE, ZQ>(offset + zc);
    }
    if (BASE == AcLogRR::BASE_2) {
      ac_fixed<zc_f_bits + 2, 2, true> zc;
      // Range reduced to 0.5 <= x < 1, -1 <= z < 1
      ac_log2_rr(x_norm_2, zc);
      ac_fixed<OFW + 1, OFW + 1, true> offset = x_exp + x_norm_exp;
      z_inter = ac_float<ZW, ZI, ZE, ZQ>(offset + zc);
    }
    
    z = z_inter;
    
    #if defined(AC_HCORDIC_H_DEBUG) && !defined(__SYNTHESIS__)
    std::string base_string = (BASE == AcLogRR::BASE_E) ? "Base e logarithm" : "Base 2 logarithm";
    std::cout << base_string << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "x_norm = " << x_norm << std::endl;
    std::cout << "x_norm_exp = " << x_norm_exp << std::endl;
    std::cout << "x_exp = " << x_exp << std::endl;
    std::cout << "z = " << z << std::endl;
    #endif
  }

//=========================================================================
// Function: ac_log2_cordic (for ac_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable log2 function for
//    ac_float datatypes.
//
//    This serves as a wrapper around the ac_log_ function for ac_float
//    types, calling it with the AcLogRR::BASE_2 enum as a template
//    parameter to select log2 computations.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_float<16, 2, 5, AC_RND> input_type;
//    typedef ac_float<32, 2, 8, AC_RND> output_type;
//    
//    typedef ac_float<50, 2, 8, AC_RND> input_2_type;
//    typedef ac_float<64, 2, 10, AC_RND> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_log2_cordic(input, output);
//      // The bitwidth of input_2 and fractional width of output_2's mantissa
//      // are so large that they would trigger static_asserts in the code if
//      // used as-is, on account of the default type widths for the temp.
//      // variables in ac_log_. To avoid this, we override the default type
//      // widths by setting the "OR_TF" template parameter to "true" and the
//      // "TF_" parameter to 30. Refer to the documentation for more details.
//      ac_log2_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 2.5;
//      output_type output;
//      input_2_type input_2 = 3.5;
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
            int AW, int AI, int AE, ac_q_mode AQ,
            int ZW, int ZI, int ZE, ac_q_mode ZQ>
  void ac_log2_cordic(const ac_float<AW, AI, AE, AQ> &x, ac_float<ZW, ZI, ZE, ZQ> &z)
  {
    ac_log_<AcLogRR::BASE_2, OR_TF, TF_, TQ>(x, z);
  }

//=========================================================================
// Function: ac_log_cordic (for ac_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable log function for
//    ac_float datatypes.
//
//    This serves as a wrapper around the ac_log_ function for ac_float
//    types, calling it with the AcLogRR::BASE_E enum as a template
//    parameter to select log computations.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_float<16, 2, 5, AC_RND> input_type;
//    typedef ac_float<32, 2, 8, AC_RND> output_type;
//    
//    typedef ac_float<50, 2, 8, AC_RND> input_2_type;
//    typedef ac_float<64, 2, 10, AC_RND> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_log_cordic(input, output);
//      // The bitwidth of input_2 and fractional width of output_2's mantissa
//      // are so large that they would trigger static_asserts in the code if
//      // used as-is, on account of the default type widths for the temp.
//      // variables in ac_log_. To avoid this, we override the default type
//      // widths by setting the "OR_TF" template parameter to "true" and the
//      // "TF_" parameter to 30. Refer to the documentation for more details.
//      ac_log_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 2.5;
//      output_type output;
//      input_2_type input_2 = 3.5;
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for x_norm_2 variable in ac_log_.
            int AW, int AI, int AE, ac_q_mode AQ,
            int ZW, int ZI, int ZE, ac_q_mode ZQ>
  void ac_log_cordic(const ac_float<AW, AI, AE, AQ> &x, ac_float<ZW, ZI, ZE, ZQ> &z)
  {
    ac_log_<AcLogRR::BASE_E, OR_TF, TF_, TQ>(x, z);
  }

//=========================================================================
// Function: ac_log_cordic (for ac_std_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable log function for
//    ac_std_float datatypes.
//
// Usage:
//    A sample testbench and its implementation look like this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_std_float<32, 8> input_type;
//    typedef ac_std_float<32, 8> output_type;
//    
//    typedef ac_std_float<64, 11> input_2_type;
//    typedef ac_std_float<64, 11> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_log_cordic(input, output);
//      // The bitwidth of input_2 and fractional width of output_2's mantissa
//      // are so large that they would trigger static_asserts in the code if
//      // used as-is, on account of the default type widths for the temp.
//      // variables in ac_log_. To avoid this, we override the default type
//      // widths by setting the "OR_TF" template parameter to "true" and the
//      // "TF_" parameter to 30. Refer to the documentation for more details.
//      ac_log_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input(2.5);
//      output_type output;
//      input_2_type input_2(3.5);
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------
  
  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for x_norm_2 variable in ac_log_.
            int AW, int AE, int ZW, int ZE>
  void ac_log_cordic(const ac_std_float<AW, AE> &x, ac_std_float<ZW, ZE> &z)
  {
    ac_float<ZW - ZE + 1, 2, ZE> z_ac_fl; // Equivalent ac_float representation for output.
    ac_log_cordic<OR_TF, TF_, TQ>(x.to_ac_float(), z_ac_fl); // Call ac_float version.
    ac_std_float<ZW, ZE> z_temp(z_ac_fl); // Convert output ac_float to ac_std_float.
    
    #ifdef AC_HCORDIC_NAN_SUPPORTED
    // If the input is +/-nan, the output will be +/-nan. This mirrors the behavior of the
    // log2 library in math.h.
    //
    // If the input is negative, output will be +nan, mirroring the math.h log2 library.
    if (x.isnan() || x.signbit()) {
      z_temp = z.nan();
      if (x.isnan()) {
        z_temp.set_signbit(x.signbit());
      }
    }
    #endif
    
    z = z_temp;
  }

//=========================================================================
// Function: ac_log2_cordic (for ac_std_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable log2 function for
//    ac_std_float datatypes.
//
// Usage:
//    A sample testbench and its implementation look like this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_std_float<32, 8> input_type;
//    typedef ac_std_float<32, 8> output_type;
//    
//    typedef ac_std_float<64, 11> input_2_type;
//    typedef ac_std_float<64, 11> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_log2_cordic(input, output);
//      // The bitwidth of input_2 and fractional width of output_2's mantissa
//      // are so large that they would trigger static_asserts in the code if
//      // used as-is, on account of the default type widths for the temp.
//      // variables in ac_log_. To avoid this, we override the default type
//      // widths by setting the "OR_TF" template parameter to "true" and the
//      // "TF_" parameter to 30. Refer to the documentation for more details.
//      ac_log2_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input(2.5);
//      output_type output;
//      input_2_type input_2(3.5);
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for x_norm_2 variable in ac_log_.
            int AW, int AE, int ZW, int ZE>
  void ac_log2_cordic(const ac_std_float<AW, AE> &x, ac_std_float<ZW, ZE> &z)
  {
    ac_float<ZW - ZE + 1, 2, ZE> z_ac_fl; // Equivalent ac_float representation for output.
    ac_log2_cordic<OR_TF, TF_, TQ>(x.to_ac_float(), z_ac_fl); // Call ac_float version.
    ac_std_float<ZW, ZE> z_temp(z_ac_fl); // Convert output ac_float to ac_std_float.
    
    #ifdef AC_HCORDIC_NAN_SUPPORTED
    // If the input is +/-nan, the output will be +/-nan. This mirrors the behavior of the
    // log2 library in math.h.
    //
    // If the input is negative, output will be +nan, mirroring the math.h log2 library.
    if (x.isnan() || x.signbit()) {
      z_temp = z.nan();
      if (x.isnan()) {
        z_temp.set_signbit(x.signbit());
      }
    }
    #endif
    
    z = z_temp;
  }

//=========================================================================
// Function: ac_log_cordic (for ac_ieee_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable log function for
//    ac_ieee_float datatypes.
//
// Usage:
//    A sample testbench and its implementation look like this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_ieee_float<binary32> input_type;
//    typedef ac_ieee_float<binary32> output_type;
//    
//    typedef ac_ieee_float<binary64> input_2_type;
//    typedef ac_ieee_float<binary64> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_log_cordic(input, output);
//      // The bitwidth of input_2 and fractional width of output_2's mantissa
//      // are so large that they would trigger static_asserts in the code if
//      // used as-is, on account of the default type widths for the temp.
//      // variables in ac_log_. To avoid this, we override the default type
//      // widths by setting the "OR_TF" template parameter to "true" and the
//      // "TF_" parameter to 30. Refer to the documentation for more details.
//      ac_log_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input(2.5);
//      output_type output;
//      input_2_type input_2(3.5);
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------
  
  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for x_norm_2 variable in ac_log_.
            ac_ieee_float_format Format, ac_ieee_float_format outFormat>
  void ac_log_cordic(const ac_ieee_float<Format> &x, ac_ieee_float<outFormat> &z)
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int ZW = T_out::width;
    const int ZE = T_out::e_width;
  
    ac_float<ZW - ZE + 1, 2, ZE> z_ac_fl; // Equivalent ac_float representation for output.
    ac_log_cordic<OR_TF, TF_, TQ>(x.to_ac_float(), z_ac_fl); // Call ac_float version.
    T_out z_temp(z_ac_fl); // Convert output ac_float to ac_ieee_float.
    
    #ifdef AC_HCORDIC_NAN_SUPPORTED
    // If the input is +/-nan, the output will be +/-nan. This mirrors the behavior of the log2
    // library in math.h. For the same reason, a negative input will result in a +nan output.
    if (x.isnan() || x.signbit()) {
      z_temp = z.nan();
      if (x.isnan()) {
        z_temp.set_signbit(x.signbit());
      }
    }
    #endif
    
    z = z_temp;
  }

//=========================================================================
// Function: ac_log2_cordic (for ac_ieee_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable log2 function for
//    ac_ieee_float datatypes.
//
// Usage:
//    A sample testbench and its implementation look like this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_ieee_float<binary32> input_type;
//    typedef ac_ieee_float<binary32> output_type;
//    
//    typedef ac_ieee_float<binary64> input_2_type;
//    typedef ac_ieee_float<binary64> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_log2_cordic(input, output);
//      // The bitwidth of input_2 and fractional width of output_2's mantissa
//      // are so large that they would trigger static_asserts in the code if
//      // used as-is, on account of the default type widths for the temp.
//      // variables in ac_log_. To avoid this, we override the default type
//      // widths by setting the "OR_TF" template parameter to "true" and the
//      // "TF_" parameter to 30. Refer to the documentation for more details.
//      ac_log2_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input(2.5);
//      output_type output;
//      input_2_type input_2(3.5);
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------
  
  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for x_norm_2 variable in ac_log_.
            ac_ieee_float_format Format, ac_ieee_float_format outFormat>
  void ac_log2_cordic(const ac_ieee_float<Format> &x, ac_ieee_float<outFormat> &z)
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int ZW = T_out::width;
    const int ZE = T_out::e_width;
  
    ac_float<ZW - ZE + 1, 2, ZE> z_ac_fl; // Equivalent ac_float representation for output.
    ac_log2_cordic<OR_TF, TF_, TQ>(x.to_ac_float(), z_ac_fl); // Call ac_float version.
    T_out z_temp(z_ac_fl); // Convert output ac_float to ac_ieee_float.
    
    #ifdef AC_HCORDIC_NAN_SUPPORTED
    // If the input is +/-nan, the output will be +/-nan. This mirrors the behavior of the log2
    // library in math.h. For the same reason, a negative input will result in a +nan output.
    if (x.isnan() || x.signbit()) {
      z_temp = z.nan();
      if (x.isnan()) {
        z_temp.set_signbit(x.signbit());
      }
    }
    #endif
    
    z = z_temp;
  }

//=========================================================================
// Function: ac_exp2_cordic (for ac_fixed)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable exp2 function for
//    ac_fixed datatypes.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_fixed<16,  8, true, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<32, 16, false, AC_RND, AC_SAT> output_type;
//    
//    // Input fractional bits must not exceed 45.
//    typedef ac_fixed<55, 10, true, AC_RND, AC_SAT> input_2_type;
//    typedef ac_fixed<65, 20, false, AC_RND, AC_SAT> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_exp2_cordic(input, output);
//      // output_2 has a large number of fractional bits, so large that
//      // it would trigger a static_assert in the code if used as-is, on
//      // account of the default fractional width for a temp. variable in
//      // ac_exp2_cordic. To avoid this, we override the default width
//      // by setting the "OR_TF" template parameter to "true" and the "TF_"
//      // parameter to 30. Refer to the documentation for more details.
//      ac_exp2_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 2.5;
//      output_type output;
//      input_2_type input_2 = 3.5;
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template <bool OR_TF = false, // Override default fractional bitwidth for temp variable?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variable.
            int AW, int AI, bool AS, ac_q_mode AQ, ac_o_mode AV,
            int ZW, int ZI, ac_q_mode ZQ, ac_o_mode ZV>
  void ac_exp2_cordic(const ac_fixed<AW,AI,AS,AQ,AV> &x, ac_fixed<ZW,ZI,false,ZQ,ZV> &z)
  {    
    // Split the input into integer (xI) and fractional (xF) parts.
    // 2^(xI + xF) = (2^xI)*(2^xF)
    // 2^xF is computed using ac_exp2_rr.
    // Multiplying with 2^xI can be done merely by left-shifting.
    
    // If (AW - AI) > 45, we'll get a compiler error further down the line, in the shift_dist function.
    static_assert(AW - AI <= 45, "Number of input fraction bits must not exceed 45.");

    ac_fixed<AC_MAX(AW - AI, 1), 0, false> xF = 0.0;
    
    if (AW > AI) {
      xF.set_slc(0, x.template slc<AC_MAX(AW - AI, 1)>(0));
    }
    
    // If OR_TF is set to true, fractional bitwidth used in temp variables (EXP2_TF) is set to TF_.
    // By default, however, EXP2_TF is set to the number of fractional bits in the output.
    const int EXP2_TF = OR_TF ? TF_ : ZW - ZI;
    
    // If EXP2_TF > 44, (EXP2_TF + 1) > 45, which will cause a compiler error further down the line, in
    // the shift_dist function.
    static_assert(EXP2_TF <= 44, "EXP2_TF must be lesser than 45. Consider overriding the default value with the OR_TF and TF_ template parameters in your function call, to satisfy these conditions. Refer to the documentation for more information.");
    
    // Caclculating 2^xF
    ac_fixed<EXP2_TF + 1, 1, false, TQ> exp2_xF;
    ac_exp2_rr(xF, exp2_xF);
    
    ac_fixed<ZW,ZI,false,ZQ,ZV> z_temp;
    // Multiply (2^xF) with (2^xI) by left-shifting.
    ac_math::ac_shift_left(exp2_xF, x.to_int(), z_temp);
    
    z = z_temp;
  }

//=========================================================================
// Function: ac_exp_cordic (for ac_fixed)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable exp function for
//    ac_fixed datatypes.
//
//    This function multiplies the input by 1/ln(2) and passes the product
//    to ac_exp2_cordic. The exp2 output of this product is the exp output
//    of the input, on account of the change of base property for exponentials.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_fixed<16,  8, true, AC_RND, AC_SAT> input_type;
//    typedef ac_fixed<32, 16, false, AC_RND, AC_SAT> output_type;
//    
//    typedef ac_fixed<32, 16, true, AC_RND, AC_SAT> input_2_type;
//    typedef ac_fixed<65, 20, false, AC_RND, AC_SAT> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_exp_cordic(input, output);
//      // output_2 has a large number of fractional bits, so large that
//      // it would trigger a static_assert in the code if used as-is, on
//      // account of the default fractional width for temp. variables in
//      // ac_exp2_cordic. To avoid this, we override the default width
//      // by setting the "OR_TF" template parameter to "true" and the
//      // "TF_" parameter to 30.
//      // Refer to the documentation for more details.
//      ac_exp_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 2.5;
//      output_type output;
//      input_2_type input_2 = 3.5;
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
            int AW, int AI, bool AS, ac_q_mode AQ, ac_o_mode AV,
            int ZW, int ZI, ac_q_mode ZQ, ac_o_mode ZV>
  void ac_exp_cordic(const ac_fixed<AW,AI,AS,AQ,AV> &x, ac_fixed<ZW,ZI,false,ZQ,ZV> &z)
  {
    // If OR_TF is set to true, fractional bitwidth used in temp variables (EXP_TF) is set to TF_.
    // By default, however, EXP_TF is set to the number of fractional bits in the output.
    const int EXP_TF = OR_TF ? TF_ : ZW - ZI;
    
    // If EXP_TF > 44, (EXP_TF + 1) > 45, which will cause a compiler error further down the line, in the
    // shift_dist function.
    static_assert(EXP_TF <= 44, "EXP_TF must be lesser than 45. Consider overriding the default value with the OR_TF and TF_ template parameters in your function call, to satisfy these conditions. Refer to the documentation for more information.");
    
    // Multiply input with 1/ln(2), and pass the product to ac_exp2_cordic.
    ac_fixed<EXP_TF + 1, 1, false> inv_ln2 = inv_ln2_function<EXP_TF + 1>();
    ac_fixed<EXP_TF + AI + 1, AI + 1, AS, TQ> x_to_exp2 = x*inv_ln2;
    ac_exp2_cordic<OR_TF, TF_, TQ>(x_to_exp2, z);
  }

  // In order to fully leverage the increased range provided by floating point types, the
  // floating point implementations for the cordic exponential cordic functions, i.e. ac_exp_cordic and
  // ac_exp2_cordic, operate differently than their fixed point counterparts. Both of them convert the
  // floating point input to an intermediate fixed point input, the fractional part of which is passed to
  // ac_exp2_rr function. The output corresponding to the fractional part is then de-normalized to
  // produce the final output.
  //
  // The floating point ac_exp_cordic function also multiplies the input with 1/ln(2) before it is
  // normalized and passed to ac_exp2_rr, so as to utilize the change-of-base property,
  // i.e. exp(x) = exp2(x/ln(2)). The fixed point version does not do that.
  //
  // Since the only difference between ac_exp2_cordic and ac_exp_cordic is the above multiplication, the
  // core code for both can be merged into a single generic function ac_exp_float_ which switches between
  // multiplying and not multiplying with 1/ln(2) based on the "base" enum.

  template <enum AcLogRR::base BASE,
            bool OR_TF, // Override default fractional bitwidth for temp variables?
            int TF_, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ, // Rounding mode allocated for temporary variable.
            int AW, int AI, int AE, ac_q_mode AQ,
            int ZW, int ZI, int ZE, ac_q_mode ZQ>
  void ac_exp_float_(const ac_float<AW,AI,AE,AQ> &x, ac_float<ZW,ZI,ZE,ZQ> &z)
  {
    const int z_exp_min = 1 << (ZE - 1);
    const int TI = ac::nbits<z_exp_min + AC_MAX(ZI - 1, ZW - ZI)>::val + 1;
    // If OR_TF is set to true, fractional bitwidth used in temp variables is set to TF_.
    const int TF = OR_TF ? TF_ : ZW - ZI;
    
    // If TF < 1, it'll cause an issue with bit-slicing operations further down the line.
    // If TF > 44, (TF + 1) > 45, which will cause a compiler error further down the line, in the
    // shift_dist function.
    static_assert(TF >= 1 && TF <= 44, "TF must be a positive integer and lesser than 45. Consider overriding the default value with the OR_TF and TF_ template parameters in your function call, to satisfy these conditions. Refer to the documentation for more information.");
    
    ac_fixed<AW, AI, true> x_mant = x.mantissa();
    ac_int<AE, true> x_exp = x.exp();
    
    ac_fixed<TF + TI, TI, true, TQ, AC_SAT> input_to_exp2_fxpt;
    if (BASE == AcLogRR::BASE_2) {
      // Convert ac_float to ac_fixed variable, through the following left-shift operation:
      // mantissa << exp. The precision of input_to_exp2_fxpt is such that it will saturate a little
      // above the input which would cause the floating point output to saturate if we were using an
      // accurate math function to calculate the exp2 output.
      ac_math::ac_shift_left(x_mant, x_exp.to_int(), input_to_exp2_fxpt);
    } else {
      // This step utilizes almost the same left-shifting operation as the one for base 2 exponentials,
      // except the input mantissa is also multiplied by 1/ln(2).
      ac_fixed<TF + 1, 1, false> inv_ln2 = inv_ln2_function<TF + 1>();
      ac_math::ac_shift_left(x_mant*inv_ln2, x_exp.to_int(), input_to_exp2_fxpt);
    }
    
    ac_fixed<TF, 0, false> input_to_exp2_fxpt_frac;
    input_to_exp2_fxpt_frac.set_slc(0, input_to_exp2_fxpt.template slc<TF>(0));
    ac_fixed<TF + 1, 1, false, TQ> output_to_fxpt_frac;
    ac_exp2_rr(input_to_exp2_fxpt_frac, output_to_fxpt_frac);
    
    ac_int<TI, true> input_to_exp2_fxpt_int = input_to_exp2_fxpt.to_int();
    ac_float<ZW, ZI, ZE, ZQ> z_temp(output_to_fxpt_frac, input_to_exp2_fxpt_int, true);
    
    z = z_temp;
    
    #if defined(AC_HCORDIC_H_DEBUG) && !defined(__SYNTHESIS__)
    std::cout << "x = " << x << std::endl;
    std::cout << "input_to_exp2_fxpt.type_name() : " << input_to_exp2_fxpt.type_name() << std::endl;
    std::cout << "input_to_exp2_fxpt = " << input_to_exp2_fxpt << std::endl;
    std::cout << "input_to_exp2_fxpt_frac = " << input_to_exp2_fxpt_frac << std::endl;
    std::cout << "input_to_exp2_fxpt_int = " << input_to_exp2_fxpt_int << std::endl;
    std::cout << "output_to_fxpt_frac = " << output_to_fxpt_frac << std::endl;
    std::cout << "z = " << z << std::endl;
    #endif
  }

//=========================================================================
// Function: ac_exp_cordic (for ac_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable exp function for
//    ac_float datatypes.
//
//    ac_exp_cordic for ac_float types is a wrapper, the heavy lifting is
//    done by ac_exp_float_, ac_exp2_rr and all the other functions called
//    by ac_exp2_rr.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_float<16, 2, 5, AC_RND> input_type;
//    typedef ac_float<32, 2, 7, AC_RND> output_type;
//    
//    typedef ac_float<32, 2,  8, AC_RND> input_2_type;
//    typedef ac_float<50, 2, 10, AC_RND> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_exp_cordic(input, output);
//      // output_2's mantissa has a large number of fractional bits, so
//      // large that it would trigger a static_assert in the code if used
//      // as-is, on account of the default fractional width for temp.
//      // variables in ac_exp_float_. To avoid this, we override the
//      // default width by setting the "OR_TF" template parameter to
//      // "true" and the "TF_" parameter to 30.
//      // Refer to the documentation for more details.
//      ac_exp_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 2.5;
//      output_type output;
//      input_2_type input_2 = 3.5;
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
            int AW, int AI, int AE, ac_q_mode AQ,
            int ZW, int ZI, int ZE, ac_q_mode ZQ>
  void ac_exp_cordic(const ac_float<AW,AI,AE,AQ> &x, ac_float<ZW,ZI,ZE,ZQ> &z)
  {
    // Reuse the BASE_E enum from the AcLogRR struct, to switch to the base e exponential.
    ac_exp_float_<AcLogRR::BASE_E, OR_TF, TF_, TQ>(x, z);
  }

//=========================================================================
// Function: ac_exp2_cordic (for ac_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable exp2 function for
//    ac_float datatypes.
//
//    ac_exp2_cordic for ac_float types is a wrapper, the heavy lifting is
//    done by ac_exp_float_, ac_exp2_rr and all the other functions called
//    by ac_exp2_rr.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_float<16, 2, 5, AC_RND> input_type;
//    typedef ac_float<32, 2, 7, AC_RND> output_type;
//    
//    typedef ac_float<32, 2,  8, AC_RND> input_2_type;
//    typedef ac_float<50, 2, 10, AC_RND> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_exp2_cordic(input, output);
//      // output_2's mantissa has a large number of fractional bits, so
//      // large that it would trigger a static_assert in the code if used
//      // as-is, on account of the default fractional width for temp.
//      // variables in ac_exp_float_. To avoid this, we override the
//      // default width by setting the "OR_TF" template parameter to
//      // "true" and the "TF_" parameter to 30.
//      // Refer to the documentation for more details.
//      ac_exp2_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input = 2.5;
//      output_type output;
//      input_2_type input_2 = 3.5;
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  // ac_exp2_cordic for ac_float types is a wrapper, the heavy lifting is done by ac_exp_float_, 
  // ac_exp2_rr and all the other functions called by ac_exp2_rr.
  
  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
            int AW, int AI, int AE, ac_q_mode AQ,
            int ZW, int ZI, int ZE, ac_q_mode ZQ>
  void ac_exp2_cordic(const ac_float<AW,AI,AE,AQ> &x, ac_float<ZW,ZI,ZE,ZQ> &z)
  {
    // Reuse the BASE_E enum from the AcLogRR struct, to switch to the base 2 exponential.
    ac_exp_float_<AcLogRR::BASE_2, OR_TF, TF_, TQ>(x, z);
  }

//=========================================================================
// Function: ac_exp_cordic (for ac_std_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable exp function for
//    ac_std_float datatypes.
//
//    ac_exp_cordic for ac_std_float types is a wrapper, the heavy lifting is
//    done by ac_exp_float_, ac_exp2_rr and all the other functions called
//    by ac_exp2_rr.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_std_float<32, 8> input_type;
//    typedef ac_std_float<32, 8> output_type;
//    
//    typedef ac_std_float<64, 11> input_2_type;
//    typedef ac_std_float<64, 11> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_exp_cordic(input, output);
//      // output_2's mantissa has a large number of fractional bits, so
//      // large that it would trigger a static_assert in the code if used
//      // as-is, on account of the default fractional width for temp.
//      // variables in ac_exp_float_. To avoid this, we override the
//      // default width by setting the "OR_TF" template parameter to
//      // "true" and the "TF_" parameter to 30.
//      // Refer to the documentation for more details.
//      ac_exp_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input(2.5);
//      output_type output;
//      input_2_type input_2(3.5);
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  // The ac_std_float implementation for ac_exp_cordic uses temporary ac_float variables which are in
  // turn passed to the ac_float implementation.
  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
            int AW, int AE, int ZW, int ZE>
  void ac_exp_cordic(const ac_std_float<AW, AE> &x, ac_std_float<ZW, ZE> &z)
  {
    ac_float<ZW - ZE + 1, 2, ZE> z_ac_fl; // Equivalent ac_float representation for output.
    ac_exp_cordic<OR_TF, TF_, TQ>(x.to_ac_float(), z_ac_fl); // Call ac_float version.
    ac_std_float<ZW, ZE> z_temp(z_ac_fl); // Convert output ac_float to ac_std_float.
    
    #ifdef AC_HCORDIC_NAN_SUPPORTED
    // If the input is +/-nan, the output will be +/-nan. This mirrors the behavior of the
    // exp library in math.h.
    if (x.isnan()) {
      z_temp = z.nan();
      z_temp.set_signbit(x.signbit());
    }
    #endif
    
    z = z_temp;
  }

//=========================================================================
// Function: ac_exp2_cordic (for ac_std_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable exp2 function for
//    ac_std_float datatypes.
//
//    ac_exp2_cordic for ac_std_float types is a wrapper, the heavy lifting is
//    done by ac_exp_float_, ac_exp2_rr and all the other functions called
//    by ac_exp2_rr.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_std_float<32, 8> input_type;
//    typedef ac_std_float<32, 8> output_type;
//    
//    typedef ac_std_float<64, 11> input_2_type;
//    typedef ac_std_float<64, 11> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_exp2_cordic(input, output);
//      // output_2's mantissa has a large number of fractional bits, so
//      // large that it would trigger a static_assert in the code if used
//      // as-is, on account of the default fractional width for temp.
//      // variables in ac_exp_float_. To avoid this, we override the
//      // default width by setting the "OR_TF" template parameter to
//      // "true" and the "TF_" parameter to 30.
//      // Refer to the documentation for more details.
//      ac_exp2_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input(2.5);
//      output_type output;
//      input_2_type input_2(3.5);
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  // The ac_std_float implementation for ac_exp2_cordic uses temporary ac_float variables which are in
  // turn passed to the ac_float implementation.
  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
            int AW, int AE, int ZW, int ZE>
  void ac_exp2_cordic(const ac_std_float<AW, AE> &x, ac_std_float<ZW, ZE> &z)
  {
    ac_float<ZW - ZE + 1, 2, ZE> z_ac_fl; // Equivalent ac_float representation for output.
    ac_exp2_cordic<OR_TF, TF_, TQ>(x.to_ac_float(), z_ac_fl); // Call ac_float version.
    ac_std_float<ZW, ZE> z_temp(z_ac_fl); // Convert output ac_float to ac_std_float.
    
    #ifdef AC_HCORDIC_NAN_SUPPORTED
    // If the input is +/-nan, the output will be +/-nan. This mirrors the behavior of the
    // exp2 library in math.h.
    if (x.isnan()) {
      z_temp = z.nan();
      z_temp.set_signbit(x.signbit());
    }
    #endif
    
    z = z_temp;
  }

//=========================================================================
// Function: ac_exp_cordic (for ac_ieee_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable exp function for
//    ac_ieee_float datatypes.
//
//    ac_exp_cordic for ac_ieee_float types is a wrapper, the heavy lifting is
//    done by ac_exp_float_, ac_exp2_rr and all the other functions called
//    by ac_exp2_rr.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_ieee_float<binary32> input_type;
//    typedef ac_ieee_float<binary32> output_type;
//    
//    typedef ac_ieee_float<binary64> input_2_type;
//    typedef ac_ieee_float<binary64> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_exp_cordic(input, output);
//      // output_2's mantissa has a large number of fractional bits, so
//      // large that it would trigger a static_assert in the code if used
//      // as-is, on account of the default fractional width for temp.
//      // variables in ac_exp_float_. To avoid this, we override the
//      // default width by setting the "OR_TF" template parameter to
//      // "true" and the "TF_" parameter to 30.
//      // Refer to the documentation for more details.
//      ac_exp_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input(2.5);
//      output_type output;
//      input_2_type input_2(3.5);
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  // The ac_ieee_float implementation for ac_exp_cordic uses temporary ac_float variables which are in
  // turn passed to the ac_float implementation.
  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
            ac_ieee_float_format Format, ac_ieee_float_format outFormat>
  void ac_exp_cordic(const ac_ieee_float<Format> &x, ac_ieee_float<outFormat> &z)
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int ZW = T_out::width;
    const int ZE = T_out::e_width;
  
    ac_float<ZW - ZE + 1, 2, ZE> z_ac_fl; // Equivalent ac_float representation for output.
    ac_exp_cordic<OR_TF, TF_, TQ>(x.to_ac_float(), z_ac_fl); // Call ac_float version.
    T_out z_temp(z_ac_fl); // Convert output ac_float to ac_ieee_float.
    
    #ifdef AC_HCORDIC_NAN_SUPPORTED
    // If the input is +/-nan, the output will be +/-nan. This mirrors the behavior of the
    // exp2 library in math.h.
    if (x.isnan()) {
      z_temp = z.nan();
      z_temp.set_signbit(x.signbit());
    }
    #endif
    
    z = z_temp;
  }

//=========================================================================
// Function: ac_exp2_cordic (for ac_ieee_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable exp2 function for
//    ac_ieee_float datatypes.
//
//    ac_exp2_cordic for ac_ieee_float types is a wrapper, the heavy lifting is
//    done by ac_exp_float_, ac_exp2_rr and all the other functions called
//    by ac_exp2_rr.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_ieee_float<binary32> input_type;
//    typedef ac_ieee_float<binary32> output_type;
//    
//    typedef ac_ieee_float<binary64> input_2_type;
//    typedef ac_ieee_float<binary64> output_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const input_type &input,
//      output_type &output,
//      const input_2_type &input_2,
//      output_2_type &output_2
//    )
//    {
//      ac_exp2_cordic(input, output);
//      // output_2's mantissa has a large number of fractional bits, so
//      // large that it would trigger a static_assert in the code if used
//      // as-is, on account of the default fractional width for temp.
//      // variables in ac_exp_float_. To avoid this, we override the
//      // default width by setting the "OR_TF" template parameter to
//      // "true" and the "TF_" parameter to 30.
//      // Refer to the documentation for more details.
//      ac_exp2_cordic<true, 30>(input_2, output_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      input_type input(2.5);
//      output_type output;
//      input_2_type input_2(3.5);
//      output_2_type output_2;
//      
//      CCS_DESIGN(project)(input, output, input_2, output_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  // The ac_ieee_float implementation for ac_exp2_cordic uses temporary ac_float variables which are in
  // turn passed to the ac_float implementation.
  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
            ac_ieee_float_format Format, ac_ieee_float_format outFormat>
  void ac_exp2_cordic(const ac_ieee_float<Format> &x, ac_ieee_float<outFormat> &z)
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int ZW = T_out::width;
    const int ZE = T_out::e_width;
    
    ac_float<ZW - ZE + 1, 2, ZE> z_ac_fl; // Equivalent ac_float representation for output.
    ac_exp2_cordic<OR_TF, TF_, TQ>(x.to_ac_float(), z_ac_fl); // Call ac_float version.
    T_out z_temp(z_ac_fl); // Convert output ac_float to ac_ieee_float.
    
    #ifdef AC_HCORDIC_NAN_SUPPORTED
    // If the input is +/-nan, the output will be +/-nan. This mirrors the behavior of the
    // exp2 library in math.h.
    if (x.isnan()) {
      z_temp = z.nan();
      z_temp.set_signbit(x.signbit());
    }
    #endif
    
    z = z_temp;
  }

//=========================================================================
// Function: ac_pow_cordic (for ac_fixed)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable pow function for
//    ac_fixed datatypes.
//
//    This implementation tries to be as general as possible without
//    taking into account the actual range of arguments 'a', and 'b'.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_fixed<11,  3, false, AC_RND, AC_SAT> a_type;
//    typedef ac_fixed<12,  4,  true, AC_RND, AC_SAT> b_type;
//    typedef ac_fixed<32, 16, false, AC_RND, AC_SAT> z_type;
//    
//    typedef ac_fixed<50,  5, false, AC_RND, AC_SAT> a_2_type;
//    typedef ac_fixed<32,  5,  true, AC_RND, AC_SAT> b_2_type;
//    typedef ac_fixed<65, 20, false, AC_RND, AC_SAT> z_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const a_type &a,
//      const b_type &b,
//      z_type &z,
//      const a_2_type &a_2,
//      const b_2_type &b_2,
//      z_2_type &z_2
//    )
//    {
//      ac_pow_cordic(a, b, z);
//      // The bitwidth of a_2, b_2 and z_2 are so large that using them
//      // would trigger static_asserts in the code on account of 
//      // the default width assignments for temp variables.
//      // To avoid this, we override the default widths by setting 
//      // the "OR_TF" template parameter to "true" and the "TF_"
//      // parameter to 30.
//      // Refer to the documentation for more details.
//      ac_pow_cordic<true, 30>(a_2, b_2, z_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      a_type a = 1.5;
//      b_type b = 2.5;
//      z_type z;
//      a_2_type a_2 = 3.5;
//      b_2_type b_2 = 4.5;
//      z_2_type z_2;
//      
//      CCS_DESIGN(project)(a, b, z, a_2, b_2, z_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------
  
  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
            int AW, int AI, ac_q_mode AQ, ac_o_mode AV,
            int BW, int BI, bool BS, ac_q_mode BQ, ac_o_mode BV,
            int ZW, int ZI, ac_q_mode ZQ, ac_o_mode ZV>
  void ac_pow_cordic(const ac_fixed<AW,AI,false,AQ,AV> &a,
                     const ac_fixed<BW,BI,BS,BQ,BV> &b,
                     ac_fixed<ZW,ZI,false,ZQ,ZV> &z)
  {
    ac_fixed<AW, 0, false, AQ, AV> a_norm;
    int a_exp = ac_normalize(a, a_norm);
    
    // By default, the number of fractional bits in a_norm_2 are set to the number of input bits.
    // The user can, however, override this default assignment by using the OR_TF and TF_ template parameters.
    const int a_norm_2_f_bits = OR_TF ? TF_ : AW;
    
    static_assert(a_norm_2_f_bits <= 44, "Number of fractional bits in a_norm_2 must be no more than 44. Consider overriding the default fractional bits allocated with the OR_TF and TF_ template parameters, so as to satisfy this condition. Refer to the documentation for more information.");
    
    ac_fixed<a_norm_2_f_bits, 0, false, TQ> a_norm_2 = a_norm;
    
    // By default, the number of fractional bits in log2_a are set to the number of output bits.
    // The user can, however, override this default assignment by using the OR_TF and TF_ template parameters.
    const int log2_a_F = OR_TF ? TF_ : ZW;
    
    static_assert(log2_a_F <= 44, "Number of fractional bits in log2_a must be no more than 44. Consider overriding the default fractional bits allocated with the OR_TF and TF_ template parameters, so as to satisfy this condition. Refer to the documentation for more information.");
    
    ac_fixed<log2_a_F + 2, 2, true> log2_a_norm;
    // Range reduced to 0.5 <= a_norm_2 < 1, -1 <= log2_a_norm < 1
    ac_log2_rr(a_norm_2, log2_a_norm);
    
    const int log2_a_I = AC_MAX(int(ac::nbits<AW-AI>::val), int(ac::nbits<AI>::val)) + 1;
    
    typedef ac_fixed<log2_a_F, log2_a_I, true> log2_a_type;
    log2_a_type log2_a = log2_a_norm + a_exp;
    
    // Find number of integer and fractional bits in product of b and log2_a.
    typedef ac_fixed<BW,BI,BS,BQ,BV> b_type;
    typedef typename ac::rt_2T<log2_a_type, b_type>::mult prod_type;
    const int prod_i_bits = prod_type::i_width;
    // Limit the number of frac. bits in the product to 45, avoiding a static_assert in ac_exp2_cordic.
    const int prod_f_bits = AC_MIN(prod_type::width - prod_i_bits, 45);
    
    // Multiply b by log2(a) and pass the product to the ac_exp2_cordic function.
    ac_fixed<prod_f_bits + prod_i_bits, prod_i_bits, true> prod = b*log2_a;
    ac_exp2_cordic<OR_TF, TF_, TQ>(prod, z);
  }

//=========================================================================
// Function: ac_pow_cordic (for ac_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable pow function for
//    ac_float datatypes.
//
//    This implementation tries to be as general as possible without
//    taking into account the actual range of arguments 'a', and 'b'.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//    
//    typedef ac_float<16, 2, 5, AC_RND> a_type;
//    typedef ac_float<16, 2, 6, AC_RND> b_type;
//    typedef ac_float<32, 2, 8, AC_RND> z_type;
//    
//    typedef ac_float<50, 2,  7, AC_RND> a_2_type;
//    typedef ac_float<51, 2,  7, AC_RND> b_2_type;
//    typedef ac_float<60, 2, 10, AC_RND> z_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const a_type &a,
//      const b_type &b,
//      z_type &z,
//      const a_2_type &a_2,
//      const b_2_type &b_2,
//      z_2_type &z_2
//    )
//    {
//      ac_pow_cordic(a, b, z);
//      // The bitwidth of a_2, b_2 and z_2 are so large that using them
//      // would trigger static_asserts in the code on account of 
//      // the default width assignments for temp variables.
//      // To avoid this, we override the default widths by setting 
//      // the "OR_TF" template parameter to "true" and the "TF_"
//      // parameter to 30.
//      // Refer to the documentation for more details.
//      ac_pow_cordic<true, 30>(a_2, b_2, z_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      a_type a = 1.5;
//      b_type b = 2.5;
//      z_type z;
//      a_2_type a_2 = 3.5;
//      b_2_type b_2 = 4.5;
//      z_2_type z_2;
//      
//      CCS_DESIGN(project)(a, b, z, a_2, b_2, z_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
            int AW, int AI, int AE, ac_q_mode AQ,
            int BW, int BI, int BE, ac_q_mode BQ,
            int ZW, int ZI, int ZE, ac_q_mode ZQ>
  void ac_pow_cordic(const ac_float<AW, AI, AE, AQ> &a,
                     const ac_float<BW, BI, BE, BQ> &b,
                     ac_float<ZW, ZI, ZE, ZQ> &z)
  {
    #ifdef ASSERT_ON_INVALID_INPUT
    AC_ASSERT(a >= 0, "Negative base input not supported.");
    #endif
    
    ac_fixed<AW - 1, AI - 1, false> a_mant = a.mantissa();
    ac_int<AE, true> a_exp = a.exp();

    ac_fixed<AW - 1, 0, false> a_norm;
    int a_norm_exp = ac_normalize(a_mant, a_norm);
    
    // Both the a_norm variable and the log2_a variable eventually have to be passed to ac_log2_rr,
    // and we might need to reduce the fractional variables in both to reduce the chances of the
    // compiler throwing an error in the shift_dist function. To do so, use the "OR_TF" and "TF_"
    // template parameters to override the default number of fractional bits assigned to both.
    
    const int a_norm_2_f_bits = OR_TF ? TF_ : AW - 1;
    
    static_assert(a_norm_2_f_bits <= 44, "Number of fractional bits in a_norm_2 must be no more than 44. Consider overriding the default fractional bits allocated with the OR_TF and TF_ template parameters, so as to satisfy this condition. Refer to the documentation for more information.");
    
    ac_fixed<a_norm_2_f_bits, 0, false, TQ> a_norm_2 = a_norm;
    
    const int log2_a_I = ac::nbits<AC_MAX(AI - 1, AW - AI - 1) + (1 << (AE - 1)) + 1>::val + 1;
    const int log2_a_F = OR_TF ? TF_ : ZW;
    
    static_assert(log2_a_F <= 44, "Number of fractional bits in log2_a must be no more than 44. Consider overriding the default fractional bits allocated with the OR_TF and TF_ template parameters, so as to satisfy this condition. Refer to the documentation for more information.");
    
    ac_fixed<log2_a_F + 2, 2, true> log2_a_norm;
    // Range reduced to 0.5 <= a_norm_2 < 1, -1 <= log2_a_norm < 1
    ac_log2_rr(a_norm_2, log2_a_norm);
    ac_fixed<log2_a_I, log2_a_I, true> offset = a_exp + a_norm_exp;
    ac_fixed<log2_a_I + log2_a_F, log2_a_I, true> log2_a = log2_a_norm + offset;
    
    const int z_exp_min = 1 << (ZE - 1);
    const int TI = ac::nbits<z_exp_min + AC_MAX(ZI - 1, ZW - ZI)>::val + 1;
    
    const int TF = OR_TF ? TF_ : ZW - ZI;
    
    // If TF < 1, it'll cause an issue with bit-slicing operations further down the line.
    // If TF > 44, (TF + 1) > 45, which will cause a compiler error further down the line, in the
    // shift_dist function.
    static_assert(TF >= 1 && TF <= 44, "TF must be a positive integer and lesser than 45. Consider overriding the default value with the OR_TF and TF_ template parameters in your function call, to satisfy these conditions. Refer to the documentation for more information.");
    
    ac_fixed<BW, BI, true> b_mant = b.mantissa();
    ac_int<BE, true> b_exp = b.exp();
    
    ac_fixed<TF + TI, TI, true, TQ, AC_SAT> input_to_exp2_fxpt;
    // Multiply the mantissa of b with log2(a) and then convert the product to an ac_fixed variable
    // through the following left-shift operation: mantissa << exp.
    // The precision of input_to_exp2_fxpt is such that it will saturate a little
    // above the input which would cause the floating point output to saturate if we were using an
    // accurate math function to calculate the exp2 output.
    ac_math::ac_shift_left(b_mant*log2_a, b_exp.to_int(), input_to_exp2_fxpt);
    
    ac_fixed<TF, 0, false> input_to_exp2_fxpt_frac;
    input_to_exp2_fxpt_frac.set_slc(0, input_to_exp2_fxpt.template slc<TF>(0));
    ac_fixed<TF + 1, 1, false, TQ> output_to_fxpt_frac;
    ac_exp2_rr(input_to_exp2_fxpt_frac, output_to_fxpt_frac);
    
    ac_int<TI, true> input_to_exp2_fxpt_int = input_to_exp2_fxpt.to_int();
    ac_float<ZW, ZI, ZE, ZQ> z_temp(output_to_fxpt_frac, input_to_exp2_fxpt_int, true);
    
    if (a == 1) {
      z = 1;
    } else if (a == 0) {
      z = (b == 0) ? 1 : 0;
    } else {
      z = z_temp;
    }
    
    #if defined(AC_HCORDIC_H_DEBUG) && !defined(__SYNTHESIS__)
    std::cout << "a = " << a << std::endl;
    std::cout << "log2_a = " << log2_a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "input_to_exp2_fxpt = " << input_to_exp2_fxpt << std::endl;
    std::cout << "z_temp = " << z_temp << std::endl;
    std::cout << "z = " << z << std::endl;
    #endif
  }

//=========================================================================
// Function: ac_pow_cordic (for ac_std_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable pow function for
//    ac_std_float datatypes.
//
//    This implementation tries to be as general as possible without
//    taking into account the actual range of arguments 'a', and 'b'.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//        
//    typedef ac_std_float<32, 8> a_type;
//    typedef ac_std_float<32, 8> b_type;
//    typedef ac_std_float<32, 8> z_type;
//        
//    typedef ac_std_float<64, 11> a_2_type;
//    typedef ac_std_float<64, 11> b_2_type;
//    typedef ac_std_float<64, 11> z_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const a_type &a,
//      const b_type &b,
//      z_type &z,
//      const a_2_type &a_2,
//      const b_2_type &b_2,
//      z_2_type &z_2
//    )
//    {
//      ac_pow_cordic(a, b, z);
//      // The bitwidth of a_2, b_2 and z_2 are so large that using them
//      // would trigger static_asserts in the code on account of 
//      // the default width assignments for temp variables.
//      // To avoid this, we override the default widths by setting 
//      // the "OR_TF" template parameter to "true" and the "TF_"
//      // parameter to 30.
//      // Refer to the documentation for more details.
//      ac_pow_cordic<true, 30>(a_2, b_2, z_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      a_type a(1.5);
//      b_type b(2.5);
//      z_type z;
//      a_2_type a_2(3.5);
//      b_2_type b_2(4.5);
//      z_2_type z_2;
//      
//      CCS_DESIGN(project)(a, b, z, a_2, b_2, z_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  // ac_std_float implementation of the ac_pow_cordic() function.
  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
            int AW, int AE, int BW, int BE, int ZW, int ZE>
  void ac_pow_cordic(const ac_std_float<AW, AE> &a,
                     const ac_std_float<BW, BE> &b,
                     ac_std_float<ZW, ZE> &z)
  {
    ac_float<ZW - ZE + 1, 2, ZE> z_ac_fl; // Equivalent ac_float representation for output.
    ac_pow_cordic<OR_TF, TF_, TQ>(a.to_ac_float(), b.to_ac_float(), z_ac_fl);
    ac_std_float<ZW, ZE> z_temp(z_ac_fl); // Convert output ac_float to ac_std_float.
    z = z_temp;
  }

//=========================================================================
// Function: ac_pow_cordic (for ac_ieee_float)
//
// Description:
//    Hyperbolic CORDIC implementation of synthesizable pow function for
//    ac_ieee_float datatypes.
//
//    This implementation tries to be as general as possible without
//    taking into account the actual range of arguments 'a', and 'b'.
//
// Usage:
//    A sample testbench and its implementation look like
//    this:
//
//    #include <ac_math/ac_hcordic.h>
//    using namespace ac_math;
//        
//    typedef ac_ieee_float<binary32> a_type;
//    typedef ac_ieee_float<binary32> b_type;
//    typedef ac_ieee_float<binary32> z_type;
//        
//    typedef ac_ieee_float<binary64> a_2_type;
//    typedef ac_ieee_float<binary64> b_2_type;
//    typedef ac_ieee_float<binary64> z_2_type;
//    
//    #pragma hls_design top
//    void project(
//      const a_type &a,
//      const b_type &b,
//      z_type &z,
//      const a_2_type &a_2,
//      const b_2_type &b_2,
//      z_2_type &z_2
//    )
//    {
//      ac_pow_cordic(a, b, z);
//      // The bitwidth of a_2, b_2 and z_2 are so large that using them
//      // would trigger static_asserts in the code on account of 
//      // the default width assignments for temp variables.
//      // To avoid this, we override the default widths by setting 
//      // the "OR_TF" template parameter to "true" and the "TF_"
//      // parameter to 30.
//      // Refer to the documentation for more details.
//      ac_pow_cordic<true, 30>(a_2, b_2, z_2);
//    }
//    
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      a_type a(1.5);
//      b_type b(2.5);
//      z_type z;
//      a_2_type a_2(3.5);
//      b_2_type b_2(4.5);
//      z_2_type z_2;
//      
//      CCS_DESIGN(project)(a, b, z, a_2, b_2, z_2);
//      
//      CCS_RETURN (0);
//    }
//    #endif
//
//-------------------------------------------------------------------------

  // ac_ieee_float implementation of the ac_pow_cordic() function.
  template <bool OR_TF = false, // Override default fractional bitwidth for temp variables?
            int TF_ = 32, // Template argument for fractional bitwidth to override with.
            ac_q_mode TQ = AC_TRN, // Rounding mode allocated for temporary variables.
            ac_ieee_float_format aFormat, ac_ieee_float_format bFormat, ac_ieee_float_format outFormat>
  void ac_pow_cordic(const ac_ieee_float<aFormat> &a,
                     const ac_ieee_float<bFormat> &b,
                     ac_ieee_float<outFormat> &z)
  {
    typedef ac_ieee_float<outFormat> T_out;
    const int ZW = T_out::width;
    const int ZE = T_out::e_width;
    ac_float<ZW - ZE + 1, 2, ZE> z_ac_fl; // Equivalent ac_float representation for output.
    ac_pow_cordic<OR_TF, TF_, TQ>(a.to_ac_float(), b.to_ac_float(), z_ac_fl);
    T_out z_temp(z_ac_fl); // Convert output ac_float to ac_ieee_float.
    z = z_temp;
  }
}

#endif // _INCLUDED_AC_HCORDIC_H_
