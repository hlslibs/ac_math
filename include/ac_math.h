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
// Revision History:
//    3.1.0  - Added ac_tanh_pwl.h
//    2.0.10 - Removed including of ac_random.h from ac_math.h header file.
//             Completed list of header file inclusions in ac_math.h

#ifndef _INCLUDED_AC_MATH_H_
#define _INCLUDED_AC_MATH_H_

#include <ac_math/ac_abs.h>
// ac_abs()

#include <ac_math/ac_arccos_cordic.h>
// ac_arccos_cordic()

#include <ac_math/ac_arcsin_cordic.h>
// ac_arcsin_cordic()

#include <ac_math/ac_atan_pwl.h>
// ac_atan_pwl()

#include <ac_math/ac_atan_pwl_ha.h>
// ac_atan_pwl_ha()

#include <ac_math/ac_atan2_cordic.h>
// ac_atan2_cordic()

#include <ac_math/ac_barrel_shift.h>
//ac_barrel_shift()

#include <ac_math/ac_div.h>
// ac_div()

#include <ac_math/ac_hcordic.h>
// ac_log_cordic()
// ac_log2_cordic()
// ac_exp_cordic()
// ac_exp2_cordic()
// ac_pow_cordic()

#include <ac_math/ac_inverse_sqrt_pwl.h>
// ac_inverse_sqrt_pwl()

#include <ac_math/ac_log_pwl.h>
// ac_log_pwl()
// ac_log2_pwl()

#include <ac_math/ac_normalize.h>
// ac_normalize()

#include <ac_math/ac_pow_pwl.h>
// ac_pow2_pwl()
// ac_exp_pwl()

#include <ac_math/ac_reciprocal_pwl.h>
// ac_reciprocal_pwl()

#include <ac_math/ac_reciprocal_pwl_ha.h>
// ac_reciprocal_pwl_ha()

#include <ac_math/ac_shift.h>
// ac_shift_left()
// ac_shift_right()

#include <ac_math/ac_sigmoid_pwl.h>
// ac_sigmoid_pwl()

#include <ac_math/ac_sincos_cordic.h>
// ac_sincos_cordic()

#include <ac_math/ac_sincos_lut.h>
// ac_sincos_lut()

#include <ac_math/ac_sqrt.h>
// ac_sqrt()

#include <ac_math/ac_sqrt_pwl.h>
// ac_sqrt_pwl()

#include <ac_math/ac_tan_pwl.h>
// ac_tan_pwl()

#include <ac_math/ac_tanh_pwl.h>
// ac_tanh_pwl()

#include <ac_math/ac_softmax_pwl.h>
// ac_softmax_pwl()

#include <ac_float_add_tree.h>
// add_tree()
// add_tree_ptr()
// exp_equalizer()
// block_add_tree_ptr()

#endif

