/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 2025.4                                              *
 *                                                                        *
 *  Release Date    : Thu Dec 11 10:23:15 PST 2025                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 2025.4.1                                            *
 *                                                                        *
 *  Copyright 2025 Siemens                                                *
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
// **************************************************************************
//   Bubble-Sort instantiation
//   ------------------------
//
//   Template parameters
//   -------------------
//   -> DataType    = pixelType       	Element data-type (8-bit unsigned in
//                                 		design_config.h).  Wider types increase the
//                                  	width of every comparator and swapper.
//
//   -> ArraySize   = NUM_ELEMENTS   	Number of items to sort.  Determines the
//                                  	loop bounds as well as the size of the local
//                                  	ac_array used inside the DUT.
//
//   -> Ascending   = true|false     	Sorting direction.  The comparison operator
//                                  	flips at compile time, so there is no run-
//                                  	time penalty for choosing descending order.
//
//
//  Example instantiation in a testbench
//  ------------------------------------
//     	ac_bubble_sort<pixelType, NUM_ELEMENTS, Ascending> dut;
//       	dut.hardwareAcceleratedSort(reference_array, hardware_array);
// **************************************************************************

#ifndef _INCLUDED_AC_BUBBLE_SORT_H_
#define _INCLUDED_AC_BUBBLE_SORT_H_

#include <ac_math/ac_sort_base.h>

#if (defined(__GNUC__) && (__cplusplus < 201103L))
#error Please use C++11 or a later standard for compilation.
#endif
#if (defined(_MSC_VER) && (_MSC_VER < 1920) && !defined(__EDG__))
#error Please use Microsoft VS 2019 or a later standard for compilation.
#endif

namespace ac_math
{

  template <typename DataType = char, uint32_t ArraySize = 16, bool Ascending = true>
  class ac_bubble_sort : public ac_sort_base<DataType>
  {
  private:
    // Add using declarations to bring base class members into scope
    using ac_sort_base<DataType>::swapper;

  public:
    ac_bubble_sort() {
#ifndef __SYNTHESIS__
      ac_sort_base<DataType>::template assertArraySize<ArraySize>();
#endif
    }

    // Top-level method
    void hardwareAcceleratedSort(
      ac_array<DataType, ArraySize> &in_array,
      ac_array<DataType, ArraySize> &out_array) {

      ac_array<DataType, ArraySize> vec;

      // 1) Copy input into a local array
      vec = in_array;

      // 2) Perform Bubble Sort
      //    - If Ascending == true, smallest elements "bubble" to the front
      //    - If Ascending == false, largest elements "bubble" to the front
      #pragma hls_unroll yes
      LOOP_BubbleExtern: for (int i = ArraySize - 2; i >= 0; i--) {
        #pragma hls_unroll yes
        #pragma hls_waive CCC
        LOOP_BubbleInt: for (int j = 0; j < ArraySize - 1; j++) {
          // Compare condition changes based on Ascending
          #pragma hls_waive CNS
          if (Ascending) {
            // Ascending: Swap if out of order
            if (vec[j] > vec[j + 1]) {
              swapper(vec[j], vec[j + 1]);
            }
          } else {
            // Descending: Swap if out of order
            if (vec[j] < vec[j + 1]) {
              swapper(vec[j], vec[j + 1]);
            }
          }

          // Break early after the i-th pass
          if (j == i) break;
        }
      }

      // 3) Copy sorted data to output
      out_array = vec;

#ifndef __SYNTHESIS__
      // 4) Simulation-only check: Verify the final array is sorted
      ac_sort_base<DataType>::template assertSortedOrder<ArraySize, Ascending>(out_array);
#endif
    }
  };

}

#endif

