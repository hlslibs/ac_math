/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 2026.2                                              *
 *                                                                        *
 *  Release Date    : Tue Jun 30 15:00:22 PDT 2026                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 2026.2.1                                            *
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
//   Brick-Sort HLS implementation
//   -----------------------------
//
//   Template parameters
//   -------------------
//   -> DataType    = pixelType        Element data-type (8-bit unsigned in
//                                  	design_config.h).  Wider types increase the
//                                  	width of every comparator and swapper.
//
//   -> ArraySize   = NUM_ELEMENTS  	Number of items to sort.  Determines the
//                                  	loop bounds as well as the size of the local
//                                  	ac_array used inside the DUT.
//
//   -> Ascending   = true|false     	Sorting direction.  The comparison operator
//                                  	flips at compile time, so there is no run-
//                                  	time penalty for choosing descending order.
//
//   C-CORE pragmas inside the DUT
//   -----------------------------
//      #pragma hls_design ccore
//      #pragma ccore_type combinational
//
//   They wrap each of the two kernel functions:
//         void sortStageEven(...);
//         void sortStageOdd (...);
//
//   -> Runtime (compile / sim)
//   		Because Catapult emits one shared RTL module
//       instead of inlining the logic in every loop iteration, large
//       instantiations (e.g. ArraySize = 64) compile and simulate
//       noticeably faster.
//
//  -> Latency
//   		Because a C-CORE is treated as a fixed black-box, Catapult can't
//  		pipeline or retime logic *across* the core boundary.  The optimisation
//  		scope is therefore smaller, so overall latency tends to increase.
//
//
//  Example instantiation in a testbench
//  ------------------------------------
//     	ac_brick_sort<pixelType, NUM_ELEMENTS, Ascending> dut;
//       	dut.hardwareAcceleratedSort(reference_array, hardware_array);
// **************************************************************************

#ifndef _INCLUDED_AC_BRICK_SORT_H_
#define _INCLUDED_AC_BRICK_SORT_H_

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
  class ac_brick_sort : public ac_sort_base<DataType>
  {
  private:
    // Add using declarations to bring base class members into scope
    using ac_sort_base<DataType>::swapper;

    // 3) Even-stage function
    //    Compare pairs (0,1), (2,3), (4,5)...
    #pragma hls_design ccore
    #pragma ccore_type combinational
    void sortStageEven(ac_array<DataType, ArraySize> &array) {
      #pragma hls_unroll yes
      LOOP_SorterEvenStage: for (unsigned int j = 0; j < ArraySize - 1; j += 2) {
        if (Ascending) {
          if (array[j] > array[j + 1]) {
            swapper(array[j], array[j + 1]);
          }
        } else {
          if (array[j] < array[j + 1]) {
            swapper(array[j], array[j + 1]);
          }
        }
      }
    }

    // 4) Odd-stage function
    //    Compare pairs (1,2), (3,4), (5,6)...
    //    Also handles wrap-around if SIZE is even and j == SIZE - 1.
    #pragma hls_design ccore
    #pragma ccore_type combinational
    void sortStageOdd(ac_array<DataType, ArraySize> &array) {
      #pragma hls_unroll yes
      LOOP_SorterOddStage: for (unsigned int j = 1; j < ArraySize; j += 2) {
        // Wrap-around swap only if SIZE is even and j == SIZE - 1
        if ((ArraySize % 2 == 0) && (j == ArraySize - 1)) {
          if (Ascending) {
            if (array[ArraySize - 1] < array[0]) {
              swapper(array[ArraySize - 1], array[0]);
            }
          } else {
            if (array[ArraySize - 1] > array[0]) {
              swapper(array[ArraySize - 1], array[0]);
            }
          }
        } else {
          // Normal adjacent comparison
          if (Ascending) {
            if (array[j] > array[j + 1]) {
              swapper(array[j], array[j + 1]);
            }
          } else {
            if (array[j] < array[j + 1]) {
              swapper(array[j], array[j + 1]);
            }
          }
        }
      }
    }

  public:
    // Constructor with member initializer list
    ac_brick_sort() {
#ifndef __SYNTHESIS__
      ac_sort_base<DataType>::template assertArraySize<ArraySize>();
#endif
    }

    // Top-level method
    void hardwareAcceleratedSort(
      ac_array<DataType, ArraySize> &in_array,
      ac_array<DataType, ArraySize> &out_array) {
      ac_array<DataType, ArraySize> vec;

      // 1) Copy input into local array
      vec = in_array;

      // 2) Determine iteration count
      // If ArraySize is odd => do ArraySize passes;
      // if ArraySize is even => do ArraySize - 1 passes.
      const uint32_t iterationCount = (ArraySize % 2 == 1) ? ArraySize : (ArraySize - 1);

      // 3) Perform brick sort by alternating calls to odd/even stages
      #pragma hls_unroll yes
      LOOP_Sort: for (unsigned int pass = 0; pass < iterationCount; pass++) {
        bool oddPass = (pass & 1) != 0;
        if (oddPass) {
          sortStageOdd(vec);
        } else {
          sortStageEven(vec);
        }
      }

      // 4) Copy sorted results to out_array
      out_array = vec;

#ifndef __SYNTHESIS__
      // 5) Simulation-only check final sorted order
      ac_sort_base<DataType>::template assertSortedOrder<ArraySize, Ascending>(out_array);
#endif
    }
  };

}

#endif

