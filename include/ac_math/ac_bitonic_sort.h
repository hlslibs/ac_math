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
// *************************************************************************
//  Bitonic Sorter HLS implementation
//  ---------------------------------
//
//  Key Components:
//  1. `ac_bitonic_sort`:
//     - For fixed array sizes that are a power of 2.
//     - Uses `swap_core` to compare and conditionally swap elements.
//     - Fully unrolled sorting stages for high-throughput pipelines.
//
//   Template parameters
//   -------------------
//   -> DataType    = pixelType        Element data-type (8-bit unsigned in
//                                  	design_config.h).  Wider types increase the
//                                  	width of every comparator and swapper.
//
//   -> ArraySize   = NUM_ELEMENTS   	Number of items to sort.  Determines the
//                                  	loop bounds as well as the size of the local
//                                  	ac_array used inside the DUT.
//
//   -> Ascending   = true|false    	Sorting direction.  The comparison operator
//                                  	flips at compile time, so there is no run-
//                                  	time penalty for choosing descending order.
//
//   -> CCore		  = true|false			Generates CCORE for conditional swapping logic
//   												if set to true, else no CCORE is generated.
//
//
//   *  C-CORE pragmas inside the DUT
//   -----------------------------
//      #pragma hls_design ccore
//      #pragma CCoretype combinational
//
//   They wrap swap functions of swap_core class :
//         void conditionalSwap(...);
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
// 	 Example instantiation in a testbench
//   ------------------------------------
//   	// Bitonic sort (power of 2)
//    	ac_bitonic_sort<pixelType, NUM_ELEMENTS, Ascending, CCore> dut;
//    	dut.hardwareAcceleratedSort(reference_array, hardware_array);
//
// *************************************************************************

#ifndef _INCLUDED_AC_BITONIC_SORT_H_
#define _INCLUDED_AC_BITONIC_SORT_H_

#include <ac_math/ac_sort_base.h>

#if (defined(__GNUC__) && (__cplusplus < 201103L))
#error Please use C++11 or a later standard for compilation.
#endif
#if (defined(_MSC_VER) && (_MSC_VER < 1920) && !defined(__EDG__))
#error Please use Microsoft VS 2019 or a later standard for compilation.
#endif

namespace ac_math
{

  // Swap function for bitonic sorter w/ #elements = 2^k
  template<typename DataType>
  class swap_core : public ac_sort_base<DataType>
  {
  public:
    // Add using declarations to bring base class members into scope
    using ac_sort_base<DataType>::swapper;

    void conditionalSwap(DataType *valA, DataType *valB, bool sortDirection) {
      // If we want an overall Ascending result:
      //    isAscending -> swap if left > right
      //    !isAscending -> swap if left < right
      //
      // If we want an overall descending result, invert that logic.
      bool cond = (sortDirection && (*valA > *valB))
                  || (!sortDirection && (*valA < *valB));

      // Perform conditional swapping of values
      if (cond) {
        swapper(*valA, *valB);
      }
    }
  };

  // Swap function for bitonic sorter w/ #elements = 2^k (w/ ccore)
  template <typename DataType>
  class swap_core_ccore : public ac_sort_base<DataType>
  {
  public:
    // Add using declarations to bring base class members into scope
    using ac_sort_base<DataType>::swapper;

    // Wrapper function to include/ exclude pragma usage for CCORE generation
    #pragma hls_design ccore
    #pragma CCoretype combinational
    void conditionalSwapCCore(DataType *valA, DataType *valB, bool sortDirection) {
      // If we want an overall Ascending result:
      //    isAscending -> swap if left > right
      //    !isAscending -> swap if left < right
      //
      // If we want an overall descending result, invert that logic.
      bool cond = (sortDirection && (*valA > *valB))
                  || (!sortDirection && (*valA < *valB));

      // Perform conditional swapping of values
      if (cond) {
        swapper(*valA, *valB);
      }
    }
  };

  // Bitonic sort for power of 2 array size
  // Need to add assertion for cases where ArraySize != 2^n
  template <typename DataType = char, uint32_t ArraySize = 32, bool Ascending = true, bool CCore = true>
  class ac_bitonic_sort : public ac_sort_base<DataType>
  {
  public:
    // Number of stages for bitonic sort
    static constexpr uint32_t NREGS = ac::log2_ceil<ArraySize>::val + 1;

    // Create an instance of swap_core
    swap_core<DataType> swap_block;
    swap_core_ccore<DataType> swap_block_ccore;

    void counterLoop(
      ac_array<DataType, ArraySize> &vec,
      ac_int<NREGS, false> shifterValue1,
      ac_int<NREGS, false> shifterValue2
    ) {
      #pragma hls_unroll yes
      LOOP_Counter: for (unsigned int cntrIndex = 0; cntrIndex < ArraySize; cntrIndex++) {
        ac_int<NREGS, false> partnerIndex = (cntrIndex ^ shifterValue2);

#ifndef __SYNTHESIS__
        // Bounds checks in simulation
        AC_ASSERT(cntrIndex < ArraySize, "Index out of range!");
        AC_ASSERT(partnerIndex < ArraySize, "Index out of range!");
#endif

        if (partnerIndex > cntrIndex) {
          // Sub-block direction from bitonic
          bool isAscending = ((shifterValue1 & cntrIndex) == 0);

          // If we want an overall Ascending result:
          //    isAscending -> swap if left > right
          //    !isAscending -> swap if left < right
          //
          // If we want an overall descending result, invert that logic.
          bool sortDir = (Ascending ? isAscending : !isAscending);

          // Use the swap_core class object
          if (CCore) {
            swap_block_ccore.conditionalSwapCCore(&vec[cntrIndex], &vec[partnerIndex], sortDir);
          } else {
            swap_block.conditionalSwap(&vec[cntrIndex], &vec[partnerIndex], sortDir);
          }
        }
      }
    }

  public:
    // Member variables for bitonic sort
    ac_int<NREGS, false> shifterValue1 = 0;
    ac_int<NREGS, false> shifterValue2 = 0;

    // Constructor with member initializer list
    ac_bitonic_sort() {
#ifndef __SYNTHESIS__
      ac_sort_base<DataType>::template assertArraySize<ArraySize>();
#endif
    }

    // Top-level method
    void hardwareAcceleratedSort(
      ac_array<DataType, ArraySize> &in_array,
      ac_array<DataType, ArraySize> &out_array
    ) {
      // Local array for in-place sorting
      ac_array<DataType, ArraySize> vec;

      // 1) Copy in_array to a local array 'vec'
      vec = in_array;

      // 2) Nested loops for bitonic compare-and_swap
      #pragma hls_unroll yes
      LOOP_ShiftReg1: for (unsigned int shiftReg1Iter = 1; shiftReg1Iter < NREGS; shiftReg1Iter++) {
        shifterValue1 = (1 << shiftReg1Iter);

        #pragma hls_unroll yes
        LOOP_ShiftReg2: for (unsigned int shiftReg2Iter = 1; shiftReg2Iter < NREGS; shiftReg2Iter++) {
          if (shiftReg2Iter <= shiftReg1Iter) {
            ac_int<NREGS, false> local_shifterValue2 = (1 << (shiftReg1Iter - shiftReg2Iter));
            shifterValue2 = local_shifterValue2;

            counterLoop(vec, shifterValue1, local_shifterValue2);
          }
        }
      }

      // 3) Copy sorted data into out_array
      out_array = vec;

#ifndef __SYNTHESIS__
      // Final assertion check
      ac_sort_base<DataType>::template assertSortedOrder<ArraySize, Ascending>(out_array);
#endif
    }
  };

}

#endif

