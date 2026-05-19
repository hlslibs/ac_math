/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 2026.2                                              *
 *                                                                        *
 *  Release Date    : Tue May 12 21:06:11 PDT 2026                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 2026.2.0                                            *
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
//  Parallel Search HLS Implementation
//  ---------------------------------
// 
//  Key Components:
//  1. `ac_parallel_search`: 
//     - For arrays with size that is multiple of 2
//     - Uses binary tree reduction for parallel comparisons
//     - Fully unrolled comparison stages for high throughput
// 
//  Template Parameters
//  ------------------
//  -> DataType    = char            Element data-type. Can be any comparable
//                                   type. Affects width of comparators.
// 
//  -> ArraySize   = 16              Number of elements to search. Determines
//                                   the number of comparison stages and local
//                                   array size inside the DUT.
// 
//  -> FindMax     = true|false      Search direction. The comparison operator
//                                   is determined at compile time, so no run-
//                                   time penalty for finding minimum.
// 
//  Example Instantiation in a Testbench
//  -----------------------------------
//  // Parallel search for maximum value
//  		ac_parallel_search<int, 32, true> searcher;
//  		searcher.hardwareAcceleratedSearch(input_array, max_value);
// 
// *************************************************************************

#ifndef _INCLUDED_AC_SEARCH_PARALLEL_H_
#define _INCLUDED_AC_SEARCH_PARALLEL_H_ 

#include <ac_math/ac_search_base.h>

#if (defined(__GNUC__) && (__cplusplus < 201103L))
#error Please use C++11 or a later standard for compilation.
#endif
#if (defined(_MSC_VER) && (_MSC_VER < 1920) && !defined(__EDG__))
#error Please use Microsoft VS 2019 or a later standard for compilation.
#endif

// Helper function to calculate log base 2 at compile time
constexpr int calculateLog2(int n) {
    return (n <= 1) ? 0 : 1 + calculateLog2(n/2);
}

template <typename DataType = char, uint32_t ArraySize = 16, bool FindMax = true>
class ac_parallel_search : public ac_search_base<DataType> {
private:
    // Add using declarations to bring base class members into scope
    using ac_search_base<DataType>::compareMax;
    using ac_search_base<DataType>::compareMin;
    
public:
    // Constructor with member initializer list
    ac_parallel_search() {
    #ifndef __SYNTHESIS__
        ac_search_base<DataType>::template assertArraySize<ArraySize>();
    #endif
    } 

    // Top-level function 
	 void hardwareAcceleratedSearch(
        ac_array<DataType, ArraySize> &in_array,
        DataType &op_value) {

        ac_array<DataType, ArraySize> temp;
        
        // Copy input to local memory
        temp = in_array;
        
        // Number of active elements in current iteration
        uint32_t activeElements = ArraySize;
        
        // Main computation loops
        #pragma hls_unroll yes 
        Loop_ComputePairs: for (unsigned int step = 0; step < calculateLog2(ArraySize); step++) {
            // Number of pairs to compare in this iteration
            int numPairs = activeElements / 2;
            
            // Perform parallel comparisons
            #pragma hls_unroll yes
            Loop_ParallelCompare: for (unsigned int i = 0; i < numPairs; i++) {
                uint32_t baseIndex = i * 2;

                if (FindMax) {
                    temp[i] = compareMax(temp[baseIndex], temp[baseIndex + 1]);
                } else {
                    temp[i] = compareMin(temp[baseIndex], temp[baseIndex + 1]);
                }
            }
            
            // Handle odd element if present
            if (activeElements & 1) {
                temp[numPairs] = temp[activeElements - 1];
                activeElements = numPairs + 1;
            } else {
                activeElements = numPairs;
            }
        }
        
        // Assign final result
        op_value = temp[0];

        #ifndef __SYNTHESIS__
        if (FindMax) {
            // Verify the max value
            ac_search_base<DataType>::template assertMaxValue<ArraySize>(in_array, op_value);
        } else {
            // Verify the min value
            ac_search_base<DataType>::template assertMinValue<ArraySize>(in_array, op_value);
        }
        #endif
    }
};

#endif // _INCLUDED_AC_SEARCH_PARALLEL_H_
