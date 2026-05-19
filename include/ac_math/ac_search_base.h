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
#ifndef _INCLUDED_AC_SEARCH_BASE_H_
#define _INCLUDED_AC_SEARCH_BASE_H_

#include <cstdint>
#include <ac_int.h>
#include <ac_array.h>
#include <mc_scverify.h>
#include <ac_assert.h>

#if (defined(__GNUC__) && (__cplusplus < 201103L))
#error Please use C++11 or a later standard for compilation.
#endif
#if (defined(_MSC_VER) && (_MSC_VER < 1920) && !defined(__EDG__))
#error Please use Microsoft VS 2019 or a later standard for compilation.
#endif

template <typename DataType>
class ac_search_base {
public:
	// Compare and return maximum of two elements
	inline DataType compareMax(DataType a, DataType b) {
		return (a > b) ? a : b;
	}
	
	// Compare and return minimum of two elements
	inline DataType compareMin(DataType a, DataType b) {
		return (a < b) ? a : b;
	}
	
	// Assert valid array size
	template <uint32_t ArraySize>
	static void assertArraySize() {
		AC_ASSERT(ArraySize > 0, "Array size must be positive");
	}
	
	// Assert maximum value is correct
	template <uint32_t ArraySize>
	static void assertMaxValue(ac_array<DataType, ArraySize> &arr, DataType maxValue) {
		for (unsigned int i = 0; i < ArraySize; i++) {
			AC_ASSERT(arr[i] <= maxValue, "Found element larger than reported maximum");
		}
		bool foundMax = false;
		for (unsigned int i = 0; i < ArraySize; i++) {
			if (arr[i] == maxValue) {
		   	foundMax = true;
		  		break;
			}
		}
		AC_ASSERT(foundMax, "Maximum value not found in array");
	}
	
	// Assert minimum value is correct
	template <uint32_t ArraySize>
	static void assertMinValue(ac_array<DataType, ArraySize> &arr, DataType minValue) {
		for (unsigned int i = 0; i < ArraySize; i++) {
			AC_ASSERT(arr[i] >= minValue, "Found element smaller than reported minimum");
		}
		bool foundMin = false;
		for (unsigned int i = 0; i < ArraySize; i++) {
			if (arr[i] == minValue) {
		   	foundMin = true;
		   	break;
			}
		}
		AC_ASSERT(foundMin, "Minimum value not found in array");
	}
};


#endif // _INCLUDED_AC_SEARCH_BASE_H_ 

