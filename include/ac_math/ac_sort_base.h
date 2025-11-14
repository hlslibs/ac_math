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

#ifndef _INCLUDED_AC_SORT_BASE_H_
#define _INCLUDED_AC_SORT_BASE_H_

#include <ac_int.h>
#include <ac_array.h>

#ifndef __SYNTHESIS__
#include <ac_assert.h>
#endif

#if (defined(__GNUC__) && (__cplusplus < 201103L))
#error Please use C++11 or a later standard for compilation.
#endif
#if (defined(_MSC_VER) && (_MSC_VER < 1920) && !defined(__EDG__))
#error Please use Microsoft VS 2019 or a later standard for compilation.
#endif

template <typename DataType>
class ac_sort_base
{
protected:
  // Swap block
  template <typename T>
  void swapper(T &a, T &b) {
    T temp = a;
    a = b;
    b = temp;
  }

#ifndef __SYNTHESIS__
  // Assertion to ensure we have a sensible array size.
  template <uint32_t ArraySize>
  static void assertArraySize() {
    AC_ASSERT(ArraySize >= 2,
              "Sorting algorithm requires ArraySize >= 2 for meaningful operation.");
  }

  // Assertion to check the final sorted order
  template <uint32_t ArraySize, bool Ascending>
  static void assertSortedOrder(const ac_array<DataType, ArraySize> &out_array) {
    LOOP_SorterCheck: for (unsigned int i = 0; i < ArraySize - 1; i++) {
      if (Ascending) {
        AC_ASSERT(out_array[i] <= out_array[i + 1],
                  "Output is not sorted in Ascending order!");
      } else {
        AC_ASSERT(out_array[i] >= out_array[i + 1],
                  "Output is not sorted in descending order!");
      }
    }
  }
#endif

};

#endif

