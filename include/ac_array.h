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
//  *************************************************************************
//  File : ac_array.h
//
//  Created on: Jun 14, 2017
//
//  Description:  Container class for multi-dimensional (up to 3-D) arrays.
//     Using this container class simplifies copy while still providing
//     the common C-style array index operator [].
//     The number of dimension sizes specified in the template parameterization
//     determines the number of dimensions of the ac_array object:
//        ac_array<int,4> x;       // a 1-D array of integers with 4 elements
//        ac_array<bool,5,6> y;    // a 2-D array of bools of size [5][6]
//        ac_array<short,3,4,5> z; // a 3-D array of size [3][4][5]
//
//     The C++ ostream left-shift operator provides pretty-print output:
//        cout << z << endl;
//     Comparison operators ==, != compare two arrays of the same size.
//     Array indexing uses []:
//        cout << "element z(2,1,1)=" << z[2][1][1] << endl;
//
//     Note: The base class 'ac_array_b' is not intended to be used
//           directly - only use ac_array.
//
//  Usage:
//    Below is a sample design and C++ testbench
//
//    #include <ac_int.h>
//    #include <ac_array.h>
//
//    typedef ac_int<8,false> dtype_t;
//
//    void design(ac_array<dtype_t,3,4> a_2Din, dtype_t scale,
//                ac_array<dtype_t,3,4> &a_2Dout)
//    {
//      for (int i=0; i<3; i++) {
//        for (int j=0; j<4; j++) {
//          a_2Dout[i][j] = a_2Din[i][j] * scale;
//        }
//      }
//    }
//
//    #ifndef __SYNTHESIS__
//    #include <mc_scverify.h>
//    CCS_MAIN(int argc, char *argv[])
//    {
//      ac_array<dtype_t,3,4> tb_in;
//      ac_array<dtype_t,3,4> tb_out;
//      dtype_t               scale;
//      for (int iteration=1; iteration<10; iteration++) {
//        // fill input 2-D array
//        for (int i=0; i<3; i++) for (int j=0; j<4; j++) { tb_in[i][j] = iteration*(i*4+j); }
//        scale = iteration;
//        // debug dump of contents
//        std::cout << "Array input:  " << tb_in << std::endl;
//        CCS_DESIGN(design)(tb_in, scale, tb_out);
//        // debug dump of output
//        std::cout << "Array output: " << tb_out << std::endl;
//      }
//      CCS_RETURN(0);
//    }
//    #endif
//
//  *************************************************************************

#ifndef _INCLUDED_AC_ARRAY_H_
#define _INCLUDED_AC_ARRAY_H_

#include <sstream>
#ifdef AC_ARRAY_ASSERT
#include <ac_assert.h>
#define AC_A_ASSERT(cond) ac::ac_assert(__FILE__, __LINE__, #cond, cond)
#else // !AC_ARRAY_ASSERT
#define AC_A_ASSERT(cond)
#endif // AC_ARRAY_ASSERT

#ifndef __SYNTHESIS__
#include <iostream>
#endif

// Forward declaration
template <typename T, unsigned D1, unsigned D2, unsigned D3> class ac_array;

//=========================================================================
// Class: ac_array_b
//
// Description: Base class for n-d array
//   This class is not directly used in designs - use ac_array<>.
//
//-------------------------------------------------------------------------


template <typename T, typename Td, unsigned D>
class ac_array_b
{
public:
  typedef T ElemType;

  ac_array_b() {}
  virtual ~ac_array_b() {}

  T &operator[] (unsigned i) {
    AC_A_ASSERT(i < D);
    return d[i];
  }

  const T &operator[] (unsigned i) const {
    AC_A_ASSERT(i < D);
    return d[i];
  }

  unsigned size() const { return D; }

  void set(const Td &ival) {
    for ( unsigned i = 0; i < D; ++i ) {
      set(i, ival);
    }
  }

  void set(unsigned i, const Td &ival) {
    AC_A_ASSERT(i < D);
    d[i] = ival;
  }

  void clearAll(bool dc=false) {
    for ( unsigned i = 0; i < D; ++i ) {
      clear(i, dc);
    }
  }

  void clear(unsigned i, bool dc=false) {
    AC_A_ASSERT(i < D);
    if ( !dc ) {
      d[i] = 0;
    } else {
      Td v;
      d[i] = v;
    }
  }

  // Relational operators - can only compare against
  // same sized arrays
  bool operator==(const ac_array_b<T,Td,D> &other) const {
    bool equal = true;
    for (int i=0; i<D; i++) {
      if (this->d[i] != other[i]) {
        equal = false;
        break;
      }
    }
    return equal;
  }

  bool operator!=(const ac_array_b<T,Td,D> &other) const {
    return !operator==(other);
  }

public: // data
  T d[D];
};

//=========================================================================
// Class: ac_array
//
//-------------------------------------------------------------------------

template <typename T, unsigned D1, unsigned D2=0, unsigned D3=0>
class ac_array : public ac_array_b< ac_array<T,D2,D3>, T, D1>
{
  typedef ac_array_b< ac_array<T,D2,D3>, T, D1> Base;
public:

  ac_array() {}
  ac_array(const T &ival) { Base::set(ival); }
  virtual ~ac_array() {}

  ac_array &operator= (const T &v) { Base::set(v); return *this; }
};

//=========================================================================
// Specialization Class: ac_array for 1 dimension (row)
//
//-------------------------------------------------------------------------

template <typename T, unsigned D1>
class ac_array<T,D1,0,0> : public ac_array_b<T,T,D1>
{
  typedef ac_array_b<T,T,D1> Base;
public:
  ac_array() {}
  ac_array(const T &ival) { Base::set(ival); }
  virtual ~ac_array() {}

  ac_array &operator= (const T &v) { Base::set(v); return *this; }
};

template <typename T>
class ac_array<T,0,0,0> : public ac_array_b<T,T,1>
{
  typedef ac_array_b<T,T,1> Base;
public:
  ac_array() {}
  ac_array(const T &ival) { Base::set(ival); }
  virtual ~ac_array() {}

  ac_array &operator= (const T &v) { Base::set(v); return *this; }
};

//=========================================================================
// Non-synthesis helper functions
//
//-------------------------------------------------------------------------
#ifndef __SYNTHESIS__

//=======================================================================
// Pretty-print with ostream operator <<
//-----------------------------------------------------------------------
template <typename T, typename Td, unsigned D>
std::ostream &operator<<(std::ostream &os, const ac_array_b<T,Td,D> &a)
{
  os << '{';
  for (int i=0; i<D; i++) { os << a[i]; if (i<D-1) os << ' '; }
  os << '}';
  return os;
}

#if defined(IEEE_1666_SYSTEMC)
template <typename T, typename Td, unsigned D>
inline void sc_trace(sc_trace_file *tf, const ac_array_b<T,Td,D> &v, const std::string &nm)
{
  for (unsigned i=0; i<D; i++) {
    std::ostringstream os;
    os << nm << "(" << i << ")";
    sc_trace(tf, v.d[i], os.str());
  }
}
#endif

#endif

#endif

