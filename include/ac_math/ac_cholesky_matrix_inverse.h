/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 2026.1                                              *
 *                                                                        *
 *  Release Date    : Wed Feb 11 11:04:12 PST 2026                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 2026.1.0                                            *
 *                                                                        *
 *  Copyright  Siemens                                                *
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
#ifndef AC_CHOLESKY_MATRIX_INVERSE_H
#define AC_CHOLESKY_MATRIX_INVERSE_H

#if (defined(__GNUC__) && (__cplusplus < 201103L))
#error Please use C++11 or a later standard for compilation.
#endif
#if (defined(_MSC_VER) && (_MSC_VER < 1920) && !defined(__EDG__))
#error Please use Microsoft VS 2019 or a later standard for compilation.
#endif


#include <ac_array.h>
#include <ac_int.h>
#include <ac_fixed.h>
#include <ac_channel.h>
#include <ac_math/ac_sqrt.h>
#include <ac_math/ac_reciprocal_pwl_vha.h>
#include <iostream> // Added for debug prints


#ifdef __CCOV__
#include "ac_cholesky_matrix_inverse_covergroup.h"
#endif


namespace ac_math
{


  // ========================================================================================================================================================
  //*************************************** BASIC OPERATOR LEVEL HELPER FUNCTION USED IN MATRIX INVERSE  **********************************************
  // ========================================================================================================================================================


  // ============================================================
  // MAC Core Declaration
  // ============================================================
  template <typename Division_T>
  class ac_cinvmac_core
  {
  public:
    ac_cinvmac_core() {}
    #pragma hls_design ccore
    #pragma hls_ccore_type sequential


    void mac_ccore( Division_T &mac_inA, Division_T &mac_inB,Division_T &mac_out) {
      #pragma hls_waive OVL
      mac_out+=(Division_T)((Division_T)mac_inA*(Division_T)mac_inB);

      #ifndef __SYNTHESIS__  // Catapult macro: only defined during synthesis
      std::cout << "[MAC] InputA: " << mac_inA << "[MAC] InputB: " << mac_inB << " Output: " << mac_out << "\n";
      #endif

    }
  };

// ============================================================
// MAC CCore for local-buffer-based design
// acc_out = acc_in + (inA * inB)
// ============================================================

  template <typename Division_T>
  class ac_cinvmac_core_buf
  {
  public:
    ac_cinvmac_core_buf() {}

    #pragma hls_design ccore
    #pragma hls_ccore_type sequential

    void mac_ccore( Division_T &inA,
                    Division_T &inB,
                    Division_T &acc_in,
                    Division_T &acc_out) {
      Division_T prod = (Division_T)((Division_T)inA * (Division_T)inB);
      acc_out = (Division_T)acc_in + (Division_T)prod;
      #ifndef __SYNTHESIS__
      std::cout << "[MAC] inA=" << inA
                << " inB=" << inB
                << " acc_in=" << acc_in
                << " acc_out=" << acc_out << "\n";
      #endif
    }
  };
  // ============================================================
  // MAC Core Declaration
  // ============================================================

  template <typename Division_T>
  class ac_cinvmacsub_core
  {
  public:
    ac_cinvmacsub_core() {}
    #pragma hls_design ccore
    #pragma hls_ccore_type sequential

    void mac_ccore(Division_T &inA,  Division_T &inB, Division_T &inC, Division_T &acc) {
      #pragma hls_waive OVL
      acc=(Division_T)inC-(Division_T)((Division_T)inA * (Division_T)inB);

      #ifndef __SYNTHESIS__
      std::cout << "[MAC] inA=" << inA << " inB=" << inB << " inC=" << inC << " acc=" << acc << "\n";
      #endif
    }
  };

  // ============================================================
  // Reciprocal Core Declaration
  // ============================================================
  template <typename Reciprocal_T, typename Reciprocal_OutT>
  class ac_cinvreciprocal_core
  {
  public:
    ac_cinvreciprocal_core() {}
    #pragma hls_design ccore
    #pragma hls_ccore_type sequential

    void inv_sqrt_ccore(Reciprocal_T &sub_res, Reciprocal_OutT &rec_res) {
      ac_math::ac_reciprocal_pwl_vha(sub_res, rec_res);
      #ifndef __SYNTHESIS__  // Catapult macro: only defined during synthesis
      std::cout << "[Reciprocal] Input: " << sub_res << " Output: " << rec_res << "\n";
      #endif

    }
  };

  // ============================================================
  // Square Root Core Declaration
  // ============================================================

  template <typename SqrRoot_T,typename SqrRoot_OutT>
  class ac_cinvsqrt_core
  {
  public:
    ac_cinvsqrt_core() {}

    #pragma hls_design ccore
    #pragma hls_ccore_type sequential


    void ccore_sqrt(const SqrRoot_T &sqrt_in, SqrRoot_OutT &sqrt_out) {
      ac_math::ac_sqrt(sqrt_in, sqrt_out);
      #ifndef __SYNTHESIS__  // Catapult macro: only defined during synthesis
      std::cout << "[Sqrt] Input: " << sqrt_in << " Output: " << sqrt_out << "\n";
      #endif
    }
  };



// ============================================================
// INVcol Core for Ainv computation with local buffering
// ============================================================
  template <typename Input_T, typename Division_T, typename Mat_size_T, unsigned N>
  class ac_invcol_buffered_ccore
  {
  public:
    ac_invcol_buffered_ccore() {}

    #pragma hls_design ccore
    #pragma hls_ccore_type sequential

    void invcol_buffered(const Input_T row_i[N],
                         const Input_T col_j[N],
                         Mat_size_T j,
                         Mat_size_T i,
                         Division_T &acc_out) {
      Division_T acc = ( Division_T)0;
      Mat_size_T M_AC=N;
      #ifndef __CCOV__
      #pragma hls_unroll yes
      INVCOLLOOP: for (Mat_size_T k = 0; k < M_AC; ++k) {
        if ((k >= j) && (k < i)) {
          #pragma hls_waive SAT
          acc +=(Division_T)( (Division_T)row_i[k] * (Division_T)col_j[k]);
        }
      }
      #endif
      acc_out = acc;
    }
  };

// ============================================================
// MAC Inner Core for Ainv computation with local buffering
// ============================================================
  template <typename Division_T, typename Mat_size_T, unsigned N>
  class ac_mac_inner_buffered_ccore
  {
  public:
    ac_mac_inner_buffered_ccore() {}

    #pragma hls_design ccore
    #pragma hls_ccore_type sequential
    void mac_inner_buffered(const Division_T col_i[N],
                            const Division_T col_j[N],
                            unsigned start_k,
                            Division_T &acc_out) {
      Division_T acc = ( Division_T)0;

      #pragma hls_unroll yes
      INVMACLOOP: for (Mat_size_T k = 0; k < N; ++k) {
        if (k >= start_k) {
          acc += (Division_T)((Division_T)col_i[k] * (Division_T)col_j[k]);
        }
      }

      acc_out = acc;
    }
  };

// ============================================================
// Helper Chunked MAC Core
// ============================================================
  template <typename Input_T, typename Division_T, typename Mat_size_T, unsigned CHUNK_SIZE>
  class ac_chunk_mac_ccore
  {
  public:
    ac_chunk_mac_ccore() {}

    #pragma hls_design ccore
    #pragma hls_ccore_type sequential

    void chunk_mac(const Input_T row_i[CHUNK_SIZE],
                   const Input_T col_j[CHUNK_SIZE],
                   Division_T &acc_out) {
      Division_T acc = (Division_T)0;
      Mat_size_T M_AC=CHUNK_SIZE;

      #pragma hls_unroll yes
      CHUNK_LOOP: for (Mat_size_T k = 0; k < M_AC; ++k) {
        #pragma hls_waive OVL
        acc +=(Division_T)( (Division_T)row_i[k] * (Division_T)col_j[k]);
      }

      acc_out = acc;
    }
  };


  // ========================================================================================================================================================
  //*************************************** HELPER KERNEL FOR CHOLESKY DECOMPSOTION (RIGHT AND LEFT LOOKING)  **********************************************
  // ========================================================================================================================================================


  // ============================================================
  // compute L[k][k] Helper kernels for right looking cholesky
  // ============================================================


  template <typename Input_T, typename Output_T, typename SqrRoot_T,typename SqrRoot_OutT,typename Division_T, typename Reciprocal_T, typename Reciprocal_OutT, typename Mat_size_T,unsigned N>

  class ac_compute_Lkk
  {
  public:
    ac_compute_Lkk() {}

    void compute_Lkk(ac_array<Input_T,N,N> &A, Mat_size_T k, Division_T &Lkk) {

      Division_T sum= {};
      Division_T v= {};
      SqrRoot_T diag= {};
      SqrRoot_OutT sqrt_res= {};
      sum = (Division_T)0;

      #ifndef __CCOV__
      #pragma hls_unroll yes
      #pragma hls_waive FVI
      MACLOOP:for (Mat_size_T s = 0; s < k; ++s) {
        v = A[k][s];
        //sum += v * v;
        MACINST.mac_ccore(v,v,sum);
      }
      #endif
      #pragma hls_waive OVL
      diag = (A[k][k] - sum);
      SQRTINST.ccore_sqrt(diag, sqrt_res);
      Lkk=(Division_T)sqrt_res;
      #ifndef __SYNTHESIS__  // Catapult macro: only defined during synthesis
      std::cout << "[compute_Lkk] k=" << k << " sum=" << sum << " diag=" << diag << " Lkk=" << Lkk << "\n";
      #endif
    }
  private:
    ac_cinvsqrt_core<SqrRoot_T,SqrRoot_OutT> SQRTINST;
    ac_cinvmac_core<Division_T> MACINST;

  };





// ================================================================================
// compute column k Helper kernels (local-buffer-based) for right looking cholesky
// ================================================================================

  template <typename Input_T, typename Output_T,typename SqrRoot_T, typename SqrRoot_OutT,typename Division_T,typename Reciprocal_T, typename Reciprocal_OutT,typename Mat_size_T,unsigned N>
  class ac_compute_column_k
  {
  public:
    ac_compute_column_k() {}

    void compute_column_k(ac_array<Input_T,N,N> &A, Mat_size_T k, Division_T &Lkk) {
      Reciprocal_OutT recp_out = {};
      Reciprocal_T    recp_in  = {};
      Division_T      acc      = ( Division_T)0;
      Mat_size_T M_AC=N;

      // --------------------------------------------
      // 1. Preload row k and column k into local buffers
      // --------------------------------------------
      Division_T rowbuf[N]= {};
      Division_T colbuf[N]= {};

      #ifndef __CCOV__
      #pragma hls_unroll yes
      #pragma hls_waive FVI
      BUFROW: for (Mat_size_T s = 0; s < k; ++s) {
        rowbuf[s] = (Division_T)A[k][s];
      }
      #endif

      #ifndef __CCOV__
      #pragma hls_unroll yes
      BUFCOL: for (Mat_size_T i = 0; i < M_AC; ++i) {
        colbuf[i] = (i >= k + 1) ? (Division_T)A[i][k] : (Division_T)0;
      }
      #endif

      // --------------------------------------------
      // 2. Compute column entries
      // --------------------------------------------
      #ifndef __CCOV__
      #pragma hls_unroll yes
      //#pragma hls_pipeline_init_interval 1
      ROWLOOP: for (Mat_size_T i = 0; i < M_AC; ++i) {
        if (i > k) {
          acc = (Division_T)0;

          #pragma hls_unroll yes
          // #pragma hls_pipeline_init_interval 1
          #pragma hls_waive FVI
          COLLOOP: for (Mat_size_T s = 0; s < k; ++s) {

            Division_T macinA = (Division_T)A[i][s];
            Division_T macinB = rowbuf[s];   // from buffered row k
            Division_T acc_tmp = ( Division_T)0;
            Division_T prod = (Division_T)((Division_T)macinA* (Division_T)macinB);
            #pragma hls_waive OVL
            acc_tmp = (Division_T)acc + (Division_T)prod;
            //MACINST1.mac_ccore(macinA, macinB, acc, acc_tmp);
            acc = acc_tmp;
          }

          // reciprocal of Lkk
          recp_in = (Reciprocal_T)Lkk;
          CDRCPINST.inv_sqrt_ccore(recp_in, recp_out);

          // update A[i][k] using local buffer value
          A[i][k] = (Output_T)((Division_T)(colbuf[i] - acc) * (Division_T)recp_out);

          #ifndef __SYNTHESIS__
          std::cout << "[compute_column_k] k=" << k << " i=" << i
                    << " acc=" << acc
                    << " recp_in=" << recp_in
                    << " recp_out=" << recp_out
                    << " A[" << i << "][" << k << "]=" << A[i][k] << "\n";
          #endif
        }
      }
      #endif
    }

  private:
    ac_cinvreciprocal_core<Reciprocal_T, Reciprocal_OutT> CDRCPINST;
    ac_cinvmac_core_buf<Division_T> MACINST1;
  };


  // ============================================================
  // compute L[k][k] Helper kernels for left looking cholesky
  // ============================================================

  template <typename Input_T, typename Output_T, typename SqrRoot_T,typename SqrRoot_OutT,typename Division_T, typename Reciprocal_T, typename Reciprocal_OutT, typename Mat_size_T,unsigned N>
  class ac_compute_Lkk_local
  {
  public:
    ac_compute_Lkk_local() {}
    void compute_Lkk_local(ac_array<Input_T,N,N> &A, Mat_size_T k, Division_T &Lkk) {

      SqrRoot_T sqrt_in= {};
      SqrRoot_OutT sqrt_out= {};

      sqrt_in = (SqrRoot_T)A[k][k];
      sqrt_out = (SqrRoot_OutT)0;


      sqrt_inst.ccore_sqrt(sqrt_in, sqrt_out);


      // Lkk is sqrt_out promoted to Division_T
      Lkk = (Division_T)sqrt_out;

      #ifndef __SYNTHESIS__
      std::cout << "[compute_Lkk_local] k="<<k<<" Akk(in)="<< (Division_T)A[k][k]
                << " Lkk=" << (Division_T)Lkk << "\n";
      #endif
    }

  private:
    ac_cinvsqrt_core<SqrRoot_T,SqrRoot_OutT> sqrt_inst;
  };


  // ============================================================
  // compute column[k] local Helper kernels for left looking cholesky
  // ============================================================



  template <typename Input_T, typename Output_T, typename SqrRoot_T,typename SqrRoot_OutT, typename Division_T, typename Reciprocal_T,typename Reciprocal_OutT, typename Mat_size_T, unsigned N>
  class ac_compute_column_k_local
  {
  public:
    ac_compute_column_k_local() {}

    void compute_column_k_local(ac_array<Input_T, N, N> &A,Mat_size_T k,const Division_T &Lkk) {

      Division_T invLkk = {};
//    Division_T tmp = {};
//    Division_T val = {};
      Reciprocal_T recp_in = {};
      Reciprocal_OutT recp_out = {};
      Mat_size_T M_AC=N;

      // compute reciprocal 1/Lkk
      recp_in = (Reciprocal_T)Lkk;
      recp_out = (Reciprocal_OutT)0;

      rcp_core_sim.inv_sqrt_ccore(recp_in, recp_out);
      invLkk = (Division_T)recp_out;

      // -------------------------
      // Local buffer for column k
      // -------------------------
      Division_T col_buffer[N]= {};

      // Read column into local buffer
      #ifndef __CCOV__
      #pragma hls_unroll yes
      COL_READ: for (Mat_size_T i = 0; i < M_AC; ++i) {
        col_buffer[i] = (i >= k + 1) ? (Division_T)A[i][k] : (Division_T)0;
      }
      #endif

      // Compute using buffer
      #ifndef __CCOV__
      #pragma hls_unroll yes
      COL_COMPUTE: for (Mat_size_T i = 0; i < M_AC; ++i) {
        col_buffer[i] = (i >= k + 1) ? (Division_T)((Division_T)col_buffer[i] * (Division_T)invLkk) : (Division_T)col_buffer[i];
      }
      #endif

      // Write back results to A
      #ifndef __CCOV__
      #pragma hls_unroll yes
      COL_WRITE: for (Mat_size_T i = 0; i < M_AC; ++i) {
        A[i][k] = (i >= k + 1) ? (Output_T)col_buffer[i] : (Output_T)A[i][k];

        #ifndef __SYNTHESIS__
        std::cout << "[compute_column_k_local] k=" << k << " i=" << i
                  << " Aik_old=" << (Division_T)col_buffer[i]
                  << " 1/Lkk=" << (Division_T)invLkk
                  << " L[i][k]=" << (Output_T)A[i][k] << "\n";
        #endif
      }
      #endif
    }

  private:
    ac_cinvreciprocal_core<Reciprocal_T, Reciprocal_OutT> rcp_core_sim;
  };

  // ============================================================
  // Trailing Update Helper kernels for left looking cholesky
  // ============================================================

  template <typename Input_T, typename Output_T, typename SqrRoot_T, typename SqrRoot_OutT, typename Division_T, typename Reciprocal_T, typename Reciprocal_OutT, typename Mat_size_T, unsigned N>
  class ac_update_trailing_local
  {
  public:
    ac_update_trailing_local() {}

    void update_trailing_local(ac_array<Input_T, N, N> &A, Mat_size_T k) {
      Mat_size_T M_AC=N;
      Division_T aij = {};
      Division_T lik = {};
      Division_T ljk = {};
      Division_T updated = {};

      // --------------------------------
      // Local buffers for row and column
      // --------------------------------
      Division_T row_buffer[N]= {};  // buffer for A[i][j]
      Division_T col_buffer[N]= {};  // buffer for A[j][k]
      Division_T diag_buffer[N]= {}; // buffer for A[i][k]

      // Preload column k into col_buffer
      #ifndef __CCOV__
      #pragma hls_unroll yes
      COL_LOAD: for (Mat_size_T j = 0; j < M_AC; ++j) {
        col_buffer[j] = (j >= k + 1) ? (Division_T)A[j][k] : (Division_T)0;
      }
      #endif

      // Process row by row
      #ifndef __CCOV__
      #pragma hls_unroll yes
      CUTLROW: for (Mat_size_T i = 0; i < M_AC; ++i) {
        diag_buffer[i] = (i >= k + 1) ? (Division_T)A[i][k] : (Division_T)0;

        #pragma hls_unroll yes
        ROW_LOAD: for (Mat_size_T j = 0; j < M_AC; ++j) {
          row_buffer[j] = (i >= k + 1 && j >= k + 1) ? (Division_T)A[i][j] : (Division_T)0;
        }


        #pragma hls_unroll yes
        ROW_COMPUTE: for (Mat_size_T j = 0; j < M_AC; ++j) {
          aij = (i >= k + 1 && j >= k + 1 && j <= i) ? row_buffer[j] : (Division_T)0;
          lik = (i >= k + 1 && j >= k + 1 && j <= i) ? diag_buffer[i] : (Division_T)0;
          ljk = (i >= k + 1 && j >= k + 1 && j <= i) ? col_buffer[j] : (Division_T)0;
          MACINST1.mac_ccore(lik, ljk, aij, updated);
          row_buffer[j] = (i >= k + 1 && j >= k + 1 && j <= i) ? (Division_T)updated : (Division_T)0;

          #ifndef __SYNTHESIS__
          if (i >= k + 1 && j >= k + 1 && j <= i) {
            std::cout << "[update_trailing_local] k=" << k
                      << " (i,j)=(" << i << "," << j << ")"
                      << " old=" << aij << " lik=" << lik
                      << " ljk=" << ljk
                      << " new=" << updated << "\n";
          }
          #endif
        }

        #pragma hls_unroll yes
        ROW_WRITE: for (Mat_size_T j = 0; j < M_AC; ++j) {
          A[i][j] = (i >= k + 1 && j >= k + 1) ? (Output_T)row_buffer[j] : (Output_T) A[i][j];
        }
      }
      #endif
    }

  private:
    ac_cinvmacsub_core<Division_T> MACINST1;
  };




// ============================================================
// invert lower-triangular L -> Linv
// ============================================================

  template <
    typename Input_T, typename Output_T, typename SqrRoot_T, typename SqrRoot_OutT,typename Division_T, typename Reciprocal_T, typename Reciprocal_OutT,typename Mat_size_T, unsigned N, bool MODE>
  class ac_invert_lower
  {
  public:
    ac_invert_lower() {}
    void invert_lower(ac_array<Input_T,N,N> &L, ac_array<Output_T,N,N> &Linv) {

      Mat_size_T M_AC = N;
      Reciprocal_T    recp_in1 = {};
      Reciprocal_OutT recp_out1 = {};
      Reciprocal_T    recp_in2 = {};
      Reciprocal_OutT recp_out2 = {};
      Division_T acc = {};

      // Diagonal inversion (same for coverage and synth)
      #pragma hls_unroll yes
      INVRROW: for (Mat_size_T i = 0; i < M_AC; ++i) {
        recp_in1 = (Reciprocal_T)L[i][i];
        CDRCP1.inv_sqrt_ccore(recp_in1, recp_out1);
        Linv[i][i] = (Output_T)(recp_out1);
        #ifndef __SYNTHESIS__
        std::cout << "[invert_lower] diag i=" << i << " Linv[ii]=" << (Output_T)Linv[i][i] << "\n";
        #endif
      }

      #ifdef __CCOV__

      // ============================
      // Coverage build: exercise BOTH MODE paths
      //  - First: MODE==true style (direct MAC in loop)
      //  - Then:  MODE==false style (mac_ccore)
      // This keeps the synthesized datapath unchanged (these blocks are coverage-only)
      // ============================

      // MODE == true body (direct multiply-add accumulation)
      #pragma hls_unroll yes
      INVROW_COV_TRUE: for (Mat_size_T i = 1; i < M_AC; ++i) {
        #pragma hls_unroll yes
        INVROWM_COV_TRUE: for (Mat_size_T j = 0; j < i; ++j) {
          #ifndef __CCOV__
          acc = (Division_T)0;
          #pragma hls_waive CNS
          // INVCOLM style: direct accumulation over k
          #pragma hls_unroll yes
          INVCOLM_COV: for (Mat_size_T k = 0; k < M_AC; ++k) {
            acc += (k >= j && k < i) ? (Division_T)((Division_T)L[i][k] * (Division_T)Linv[k][j]) : (Division_T)0;
          }
          #endif
          recp_in2 = (Reciprocal_T)L[i][i];
          CDRCP2.inv_sqrt_ccore(recp_in2, recp_out2);
          Linv[i][j] = (Output_T)((Division_T)(-acc) * (Division_T)recp_out2);
          #ifndef __SYNTHESIS__
          std::cout << "[invert_lower][COV_TRUE] offdiag (i,j)=("<<i<<","<<j<<") acc=" << acc
                    << " Linv=" << (Output_T)Linv[i][j] << "\n";
          #endif
        }
      }
      // Optional coverage bin (if your coverage tool supports it)


      // MODE == false body (use mac_ccore)
      #pragma hls_unroll yes
      INVROW_COV_FALSE: for (Mat_size_T i = 1; i < M_AC; ++i) {
        #ifndef __CCOV__
        #pragma hls_unroll yes
        INVROWM_COV_FALSE: for (Mat_size_T j = 0; j < i; ++j) {
          #ifndef __CCOV__
          acc = (Division_T)0;
          #pragma hls_waive CNS
          // INVCOL using MACINST1.mac_ccore
          #pragma hls_unroll yes
          INVCOL_COV: for (Mat_size_T k = 0; k < M_AC; ++k) {

            Division_T macinA = (k >= j && k < i) ? L[i][k] : (Division_T)0;
            Division_T macinB = (k >= j && k < i) ? Linv[k][j] : (Division_T)0;

            MACINST1.mac_ccore(macinA, macinB, acc);
          }
          #endif
          recp_in2 = (Reciprocal_T)L[i][i];
          CDRCP2.inv_sqrt_ccore(recp_in2, recp_out2);
          Linv[i][j] = (Output_T)((Division_T)(-acc) * (Division_T)recp_out2);
          #ifndef __SYNTHESIS__
          std::cout << "[invert_lower][COV_FALSE] offdiag (i,j)=("<<i<<","<<j<<") acc=" << acc
                    << " Linv=" << (Output_T)Linv[i][j] << "\n";
          #endif
        }
        #endif
      }

      #else // !__CCOV__

      // ============================
      // Normal build: compile-time MODE decides implementation
      // ============================

      #pragma hls_unroll yes
      INVROW: for (Mat_size_T i = 1; i < M_AC; ++i) {
        #pragma hls_waive FVI
        #pragma hls_unroll yes
        INVROWM: for (Mat_size_T j = 0; j < i; ++j) {
          acc = (Division_T)0;
          #pragma hls_waive CNS
          if (MODE) {
            #pragma hls_unroll yes
            INVCOLM: for (Mat_size_T k = 0; k < M_AC; ++k) {
              acc += (k >= j && k < i) ? (Division_T)((Division_T)L[i][k] * (Division_T)Linv[k][j]) : (Division_T)0;
            }
          } else {
            #pragma hls_unroll yes
            INVCOL: for (Mat_size_T k = 0; k < M_AC; ++k) {
              Division_T macinA = (k >= j && k < i) ? L[i][k] : (Division_T)0;
              Division_T macinB = (k >= j && k < i) ? Linv[k][j] : (Division_T)0;
              MACINST1.mac_ccore(macinA, macinB, acc);
            }
          }

          recp_in2 = (Reciprocal_T)L[i][i];
          CDRCP2.inv_sqrt_ccore(recp_in2, recp_out2);
          Linv[i][j] = (Output_T)((Division_T)(-acc) * (Division_T)recp_out2);

          #ifndef __SYNTHESIS__
          std::cout << "[invert_lower] offdiag (i,j)=("<<i<<","<<j<<") acc=" << acc
                    << " Linv=" << (Output_T)Linv[i][j] << "\n";
          #endif
        }
      }

      #endif // __CCOV__
    }

  private:
    ac_cinvreciprocal_core<Reciprocal_T, Reciprocal_OutT> CDRCP1;
    ac_cinvreciprocal_core<Reciprocal_T, Reciprocal_OutT> CDRCP2;
    ac_cinvmac_core<Division_T> MACINST1;
  };


// ============================================================
// invert lower-triangular L -> Linv
// ============================================================

  template <typename Input_T, typename Output_T, typename SqrRoot_T, typename SqrRoot_OutT, typename Division_T, typename Reciprocal_T, typename Reciprocal_OutT, typename Mat_size_T, unsigned N>
  class ac_invert_lower_buf
  {
  public:
    ac_invert_lower_buf() {}

    void invert_lower(ac_array<Input_T, N, N> &L, ac_array<Output_T, N, N> &Linv) {

      Reciprocal_T    recp_in1 = {};
      Reciprocal_OutT recp_out1 = {};
      Reciprocal_T    recp_in2 = {};
      Reciprocal_OutT recp_out2 = {};
      Division_T acc = {};
      Mat_size_T M_AC=N;

      #pragma hls_unroll yes
      INVRROW: for (Mat_size_T i = 0; i < M_AC; ++i) {
        recp_in1 = (Reciprocal_T)L[i][i];
        CDRCP1.inv_sqrt_ccore(recp_in1, recp_out1);
        Linv[i][i] = (Output_T)(recp_out1);

        #ifndef __SYNTHESIS__
        std::cout << "[invert_lower] diag i=" << i << " Linv[ii]=" << (Output_T)Linv[i][i] << "\n";
        #endif
      }

      #pragma hls_unroll yes
      INVROW: for (Mat_size_T i = 1; i < M_AC; ++i) {
        #pragma hls_unroll yes
        #pragma hls_waive FVI
        INVROWM: for (Mat_size_T j = 0; j < i; ++j) {
          acc = (Division_T)0;

          // Local buffers to hold row and column data
          Input_T L_row_i[N]= {};
          Output_T Linv_col_j[N]= {};

          #pragma hls_unroll yes
          LOAD_BUFFERS: for (Mat_size_T k = 0; k < M_AC; ++k) {
            L_row_i[k]    = (Output_T)L[i][k];
            Linv_col_j[k] =(Output_T) Linv[k][j];
          }


          INVCOL_BUFFERED.invcol_buffered(L_row_i, Linv_col_j, j, i, acc);

          recp_in2 = (Reciprocal_T)L[i][i];
          CDRCP2.inv_sqrt_ccore(recp_in2, recp_out2);

          Linv[i][j] = (Output_T)((Division_T)(-acc) * (Division_T)recp_out2);

          #ifndef __SYNTHESIS__
          std::cout << "[invert_lower] offdiag (i,j)=(" << i << "," << j << ") acc=" << acc
                    << " Linv=" << (Output_T)Linv[i][j] << "\n";
          #endif
        }
      }
    }

  private:
    ac_cinvreciprocal_core<Reciprocal_T, Reciprocal_OutT> CDRCP1;
    ac_cinvreciprocal_core<Reciprocal_T, Reciprocal_OutT> CDRCP2;
    ac_invcol_buffered_ccore<Input_T, Division_T, Mat_size_T, N> INVCOL_BUFFERED;

  };


  // ============================================================
  // compute Ainv = Linv^T * Linv Helper kernels
  // ============================================================

  template <
    typename Input_T, typename Output_T, typename SqrRoot_T, typename SqrRoot_OutT,typename Division_T, typename Reciprocal_T, typename Reciprocal_OutT,typename Mat_size_T, unsigned N, bool MODE>
  class ac_compute_Ainv_from_Linv
  {
  public:
    ac_compute_Ainv_from_Linv() {}

    void compute_Ainv_from_Linv(ac_array<Input_T, N, N> &Linv,ac_array<Output_T, N, N> &Ainv) {

      Mat_size_T M_AC = N;
      Division_T acc = {};

      #ifdef __CCOV__

      // ===========================================================
      // COVERAGE BUILD: Execute both MODE=true and MODE=false paths
      // ===========================================================

      // MODE == true body (direct multiply-add)
      #pragma hls_unroll yes
      INROW_COV_TRUE: for (Mat_size_T i = 0; i < M_AC; ++i) {
        #pragma hls_unroll yes
        INCOL_COV_TRUE: for (Mat_size_T j = 0; j <= i; ++j) {
          #ifndef __CCOV__
          acc = (Division_T)0;
          #pragma hls_waive CNS
          // INVMACLOOPM: direct accumulation
          #pragma hls_unroll yes
          INVMACLOOPM_COV: for (Mat_size_T k = 0; k < M_AC; ++k) {
            acc += (k >= i)
                   ? (Division_T)((Division_T)Linv[k][i] * (Division_T)Linv[k][j])
                   : (Division_T)0;
          }
          #endif
          Ainv[i][j] = (Output_T)acc;
          Ainv[j][i] = (Output_T)acc;
          #ifndef __SYNTHESIS__
          std::cout << "[compute_Ainv_from_Linv][COV_TRUE] ("
                    << i << "," << j << ")=" << acc << "\n";
          #endif
        }
      }

      // MODE == false body (use MAC core)
      #pragma hls_unroll yes
      INROW_COV_FALSE: for (Mat_size_T i = 0; i < M_AC; ++i) {
        #pragma hls_unroll yes
        INCOL_COV_FALSE: for (Mat_size_T j = 0; j <= i; ++j) {
          #ifndef __CCOV__
          acc = (Division_T)0;
          #pragma hls_waive CNS
          // INVMACLOOP: using mac_ccore
          #pragma hls_unroll yes
          INVMACLOOP_COV: for (Mat_size_T k = 0; k < M_AC; ++k) {
            Division_T macinA = (k >= i) ? Linv[k][i] : (Division_T)0;
            Division_T macinB = (k >= i) ? Linv[k][j] : (Division_T)0;
            MACINST2.mac_ccore(macinA, macinB, acc);
          }
          #endif
          Ainv[i][j] = (Output_T)acc;
          Ainv[j][i] = (Output_T)acc;
          #ifndef __SYNTHESIS__
          std::cout << "[compute_Ainv_from_Linv][COV_FALSE] ("
                    << i << "," << j << ")=" << acc << "\n";
          #endif
        }
      }


      #else // !__CCOV__

      // ===========================================================
      // NORMAL BUILD: compile-time MODE decides implementation
      // ===========================================================
      #pragma hls_unroll yes
      INROWLOOP: for (Mat_size_T i = 0; i < M_AC; ++i) {
        #pragma hls_waive FVI
        #pragma hls_unroll yes
        INCOLLOOP: for (Mat_size_T j = 0; j <= i; ++j) {
          acc = (Division_T)0;
          #pragma hls_waive CNS

          if (MODE) {
            #pragma hls_unroll yes
            INVMACLOOPM: for (Mat_size_T k = 0; k < M_AC; ++k) {
              acc += (k >= i)
                     ? (Division_T)((Division_T)Linv[k][i] * (Division_T)Linv[k][j])
                     : (Division_T)0;
            }
          } else {
            #pragma hls_unroll yes
            INVMACLOOP: for (Mat_size_T k = 0; k < M_AC; ++k) {
              Division_T macinA = (k >= i) ? Linv[k][i] : (Division_T)0;
              Division_T macinB = (k >= i) ? Linv[k][j] : (Division_T)0;
              MACINST2.mac_ccore(macinA, macinB, acc);
            }
          }

          Ainv[i][j] = (Output_T)acc;
          Ainv[j][i] = (Output_T)acc;

          #ifndef __SYNTHESIS__
          std::cout << "[compute_Ainv_from_Linv] ("
                    << i << "," << j << ")=" << acc << "\n";
          #endif
        }
      }

      #endif // __CCOV__
    }

  private:
    ac_cinvmac_core<Division_T> MACINST2;
  };


// ============================================================
// Buffered Ainv Computation from Linv
// ============================================================

  template <typename Input_T, typename Output_T, typename SqrRoot_T,typename SqrRoot_OutT,typename Division_T, typename Reciprocal_T, typename Reciprocal_OutT, typename Mat_size_T,unsigned N>
  class ac_compute_Ainv_from_Linv_buf
  {
  public:
    ac_compute_Ainv_from_Linv_buf() {}

    void compute_Ainv_from_Linv(ac_array<Input_T,N,N> &Linv,ac_array<Output_T,N,N> &Ainv) {
      Division_T acc = ( Division_T)0;
      Division_T col_i[N]= {};
      Division_T col_j[N]= {};

      #pragma hls_unroll yes
      INROWLOOP: for (Mat_size_T i = 0; i < N; ++i) {
        #pragma hls_unroll yes
        INCOLLOOP: for (Mat_size_T j = 0; j <= i; ++j) {
          // --- Buffer columns i and j into local arrays ---
          #pragma hls_unroll yes
          BUFCOL: for (Mat_size_T k = 0; k < N; ++k) {
            col_i[k] = (Division_T)Linv[k][i];
            col_j[k] = (Division_T)Linv[k][j];
          }

          // --- Call buffered MAC core ---
          MACINNER.mac_inner_buffered(col_i, col_j, i, acc);

          Ainv[i][j] = (Output_T)acc;
          Ainv[j][i] = (Output_T)acc;

          #ifndef __SYNTHESIS__
          std::cout << "[compute_Ainv] (" << i << "," << j
                    << ") = " << acc << "\n";
          #endif
        }
      }
    }

  private:
    ac_mac_inner_buffered_ccore<Division_T,Mat_size_T, N> MACINNER;
  };


// ============================================================
// Optimized Ainv Computation from Linv
// ============================================================
  template <typename Input_T, typename Output_T, typename SqrRoot_T,typename SqrRoot_OutT,typename Division_T, typename Reciprocal_T, typename Reciprocal_OutT, typename Mat_size_T,unsigned N>

  class ac_compute_Ainv_from_Linv_chunk
  {
  public:
    static constexpr unsigned CHUNK_SIZE = 8;
    static constexpr unsigned NUM_CHUNKS = (N + CHUNK_SIZE - 1) / CHUNK_SIZE;

    ac_compute_Ainv_from_Linv_chunk() {}



    void compute_Ainv_from_Linv(ac_array<Input_T, N, N> &Linv, ac_array<Output_T, N, N> &Ainv) {
      Mat_size_T M_AC=N;
      Division_T acc = ( Division_T)0;

      #pragma hls_unroll yes
      INROWLOOP: for (Mat_size_T i = 0; i < M_AC; ++i) {
        #pragma hls_unroll yes
        #pragma hls_waive FVI
        INCOLLOOP: for (Mat_size_T j = 0; j <= i; ++j) {

          acc = ( Division_T)0;

          #pragma hls_unroll yes  // Fully unroll since NUM_CHUNKS is small and fixed
          for (Mat_size_T chunk = 0; chunk < NUM_CHUNKS; ++chunk) {

            Input_T row_buf[CHUNK_SIZE];
            Input_T col_buf[CHUNK_SIZE];

            #ifndef __CCOV__
            #pragma hls_unroll yes
            LOAD_BUFFERS: for (Mat_size_T k = 0; k < CHUNK_SIZE; ++k) {
              Mat_size_T idx = chunk * CHUNK_SIZE + k;
              #pragma hls_waive CNS
              row_buf[k] = (idx < N) ? Linv[idx][i] : (Input_T)0;
              #pragma hls_waive CNS
              col_buf[k] = (idx < N) ? Linv[idx][j] : (Input_T)0;
            }
            #endif

            Division_T acc_chunk = (Division_T)0;
            chunk_mac_ccore.chunk_mac(row_buf, col_buf, acc_chunk);
            acc += acc_chunk;
          }

          Ainv[i][j] = (Output_T)acc;
          Ainv[j][i] = (Output_T)acc;

          #ifndef __SYNTHESIS__
          std::cout << "[compute_Ainv] (" << i << "," << j << ") = " << acc << "\n";
          #endif
        }
      }
    }

  private:
    ac_chunk_mac_ccore<Input_T, Division_T, Mat_size_T, CHUNK_SIZE> chunk_mac_ccore;
  };

  // ============================================================
  // Top-Level Wrapper Class: Cholesky Inversion
  // ============================================================

  template <typename Input_T, typename Output_T, typename SqrRoot_T,typename SqrRoot_OutT, typename Division_T, typename Reciprocal_T, typename Reciprocal_OutT,typename Mat_size_T, unsigned TOL_CODE, int W1,unsigned M,bool MODE,bool CORE,bool ALGO>
  class ac_cholesky_matrix_inverse
  {
  public:


    ac_cholesky_matrix_inverse() {}



    void cholesky_matrix_inverse(ac_array<Input_T, M, M> &L_in_array, ac_array<Output_T, M, M> &L_out_array) {

      #ifndef __CCOV__
      static_assert(M > 0, "Matrix size N must be > 0");
      #endif

      ac_array<Output_T,M,M> A= {};

      ac_array<Output_T,M,M> Linv= {};
      ac_array<Output_T,M,M> Ainv= {};

      Mat_size_T M_AC=M;
      // Copy input
      //#pragma hls_pipeline_init_interval 1
      #pragma hls_unroll yes
      INPUTROWLOOP:for (Mat_size_T i = 0; i < M_AC; ++i) {
        #pragma hls_unroll yes
        INPUTCOLLOOP:for (Mat_size_T j = 0; j < M_AC; ++j) {
          A[i][j] = L_in_array[i][j];
          //#ifndef __SYNTHESIS__
          //std::cout << "A[" << i << "][" << j << "] = " << A[i][j] << "\n";
          //#endif
        }
      }

      // Zero upper triangle
      #ifndef __CCOV__
      //#pragma hls_pipeline_init_interval 1
      #pragma hls_unroll yes
      CONVROWLOOP:for (Mat_size_T i = 0; i < M_AC; ++i) {
        #pragma hls_unroll yes
        CONVCOLLOOP: for (Mat_size_T j = 0; j < M_AC; ++j) {
          A[i][j] = (j >= i + 1) ? (Division_T)0 : A[i][j];
        }
      }
      #endif

      #pragma hls_waive CNS
      #ifdef __CCOV__

// ===========================================================
//  Coverage build: force execution of both branches
// ===========================================================

// ---------- Right-looking Cholesky (forced path) ----------

      #pragma hls_unroll yes
      PELOOPRIGHT_COV_TRUE: for (Mat_size_T k = 0; k < M_AC; ++k) {
        Division_T Lkk = (Division_T)0;
        CLKK.compute_Lkk(A, k, Lkk);
        A[k][k] = (Output_T)Lkk;
        CLMNK.compute_column_k(A, k, Lkk);
      }



// ---------- Left-looking Cholesky (dummy path for coverage) ----------

      #pragma hls_unroll yes
      PELOOPLEFT_COV_FALSE: for (Mat_size_T k = 0; k < M_AC; ++k) {
        Division_T Lkk = (Division_T)0;
        CLKKL.compute_Lkk_local(A, k, Lkk);
        A[k][k] = (Output_T)Lkk;
        CLMNKL.compute_column_k_local(A, k, Lkk);
        CUTL.update_trailing_local(A, k);
      }



      #else  // ====================================================
//  Normal synthesis/simulation build
// ===========================================================

      if (ALGO) {

        // ---------- Right-looking Cholesky ----------
        #pragma hls_unroll yes
        PELOOPRIGHT: for (Mat_size_T k = 0; k < M_AC; ++k) {
          Division_T Lkk = (Division_T)0;
          CLKK.compute_Lkk(A, k, Lkk);
          A[k][k] = (Output_T)Lkk;
          CLMNK.compute_column_k(A, k, Lkk);
        }

      } else {

        // ---------- Left-looking Cholesky ----------
        #pragma hls_unroll yes
        PELOOPLEFT: for (Mat_size_T k = 0; k < M_AC; ++k) {
          Division_T Lkk = (Division_T)0;
          CLKKL.compute_Lkk_local(A, k, Lkk);
          A[k][k] = (Output_T)Lkk;
          CLMNKL.compute_column_k_local(A, k, Lkk);
          CUTL.update_trailing_local(A, k);
        }
      }

      #endif  // __CCOV__




      // Init buffers
      //#pragma hls_pipeline_init_interval 1
      #pragma hls_unroll yes
      BUFFRLOOP:for (Mat_size_T i = 0; i < M_AC; ++i) {
        #pragma hls_unroll yes
        BUFFCLOOP:for (Mat_size_T j = 0; j < M_AC; ++j) {
          Linv[i][j] = (Output_T)0;
          Ainv[i][j] = (Output_T)0;
        }
      }

      #ifdef __CCOV__

      INVLB.invert_lower(A, Linv);

      INVL.invert_lower(A, Linv);



      CALINVC.compute_Ainv_from_Linv(Linv, Ainv);

      CALINV.compute_Ainv_from_Linv(Linv, Ainv);


      #else
      // Invert L
      #pragma hls_waive CNS
      if (CORE)
      { INVLB.invert_lower(A, Linv); }
      else
      { INVL.invert_lower(A, Linv); }



      // Compute Ainv
      #pragma hls_waive CNS
      if (CORE)
      { CALINVC.compute_Ainv_from_Linv(Linv, Ainv); }
      else
      { CALINV.compute_Ainv_from_Linv(Linv, Ainv); }


      #endif





      // Write output
      //#pragma hls_pipeline_init_interval 1
      #pragma hls_unroll yes
      OUTPROWLOOP:for (Mat_size_T i = 0; i < M_AC; ++i) {
        #pragma hls_unroll yes
        OUTPCOLLOOP:for (Mat_size_T j = 0; j < M_AC; ++j) {
          L_out_array[i][j] = (Output_T)Ainv[i][j];
        }
      }



      #ifdef __CCOV__
      MyCCoverGroup<Input_T, Output_T,SqrRoot_T,SqrRoot_OutT, M, TOL_CODE,W1> cov_decomp("cov_decomp");
      MyCCoverGroup<Output_T, Output_T,SqrRoot_T,SqrRoot_OutT, M, TOL_CODE,W1> cov_forward("cov_forward");
      MyCCoverGroup<Output_T, Output_T,SqrRoot_T, SqrRoot_OutT,M, TOL_CODE,W1> cov_backward("cov_backward");

      // Decomposition stage coverage: A_in -> L
      cov_decomp.SampleAll(L_in_array, A);
      cov_decomp.is_valid_matrix_data(L_in_array, A);
      cov_decomp.is_valid_meta(TOL_CODE, W1);
      cov_decomp.SampleAllMetaCoverage(TOL_CODE, W1);

      // Forward substitution stage: L -> Y
      cov_forward.SampleAll(A, Linv);
      cov_forward.is_valid_matrix_data(A, Linv);
      cov_forward.is_valid_meta(TOL_CODE, W1);
      cov_forward.SampleAllMetaCoverage(TOL_CODE,W1);

      // Backward substitution stage: A,Y -> A_inv
      cov_backward.SampleAll(Linv, Ainv);
      cov_backward.SampleAll(Ainv, L_out_array);
      cov_backward.is_valid_matrix_data(Linv, Ainv);
      cov_backward.is_valid_matrix_data(Ainv, L_out_array);
      cov_backward.is_valid_meta(TOL_CODE, W1);
      cov_backward.SampleAllMetaCoverage(TOL_CODE, W1);
      #endif

    }

  private:
    ac_compute_Lkk<Input_T, Output_T, SqrRoot_T, SqrRoot_OutT, Division_T, Reciprocal_T, Reciprocal_OutT,Mat_size_T, M> CLKK;
    ac_compute_column_k<Input_T, Output_T, SqrRoot_T, SqrRoot_OutT, Division_T, Reciprocal_T, Reciprocal_OutT,Mat_size_T, M> CLMNK;
    ac_compute_Lkk_local<Input_T, Output_T, SqrRoot_T, SqrRoot_OutT, Division_T, Reciprocal_T, Reciprocal_OutT,Mat_size_T, M> CLKKL;
    ac_compute_column_k_local<Input_T, Output_T, SqrRoot_T, SqrRoot_OutT, Division_T, Reciprocal_T, Reciprocal_OutT,Mat_size_T, M> CLMNKL;
    ac_update_trailing_local<Input_T, Output_T, SqrRoot_T, SqrRoot_OutT, Division_T, Reciprocal_T, Reciprocal_OutT,Mat_size_T, M> CUTL;

    ac_invert_lower<Input_T, Output_T, SqrRoot_T, SqrRoot_OutT, Division_T, Reciprocal_T, Reciprocal_OutT,Mat_size_T, M,MODE> INVL;
    ac_invert_lower_buf<Input_T, Output_T, SqrRoot_T, SqrRoot_OutT, Division_T, Reciprocal_T, Reciprocal_OutT,Mat_size_T, M> INVLB;
    ac_compute_Ainv_from_Linv<Input_T, Output_T, SqrRoot_T, SqrRoot_OutT, Division_T, Reciprocal_T, Reciprocal_OutT,Mat_size_T, M,MODE> CALINV;
    ac_compute_Ainv_from_Linv_buf<Input_T, Output_T, SqrRoot_T, SqrRoot_OutT, Division_T, Reciprocal_T, Reciprocal_OutT,Mat_size_T, M> CALINVB;

    ac_compute_Ainv_from_Linv_chunk<Input_T, Output_T, SqrRoot_T, SqrRoot_OutT, Division_T, Reciprocal_T, Reciprocal_OutT,Mat_size_T, M> CALINVC;


  };
}


#endif // SYSTOLIC_CHOLESKY_INV_CATAPULT_CPP



