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
#ifndef _INCLUDED_AC_CHOLESKY_INV_1D_SYSTOLIC_RL_H_
#define _INCLUDED_AC_CHOLESKY_INV_1D_SYSTOLIC_RL_H_

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


#ifdef __CCOV__
#include <ac_cholesky_inv_1d_systolic_rl_covergroup.h>
#endif


namespace ac_math
{

  // ============================================================
  // Reciprocal Core Declaration
  // ============================================================
  template <typename Reciprocal_T, typename Reciprocal_OutT>
  class ac_myreciprocal_core
  {
  public:
    ac_myreciprocal_core() {}
    #pragma hls_design ccore
    #pragma hls_ccore_type sequential

    void inv_sqrt_ccore(Reciprocal_T &sub_res, Reciprocal_OutT &rec_res) {
      ac_math::ac_reciprocal_pwl_vha(sub_res, rec_res);
    }
  };

  // ============================================================
  // Square Root Core Declaration
  // ============================================================

  template <typename SqrRoot_T>
  class ac_mysqrt_core
  {
  public:
    ac_mysqrt_core() {}

    #pragma hls_design ccore
    #pragma hls_ccore_type sequential

    void ccore_sqrt(const SqrRoot_T &sqrt_in, SqrRoot_T &sqrt_out) {
      ac_math::ac_sqrt(sqrt_in, sqrt_out);
    }
  };

  // ============================================================
  // Cholesky Decomposition Stage
  // ============================================================

  template <typename Input_T, typename Output_T, typename SqrRoot_T, typename Division_T, typename Reciprocal_T, typename Reciprocal_OutT, int M>
  class ac_cholesky_decomposition
  {
  public:
    ac_cholesky_decomposition() {}

    void choldecomp(ac_array<Input_T, M, M> &L_in_array, ac_array<Output_T, M, M> &L_out_array) {
      Division_T L[M][M] = {};
      Division_T A_local[M][M] = {};
      Division_T div_res_in, div_res_out = {};
      Division_T diag = {};
      SqrRoot_T sqrt_res_in, sqrt_res = {};
      Reciprocal_T recp_in;
      Reciprocal_OutT recp_out;

      #pragma hls_unroll yes
      LOOP_ReadRow:
      for (int i = 0; i < M; ++i) {
        #pragma hls_unroll yes
        LOOP_ReadCol:
        for (int j = 0; j < M; ++j) {
          A_local[i][j] = L_in_array[i][j];
        }
      }

      #pragma hls_unroll yes
      LOOP_ProcessingElement:
      for (int k = 0; k < M; ++k) {
        sqrt_res_in = (SqrRoot_T)(A_local[k][k]);

        // In the FPGA synthesis flow, the design bypasses the CCore implementation
        // and instead uses a custom RTL-optimized version to enhance QoR metrics
        // such as area utilization, timing closure, and throughput
        #ifdef FPGA
        ac_math::ac_sqrt(sqrt_res_in, sqrt_res);
        #else
        SQRTINST.ccore_sqrt(sqrt_res_in, sqrt_res);
        #endif

        L[k][k] = sqrt_res;
        recp_in = sqrt_res;

        // In the FPGA synthesis flow, the design bypasses the CCore implementation
        // and instead uses a custom RTL-optimized version to enhance QoR metrics
        // such as area utilization, timing closure, and throughput
        #ifdef FPGA
        ac_math::ac_reciprocal_pwl_vha(recp_in, recp_out);
        #else
        CDRCPINST.inv_sqrt_ccore(recp_in, recp_out);
        #endif


        diag = (Division_T)(recp_out);

        // In C design check (simulation), the loop is fully unrolled by factor M
        // to accelerate functional verification. For synthesis targeting FPGA or ASIC,
        // the loop remains rolled to maintain the systolic architecture, which
        // enables efficient pipelining and resource sharing in the final hardware.

        #ifndef __SYNTHESIS__
        #pragma hls_unroll factor = M
        #endif

        #ifdef __SYNTHESIS__
        #ifdef FPGA
        #pragma hls_pipeline_init_interval 1
        #else
        #pragma hls_unroll yes
        #endif
        #endif

        // 3. Column update: L[i][k] = A[i][k] / L[k][k]
        #pragma hls_waive FVI
        #pragma hls_waive CNS
        #pragma hls_waive DBZ

        //CCOVERAGE: is blocked by the pragma intenetionally as only the lower triangular part is updated
        //Upper triangular part is not updated therefore brach coverage will alwyasbe 50%

        #ifndef __CCOV__
        LOOP_DivisionOffDiag:
        for (int i = k + 1; i < M; ++i) {
          div_res_in = A_local[i][k];
          div_res_out = (Division_T)(div_res_in * diag);    // Replace division with multiplication by reciprocal
          L[i][k] = div_res_out;
        }
        #endif

        // 4. Trailing submatrix update:
        // A[i][j] -= L[i][k] * L[j][k]

        // In C design check (simulation), the loop is fully unrolled by factor M
        // to accelerate functional verification. For synthesis targeting FPGA or ASIC,
        // the loop remains rolled to maintain the systolic architecture, which
        // enables efficient pipelining and resource sharing in the final hardware.
        #ifndef __SYNTHESIS__
        #pragma hls_unroll factor = M
        #endif
        #ifdef __SYNTHESIS__
        #ifdef FPGA
        #pragma hls_pipeline_init_interval 1
        #else
        #pragma hls_unroll yes
        #endif
        #endif

        #pragma hls_waive FVI
        LOOP_UpdateOffDiagRow:
        for (int i = k + 1; i < M; ++i) {

          #ifndef __SYNTHESIS__
          #pragma hls_unroll factor = M
          #endif
          #ifdef __SYNTHESIS__
          #pragma hls_unroll yes
          #endif


          #pragma hls_waive NCO
          #pragma hls_waive FVI
          LOOP_UpdateOffDiagCol:
          for (int j = k + 1; j < M; ++j) {
            A_local[i][j] = (j <= i) ? (Division_T)(A_local[i][j] - (L[i][k] * L[j][k])) : (Division_T)(0.0);
          }
        }
      }

      #ifdef FPGA
      #pragma hls_pipeline_init_interval 1
      #else
      #pragma hls_unroll yes
      #endif
      //CCOVERAGE: is blocked by the pragma intenetionally as only the lower triangular part is updated
      //Upper triangular part is not updated therefore brach coverage will alwyasbe 50%

      #ifndef __CCOV__
      LOOP_LowerMatOutRow:
      for (int i = 0; i < M; i++) {
        #ifdef FPGA
        #pragma hls_pipeline_init_interval 1
        #else
        #pragma hls_unroll yes
        #endif
        LOOP_LowerMatOutCol:
        for (int j = 0; j < M; j++) {
          L_out_array[i][j] = (j <= i) ? (Output_T)(L[i][j]) : (Output_T)(0);
        }
      }
      #endif

    }

  private:
    ac_mysqrt_core<SqrRoot_T> SQRTINST;
    ac_myreciprocal_core<Reciprocal_T, Reciprocal_OutT> CDRCPINST;
  };

  // ============================================================
  // Forward Substitution Stage
  // ============================================================

  template <typename Input_T, typename Output_T, typename SqrRoot_T, typename Division_T, typename Reciprocal_T, typename Reciprocal_OutT, int M>
  class ac_forward_substitution
  {
  public:
    ac_forward_substitution() {}

    void forwdsub(ac_array<Output_T, M, M> &L_chol_array, ac_array<Output_T, M, M> &fs_out_array) {
      Division_T fsint_result = {};
      Reciprocal_T rcp_in = {};
      Reciprocal_OutT rcp_out = {};
      Division_T rcp = {};
      Division_T sum = {};

      #pragma hls_pipeline_init_interval 1
      LOOP_FSDataCol:
      for (int col = 0; col < M; ++col) {
        #pragma hls_waive CNS
        #pragma hls_waive DBZ
        #pragma hls_unroll yes
        LOOP_FSDataOuterRow:
        for (int row = 0; row < M; ++row) {
          sum = (row == col) ? (Division_T)(1) : (Division_T)(0);

          #pragma hls_unroll yes
          #pragma hls_waive SAT

          #ifndef __CCOV__

          LOOP_FSDataInnerRow:
          for (int k = 0; k < M; ++k) {
            sum = (k < row) ? (Division_T)(sum - (fs_out_array[k][col] * L_chol_array[row][k])) : sum;
          }
          #endif

          rcp_in = (Reciprocal_T)(L_chol_array[row][row]);

          // In the FPGA synthesis flow, the design bypasses the CCore implementation
          // and instead uses a custom RTL-optimized version to enhance QoR metrics
          // such as area utilization, timing closure, and throughput

          #ifdef FPGA
          ac_math::ac_reciprocal_pwl_vha(rcp_in, rcp_out);
          #else
          FSRCPINST.inv_sqrt_ccore(rcp_in, rcp_out);
          #endif
          rcp = (Division_T)(rcp_out);
          fsint_result = (Division_T)(sum * rcp);
          fs_out_array[row][col] = (Output_T)(fsint_result);
        }
      }
    }

  private:
    ac_myreciprocal_core<Reciprocal_T, Reciprocal_OutT> FSRCPINST;
  };

  // ============================================================
  // Backward Substitution Stage
  // ============================================================

  template <typename Input_T, typename Output_T, typename SqrRoot_T, typename Division_T, typename Reciprocal_T, typename Reciprocal_OutT, int M>
  class ac_backward_substitution
  {
  public:
    ac_backward_substitution() {}

    void backdsub(ac_array<Output_T, M, M> &L_chol_array, ac_array<Output_T, M, M> &fs_out_array, ac_array<Output_T, M, M> &bs_out_array) {
      Division_T sum = {};
      Division_T bsint_result = {};
      Reciprocal_T rcp_in = {};
      Reciprocal_OutT rcp_out = {};
      Division_T rcp = {};

      #pragma hls_pipeline_init_interval 1
      LOOP_BKDSDataCol:
      for (int col = 0; col < M; ++col) {
        #ifndef __SYNTHESIS__
        #pragma hls_unroll factor = M
        #endif
        #ifdef __SYNTHESIS__
        #pragma hls_unroll yes
        #endif

        #pragma hls_waive CNS
        #pragma hls_waive DBZ
        LOOP_BKDSDataOuterRow:
        for (int row = M - 1; row >= 0; --row) {
          sum = fs_out_array[row][col];

          #ifndef __SYNTHESIS__
          #pragma hls_unroll factor = M
          #endif
          #ifdef __SYNTHESIS__
          #pragma hls_unroll yes
          #endif

          #pragma hls_waive FVI

          //CCOVERAGE: is blocked by the pragma intenetionally as only the lower triangular part is updated
          //Upper triangular part is not updated therefore brach coverage will alwyasbe 50%
          #ifndef __CCOV__
          LOOP_BKDSDataInnerRow:
          for (int k = row + 1; k < M; ++k) {
            sum = (Division_T)(sum - (L_chol_array[k][row] * bs_out_array[k][col]));
          }
          #endif

          rcp_in = (Reciprocal_T)(L_chol_array[row][row]);

          // In the FPGA synthesis flow, the design bypasses the CCore implementation
          // and instead uses a custom RTL-optimized version to enhance QoR metrics
          // such as area utilization, timing closure, and throughput

          #ifdef FPGA
          ac_math::ac_reciprocal_pwl_vha(rcp_in, rcp_out);
          #else
          BKDRCPINST.inv_sqrt_ccore(rcp_in, rcp_out);
          #endif

          rcp = (Division_T)(rcp_out);

          bsint_result = (Division_T)(sum * rcp);
          bs_out_array[row][col] = (Output_T)(bsint_result);
        }
      }
    }

  private:
    ac_myreciprocal_core<Reciprocal_T, Reciprocal_OutT> BKDRCPINST;
  };

  // ============================================================
  // Top-Level Wrapper Class: Cholesky Inversion
  // ============================================================

  template <typename Input_T, typename Output_T, typename SqrRoot_T, typename Division_T, typename Reciprocal_T, typename Reciprocal_OutT, int TOL_CODE, int W1,int M>
  class ac_cholesky_inv_1d_systolic_rl
  {
  public:
    ac_cholesky_inv_1d_systolic_rl() {}

    void cholinv(ac_array<Input_T, M, M> &A_in_array, ac_array<Output_T, M, M> &Ainv_out_array) {
      ac_array<Output_T, M, M> L, Y= {};

      CDINST.choldecomp(A_in_array, L);
      FSINST.forwdsub(L, Y);
      BSINST.backdsub(L, Y, Ainv_out_array);



      #ifdef __CCOV__
      MyCCoverGroup<Input_T, Output_T,SqrRoot_T, M, TOL_CODE,W1> cov_decomp("cov_decomp");
      MyCCoverGroup<Output_T, Output_T,SqrRoot_T, M, TOL_CODE,W1> cov_forward("cov_forward");
      MyCCoverGroup<Output_T, Output_T,SqrRoot_T, M, TOL_CODE,W1> cov_backward("cov_backward");

      // Decomposition stage coverage: A_in -> L
      cov_decomp.SampleAll(A_in_array, L);
      cov_decomp.is_valid_matrix_data(A_in_array, L);
      cov_decomp.is_valid_meta(TOL_CODE, W1);
      cov_decomp.SampleAllMetaCoverage(TOL_CODE, W1);

      // Forward substitution stage: L -> Y
      cov_forward.SampleAll(L, Y);
      cov_forward.is_valid_matrix_data(L, Y);
      cov_forward.is_valid_meta(TOL_CODE, W1);
      cov_forward.SampleAllMetaCoverage(TOL_CODE,W1);

      // Backward substitution stage: A,Y -> A_inv
      cov_backward.SampleAll(Y, Ainv_out_array);
      cov_backward.SampleAll(L, Ainv_out_array);
      cov_backward.is_valid_matrix_data(Y, Ainv_out_array);
      cov_backward.is_valid_matrix_data(L, Ainv_out_array);
      cov_backward.is_valid_meta(TOL_CODE, W1);
      cov_backward.SampleAllMetaCoverage(TOL_CODE, W1);
      #endif

    }

  private:
    ac_cholesky_decomposition<Input_T, Output_T, SqrRoot_T, Division_T, Reciprocal_T, Reciprocal_OutT, M> CDINST;
    ac_forward_substitution<Input_T, Output_T, SqrRoot_T, Division_T, Reciprocal_T, Reciprocal_OutT, M> FSINST;
    ac_backward_substitution<Input_T, Output_T, SqrRoot_T, Division_T, Reciprocal_T, Reciprocal_OutT, M> BSINST;
  };

}

#endif // _INCLUDED_AC_CHOLESKY_INV_1D_SYSTOLIC_RL_H_
