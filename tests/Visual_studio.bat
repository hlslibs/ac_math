::This is a windows .bat file which runs the ac_ipl unit tests using the visual studio developer command prompt 

::Invoke the visual studio developer command prompt
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\Tools\VsDevCmd.bat"

::Set the path for Mgc_home
set MGC_HOME=C:\MGCNoScan\abeemana\sb\sif\ixn\Mgc_home

::All the unit tests
set SOURCES_CPP=rtest_ac_div;^
  rtest_ac_atan_pwl;^
  rtest_ac_atan_pwl_ha;^
  rtest_ac_atan_pwl_vha;^
  rtest_ac_barrel_shift;^
  rtest_ac_tan_pwl;^
  rtest_ac_sigmoid_pwl;^
  rtest_ac_tanh_pwl;^
  rtest_ac_cholinv;^
  rtest_ac_exp_cordic;^
  rtest_ac_exp2_cordic;^
  rtest_ac_pow2_pwl;^
  rtest_ac_exp_pwl;^
  rtest_ac_determinant;^
  rtest_ac_chol_d;^
  rtest_ac_qrd;^
  rtest_ac_abs;^
  rtest_ac_arccos_cordic;^
  rtest_ac_arcsin_cordic;^
  rtest_ac_atan2_cordic;^
  rtest_ac_exp_pwl;^
  rtest_ac_inverse_sqrt_pwl;^
  rtest_ac_inverse_sqrt_pwl_vha;^
  rtest_ac_log_cordic;^
  rtest_ac_log2_cordic;^
  rtest_ac_log2_pwl;^
  rtest_ac_log_pwl;^
  rtest_ac_array;^
  rtest_ac_matrix;^
  rtest_ac_matrixmul;^
  rtest_ac_normalize;^
  rtest_ac_pow2_pwl;^
  rtest_ac_pow_cordic;^
  rtest_ac_pow_pwl;^
  rtest_ac_reciprocal_pwl;^
  rtest_ac_reciprocal_pwl_ha;^
  rtest_ac_reciprocal_pwl_vha;^
  rtest_ac_shift;^
  rtest_ac_sincos_cordic;^
  rtest_ac_sincos_lut;^
  rtest_ac_softmax_pwl;^
  rtest_ac_sqrt;^
  rtest_ac_sqrt_pwl;^
  rtest_ac_leading;^
  rtest_signed_ac_div_v2;^
  rtest_unsigned_ac_div_v2

::Compile and execute the unit tests
(for %%a in (%SOURCES_CPP%) do ( 
   cl /EHsc /I%MGC_HOME%\shared\include %%a.cpp
   %%a.exe
))
