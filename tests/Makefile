
ifeq "$(AC_TYPES_INC)" ""
  $(error Error: The AC_TYPES_INC variable was not set. Cannot locate AC Datatypes installation)
endif

AC_MATH_INC = ../include

ifeq "$(GCC_EXEC)" ""
GCC_EXEC = g++
endif

DIE        = || $(RM) $@
RM        := rm -f
ECHO      := echo

SOURCES_CPP = \
  rtest_ac_log2_pwl.cpp \
  rtest_ac_log_pwl.cpp \
  rtest_ac_sincos_lut.cpp \
  rtest_ac_normalize.cpp \
  rtest_ac_reciprocal_pwl.cpp \
  rtest_ac_exp_pwl.cpp \
  rtest_ac_pow2_pwl.cpp \
  rtest_ac_sqrt_pwl.cpp \
  rtest_ac_inverse_sqrt_pwl.cpp

OBJS = $(SOURCES_CPP:.cpp=.o)

%.o: %.cpp
	-@$(ECHO) "------------------------------ Testcase $< ----------------------------------"
	@$(GCC_EXEC) -g -std=c++11 -I$(AC_TYPES_INC) -I$(AC_MATH_INC) $< -o $@ $(DIE)
	@./$@ $(DIE)

all: $(OBJS)

.PHONY: clean
clean:
	-@$(RM) $(OBJS)
