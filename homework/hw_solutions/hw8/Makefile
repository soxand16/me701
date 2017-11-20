SWIG_INCLUDES=-I/Users/Sam/opt/anaconda/include/python3.6m -I/Users/Sam/opt/anaconda/lib/python3.6/site-packages/numpy/core/include

# F2PY
f2py:
	f2py -c -m oxandale_samuel_hw8_p1_fortran oxandale_samuel_hw8_p1.f90



# SWIG
swig:
	swig -python -c++ oxandale_samuel_hw8_p1_cpp.i
	g++ -fPIC -c oxandale_samuel_hw8_p1.cpp oxandale_samuel_hw8_p1_cpp_wrap.cxx ${SWIG_INCLUDES}
	g++ -shared oxandale_samuel_hw8_p1.o oxandale_samuel_hw8_p1_cpp_wrap.o ${SWIG_INCLUDES} -o _oxandale_samuel_hw8_p1_cpp.so

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .f90 .o
.PHONY: clean neat

%.o: %.f90
	$(F90) $(F90FLAGS) $(INCLUDE) -c $<