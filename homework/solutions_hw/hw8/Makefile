SWIG_INCLUDES=-I/Users/Sam/opt/anaconda/include/python3.6m -I/Users/Sam/opt/anaconda/lib/python3.6/site-packages/numpy/core/include

# F2PY
f2py:
	f2py -c -m oxandale_samuel_hw8_p1_fortran oxandale_samuel_hw8_p1.f90



# SWIG
swig:
	swig -python -c++ oxandale_samuel_hw8_p1_cpp.i
	g++ -fPIC -c oxandale_samuel_hw8_p1.cpp oxandale_samuel_hw8_p1_cpp_wrap.cxx ${SWIG_INCLUDES}
	g++ -shared oxandale_samuel_hw8_p1.o oxandale_samuel_hw8_p1_cpp_wrap.o ${SWIG_INCLUDES} -o _oxandale_samuel_hw8_p1_cpp.so

# Integration compilation
integrate:
	g++ oxandale_samuel_hw8_p2.cpp -o oxandale_samuel_hw8_p2


# Remove objects, module files, and exectables
clean:
	@rm -f *.o *.mod *-x *wrap.cxx *.so *.pyc

# Same, but leave the shared objects
neat:
	@rm -f *.o *.mod *-x *wrap.cxx *.pyc

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .f90 .o
.PHONY: clean neat

%.o: %.f90
	$(F90) $(F90FLAGS) $(INCLUDE) -c $<