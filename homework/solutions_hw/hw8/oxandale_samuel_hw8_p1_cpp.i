// The syntax is %module module_name, i.e., the name
// you want to import.  All of the C++ includes you need
// should be included here.
%module oxandale_samuel_hw8_p1_cpp
%{
#include "oxandale_samuel_hw8_p1_cpp.hh"
    %}


// this stuff is boiler plate:
%{
#define SWIG_FILE_WITH_INIT
    %}
%include "numpy.i"
%init %{
    import_array();
    %}

// The following are examples of a SWIG "typemap".  This is an advanced
// feature that I've only ever tinkered with seriously for one project, but
// I've **used** the labor of others liberally.

// This is what we want for input arrays that should not get changed
%apply (double* IN_ARRAY1, int DIM1) {(double* vec1, int len1),
    (double* vec2, int len2)}
%rename (interpolate) my_interpolate;
%exception interpolate {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%inline %{
    double interpolate(double x_new, double* vec1, int len1, double* vec2, int len2, int order) {
        if (len1 != len2) {
            PyErr_Format(PyExc_ValueError,
                         "Arrays of lengths (%d,%d) given",
                         len1, len2);
            return 0.0;
        }
        return interpolate(x_new, vec1, vec2, len1, order);
    }
    %}

// This is what we want for an array produced by the function.
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* a, int n)}

// This is what we want when we need the raw data already stored in some
// C++ object (e.g., a class, etc.).  BE CAREFUL.
//%apply (double** ARGOUTVIEW_ARRAY1, int *DIM1) {(double** a, int *n)}

// INCLUDE THE HEADER FILE WITH THE ACTUAL DECLARATIONS
%include "oxandale_samuel_hw8_p1_cpp.hh"

