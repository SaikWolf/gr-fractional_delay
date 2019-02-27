/* -*- c++ -*- */

#define FRACTIONAL_DELAY_API

%include "gnuradio.i"			// the common stuff
%include "std_vector.i"

//load generated python docstrings
%include "fractional_delay_swig_doc.i"

%{
#define SWIGPY_SLICE_ARG(obj) ((PySliceObject*) (obj))
#include "fractional_delay/tap_gen.h"
%}

%rename(__assign__) *::operator=;

namespace std{
  %template(FloatVector) vector<float>;
  %template(DoubleVector) vector<double>;
  %template(SizeVector) vector<size_t>;
}

%include "fractional_delay/tap_gen.h"

