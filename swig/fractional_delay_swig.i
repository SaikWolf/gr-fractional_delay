/* -*- c++ -*- */

#define FRACTIONAL_DELAY_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "fractional_delay_swig_doc.i"

%{
#include "fractional_delay/tap_gen.h"
#include "fractional_delay/fd_fft_cc.h"
%}

%include "fractional_delay/tap_gen.h"
%include "fractional_delay/fd_fft_cc.h"
GR_SWIG_BLOCK_MAGIC2(fractional_delay, fd_fft_cc);
