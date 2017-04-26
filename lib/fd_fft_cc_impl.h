/* -*- c++ -*- */
/*
 * Copyright 2017 <+YOU OR YOUR COMPANY+>.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_FRACTIONAL_DELAY_FD_FFT_CC_IMPL_H
#define INCLUDED_FRACTIONAL_DELAY_FD_FFT_CC_IMPL_H

#include <fractional_delay/fd_fft_cc.h>
#include <boost/math/special_functions/sinc.hpp>
#include <gnuradio/filter/fir_filter.h>
#include <vector>

namespace gr {
  namespace fractional_delay {

    class fd_fft_cc_impl : public fd_fft_cc
    {
     private:
      typedef std::complex<float> complexf;
      float  d_fd;//fractional_delay
      gr::filter::firdes::win_type d_wt;
      bool d_updated;

      //kernels
      std::vector<float>    d_window;

      std::vector<float> d_proto;
      std::vector<float> d_taps;
      gr::filter::kernel::fir_filter_ccf* d_filt;

      void install_taps();
      void gen_proto();

     public:
      fd_fft_cc_impl(float fd, int wt);
      ~fd_fft_cc_impl();

      // Where all the action really happens
      int work(int noutput_items,
                gr_vector_const_void_star &input_items,
                gr_vector_void_star &output_items);

      void set_taps(float fd, int wt);
      void set_fd(float fd);
      void set_window(int wt);
    };

  } // namespace fractional_delay
} // namespace gr

#endif /* INCLUDED_FRACTIONAL_DELAY_FD_FFT_CC_IMPL_H */
