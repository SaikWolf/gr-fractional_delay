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
#include <gnuradio/filter/fir_filter.h>
#include <gnuradio/filter/firdes.h>
#include <gnuradio/fft/fft.h>
#include <volk/volk.h>
#include <fftw3.h>
#include <vector>

namespace gr {
  namespace fractional_delay {

    class fd_fft_cc_impl : public fd_fft_cc
    {
     private:
      typedef std::complex<float> complexf;
      size_t d_block_size;
      size_t d_fft_size;
      float  d_fd;//fractional_delay

      //state machine needs
      size_t d_state;
      size_t d_bip;//buff input pointer
      size_t d_bop;//buff_output pointer

      //kernels
      std::vector<float>    d_window;
      std::vector<float>    d_window_inv;
      gr::fft::fft_complex* d_fft;
      gr::fft::fft_complex* d_ifft;
      std::vector<complexf> d_rotator;

      std::vector<gr::filter::kernel::fir_filter_ccf*> d_interp;
      gr::filter::kernel::fir_filter_ccf* d_decim;

      std::vector<float> d_taps_interp;
      std::vector<float> d_taps_decim;

      void install_taps();

      //volk things
      int d_align;

     public:
      fd_fft_cc_impl(size_t block_size, float fd,
                     const std::vector<float> &window);
      ~fd_fft_cc_impl();

      // Where all the action really happens
      int work(int noutput_items,
                gr_vector_const_void_star &input_items,
                gr_vector_void_star &output_items);
    };

  } // namespace fractional_delay
} // namespace gr

#endif /* INCLUDED_FRACTIONAL_DELAY_FD_FFT_CC_IMPL_H */
