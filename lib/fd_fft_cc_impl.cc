/* -*- c++ -*- */
/*
 * Copyright 2017 Bill Clark.
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "fd_fft_cc_impl.h"
#include <stdio.h>

// typdef notes
// complexf == std::complex<float>

namespace gr {
  namespace fractional_delay {

    fd_fft_cc::sptr
    fd_fft_cc::make(float fd, int wt)
    {
      return gnuradio::get_initial_sptr
        (new fd_fft_cc_impl(fd, wt));
    }

    /*
     * The private constructor
     */
    fd_fft_cc_impl::fd_fft_cc_impl(float fd, int wt)
      : gr::sync_block("fd_fft_cc",
              gr::io_signature::make(1, 1, sizeof(complexf)),
              gr::io_signature::make(1, 1, sizeof(complexf))),
        d_fd(fd)
    {
      if((d_fd > 0.5)||(d_fd < -0.5)){
        throw std::runtime_error(
          "fd_fft_cc: fractional delay must be within [-0.5,0.5]\n");
      }
      gen_proto();
      set_window(wt);

      std::vector<float> dummy_taps;
      d_filt = new gr::filter::kernel::fir_filter_ccf(1,dummy_taps);
      install_taps();
    }

    /*
     * Our virtual destructor.
     */
    fd_fft_cc_impl::~fd_fft_cc_impl()
    {
      delete d_filt;
    }

    void
    fd_fft_cc_impl::gen_proto(){
      d_proto = gr::filter::firdes::low_pass_2(1,1,0.5,0.05,61,
                        gr::filter::firdes::WIN_BLACKMAN_hARRIS);
      d_updated = true;
    }

    void
    fd_fft_cc_impl::install_taps(){
      if((d_fd>-1.19e-07)&&(d_fd<1.19e-07)){
        d_taps = std::vector<float>(d_proto.begin(), d_proto.end());
      }
      else{
        d_taps = std::vector<float>(d_proto.size(),0.);
        if(d_window.size()){
          for(int idx = 0; idx < d_proto.size(); idx++){
            for(int ind = 0; ind < d_proto.size(); ind++){
              d_taps[idx] += d_proto[ind]*d_window[ind]*
                  boost::math::sinc_pi(M_PI*(float(idx)-float(ind) - d_fd));
            }
          }
        }
        else{
          for(int idx = 0; idx < d_proto.size(); idx++){
            for(int ind = 0; ind < d_proto.size(); ind++){
              d_taps[idx] += d_proto[ind]*
                  boost::math::sinc_pi(M_PI*(float(idx)-float(ind) - d_fd));
            }
          }
        }
      }
      d_filt->set_taps(d_taps);
      set_history(d_taps.size());
      d_updated = false;
    }

    void
    fd_fft_cc_impl::set_taps(float fd, int wt) {
      d_fd = fd;
      set_window(wt);
    }

    void
    fd_fft_cc_impl::set_fd(float fd) {
      d_fd = fd;
      d_updated = true;
    }

    void
    fd_fft_cc_impl::set_window(int wt) {
      switch(wt){
        case 0:
          d_wt = gr::filter::firdes::WIN_HAMMING;
          break;
        case 1:
          d_wt = gr::filter::firdes::WIN_HANN;
          break;
        case 2:
          d_wt = gr::filter::firdes::WIN_BLACKMAN;
          break;
        case 3:
          d_wt = gr::filter::firdes::WIN_RECTANGULAR;
          break;
        case 5:
          d_wt = gr::filter::firdes::WIN_BLACKMAN_hARRIS;
          break;
        case 6:
          d_wt = gr::filter::firdes::WIN_BARTLETT;
          break;
        case 7:
          d_wt = gr::filter::firdes::WIN_FLATTOP;
          break;
        default:
          d_wt = gr::filter::firdes::WIN_NONE;
          break;
      }
      if(d_wt != gr::filter::firdes::WIN_NONE){
        d_window = gr::filter::firdes::window(d_wt,d_proto.size(),6.76);
      }
      else{
        d_window = std::vector<float>(0);
      }
      d_updated = true;
    }

    int
    fd_fft_cc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      if(d_updated){
        install_taps();
        //proto is constant length
      }

      const complexf *in = (const complexf *) input_items[0];
      complexf *out = (complexf *) output_items[0];

      d_filt->filterN(&out[0], &in[0], noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace fractional_delay */
} /* namespace gr */
