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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "fd_fft_cc_impl.h"

// typdef notes
// complexf == std::complex<float>

namespace gr {
  namespace fractional_delay {

    fd_fft_cc::sptr
    fd_fft_cc::make(size_t block_size, float fd,
                    const std::vector<float> &window)
    {
      return gnuradio::get_initial_sptr
        (new fd_fft_cc_impl(block_size, fd, window));
    }

    /*
     * The private constructor
     */
    fd_fft_cc_impl::fd_fft_cc_impl( size_t block_size, float fd,
                                    const std::vector<float> &window)
      : gr::sync_block("fd_fft_cc",
              gr::io_signature::make(1, 1, sizeof(complexf)),
              gr::io_signature::make(1, 1, sizeof(complexf))),
        d_block_size(block_size),
        d_fft_size(4*block_size),
        d_fd(fd),
        d_interp(4),
        d_state(0),
        d_bip(0),
        d_bop(0)
    {
      if((d_fd > 0.5)||(d_fd < -0.5)){
        throw std::runtime_error(
          "fd_fft_cc: fractional delay must be within [-0.5,0.5]\n");
      }
      if(!((window.size() == d_block_size*4)||(window.size()==0))){
        throw std::runtime_error(
          "fd_fft_cc: len(window) is expected to be 4*block_size or empty\n");
      }
      d_window = std::vector<float>(window.begin(), window.end());
      d_window_inv = std::vector<float>(d_window.size());
      for(size_t idx = 0; idx < d_window.size(); idx++){
        d_window_inv[idx] = (std::abs(d_window[idx])>1e-6) ? 1./d_window[idx] : 1.;
      }
      d_fft = new gr::fft::fft_complex(d_fft_size,true);
      d_ifft = new gr::fft::fft_complex(d_fft_size,false);
      d_rotator = std::vector<complexf>(d_fft_size);
      std::vector<float> dummy_taps;
      for(size_t idx = 0; idx < 4; idx++){
        d_interp[idx] = new gr::filter::kernel::fir_filter_ccf(1,dummy_taps);
      }
      install_taps();
      set_output_multiple(d_block_size);

      d_align = volk_get_alignment();
    }

    /*
     * Our virtual destructor.
     */
    fd_fft_cc_impl::~fd_fft_cc_impl()
    {
    }

    void
    fd_fft_cc::install_taps(){
      d_taps_interp = gr::filter::firdes::low_pass_2(2,4,0.5,0.05,61,5);
      int nt = d_taps_interp.size()/4;

      std::vector< std::vector<float> > xtaps(4);
      for(size_t idx = 0; idx < 4; idx++){
        xtaps[idx].resize(nt);
      }

      for(size_t idx = 0; idx < 4; idx++){
        xtaps[idx%4][idx/4] = d_taps_interp[idx];
      }
      for(size_t idx = 0; idx < 4; idx++){
        d_interp[idx]->set_taps(xtaps[idx]);
      }

      d_taps_decim = std::vector<float>(d_taps)
      d_decim = new gr::filter::fir_filter_ccf(4,d_taps_decim);

      for(size_t idx = 0; idx < d_fft_size; idx++){
        d_rotator[idx] = std::exp(complexf(0.,
          2*M_PI*float(idx)*(4*d_fd)/float(d_fft_size)));
      }

      set_history(nt);
    }

    int
    fd_fft_cc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const complexf *in = (const copmlexf *) input_items[0];
      complexf *out = (complexf *) output_items[0];

      complexf* fft_in = d_fft->get_inbuf();
      complexf* fft_out = d_fft->get_outbuf();
      complexf* ifft_in = d_ifft->get_inbuf();
      complexf* ifft_out = d_ifft->get_outbuf();

      ////// noutput_items % d_block_size == 0 ? TRUE

      size_t oo(0),ii(0);
      while(oo < noutput_items){
        if(d_state == 0){
          while(d_bip < d_fft_size){
            for(size_t idx = 0; idx < 4; idx++){
              fft_in[d_bip+idx] = d_interp[idx]->filter(&in[ii]);
            }
            ii++;
            d_bip+=4;
          }
          d_state = 1;
          d_bip = 0;
        }
        if(d_state == 1){
          if(d_window.size()){
            volk_32fc_32f_multiply_32fc(&fft_in[0], &fft_in[0],
                                        &d_window[0], d_fft_size);
          }
          d_fft->execute();
          d_state = 2;
        }
        if(d_state == 2){
          volk_32fc_x2_multiply_32fc(&ifft_in[0], &fft_out[0],
                                    &d_rotator[0], d_fft_size);
          d_state = 3;
        }
        if(d_state == 3){
          d_ifft->execute();
          if(d_window.size()){
            volk_32fc_32f_multiply_32fc(&ifft_out[0], &ifft_out[0],
                                        &d_window_inv[0], d_fft_size)
          }
          d_state = 4;
        }
        if(d_state == 4){
          d_decim->filterNdec(&out[oo], &ifft_out[0], d_block_size);
          oo += d_block_size;
          d_state = 0;
        }
      }

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace fractional_delay */
} /* namespace gr */
