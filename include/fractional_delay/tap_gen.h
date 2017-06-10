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


#ifndef INCLUDED_FRACTIONAL_DELAY_TAP_GEN_H
#define INCLUDED_FRACTIONAL_DELAY_TAP_GEN_H

#include <fractional_delay/api.h>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <gnuradio/gr_complex.h>
#include <gnuradio/fft/window.h>
#include <gnuradio/fft/fft.h>
#include <boost/math/special_functions/sinc.hpp>

namespace gr {
  namespace fractional_delay {

    class FRACTIONAL_DELAY_API tap_gen {
    public:
      enum method_type {
        SINC_INTERP = 0,
        SINC_SPLINE = 1,
        GLS_APPROX = 2,
        LAGRANGE = 3,
        EQUIRIPPLE = 4,
        AUG_SINC = 2,
        AUG_DIRICHLET = 6,
        AUG_LMS = 7,
      };

      enum win_type {
        WIN_HAMMING = 0,         //!< Hamming window; max attenuation 53 dB
        WIN_HANN = 1,            //!< Hann window; max attenuation 44 dB
        WIN_BLACKMAN = 2,        //!< Blackman window; max attenuation 74 dB
        WIN_RECTANGULAR = 3,     //!< Basic rectangular window; max attenuation 21 dB
        WIN_KAISER = 4,          //!< Kaiser window; max attenuation see window::max_attenuation
        WIN_BLACKMAN_hARRIS = 5, //!< Blackman-harris window; max attenuation 92 dB
        WIN_BLACKMAN_HARRIS = 5, //!< alias to WIN_BLACKMAN_hARRIS for capitalization consistency
        WIN_BARTLETT = 6,        //!< Barlett (triangular) window; max attenuation 26 dB
        WIN_FLATTOP = 7,         //!< flat top window; useful in FFTs; max attenuation 93 dB
      };

      static std::vector<float>
      sinc_interp(double gain, double sampling_freq, double pass_freq,
        int ntaps, double fractional_delay, win_type wt, double beta = 6.76);

      static std::vector<float>
      sinc_spline(double gain, double sampling_freq, double pass_freq,
        double cutoff_freq, int spline_order, int ntaps,
        double fractional_delay);

      static std::vector<float>
      gls_approx(double gain, double sampling_freq, double pass_freq,
        int ntaps, double fractional_delay);

      static std::vector<float>
      lagrange_interp(double gain, int ntaps, double fractional_delay);

      static std::vector<float>
      augment_sinc(double gain,  double sampling_freq, double pass_freq,
                   double fractional_delay, double interp,
                   const std::vector<float> &proto);

      static std::vector<float>
      augment_dirichlet(double gain, double sampling_freq, double pass_freq,
                        double fractional_delay, double interp,
                        const std::vector<float> &proto);

      static std::vector<float>
      augment_lms(double gain, double sampling_freq, double pass_freq,
                  double fractional_delay, double interp,
                  const std::vector<float> &proto, int resolution);


      static std::vector<double>
      get_aug_lms_Pinv_matrix(double sampling_freq, double pass_freq,
                              const std::vector<float> &proto, int resolution);

      static std::vector<double>
      get_aug_lms_frac_delay_p(double sampling_freq, double pass_freq,
                           double fractional_delay, double interp,
                           const std::vector<float> &proto, int resolution);

      static std::vector<float>
      augment_lms_updatable(double gain,
                            const std::vector<float> &proto,
                            const std::vector<double> &Pinv_raw,
                            const std::vector<double> &p_raw);


    };
  }
}
#endif /* INCLUDED_FRACTIONAL_DELAY_TAP_GEN_H */
