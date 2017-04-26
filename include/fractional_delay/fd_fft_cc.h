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


#ifndef INCLUDED_FRACTIONAL_DELAY_FD_FFT_CC_H
#define INCLUDED_FRACTIONAL_DELAY_FD_FFT_CC_H

#include <fractional_delay/api.h>
#include <gnuradio/sync_block.h>
#include <gnuradio/filter/firdes.h>

namespace gr {
  namespace fractional_delay {

    /*!
     * \brief Fractional Delay FIR filter
     * Utilizes the FFT/Sinc interpolation approach.
     * For now, the fractional delay comes at the
     * cost of a 27 sample delay as well, purely for
     * implementation ease. It is within reason for the
     * the delay to >= 5 samples (worse performance as
     * delay approaches 5 samples however).
     * \ingroup fractional_delay
     *
     */
    class FRACTIONAL_DELAY_API fd_fft_cc : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<fd_fft_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of fractional_delay::fd_fft_cc.
       *
       * To avoid accidental use of raw pointers, fractional_delay::fd_fft_cc's
       * constructor is in a private implementation
       * class. fractional_delay::fd_fft_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(float fd, int wt);

      /*!
       * \brief Update both the fractional_delay and window type.
       *
       * Update both the fractional_delay and window type.
       * (Works on the fly since the prototype filter doesn't change)
       */
      virtual void set_taps(float fd, int wt)=0;

      /*!
       * \brief Update the fractional_delay.
       *
       * Update the fractional_delay.
       * (Works on the fly since the prototype filter doesn't change)
       */
      virtual void set_fd(float fd)=0;

      /*!
       * \brief Update the window type.
       *
       * Update window type.
       * (Works on the fly since the prototype filter doesn't change)
       * (Little net effect in current implementation has been seen)
       */
      virtual void set_window(int wt)=0;
    };

  } // namespace fractional_delay
} // namespace gr

#endif /* INCLUDED_FRACTIONAL_DELAY_FD_FFT_CC_H */
