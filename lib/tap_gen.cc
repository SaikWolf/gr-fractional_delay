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

 #include <fractional_delay/tap_gen.h>
 #include <iostream>

 namespace gr{
   namespace fractional_delay{

     typedef gr_complex complexf;
     typedef gr_complexd complexd;

    std::vector<float>
    tap_gen::sinc_interp(double gain, double sampling_freq, double pass_freq,
                         int ntaps, double fractional_delay,
                         win_type wt, double beta)
    {
      double alpha = pass_freq/sampling_freq;
      if(alpha > .5) alpha = .5;
      gr::fft::window::win_type twt = static_cast<gr::fft::window::win_type>(wt);
      std::vector<float> window = gr::fft::window::build(twt, ntaps, beta);
      std::vector<float> taps(ntaps);
      double center;
      double fd = std::fmod(fractional_delay,1.0);
      if(ntaps%2){
        //odd
        center = double(ntaps-1)/2.;
      }
      else{
        //even
        center = double(ntaps)/2.;
      }
      if(fractional_delay < 0.) center -= 1.;

      double m2 = 0.;
      for(size_t n = 0; n < ntaps; n++){
        taps[n] = window[n]*2.*alpha*
                  boost::math::sinc_pi(2.*M_PI*alpha*(double(n-(center+fd))));
        m2 += taps[n]*taps[n];
      }
      //double scale = std::sqrt(gain*ntaps/m2);
      double scale = std::sqrt(gain/m2);
      for(size_t n = 0; n < ntaps; n++){
        taps[n] *= scale;
      }
      return taps;
    }

    std::vector<float>
    tap_gen::sinc_spline(double gain, double sampling_freq, double pass_freq,
                         double cutoff_freq, int spline_order, int ntaps,
                         double fractional_delay)
    {
      double fp = pass_freq/sampling_freq;
      double fc = cutoff_freq/sampling_freq;
      if(fp > 0.5) fp = 0.5;
      if(fc < fp) fc = fp;
      double df = fc-fp;
      double sf = fc+fp;
      double center;
      double fd = std::fmod(fractional_delay,1.0);
      if(ntaps%2){
        //odd
        center = double(ntaps-1)/2.;
      }
      else{
        //even
        center = double(ntaps)/2.;
      }
      if(fractional_delay < 0) center -= 1.;
      double P(spline_order);

      std::vector<float> taps(ntaps);
      double m2 = 0.;
      double sc = 0.;
      for(size_t n = 0; n < ntaps; n++){
        taps[n] = boost::math::sinc_pi(M_PI*df*(double(n)-(center+fd))/(2.*P));
        sc = double(n)-(center+fd);
        if(std::abs(sc) < 3.0e-16) sc = 2*P*sf;//sin(x) approx x
        else sc = std::sin(M_PI*sc*sf)/(M_PI*sc/(2*P));
        taps[n] = std::pow(taps[n],P)*sc;
        m2 += taps[n]*taps[n];
      }
      //double scale = std::sqrt(gain*ntaps/m2);
      double scale = std::sqrt(gain/m2);
      for(size_t n = 0; n < ntaps; n++){
        taps[n] *= scale;
      }
      return taps;
    }

    std::vector<float>
    tap_gen::gls_approx(double gain, double sampling_freq, double pass_freq,
                        int ntaps, double fractional_delay)
    {
      double alpha = pass_freq/sampling_freq;
      if(alpha > 0.5) alpha = 0.5;

      double center;
      double fd = std::fmod(fractional_delay,1.);
      if(ntaps%2){
        //odd
        center = double(ntaps-1)/2.;
      }
      else{
        //even
        center = double(ntaps)/2.;
      }
      if(fractional_delay < 0.) center -= 1.;

      std::vector<double> P_raw(ntaps*ntaps);
      std::vector<double> p_raw(ntaps);

      for(size_t k = 0; k < ntaps; k++){
        p_raw[k] = 2.*alpha*boost::math::sinc_pi(2.*M_PI*alpha*(double(k)-(center+fd)));
        for(size_t l = 0; l < ntaps; l++){
          P_raw[k*ntaps+l] = 2.*alpha*
              boost::math::sinc_pi(2.*M_PI*alpha*double(k-l));
        }
      }

      Eigen::Map< Eigen::MatrixXd > mapP( &P_raw[0], ntaps, ntaps );
      Eigen::Map< Eigen::VectorXd > mapp( &p_raw[0], ntaps );
      Eigen::MatrixXd P = mapP;
      Eigen::VectorXd p = mapp;

      Eigen::VectorXd x = P.colPivHouseholderQr().solve(p);

      std::vector<double> x_raw( &x(0), &x(0) + ntaps );
      std::vector<float> taps( x_raw.begin(), x_raw.end() );
      double m2 = 0.;
      for(size_t n = 0; n < ntaps; n++){
        m2 += taps[n]*taps[n];
      }
      //double scale = std::sqrt(gain*ntaps/m2);
      double scale = std::sqrt(gain/m2);
      for(size_t n = 0; n < ntaps; n++){
        taps[n] *= scale;
      }
      return taps;
    }

    std::vector<float>
    tap_gen::lagrange_interp(double gain, int ntaps, double fractional_delay)
    {
      int N = ntaps-1;

      double center;
      double fd = std::fmod(fractional_delay,1.0);
      if(N%2){//odd
        center = double(ntaps-1)/2.;
      }
      else{//even
        center = double(ntaps)/2.;
      }
      if(fractional_delay < 0.) center -= 1;

      double D = center + fd;
      std::vector<float> taps(ntaps,1.);

      for(int n = 0; n < ntaps; n++){
        for(int k = 0; k < ntaps; k++){
          taps[n] = (n==k) ? taps[n] : taps[n]*double(D-k)/double(n-k);
        }
      }
      return taps;
    }

    std::vector<float>
    tap_gen::augment_sinc(double gain, double sampling_freq, double pass_freq,
                          double fractional_delay, double interp,
                          const std::vector<float> &proto)
    {
      int ntaps = proto.size();
      int N = ntaps - 1;
      if(N){
        double alpha = pass_freq/sampling_freq;
        double center;
        double fd = std::fmod(fractional_delay,1.0)*interp;
        if(fractional_delay < 0.) fd -= interp;

        std::vector<float> prototaps( proto.begin(), proto.end() );
        double p2(0.);
        for(size_t n = 0; n < ntaps; n++){
          p2 += prototaps[n]*prototaps[n];
        }
        double scaleP = std::sqrt(1./p2);
        for(size_t n = 0; n < ntaps; n++){
          prototaps[n] *= scaleP;
        }

        std::vector<float> taps(ntaps,0.);
        for(int n = 0; n < ntaps; n++){
          for(int k = 0; k < ntaps; k++){
            taps[n] += 2.*alpha*prototaps[k]*
              boost::math::sinc_pi(2.*alpha*M_PI*(double(n-k)-fd));
          }
        }

        double t2 = 0.;
        for(size_t idx = 0; idx < taps.size(); idx++){
          t2 += taps[idx]*taps[idx];
        }
        double scaleT = std::sqrt(1./t2);
        for(size_t idx = 0; idx < taps.size(); idx++){
          taps[idx] *= scaleT;
        }

        return taps;
      }
      else{
        std::vector<float> taps(0);
        return taps;
      }
    }

    std::vector<float>
    tap_gen::augment_dirichlet(double gain, double sampling_freq, double pass_freq,
                               double fractional_delay, double interp,
                               const std::vector<float> &proto)
    {
      int ntaps = proto.size();
      int N = ntaps - 1;
      if(N){
        double alpha = pass_freq/sampling_freq;
        double fd = std::fmod(fractional_delay,1.0)*interp;
        if(fractional_delay < 0.) fd -= interp;

        std::vector<float> prototaps( proto.begin(), proto.end() );
        double p2(0.);
        for(size_t n = 0; n < ntaps; n++){
          p2 += prototaps[n]*prototaps[n];
        }
        double scaleP = std::sqrt(1./p2);
        for(size_t n = 0; n < ntaps; n++){
          prototaps[n] *= scaleP;
        }

        std::vector<float> taps(ntaps,0.);
        double denom, weight;
        for(int n = 0; n < ntaps; n++){
          for(int k = 0; k < ntaps; k++){
            denom = ntaps*std::sin(2.*alpha*M_PI*(double(n-k)-fd)/double(ntaps));
            if((denom < 2.22e-16)&&(denom > -2.22e-16)) weight = 1.;
            else weight = std::sin(2.*alpha*M_PI*(double(n-k)-fd))/denom;
            taps[n] += 2.*alpha*prototaps[k]*weight;
          }
        }

        double t2 = 0.;
        for(size_t idx = 0; idx < taps.size(); idx++){
          t2 += taps[idx]*taps[idx];
        }
        double scaleT = std::sqrt(1./t2);
        for(size_t idx = 0; idx < taps.size(); idx++){
          taps[idx] *= scaleT;
        }

        return taps;
      }
      else{
        std::vector<float> taps(0);
        return taps;
      }
    }

    std::vector<float>
    tap_gen::augment_lms(double gain, double sampling_freq, double pass_freq,
                         double fractional_delay, double interp,
                         const std::vector<float> &proto, int resolution)
    {
      int ntaps = proto.size();
      int N = ntaps-1;
      if(N){
        double alpha = pass_freq/sampling_freq;
        double center;
        double fd = std::fmod(fractional_delay,1.0)*interp;
        if(N%2){//odd
          center = double(ntaps-1)/2.;
        }
        else{//even
          center = double(ntaps)/2.;
        }
        if(fractional_delay < 0.) center -= interp;

        double D = center + fd;


        std::vector<float> prototaps( proto.begin(), proto.end() );
        double p2(0.);
        for(size_t n = 0; n < ntaps; n++){
          p2 += prototaps[n]*prototaps[n];
        }
        double scaleP = std::sqrt(1./p2);
        for(size_t n = 0; n < ntaps; n++){
          prototaps[n] *= scaleP;
        }

        int nbins = resolution*ntaps;
        gr::fft::fft_complex to_spec(nbins, true);

        complexf *ib1 = to_spec.get_inbuf();
        complexf *ob1 = to_spec.get_outbuf();

        memset( ib1, 0, nbins*sizeof(complexf) );
        memset( ob1, 0, nbins*sizeof(complexf) );

        std::vector<double> mag_spec0(nbins);
        std::vector<double> mag_specA(nbins);

        for(size_t idx = 0; idx < ntaps; idx++){
          ib1[idx] = complexf(prototaps[idx],0.);
        }
        to_spec.execute();

        for(size_t idx = 0; idx < nbins; idx++){
          mag_spec0[idx] = double((ob1[idx]*std::conj(ob1[idx])).real());
        }

        ///////////fft shift needed
        size_t pointA = int(std::ceil(float(nbins)/2.));
        memcpy( &mag_specA[0], &mag_spec0[pointA], (nbins-pointA)*sizeof(double) );
        memcpy( &mag_specA[pointA], &mag_spec0[0], (pointA)*sizeof(double) );

        double incrmt = (2*M_PI - 1./double(nbins))/double(nbins-1);
        std::vector<double> omegaA(nbins,0.);
        for(size_t idx = 0; idx < nbins; idx++){
          omegaA[idx] = -M_PI + double(idx)*(incrmt);
        }

        int start_at = int(std::max(0.,
                      round(double(nbins)*(1.-2*alpha)/2.)));
        int end_at = int(std::min(double(nbins),
                      round(double(nbins)*(1.-(1.-2*alpha)/2.))));//one shy

        std::vector<double> omegaB( omegaA.begin()+start_at, omegaA.begin()+end_at );
        std::vector<double> mag_specB( mag_specA.begin()+start_at, mag_specA.begin()+end_at );

        int chk = omegaB.size();
        std::vector<complexd> Hid(chk);
        std::vector<double> c(chk);
        std::vector<double> s(chk);
        std::vector<complexd> ek(chk);
        std::vector<complexd> el(chk);
        for(size_t idx = 0; idx < chk; idx++){
          Hid[idx] = std::exp(complexd(0.,-omegaB[idx]*D));
        }

        std::vector<double> p_raw(ntaps);
        std::vector<double> P_raw(ntaps*ntaps);

        double Ckl(0.),pk(0.),scale_w(1./double(chk));
        for(int k = 0; k < ntaps; k++){
          pk = 0.;
          for(size_t w = 0; w < chk; w++){
            ek[w] = std::exp(complexd(0.,-omegaB[w]*k));
            c[w] = std::cos(k*omegaB[w]);
            s[w] = std::sin(k*omegaB[w]);
            pk += mag_specB[w]*(Hid[w].real()*c[w]-Hid[w].imag()*s[w]);
          }
          p_raw[k] = pk*scale_w;
          for(int l = 0; l < ntaps; l++){
            Ckl = 0.;
            for(size_t w = 0; w < chk; w++){
              el[w] = std::exp(complexd(0.,omegaB[w]*l));
              Ckl += mag_specB[w]*(ek[w]*el[w]).real();
            }
            P_raw[k*ntaps+l] = Ckl*scale_w;
          }
        }

        Eigen::Map< Eigen::MatrixXd > mapP( &P_raw[0], ntaps, ntaps );
        Eigen::Map< Eigen::VectorXd > mapp( &p_raw[0], ntaps );
        Eigen::MatrixXd P = mapP;
        Eigen::VectorXd p = mapp;

        Eigen::VectorXd x = P.colPivHouseholderQr().solve(p);

        std::vector<double> x_raw( x.data(), x.data()+x.size() );
        std::vector<float> augmentor( x_raw.begin(), x_raw.end() );
        std::vector<float> taps( 2*ntaps-1, 0. );
        double m2(0.);
        for(size_t n = 0; n < ntaps; n++){
          m2 += augmentor[n]*augmentor[n];
        }

        double scaleM = std::sqrt(1./m2);
        for(size_t n = 0; n < ntaps; n++){
          augmentor[n] *= scaleM;
        }
        double t2(0.);
        for(int idx = 0; idx < taps.size(); idx++){
          if(idx < ntaps){
            for(int m = 0; m <= idx; m++){
              taps[idx] += prototaps[m] * augmentor[idx-m];
            }
          }
          else{
            for(int m = idx-(ntaps-1); m < ntaps; m++){
              taps[idx] += prototaps[m] * augmentor[idx-m];
            }
          }
          t2 += taps[idx]*taps[idx];
        }
        double scaleT = std::sqrt(gain/t2);
        for(size_t idx = 0; idx < taps.size(); idx++){
          taps[idx] *= scaleT;
        }
        return taps;
      }
      else{
        std::vector<float> taps(0);
        return taps;
      }
    }

   }
 }
