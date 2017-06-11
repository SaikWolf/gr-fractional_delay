#!/usr/bin/python

from . import tap_gen as TG
from gnuradio.fft import window as W
from numpy import sinc

class tap_generator:
  def __init__(self):
    pass

  @staticmethod
  def gen_sinc_frac_delay(ntaps=21,fd=0.):
    a = -(ntaps-1)/2 if ntaps%2 else -(ntaps)/2
    b = (ntaps-1)/2 if ntaps%2 else (ntaps)/2-1
    r = range(a,b)
    taps = [sinc(x-fd) for x in r]
    return taps

  @staticmethod
  def gen_sinc_interp(gain=1., rate_in=1., rate_out=2.0, ntaps=21, fd=0., win=W.WIN_RECTANGULAR, beta=6.76):
    interp = float(rate_out)/rate_in
    sf = 1.
    pf = sf/interp/2.
    taps = TG.sinc_interp(gain,sf,pf,ntaps,fd,win,beta)
    return list(taps)

  @staticmethod
  def gen_sinc_spline(gain=1., rate_in=1.,rate_out=2.0,xper=20.,order=3,ntaps=21,fd=0.):
    interp = float(rate_out)/rate_in
    sf = 1.
    pf = sf/interp/2.
    xf = pf*(1.+xper/100.) if pf*(1.+xper/100.)<0.5 else 0.5
    print sf,pf,xf
    taps = TG.sinc_spline(gain,sf,pf,xf,order,ntaps,fd)
    return list(taps)

  @staticmethod
  def gen_lagrange_interp(gain=1.,ntaps=21,fd=0.):
    taps = TG.lagrange_interp(gain,ntaps,fd)
    return taps

  @staticmethod
  def gen_gls_frac_delay(gain=1.,xper=30.,ntaps=21,fd=0.):
    taps = TG.gls_approx(gain,1.,xper/100.,ntaps,fd)
    return list(taps)

  @staticmethod
  def gen_gls_updatable(gain=1.,xper=30.,ntaps=21,fd=0.,Pinv=None):
    ret_pinv = False
    if(Pinv is None):
      ret_pinv = True
      Pinv = TG.get_gls_Pinv_matrix(1.,xper/100.,ntaps)

    p = TG.get_gls_frac_delay_p(1.,xper/100.,ntaps,fd)
    taps = TG.gls_updatable(gain,ntaps,Pinv,p)
    if(ret_pinv):
      return (list(taps),list(Pinv))
    else:
      return list(taps)





