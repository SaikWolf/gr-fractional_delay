#!/bin/usr/python

from . import tap_generator as TG
from . import tap_augmentor as TA

class fd_struct:
  def __init__(self):
    self.pinv = None
    self.mode = 0
    self.gain = 1.
    self.rate_in = 1.
    self.rate_out = 2.
    self.xper = 20.
    self.ntaps = 21
    self.fd = 0.
    self.order = 3
    self.win = 0
    self.proto = [1.,]
    self.res = 2
    self.beta = 6.76

  def get_taps(self,fd=None):
    if(fd is not None):
      self.update(fd=fd)

    if(self.mode == 0):
      taps = TG.gen_sinc_frac_delay(self.ntaps,self.fd)
    elif(self.mode == 1):
      taps = TG.gen_sinc_interp(self.gain, self.rate_in, self.rate_out, self.ntaps, self.fd, self.win, self.beta)
    elif(self.mode == 2):
      taps = TG.gen_sinc_spine(self.gain, self.rate_in, self.rate_out, self.xper, self.order, self.ntaps, slef.fd)
    elif(self.mode == 3):
      taps = TG.gen_lagrange_interp(self.gain, self.ntaps, self.fd)
    elif(self.mode == 4):
      taps = TG.gen_gls_frac_delay(self.gain,self.xper,self.ntaps,self.fd)
    elif(self.mode == 5):
      if(self.pinv is None):
        (taps,self.pinv) = TG.gen_gls_updatable(self.gain,self.xper,self.ntaps,self.fd,self.pinv)
      else:
        taps = TG.gen_gls_updatable(self.gain,self.xper,self.ntaps,self.fd,self.pinv)
    elif(self.mode == 6):
      taps = TA.aug_sinc(self.gain,self.rate_in,self.rate_out,self.fd,self.proto)
    elif(self.mode == 7):
      taps = TA.aug_dirichlet(self.gain,self.rate_in,self.rate_out,self.fd,self.proto)
    elif(self.mode == 8):
      taps = TA.aug_lms(self.gain,self.rate_in,self.rate_out,self.xper,self.fd,self.proto,self.res)
    elif(self.mode == 9):
        if(self.pinv is None):
          (taps,self.pinv) = TA.aug_lms_updatable(self.gain,self.rate_in,self.rate_out,self.fd,self.proto,self.res,self.pinv)
        else:
          taps = TA.aug_lms_updatable(self.gain,self.rate_in,self.rate_out,self.fd,self.proto,self.res,self.pinv)
    else:
      taps = [-1.]
    return taps

  def update(self, mode=None, gain=None, rate_in=None, rate_out=None, xper=None, ntaps=None, fd=None, order=None, win=None, proto=None, res=None, beta=None):
    self.mode     = mode      if mode     is not None else self.mode
    self.gain     = gain      if gain     is not None else self.gain
    self.rate_in  = rate_in   if rate_in  is not None else self.rate_in
    self.rate_out = rate_out  if rate_out is not None else self.rate_out
    self.xper     = xper      if xper     is not None else self.xper
    self.ntaps    = ntaps     if ntaps    is not None else self.ntaps
    self.fd       = fd        if fd       is not None else self.fd
    self.order    = order     if order    is not None else self.order
    self.win      = win       if win      is not None else self.win
    self.proto    = proto     if proto    is not None else self.proto
    self.res      = res       if res      is not None else self.res
    self.beta     = beta      if beta     is not None else self.beta
    if((proto is not None) or (res is not None) or (rate_in is not None) or (rate_out is not None) or (ntaps is not None)):
      self.pinv = None





































