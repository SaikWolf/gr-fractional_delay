#!/usr/bin/python

from . import tap_gen as TG
from gnuradio.fft import window as W
from gnuradio.filter import firdes as FD
from numpy import sinc
from numpy import convolve as conv
from math import sqrt

class tap_augmentor:
  def __init__(self):
    pass

  def __create_def_proto(alpha):
    ntaps = 21
    a = -(ntaps-1)/2 if ntaps%2 else -(ntaps)/2
    b = (ntaps-1)/2 if ntaps%2 else (ntaps)/2-1
    b = b+1
    r = range(a,b)
    taps = [sinc(alpha*x) for x in r]
    return taps

  @staticmethod
  def aug_sinc(gain=1.,rate_in=1.,rate_out=2.,fd=0.,proto=__create_def_proto(0.5)):
    interp = float(rate_out)/rate_in
    sf = 1.
    pf = sf/interp/2.
    taps = TG.augment_sinc(gain,sf,pf,fd,interp,proto)
    return list(taps)

  @staticmethod
  def aug_dirichlet(gain=1.,rate_in=1.,rate_out=2.,fd=0.,proto=__create_def_proto(0.5)):
    interp = float(rate_out)/rate_in
    sf = 1.
    pf = sf/interp/2.
    taps = TG.augment_dirichlet(gain,sf,pf,fd,interp,proto)
    return list(taps)

  @staticmethod
  def aug_lms(gain=1.,rate_in=1.,rate_out=2.,xper=20.,fd=0.,proto=__create_def_proto(0.5),res=2):
    interp = float(rate_out)/rate_in
    sf = 1.
    pf = sf/interp/2.
    xf = pf*(1+xper/100.) if pf*(1+xper/100.)<0.5 else 0.5
    taps = TG.augment_lms(gain,sf,xf,fd,interp,proto,res)
    return list(taps)

  @staticmethod
  def aug_lms_updatable(gain=1.,rate_in=1.,rate_out=2.,xper=20.,fd=0.,proto=__create_def_proto(0.5),res=2,Pinv=None):
    interp = float(rate_out)/rate_in
    sf = 1.
    pf = sf/interp/2.
    xf = pf*(1+xper/100.) if pf*(1+xper/100.)<0.5 else 0.5

    ret_pinv = False

    if(Pinv is None):
      ret_pinv = True
      Pinv = TG.get_aug_lms_Pinv_matrix(sf,xf,proto,res)

    p = TG.get_aug_lms_frac_delay_p(sf, xf, fd, interp, proto, res)

    taps = TG.augment_lms_updatable(gain,proto,Pinv,p)

    if(ret_pinv):
      return (list(taps),list(Pinv))
    else:
      return list(taps)



  @staticmethod
  def aug_lms_dsp(gain=1.,interp_proto=2.,interp_base=1.,fd=0.,proto=__create_def_proto(0.5),res=2,Pinv=None):
    proto_len = len(proto)
    intp_mask = FD.low_pass_2(sqrt(2),2,.5,.1,120,5)
    new_proto = (conv([x for y in [[x,0] for x in proto] for x in y],intp_mask))[0:-1]

    ret_pinv = False
    if(Pinv is None):
      ret_pinv = True
      Pinv = TG.get_aug_lms_Pinv_matrix(1,0.3,new_proto,res)

    p = TG.get_aug_lms_frac_delay_p(1,.3,fd,2*interp_proto*interp_base,new_proto,res)
    fd_taps = TG.augment_lms_updatable(1,new_proto,Pinv,p)
    ext_taps = conv(fd_taps,intp_mask)
    branched_taps = ext_taps[0::2]
    branch_len = len(branched_taps)

    chk1 = (proto_len-1)/2 if proto_len%2 else proto_len/2
    chk2 = (branch_len-1)/2 if branch_len%2 else branch_len/2

    aug_taps = branched_taps[chk2-chk1:chk2+chk1+1]

    tw = sum([x*x for x in aug_taps])
    taps = [x/sqrt(tw) for x in aug_taps]

    if(ret_pinv):
      return (list(taps), list(Pinv))
    else:
      return list(taps)









