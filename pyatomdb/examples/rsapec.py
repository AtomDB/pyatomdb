# import relevant models
# inital import creates rsapec, rsvapec, rsvvapec, analagous to apec, vapec, vvapec but resonance scattered

from xspec import *
import rsapec_xspec
import xspec
import matplotlib.pyplot as plt


#import the response and spectrum file of Hitomi observation of Perseus core
s = Spectrum('ah100040030sxs_p0px1010_cl2_HP.pi')
s.response = 'ah100040030_sxs_cl2_HP_small.rmf'
s.response.arf = 'ah100040030sxs_p0px1010_ptsrc.arf'

#set the energy range for the fitting
s.ignore('**-6.0,7.0-**')


# declare a new model and initiate the model parameters
m1 = xspec.Model('tbabs*pyapecrs')
m1.TBabs.nH=0.138
m1.TBabs.nH.frozen=True
m1.pyapecrs.kT=4.06
m1.pyapecrs.nL=8.05E+21
m1.pyapecrs.Velocity=175
m1.pyapecrs.Abundanc=0.39
m1.pyapecrs.Redshift=1.72391E-02
m1.pyapecrs.Redshift.frozen=True
m1.pyapecrs.norm=0.97
Fit.statMethod = "cstat"


# let's overplot the model and the spectrum 

Plot.device='/null'
Plot.xAxis = "keV"
Plot.setRebin(5,10)
Plot('data')
chans = Plot.x()
rates = Plot.y()
xErrs = Plot.xErr()
yErrs = Plot.yErr()
folded= xspec.Plot.model(1)
plt.xlim(6.5,6.6)
plt.plot(chans,folded, color='orange', linewidth=2, label='rsapec model')
plt.errorbar(chans, rates, xerr=xErrs,yerr=yErrs, label='Hitomi observation', fmt=",", color='blue')
plt.xlabel('Energy(keV)')
plt.ylabel('counts/s/keV')
plt.legend()
plt.show()

# save image files
plt.savefig('Rsapec_hitomi.pdf')
plt.savefig('Rsapec_hitomi.svg')
