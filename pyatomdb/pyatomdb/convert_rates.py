import sys
import read_spex_urdam
import numpy
import scipy.special
import pylab
import pyatomdb






mec2 = 510.99895000



g7a1 = 3.480230906913262     # sqrt(pi)*(gamma + ln(4))
g7a2 = 1.772453850905516     # sqrt(pi)
g7a3 = 1.360544217687075E-02 # 2/147
g7a4 = 0.4444444444444444    # 4/9
g7c1 = -4881/8.
g7c2 =  1689/16.
p7 = numpy.array([1.000224,-0.113011,1.851039,0.0197311,0.921832,2.651957])
gamma = 0.577215664901532860606512    #Euler's constant
a = numpy.array([0.999610841,3.50020361,-0.247885719,0.0100539168,1.39075390e-3,1.84193516,4.64044905])
r0 = 4.783995473666830e-10   #=2*sqrt(2/pi)*c*1e-18



def eval_udi_shell(coeff, kT):
  eion = coeff['I_keV']
  y = eion/kT
  lam = eion/mec2
  en1 = numpy.exp(y) * scipy.special.exp1( y)
  
  g=numpy.zeros((8, len(kT)))
  g[0,:] = 1/y
  g[1,:] = en1
  g[2,:] = scipy.special.expn(2, y)*numpy.exp(y)
  g[3,:] = en1/y
  g[4,:] = (1+en1)/y**2
  g[5,:] = (3+y+2*en1)/y**3
  
  
  k = numpy.where(y<0.6)[0]
  if len(k) > 0:
    yy = y[k]
    g[6,k] = numpy.exp(yy) * \
          ((((yy / 486.0 - g7a3) * yy + 0.08) * yy - g7a4) * yy + 4 \
	        - (g7a2 * numpy.log(yy) + g7a1) / numpy.sqrt(yy))
  
  k = numpy.where((y>=0.6) & (y<=20.0))[0]
  if len(k) > 0:
    yy=y[k]
    g[6,k] = (p7[0] + (p7[1] + p7[3] * numpy.log(yy)) /\
              numpy.sqrt(yy) + p7[2] / yy) / (yy + p7[4]) / (yy + p7[5])
              
  k = numpy.where(y>20.0)[0]
  if len(k) > 0:
    yy=y[k]
    g[6,k] = (((((g7c1 / yy + g7c2) / yy - 22.0) /\
                 yy + 5.75) / yy - 2) / yy + 1)/yy**2


  k = numpy.where(y<0.5)[0]
  if len(k) > 0:
    yy = y[k]
    g[7,k] = (((((-yy / 3000.0 - (1.0 / 384.0)) * yy - (1.0 / 54.0)) * yy +\
             0.125) * yy - 1.) * yy + 0.989056 + \
              (numpy.log(yy) / 2.0 + gamma) * numpy.log(yy)) * numpy.exp(yy)
  

  k = numpy.where(((y>=0.5) & (y<=20.0)))[0]
  if len(k) > 0:
    yy = y[k]
#    print(yy)
    g[7,k] = ((((a[4] / yy + a[3]) / yy + a[2]) / yy + a[1]) / yy + a[0]) /\
               (yy + a[5]) / (yy + a[6])
  k = numpy.where(y>20.0)[0]
  if len(k) > 0:
    yy = y[k]
    g[7,k] = ((((((13068 / yy - 1764) / yy + 274) / yy - 50) / yy +\
               11) / yy - 3) / yy + 1) / yy**2

  
  shlout =  coeff['A'] * 1e24*(g[0,:] - g[1,:]) + \
	    coeff['B'] * 1e24* (g[0,:] - 2 * g[1,:] + g[2,:]) + \
	    coeff['C'] * 1e24* (g[3,:] + 1.5 * lam * g[4,:] + 0.25 * lam**2 * g[5,:]) + \
	    coeff['D'] * 1e24* g[6,:] +\
	    coeff['E'] * 1e24* g[7,:]

  shlout *= r0 * numpy.exp(-y) / eion**2 * y**1.5 * numpy.sqrt(lam)
  
  return(shlout)

def eval_udi(coeff, kT):
  rate = numpy.zeros(len(kT))
  for i in range(len(coeff)):
#    print(coeff[i])
    rate += eval_udi_shell(coeff[i], kT)

  return(rate)
    
def eval_uea_shell(coeff, kT, altmethod=False):
  eion = coeff['I_keV']
  y = eion/kT
  #print(eion)
  #print(kT)
  #print(y)
#  zzz=input()
  lam = eion/mec2
  exp1=scipy.special.exp1( y)
  en1 = numpy.exp(y) * exp1
  eminusy = numpy.exp(-y)
  
#dtype=[('I_keV', '<f8'), ('A', '<f8'), ('B', '<f8'), ('C', '<f8'), ('D', '<f8'), ('E', '<f8')])
  m1 = (1/y)* eminusy
  m2 = exp1
  m3 = eminusy- y*exp1
  m4 = (1-y)*eminusy/2 + y**2*exp1/2
  m5 = exp1/y
  
  c_ea = m1 * coeff['A']* 1e24 +\
         m2 * coeff['B']* 1e24 +\
         m3 * coeff['C']* 1e24 +\
         m4 * 2 * coeff['D']* 1e24 +\
         m5 * coeff['E']* 1e24
         
  c_ea *=eminusy*y**1.5*numpy.sqrt(lam)/eion**2
  
  if altmethod:
    #print("ALTMETHOD")
  
#    m1 =  (1/y)* eminusy
#    m2 = scipy.special.expn(1, y)
#    m3 = scipy.special.expn(2, y)
#    m4 = scipy.special.expn(3, y)
#    m5 = scipy.special.expn(1, y)/y
#    c_ea = m1 * coeff['A']* 1e24 +\
#           m2 * coeff['B']* 1e24 +\
#           m3 * coeff['C']* 1e24 +\
#           m4 * 2 * coeff['D']* 1e24 +\
#           m5 * coeff['E']* 1e24


    em1 = scipy.special.expn(1, y)*numpy.exp(y)
    em2 = scipy.special.expn(2, y)*numpy.exp(y)
    em3 = scipy.special.expn(3, y)*numpy.exp(y)
    c_ea = coeff['A']*1e24 + y * (coeff['B']*1e24*em1 + coeff['C']*1e24 * em2 + 2*coeff['D']*1e24 * em3) + coeff['E']*1e24*em1
    
    
#    c_ea *=eminusy*y**1.5*numpy.sqrt(lam)/eion**2
#    print(r0)
#    r0tmp = 2*numpy.sqrt(2/numpy.pi) * kT**-1.5 * mec2**-0.5
#    c_ea *= numpy.exp(-y)
#    print(r0)
    c_ea *=  numpy.exp(-y)/eion**2 * y**1.5 * lam**0.5
  else:
    #print("NORMMETHOD")
    r0tmp=r0*1

  
  return(r0*c_ea)

def eval_uea(coeff, kT, altmethod=False):
  rate = numpy.zeros(len(kT))
  for i in range(len(coeff)):
#    print(coeff[i])
    rate += eval_uea_shell(coeff[i], kT, altmethod=altmethod)

  return(rate)
    
#print("UDICOEFF", udicoeff)
#print("UEACOEFF", ueacoeff)


#udirates = eval_udi(udicoeff, kT)
#print("UDI:", udirates)


#uearates = eval_uea(ueacoeff, kT)
#print("UEA:", uearates)

#altuearates = eval_uea(ueacoeff, kT, altmethod=True)
#print("altUEA:", altuearates)

'''
tmpx = numpy.array([1000000,
1988011.86625617,
2879310.74238541,
5941731.78504284,
10267446.547635,
21557422.7223439,
45842133.5603972,
105420458.608543,
231012970.008316,
442838182.117041,
988238251.65563])

tmpy = numpy.array([2.05891626998812E-12,
5.67713295443297E-11,
1.47887064686471E-10,
3.244537012862E-10,
3.98138388383456E-10,
4.66104062060517E-10,
5.12091848890409E-10,
5.0491481852928E-10,
4.58491741169489E-10,
4.11468691505062E-10,
3.57305196573902E-10])


tmpx2 = numpy.array([1165270.83668186,
1467443.08078452,
1847065.14667663,
2324894.0287644,
2926069.14767006,
3683032.35075922,
4635819.11847169,
7344603.90872856,
9244625.6168961,
14645078.3547245,
23202436.7169887,
36756653.8465109,
58234166.7381248,
116118438.162252,
231560261.648564,
461729212.061864,
920767957.681571])

tmpy2 = numpy.array([8.68842218904234E-11,
2.88790116474045E-10,
7.18896052857088E-10,
1.415448007261E-09,
2.31644539470288E-09,
3.27495378854931E-09,
4.14837899614404E-09,
5.3545708850125E-09,
5.70301095313921E-09,
6.10995202362257E-09,
6.32641738340257E-09,
6.41181672408249E-09,
6.34579375062195E-09,
5.95387800755902E-09,
5.316915999644E-09,
4.60509452778059E-09,
3.97639196729588E-09])/10

#di, ea, rr, dr = pyatomdb.atomdb.get_ionrec_rate(priyanka_data['Temperature'], Te_unit='K', separate=True, Z=26, z1=16, extrap=True)

#print(kT*11604.5*1000, udirates)



fig = pylab.figure()
fig.show()
ax = fig.add_subplot(111)
#ax.loglog(kT*11604.5*1000, uearates, '+', label='EA')
#ax.loglog(kT*11604.5*1000, altuearates, '+', label='ALTEA')
#ax.loglog(kT*11604.5*1000, udirates, '^', label='DI')
ax.loglog(kT*11604.5*1000, udirates+uearates, '*', label='EA+DI')
ax.loglog(kT*11604.5*1000, udirates+altuearates, 'o', label='ALTEA+DI')
#ax.loglog(kT*11604.5*1000, priyanka_data['Total'], 'x', label='PC')
#ax.loglog(kT*11604.5*1000, priyanka_data['C_DI'], 'x', label='PC_DI')
#ax.loglog(kT*11604.5*1000, priyanka_data['C_EA'], 'x', label='PC_EA')
#ax.loglog(kT*11604.5*1000, di+ea, ':s', label='Bryans')
#ax.loglog(tmpx, tmpy, '*', label='paper')
#ax.loglog(tmpx2, tmpy2, '*', label='paper')
#ax.legend(loc=0)
pylab.draw()
zzz=input()
'''

