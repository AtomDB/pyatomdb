
import numpy
import xspec
import pyatomdb


# These are globals to hold the model data
kappamodelobject = pyatomdb.spectrum.KappaSession(
                   elements=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                             13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                             23, 24, 25, 26, 28])
kappavmodelobject = False
kappavvmodelobject = False


kappamodelobject.abund_xspectoatomdb = {}
kappamodelobject.abund_xspectoatomdb['angr'] = 'AG89'
kappamodelobject.abund_xspectoatomdb['aspl'] = 'Asplund09'
kappamodelobject.abund_xspectoatomdb['feld'] = 'Feldman'
kappamodelobject.abund_xspectoatomdb['aneb'] = 'AE82'
kappamodelobject.abund_xspectoatomdb['grsa'] = 'GS98'
kappamodelobject.abund_xspectoatomdb['wilm'] = 'Wilm00'
kappamodelobject.abund_xspectoatomdb['lodd'] = 'Lodd03'

pykappaInfo = ("kT            \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
             "kappa         \"\"      3.0 2.0 2.0 100. 1000. 0.01",
             "abund         \"\"      1.0 0.0 0.0 10.0 10.0 0.01")

pyvkappaInfo = ("kT            \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
              "kappa         \"\"      3.0 2.0 2.0 100. 1000. 0.01",
              "H             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "He            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "C             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "N             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "O             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "F             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ne            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Mg            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Al            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Si            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "S             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ar            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ca            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Fe            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ni            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01")


pyvvkappaInfo = ("kT            \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
               "kappa         \"\"      3.0 2.0 2.0 100. 1000. 0.01",
               "H             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "He            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Li            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Be            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "B             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "C             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "N             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "O             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "F             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Ne            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Na            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Mg            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Al            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Si            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "P             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "S             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Cl            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Ar            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "K             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Ca            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Sc            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Ti            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "V             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Cr            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Mn            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Fe            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
               "Ni            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01")


def pykappa(engs, params, flux):
  """
  Kaapa model for data

  PARAMETERS
  ----------
  engs : list[float]
    The energy bin edges (from xspec)
  params : list[float]
    The parameter list. See kappaInfo for definition
  flux : list[float]
    The array to fill with return values

  RETURNS
  -------
  None
    Fills out the flux array with photon cm3 s-1 bin-1 x 1e14

  USAGE
  -----
    # load the model into XSPEC
    xspec.AllModels.addPyMod(kappa, kappaInfo, 'add')
    # make a model
    m = xspec.Model('kappa')
  """

  # This is the call that will return everything. So set everything!
  ebins = numpy.array(engs)
  kappamodelobject.set_response(ebins, raw=True)
  # kappa model has the 14 main elements
  elarray = kappamodelobject.elements
  abund = numpy.zeros(len(elarray))

  
  if len(params)==4:
    # kappa case
    elarray=[2,6,7,8,10,12,13,14,16,18,20,26,28]
    abund = float(params[2])
  elif len(params)==18:
    # vkappa case
    elarray=[1,2,6,7,8,9,10,12,13,14,16,18,20,26,28]
    abund = numpy.array(params[2:17])
  elif len(params)==30:
    # vvkappa case
    elarray=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,\
             20,21,22,23,24,25,26,28]
    abund = numpy.array(params[2:29])
    

  T = float(params[0])
  k = float(params[1])

  # correct abundance set if required
  if kappamodelobject.abund_xspectoatomdb[xspec.Xset.abund.split(':')[0]] != \
     kappamodelobject.abundset:
    kappamodelobject.set_abundset(kappamodelobject.abund_xspectoatomdb[
     xspec.Xset.abund.split(':')[0]])

  # set the abundance

  kappamodelobject.set_abund(elarray, abund)

  spec = kappamodelobject.return_spectrum(T, k)

  flux[:] = spec*1e14


def pyvkappa(engs, params, flux):
  pykappa(engs,params,flux)

def pyvvkappa(engs, params, flux):
  pykappa(engs,params,flux)


#def pyvkappa(engs, params, flux):
#
#  """
#  Kaapa model for data
#
#  PARAMETERS
#  ----------
#  engs : list[float]
#    The energy bin edges (from xspec)
#  params : list[float]
#    The parameter list. See vkappaInfo for definition
#  flux : list[float]
#    The array to fill with return values
#
#  RETURNS
#  -------
#  None
#    Fills out the flux array with photon cm3 s-1 bin-1 x 1e14
#
#  USAGE
#  -----
#    # load the model into XSPEC
#    xspec.AllModels.addPyMod(vkappa, vkappaInfo, 'add')
#    # make a model
#    m = xspec.Model('vkappa')
#  """
#
#  # This is the call that will return everything. So set everything!
#  ebins = numpy.array(engs)
#  kappamodelobject.set_response(ebins, raw=True)
#
#  # kappa model has the 14 main elements
##  kappamodelobject.elements = [1,2,6,7,8,10,12,13,14,16,18,20,26,28]
#
#  abundset = xspec.Xset.abund
#
#  T = params[0]
#  k = params[1]
#
#  # correct abundance set if required
#  if kappamodelobject.abund_xspectoatomdb[xspec.Xset.abund] != \
#     kappamodelobject.abundset:
#    kappamodelobject.set_abundset(kappamodelobject.abund_xspectoatomdb[
#     xspec.Xset.abund])
#
#  # set the abundance
#  elarray = kappamodelobject.elements
#  abund = numpy.zeros(len(elarray))
#  abund[numpy.array([1, 2, 6, 7, 8, 10, 12, 13, 14, 16, 18,
#                     20, 26, 27])-1] = params[2:16]
#  kappamodelobject.set_abund(elarray, abund)
#
#  spec = kappamodelobject.return_spectrum(T, k)
#
#  flux[:] = spec*1e14
#
#
#def vvkappa(engs, params, flux):
#
#  """
#  Kaapa model for data
#
#  PARAMETERS
#  ----------
#  engs : list[float]
#    The energy bin edges (from xspec)
#  params : list[float]
#    The parameter list. See vvkappaInfo for definition
#  flux : list[float]
#    The array to fill with return values
#
#  RETURNS
#  -------
#  None
#    Fills out the flux array with photon cm3 s-1 bin-1 x 1e14
#
#  USAGE
#  -----
#    # load the model into XSPEC
#    xspec.AllModels.addPyMod(vvkappa, vvkappaInfo, 'add')
#    # make a model
#    m = xspec.Model('vvkappa')
#  """
#
#  # This is the call that will return everything. So set everything!
#  ebins = numpy.array(engs)
#  kappamodelobject.set_response(ebins, raw=True)
#
#  # kappa model has the 14 main elements
#  Zlist = list(range(1, 27))
#  Zlist.append(28)
#  kappamodelobject.elements = Zlist
#
#  abundset = xspec.Xset.abund
#
#  T = params[0]
#  k = params[1]
#  abund = {}
#  abund[1] = params[2]
#  abund[2] = params[3]
#  abund[3] = params[4]
#  abund[4] = params[5]
#  abund[5] = params[6]
#  abund[6] = params[7]
#  abund[7] = params[8]
#  abund[8] = params[9]
#  abund[9] = params[10]
#  abund[10] = params[11]
#  abund[11] = params[12]
#  abund[12] = params[13]
#  abund[13] = params[14]
#  abund[14] = params[15]
#  abund[15] = params[16]
#  abund[16] = params[17]
#  abund[17] = params[18]
#  abund[18] = params[19]
#  abund[19] = params[20]
#  abund[20] = params[21]
#  abund[21] = params[22]
#  abund[22] = params[23]
#  abund[23] = params[24]
#  abund[24] = params[25]
#  abund[25] = params[26]
#  abund[26] = params[27]
#  abund[28] = params[28]
#
#  # correct abundance set if required
#  if kappamodelobject.abund_xspectoatomdb[xspec.Xset.abund] != \
#     kappamodelobject.abundset:
#    kappamodelobject.set_abundset(kappamodelobject.abund_xspectoatomdb[
#                                                      xspec.Xset.abund])
#
#  # set the abundance
#  elarray = kappamodelobject.elements
#  kappamodelobject.set_abund(elarray, params[2:29])
#
#  spec = kappamodelobject.return_spectrum(T, k)
#
#  flux[:] = spec*1e14
#

xspec.AllModels.addPyMod(pykappa, pykappaInfo, 'add')
xspec.AllModels.addPyMod(pyvkappa, pyvkappaInfo, 'add')
xspec.AllModels.addPyMod(pyvvkappa, pyvvkappaInfo, 'add')
print("")
print("Models kappa, vkappa, and vvkappa ready for use")
