import numpy, xspec
import pyatomdb


#initialize CX object
cie_rs = pyatomdb.spectrum.CIESession_RS()
cie_rs.set_broadening(True)

# These are the definitions XSPEC uses for the inputs



pyrsapecInfo = ("kT   \"keV\"   4.0 0.00862 0.00862 86. 86. 0.01",
            "Abundanc   \"\" 0.6 0.01 0.2 100. 1000. 0.01",
            "*Redshift    \"\"     0.0",\
            "Velocity \"km/s\" 150.0 -10.0 0. 0. 1e6 1e6",
            "nL   \"cm^-2\" 1e21 1e20 1e20 1e23 1e23 1e18")



pyrsvapecInfo = ("kT   \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
            "H   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "He   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "C   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "N   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "O   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Ne   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Mg   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Al   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Si   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "S   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Ar   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Ca   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Fe   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Ni   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "*Redshift    \"\"     0.0",\
            "Velocity \"km/s\" 0.0 -10.0 0. 0. 1e6 1e6",
            "nL   \"cm^-2\" 1e21 1e20 1e20 1e23 1e23 1e18")


pyrsvvapecInfo = ("kT   \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
            "H   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "He   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Li   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Be   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "B    \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "C   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "N   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "O   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "F   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Ne   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Na   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Mg   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Al   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Si   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "P   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "S   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Cl  \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Ar   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "K   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Ca   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Sc   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Ti   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "V    \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Cr   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Mn   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Fe   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Co   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Ni   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Cu   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "Zn   \"\" 1.0 0.01 0.2 100. 1000. 0.01",
            "*Redshift    \"\"     0.0",\
            "Velocity \"km/s\" 0.0 -10.0 0. 0. 1e6 1e6",\
            "nL   \"cm^-2\" 1e21 1e20 1e20 1e23 1e23 1e18")


#"Density   \"cm^-3\" 0.03 1e-3 1e-3 0.05 0.05 1e-3")

def pyapecrs(engs, params, flux):

  """
  ACX2 model for data

  PARAMETERS
  ----------
  engs : list[float]
    The energy bin edges (from xspec)
  params : list[float]
    The parameter list. See acx2Info for definition
  flux : list[float]
    The array to fill with return values

  RETURNS
  -------
  None
    Fills out the flux array with photon cm3 s-1 bin-1, x1e10

  USAGE
  -----
    # load the model into XSPEC
    xspec.AllModels.addPyMod(pyapec, pyapecInfo, 'add')
    # make a model
    m = xspec.Model('pyapec')
  """

  # This is the call that will return everything. So set everything!
  ebins = numpy.array(engs)
  #print(ebins)
  #ebins = numpy.linspace(6.0,7.0,10000)

  elements = []
  for i in range(1,31): elements.append(i)

  if len(params)==6:
    # apec model
    abund=numpy.array(params[1])
    elements=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,\
              16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
    eloffset = 2

  elif len(params)==19:
    # apec model
    elements=[1,2,6,7,8,10,12,13,14,16,18,20,26,28]
    abund=params[1:15]
    eloffset = 15

  elif len(params)==35:
    # apec model
    elements=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,\
              16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
    abund=params[1:31]
    eloffset = 31

  
  redshift = float(params[eloffset+0])
  
  broadening = float(params[eloffset+1])

  nL_integral = float(params[eloffset+2])
  
  
  # set energy bins, accounting for redshift
  cie_rs.set_response(ebins*(1.0+redshift), raw=True)


  # set abundance vector
#  abund = numpy.array(params[1])
  cie_rs.set_abund(elements, abund)
  cie_rs.set_eebrems(False)
  
  # set the response
  cie_rs.set_broadening(True, velocity_broadening=broadening)


  
  # get the spectrum
  spec = cie_rs.return_spectrum(params[0],  N_e= nL_integral, Ab=abund,  log_interp=False)

  # return the flux.
  flux[:] = spec*1e14


def pyvapecrs(engs, params, flux):
  pyapecrs(engs, params, flux)
  

def pyvvapecrs(engs, params, flux):
  pyapecrs(engs, params, flux)




# this is how to import the models into pyxspec.
xspec.AllModels.addPyMod(pyapecrs, pyrsapecInfo, 'add')
xspec.AllModels.addPyMod(pyvapecrs, pyrsvapecInfo, 'add')
xspec.AllModels.addPyMod(pyvvapecrs, pyrsvvapecInfo, 'add')