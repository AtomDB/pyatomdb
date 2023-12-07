import pyatomdb

"""
This code will produce a list of lines in a given wavelength range at a
given temperature. It also shows the use of an NEI version, where you 
have to additionally specify the initial ionization temperature (or the
ionization fraction directly) and the elapsed Ne*t.

The results of the list_lines codes are numpy arrays which can be sorted any
way you wish. You can, of course, extract the lines easily at this point. There
is also a print_lines routine for a fixed format output.

Parameters
----------
none

Returns
-------
none

"""

# Adam Foster 2015-12-02
# version 0.1

#specify wavelength range, in Angstroms
wl = [8.0,9.0]

# electron temperature in K
Te = 1e7

# get equilibrium line list
CIE = pyatomdb.spectrum.CIESession()

res = CIE.return_linelist(Te, wl, specunit = 'A', teunit='K')

#trim off the weaker lines
res = res[res['Epsilon']> 1e-18]

# reprocess lines for printing
print("Unsorted line list:")
print(CIE.format_linelist(res))

# re-sort lines
# for more information, look up numpy.sort: res is a numpy array.
# http://docs.scipy.org/doc/numpy/reference/generated/numpy.sort.html

res.sort(order=['Epsilon'])
res=res[::-1]
print("Top 10 sorted by Emissivity:")
print(CIE.format_linelist(res[:10]))

# re-sort by element, ion then emissivity
res.sort(order=['Element','Ion','Epsilon'])
print("Top 10 sorted by Element, Ion, Emissivity:")
print(CIE.format_linelist(res[:10]))


# now do an NEI version. This is slow at the moment, but functional.
Te_init = 1e4
tau = 1e11
NEI = pyatomdb.spectrum.NEISession()
res_nei = NEI.return_linelist(Te, tau, wl, teunit='K', \
                              specunit='A',\
                              init_pop=Te_init)
res_nei = res_nei[res_nei['Epsilon']> 1e-19]
print("NEI linelist:")
print(NEI.format_linelist(res_nei[:10]))
                        






