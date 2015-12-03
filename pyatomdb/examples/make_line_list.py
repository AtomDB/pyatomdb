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

res = pyatomdb.spectrum.list_lines(wl,Te=Te, teunit='K', minepsilon=1e-18)

# reprocess lines for printing
print "Unsorted line list:"
pyatomdb.spectrum.print_lines(res)

# re-sort lines, for a giggle
# for more information, look up numpy.sort: res is a numpy array.
# http://docs.scipy.org/doc/numpy/reference/generated/numpy.sort.html

res.sort(order=['Epsilon'])
print "sorted by Emissivity:"
pyatomdb.spectrum.print_lines(res)

# re-sort by element, ion then emissivity
res.sort(order=['Element','Ion','Epsilon'])
print "sorted by Element, Ion, Emissivity:"
pyatomdb.spectrum.print_lines(res)


# now do an NEI version. This is slow at the moment, but functional.
Te_init = 1e4
tau = 1e11
res_nei = pyatomdb.spectrum.list_nei_lines(wl,Te=Te, teunit='K', \
                                           minepsilon=1e-18,\
                                           Te_init=Te_init,\
                                           tau = tau)
print "NEI linelist (this takes a while):"
pyatomdb.spectrum.print_lines(res_nei)
                        






