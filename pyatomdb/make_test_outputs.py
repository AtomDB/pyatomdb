import pyatomdb

pyatomdbversion=pyatomdb.__version__

PYADBVERS = pyatomdbversion.replace('.','')

ADBVERS = '309'

s = pyatomdb.spectrum.CIESession()

ebins=pyatomdb.numpy.linspace(0.1,10,1000)
s.set_response(ebins*1.0, raw=True)
y = s.return_spectrum(1.0)

pyatomdb.numpy.save("tests/testdata/test_set_response_spec1_%s.npy"%(ADBVERS),y)

     
s = pyatomdb.spectrum.CIESession()
t = {}
s.set_response('tests/testdata/aciss_heg1_cy19.grmf',\
               arf='tests/testdata/aciss_heg1_cy19.garf')

t['spec'] = s.return_spectrum(1.0)

s.set_broadening(True)
t['tbroadspec'] = s.return_spectrum(1.0)

s.set_response('tests/testdata/aciss_heg1_cy19.grmf',arf=False)
t['tbroadspec_noarf'] = s.return_spectrum(1.0)

pyatomdb.numpy.savez("tests/testdata/test_set_response_spec2_%s.npz"%(ADBVERS),\
                     spec=t['spec'],\
                     tbroadspec=t['tbroadspec'],\
                     tbroadspec_noarf=t['tbroadspec_noarf'])

