import unittest, warnings
import pyatomdb

ATOMDBVERSION='309'



class TestSpectrum(unittest.TestCase):


    def test_convert_temp(self):
      """ Test temperature conversion in test_convert_temp"""
      T={}
      T['keV'] = 2.0

      Unitsin=['K','keV','eV']
      Unitsout=['K','keV','eV']
      T['K'] = T['keV']/pyatomdb.const.KBOLTZ
      T['eV'] = T['keV']*1000
      for Uin in Unitsin:
        for Uout in Unitsout:
          c = T[Uout]
          i = T[Uin]
          pyatomdb.numpy.testing.assert_allclose(c, pyatomdb.spectrum.convert_temp(i,Uin,Uout))


      self.assertRaises(pyatomdb.util.UnitsError, pyatomdb.spectrum.convert_temp, 1.0,'krf','kev')

      self.assertRaises(pyatomdb.util.UnitsError, pyatomdb.spectrum.convert_temp, 1.0,'kev','krh')



class TestSpectrumCIESession(unittest.TestCase):

    def setUp(self):
      warnings.simplefilter("ignore", ResourceWarning)

      self.s = pyatomdb.spectrum.CIESession()

    def test_init(self):
      """ Testing initialization of the CIESession"""
      self.assertEqual(self.s.broaden_limit, 1e-18)


    def test_set_raw_response(self):
      """ Testing dummy (raw) responses"""

      self.assertRaises(pyatomdb.util.ReadyError, self.s.return_spectrum,1.0)

      ebins=pyatomdb.numpy.linspace(0.1,10,1000)
      self.s.set_response(ebins*1.0, raw=True)
      pyatomdb.numpy.testing.assert_allclose(self.s.specbins, ebins)

      # now test spectrum calculation

      # load the correct result
      testdat = pyatomdb.numpy.load("tests/testdata/test_set_response_spec1_%s.npy"%(ATOMDBVERSION))

      t = self.s.return_spectrum(1.0)
      pyatomdb.numpy.testing.assert_allclose(t, testdat)


    def test_set_response(self):
      """ Testing real responses and spectral broadening"""
      # test with real response
      self.s.set_response('tests/testdata/aciss_heg1_cy19.grmf',\
                          arf='tests/testdata/aciss_heg1_cy19.garf')

      # load the correct result
      testdat = pyatomdb.numpy.load("tests/testdata/test_set_response_spec2_%s.npz"%(ATOMDBVERSION))

      # calc spectrum
      t = self.s.return_spectrum(1.0)
      pyatomdb.numpy.testing.assert_allclose(t, testdat['spec'])


      # turn on thermal broadening
      self.s.set_broadening(True)
      t = self.s.return_spectrum(1.0)
      pyatomdb.numpy.testing.assert_allclose(t, testdat['tbroadspec'])

    def test_nearest_spectrum(self):
      """ Testing real responses with the nearest spectrum"""
      # test with real response
      self.s.set_response('tests/testdata/aciss_heg1_cy19.grmf',\
                          arf='tests/testdata/aciss_heg1_cy19.garf')

      # load the correct result
      testdat = pyatomdb.numpy.load("tests/testdata/test_set_response_spec2_%s.npz"%(ATOMDBVERSION))

      # calc spectrum
      t = self.s.return_spectrum(1.0)
      pyatomdb.numpy.testing.assert_allclose(t, testdat['spec'])


      # turn on thermal broadening
      self.s.set_broadening(True)
      t = self.s.return_spectrum(1.0)
      pyatomdb.numpy.testing.assert_allclose(t, testdat['tbroadspec'])


if __name__ == '__main__':
    print("WHEE")
    unittest.main()
    print("WHOO")
