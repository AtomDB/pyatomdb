import pyatomdb
import class_apec
import atomdb,os
import numpy
Abund= atomdb.get_abundance(abundfile=False, abundset='AG89', element=[-1], datacache=False, settings=False, show=False)
settings=class_apec.parse_par_file(os.path.expandvars('$ATOMDB/apec.par'))
d=class_apec.element(26, 40000000, 1, 0, 1.0, 'None', Abund)
ionfrac =  d.calc_elem_ionbal_delta(26, 40000000, 1.0, 'None')
a=class_apec.variableapec(26,25,10000000,1,[2,5],[1,1],[[0.8,0.8],[0.8,0.8]],['LA','EC'],[2,5],[1,1],0,Abund,settings,ionfrac)
a.run_apec_ion_delta_single_exc(settings, 26,25, 10000000,1,[2,5],[1,1], [[0.8,0.8],[0.8,0.8]],['LA','EC'],[2,5],[1,1],ionfrac, Abund)
