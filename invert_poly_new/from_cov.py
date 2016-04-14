from cov_T import *
from zeroth_fit_2var import *
from zeroth_fit import *
deca=3
#cov_T(File, molec,n_step,ene_min,ene_maxt_start,t_stop,order)

#cov_T('15N2_100K_0.16','n2',30,670,1650,25-deca,55,1)
#cov_T('15N2_100K_0.37','n2',30,670,1650,25-deca,55,1)
#cov_T('15N2_100K_0.72','n2',30,670,1650,29-deca,55,1)
#cov_T('15N2_100K_1.4','n2',30,670,1650,29.8-deca,55,1)#
#cov_T('15N2_14K','n2',55,1100,1900,35-deca,65,1)

#cov_T('13CO_100K_0.19','co',30,670,1650,25-deca,58,1)
#cov_T('13CO_100K_0.33','co',30,670,1650,25-deca,58,1)
#cov_T('13CO_100K_0.75','co',30,670,1650,31.8-deca,58,1)
#cov_T('13CO_100K_1.26','co',30,670,1650,32.5-deca,58,1)
#cov_T('13CO_14K','co',35,1300,2000,32-deca,63,1)



#Good prams N2

#print '------------------  15N2 2var fits ---------------------------'
zeroth_fit_2var(11,767,7e11,23.1,27.6)

#print '------------------ Zeroth order 15N2 1var---------------------------'
#zeroth_fit('15N2_dec',deca,'n2',4.5,767,7e11,24.55,27.9)


#print '------------------  13CO 2var fits ---------------------------'
#zeroth_fit_2var('13CO_bare_4.95',deca,'co',4.95,855,7e11,26.41,30.8)

#print '------------------ Zeroth order 13CO 1var---------------------------'
#zeroth_fit('13CO_bare_4.95',deca,'co',4.95,855,7e11,26.41,30.8)


#zeroth_fit('CO2_TPD',deca,'co2',5.6,2767.8,70,77)
#zeroth_fit_2var('H2O_TPD',deca,'h2o',4.5,5773,1e15,128,145)
#zeroth_fit('13CO_test',6.0,837,22.5,27.8)