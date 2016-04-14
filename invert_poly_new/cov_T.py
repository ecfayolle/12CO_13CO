import numpy as np
import math
import matplotlib.pyplot as plt
import scipy
from matplotlib import rc, font_manager
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.integrate import odeint 
from scipy import integrate
from scipy.integrate import quad
from scipy import linalg
from scipy import optimize
from scipy import ndimage as nd
import csv

deca=3
#cov_T('15N2_100K_0.16','n2',30,26,55,1)
#cov_T('15N2_100K_0.37','n2',30,26,55,1)
#cov_T('15N2_100K_0.72','n2',30,29,55,1)
#cov_T('15N2_100K_1,4','n2',30,29,55,1)
#cov_T('15N2_14K','n2',30,26,55,1)


def cov_T(File, molec,n_step,ene_min,ene_max,t_start,t_stop,order):

	sizeOfFont = 16
	fontProperties = {'family':'sans-serif','sans-serif':['cm'],
	    'weight' : 'normal', 'size' : sizeOfFont}
	#rc('text', usetex=True)
	rc('font',**fontProperties)
	ticks_font = font_manager.FontProperties(family='cm', style='normal',
	    size=sizeOfFont, weight='normal', stretch='normal')

	if molec =='n2':
		m_mol=30.0e-3
		E_zero=763

	if molec =='co':
		m_mol=29.0e-3
		E_zero=859
		
	print 'facto'

		
	filenm=File+'_T_molec_dose.txt'

	T_tot, N_tot = np.loadtxt(filenm,unpack=True)
	T_tot=T_tot-deca


	h_rate=1/60.0 #Kmin-1
	v_pump=0.00000

	T=T_tot[np.where((T_tot>t_start) & (T_tot<t_stop))]
	N=N_tot[np.where((T_tot>t_start) & (T_tot<t_stop))]

	N_max=scipy.integrate.trapz(N,x=T)/1.0e15
	print 'coverage',N_max
	#y= scipy.integrate.cumtrapz(N, x=T)
	#covt=N_max-y

	#E_fit = np.linspace(ene_min, ene_max, n_ene)
	E_fit = np.arange(ene_min, ene_max, n_step)

	A = np.zeros((len(T),len(E_fit)),np.float32)
	b = N/1.0e15

	success = 0
	while (not success):

		for i in np.arange(len(E_fit)):
	  		nu_des=np.sqrt(2.0*1e19*8.31*E_fit[i]/(np.pi*np.pi*m_mol))
	  		print nu_des
	  		def funct(x,t):
	  			f_ml = -(nu_des/h_rate)*(x[0]**order)*np.exp(-E_fit[i]/t)
	  			f_g= (nu_des/h_rate)*(x[0]**order)*np.exp(-E_fit[i]/t)-v_pump*x[1]
	  			return [f_ml,f_g]
		  	F_pure = odeint(funct,[N_max,0.00],T)
		  	N_ml = np.clip(F_pure[:,0],0.0,1)
		  	N_g = F_pure[:,1]
		  	k_d_ml = (nu_des/h_rate)*np.exp(-E_fit[i]/T)*(N_ml**order)
		  	print "Ns",N_ml
		  	#if order==0:
		  	#	k_d_ml[np.where(k_d_ml>np.max(N))]=0.0
		  	A[:,i] = k_d_ml
                #A[:,i] = np.clip(k_d_ml,-1e15,1e15)
        
		x,resid = optimize.nnls(A,b)
		nx = nd.gaussian_filter(x,1.5)
		bfit = np.dot(A,x)
		resid = np.sum((bfit-b)*(bfit-b),axis=0)/np.sum(b*b,axis=0)
		print ("resid",resid)
		if resid < 0.05:
			success = 1


	k_d_z = (np.sqrt(2.0*1e19*8.31*E_zero/(np.pi*np.pi*m_mol))/h_rate)*np.exp(-E_zero/T_tot)


	#T=T_tot[np.where((T_tot>t_start) & (T_tot<t_stop))]


	#print(E_fit[np.where(x==np.max(x))])
	mefit = np.sum(nx*E_fit,axis=0)/np.sum(nx,axis=0)
	print 'mean E_fit',mefit
	print 'width E_fit',(np.sum(nx*E_fit*E_fit,axis=0)/np.sum(nx,axis=0)-mefit**2)**0.5


	##### plot the TPD fit

	plot_tpd=File+'_tpdfit.pdf'
	plt.plot(T_tot,N_tot*1e-13,color='k')
	plt.plot(T,bfit*1e2,color='r',linestyle='--',linewidth=3.0)
	#plt.plot(T_tot,cfit*1e2,color='r',linestyle='--',linewidth=3.0)
	#plt.plot(T_tot,1e-13*np.clip(k_d_z,-1,1e15),color='r',linestyle='--',linewidth=3.0)
	plt.xlim(15,55)
	plt.ylim(-np.max(N_tot*1e-13)*0.01,np.max(N_tot*1e-13)*1.1)
	plt.xlabel(r'Temperature / K')
	plt.ylabel(r' Desorption rate / 10$^{13}$ molecules.K$^{-1}$')
	print 'TPD integ',N_max
	print 'Fit integ',scipy.integrate.trapz(bfit, x=T)
	plt.savefig(plot_tpd)
	plt.close('all')


	##### plot the energy distribution
	plot_ene=File +'_enedistri.pdf'
	file_ene=File +'_enedistri.csv'
	file_fit=File+'_fit.csv'
	plt.scatter(E_fit,x,color='k')
	plt.plot(E_fit,nx,color='k',linestyle='--',linewidth=2.0)
	plt.xlabel(r'Desorption energy / K')
	plt.ylabel(r'Fractional coverage')
	plt.savefig(plot_ene)
	plt.close('all')

	with open(file_ene, 'w') as f:
		writer = csv.writer(f)
		writer.writerows(zip(E_fit,x,nx))

	with open(file_fit, 'w') as f:
		writer = csv.writer(f)
		writer.writerows(zip(T,bfit*1e15))

	np.savetxt(File+'_sol.txt', (mefit,(np.sum(nx*E_fit*E_fit,axis=0)/np.sum(nx,axis=0)-mefit**2)**0.5))

