import numpy as np
import math
import matplotlib 
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
from scipy import stats
import csv
from scipy.odr import ODR, Model, Data, RealData

def zeroth_fit_2var(Nml,Eth,nu_des_ini,t_start,t_stop):
	
	deca=0
	sizeOfFont = 16
	fontProperties = {'family':'sans-serif',
	    'weight' : 'normal', 'size' : sizeOfFont}
	#rc('text', usetex=True)
	rc('font',**fontProperties)
	ticks_font = font_manager.FontProperties(family='cm', style='normal',
	    size=sizeOfFont, weight='normal', stretch='normal')
	filenm='12CO_13CO_1_1_fit.txt'


	#if molec =='n2':
	#	m_mol=30.0e-3
	#	xtex=0.042

	#if molec =='co':
	m_mol=28.0e-3
	xtex=0.0375

	#if molec =='co2':
	#	m_mol=44.0e-3
	#if molec =='h2o':
	#	m_mol=18.0e-3

	errort=0.1

	T_tot, N_tot,N_13co = np.loadtxt(filenm,unpack=True,skiprows=1)
	N_tot=N_13co
	T_ini=T_tot
	T_tot=T_tot
	T_tot_err=T_tot*0.0+errort
	N_tot=np.clip(N_tot,1e-12,1e30)
	N_tot_err=1e12+N_tot*0

	h_rate=1/60.0 #Kmin-1

	T=T_tot[np.where((T_tot>t_start) & (T_tot<t_stop))]
	N=N_tot[np.where((T_tot>t_start) & (T_tot<t_stop))]
	N_err=1e10+N*0
	T_err=T*0.0+errort
	
	x=1.0/T

	x_err=T_err/(T**2)
	
	#print x_err
	y=np.log(N/(1e15))
	y_err=N_err/N
	#y=np.log(N)
	p=[0.,0.]
	#nu_des_ini=7e11#np.sqrt(2.0*1e19*8.31*Eth/(np.pi*np.pi*m_mol))

	#slope, intercept, r_value, p_value, std_err=scipy.stats.linregress(x, y)
	def func(p, t):
		return p[1] +t*p[0]
	#def func (p,t):
	#	return np.log(1e12/h_rate)+t*p[0]
	#def func(p, t):
	#	return np.log(np.sqrt(2.0*1e19*8.31*-1*p[0]/(np.pi*np.pi*m_mol))/h_rate) +t*p[0]
	#p[1]=np.log(1e12/h_rate)
	
	#p,pcov=curve_fit(func, x, y)
	#err = np.sqrt(np.diag(pcov))
	data = RealData(x, y, sx=x_err,sy=y_err)#, sy=np.log(0.001))
	
	model = Model(func)

	#odr = ODR(data, model, beta0=[-800,1e12])
	odr = ODR(data, model, beta0=[-750,7e11])
	odr.set_job(fit_type=2)
	output = odr.run()

	p=output.beta
	err=output.sd_beta

	E_des=-p[0]
	nu_des=h_rate*np.exp(p[1])
	nu_theo=h_rate*np.exp(np.log(np.sqrt(2.0*1e19*8.31*-1*p[0]/(np.pi*np.pi*m_mol))/h_rate))
	#p[1]=np.log(np.sqrt(2.0*1e19*8.31*-1*p[0]/(np.pi*np.pi*m_mol))/h_rate)
	#err[1]=0.0
	
	#nu_theo=np.sqrt(2.0*1e19*8.31*Eth/(np.pi*np.pi*m_mol))
	print 'E_des fit', p[0]*-1, '+-',err[0]
	print 'E_des litt', Eth
	print 'nu fit e12', 1e-12*h_rate*np.exp(p[1]), '+-',err[1]*1e-12*h_rate*np.exp(p[1])

	print 'nu litt e12', 1e-12*nu_des_ini
	print 'nu theo e12', 1e-12*nu_theo


	#nu_des=2e11
	#print "nu e12n", nu_des*1e-12
	# def funct(x,t):
	# 	if x[0]>0.0:
	#   		f_ml = -(nu_des/h_rate)*np.exp(-E_des/t)
	#   		f_g= f_ml-x[1]
	#   	if x[0]<0.0:
	#   		f_ml=0
	#   		f_g= f_ml-x[1]
	#   	return [f_ml,f_g]
		
	# F_pure = odeint(funct,[Nml,0.00],T_tot)
	# N_ml = np.clip(F_pure[:,0],0.0,10)
	# N_g = -F_pure[:,1]

	# k_d_ml = (nu_des/h_rate)*np.exp(-E_des/T_tot)
	

	# k_d_mode=(nu_des_ini/h_rate)*np.exp(-Eth/T_tot)
		#print "Ns",N_ml
		  	#if order==0:
		  	#	k_d_ml[np.where(k_d_ml>np.max(N))]=0.0
		#A[:,i] = k_d_ml
                #A[:,i] = np.clip(k_d_ml,-1e15,1e15)
        
		#x,resid = optimize.nnls(A,b)

	#slope, intercept, r_value, p_value, std_err=scipy.stats.linregress(x, y)

	plot_tpd='12CO_reglinfit_2var.pdf'
	plt.errorbar(x,y,xerr=x_err,yerr=y_err,color='k')
	#plt.errorbar(a,b,yerr=c, linestyle="None")
	plt.plot(x,p[1]+p[0]*x,color='r',linestyle='--',linewidth=3.0)
	#plt.plot(x,p[1]+p[0]*x,color='r',linestyle='--',linewidth=3.0)

	plt.text(xtex,0.0,'-Slope: '+str(round(E_des))+' +- '+str(round(err[0])))
	plt.text(xtex,-0.5,'exp(Intercept)*h_rate: '+ str(round(1e-12*h_rate*np.exp(p[1]),2))+ '+/-'+str(round(err[1]*1e-12*h_rate*np.exp(p[1]),2))+ 'e12')
	#plt.plot(x,p[0]+(p[1]-err[1])*x,color='r',linestyle='--',linewidth=3.0)
	#plt.plot(x,p[0]+(p[1]+err[1])*x,color='r',linestyle='--',linewidth=3.0)
	#plt.plot(x,np.log(1e12/h_rate)-Eth*x,color='b',linestyle='--',linewidth=3.0)
	#plt.xlim(1/(t_stop+5),1/(t_start-5))
	plt.savefig(plot_tpd)
	plt.close('all')
	
	#ind=np.where(N_tot<np.max(N_tot))

	#plot_tpd='12CO_tpdfit_2var.pdf'
	#plt.errorbar(T_tot[ind],N_tot[ind]/(1e15),xerr=T_tot_err[ind],yerr=N_tot_err[ind]/1e15,color='k')
	#plt.plot([t_start-deca,t_start-deca],[0,np.max(N_tot)/1e15],linestyle='--',color='g')
	#plt.plot([t_stop-deca,t_stop-deca],[0,np.max(N_tot)/1e15],linestyle='--',color='g')
	#plt.plot(T,-np.exp(intercept)*np.exp(-1*slope/T),color='r',linestyle='--',linewidth=3.0)
	#plt.plot(T,(60*1.5e12/5.3e15)*np.exp(-950.0/T),color='b',linestyle='--',linewidth=3.0)
	#plt.plot(T_tot,np.clip(k_d_ml,0,2*np.max(N_tot/1e15)),color='r')
	#plt.plot(T_tot,np.clip(k_d_mode,0,2*np.max(N_tot/1e15)),color='r')
	
	#plt.xlim(15,50)
	#plt.ylim(-0.005,np.max(N_tot)/1e15+0.1)
	#plt.savefig(plot_tpd)
	plt.close('all')