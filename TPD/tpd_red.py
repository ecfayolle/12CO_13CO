import numpy as np
import csv
import xlrd
import time
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from scipy import interpolate as interp
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
import scipy

def savitzky_golay(y, window_size, order, deriv=0, rate=1):

    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

file_temp=['TPD_12CO_13CO_1_1_1K_03242016.xls','TPD_12CO_13CO_1_1_1K_03242016.xls','TPD_12CO_13CO_1_1_1K_03242016.xls',
			'TPD_12CO_13COonH2O_1_1_1K_03242016.xls','TPD_12CO_13COonH2O_1_1_1K_03242016.xls','TPD_12CO_13COonH2O_1_1_1K_03242016.xls',
			'TPD_12CO_13CO_onH2O_03252016.xls','TPD_12CO_13CO_onH2O_03252016.xls','TPD_12CO_13CO_onH2O_03252016.xls','TPD_12CO_13CO_onH2O_03252016.xls',
			'TPD_12CO_13CO_10_1_onH2O_03292016.xls','TPD_12CO_13CO_10_1_onH2O_03292016.xls','TPD_12CO_13CO_10_1_onH2O_03292016.xls','TPD_12CO_13CO_10_1_onH2O_03292016.xls',
			'TPD2_12CO_13CO_10_1_onH2O_03292016.xls','TPD2_12CO_13CO_10_1_onH2O_03292016.xls','TPD2_12CO_13CO_10_1_onH2O_03292016.xls','TPD2_12CO_13CO_10_1_onH2O_03292016.xls',
			'TPD_12CO_13C18O_H2O_04072016.xls','TPD_12CO_13C18O_H2O_04072016.xls','TPD2_12CO_13C18O_H2O_04072016.xls','TPD2_12CO_13C18O_H2O_04072016.xls']

file_mass=['TPD_12CO_13CO_1_1_1K_03242016.csv','TPD_12CO_13CO_1_1_1K_03242016.csv','TPD_12CO_13CO_1_1_1K_03242016.csv',
			'TPD_12CO_13COonH2O_1_1_1K_03242016.csv','TPD_12CO_13COonH2O_1_1_1K_03242016.csv','TPD_12CO_13COonH2O_1_1_1K_03242016.csv',
			'TPD_12CO_13CO_onH2O_03252016.csv','TPD_12CO_13CO_onH2O_03252016.csv','TPD_12CO_13CO_onH2O_03252016.csv','TPD_12CO_13CO_onH2O_03252016.csv',
			'TPD_12CO_13CO_10_1_onH2O_03292016.csv','TPD_12CO_13CO_10_1_onH2O_03292016.csv','TPD_12CO_13CO_10_1_onH2O_03292016.csv','TPD_12CO_13CO_10_1_onH2O_03292016.csv',
			'TPD2_12CO_13CO_10_1_onH2O_03292016.csv','TPD2_12CO_13CO_10_1_onH2O_03292016.csv','TPD2_12CO_13CO_10_1_onH2O_03292016.csv','TPD2_12CO_13CO_10_1_onH2O_03292016.csv',
			'TPD_12CO_13C18O_H2O_04072016.csv','TPD_12CO_13C18O_H2O_04072016.csv','TPD2_12CO_13C18O_H2O_04072016.csv','TPD2_12CO_13C18O_H2O_04072016.csv']
#n_col=[8,7,9,5]
n_col=[6,7,8,6,7,8,6,7,8,9,6,7,8,9,6,7,8,9,6,9,6,9]

t_in=np.zeros(len(file_temp))+21.
t_out=np.zeros(len(file_temp))+41.



for i in np.arange(len(file_temp)):
#for i in [1]:
	file_out='TPD_results/'+file_mass[i]+'_'+str(n_col[i])+'.dat'
	print file_out
	time_i,m_i=np.genfromtxt('Raw_data/'+file_mass[i],delimiter=',',unpack=True,skip_header=42,usecols=[1,n_col[i]])

	with open(file_mass[i], 'r') as infile:
	    data = infile.read()
	infile.close()
	my_list = data.splitlines()
	t_ini=time.strptime(str(my_list[2].split(',')[3]),"%I:%M:%S %p")
	t_ini_i=1e3*(t_ini.tm_hour*3600+t_ini.tm_min*60+t_ini.tm_sec)
	time_i=time_i+t_ini_i
	print t_ini_i

	book = xlrd.open_workbook(file_temp[i])
	sh = book.sheet_by_index(0)
	t_ini=time.strptime(str(sh.cell_value(1,1)),"%a %b %d %H:%M:%S %Z %Y")
	t_ini_t=1e3*(t_ini.tm_hour*3600+t_ini.tm_min*60+t_ini.tm_sec)
	time_t=np.arange(4,sh.nrows,dtype=np.float64)
	temp_t=np.arange(4,sh.nrows,dtype=np.float64)
	for rx in np.arange(4,sh.nrows):
	    time_t[rx-4]=np.float(sh.cell_value(rx,0))+t_ini_t
	    temp_t[rx-4]=np.float(sh.cell_value(rx,2))-4.5
	print t_ini_i,t_ini_t
	

	if t_ini_t > t_ini_i:
		print t_ini_t, t_ini_i
		print 'timt',time_t[0:4]
		coeff = np.polyfit( time_t[6:628], temp_t[6:628], 1 )
		print t_ini_i
		temp_ex=coeff[0]*(t_ini_i)+coeff[1]
		print temp_ex
		temp_t=np.insert(temp_t,0,temp_ex)
		time_t=np.insert(time_t,0,t_ini_i)
		print 'timt',time_t[0:4]

	if np.max(time_t) < np.max(time_i):
		#print t_ini_t, t_ini_i
		#print 'timt',time_t[0:4]
		coeff = np.polyfit( time_t[-20:], temp_t[-20:], 1 )
		#print t_ini_i
		temp_ex=coeff[0]*np.max(time_i)+coeff[1]
		print temp_ex
		temp_t=np.append(temp_t,temp_ex)
		time_t=np.append(time_t,np.max(time_i))
		print 'timt',time_t[-5:]

	interpfunc = interpolate.interp1d(time_t,temp_t, kind='linear')
	temp_i=interpfunc(time_i)

	ind=np.append(np.where(temp_i<t_in[i]),np.where(temp_i>t_out[i]))#+np.where(temp_i>t_out[i])
	#print ind
	#print np.shape(ind)
	#print len(temp_i)
	x_i=temp_i[ind]
	#print x_i
	y_i=m_i[ind]

	#spl = savitzky_golay(y_i,26,4)
	#s = interpolate.UnivariateSpline(x_i, y_i)(temp_i)
	#s=np.poly1d(np.polyfit(x_i, y_i, 7))(temp_i)
	#s=interp1d(x_i, y_i, kind='cubic')(temp_i)
	#func = scipy.interpolate.splrep(x_i, y_i, s=9)
	#yvals = scipy.interpolate.splev(temp_i, func, der=0)


	#m_i_bsln=m_i

	#plt.plot(temp_i,spl(temp_i),color='g')
	#plt.plot(x_i,spl,color='g')
	plt.plot(temp_i,m_i,linestyle='-',color='b')
	plt.scatter(x_i,y_i,color='r')
	#plt.plot(temp_i,s,color='r')
	print 'OK',i




	with open(file_out, 'w') as f:
		writer = csv.writer(f, delimiter='\t')
		writer.writerows(zip(temp_i,m_i))



#plt.xlim([20,65])
#plt.ylim([-1e4,4e5])
plt.show()