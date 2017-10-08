from scipy.optimize import leastsq
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants 
import scipy
from scipy.optimize import curve_fit

#---------------------------------------------------------------------------Constants-----------------------------------------------------------------
Ry=constants.physical_constants["Rydberg constant times hc in J"]
au=constants.physical_constants["Bohr radius"]

conv=(np.power(10,9)*np.power(au[0],3)/Ry[0]) #conversion factor

#-------------------------------------------------------------------Data + second order fit-----------------------------------------------------------
calc={}
xdata=[]
ydata=[]

name=raw_input("Name of the file: ")
n=input("Number of atoms in the cell: ")

f=open(name, "r")
lines = f.readlines()
f.close()

for i in lines:
	v=i.strip().split()
	calc[float(v[0])**3/n]=float(v[1])
	xdata.append(float(v[0])**3/n)
	ydata.append(float(v[1]))

z = np.polyfit(xdata,ydata, 2) #fitted coeficients
p = np.poly1d(z) #fitted function

xfit=np.linspace(min(xdata), max(xdata), 100)

#----------------------------------------------------------------------Murnaghan function--------------------------------------------------------------------------
def Murnaghan(vols, E0, B0, dB, V0):
	E = E0 -B0*V0*conv/(dB-1) + (B0*vols*conv/dB) * (1 + ( (V0/vols)**dB )/(dB-1) )
	return E


#----------------------------------------------------------------------Murnaghan parameters------------------------------------------------------------------------

#np.array for the x and y data
vols=np.array(xdata)           
energies=np.array(ydata)

#Parameters obtained from the polynomial fit
V0=-z[1]/(2*z[0])
B0=-z[1]
E0=z[0]*V0**2+z[1]*V0+z[2]
dB=2 #assumption
p0 = scipy.array([E0, B0, dB, V0])

popt, pcov = curve_fit(Murnaghan, vols, energies, p0, method='trf') #fit for the Murnaghan equation

# 'trf': when the number of observations is less than the number of variables

print "Final parameters (E0, B0, dB, V0): ", popt


#-------------------------------------------------------------------------Graphic----------------------------------------------------------------------------------
plt.figure()

plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)

plt.xlabel('Volume (a.u. ^3)')
plt.ylabel('Energia total (Ry)')

plt.plot(*zip(*sorted(calc.items())), linestyle=':', linewidth=1.0 , color='b', marker='o', markersize=3.0, label="data")
plt.plot(xfit, p(xfit), linestyle='-', linewidth=1.0 , color='r', marker='', markersize=3.0, label="polynomial fit")
plt.plot(xfit, Murnaghan(xfit, *popt), linestyle='-', linewidth=1.0 , color='g', marker='', markersize=3.0,label="Murnaghan fit")

plt.gcf().subplots_adjust(bottom=0.15)
plt.tight_layout()

plt.legend(loc=1)

plt.savefig('fit.png')

plt.show()
