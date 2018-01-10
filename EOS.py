import numpy as np
import matplotlib.pyplot as plt
from scipy import constants 
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from scipy.optimize import minimize_scalar

#--------------------------------------------------------------------------------Constants--------------------------------------------------------------------------------
Ry=constants.physical_constants["Rydberg constant times hc in J"] #Rydberg (energy unit)
au=constants.physical_constants["Bohr radius"] # Bohr radius 

conv=(1e9*np.power(au[0],3)/Ry[0]) #conversion factor for bulk modulus in GPa


#--------------------------------------------------------------------------------EOS class start--------------------------------------------------------------------------------
class eos:
	"""
		Curve fit of energies and volumes from equations of state. Options: Murnaghan and Birch-Murnaghan.
	
															"""
	def __init__(self, name_file):
		self.file=name_file #name of the archive with the lattice parameters and energies
		

	def select(self):
		"""
			Read the files with information about the lattice parameters and energies.
													"""
		f=open(self.file, "r")
		lines = f.readlines()
		f.close()
		xdata=[] 
		ydata=[]
		calc={}
		self.n_atoms=int(lines[4].strip().split()[1]) #number of atoms in the unit cell
		for l in range(10, len(lines)):
			value=lines[l].strip().split()
			calc[float(value[0])**3/self.n_atoms]=float(value[1])
			xdata.append(float(value[0])**3/self.n_atoms)
			ydata.append(float(value[1]))
		self.volumes=np.array(xdata)
		self.energies=np.array(ydata)
		self.calc=calc
		return self.volumes, self.energies, self.calc, self.n_atoms

	
	def fit_2(self):
		"""
			Second order polinomial fit. It is the first guess to fit the curve with equations of state.
															"""
		self.fit_parameters=np.polyfit(self.volumes,self.energies, 2) #fitted coeficients
		self.fit_function=np.poly1d(self.fit_parameters) #fitted function
		V0=-self.fit_parameters[1]/(2*self.fit_parameters[0])
		B0=-self.fit_parameters[1]
		E0=self.fit_parameters[0]*V0**2+self.fit_parameters[1]*V0+self.fit_parameters[2]
		dB0=2. #assumption
		self.fit_guess=(E0, B0, dB0, V0)
		return self.fit_parameters, self.fit_function, self.fit_guess

		
	def Murnaghan(self, V, *args):
		"""
			Definition of the Murnaghan equation of state.
									"""
		E0, B0, dB0, V0 = args
		return E0 - B0*V0*conv/(dB0-1.) + (B0*V*conv/dB0) * (1. + ( (V0/V)**dB0 )/(dB0-1) )


	def fit(self, method='Murnaghan'):
		"""
			Choose of the equation of state: Murnaghan or Birch-Murnaghan. In this case, or method='Murnaghan' (standard) or method='Birch-Murnaghan'.
																			"""
		if method=="Murnaghan":
			popt, pcov=curve_fit(self.Murnaghan, self.volumes, self.energies, p0=self.fit_guess, method='trf') #fit for the Murnaghan equation
			self.E0=popt[0]
			self.B0=popt[1]
			self.dB0=popt[2]
			self.V0=popt[3]
			self.m_parameters=popt
			return self.E0, self.B0, self.dB0, self.V0, self.m_parameters
		elif method=="Birch-Murnaghan":
			print "Not implemented"
		else: 
			print "Just the Murnaghan and Birch-Murnaghan equations of state area available."
	

	def min(self, method='Murnaghan'):
		if method=='Murnaghan':
			sol= minimize_scalar(self.Murnaghan, args=tuple(self.m_parameters), method="brent")
			self.min_eos=sol.x
			return self.min_eos
		else:
			return "ERROR."
		

