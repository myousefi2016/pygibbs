import numpy as np
import matplotlib.pyplot as plt
from scipy import constants 
from scipy.integrate import quad


#-----------------------------------------------------------------------Constants-----------------------------------------------------------------------------------------------
hbar=constants.hbar
k=constants.k
pi=constants.pi
N=constants.Avogadro
u=constants.physical_constants['unified atomic mass unit'] 
Ry=constants.physical_constants["Rydberg constant times hc in J"]


#--------------------------------------------------------------------------------Debye class start--------------------------------------------------------------------------------
class Debye:
	"""
		Implementation of the Debye model and of the Debye-Gruneisen model for the calculation of the Debye temperature and specific heat for constant volume.
	
																								"""
	
	def __init__(self, name_file):
		self.name_file=name_file #file with the elastic constants
	
	
	def select_info(self):
		"""
			Read the files with information about the material.
																		"""
		f=open(self.name_file, "r") 
		lines= f.readlines()
		f.close()
		self.rho= 1000*float(lines[5].strip().split()[1]) #density
		self.M= constants.m_u*float(lines[6].strip().split()[1]) #atomic mass
		self.Vm=float(N*self.M/self.rho) #molar volume
		return self.Vm, self.M, self.rho
	
	
	def select_constants(self):
		"""
			Read the files with information about the elastic constants for cubic materials (3 independent elastic constants).
																		"""
		f=open(self.name_file, "r") 
		lines= f.readlines()
		f.close()
		constants=lines[8].strip().split()
		if len(constants)==3:
			"""
				Cubic material. 
						"""
			self.C11=float(constants[0])
			self.C12=float(constants[1])
			self.C44=float(constants[2])
			return self.C11, self.C12, self.C44 #GPa units 
		else:
			return "Non cubic material is not supported."


	def VRH(self):
		"""
			Calculation of the elastic moduli and Poisson coeficient.
											"""
		self.B=1./3*(self.C11 + 2*self.C12)
		Sv=1./5*(self.C11 - self.C12 + 3*self.C44)
		Sr=5.*(self.C11 - self.C12)*self.C44/(4*self.C44 + 3*(self.C11 - self.C12))
		nuv=1./2*(3*self.B - 2*Sv)/(3*self.B + Sv)
		nur=1./2*(3*self.B - 2*Sr)/(3*self.B + Sr)
		self.nu=(nuv + nur)/2
		self.S=(Sv + Sr)/2
		self.L=3*self.B*(1 - self.nu)/(1 + self.nu)
		return self.nu, self.S, self.L, self.B

	def sound_v(self):
		"""
			Calculation of the Debye sound velocity.
									"""
		vl=np.sqrt(1e9*self.L/self.rho)
		vs=np.sqrt(1e9*self.S/self.rho)
		self.vd=(1./3*(1./vl**3 + 2./vs**3))**(-1./3)
		return self.vd

	def debye_temp(self):
		"""
			Calculation of the Debye temperature from Debye model.
										"""
		self.theta0=self.vd*hbar/k*(6.*N*pi**2/self.Vm)**(1./3)
		return self.theta0

	def integrand(self,x):
		return np.power(x,4)*np.exp(x)/np.power(np.exp(x)-1.,2)

	def c_v0(self, T):
		"""
			Specific heat at constant volume with Debye model.
										"""
		self.T=T
		I=quad(self.integrand, 0, self.theta0/self.T)
		self.cv0 =I[0]*9.*N*k*(self.T/self.theta0)**3
		return self.cv0


