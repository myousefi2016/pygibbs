import numpy as np
from scipy import constants 
from scipy.optimize import curve_fit
from scipy.optimize import minimize_scalar
from scipy.optimize import leastsq
from scipy.integrate import quad
from EOS import eos
#-----------------------------------------------------------------------------------------Constants----------------------------------------------------------------------------------
Ry=constants.physical_constants["Rydberg constant times hc in J"]
au=constants.physical_constants["Bohr radius"]
hbar=constants.hbar
k=constants.k
pi=constants.pi
N=constants.Avogadro
u=constants.physical_constants['unified atomic mass unit'] 
R=constants.R

conv=(1e9*np.power(au[0],3)/Ry[0]) #conversion factor for bulk modulus in GPa

#--------------------------------------------------------------------------------Helmholtz class start--------------------------------------------------------------------------------

class Gruneisen_Helmholtz:


	def __init__(self, parameters, T, theta, Vm):
		self.eos_parameters=parameters #fit parameters for Murnaghan or Birch murnaghan equation of state (EOS) (E0, B0, dB0, V0)
		self.T=T #temperature (debye) 
		self.Vm=Vm #molar volume (debye-gruneisen)
		self.E0=parameters[0]
		self.B0=parameters[1]
		self.dB0=parameters[2]
		self.V0=parameters[3]
		self.theta0=theta


	def integrand(self, x):
		return x**3/(np.exp(x)-1.)


	def gru_parameter(self):
		"""
			Gruneisen parameter.
									"""
		#if self.T<=self.theta0:
		self.g= -1. + (1. + self.dB0)/2
		#else:
		#	self.g= -2./3 + (1. + self.dB0)/2
		return self.g


	def dg_temp(self, V):
		"""
			Debye temperature from Debye-Gruneisen model.
													"""
		self.TD = self.theta0*(self.V0/V)**self.g
		return self.TD

	
	def D_function(self, V):
		"""
			Debye function.
													"""
		I=quad(self.integrand, 0, self.dg_temp(V)/self.T)
		self.Dfunc = I[0]*3*(self.dg_temp(V)/self.T)**(-3)
		return self.Dfunc	


	def Fvib(self, V):
		"""
			Vibrational Helmholtz free energy.
									"""
		self.F_vib= 9./8*k/Ry[0]*self.dg_temp(V) + k/Ry[0]*self.T*( 3*np.log(1-np.exp(-self.dg_temp(V)/self.T )) - self.D_function(V))
		return self.F_vib 
	
	def Murnaghan(self, V):
		"""
			Definition of the Murnaghan equation of state.
									"""
		self.E= self.E0 - self.B0*self.V0*conv/(self.dB0-1.) + (self.B0*V*conv/self.dB0) * (1. + ( (self.V0/V)**self.dB0 )/(self.dB0-1) )
		return self.E 

	def F(self, V, method='Murnaghan'):
		"""
			Helmholtz free energy.
									"""
		self.HF = self.Fvib(V) + self.Murnaghan(V)
		return self.HF	

	def minimum(self):
		"""
			Minimization of the Helmholtz free energy.
									"""
		F_min=minimize_scalar(self.F, bracket=(0.9*self.V0, self.V0), method="brent")
		self.V=F_min.x
		return self.V

	
	def bulk_m(self):
		"""
			Calculation of the bulk modulus for the temperature T given.
									"""
		B=self.B0*1e9*(self.V0/self.V)**self.dB0 + 1/self.Vm*(9./8*self.g*(1+self.g)*R*self.dg_temp(self.V) + 3*(1-3*self.g)*R*self.T*self.D_function(self.V) + 9.*self.g**2*R*self.dg_temp(self.V)/(np.exp(self.dg_temp(self.V)/self.T)-1) )
		self.B=B/1e9
		return self.B


	def integrand2(self,x):
		return np.power(x,4)*np.exp(x)/np.power(np.exp(x)-1.,2)


	def cv(self):
		"""
			Specific heat for constant volume from Debye-Gruneisen model.
													"""
		I=quad(self.integrand2, 0, self.TD/self.T)
		self.Cv =I[0]*9.*N*k*(self.T/self.TD)**3
		return self.Cv


	def alpha(self):
		"""
			Thermal expansion coefficient.
									"""
		self.a=self.g*self.Cv/(self.V*au[0]**3*self.B*1e9)/N
		return self.a

	
	def cp(self):
		"""
			Calculation of the specific heat for constant pressure.
										"""
		self.c_p=self.Cv + (self.a)**2*self.B*1e9*self.T*self.V*au[0]**3*N
		return self.c_p


















