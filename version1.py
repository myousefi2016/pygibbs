from scipy import constants 
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


#-----------------------------------------------Physical Constants-----------------------------------------------------
hbar=constants.hbar
k=constants.k
pi=constants.pi
N=constants.Avogadro
u=constants.physical_constants['unified atomic mass unit'] 


#-------------------------------------------Calculating the elastic moduli---------------------------------------------
option = input("Do you want to insert a text file (1) with the elastic constants or each constant (2)? ")

if option==1:
	name = raw_input("Insert the name of the text file: ")                     #txt file with the elastic constants
	f=open(name, "r")
	lines = f.readlines()
	f.close()
	C=[]
	for i in lines:
		I=i.strip().split("\n")
		C.append(float(I[0]))
else: 
	print "Insert the three elastic constants:"                                #tuple with the elastic constants
	C=input("c11, c12, c44: ")


def voigt(C):                                                                      #Voigt model
	Bv=1.0/3*(C[0]+2*C[1])
	Sv=1.0/5*(C[0]-C[1]+3*C[2])
	nu_v=(3*Bv-2*Sv)/(2.0*(3*Bv+Sv))
	voigt=(Bv, Sv, nu_v)
	return voigt
modulusV=voigt(C)


def reuss(C):                                                                      #Reuss model
	Br=1.0/3*(C[0]+2*C[1])
	Sr=(5.0*(C[0]-C[1])*C[2])/(4*C[2]+3*(C[0]-C[1]))
	nu_r=(3*Br-2*Sr)/(2.0*(3*Br+Sr))
	reuss=(Br, Sr, nu_r)
	return reuss
modulusR=reuss(C)


def VRH(reuss, voigt):                                                             #Hill model
	B=(modulusV[0]+modulusR[0])/2
	nu=(modulusV[2]+modulusR[2])/2
	L=3.0*B*(1-nu)/(1+nu)
	S=(3.0*(1-2*nu)*B)/(2*(1+nu))
	VRH=(nu, B, L, S)
	return VRH
modulus= VRH(modulusV, modulusR)


print "nu=",modulus[0]                                                             #Elastic moduli calculated
print "L=", modulus[2],"GPa"
print "S=", modulus[3],"GPa"
print "B=", modulus[1],"GPa"


#--------------------------------------Calculating the Debye temperature-----------------------------------------------------
V = input("Insert the molar volume [m^3/mol]: ")
M = constants.m_u*input("Insert the atomic weight [u]: ")

rho=float(N*M/V)                                                                   #density
vl=np.sqrt(10**9*modulus[2]/rho)                                                   #sound velocities
vs=np.sqrt(10**9*modulus[3]/rho)
const0=1.0/vl**3+2.0/vs**3
vd=np.power(3.0/const0,1.0/3)

def debye1(rho, V, vd):                                                            #Debye temperature
	const = (1.0/V)*6*N*np.power(pi,2)
	debye = (hbar/k)*vd*np.power(const, 1.0/3)
	return debye
theta=debye1(rho, V, vd)

print "The Debye temperature is:", theta


#----------------------------------------Calculating the specific heat---------------------------------------------------------
T=input("Insert the maximum temperature [K]: ")                                    #Range of temperatures
temp=np.linspace(1,T, num=T)


def integrand(x):                                                                  #Defining the integrand
	return np.power(x,4)*np.exp(x)/np.power(np.exp(x)-1,2)


def heat(temp, theta):                                                             #Specific heat
	c_v={}
	for i in temp:
		r1=theta/i
		r2=i/theta
		I=quad(integrand, 0, r1)
		c=I[0]*9.0*N*k*np.power(r2,3)
		c_v[r2]=c
	return c_v
c_v=heat(temp,theta)

                                                                  
plt.figure(1)                                                                      #Plot of the specific heat versus T/theta
plt.plot(*zip(*sorted(c_v.items())), linestyle='-', linewidth=1.0, color='r', marker='.', markersize=3.0)
plt.savefig('c_v.png')
plt.show()



















