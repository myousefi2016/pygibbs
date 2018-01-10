import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import interface
import EOS
import debye
import helmholtz_gruneisen 

#-------------------------------------------------------------------------Reading variables----------------------------------------------------------------------------------
T=interface.T
name=interface.name_file.strip().split('.')

#---------------------------------------------------------------------------Reading file-------------------------------------------------------------------------------------
f_data=open(interface.name_file, 'r')
lines=f_data.readlines()
f_data.close()

#-------------------------------------------------------------Creating the directory and the file--------------------------------------------------------------------------------
shutil.rmtree(name[0])
os.mkdir(name[0])
f = open('./'+name[0]+'/'+name[0]+'.txt', 'w')

#-------------------------------------------------------------------------Basic information-----------------------------------------------------------------------------------
f.write("1. INFORMATIONS:" +"\n\n")
f.write("	1.1. Material: "+name[0] +"\n")
for i in range(1,7):
	f.write( '%45s' % (lines[i]))
f.write("\n")
f.write("	1.2. Method: Murnaghan equation of state."+"\n")
C=lines[8].strip().split()
f.write("	1.3. Elastic constants:"+"\n")
f.write("			C11=" + str(C[0]) + " GPa" + '\n' )
f.write("			C12="+ str(C[1] ) + " GPa"+ '\n' )
f.write("			C44="+ str(C[2] ) + " GPa"+ '\n' )

#Implementation of the equation of state (EOS) class
f.write('\n')

f.write("2. RESULTS FROM THE FIT OF MURNAGHAN EQUATION OF STATE:" + '\n\n')
material_eos = EOS.eos(interface.name_file)
f.write("	2.1. Experimental data:"+"\n")
f.write('%18s  %15s  %12s' % ('a (au)', 'Volume (au^3)', 'Energy (Ry)') + '\n' )
for j in range(10, len(lines)):
	value=lines[j].strip().split()
	f.write('%18s  %15s  %12s' % (str(value[0]), str(float(value[0])**3/material_eos.select()[3]), str(value[1])) + '\n' )
f.write('\n')
coef_pol=material_eos.fit_2()[0]
f.write( "	2.2. Second order fit: " + '\n')
f.write("			a=" + str(coef_pol[0]) + '\n' )
f.write("			b="+ str(coef_pol[1] ) + '\n' )
f.write("			c="+ str(coef_pol[2] ) + '\n\n' )
coef_guess=material_eos.fit_2()[2]
f.write("	2.3. Initial guess for Murnaghan parameters:" + '\n')
f.write("			E0=" + str(coef_guess[0]) + " Ry" + '\n' )
f.write("			B0="+ str(coef_guess[1] ) + " GPa" + '\n' )
f.write("			dB0="+ str(coef_guess[2] ) + '\n' )
f.write("			V0="+ str(coef_guess[3] ) + " au^3" + '\n\n' )
coef_eos=material_eos.fit()[4] 
f.write("	2.4. Murnaghan parameters:" + '\n')
f.write("			E0=" + str(coef_eos[0]) + " Ry" + '\n' )
f.write("			B0="+ str(coef_eos[1] ) + " GPa"+ '\n' )
f.write("			dB0="+ str(coef_eos[2] ) + '\n' )
f.write("			V0="+ str(coef_eos[3] ) + " au^3"+ '\n\n' )
f.write("	2.5. Equilibrium volume from Murnaghan equation: "+str(material_eos.min()) +' au^3'+ '\n')


#Implementation of the Debye model from Debye class
f.write('\n')
f.write("3. RESULTS FROM THE DEBYE MODEL:" + '\n\n')
material_debye=debye.Debye(interface.name_file)
material_debye.select_info()
material_debye.select_constants()
moduli=material_debye.VRH()
f.write("	3.1. Poisson coefficient and the elastic moduli:" +'\n')
f.write("				nu=" + str(moduli[0]) + '\n' )
f.write("				S="+ str(moduli[1] ) +" GPa" + '\n' )
f.write("				L="+ str(moduli[2] ) +" GPa" + '\n' )
f.write("				B="+ str(moduli[3] ) +" GPa" + '\n\n' )
material_debye.sound_v()
f.write("	3.2. Debye temperature: "+ str(material_debye.debye_temp())+' K' +'\n\n')
f.write("	3.3. Specific heat (constant volume): "+ str(material_debye.c_v0(T)) +' J/molK' +'\n')

#Implementation of the Debye-Gruneisen model 
f.write('\n\n')
f.write("4. RESULTS FROM THE DEBYE-GRUNEISEN MODEL AND HELMHOLTZ FREE ENERGY:" + '\n\n')
material_HG=helmholtz_gruneisen.Gruneisen_Helmholtz(material_eos.fit()[4], T, material_debye.debye_temp(), material_debye.select_info()[0])
f.write("	4.1. Gruneisen parameter: "+str(material_HG.gru_parameter()) +'\n\n')
f.write("	4.2. Minimum volume from Helmholtz free energy: " +str(material_HG.minimum())+ ' au^3'+'\n\n')
f.write("	4.3. Helmholtz free energy: " +str(material_HG.F(material_HG.minimum()))+ ' Ry'+'\n\n')
f.write("	4.4. Debye temperature:"+str( material_HG.TD ) + ' K'+'\n\n')
f.write("	4.5. Specific heat (constant volume):"+str( material_HG.cv())+' J/molK' +'\n\n')
f.write("	4.6. Bulk modulus: " +str(material_HG.bulk_m())+ ' GPa'+'\n\n')
f.write("	4.7. Expansion coefficient: " +str(material_HG.alpha())+ ' K^-1'+'\n\n')
f.write("	4.8. Specific heat (constant pressure): " +str(material_HG.cp())+' J/molK' +'\n\n')



#-----------------------------------------------------Plots-------------------------------------------------------------------------------
theta0=material_debye.debye_temp()
#Debye class
cv0={}
for j in range(1, int(2*theta0)):
	cv0[j/theta0]=material_debye.c_v0(j)
plt.figure(1)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel(r'T/$\theta _D$', fontsize=18)
plt.ylabel('specific heat (J/mol.K)', fontsize=18)
plt.plot(*zip(*sorted(cv0.items())), linestyle='-', linewidth=1.0 , color='b', marker='', markersize=1.0, label="Data")
plt.gcf().subplots_adjust(bottom=0.15)
plt.tight_layout()
plt.savefig( './'+name[0].strip().split('.')[0]+'/'+name[0].strip().split('.')[0] + '_cv_0.png')

#EOS class
xfit=np.linspace(min(material_eos.volumes), max(material_eos.volumes), 100)		
plt.figure(2) 
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel(r'Volume (au$ ^3$)', fontsize=18)
plt.ylabel('Total Energy (Ry)', fontsize=18)
plt.plot(*zip(*sorted(material_eos.calc.items())), linestyle=':', linewidth=1.0 , color='b', marker='o', markersize=3.0, label="Data")
plt.plot(xfit, material_eos.fit_function(xfit), linestyle='-', linewidth=1.0 , color='r', marker='', markersize=3.0, label="Polynomial fit")
plt.plot(xfit, material_eos.Murnaghan(xfit, *material_eos.m_parameters), linestyle='-', linewidth=1.0 , color='g', marker='', markersize=3.0,label="Murnaghan fit")
plt.gcf().subplots_adjust(bottom=0.15)
plt.tight_layout()
plt.legend(loc=1)
plt.savefig( './'+name[0].strip().split('.')[0]+'/'+name[0].strip().split('.')[0] + '_fit.png')

#Gruneisen-Helmholtz class
F={}
cp={}
cv={}
for k in range(1, int(2*theta0)):
	material_HG2=helmholtz_gruneisen.Gruneisen_Helmholtz(material_eos.fit()[4], k, material_debye.debye_temp(), material_debye.select_info()[0])
	material_HG2.gru_parameter()
	material_HG2.minimum()
	material_HG2.TD
	material_HG2.Dfunc
	material_HG2.Fvib(material_HG2.minimum())
	material_HG2.Murnaghan(material_HG2.minimum())
	material_HG2.F(material_HG2.minimum())
	#F[material_HG2.minimum()]= material_HG2.F(material_HG2.minimum())
	material_HG2.bulk_m()
	cv[k/theta0]=material_HG2.cv()
	material_HG2.alpha()
	cp[k/theta0]= material_HG2.cp()

V=np.linspace(0.85*material_HG.minimum(), 1.2*material_HG.minimum(), num=100)
for l in V:
	F[l]=material_HG.F(l)

plt.figure(3)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel(r'V (au$^3$)', fontsize=18)
plt.ylabel('Energy (Ry)', fontsize=18)
plt.plot(V, material_HG.Murnaghan(V), linestyle='-', linewidth=1.0 , color='r', marker='', markersize=2.0, label="Murnaghan")
plt.plot(*zip(*sorted(F.items())), linestyle='-', linewidth=1.0 , color='b', marker='', markersize=2.0, label="Helmholtz free energy")
plt.gcf().subplots_adjust(bottom=0.15)
plt.tight_layout()
plt.legend(loc=2)
plt.savefig( './'+name[0].strip().split('.')[0]+'/'+name[0].strip().split('.')[0] + '_F.png')

plt.figure(4)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel(r'T/$\theta _D$', fontsize=18)
plt.ylabel('specific heat (J/mol.K)', fontsize=18)
plt.plot(*zip(*sorted(cp.items())), linestyle='-', linewidth=1.0 , color='b', marker='', markersize=3.0, label="Constant pressure")
plt.plot(*zip(*sorted(cv.items())), linestyle='-', linewidth=1.0 , color='r', marker='', markersize=3.0, label="Constant volume")
plt.gcf().subplots_adjust(bottom=0.15)
plt.tight_layout()
plt.legend(loc=4)
plt.savefig( './'+name[0].strip().split('.')[0]+'/'+name[0].strip().split('.')[0] + '_cv_cp.png')


f2=open("diamante.dat", "r")
lines2 = f2.readlines()
f2.close()

cvexp={}

for CV in lines2[1:]:
	v=CV.strip().split()
	cvexp[float(v[0])/theta0]=float(v[1])

plt.figure(5)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel(r'T/$\theta _D$', fontsize=18)
plt.ylabel('specific heat (J/mol.K)', fontsize=18)
plt.plot(*zip(*sorted(cp.items())), linestyle='-', linewidth=1.0 , color='b', marker='', markersize=3.0, label="Constant pressure")
plt.plot(*zip(*sorted(cvexp.items())), linestyle='-', linewidth=1.0 , color='g', marker='', markersize=3.0, label="Experimental data")
plt.gcf().subplots_adjust(bottom=0.15)
plt.tight_layout()
plt.legend(loc=4)
plt.savefig( './'+name[0].strip().split('.')[0]+'/'+name[0].strip().split('.')[0] + '_exp.png')
