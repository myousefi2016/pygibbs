import matplotlib.pyplot as plt
import numpy as np
import interface
import debye
import EOS
import helmholtz_gruneisen 

#-------------------------------------------------------------------------Reading variables----------------------------------------------------------------------------------
T=interface.T
name=interface.name_file.strip().split('.')


#Debye class
material_debye=debye.Debye(interface.name_file)
material_debye.select_info()
material_debye.select_constants()
material_debye.VRH()
material_debye.sound_v()
theta0=int(material_debye.debye_temp())
cv0={}
for i in range(1, 2*theta0 + 1):
	cv0[i/theta0]=material_debye.c_v0(i)

plt.figure(1)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel(r'T/$\theta ^0$')
plt.ylabel('heat capacity (J/mol.K)')
plt.plot(*zip(*sorted(cv0.items())), linestyle='', linewidth=1.0 , color='r', marker='o', markersize=1.0, label="Data")
plt.gcf().subplots_adjust(bottom=0.15)
plt.tight_layout()
plt.savefig( './'+name[0].strip().split('.')[0]+'/'+name[0].strip().split('.')[0] + '_cv_0.png')


#EOS class
material_eos = EOS.eos(interface.name_file)
material_eos.select()
material_eos.fit_2()
material_eos.fit()
material_eos.min()
xfit=np.linspace(min(material_eos.volumes), max(material_eos.volumes), 100)		
plt.figure(2) 
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel(r'Volume (a.u.$ ^3$)')
plt.ylabel('Total Energy (Ry)')
plt.plot(*zip(*sorted(material_eos.calc.items())), linestyle=':', linewidth=1.0 , color='b', marker='o', markersize=3.0, label="Data")
plt.plot(xfit, material_eos.fit_function(xfit), linestyle='-', linewidth=1.0 , color='r', marker='', markersize=3.0, label="Polynomial fit")
plt.plot(xfit, material_eos.Murnaghan(xfit, *material_eos.m_parameters), linestyle='-', linewidth=1.0 , color='g', marker='', markersize=3.0,label="Murnaghan fit")
plt.gcf().subplots_adjust(bottom=0.15)
plt.tight_layout()
plt.legend(loc=1)
plt.savefig( './'+name[0].strip().split('.')[0]+'/'+name[0].strip().split('.')[0] + '_fit.png')
