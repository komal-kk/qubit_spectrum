#!/usr/bin/env python
# coding: utf-8

# In[24]:


from qutip import*

import matplotlib.pyplot as plt
import numpy as np
import math
import warnings
import types
import scipy.sparse as sp
import scipy.linalg as la
from scipy.interpolate import BSpline, make_interp_spline
import qutip.settings as settings
from qutip import __version__
from qutip.fastsparse import fast_csr_matrix, fast_identity

from qutip.sparse import (sp_eigs, sp_expm, sp_fro_norm, sp_max_norm,
                          sp_one_norm, sp_L2_norm)
from qutip.dimensions import type_from_dims, enumerate_flat, collapse_dims_super
from qutip.cy.spmath import (zcsr_transpose, zcsr_adjoint)
import scipy.integrate as integrate

from scipy.integrate import dblquad



nn=13
h=6.63*10**(-34)
ghz=10**(9)
kb=1.38*10**(-23)
T=17*10**(-3)
bt=1/(kb*T)
e=1.6*10**(-19)

#ecj=12.52
ec=3.088
ej=17.45
#el=0.038

#cc=e**2/(2*ec*h*ghz)
#cj=e**2/(2*ecj*h*ghz)
phi0=h/(2*np.pi*e)
ic=ej/phi0
ctot=(2*e)**2/ec

ecsig= e**2/(2*(ctot)*h*ghz)

xi=(2*ec/ej)**(0.25)
#xi2=0.5*cc/(cc+cj)
#csig=cc+cj
#u=0.235


iota=complex(0,1)

a1=destroy(nn,(1-nn)/2)
c1=create(nn,(1-nn)/2)
#no4=tensor([identity(nn),identity(nn),identity(nn)])
no4=identity(nn)
a=iota*(c1-a1)
cp=a.eigenstates()
cp1=0
sp1=0

'''
nmat=0

print(cp[1])
for i in range(0,nn):
	nmat+=cp[0][i]*(cp[1][i]*(cp[1][i]).dag())
'''
#print(nmat,a)
for i in range(nn-1):
    cp1+=(1/2)*(cp[1][i]*(cp[1][i+1]).dag()+cp[1][i+1]*(cp[1][i]).dag())
    sp1+=(1/2)*(-iota)*(cp[1][i]*(cp[1][i+1]).dag()-cp[1][i+1]*(cp[1][i]).dag())
#count=0
#print(cp1)
phi1=(np.pi/2.0-(cp1+cp1**3/6.0+3*cp1**5/40.0+15*cp1**7/(7*48))-35*(cp1**9)/1152)

hamtot0=[]
hamtot1=[]
hamtot2=[]
e10=[]
'''
hamtot3=[]
hamtot4=[]
hamtot5=[]
hamtot6=[]
hamtot7=[]
hamtot8=[]
hamtot9=[]
hamtot10=[]
'''
hamtot=[]
for ng in np.arange(-1,1.5,0.5):
    ham=(1/2)*ec*(a-ng*no4)*(a-ng*no4)-ej*cp1
    hamtot.append(ham.eigenenergies())
    hamtot0.append((ham.eigenenergies())[0])
    hamtot1.append((ham.eigenenergies())[1])
    hamtot2.append((ham.eigenenergies())[2])
    e10.append((ham.eigenenergies())[1]-(ham.eigenenergies())[0])
    '''
    hamtot3.append((ham.eigenenergies())[3])
    hamtot4.append((ham.eigenenergies())[4])
    hamtot5.append((ham.eigenenergies())[5])
    hamtot6.append((ham.eigenenergies())[6])
    hamtot7.append((ham.eigenenergies())[7])
    hamtot8.append((ham.eigenenergies())[8])
    hamtot9.append((ham.eigenenergies())[9])
    hamtot10.append((ham.eigenenergies())[10])
    '''
    
tval1 = np.arange(-1,1.5,0.5)
#print(xi,hamtot1[0])
#print(hamtot1)
#Energy2 = (ham).eigenenergies()
#eigen=(ham).eigenstates()
''''''
xnew = np.linspace(tval1.min(), tval1.max(),50)
spl0 = make_interp_spline(tval1, hamtot0, k=3, bc_type='natural')  # type: BSpline
power_smooth0 = spl0(xnew)
spl1 = make_interp_spline(tval1, hamtot1, k=3, bc_type='natural')  # type: BSpline
power_smooth1 = spl1(xnew)
spl2 = make_interp_spline(tval1, hamtot2, k=3, bc_type='natural')  # type: BSpline
power_smooth2 = spl2(xnew)
spl10= make_interp_spline(tval1, e10, k=3, bc_type='natural')  # type: BSpline
power_smooth10 = spl10(xnew)
#print(hamtot)
plt.plot(xnew,power_smooth0,marker="o",markersize=1,linewidth=1)
plt.plot(xnew,power_smooth1,marker="o",markersize=1,linewidth=1)
plt.plot(xnew,power_smooth2,marker="o",markersize=1,linewidth=1)
plt.plot(tval1,hamtot0,marker="o",markersize=1,linewidth=1,color='green')
plt.plot(tval1,hamtot1,marker="o",markersize=1,linewidth=1,color='green')
plt.plot(tval1,hamtot2,marker="o",markersize=1,linewidth=1,color='green')
'''
plt.plot(tval1,hamtot0,marker="o",markersize=1,linewidth=1,color='green')
plt.plot(tval1,hamtot1,marker="o",markersize=1,linewidth=1,color='green')
plt.plot(tval1,hamtot2,marker="o",markersize=1,linewidth=1,color='green')
#plt.plot(xnew,power_smooth3,marker="o",markersize=1,linewidth=1)
plt.plot(tval1,hamtot3,marker="o",markersize=1,linewidth=1,color='green')
plt.plot(tval1,hamtot4,marker="o",markersize=1,linewidth=1,color='green')
plt.plot(tval1,hamtot5,marker="o",markersize=1,linewidth=1,color='green')
plt.plot(tval1,hamtot6,marker="o",markersize=1,linewidth=1,color='green')
plt.plot(tval1,hamtot7,marker="o",markersize=1,linewidth=1,color='green')
plt.plot(tval1,hamtot8,marker="o",markersize=1,linewidth=1,color='green')
plt.plot(tval1,hamtot9,marker="o",markersize=1,linewidth=1,color='green')
plt.plot(tval1,hamtot10,marker="o",markersize=1,linewidth=1,color='green')
'''
plt.xlabel("$n_g$",fontsize=18)
plt.ylabel("Energy (GHz)",fontsize=18)

plt.show()

plt.plot(xnew,power_smooth10,marker="o",markersize=1,linewidth=1)
plt.xlabel("$n_g$",fontsize=18)
plt.ylabel("Transition frequency (MHz)",fontsize=18)
plt.show()


# In[ ]:




