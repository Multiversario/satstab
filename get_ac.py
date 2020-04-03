import numpy as np
import sys
from scipy.interpolate import griddata
from scipy.interpolate import interp2d

a_p = float(sys.arg[1])
m_p = float(sys.arg[2])
m_star = float(sys.arg[3]) 

R_H = a_p*(m_p/3.*m_star)**{1/3.}

e_p = float(sys.arg[4])
e_sat = float(sys.arg[5])
sub = int(sys.arg[6])# 0 = moon, 1 = submoon

if sub == 0:
  data = np.genfromtxt("Contour_moon.txt",delimiter=',',comments='#')
else:
    data = np.genfromtxt("Contour_submoon.txt",delimiter=',',comments='#')


X = data[:,0]
Y = data[:,1]
Z = data[:,2]

xi = np.arange(0.0,0.51,0.01)
yi = np.arange(0.0,0.51,0.01)
zi = griddata((X,Y),Z,(xi[None,:],yi[:,None]),method = 'linear',fill_value=0)

f = interp2d(xi, yi, zi, kind='linear')

print ("a_c = %1.3f AU" % (f(e_p,e_sat)[0] * R_H ))
