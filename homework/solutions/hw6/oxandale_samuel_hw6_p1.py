import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
import scipy as sp
from scipy import optimize

plt.ioff()
exp = np.exp

def Q(rho_e, rho_h) :
    
    return rho_e + rho_e**2*(exp(-1.0/rho_e)-1.0) + \
          rho_h + rho_h**2*(exp(-1.0/rho_h)-1.0)
          
          
def sig_Q(rho_e, rho_h) :
    
    a = rho_e**2 + 2.*rho_e**3*(exp(-1.0/rho_e)-1) + \
       0.5*rho_e**3*(1-exp(-2.0/rho_e))
       
    b = rho_h**2 + 2.*rho_h**3*(exp(-1.0/rho_h)-1) + \
       0.5*rho_h**3*(1-exp(-2.0/rho_h))
       
    c = 2.*rho_e*rho_h + 2.*rho_e**2*rho_h*(exp(-1.0/rho_e)-1) + \
       2.*rho_h**2*rho_e*(exp(-1.0/rho_h)-1)
       
    d = 2.*(rho_e*rho_h)**2/(rho_e-rho_h)*(exp(-1.0/rho_e)-exp(-1.0/rho_h))
    
    return np.sqrt( a+b+c+d-Q(rho_e,rho_h)**2)


def R(rho_e, rho_h) :
    
    return 100*sig_Q(rho_e, rho_h)/Q(rho_e, rho_h)


n=100
H = np.logspace(-2, 2, n) 
E = np.logspace(-2, 2, n)

# Contour levels
levels = np.array([0.1, 0.2, 0.5, 1., 2., 5., 10., 15., 20., 30., 40.])

H, E = np.meshgrid(H, E, sparse=False, indexing='ij')

res = R(E, H)

# Init figure
plt.figure(1, figsize=(8,8))
plt.contour(np.log10(E),np.log10(H), res, levels, colors='k')

# Axis labels
plt.ylabel('Hole Extraction Factor')
plt.xlabel('Electron Extaction Factor')

# Axis major ticks
plt.xticks([-1, 0, 1, 2], [0.1, 1, 10, 100])
plt.yticks([-1, 0, 1, 2], [0.1, 1, 10, 100])

# Set locations of minor ticks
minor_ticks = np.array([])

for i in np.arange(-2, 2) :
    minor_ticks = np.append(minor_ticks, np.arange(10.**i, 10.**(i+1), 10.**i))
    
minor_ticks = np.log10(minor_ticks)
    
plt.gca().yaxis.set_minor_locator(FixedLocator(minor_ticks))
plt.gca().xaxis.set_minor_locator(FixedLocator(minor_ticks))

# Contour line labels
font = 8
for x in [40., 30., 20., 15.] :
    xloc = np.log10(sp.optimize.newton(lambda rho_e: x - R(rho_e, 100), 0.1))
    plt.text(xloc, 2.01, '{}%'.format(int(x)), fontsize=font, horizontalalignment='center')

for y in levels[:3] :
    yloc = np.log10(sp.optimize.newton(lambda rho_h: y - R(100, rho_h), 100-100*y))
    plt.text(2.02, yloc, '{}%'.format(y), fontsize=font, verticalalignment='center')
    
for y in levels[4:6] :
    yloc = np.log10(sp.optimize.newton(lambda rho_h: y - R(100, rho_h), 5))
    plt.text(2.02, yloc, '{}%'.format(y), fontsize=font, verticalalignment='center')
    
for y in levels[7:] :
    yloc = np.log10(sp.optimize.newton(lambda rho_h: y - R(100, rho_h), 0.01))
    plt.text(2.02, yloc, '{}%'.format(y), fontsize=font, verticalalignment='center')

plt.savefig(’new_contour.png’)