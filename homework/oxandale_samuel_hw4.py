#!/usr/bin/env python3

import re
import matplotlib.pyplot as plt
import scipy as sp
from scipy import stats
from scipy import integrate
import numpy as np

plots = True

# =============================================================================
# Problem 1
# =============================================================================

matches = ['pit', 'spot', 'spate', 'slap two', 'respite']
not_matches = ['pt', 'Pot', 'peat', 'part']

pattern = '\w*p[\w\s]t\w*'

for s in matches :
    if not re.match(pattern, s):
        print(s + ' does not match the pattern and should.')
        
for s in not_matches :
    if re.match(pattern, s) :
        print(s + 'matches the pattern and should not.')
        
        
# =============================================================================
# Problem 2
# =============================================================================

f = open('regex_sample_mcnp.txt', 'r')

s = f.read()

f.close()

tally = re.findall('1tally\s+(\d+)', s)

energy = re.findall('energy\s+\n\s+([0-9E\+-\.]{10}).+\n\s+([0-9E\+-\.]{10}).+\n\s+([0-9E\+-\.]{10}).+\n', s)
energy = [[float(i) for i in j] for j in energy]

value = re.findall('energy\s+\n(?:\s+.+([0-9E\+-\.]{11}).+\n)(?:\s+.+([0-9E\+-\.]{11}).+\n)(?:\s+.+([0-9E\+-\.]{11}).+\n)', s)
value = [[float(i) for i in j] for j in value]

sigma = re.findall('energy\s+\n(?:.+([0-9\.]{6})\n)(?:.+([0-9\.]{6})\n)(?:.+([0-9\.]{6})\n)', s)
sigma = [[float(i) for i in j] for j in sigma]

    

d = {}
for i in range(len(tally)) :
    d[int(tally[i])] = {'energy': energy[i], 'value': value[i], 'sigma': sigma[i]}
    
    
# =============================================================================
# Problem 3
# =============================================================================

# *****************************************************************************
#    Read in Data
# *****************************************************************************    
f = open('splits.txt', 'r')

s = f.read()

f.close()
# Splits are given in the following format :
# Distance, km  :   Total time, (h):mm:ss, : Time from last split, m:ss
# i.e. 02 : 5:41	(2:51)  --->  At 2km, 5:41 has elapsed in the race, 2:51 since last split.
splits = re.findall('(\d+(?:.\d+)?)\s+:\s+(\d?:?\d+:\d+)\s+\((\d+(?::|.)\d+)\)', s)

distance = []
time = []
split = []

for data in splits :
    # Adds value to distance array
    distance.append(float(data[0]))
    
    # Adds value to time array
    t_total_split = data[1].split(':')
    # If format is mm:ss...
    if len(t_total_split) == 2 :
        # Convert to seconds
        t_total = 60*float(t_total_split[0]) + float(t_total_split[1])
    # Format is h:mm:ss...
    else :
        # Converts to seconds
        t_total = 60*60*float(t_total_split[0]) + 60*float(t_total_split[1]) + float(t_total_split[2])
    time.append(t_total)
    
    # Adds value to split array
    t_total_split = data[2].split(':')
    # If format is m:ss
    if len(t_total_split) == 2 :
        # Convert to seconds
        t_total = 60*float(t_total_split[0]) + float(t_total_split[1])
    # Format is ss.s
    else :
        t_total = float(t_total_split[0])
    split.append(t_total)
 
distance = np.array(distance)
time = np.array(time)
split = np.array(split)    
# *****************************************************************************
#  Least Squares
# ***************************************************************************** 
    
m1, b1, r, p, sderr = sp.stats.linregress(distance, time)

# Last value thrown out as it is not the split of a full km 
coeff = np.polyfit(distance[:-1], split[:-1], 4)  

time_fit = m1*distance + b1

split_fit_poly = np.poly1d(coeff)(distance)

# *****************************************************************************
#    Minimax
# ***************************************************************************** 

def fun(coeff, x, y_act, order):
    a, b, c, d, e = coeff
    # There were no existing models of the data that I could fit my equation
    # to, but had I found one, I could have put the general form of it here
    # and used this method for least squares and minimax, adjusting the order
    # to 2 and inf respectively.
    y = a*x**4 + b*x**3 + c*x**2 + d*x**1 + e
    return np.linalg.norm(y - y_act, order)

coeff2 = sp.optimize.minimize(fun, [1,-1,1,-1,170] , (distance[:-1], split[:-1], np.inf))

split_fit_minimax = np.poly1d(coeff2.x)(distance)

#if plots :
#    plt.figure(-1)
#    plt.plot(distance, np.array(time)/60, 'r.', distance, time_fit/60, '-b')
#    plt.title("Time Elapsed vs. Distance for \nEliud Kipchoge's Sub-2 Marathon Attempt")
#    plt.xlabel('Distance, km')
#    plt.ylabel('Time Elapsed, min')
#    plt.legend(['Actual Data', 'Least Squares Fit'])

# *****************************************************************************
#    Plots
# ***************************************************************************** 

if plots :
    plt.figure(0)
    # Last value thrown out as it is not the split of a full km
    plt.plot(distance[:-1], split[:-1], 'r.', distance[:-1], split_fit_poly[:-1], '-b', distance[:-1], split_fit_minimax[:-1], '--r')
    plt.title("Splits for each km for \nEliud Kipchoge's Sub-2 Marathon Attempt")
    plt.xlabel('Distance, km')
    plt.ylabel('Time Elapsed, s')
    plt.legend(['Actual Data', 'Least Squares Fit', 'Minimax Fit'])


# =============================================================================
# Problem 4
# =============================================================================

# 4.1
def fun1(y, t0):
    yp = y + 1
    return yp

t1 = np.linspace(0,10, 101)
y1_0 = 1

y1 = sp.integrate.odeint(fun1, y1_0, t1)

if plots :
    plt.figure(1)
    plt.plot(t1, y1)
    plt.title('Problem 4.1')
    plt.legend('y1') 

# 4.2
def fun2(y, t0) :
    yp = y[1]
    ypp = y[2]
    yppp = y[0] - y[1]
    
    return [yp, ypp, yppp]
    
t2 = t1
y2_0 = np.ones(3)

y2 = sp.integrate.odeint(fun2, y2_0, t2)

if plots :
    plt.figure(2)
    plt.plot(t2, y2.T[0])
    plt.title('Problem 4.2')
    plt.legend(["y"])

# 4.3
def fun3(y, t0):
    y, z = y[0], y[1]
    yp = 1000*y + 1
    zp = 0.0001*z + 1000*y
    return [yp, zp]

t3 = t1
y3_0 = np.ones(2)

y3 = sp.integrate.odeint(fun3, y3_0, t3)

if plots :
    plt.figure(3)
    plt.plot(t2, y2.T[0], t2, y2.T[1])
    plt.title('Problem 4.3')
    plt.legend(["y", 'z'])

# 4.4

def forward_euler(p, q, y0=0.0, N=25, t_max=10.0) :
    """ Perform forward Euler on equation of the form
        y' = p*y**2 + q
    
    Inputs:
        p : coefficient of y(x)
        q : right hand side of differential equation
        y0 : initial value (y0 = y(0))
        N : number of points
        t_max : stopping point for approximation
        
    Returns :
        y : np.array of the approximate y values coressponding to t values
        t : np.array of the t values corresponding to y values
    """   
    t = np.linspace(0, t_max, N)
    delta = t[1]-t[0]
    y = np.zeros(N)
    y[0] = y0
    p = p*np.ones(N)
    q = q*np.ones(N)
    for i in range(1, N) :
        y[i] = y[i-1] + delta*(p[i-1]*y[i-1]**2 + q[i-1])
    return y, t

y4, t4 = forward_euler(1, 1, 1, 11, t_max=10)

if plots :
    plt.figure(4)
    plt.plot(t4, y4)
    plt.title('Problem 4.4')
    plt.yscale('log')

# 4.5

# Shooting method is used. The ode is solved as an ivp with a guess of an
# initial derivative. Newtons method is then used to find the initial
# dervivative that results in an equation that satisfies the second boundary
# condition.

t5 = t1

def solve_second_order_ivp(u0, t, p, q, f) :
    # derivative function
    def u_prime(u, t, p, q, f) :
        y, z = u[0], u[1]
        y_prime = z
        z_prime = f(t)-p(t)*z-q(t)*y
        return sp.array([y_prime, z_prime])
    u = sp.integrate.odeint(u_prime, u0, t, args=(p, q, f))
    return u.T[0] # return just y(x)
 
def r(d, a, b, y_a, y_b, p, q, f) :
    y = solve_second_order_ivp(sp.array([y_a, d]), [a, b], p, q, f)
    return y[-1] - y_b

p, q, f = lambda x: 0, lambda x: 1., lambda x: -1.
d_sol = sp.optimize.newton(r, x0=-1, fprime=None, args=(0.0, 10.0, 1., 1., p, q, f))

y5 = solve_second_order_ivp([1., d_sol], t5, p, q, f)

if plots :
    plt.figure(5)
    plt.plot(t5,y5)
    plt.title('Problem 4.5')
    