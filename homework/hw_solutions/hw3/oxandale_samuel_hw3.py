#!/usr/bin/env python3

import numpy as np
import decimal as dec

# =============================================================================
# me701, hw3, Sam Oxandale
# =============================================================================


# Problem 1

a = 'hello'
b = 'world'
c = a + ' ' + b
d = c.title()
e, f = d.split(' ')
g, h, i = 123, 3.141592653589793, 6.022e23
j = '|'.join([str(g), str(round(h,5)), ' ' + str(round(i,3))])
k = '..'.join(str(x) for x in range(5))


# Problem 2

"""
1. This occurs because both variables point to the same memory on the computer.
    They are different names for the same spot, so when you change either 
    variable, the value at that memory location changes, changing both 
    variables.
    
2. A second list could be created by the following code :
    
    >>>>> b = a.copy()
    
    ---or---
    
    >>>>> b = list(a)
"""


# Problem 3

powers = [2**i for i in range(20)]
xpoints = [1,2,-1]
ypoints = [8,4,3,0]
zpoints = [0,-1]

points = [(x, y, z) for x in xpoints for y in ypoints for z in zpoints]

# Problem 4

def decimal_to_binary(x, n) :
    """Converts x to binary to n bits after the decimal
    
    Arguments :
        x: the number to be converted
        n: number of bits to be used to the right of the decimal
        
    Returns :
        bin_num: the binary equivelent of x to n-bits 
    """
    # Checks that values are of proper type
    assert type(x) == float, 'x must be a float'
    assert type(n) == int, 'n must be an int'
    
    # Splits x into the whole number and decimal portion
    x_split = (x//1, x % 1)
    x_int = x_split[0]
    x_dec = x_split[1]
    
    # Creates list of decimal portion of binary number
    dec_list = []
    while n > 0 :
        x_dec = x_dec*2
        # 0 if x_dec < 1
        dec_list.append(str(int(x_dec)))
        x_dec = x_dec % 1
        n = n - 1
    # Creates list of integer portion of binary number
    int_list = []
    while x_int > 0:
        int_remainder = x_int % 2
        int_list.append(str(int(int_remainder)))
        x_int = (x_int - int_remainder)/2
    int_list.reverse()
    bin_dec = ''.join(dec_list)
    bin_int = ''.join(int_list)
    bin_num = str(bin_int) + '.' + str(bin_dec)
    return bin_num
    
def binary_to_decimal(i) :
    """Converts x to binary to n bits after the decimal
    
    Arguments :
        i: the binary number to be converted to decimal
        
    Returns :
        dec: the decimal equivelent of i
    """
    # Splits x into the whole number and decimal portion
    i_split = i.split('.')
    i_int_list = list(i_split[0])
    i_dec_list = list(i_split[1])
    
    # Initializes values
    base_int = 0
    base_dec = 0
    
    # Calculates the whole number portion
    n = 0
    while n < len(i_int_list):
        power = len(i_int_list) - 1 - n
        base_int = base_int + int(i_int_list[n])*(2**power)
        n += 1
    
    # Calculates the decimal portion
    n = 0
    while n < len(i_dec_list):
        base_dec = base_dec + int(i_dec_list[n])*(2**(-1-n))
        n += 1
        
    dec = float(str(int(base_int)) + str(base_dec)[1:])
    return dec


# Problem 5
    
# option 1
def opt1(a) :
    s=0
    for i in range(len(a)) :
        s += a[i]
    return s

# option 2 s=0
def opt2(a) :
    a = sorted(a)
    s = 0
    for i in range(len(a)) :
        s += a[i]
    return s

# option 3
def sumr(a):
    if len(a) <= 2:
        return sum(a)
    else:
        return sumr(a[:len(a)//2]) + sumr(a[len(a)//2:])

n = 100000
array = np.random.rand(n) 

print(opt1(array)==sum(array))

print('   opt1 = %.50f' % opt1(array))
print('   opt2 = %.50f' % opt2(array))
print('   sumr = %.50f' % sumr(array))
print('    sum = %.50f' % sum(array))
print(' np.sum = %.50f' % np.sum(array))
print('decimal = %.50f' % sum(map(dec.Decimal, array)))

print('\nopt1 and sum are equivelent and are accurate to about 8 or 9 places.')
print('opt2 has the same accuracy but returns a slightly different value.')
print('sumr, np.sum, and decimal return the same value and are most accurate.')

