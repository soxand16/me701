#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

double f(double *x, int n)
{
    /*
    Calculates f(x1, x2, ... xn) where f = (x1 + x2 + ... + xn)**2
      
     Arguments:
        *x - array of x-values to evaluate function at
        n - number of x-values
     
     Returns:
        value - f(x1, x2, ... xn)
    */
    
    double value;
    value = 0.0;
    
    // Sums all of the x values
    for ( int i = 0; i < n; i = i + 1)
    {
        value = value + x[i];
    }
    
    value = value*value;
    
    return value;
}

double factorial(double n)
{
    if (n<0)
    {
        return 0;
    }
    if (n==0)
    {
        return 1;
    }
    double factorial = 1.0;
    for ( int i = 2; i <= n; i = i + 1)
    {
        factorial *= i;
    }
    return factorial;
}

double monte_carlo(int n, int m)
{
    double x[n];
    double f_val = 0.0;
    double f_test;
    // Function max is where all x=1, therefore it is number of elements squared
    // This is also the area of evaluation
    double f_max = n*n;
    double total = 0.0;
    
    for ( int j = 0; j < pow(m, n); j = j + 1)
    {
        for ( int i = 0; i < n; i = i + 1)
        {
            // Gets random x-values
            x[i] = ((double)rand()/RAND_MAX);
            
        }
        // Evaluates function for random x-values
        f_val = f_val + f(x, n);
        total = total + 1;
    }
    return f_val / total;
}

double midpoint(int n, int m)
{
    double x_vals[m];
    double x[n];
    int index[n];
    double I = 0.0;
    
    // Gets values for different x's
    for ( int i = 0; i < m; i = i + 1)
    {
        x_vals[i] = (i+0.5)/m;
    }
    
    // Initializes indexes
    for ( int i = 0; i < n; i = i + 1)
    {
        index[i] = 0;
    }
    
    // Counts up the indexes to get all of the combos
    while (1)
    {
        index[n-1] = index[n-1] + 1;
        for ( int i = n - 1; i > 0; i = i - 1)
        {
            if (index[i] == m)
            {
                index[i] = 0;
                index[i-1] = index[i-1] + 1;
            }
        }
        
        if (index[0] == m)
        {
            break;
        }
        
        for ( int i = 0; i < n; i = i + 1)
        {
            x[i] = x_vals[index[i]];
        }
        
        I = I + f(x,n)/m;
        

    }
    
    return I;
}

double I_ref(int n, int m)
{
    double bi_coeff;
    if (n < 2)
    {
        bi_coeff = 0;
    }
    else
    {
    bi_coeff = factorial(n)/(factorial(n-2)*2.0);
    }
    return bi_coeff/2 + n/3.0;
}

int main( int argc, char *argv[] )
{
    int n, m;
    if ( argc != 3 )
    {
        cout<<"usage: "<< argv[0] <<" <n> <m>\n"<< endl;
    }
    else
    {
        n = atoi(argv[1]);
        m = atoi(argv[2]);
        cout<<"monte_carlo = "<<monte_carlo(n,m)<<endl;
        cout<<"midpoint = "<<midpoint(n,m)<<endl;
        cout<<"I_ref = "<<I_ref(n,m)<<endl;
    }
    return 0;
}
