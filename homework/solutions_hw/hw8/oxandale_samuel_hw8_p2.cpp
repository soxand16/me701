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
        cout<<monte_carlo(n,m)<<endl;
    }
    return 0;
}
