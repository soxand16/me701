#include <stdio.h>
#include <iostream>
using namespace std;

double interpolate(double x_new, double *x, double *y, int n, int order)
{
    /*
    Interpolates between points contained in x and y at the location of x_new using
    n points and using interpolation of order order
     
     Arguments:
        x_new - new x-point
        *x - array of known x-points
        *y - array of known y-points
        n - number of points passed in *x and *y
        order - order of interpolation to be used
     
     Returns:
        y_new - new y-point
    */
    double y_new;
    double coeff[order+1];
    double x_test;
    int index;
    int start;
    
    // Check that enough points are passed to perform an interpolation
    // of the given order
    if (n <= order)
    {
        cout << "Number of points must be greater than order" << endl;
        return 0.0;
    }
    
    index = -1;
    
    // Find the index of the first x value that is greater than x_new
    do
    {
        index += 1;
        x_test = x[index];
        
    }while(x_new>x_test and index<n);
    
    // Find the starting point for interpolation values
    start = index - order/2 - order%2;
    
    // If the starting point is less than 0, make it 0
    if (start<0)
    {
        start = 0;
    }
    // If end point is out of the range of points, set
    // start so the last point used will be the last point
    if (start+order+1>n-1)
    {
        start = n - 1 - order;
    }
    
    y_new = 0;
    
    // Fill up the coefficient array
    for ( int i = 0; i < order + 1; i = i + 1)
    {
        coeff[i] = y[i+start];
        for ( int j = 0; j < order + 1; j = j + 1)
        {
            if (i != j)
            {
                coeff[i] = coeff[i] * (x_new - x[j+start]) / (x[i+start] - x[j+start]);
            }
        }
        // Add that term to y_new
        y_new = y_new + coeff[i];
    }
    
    return y_new;
}


int main()
{
    double x_new = 1;
    int n = 4 ;
    double x[4] = {0,2,3,4};
    double y[4] = {0,4,9,10};
    int order = 1;
    double y_new ;
    
    y_new = interpolate(x_new, x, y, n, order);
    
    cout << y_new << endl;
}

