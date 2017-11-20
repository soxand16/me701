program main

    double precision :: x_new
    integer :: n=4
    double precision, dimension(1:4) :: x, y
    integer :: order
    double precision :: y_new

    x_new = 1
    n = 4
    x = (/ 0,2,3,4 /)
    y = (/ 0,4,9,10 /)
    order = 2

    y_new = interpolate(x_new, x, y, n, order)
    print *, y_new


CONTAINS
    double precision function interpolate(x_new, x, y, n, order)
        ! Interpolates between points contained in x and y at the location of x_new using
        ! n points and using interpolation of order order
        !
        ! Arguments:
        !   x_new - new x-point
        !   x - array of known x-points
        !   y - array of known y-points
        !   n - number of points passed in *x and *y
        !   order - order of interpolation to be used
        !
        ! Returns:
        !   y_new - new y-point
        double precision, intent(in) :: x_new
        double precision, dimension(:), intent(in) :: x, y
        integer, intent(in) :: n, order

        double precision :: y_new
        double precision, dimension(n) :: coeff
        double precision :: x_test
        integer :: index
        integer :: start

        ! Check that enough points are passed to perform an interpolation
        ! of the given order
        if (n <= order) then
            print *, "Number of points must be greater than order"
            stop
        end if

        ! Find the index of the first x value that is greater than x_new
        index = 0
        x_test = x(index)
        do while (x_new>x_test .AND. index<n)
            index = index + 1
            x_test = x(index)
        end do

        ! Find the starting point for interpolation values
        start = index - order/2 - mod(order, 2)

        ! If the starting point is less than 0, make it 0
        if (start<0) then
            start = 0
        end if
        ! If end point is out of the range of points, set
        ! start so the last point used will be the last point
        if (start+order+1>n-1) then
            start = n - 1 - order
        end if

        y_new = 0

        ! Fill up the coefficient array
        DO i = 1, order + 1, 1
            coeff(i) = y(i+start)
            DO j = 1, order + 1, 1
                if (i /= j) then
                    coeff(i) = coeff(i) * (x_new - x(j+start)) / (x(i+start) - x(j+start))
                end if
            END DO
            ! Add that term to y_new
            y_new = y_new + coeff(i)
        END DO
        interpolate =  y_new

    end function interpolate
end program main
