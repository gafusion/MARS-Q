program polynomial_integral_oracle
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer, parameter :: n = 7
    real(dp) :: x(n), y(n), m(n), actual(n), expected(n), coeff(0:3)
    real(dp) :: at(10), cipr(n), cidq(n), exact_m(n), error, anchor
    integer :: degree, direction, k, i
    character(len=32) :: arg
    external :: native_integrate, native_spline
    call get_command_argument(1, arg)
    read(arg, *) degree
    call get_command_argument(2, arg)
    read(arg, *) direction
    x = [-0.70_dp, -0.53_dp, -0.31_dp, -0.20_dp, -0.07_dp, -0.01_dp, 0.0_dp]
    coeff = 0.0_dp
    coeff(0) = 0.2_dp
    do k = 1, degree
        coeff(k) = 0.2_dp + 0.1_dp*k
    end do
    y = 0.0_dp
    exact_m = 0.0_dp
    do k = 0, degree
        y = y + coeff(k)*x**k
        if (k >= 2) exact_m = exact_m + k*(k - 1)*coeff(k)*x**(k - 2)
    end do
    call native_spline(n, x, y, m)
    if (.not. all(abs(m - exact_m) <= 2.0e-11_dp)) error stop 'Native M is not d2/dpsi2'
    at = 0.0_dp
    cipr = 1.0_dp
    cidq = 1.0_dp
    actual = -999.0_dp
    call native_integrate(n, direction, 1, 0, at, x, y, m, actual, cipr, cidq)
    anchor = x(n)
    if (direction == 1) anchor = x(1)
    expected = 0.5_dp
    do i = 1, n
        do k = 0, degree
            expected(i) = expected(i) + coeff(k)/(k + 1)*(x(i)**(k + 1) - &
                anchor**(k + 1))
        end do
    end do
    error = maxval(abs(actual - expected))
    print '(a,i0,a,i0,a,es17.9)', 'degree=', degree, ' direction=', direction, &
        ' independent primitive max error=', error
    if (.not. all(abs(actual - expected) <= 3.0e-13_dp)) error stop 'Native cubic spline integral disagrees'
end program
