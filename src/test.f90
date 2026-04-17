program test
    use MOD_Select_Kind, only: pv
    use MOD_Roots
    implicit none

    real(pv) :: bounds(2), root
    integer :: status, niter
    procedure(sv_function), pointer :: f_pointer => null()
    procedure(sv_function), pointer :: df_pointer => null()
    type(fun_zero_options) :: options

    ! -------------------------------------------------------
    ! Test 1: fun_zero (Brent / bisection) on a non-trivial function.
    ! -------------------------------------------------------
    print*, "============================================"
    print*, "  Test 1: fun_zero"

    bounds = [-2.0_pv, 2.0_pv]
    f_pointer => my_function

    options%show_iterations = .true.
    options%show_results    = .true.
    options%maxiter         = 100
    options%tolerance       = 1.0E-10_pv

    call fun_zero(f_pointer, bounds, status, root, niter = niter, f0_options = options)

    if (status == 0) then
        write(*,'(A, ES15.8, A, I0)') &
            "  Verified root: ", root, "  iterations: ", niter
        write(*,'(A, ES10.3)') &
            "  Residual |f(root)|: ", abs(my_function(root))
    else
        print*, "  fun_zero returned status:", status
    end if

    ! -------------------------------------------------------
    ! Test 2: newton_raphson on the same function.
    ! -------------------------------------------------------
    print*, " "
    print*, "============================================"
    print*, "  Test 2: newton_raphson"

    f_pointer  => my_function
    df_pointer => my_function_deriv

    options%show_iterations = .true.
    options%show_results    = .true.
    options%maxiter         = 100
    options%tolerance       = 1.0E-10_pv

    call newton_raphson(f_pointer, df_pointer, x0 = -1.0_pv, &
                        status = status, root = root, niter = niter, &
                        nr_options = options)

    if (status == 0) then
        write(*,'(A, ES15.8, A, I0)') &
            "  Verified root: ", root, "  iterations: ", niter
        write(*,'(A, ES10.3)') &
            "  Residual |f(root)|: ", abs(my_function(root))
    else
        print*, "  newton_raphson returned status:", status
    end if

    print*, "============================================"

contains

    function my_function(x) result(y)
        real(pv), intent(in) :: x
        real(pv) :: y
        real(pv), parameter :: pi = acos(-1.0_pv)   ! Full machine-precision pi.
        y = -exp(-x) + 2.0_pv * x**2 / exp(4.0_pv * x) - pi
    end function my_function

    function my_function_deriv(x) result(dy)
        ! Analytic derivative of my_function for use with newton_raphson.
        real(pv), intent(in) :: x
        real(pv) :: dy
        ! d/dx [ -exp(-x) ] = exp(-x)
        ! d/dx [ 2*x^2 * exp(-4x) ] = (4x - 8x^2) * exp(-4x)
        dy = exp(-x) + (4.0_pv * x - 8.0_pv * x**2) * exp(-4.0_pv * x)
    end function my_function_deriv

end program test