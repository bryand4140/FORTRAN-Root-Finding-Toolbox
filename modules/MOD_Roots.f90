module MOD_Roots
    use MOD_Select_Kind, only: pv
    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
    implicit none

    private
    public :: fun_zero
    public :: sv_function

    abstract interface
        function sv_function(x) result(y)
            ! Abstract type for a generic single variable function
            import pv
            real(pv), intent(in) :: x
            real(pv)             :: y
        end function sv_function
    end interface



    type, public :: fun_zero_options
        ! Optional inputs for fun_zero.
        real(pv) :: tolerance
        integer  :: maxiter

        logical :: show_iterations !Show the iteration number and the current root estimate.
        logical :: show_results    !Show the root and the status when fun_zero finishes.

    end type fun_zero_options

contains
  
subroutine fun_zero(f_pointer, bounds, status, root, niter, f0_options)
    ! General Description: Finds a zero of f(x) in the interval [a, b] such that f(a)*f(b) < 0.
    ! Inputs:
    !   f_pointer  - A pointer to the function f(x) to be evaluated.
    !   bounds     - A 2-element array containing the lower and upper bounds of the interval [a, b].
    ! Optional Inputs:
    !   f0_options - A derived type containing additional options:
    !                .tolerance, .maxiter, .show_iterations, .show_results.
    ! Outputs:
    !   status     - An integer indicating the status of the root finding process.
    !                0 - The root was found.
    !                1 - The interval does not bracket a root.
    !                2 - The maximum number of iterations was reached.
    !   root       - The computed (approximate) root of the function.
    ! Optional Outputs:
    !   niter      - The number of iterations required to find the root.

    implicit none
    procedure(sv_function), pointer :: f_pointer   
    real(pv), intent(in) :: bounds(2)
    integer, intent(out) :: status 
    real(pv), intent(out) :: root
    integer, intent(out), optional :: niter
    type(fun_zero_options), intent(inout), optional :: f0_options

    ! Local variables.
    real(pv) :: a_local, b_local, c_local, d, e
    real(pv) :: fa, fb, fc, s, m, tol_act, tol_local
    integer :: iter, maxit_local
    real(pv), parameter :: eps = 2.2204460492503131e-16
    ! Temporary variables for swapping.
    real(pv) :: temp, temp_val

    real(pv) :: current_f_value
    type(fun_zero_options) :: default_options

    logical :: show_iter_local, show_res_local

    !-------------------------------------------------------
    ! Set tolerance and maximum iterations using the options type.
    ! If f0_options is not provided, use defaults.
    !
    call init_fun_zero_options(default_options)

    if (present(f0_options)) then

        !Check to see if the user has provided a value for the tolerance.
        if(f0_options%tolerance <= 0.0_pv) then
            print*, "Warning: Tolerance is not assigned within f0_options."
            print*, "Using default tolerance value."
            tol_local = default_options%tolerance
        else
            tol_local = f0_options%tolerance
        end if

        !Check to see if the user has provided a value for the maximum 
        !number of iterations.
        if(f0_options%maxiter <= 0) then
            print*, "Warning: Maximum number of iterations is not assigned within f0_options."
            print*, "Using default maximum number of iterations."
            maxit_local = default_options%maxiter
        else
            maxit_local = f0_options%maxiter
        end if

        !Check to see if the user has provided a value for the show_iterations flag.
        if (.not. f0_options%show_iterations) then
            !User has not provided a value for show_iterations -- Use the default value.
            show_iter_local = default_options%show_iterations
        else
            show_iter_local = f0_options%show_iterations
        end if

        !Check to see if the user has provided a value for the show_results flag.
        if (.not. f0_options%show_results) then
            !User has not provided a value for show_results -- Use the default value.
            show_res_local = default_options%show_results
        else
            show_res_local = f0_options%show_results
        end if

        !Check to see if the user has entered a low value for the maximum 
        !number of iterations.
        if(f0_options%maxiter < 10) then
            print*, "Warning: Maximum number of iterations is set to a low value."
            print*, "This may affect the accuracy of the root finding process."
        end if

    else !If f0_options is not provided, use default values.
        tol_local       = default_options%tolerance
        maxit_local     = default_options%maxiter
        show_iter_local = default_options%show_iterations
        show_res_local  = default_options%show_results
    end if

    ! Initialize the bracketing interval.
    a_local = bounds(1)
    b_local = bounds(2)
    fa = f_pointer(a_local)
    fb = f_pointer(b_local)

    !Check to see if fa or fb is NaN.
    if (IEEE_IS_NAN(fa)) then
        write(*,"(A, ES12.4)") "Error: f(x) is undefined at the bound x =", a_local
        print*, " " !Spacer
        print*, " " !Spacer
        !Stop the program.
        stop
    end if

    if (IEEE_IS_NAN(fb)) then
        write(*,"(A, ES12.4)") "Error: f(x) is undefined at the bound x =", b_local
        print*, " " !Spacer
        print*, " " !Spacer
        !Stop the program.
        stop
    end if

    ! Check that the initial interval brackets a root.
    if (fa * fb >= 0.0d0) then
        print*, "Error: The interval [", a_local, ", ", b_local, "] does not bracket a root."
        print*, "Status: 1"
        status = 1
        root = a_local
        if (present(niter)) then
            niter = 0
        end if
        return
    end if

    c_local = a_local
    fc = fa
    d = b_local - a_local
    e = d

    ! Main iteration loop.
    do iter = 1, maxit_local
        if (fb * fc > 0.0d0) then
            ! Rename a_local and c_local so that b_local is the best approximation.
            c_local = a_local
            fc      = fa
            d       = b_local - a_local
            e       = d
        end if

        if (abs(fc) < abs(fb)) then
            ! Swap b_local and c_local (and corresponding function values)
            temp      = b_local
            b_local   = c_local
            c_local   = temp

            temp_val  = fb
            fb        = fc
            fc        = temp_val
        end if

        tol_act = 2.0d0 * eps * abs(b_local) + tol_local / 2.0d0
        m = (c_local - b_local) / 2.0d0

        ! Optionally show iteration details.
        if (show_iter_local) then
            current_f_value = f_pointer(b_local)

            !Check to see if current_f_value is NaN.
            if (IEEE_IS_NAN(current_f_value)) then
                write(*,"(A, ES12.4)") "Error: f(x) is NaN at x =", b_local
                print*, " " !Spacer
                print*, " " !Spacer
                !Stop the program.
                stop
            else
                write(*,"(A, I5, A, ES12.4)") "Iteration:", iter, "     f(x) =", current_f_value
            end if
        end if


        ! Check for convergence.
        if (abs(m) <= tol_act .or. fb == 0.0d0) then
            root = b_local
            if (present(niter)) then
                niter = iter
            end if
            status = 0
            if (show_res_local) then
                print*, " " !Spacer
                write(*,"(A, ES12.5, A, I10)") "Root found: ", root, " after ", iter, " iterations."
                print*,"Status: ", status
                print*, " " 
            end if
            return
        end if

        ! Decide whether to use interpolation or bisection.
        if (abs(e) >= tol_act .and. abs(fa) > abs(fb)) then
            ! Attempt inverse quadratic interpolation (or secant if only two points available).
            s = b_local - fb * (b_local - a_local) / (fb - fa)
            ! Check conditions for interpolation to be acceptable.
            if ((s < (3.0d0 * a_local + b_local) / 4.0d0) .or. (s > b_local) .or. &
                (abs(s - b_local) >= abs(m) / 2.0d0)) then
                ! If not acceptable, use bisection.
                s = b_local + m
                e = d
                d = m
            else
                d = s - b_local
            end if
        else
            ! Use bisection.
            s = b_local + m
            e = d
            d = m
        end if

        ! Prepare for next iteration.
        a_local = b_local
        fa = fb
        if (abs(d) > tol_act) then
            b_local = b_local + d
        else
            b_local = b_local + sign(tol_act, m)
        end if
        fb = f_pointer(b_local)
    end do

    ! Maximum iterations reached.
    root = b_local
    if (present(niter)) then
        niter = iter
    end if
    status = 2
    if (present(f0_options)) then
        if (f0_options%show_results) then
            write(*,*) "Maximum iterations reached. Approximate root: ", root
        end if
    end if
end subroutine fun_zero


subroutine init_fun_zero_options(f0_options)
    ! General Description: Initializes the fun_zero_options type.
    ! Outputs:
    !   f0_options - The fun_zero_options type to be initialized.

    type(fun_zero_options), intent(out) :: f0_options
    f0_options%tolerance       = 1.0e-6
    f0_options%maxiter         = 500
    f0_options%show_iterations = .false.
    f0_options%show_results    = .true.
end subroutine init_fun_zero_options


  
end module MOD_Roots
