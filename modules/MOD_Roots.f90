module MOD_Roots
    use MOD_Select_Kind, only: pv
    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
    implicit none

    private
    public :: fun_zero
    public :: sv_function
    public :: poly_roots
    

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


SUBROUTINE QR_Decomposition(A, Q, R, n, status)
    !General Description: This subroutine performs a QR decomposition of a matrix A
    !                     using the Gram-Schmidt process. The QR decomposition is
    !                     a factorization of a matrix A into a product A = QR of an
    !                     orthogonal matrix Q and an upper triangular matrix R.
    !                     The QR decomposition is often used to solve the linear
    !                     least squares problem, and is the basis for a particular
    !                     eigenvalue algorithm.

    !Input Parameters:
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n ! Size of the matrix A
    REAL(pv), DIMENSION(n,n), INTENT(IN) :: A ! Matrix to be decomposed
    REAL(pv), DIMENSION(n,n), INTENT(OUT) :: Q, R ! Q and R matrices

    !Local Variables:
    REAL(pv), DIMENSION(n) :: u, e ! Temporary vectors
    INTEGER :: i, j ! Loop indices
    INTEGER, INTENT(OUT) :: status ! Status of the decomposition

    ! Initialize status to success
    status = 0

    ! Perform Error Checking
    IF (SIZE(A,1) /= SIZE(A,2)) THEN
        PRINT*, "ERROR: Matrix A must be square"
        status = 1 ! Indicate error
        RETURN
    END IF

    DO i = 1, n
        u = A(:,i)
        DO j = 1, i-1
            R(j,i) = DOT_PRODUCT(Q(:,j), A(:,i))
            u = u - R(j,i) * Q(:,j)
        END DO
        IF (NORM2(u) == 0.0_pv) THEN
            PRINT*, "Error: Matrix is singular or nearly singular"
            status = 2 ! Indicate error
            RETURN
        END IF
        R(i,i) = NORM2(u)
        e = u / R(i,i)
        Q(:,i) = e
    END DO
END SUBROUTINE QR_Decomposition


subroutine poly_roots(coeff, deg, roots, status)
    ! Inputs:
    !   coeff  - Real array of polynomial coefficients
    !            (e.g., coeff(0) = constant term, coeff(deg) = leading coefficient)
    !   deg    - Degree of the polynomial (deg >= 1)
    ! Outputs:
    !   roots  - Complex array of roots (length = deg)
    !   status - Status flag: 0 if successful, nonzero otherwise.

    implicit none
    integer, intent(in)      :: deg
    real(pv), intent(in)     :: coeff(0:deg)
    complex(pv), intent(out) :: roots(deg)
    integer, intent(out)     :: status

    ! Local variables.
    integer :: i, j, maxiter, iter, n
    real(pv) :: tol
    real(pv), allocatable :: comp_mat(:,:)
    real(pv), allocatable :: Q(:,:), R(:,:)
    integer :: qr_status
    real(pv) :: offdiag_norm

    ! Parameters for convergence.
    tol = 1.0e-8_pv
    maxiter = 1000
    n = deg  ! The companion matrix is of size n x n

    ! Allocate the companion matrix and Q, R matrices.
    allocate(comp_mat(n, n))
    allocate(Q(n, n))
    allocate(R(n, n))

    ! Form the companion matrix.
    ! First, normalize coefficients (make polynomial monic)
    if (abs(coeff(deg)) < tol) then
       status = 1
       print*, "Error: Leading coefficient is zero."
       return
    end if

    ! Set the companion matrix to zero.
    comp_mat = 0.0_pv

    ! Fill in the subdiagonal with ones.
    do i = 1, n-1
       comp_mat(i+1, i) = 1.0_pv
    end do

    ! Fill in the first row: negative of normalized coefficients.
    ! Note: if we assume coeff(0) is the constant term and coeff(deg) is the leading coefficient,
    ! then the first row becomes: -coeff(deg-1:0)/coeff(deg)
    do j = 1, n
       comp_mat(1, j) = -coeff(j-1) / coeff(deg)
    end do

    ! QR Iteration loop.
    iter = 0
    do while (iter < maxiter)
       iter = iter + 1

       ! Call your QR decomposition subroutine.
       call QR_Decomposition(comp_mat, Q, R, n, qr_status)
       if (qr_status /= 0) then
          status = 2
          print*, "Error in QR decomposition."
          return
       end if

       ! Update comp_mat = R * Q.
       comp_mat = matmul(R, Q)

       ! Check for convergence.
       ! For a nearly upper triangular matrix, the subdiagonal entries should be small.
       offdiag_norm = 0.0_pv
       do i = 2, n
          offdiag_norm = offdiag_norm + abs(comp_mat(i, i-1))**2
       end do
       offdiag_norm = sqrt(offdiag_norm)
       if (offdiag_norm < tol) exit
    end do

    if (iter >= maxiter) then
       status = 3
       print*, "QR iteration did not converge."
    else
       status = 0
    end if

    ! Extract the eigenvalues from the (nearly) upper triangular matrix.
    do i = 1, n
       ! This simple version ignores the possibility of 2x2 blocks for complex pairs.
       roots(i) = cmplx(comp_mat(i,i), 0.0_pv, kind=pv)
    end do

    ! Deallocate local arrays.
    deallocate(comp_mat, Q, R)
end subroutine poly_roots


subroutine qr_algeq_solver(n, c, zr, zi, istatus, detil)

    implicit none

    integer, intent(in)   :: n !! degree of the monic polynomial.
    real(pv), intent(in)  :: c(n + 1) !! coefficients of the polynomial. in order of decreasing powers.
    real(pv), intent(out) :: zr(n) !! real part of output roots
    real(pv), intent(out) :: zi(n) !! imaginary part of output roots
    integer, intent(out)  :: istatus !! return code:
                                    !!
                                    !! * -1 : degree <= 0
                                    !! * -2 : leading coefficient `c(1)` is zero
                                    !! * 0 : success
                                    !! * otherwise, the return code from `hqr_eigen_hessenberg`
    real(pv), intent(out), optional :: detil !! accuracy hint.

    real(pv), allocatable :: a(:, :) !! work matrix
    integer, allocatable :: cnt(:) !! work area for counting the qr-iterations
    integer :: i, iter
    real(pv) :: afnorm

    ! check for invalid arguments
    if (n <= 0) then
        istatus = -1
        return
    end if
    if (c(1) == 0.0_pv) then
        ! leading coefficient is zero.
        istatus = -2
        return
    end if

    allocate (a(n, n))
    allocate (cnt(n))

    ! build the companion matrix a.
    call build_companion(n, a, c)

    ! balancing the a itself.
    call balance_companion(n, a)

    ! qr iterations from a.
    call hqr_eigen_hessenberg(n, a, zr, zi, cnt, istatus)
    if (istatus /= 0) then
        write (*, '(A,1X,I4)') 'abnormal return from hqr_eigen_hessenberg. istatus=', istatus
        if (istatus == 1) write (*, '(A)') 'matrix is completely zero.'
        if (istatus == 2) write (*, '(A)') 'qr iteration did not converge.'
        if (istatus > 3) write (*, '(A)') 'arguments violate the condition.'
        return
    end if

    if (present(detil)) then

        ! compute the frobenius norm of the balanced companion matrix a.
        afnorm = frobenius_norm_companion(n, a)

        ! count the total qr iteration.
        iter = 0
        do i = 1, n
            if (cnt(i) > 0) iter = iter + cnt(i)
        end do

        ! calculate the accuracy hint.
        detil = eps*n*iter*afnorm

    end if
end subroutine qr_algeq_solver


subroutine build_companion(n, a, c)

    !!  build the companion matrix of the polynomial.
    !!  (this was modified to allow for non-monic polynomials)

    implicit none

    integer, intent(in) :: n
    real(pv), intent(out) :: a(n, n)
    real(pv), intent(in) :: c(n + 1) !! coefficients in order of decreasing powers

    integer :: i !! counter

    ! create the companion matrix
    a = 0.0_pv
    do i = 1, n - 1
        a(i + 1, i) = 1.0_pv
    end do
    do i = n, 1, -1
        a(n - i + 1, n) = -c(i + 1)/c(1)
    end do

end subroutine build_companion

subroutine balance_companion(n, a)

    !!  blancing the unsymmetric matrix `a`.
    !!
    !!  this fortran code is based on the algol code "balance" from paper:
    !!   "balancing a matrix for calculation of eigenvalues and eigenvectors"
    !!   by b.n.parlett and c.reinsch, numer. math. 13, 293-304(1969).
    !!
    !!  note: the only non-zero elements of the companion matrix are touched.

    implicit none

    integer, intent(in) :: n
    real(pv), intent(inout) :: a(n, n)

    integer, parameter :: b = radix(1.0_pv) !! base of the floating point representation on the machine
    integer, parameter :: b2 = b**2

    integer :: i, j
    real(pv) :: c, f, g, r, s
    logical :: noconv

    if (n <= 1) return ! do nothing

    ! iteration:
    main: do
        noconv = .false.
        do i = 1, n
            ! touch only non-zero elements of companion.
            if (i /= n) then
                c = abs(a(i + 1, i))
            else
                c = 0.0_pv
                do j = 1, n - 1
                    c = c + abs(a(j, n))
                end do
            end if
            if (i == 1) then
                r = abs(a(1, n))
            elseif (i /= n) then
                r = abs(a(i, i - 1)) + abs(a(i, n))
            else
                r = abs(a(i, i - 1))
            end if

            if (c /= 0.0_pv .and. r /= 0.0_pv) then

                g = r/b
                f = 1.0_pv
                s = c + r
                do
                    if (c >= g) exit
                    f = f*b
                    c = c*b2
                end do
                g = r*b
                do
                    if (c < g) exit
                    f = f/b
                    c = c/b2
                end do
                if ((c + r)/f < 0.95_pv*s) then
                    g = 1.0_pv/f
                    noconv = .true.
                    ! touch only non-zero elements of companion.
                    if (i == 1) then
                        a(1, n) = a(1, n)*g
                    else
                        a(i, i - 1) = a(i, i - 1)*g
                        a(i, n) = a(i, n)*g
                    end if
                    if (i /= n) then
                        a(i + 1, i) = a(i + 1, i)*f
                    else
                        do j = 1, n
                            a(j, i) = a(j, i)*f
                        end do
                    end if
                end if
            end if
        end do
        if (noconv) cycle main
        exit main
    end do main

end subroutine balance_companion


function frobenius_norm_companion(n, a) result(afnorm)

    !!  calculate the frobenius norm of the companion-like matrix.
    !!  note: the only non-zero elements of the companion matrix are touched.

    implicit none

    integer, intent(in) :: n
    real(pv), intent(in) :: a(n, n)
    real(pv) :: afnorm

    integer :: i, j
    real(pv) :: sum

    sum = 0.0_pv
    do j = 1, n - 1
        sum = sum + a(j + 1, j)**2
    end do
    do i = 1, n
        sum = sum + a(i, n)**2
    end do
    afnorm = sqrt(sum)

end function frobenius_norm_companion


subroutine hqr_eigen_hessenberg(n0, h, wr, wi, cnt, istatus)

    !!  eigenvalue computation by the householder qr method
    !!  for the real hessenberg matrix.
    !!
    !! this fortran code is based on the algol code "hqr" from the paper:
    !!       "the qr algorithm for real hessenberg matrices"
    !!       by r.s.martin, g.peters and j.h.wilkinson,
    !!       numer. math. 14, 219-231(1970).
    !!
    !! comment: finds the eigenvalues of a real upper hessenberg matrix,
    !!          h, stored in the array h(1:n0,1:n0), and stores the real
    !!          parts in the array wr(1:n0) and the imaginary parts in the
    !!          array wi(1:n0).
    !!          the procedure fails if any eigenvalue takes more than
    !!          `maxiter` iterations.

    implicit none

    integer, intent(in) :: n0
    real(pv), intent(inout) :: h(n0, n0)
    real(pv), intent(out) :: wr(n0)
    real(pv), intent(out) :: wi(n0)
    integer, intent(inout) :: cnt(n0)
    integer, intent(out) :: istatus

    integer :: i, j, k, l, m, na, its, n
    real(pv) :: p, q, r, s, t, w, x, y, z
    logical :: notlast, found

    maxiter = 500 ! maximum number of iterations for the qr algorithm

    ! note: n is changing in this subroutine.
    n = n0
    istatus = 0
    t = 0.0_pv

    main: do

        if (n == 0) return

        its = 0
        na = n - 1

        do

            ! look for single small sub-diagonal element
            found = .false.
            do l = n, 2, -1
                if (abs(h(l, l - 1)) <= eps*(abs(h(l - 1, l - 1)) + abs(h(l, l)))) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) l = 1

            x = h(n, n)
            if (l == n) then
                ! one root found
                wr(n) = x + t
                wi(n) = 0.0_pv
                cnt(n) = its
                n = na
                cycle main
            else
                y = h(na, na)
                w = h(n, na)*h(na, n)
                if (l == na) then
                    ! comment: two roots found
                    p = (y - x)/2
                    q = p**2 + w
                    y = sqrt(abs(q))
                    cnt(n) = -its
                    cnt(na) = its
                    x = x + t
                    if (q > 0.0_pv) then
                        ! real pair
                        if (p < 0.0_pv) y = -y
                        y = p + y
                        wr(na) = x + y
                        wr(n) = x - w/y
                        wi(na) = 0.0_pv
                        wi(n) = 0.0_pv
                    else
                        ! complex pair
                        wr(na) = x + p
                        wr(n) = x + p
                        wi(na) = y
                        wi(n) = -y
                    end if
                    n = n - 2
                    cycle main
                else
                    if (its == maxiter) then ! 30 for the original double precision code
                        istatus = 1
                        return
                    end if
                    if (its == 10 .or. its == 20) then
                        ! form exceptional shift
                        t = t + x
                        do i = 1, n
                            h(i, i) = h(i, i) - x
                        end do
                        s = abs(h(n, na)) + abs(h(na, n - 2))
                        y = 0.75_pv*s
                        x = y
                        w = -0.4375_pv*s**2
                    end if
                    its = its + 1
                    ! look for two consecutive small sub-diagonal elements
                    do m = n - 2, l, -1
                        z = h(m, m)
                        r = x - z
                        s = y - z
                        p = (r*s - w)/h(m + 1, m) + h(m, m + 1)
                        q = h(m + 1, m + 1) - z - r - s
                        r = h(m + 2, m + 1)
                        s = abs(p) + abs(q) + abs(r)
                        p = p/s
                        q = q/s
                        r = r/s
                        if (m == l) exit
                        if (abs(h(m, m - 1))*(abs(q) + abs(r)) <= eps*abs(p) &
                            *(abs(h(m - 1, m - 1)) + abs(z) + abs(h(m + 1, m + 1)))) exit
                    end do

                    do i = m + 2, n
                        h(i, i - 2) = 0.0_pv
                    end do
                    do i = m + 3, n
                        h(i, i - 3) = 0.0_pv
                    end do
                    ! double qr step involving rows l to n and columns m to n
                    do k = m, na
                        notlast = (k /= na)
                        if (k /= m) then
                            p = h(k, k - 1)
                            q = h(k + 1, k - 1)
                            if (notlast) then
                                r = h(k + 2, k - 1)
                            else
                                r = 0.0_pv
                            end if
                            x = abs(p) + abs(q) + abs(r)
                            if (x == 0.0_pv) cycle
                            p = p/x
                            q = q/x
                            r = r/x
                        end if
                        s = sqrt(p**2 + q**2 + r**2)
                        if (p < 0.0_pv) s = -s
                        if (k /= m) then
                            h(k, k - 1) = -s*x
                        elseif (l /= m) then
                            h(k, k - 1) = -h(k, k - 1)
                        end if
                        p = p + s
                        x = p/s
                        y = q/s
                        z = r/s
                        q = q/p
                        r = r/p
                        ! row modification
                        do j = k, n
                            p = h(k, j) + q*h(k + 1, j)
                            if (notlast) then
                                p = p + r*h(k + 2, j)
                                h(k + 2, j) = h(k + 2, j) - p*z
                            end if
                            h(k + 1, j) = h(k + 1, j) - p*y
                            h(k, j) = h(k, j) - p*x
                        end do
                        if (k + 3 < n) then
                            j = k + 3
                        else
                            j = n
                        end if
                        ! column modification;
                        do i = l, j
                            p = x*h(i, k) + y*h(i, k + 1)
                            if (notlast) then
                                p = p + z*h(i, k + 2)
                                h(i, k + 2) = h(i, k + 2) - p*r
                            end if
                            h(i, k + 1) = h(i, k + 1) - p*q
                            h(i, k) = h(i, k) - p
                        end do
                    end do
                    cycle
                end if
            end if

        end do

    end do main

end subroutine hqr_eigen_hessenberg
  
!-----------------------------------------------------------------------------------------------
end module MOD_Roots
