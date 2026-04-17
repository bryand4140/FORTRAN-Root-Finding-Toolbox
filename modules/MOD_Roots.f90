module MOD_Roots
    use MOD_Select_Kind, only: pv
    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
    implicit none

    private
    public :: fun_zero
    public :: sv_function
    public :: poly_roots
    public :: newton_raphson
    public :: init_fun_zero_options
    

    abstract interface
        function sv_function(x) result(y)
            ! Abstract type for a generic single variable function
            import pv
            real(pv), intent(in) :: x
            real(pv)             :: y
        end function sv_function
    end interface


    type, public :: fun_zero_options
        ! Optional inputs for fun_zero / newton_raphson.
        ! Default values are provided so that an uninitialized struct behaves sensibly.
        real(pv) :: tolerance       = 1.0e-6_pv
        integer  :: maxiter         = 500

        logical :: show_iterations  = .false. !Show the iteration number and the current root estimate.
        logical :: show_results     = .true.  !Show the root and the status when the solver finishes.

    end type fun_zero_options
contains
  
subroutine fun_zero(f_pointer, bounds, status, root, niter, f0_options)
    ! General Description: Finds a zero of f(x) in the interval [a, b] such that f(a)*f(b) < 0.
    !   Uses a hybrid Brent / bisection algorithm for robust, fast convergence.
    ! Inputs:
    !   f_pointer  - A pointer to the function f(x) to be evaluated.
    !   bounds     - A 2-element array containing the lower and upper bounds of the interval [a, b].
    ! Optional Inputs:
    !   f0_options - A derived type containing additional options:
    !                .tolerance, .maxiter, .show_iterations, .show_results.
    !                All fields have sensible defaults; see fun_zero_options type definition.
    ! Outputs:
    !   status     - An integer indicating the status of the root finding process.
    !                0 - The root was found.
    !                1 - The interval does not bracket a root.
    !                2 - The maximum number of iterations was reached.
    !                3 - f(x) is NaN at one of the initial bounds.
    !   root       - The computed (approximate) root of the function.
    ! Optional Outputs:
    !   niter      - The number of iterations required to find the root.

    implicit none
    procedure(sv_function), pointer :: f_pointer
    real(pv), intent(in) :: bounds(2)
    integer, intent(out) :: status
    real(pv), intent(out) :: root
    integer, intent(out), optional :: niter
    type(fun_zero_options), intent(in), optional :: f0_options

    ! Local variables.
    real(pv) :: a_local, b_local, c_local, d, e
    real(pv) :: fa, fb, fc, s, m, tol_act, tol_local
    integer :: iter, maxit_local
    ! Use machine epsilon for the convergence tolerance term (Brent's algorithm).
    real(pv), parameter :: eps = epsilon(1.0_pv)
    ! Temporary variables for swapping.
    real(pv) :: temp, temp_val
    type(fun_zero_options) :: opts

    logical :: show_iter_local, show_res_local

    !-------------------------------------------------------
    ! Merge user-supplied options with defaults.
    ! Because fun_zero_options fields carry default values, opts is already
    ! initialised correctly via init_fun_zero_options; user overrides are applied below.
    call init_fun_zero_options(opts)
    if (present(f0_options)) then
        if (f0_options%tolerance > 0.0_pv) opts%tolerance = f0_options%tolerance
        if (f0_options%maxiter   > 0)       opts%maxiter   = f0_options%maxiter
        ! Use the caller's logical flags directly – .false. is a valid explicit choice.
        opts%show_iterations = f0_options%show_iterations
        opts%show_results    = f0_options%show_results

        if (opts%tolerance <= 0.0_pv) then
            print*, "Warning: Tolerance is not assigned within f0_options."
            print*, "Using default tolerance value."
            opts%tolerance = 1.0e-6_pv
        end if
        if (opts%maxiter < 10) then
            print*, "Warning: Maximum number of iterations is set to a low value."
            print*, "This may affect the accuracy of the root finding process."
        end if
    end if

    tol_local       = opts%tolerance
    maxit_local     = opts%maxiter
    show_iter_local = opts%show_iterations
    show_res_local  = opts%show_results

    ! Initialize the bracketing interval.
    a_local = bounds(1)
    b_local = bounds(2)
    fa = f_pointer(a_local)
    fb = f_pointer(b_local)

    ! Check for NaN at the bounds – return status 3 instead of stopping the program.
    if (IEEE_IS_NAN(fa)) then
        write(*,"(A, ES12.4)") "Error: f(x) is undefined at the bound x =", a_local
        status = 3
        root   = a_local
        if (present(niter)) niter = 0
        return
    end if

    if (IEEE_IS_NAN(fb)) then
        write(*,"(A, ES12.4)") "Error: f(x) is undefined at the bound x =", b_local
        status = 3
        root   = b_local
        if (present(niter)) niter = 0
        return
    end if

    ! Check that the initial interval brackets a root.
    if (fa * fb >= 0.0_pv) then
        print*, "Error: The interval [", a_local, ", ", b_local, "] does not bracket a root."
        status = 1
        root   = a_local
        if (present(niter)) niter = 0
        return
    end if

    c_local = a_local
    fc = fa
    d = b_local - a_local
    e = d

    ! Main iteration loop.
    do iter = 1, maxit_local
        if (fb * fc > 0.0_pv) then
            ! Rename a_local and c_local so that b_local is the best approximation.
            c_local = a_local
            fc      = fa
            d       = b_local - a_local
            e       = d
        end if

        if (abs(fc) < abs(fb)) then
            ! Swap b_local and c_local (and corresponding function values).
            temp      = b_local
            b_local   = c_local
            c_local   = temp

            temp_val  = fb
            fb        = fc
            fc        = temp_val
        end if

        tol_act = 2.0_pv * eps * abs(b_local) + tol_local / 2.0_pv
        m = (c_local - b_local) / 2.0_pv

        ! Optionally show iteration details.
        ! fb already equals f(b_local) from the previous step; no extra evaluation needed.
        if (show_iter_local) then
            if (IEEE_IS_NAN(fb)) then
                write(*,"(A, ES12.4)") "Error: f(x) is NaN at x =", b_local
                status = 3
                root   = b_local
                if (present(niter)) niter = iter
                return
            end if
            write(*,"(A, I5, A, ES12.4, A, ES12.4)") &
                "Iteration:", iter, "    x =", b_local, "    f(x) =", fb
        end if

        ! Check for convergence.
        if (abs(m) <= tol_act .or. fb == 0.0_pv) then
            root = b_local
            if (present(niter)) niter = iter
            status = 0
            if (show_res_local) then
                print*, " "
                write(*,"(A, ES12.5, A, I0, A)") &
                    "Root found: ", root, " after ", iter, " iterations."
                print*, "Status: ", status
                print*, " "
            end if
            return
        end if

        ! Decide whether to use interpolation or bisection.
        if (abs(e) >= tol_act .and. abs(fa) > abs(fb)) then
            ! Attempt secant interpolation.
            s = b_local - fb * (b_local - a_local) / (fb - fa)
            ! Check conditions for the interpolation step to be acceptable.
            if ((s < (3.0_pv * a_local + b_local) / 4.0_pv) .or. (s > b_local) .or. &
                (abs(s - b_local) >= abs(m) / 2.0_pv)) then
                ! If not acceptable, fall back to bisection.
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
    if (present(niter)) niter = iter
    status = 2
    if (show_res_local) then
        write(*,"(A, ES12.5)") &
            "Warning: Maximum iterations reached. Approximate root: ", root
    end if
end subroutine fun_zero


subroutine init_fun_zero_options(f0_options)
    ! General Description: Initialises a fun_zero_options variable to its default values.
    !   Calling this is optional when the type fields already carry default values
    !   (Fortran 2003 component initialisation), but is useful for explicit resets.
    ! Outputs:
    !   f0_options - The fun_zero_options variable to initialise.

    type(fun_zero_options), intent(out) :: f0_options
    f0_options%tolerance       = 1.0e-6_pv
    f0_options%maxiter         = 500
    f0_options%show_iterations = .false.
    f0_options%show_results    = .true.
end subroutine init_fun_zero_options


subroutine poly_roots(degree, coeffs, roots, status_flag, accuracy_hint)
    !********************************************************************
    ! General Description:
    !   Solves a monic algebraic equation (i.e., finds the roots of a polynomial)
    !   by constructing its companion matrix and computing its eigenvalues using
    !   the Householder QR algorithm.
    !
    !   The polynomial is defined as:
    !      coeffs(1)*x^(degree) + coeffs(2)*x^(degree-1) + ... + coeffs(degree)*x + coeffs(degree+1)
    !
    ! Inputs:
    !   degree      - An integer representing the degree of the monic polynomial.
    !   coeffs      - A real array of length (degree + 1) containing the polynomial
    !                 coefficients in order of decreasing powers. Note that coeffs(1)
    !                 must be nonzero.
    !
    ! Outputs:
    !   roots       - A vector of length (degree) containing the computed roots of the polynomial.
    !   status_flag - An integer return code indicating:
    !                   - -1: Invalid polynomial degree (<= 0)
    !                   - -2: Leading coefficient (coeffs(1)) is zero
    !                   -  0: Success
    !                   - >0: Error code returned from the QR eigenvalue solver.
    !
    ! Optional Outputs:
    !   accuracy_hint - A real value providing an accuracy estimate based on machine
    !                   precision, the degree of the polynomial, the total number of
    !                   QR iterations, and the Frobenius norm of the balanced companion
    !                   matrix.
    !
    ! Notes:
    !   This subroutine uses several helper routines:
    !     - build_companion_matrix: Constructs the companion matrix from the polynomial.
    !     - balance_companion: Balances the companion matrix to improve numerical stability.
    !     - hqr_eigen_hessenberg: Computes the eigenvalues of a Hessenberg matrix.
    !     - frobenius_norm_companion: Computes the Frobenius norm of the companion matrix.
    !********************************************************************
    implicit none

    ! Input arguments:
    integer, intent(in) :: degree
    real(pv), intent(in) :: coeffs(degree + 1)

    ! Output arguments:
    Complex(pv), intent(out) :: roots(degree)  ! Roots of the polynomial.
    integer, intent(out) :: status_flag        ! Return code (see header for details).

    ! Optional output:
    real(pv), intent(out), optional :: accuracy_hint   ! Accuracy estimate (if requested).

    ! Machine parameters:
    real(pv), parameter :: eps = epsilon(1.0_pv)

    ! Local working arrays:
    real(pv), allocatable :: companion_matrix(:,:)  ! Companion matrix (degree x degree)
    integer, allocatable  :: iteration_counts(:)    ! QR iteration counts for each eigenvalue.

    ! Local variables:
    integer :: i, total_iterations
    real(pv) :: frob_norm   ! Frobenius norm of the balanced companion matrix.
    real(pv) :: roots_real(degree), roots_imag(degree)  ! Real and imaginary parts of the roots.

    !--------------------------------------------------------------------
    ! Check for invalid input arguments.
    !--------------------------------------------------------------------
    if (degree <= 0) then
        ! Return -1 to distinguish from QR solver error codes (which are positive).
        status_flag = -1
        return
    end if
    if (coeffs(1) == 0.0_pv) then
        ! Leading coefficient is zero: cannot form a proper companion matrix.
        ! Return -2 to distinguish from QR solver error codes.
        status_flag = -2
        return
    end if

    !--------------------------------------------------------------------
    ! Allocate workspace arrays.
    !--------------------------------------------------------------------
    allocate(companion_matrix(degree, degree))
    allocate(iteration_counts(degree))

    !--------------------------------------------------------------------
    ! Build the companion matrix for the given polynomial.
    !--------------------------------------------------------------------
    call build_companion_matrix(degree, companion_matrix, coeffs)

    !--------------------------------------------------------------------
    ! Balance the companion matrix to improve numerical stability.
    !--------------------------------------------------------------------
    call balance_companion(degree, companion_matrix)

    !--------------------------------------------------------------------
    ! Compute the eigenvalues of the balanced companion matrix using the
    ! Householder QR algorithm. The real parts are stored in roots_real and
    ! the imaginary parts in roots_imag.
    ! The iteration_counts array stores the number of iterations taken for each root.
    !--------------------------------------------------------------------
    call hqr_eigen_hessenberg(degree, companion_matrix, roots_real, roots_imag, iteration_counts, status_flag)
    if (status_flag /= 0) then
        write(*, '(A,1X,I4)') 'Abnormal return from hqr_eigen_hessenberg. status_flag = ', status_flag
        if (status_flag == 1) write(*, '(A)') 'QR iteration did not converge.'
        deallocate(companion_matrix, iteration_counts)
        return
    else
        ! Store results in the complex roots array.
        do i = 1, degree
            roots(i) = cmplx(roots_real(i), roots_imag(i), kind=pv)
        end do
    end if

    !--------------------------------------------------------------------
    ! If an accuracy hint is requested, compute it.
    !   The hint is based on machine precision, the degree, the total number
    !   of QR iterations, and the Frobenius norm of the balanced companion matrix.
    !--------------------------------------------------------------------
    if (present(accuracy_hint)) then
        frob_norm = frobenius_norm_companion(degree, companion_matrix)

        total_iterations = 0
        do i = 1, degree
            if (iteration_counts(i) > 0) total_iterations = total_iterations + iteration_counts(i)
        end do

        accuracy_hint = eps * real(degree, pv) * real(total_iterations, pv) * frob_norm
    end if

    deallocate(companion_matrix, iteration_counts)

end subroutine poly_roots


subroutine build_companion_matrix(poly_degree, companion_matrix, coeffs)
    !********************************************************************
    ! General Description:
    !   Constructs the companion matrix for a polynomial. The eigenvalues
    !   of this companion matrix are the roots of the polynomial.
    !
    !   The polynomial is represented in the following form:
    !     P(x) = coeffs(1)*x^(poly_degree) + coeffs(2)*x^(poly_degree-1) + 
    !            ... + coeffs(poly_degree)*x + coeffs(poly_degree+1)
    !
    !   Note: The polynomial does not need to be monic. The companion matrix
    !         is built using normalized coefficients (dividing by coeffs(1)).
    !
    ! Inputs:
    !   poly_degree       - An integer representing the degree of the polynomial.
    !   coeffs            - A real array of length poly_degree+1 containing the
    !                       polynomial coefficients in decreasing order of power.
    !
    ! Outputs:
    !   companion_matrix  - A real matrix of dimensions (poly_degree, poly_degree)
    !                       that is the companion matrix of the polynomial.
    !********************************************************************
    implicit none

    integer, intent(in) :: poly_degree
    real(pv), intent(in) :: coeffs(poly_degree + 1)   ! Polynomial coefficients: highest to lowest degree.
    real(pv), intent(out) :: companion_matrix(poly_degree, poly_degree)

    integer :: i   ! Loop counter

    !--------------------------------------------------------------------
    ! Initialize the companion matrix to all zeros.
    !--------------------------------------------------------------------
    companion_matrix = 0.0_pv

    !--------------------------------------------------------------------
    ! Fill the subdiagonal with ones.
    ! The subdiagonal (just below the main diagonal) helps shift the basis
    ! so that the matrix's characteristic polynomial matches the input polynomial.
    !--------------------------------------------------------------------
    do i = 1, poly_degree - 1
        companion_matrix(i + 1, i) = 1.0_pv
    end do

    !--------------------------------------------------------------------
    ! Fill the last column with the normalized negative coefficients.
    ! The coefficients (except for the leading coefficient) are inserted in
    ! reverse order:
    !   - For i = poly_degree, companion_matrix(1, poly_degree) = -coeffs(poly_degree+1)/coeffs(1)
    !   - For i = poly_degree-1, companion_matrix(2, poly_degree) = -coeffs(poly_degree)/coeffs(1)
    !   - ...
    !   - For i = 1, companion_matrix(poly_degree, poly_degree) = -coeffs(2)/coeffs(1)
    !
    ! This arrangement ensures that the companion matrix accurately represents the
    ! characteristic polynomial of the original polynomial.
    !--------------------------------------------------------------------
    do i = poly_degree, 1, -1
        companion_matrix(poly_degree - i + 1, poly_degree) = -coeffs(i + 1) / coeffs(1)
    end do

end subroutine build_companion_matrix


subroutine balance_companion(matrix_dim, comp_matrix)
    !********************************************************************
    ! General Description:
    !   Balances an unsymmetric companion matrix by scaling selected rows
    !   and columns. This helps improve the accuracy and stability of
    !   eigenvalue and eigenvector computations.
    !
    !   This Fortran code is adapted from the Algol code "balance" in the
    !   paper:
    !       "Balancing a Matrix for the Calculation of Eigenvalues and
    !        Eigenvectors"
    !   by B.N. Parlett and C. Reinsch, Numerische Mathematik, 13, 293-304 (1969).
    !
    !   Note: Only the nonzero elements of the companion matrix are modified.
    !
    ! Inputs:
    !   matrix_dim  - Integer representing the dimension of the square matrix.
    !
    ! In/Out:
    !   comp_matrix - A real matrix of dimensions (matrix_dim, matrix_dim)
    !                 containing the companion matrix. This matrix is modified
    !                 in-place to achieve a balanced form.
    !********************************************************************
    implicit none

    ! Input/Output variables:
    integer, intent(in) :: matrix_dim
    real(pv), intent(inout) :: comp_matrix(matrix_dim, matrix_dim)

    ! Machine parameters:
    integer, parameter :: base = radix(1.0_pv)    ! Base (radix) of the floating-point system.
    integer, parameter :: base_squared = base**2  ! Square of the base.

    ! Local variables:
    integer :: i, j
    real(pv) :: colMagnitude      ! Magnitude from "column" nonzero elements.
    real(pv) :: rowMagnitude      ! Magnitude from "row" nonzero elements.
    real(pv) :: scalingFactor     ! Factor by which the matrix is scaled.
    real(pv) :: reciprocalScaling ! Reciprocal of scalingFactor.
    real(pv) :: totalMagnitude    ! Sum of rowMagnitude and colMagnitude.
    logical :: scaling_applied    ! Flag indicating if any scaling was done in an iteration.

    ! If the matrix is 1x1 or smaller, no balancing is needed.
    if (matrix_dim <= 1) return

    ! Begin the iterative balancing process.
    do  ! Main iteration loop.
        scaling_applied = .false.

        ! Loop over each index (corresponding to a row/column) in the matrix.
        do i = 1, matrix_dim

            !-------------------------------------------------------------
            ! Compute the "column magnitude" (colMagnitude) for the current index.
            ! For rows 1 to matrix_dim-1, the only nonzero element is the
            ! subdiagonal element at (i+1, i). For the last row (i = matrix_dim),
            ! the nonzero elements lie in the last column (rows 1 to matrix_dim-1).
            !-------------------------------------------------------------
            if (i /= matrix_dim) then
                colMagnitude = abs(comp_matrix(i + 1, i))
            else
                colMagnitude = 0.0_pv
                do j = 1, matrix_dim - 1
                    colMagnitude = colMagnitude + abs(comp_matrix(j, matrix_dim))
                end do
            end if

            !-------------------------------------------------------------
            ! Compute the "row magnitude" (rowMagnitude) for the current index.
            ! For the first row, the nonzero element is in the last column.
            ! For rows 2 to matrix_dim-1, the nonzero elements are the one
            ! immediately to the left (subdiagonal) and the one in the last column.
            ! For the last row, the only nonzero element is immediately to the left.
            !-------------------------------------------------------------
            if (i == 1) then
                rowMagnitude = abs(comp_matrix(1, matrix_dim))
            elseif (i /= matrix_dim) then
                rowMagnitude = abs(comp_matrix(i, i - 1)) + abs(comp_matrix(i, matrix_dim))
            else
                rowMagnitude = abs(comp_matrix(matrix_dim, matrix_dim - 1))
            end if

            ! Only proceed if both magnitudes are nonzero.
            if (colMagnitude /= 0.0_pv .and. rowMagnitude /= 0.0_pv) then

                !-------------------------------------------------------------
                ! Determine an appropriate scaling factor to balance the row and
                ! column magnitudes.
                !
                ! Initialize scalingFactor to 1.0 and record the original total
                ! magnitude.
                !-------------------------------------------------------------
                scalingFactor = 1.0_pv
                totalMagnitude = colMagnitude + rowMagnitude

                ! Increase scalingFactor until the column magnitude is no longer too small.
                do
                    if (colMagnitude >= rowMagnitude / base) exit
                    scalingFactor = scalingFactor * base
                    colMagnitude = colMagnitude * base_squared
                end do

                ! Decrease scalingFactor until the column magnitude is no longer too large.
                do
                    if (colMagnitude < rowMagnitude * base) exit
                    scalingFactor = scalingFactor / base
                    colMagnitude = colMagnitude / base_squared
                end do

                !-------------------------------------------------------------
                ! Check if the scaling improves the balance by more than 5%.
                ! If so, apply the scaling. The test compares the scaled sum
                ! of magnitudes with 0.95 times the original total magnitude.
                !-------------------------------------------------------------
                if ((colMagnitude + rowMagnitude) / scalingFactor < 0.95_pv * totalMagnitude) then
                    reciprocalScaling = 1.0_pv / scalingFactor
                    scaling_applied = .true.

                    ! Scale the nonzero elements in the current row:
                    ! For the first row, only the last column element is scaled.
                    ! For subsequent rows, scale the subdiagonal element (left) and
                    ! the element in the last column.
                    if (i == 1) then
                        comp_matrix(1, matrix_dim) = comp_matrix(1, matrix_dim) * reciprocalScaling
                    else
                        comp_matrix(i, i - 1) = comp_matrix(i, i - 1) * reciprocalScaling
                        comp_matrix(i, matrix_dim) = comp_matrix(i, matrix_dim) * reciprocalScaling
                    end if

                    ! Scale the nonzero elements in the corresponding column:
                    ! For rows other than the last, scale the subdiagonal element below.
                    ! For the last row, scale every element in column i.
                    if (i /= matrix_dim) then
                        comp_matrix(i + 1, i) = comp_matrix(i + 1, i) * scalingFactor
                    else
                        do j = 1, matrix_dim
                            comp_matrix(j, i) = comp_matrix(j, i) * scalingFactor
                        end do
                    end if
                end if
            end if
        end do  ! End loop over i.

        ! If any scaling was applied in this iteration, repeat the process.
        if (scaling_applied) cycle
        exit  ! Exit the main loop when no further scaling is needed.
    end do  ! End of the main iteration loop.

end subroutine balance_companion


function frobenius_norm_companion(matrix_dim, comp_matrix) result(frobenius_norm)
    !********************************************************************
    ! General Description:
    !   Computes the Frobenius norm of a companion-like matrix.
    !
    !   The companion matrix is assumed to have nonzero entries only in:
    !     - The subdiagonal positions (i.e., element (j+1, j) for j = 1 to matrix_dim-1)
    !     - The last column (i.e., element (i, matrix_dim) for i = 1 to matrix_dim)
    !
    !   The Frobenius norm is defined as:
    !       ||A||_F = sqrt( sum_i sum_j (A(i,j)**2) )
    !
    !   Since only the nonzero entries are stored in the above locations,
    !   we only sum their squares.
    !
    ! Inputs:
    !   matrix_dim    - Integer representing the dimension (order) of the matrix.
    !   comp_matrix   - Real matrix of size (matrix_dim, matrix_dim) that represents
    !                   the companion matrix.
    !
    ! Output:
    !   frobenius_norm - The computed Frobenius norm of the companion matrix.
    !********************************************************************
    implicit none

    ! Input arguments:
    integer, intent(in) :: matrix_dim
    real(pv), intent(in) :: comp_matrix(matrix_dim, matrix_dim)

    ! Function result:
    real(pv) :: frobenius_norm

    ! Local variables:
    integer :: i, j
    real(pv) :: sumSquares   ! Accumulator for the sum of squares

    ! Initialize the sum of squares to zero.
    sumSquares = 0.0_pv

    !--------------------------------------------------------------------
    ! Sum the squares of the subdiagonal elements.
    ! These are located at positions (j+1, j) for j = 1 to matrix_dim-1.
    !--------------------------------------------------------------------
    do j = 1, matrix_dim - 1
        sumSquares = sumSquares + comp_matrix(j + 1, j)**2
    end do

    !--------------------------------------------------------------------
    ! Sum the squares of the elements in the last column.
    ! These are located at positions (i, matrix_dim) for i = 1 to matrix_dim.
    !--------------------------------------------------------------------
    do i = 1, matrix_dim
        sumSquares = sumSquares + comp_matrix(i, matrix_dim)**2
    end do

    ! Compute the Frobenius norm as the square root of the total sum.
    frobenius_norm = sqrt(sumSquares)

end function frobenius_norm_companion


subroutine hqr_eigen_hessenberg(order, hessenberg, eigen_real, eigen_imag, iteration_count, status_flag)
    !********************************************************************
    ! General Description:
    !   Computes the eigenvalues of a real upper Hessenberg matrix using the
    !   Householder QR algorithm. This routine is adapted from the Algol code
    !   "hqr" in the paper:
    !
    !       "The QR Algorithm for Real Hessenberg Matrices"
    !       by R.S. Martin, G. Peters and J.H. Wilkinson,
    !       Numer. Math. 14, 219-231 (1970).
    !
    !   The algorithm reduces the matrix and then finds its eigenvalues. The
    !   real parts of the eigenvalues are stored in eigen_real and the imaginary
    !   parts in eigen_imag. The routine also records (in iteration_count) the
    !   number of iterations required for convergence for each eigenvalue.
    !
    ! Inputs:
    !   order       - The order (dimension) of the Hessenberg matrix.
    !   hessenberg  - A real matrix (order x order) in upper Hessenberg form.
    !
    ! In/Out:
    !   iteration_count - An integer array of length order that records the number
    !                     of iterations used for each eigenvalue.
    !
    ! Outputs:
    !   eigen_real  - Array (length order) containing the real parts of the computed
    !                 eigenvalues.
    !   eigen_imag  - Array (length order) containing the imaginary parts of the computed
    !                 eigenvalues.
    !   status_flag - Integer flag: 0 indicates success; 1 indicates failure to converge.
    !********************************************************************
    implicit none

    ! Arguments:
    integer, intent(in) :: order
    real(pv), intent(inout) :: hessenberg(order, order)
    real(pv), intent(out) :: eigen_real(order)
    real(pv), intent(out) :: eigen_imag(order)
    integer, intent(inout) :: iteration_count(order)
    integer, intent(out) :: status_flag

    ! Local variables:
    integer :: i, j, k, l, m
    integer :: one_less_order, current_order, num_iterations
    real(pv) :: eigen_approx        ! Approximation taken from the (current_order,current_order) diagonal.
    real(pv) :: next_diag           ! hessenberg(one_less_order, one_less_order)
    real(pv) :: offdiag_product     ! Product: hessenberg(current_order, one_less_order)*hessenberg(one_less_order, current_order)
    real(pv) :: shift_accum         ! Accumulates shifts over iterations.
    real(pv) :: subdiag_sum         ! Sum of magnitudes of selected subdiagonal elements.
    real(pv) :: temp_p = 0.0_pv, temp_q = 0.0_pv, temp_r = 0.0_pv  ! Temporary variables for QR computations.
    real(pv) :: temp_val              ! Temporary values used during exceptional shifts.
    real(pv) :: temp_r_old              ! Saves the first-column HQR vector element before temp_r is overwritten.
    real(pv) :: p_save, q_save, r_save ! Householder vector from the m-search, used at k==m in the QR step.
    real(pv) :: x, y, z           ! Used in various computations (e.g., shift approximations)
    real(pv) :: s                 ! Used for computing a norm in the QR step.
    logical :: small_subdiag_found, not_last

    ! Maximum iterations allowed per eigenvalue.
    integer, parameter :: max_iterations = 500
    ! Use machine epsilon for the subdiagonal negligibility test.
    real(pv), parameter :: eps = epsilon(1.0_pv)

    ! Initialization.
    current_order = order
    status_flag = 0
    shift_accum = 0.0_pv

    main_loop: do
        ! If the current submatrix is empty, all eigenvalues have been found.
        if (current_order == 0) return

        num_iterations = 0
        one_less_order = current_order - 1

        !--------------------------------------------------------------------
        ! Search for a small subdiagonal element. We look for an index l (2<=l<=current_order)
        ! such that the subdiagonal element hessenberg(l, l-1) is negligible relative to
        ! the sum |hessenberg(l-1, l-1)| + |hessenberg(l, l)|.
        !--------------------------------------------------------------------
        small_subdiag_found = .false.
        do l = current_order, 2, -1
            if ( abs(hessenberg(l, l-1)) <= eps * ( abs(hessenberg(l-1, l-1)) + abs(hessenberg(l, l)) ) ) then
                small_subdiag_found = .true.
                exit
            end if
        end do
        if (.not. small_subdiag_found) then
            l = 1
        end if

        ! Set eigen_approx to the bottom-right diagonal element.
        eigen_approx = hessenberg(current_order, current_order)

        if (l == current_order) then
            !----------------------------------------------------------------
            ! A single eigenvalue has converged.
            !----------------------------------------------------------------
            eigen_real(current_order) = eigen_approx + shift_accum
            eigen_imag(current_order) = 0.0_pv
            iteration_count(current_order) = num_iterations
            current_order = one_less_order
            cycle main_loop
        else
            !----------------------------------------------------------------
            ! If l == one_less_order, then a 2x2 block is isolated and yields
            ! a pair of eigenvalues.
            !----------------------------------------------------------------
            next_diag = hessenberg(one_less_order, one_less_order)
            offdiag_product = hessenberg(current_order, one_less_order) * &
                              hessenberg(one_less_order, current_order)
            if (l == one_less_order) then
                ! Compute the eigenvalues from the 2x2 block.
                temp_p = (next_diag - eigen_approx) / 2.0_pv
                temp_q = temp_p**2 + offdiag_product
                next_diag = sqrt( abs(temp_q) )
                iteration_count(current_order) = -num_iterations
                iteration_count(one_less_order) = num_iterations
                eigen_approx = eigen_approx + shift_accum
                if (temp_q > 0.0_pv) then
                    !----------------------------------------------------------------
                    ! Real eigenvalue pair.
                    if (temp_p < 0.0_pv) then
                        next_diag = -next_diag
                    end if
                    next_diag = temp_p + next_diag
                    eigen_real(one_less_order) = eigen_approx + next_diag
                    eigen_real(current_order) = eigen_approx - offdiag_product / next_diag
                    eigen_imag(one_less_order) = 0.0_pv
                    eigen_imag(current_order) = 0.0_pv
                else
                    !----------------------------------------------------------------
                    ! Complex conjugate eigenvalue pair.
                    eigen_real(one_less_order) = eigen_approx + temp_p
                    eigen_real(current_order) = eigen_approx + temp_p
                    eigen_imag(one_less_order) = next_diag
                    eigen_imag(current_order) = -next_diag
                end if
                current_order = current_order - 2
                cycle main_loop
            else
                !----------------------------------------------------------------
                ! No eigenvalue has converged yet.
                ! If we have reached the maximum number of iterations, then the
                ! algorithm fails.
                !----------------------------------------------------------------
                if (num_iterations == max_iterations) then
                    status_flag = 1
                    return
                end if

                !----------------------------------------------------------------
                ! Exceptional shift: after 10 or 20 iterations, modify the shift.
                !----------------------------------------------------------------
                if (num_iterations == 10 .or. num_iterations == 20) then
                    shift_accum = shift_accum + eigen_approx
                    do i = 1, current_order
                        hessenberg(i, i) = hessenberg(i, i) - eigen_approx
                    end do
                    subdiag_sum = abs(hessenberg(current_order, one_less_order)) + &
     &                             abs(hessenberg(one_less_order, (current_order - 1)))
                    next_diag = 0.75_pv * subdiag_sum
                    eigen_approx = next_diag
                    offdiag_product = -0.4375_pv * subdiag_sum**2
                end if

                num_iterations = num_iterations + 1

                !----------------------------------------------------------------
                ! Look for two consecutive small subdiagonal elements by scanning
                ! upward from index (current_order-2) down to l.
                !----------------------------------------------------------------
                do m = current_order - 2, l, -1
                    z = hessenberg(m, m)
                    temp_p = eigen_approx - z
                    temp_q = next_diag - z
                    ! Compute the first component of the Householder starting vector.
                    ! Save it in temp_r_old before temp_r is reused below.
                    temp_r_old = ( temp_p * temp_q - offdiag_product ) / hessenberg(m+1, m) + &
                                 hessenberg(m, m+1)
                    temp_val = hessenberg(m+1, m+1) - z - temp_p - temp_q
                    ! Third component of the Householder starting vector.
                    temp_r = hessenberg(m+2, m+1)
                    ! Normalise all three components together (previously temp_r_old was
                    ! lost and temp_r appeared twice, giving a wrong normalisation).
                    subdiag_sum = abs(temp_r_old) + abs(temp_val) + abs(temp_r)
                    temp_p = temp_r_old / subdiag_sum
                    temp_q = temp_val   / subdiag_sum
                    temp_r = temp_r     / subdiag_sum
                    if (m == l) exit
                    if ( abs(hessenberg(m, m-1))*(abs(temp_q)+abs(temp_r)) <= &
                         eps*abs(temp_p)*( abs(hessenberg(m-1, m-1)) + abs(z) + &
                         abs(hessenberg(m+1, m+1)) ) ) then
                        exit
                    end if
                end do

                ! Save the Householder vector from the m-search for use at k==m
                ! in the QR step below (avoids the out-of-bounds read when m==1).
                p_save = temp_p
                q_save = temp_q
                r_save = temp_r

                !----------------------------------------------------------------
                ! Set to zero any elements that have become negligibly small,
                ! thus preserving the Hessenberg structure.
                !----------------------------------------------------------------
                do i = m+2, current_order
                    hessenberg(i, i-2) = 0.0_pv
                end do
                do i = m+3, current_order
                    hessenberg(i, i-3) = 0.0_pv
                end do

                !----------------------------------------------------------------
                ! Perform a double QR step on the submatrix spanning rows l to
                ! current_order and columns m to current_order.
                !----------------------------------------------------------------
                do k = m, one_less_order
                    not_last = (k /= one_less_order)
                    if (k /= m) then
                        temp_p = hessenberg(k, k-1)
                        temp_q = hessenberg(k+1, k-1)
                        if (not_last) then
                            temp_r = hessenberg(k+2, k-1)
                        else
                            temp_r = 0.0_pv
                        end if
                        x = abs(temp_p) + abs(temp_q) + abs(temp_r)
                        if (x == 0.0_pv) cycle
                        temp_p = temp_p / x
                        temp_q = temp_q / x
                        temp_r = temp_r / x
                    else
                        ! Use the Householder vector saved from the m-search.
                        ! Reading hessenberg(k, k-1) at k==m (possibly m==1, i.e., column 0)
                        ! would be an out-of-bounds access and gives the wrong values.
                        temp_p = p_save
                        temp_q = q_save
                        temp_r = r_save
                        x = abs(temp_p) + abs(temp_q) + abs(temp_r)
                        if (x == 0.0_pv) cycle
                        temp_p = temp_p / x
                        temp_q = temp_q / x
                        temp_r = temp_r / x
                    end if
                    s = sqrt( temp_p**2 + temp_q**2 + temp_r**2 )
                    if (temp_p < 0.0_pv) s = -s
                    if (k /= m) then
                        hessenberg(k, k-1) = -s * x
                    elseif (l /= m) then
                        hessenberg(k, k-1) = -hessenberg(k, k-1)
                    end if
                    temp_p = temp_p + s
                    x = temp_p / s
                    y = temp_q / s
                    z = temp_r / s
                    temp_q = temp_q / temp_p
                    temp_r = temp_r / temp_p

                    !-----------------------
                    ! Row modification:
                    !-----------------------
                    do j = k, current_order
                        temp_val = hessenberg(k, j) + temp_q * hessenberg(k+1, j)
                        if (not_last) then
                            temp_val = temp_val + temp_r * hessenberg(k+2, j)
                            hessenberg(k+2, j) = hessenberg(k+2, j) - temp_val * z
                        end if
                        hessenberg(k+1, j) = hessenberg(k+1, j) - temp_val * y
                        hessenberg(k, j)   = hessenberg(k, j)   - temp_val * x
                    end do

                    !-----------------------
                    ! Column modification:
                    !-----------------------
                    if (k + 3 < current_order) then
                        j = k + 3
                    else
                        j = current_order
                    end if
                    do i = l, j
                        temp_val = x * hessenberg(i, k) + y * hessenberg(i, k+1)
                        if (not_last) then
                            temp_val = temp_val + z * hessenberg(i, k+2)
                            hessenberg(i, k+2) = hessenberg(i, k+2) - temp_val * temp_r
                        end if
                        hessenberg(i, k+1) = hessenberg(i, k+1) - temp_val * temp_q
                        hessenberg(i, k)   = hessenberg(i, k)   - temp_val
                    end do
                end do

                cycle  ! End of this QR step; restart main loop.
            end if
        end if
    end do main_loop

end subroutine hqr_eigen_hessenberg
  

!-----------------------------------------------------------------------------------------------
subroutine newton_raphson(f_pointer, df_pointer, x0, status, root, niter, nr_options)
    !********************************************************************
    ! General Description:
    !   Finds a root of f(x) using Newton-Raphson iteration:
    !       x_{n+1} = x_n - f(x_n) / f'(x_n)
    !
    !   This method converges quadratically near a simple root but requires
    !   an analytic (or numerical) derivative and a good initial guess.
    !
    ! Inputs:
    !   f_pointer   - Procedure pointer to f(x).
    !   df_pointer  - Procedure pointer to f'(x).
    !   x0          - Initial guess for the root.
    ! Optional Inputs:
    !   nr_options  - A fun_zero_options struct controlling tolerance, maxiter,
    !                 show_iterations, and show_results. Defaults are applied when absent.
    ! Outputs:
    !   status      - Integer return code:
    !                   0 - Root found within tolerance.
    !                   1 - Derivative is zero (division by zero would occur).
    !                   2 - Maximum iterations reached without convergence.
    !                   3 - f(x) or f'(x) is NaN during iteration.
    !   root        - The computed (approximate) root.
    ! Optional Outputs:
    !   niter       - Number of iterations performed.
    !********************************************************************
    implicit none

    procedure(sv_function), pointer, intent(in) :: f_pointer
    procedure(sv_function), pointer, intent(in) :: df_pointer
    real(pv), intent(in)  :: x0
    integer,  intent(out) :: status
    real(pv), intent(out) :: root
    integer,  intent(out), optional :: niter
    type(fun_zero_options), intent(in), optional :: nr_options

    ! Local variables.
    real(pv) :: x_cur, fx, dfx, dx
    integer  :: iter, maxit_local
    real(pv) :: tol_local
    logical  :: show_iter_local, show_res_local
    type(fun_zero_options) :: opts

    !-- Merge options with defaults -------------------------------------------
    call init_fun_zero_options(opts)
    if (present(nr_options)) then
        if (nr_options%tolerance > 0.0_pv) opts%tolerance = nr_options%tolerance
        if (nr_options%maxiter   > 0)       opts%maxiter   = nr_options%maxiter
        opts%show_iterations = nr_options%show_iterations
        opts%show_results    = nr_options%show_results
    end if

    tol_local       = opts%tolerance
    maxit_local     = opts%maxiter
    show_iter_local = opts%show_iterations
    show_res_local  = opts%show_results

    !-- Iteration -------------------------------------------------------------
    x_cur = x0

    do iter = 1, maxit_local
        fx  = f_pointer(x_cur)
        dfx = df_pointer(x_cur)

        ! Guard against NaN in function or derivative evaluations.
        if (IEEE_IS_NAN(fx) .or. IEEE_IS_NAN(dfx)) then
            write(*,"(A, ES12.4)") "Error: NaN encountered at x =", x_cur
            status = 3
            root = x_cur
            if (present(niter)) niter = iter
            return
        end if

        ! Guard against zero derivative (method is undefined).
        if (dfx == 0.0_pv) then
            write(*,"(A, ES12.4)") "Error: Derivative is zero at x =", x_cur
            status = 1
            root = x_cur
            if (present(niter)) niter = iter
            return
        end if

        dx    = fx / dfx
        x_cur = x_cur - dx

        if (show_iter_local) then
            write(*,"(A, I5, A, ES12.4, A, ES12.4)") &
                "Iteration:", iter, "    x =", x_cur, "    f(x) =", fx
        end if

        ! Convergence test: the correction |dx| is smaller than the tolerance.
        if (abs(dx) <= tol_local + epsilon(1.0_pv) * abs(x_cur)) then
            root = x_cur
            status = 0
            if (present(niter)) niter = iter
            if (show_res_local) then
                print*, " "
                write(*,"(A, ES12.5, A, I0, A)") &
                    "Root found: ", root, " after ", iter, " iterations."
                print*, "Status: ", status
                print*, " "
            end if
            return
        end if
    end do

    ! Maximum iterations reached.
    root = x_cur
    status = 2
    if (present(niter)) niter = iter
    if (show_res_local) then
        write(*,"(A, ES12.5)") &
            "Warning: Maximum iterations reached. Approximate root: ", root
    end if

end subroutine newton_raphson

!-----------------------------------------------------------------------------------------------
end module MOD_Roots
