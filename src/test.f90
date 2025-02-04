program test
    use MOD_Select_Kind, only: pv
    use MOD_Roots
    implicit none

    real(pv) :: bounds(2), root
    integer :: status
    procedure(sv_function), pointer :: f_pointer => null()
    type(fun_zero_options) :: options

    real(pv) :: coeff(3)
    complex(pv) :: pzeros(2)
    integer :: deg



    bounds = [-2.0, 2.0]
    f_pointer => f

    options%show_iterations = .true.
    options%show_results    = .true.
    options%maxiter         = 100
    options%tolerance       = 1.0E-6_pv

    !call fun_zero(f_pointer, bounds, status, root, f0_options = options)
    coeff = [1.0, -1.0, -3.0]
    deg   = 2
    call poly_roots(coeff, deg, pzeros, status)
    if(status == 0) then
        print *, 'Roots: ', pzeros
    else
        print *, 'No roots found'
    end if


    contains

    function f(x) result(y)
        real(pv), intent(in) :: x
        real(pv) :: y
        real(pv) :: pi

        pi = 3.141592653_pv
        y = -exp(-x) + 2.0 * x**2 / (exp(4.0*x)) - pi
        !y = -log(1.0/x) + exp(-1.0 / x)
    end function f


end program test