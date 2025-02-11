program test
    use MOD_Select_Kind, only: pv
    use MOD_Roots
    implicit none

    real(pv) :: bounds(2), root
    integer :: status
    procedure(sv_function), pointer :: f_pointer => null()
    type(fun_zero_options) :: options

    bounds = [-2.0, 2.0]
    f_pointer => my_function  ! Pointer to the function to be solved

    options%show_iterations = .true.
    options%show_results    = .true.
    options%maxiter         = 100
    options%tolerance       = 1.0E-6_pv

    call fun_zero(f_pointer, bounds, status, root, f0_options = options)


    contains

    function my_function(x) result(y)
        real(pv), intent(in) :: x
        real(pv) :: y
        real(pv) :: pi

        pi = 3.141592653_pv
        y = -exp(-x) + 2.0 * x**2 / (exp(4.0*x)) - pi
    end function my_function


end program test