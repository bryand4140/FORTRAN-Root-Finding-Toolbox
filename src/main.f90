program main
    use MOD_Select_Kind, only: pv
    use MOD_Roots
    implicit none

    integer, parameter :: degree = 4
    real(pv) :: coeffs(degree + 1)
    integer :: status_flag, i
    real(pv) :: roots_real(degree), roots_imag(degree)
    complex(pv),allocatable :: roots(:)


    print*, "============================================"
    print*, "             ---> Main < ---"

    allocate(roots(degree))
    coeffs = [1.0, 2.0, 3.0, 4.0, 5.0]



    call poly_roots(degree, coeffs, roots, status_flag)

    if(status_flag == 0) then
        do i = 1,size(roots)
            if (aimag(roots(i)) == 0.0) then
                write(*,'(F10.4)') real(roots(i))
            else
                write(*,'(F10.4, 1X, F10.4, 1X, A)') real(roots(i)), aimag(roots(i)), 'i'
            end if
        end do
    else
        print*, "Error: ", status_flag
    end if






    print*, "============================================"
    print*, " " !Spacer


end program main

