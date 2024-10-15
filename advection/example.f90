! example.f90
module example_module
    use iso_c_binding, only: C_INT
    implicit none
contains
    subroutine fortran_hello() bind(C, name="fortran_hello")
        print *, 'Hello from Fortran! is this compiled on debug'
    end subroutine fortran_hello

    subroutine add_numbers(a, b, result) bind(C, name="add_numbers")
        ! Input arguments
        integer(C_INT), intent(in) :: a, b
        ! Output argument
        integer(C_INT), intent(out) :: result
        ! Perform the addition
        result = a + b
    end subroutine add_numbers
    
end module example_module