
! This file handles some input and output functions for the solver
! the big one being  console output to the GUI which can be seen below


    module solver_flags
        use iso_c_binding, only: C_INT8_T
        implicit none
        INTEGER(C_INT8_T) :: stopit
    end module solver_flags

    subroutine set_stopit_flag() bind(C, name="set_stopit_flag")
        use iso_c_binding
        use solver_flags
        implicit none
        ! Update the stopit flag in Fortran (use 1/0 integer)
        stopit = 1
        if (stopit /= 0) then
            print *, "Stopit flag set to true in Fortran."
        else
            print *, "Stopit flag set to false in Fortran."
        end if
    end subroutine set_stopit_flag

    module io_module
        use iso_c_binding
        use types
        use conversion
        implicit none
    
        interface
            subroutine emit_console_signal(text, length) bind(C, name="emit_console_signal")
                use iso_c_binding, only: c_char, c_int
                character(kind=c_char), intent(in) :: text(*)  ! Pass the string as an array
                integer(c_int), intent(in), value :: length                ! Length of the string
            end subroutine emit_console_signal

            subroutine emit_grid_signal(g) bind(C, name="emit_grid_signal")
                use types
                type(t_grid_c), intent(in) :: g
            end subroutine emit_grid_signal

            subroutine emit_conv_point_signal(cp) bind(C, name="emit_conv_point_signal")
                use types
                type(t_conv_point), intent(in) :: cp
            end subroutine emit_conv_point_signal

            subroutine emit_grid_vector_signal(g, length) bind(C, name="emit_grid_vector_signal")
                use iso_c_binding, only: c_ptr, c_int
                use types
                ! g is a pointer to an array of t_grid_c
                type(c_ptr), intent(in), value :: g
                integer(c_int), intent(in), value :: length
            end subroutine emit_grid_vector_signal

            subroutine emit_mesh(mesh, length) bind(C, name="emit_mesh")
                use iso_c_binding, only: c_ptr, c_int
                use types
                ! mesh is a pointer to an array of cell
                type(c_ptr), intent(in), value :: mesh
                integer(c_int), intent(in), value :: length
            end subroutine emit_mesh
            
        end interface
    
    contains
    
        subroutine write_to_qt(text)
            character(len=*), intent(in) :: text
            character(kind=c_char), allocatable :: c_text(:)
            integer :: length

            ! still want it command line outout
            write(6, *) trim(text)
    
            ! trim and get the length of the string
            length = len_trim(text)
    
            ! allocate and copy fortran string to be C-compatible
            allocate(c_text(0:length-1))
            c_text = transfer(trim(text), c_text)
    
            call emit_console_signal(c_text, length)    
            ! deallocate for safety
            deallocate(c_text)
        end subroutine write_to_qt

        subroutine grid_to_qt(g)

            type(t_grid), intent(inout) :: g
            type(t_grid_c) :: g_c

            call grid_to_c(g, g_c)
            call emit_grid_signal(g_c)

            ! testing has shown that this delay is required for
            ! the main thread to copy the data before fortran 
            ! steams ahead and overwrites the data
            ! call sleepqq(20) ! no longer required as signal is blocking now

        end subroutine grid_to_qt

        subroutine conv_point_to_qt(cp)

            type(t_conv_point), intent(in) :: cp
            call emit_conv_point_signal(cp)

        end subroutine conv_point_to_qt

        subroutine grids_to_qt(g, length)

            use iso_c_binding, only: c_loc

            type(t_grid), allocatable, target, intent(in) :: g(:)
            type(t_grid), pointer :: g_ptr(:)
            integer, intent(in) :: length
            type(t_grid_c), pointer :: g_c(:)
            integer :: i

            allocate(g_c(length))

            g_ptr => g

            do i = 1, length
                call grid_to_c(g_ptr(i), g_c(i))
            end do

            !do i = 1, g(1)%nj
            !    write(6,*) 'x(ni,1) = ', g(1)%x(g(1)%ni, i)
            !end do

            call emit_grid_vector_signal(c_loc(g_c), length)

        end subroutine grids_to_qt
    
    end module io_module
