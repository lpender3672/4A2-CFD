
! This file handles some input and output functions for the solver
! the big one being  console output to the GUI which can be seen below

    module io_module
        use iso_c_binding
        use types
        use conversion
        implicit none
    
        interface
            subroutine emit_console_signal(text, length) bind(C, name="emit_console_signal")
                use iso_c_binding, only: c_char, c_int
                character(kind=c_char), intent(in) :: text(*)  ! Pass the string as an array
                integer(c_int), value :: length                ! Length of the string
            end subroutine emit_console_signal

            subroutine emit_grid_signal(g) bind(C, name="emit_grid_signal")
                use types
                type(t_grid_c), intent(in) :: g
            end subroutine emit_grid_signal
        end interface
    
    contains
    
        subroutine write_to_qt(text)
            character(len=*), intent(in) :: text
            character(kind=c_char), allocatable :: c_text(:)
            integer :: length

            ! still want it command line outout
            write(6, *) text
    
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

            type(t_grid), intent(in) :: g
            type(t_grid_c) :: g_c

            call grid_to_c(g, g_c)
            call emit_grid_signal(g_c)

        end subroutine grid_to_qt
    
    end module io_module
