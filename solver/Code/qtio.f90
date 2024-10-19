
! This file handles some input and output functions for the solver
! the big one being  console output to the GUI which can be seen below

    module io_module
        use iso_c_binding
        implicit none
    
        interface
            subroutine qt_console_write(text, length) bind(C, name="qt_console_write")
                use iso_c_binding, only: c_char, c_int
                character(kind=c_char), intent(in) :: text(*)  ! Pass the string as an array
                integer(c_int), value :: length                ! Length of the string
            end subroutine qt_console_write
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
    
            call qt_console_write(c_text, length)    
            ! deallocate for safety
            deallocate(c_text)
        end subroutine write_to_qt
    
    end module io_module
