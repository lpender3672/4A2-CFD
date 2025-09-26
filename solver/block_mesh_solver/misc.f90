
module debug
    use io_module
    use types
    implicit none
    contains
    
    subroutine find_NaN(x, msg_bfr)

        real, intent(in), dimension(:,:) :: x
        character(len=64), intent(out) :: msg_bfr

        integer :: i, j, ni, nj

        ni = size(x,1)
        nj = size(x,2)

        do i = 1, ni
            do j = 1, nj
                if(isnan(x(i,j))) then
                    write(msg_bfr, *) 'NaN found at i = ', i, ', j = ', j
                    return
                end if
            end do
        end do
  end subroutine find_NaN

  subroutine find_all_NaN(g, msg_bfr)

        type(t_grid), intent(in) :: g
        character(len=64), intent(out) :: msg_bfr
        
        write(msg_bfr,*) 'ro'
        call write_to_qt(msg_bfr)
        call find_NaN(g%ro,msg_bfr)
        call write_to_qt(msg_bfr)
        write(msg_bfr,*) 'roe'
        call write_to_qt(msg_bfr)
        call find_NaN(g%roe,msg_bfr)
        call write_to_qt(msg_bfr)
        write(msg_bfr,*) 'rovx'
        call write_to_qt(msg_bfr)
        call find_NaN(g%rovx,msg_bfr)
        call write_to_qt(msg_bfr)
        write(msg_bfr,*) 'rovy'
        call write_to_qt(msg_bfr)
        call find_NaN(g%rovy,msg_bfr)
        call write_to_qt(msg_bfr)
        
  end subroutine find_all_NaN
end module debug