
module cfill_types

    use iso_c_binding, only: c_double, c_int, c_long_long
    implicit none
    
    integer, parameter :: rk = 8
    integer, parameter :: i4 = selected_int_kind(9)
    integer, parameter :: i8 = selected_int_kind(18)

    type :: cell
        real(rk)    :: x, y, z
        integer     :: level
        integer(i8) :: key
        integer(i4) :: id
    end type cell


    type, bind(c) :: cell2d

        real(c_double)    :: xmin, xmax, ymin, ymax
        integer(c_int)    :: level
        integer(c_long_long) :: key
        integer(c_int)    :: id
    end type cell2d

end module cfill_types