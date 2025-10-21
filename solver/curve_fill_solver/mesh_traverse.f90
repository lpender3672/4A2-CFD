module mesh_traverse

    implicit none
    contains
    
end module mesh_traverse

module mesh_alloc
    use mesh_utils
    use types
    implicit none

    contains

    subroutine calc_ncells(n, m, airfoil_poly, ILOD, extra_global_levels, ncells)

        implicit none
        integer, intent(in) :: n, m
        real(8), intent(in) :: airfoil_poly(:, :)       ! (Nf,2)
        integer, intent(in) :: ILOD(n,m)
        integer, intent(in), optional :: extra_global_levels
        integer, intent(out) :: ncells

        ncells = 1
        ! traverse incrementing ncells

    end subroutine calc_ncells


    subroutine calc_nindex(n, m, cells, nindex)

        implicit none
        integer, intent(in) :: n, m
        type(cell2d), intent(inout) :: cells(:)
        integer, intent(out) :: nindex

        ! traverse neighbouring cells incrementing nindex

        ! this will fill up neigh_index of all neighbours
        ! then theres two arrays of length ncells neigh_start and neigh_count

    end subroutine calc_nindex

end module mesh_alloc