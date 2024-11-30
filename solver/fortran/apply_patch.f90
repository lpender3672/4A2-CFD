

module patches
    
        use types
        implicit none
    
        contains

subroutine apply_patch(gs, p)
!   Apply single patch to the grids

    use types
    implicit none

    type(t_grid), allocatable, intent(inout) :: gs(:)
    type(t_match), intent(in) :: p

    integer :: q
    real :: ro_avg, roe_avg, rovx_avg, rovy_avg

    type(t_grid) :: g1, g2

    ! avg primary flow variables between gs(p%n_1) and gs(p%n_2)
    ! dont know if patch is in i or j direction
    ! luckily patch contains a list of i and j values which is inefficient
    ! for the square grid boundaries, but i guess works for non-square grids
    ! but non square grids would not be supported by the implemented square cell

    ! anyway, using the list of indicies

    do q = 1, p%nk

        ro_avg = 0.5 * (gs(p%n_1)%ro(p%i_1(q), p%j_1(q)) + gs(p%n_2)%ro(p%i_2(q), p%j_2(q)))
        gs(p%n_1)%ro(p%i_1(q), p%j_1(q)) = ro_avg
        gs(p%n_2)%ro(p%i_2(q), p%j_2(q)) = ro_avg

        roe_avg = 0.5 * (gs(p%n_1)%roe(p%i_1(q), p%j_1(q)) + gs(p%n_2)%roe(p%i_2(q), p%j_2(q)))
        gs(p%n_1)%roe(p%i_1(q), p%j_1(q)) = roe_avg
        gs(p%n_2)%roe(p%i_2(q), p%j_2(q)) = roe_avg

        rovx_avg = 0.5 * (gs(p%n_1)%rovx(p%i_1(q), p%j_1(q)) + gs(p%n_2)%rovx(p%i_2(q), p%j_2(q)))
        gs(p%n_1)%rovx(p%i_1(q), p%j_1(q)) = rovx_avg
        gs(p%n_2)%rovx(p%i_2(q), p%j_2(q)) = rovx_avg

        rovy_avg = 0.5 * (gs(p%n_1)%rovy(p%i_1(q), p%j_1(q)) + gs(p%n_2)%rovy(p%i_2(q), p%j_2(q)))
        gs(p%n_1)%rovy(p%i_1(q), p%j_1(q)) = rovy_avg
        gs(p%n_2)%rovy(p%i_2(q), p%j_2(q)) = rovy_avg

    end do

    ! do the above code using array indexing
    ! this should be more efficient

!   allocate(ro_avgs(p%nk))
!   allocate(roe_avgs(p%nk))
!   allocate(rovx_avgs(p%nk))
!   allocate(rovy_avgs(p%nk))
!
!   ro_avgs(:) = 0.5 * (gs(p%n_1)%ro(p%i_1(:), p%j_1(:)) + gs(p%n_2)%ro(p%i_2(:), p%j_2(:)))
!   gs(p%n_1)%ro(p%i_1, p%j_1) = ro_avgs
!   gs(p%n_2)%ro(p%i_2, p%j_2) = ro_avgs
!
!   roe_avgs(:) = 0.5 * (gs(p%n_1)%roe(p%i_1(:), p%j_1(:)) + gs(p%n_2)%roe(p%i_2(:), p%j_2(:)))
!   gs(p%n_1)%roe(p%i_1, p%j_1) = roe_avgs
!   gs(p%n_2)%roe(p%i_2, p%j_2) = roe_avgs
!
!   rovx_avgs = 0.5 * (gs(p%n_1)%rovx(p%i_1, p%j_1) + gs(p%n_2)%rovx(p%i_2, p%j_2))
!   gs(p%n_1)%rovx(p%i_1, p%j_1) = rovx_avgs
!   gs(p%n_2)%rovx(p%i_2, p%j_2) = rovx_avgs
!
!   rovy_avgs = 0.5 * (gs(p%n_1)%rovy(p%i_1, p%j_1) + gs(p%n_2)%rovy(p%i_2, p%j_2))
!   gs(p%n_1)%rovy(p%i_1, p%j_1) = rovy_avgs
!   gs(p%n_2)%rovy(p%i_2, p%j_2) = rovy_avgs


end subroutine apply_patch

end module patches