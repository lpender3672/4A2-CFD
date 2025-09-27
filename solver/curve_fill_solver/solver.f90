
subroutine curve_fill_solver() bind(C, name="curve_fill_solver")
  USE SFC_REORDER
      IMPLICIT NONE
      INTEGER :: NC, B
      REAL*8, ALLOCATABLE :: XC(:), YC(:), ZC(:)
      INTEGER*8, ALLOCATABLE :: KEYS(:)
      INTEGER*4, ALLOCATABLE :: P(:)
      REAL*8 :: BB(6)

      NC = 1000
      B  = 19
      ALLOCATE(XC(NC),YC(NC),ZC(NC),KEYS(NC),P(NC))

      ! fill XC,YC,ZC with cell centroids
      ! compute bounding box
      BB(1)=MINVAL(XC); BB(2)=MAXVAL(XC)
      BB(3)=MINVAL(YC); BB(4)=MAXVAL(YC)
      BB(5)=MINVAL(ZC); BB(6)=MAXVAL(ZC)

      CALL BUILD_KEYS_3D(XC,YC,ZC,BB,B,KEYS,NC)
      CALL ARGSORT_KEYS(KEYS,P,NC)

      ! now reorder XC,YC,ZC, U(:,:), etc using permutation P
end subroutine curve_fill_solver

