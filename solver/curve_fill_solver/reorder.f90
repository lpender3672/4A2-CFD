      MODULE SFC_REORDER
      IMPLICIT NONE
      INTEGER, PARAMETER :: I4 = 4
      INTEGER, PARAMETER :: I8 = 8
      INTEGER, PARAMETER :: RK = 8
      CONTAINS

      INTEGER*4 FUNCTION QUANTIZE(V, VMIN, VMAX, B)
      REAL*8, INTENT(IN) :: V, VMIN, VMAX
      INTEGER,  INTENT(IN) :: B
      REAL*8 :: T
      IF (VMAX .LE. VMIN) THEN
         QUANTIZE = 0
         RETURN
      END IF
      T = (V - VMIN) / (VMAX - VMIN)
      IF (T .LT. 0.0D0) T = 0.0D0
      IF (T .GT. 1.0D0) T = 1.0D0
      QUANTIZE = INT( T * DBLE(2**B - 1) )
      END FUNCTION QUANTIZE

      INTEGER*8 FUNCTION MORTON3(IX,IY,IZ,B)
      INTEGER*4, INTENT(IN) :: IX,IY,IZ
      INTEGER,   INTENT(IN) :: B
      INTEGER :: I
      INTEGER*8 :: K
      MORTON3 = 0_8
      DO I = 0, B-1
         IF (BTEST(IX,I)) THEN
            K = 1_8
         ELSE
            K = 0_8
         END IF
         MORTON3 = IOR(MORTON3, ISHFT(K, 3*I+2))
         IF (BTEST(IY,I)) THEN
            K = 1_8
         ELSE
            K = 0_8
         END IF
         MORTON3 = IOR(MORTON3, ISHFT(K, 3*I+1))
         IF (BTEST(IZ,I)) THEN
            K = 1_8
         ELSE
            K = 0_8
         END IF
         MORTON3 = IOR(MORTON3, ISHFT(K, 3*I))
      END DO
      END FUNCTION MORTON3

      SUBROUTINE BUILD_KEYS_3D(XC,YC,ZC,BB,B,KEYS,NC)
      REAL*8, INTENT(IN) :: XC(NC), YC(NC), ZC(NC)
      REAL*8, INTENT(IN) :: BB(6)
      INTEGER, INTENT(IN) :: B, NC
      INTEGER*8, INTENT(OUT) :: KEYS(NC)
      INTEGER :: C
      INTEGER*4 :: IX,IY,IZ
      DO C = 1, NC
         IX = QUANTIZE(XC(C), BB(1), BB(2), B)
         IY = QUANTIZE(YC(C), BB(3), BB(4), B)
         IZ = QUANTIZE(ZC(C), BB(5), BB(6), B)
         KEYS(C) = MORTON3(IX,IY,IZ,B)
      END DO
      END SUBROUTINE BUILD_KEYS_3D

      SUBROUTINE ARGSORT_KEYS(KEYS,P,N)
      INTEGER, INTENT(IN) :: N
      INTEGER*8, INTENT(IN)  :: KEYS(N)
      INTEGER*4, INTENT(OUT) :: P(N)
      INTEGER :: I
      DO I=1,N
         P(I)=I
      END DO
      CALL MERGESORT(KEYS,P,1,N)
      END SUBROUTINE ARGSORT_KEYS

      RECURSIVE SUBROUTINE MERGESORT(K,IDX,L,R)
      INTEGER*8, INTENT(IN) :: K(:)
      INTEGER*4, INTENT(INOUT) :: IDX(:)
      INTEGER, INTENT(IN) :: L,R
      INTEGER :: M
      IF (L .GE. R) RETURN
      M = (L+R)/2
      CALL MERGESORT(K,IDX,L,M)
      CALL MERGESORT(K,IDX,M+1,R)
      CALL MERGE(K,IDX,L,M,R)
      END SUBROUTINE MERGESORT

      SUBROUTINE MERGE(K,IDX,L,M,R)
      INTEGER*8, INTENT(IN) :: K(:)
      INTEGER*4, INTENT(INOUT) :: IDX(:)
      INTEGER, INTENT(IN) :: L,M,R
      INTEGER :: N1,N2,I,J,POS,T
      INTEGER*4, ALLOCATABLE :: LIDX(:), RIDX(:)
      N1 = M-L+1
      N2 = R-M
      ALLOCATE(LIDX(N1),RIDX(N2))
      LIDX = IDX(L:M)
      RIDX = IDX(M+1:R)
      I=1; J=1; POS=L
      DO WHILE (I<=N1 .AND. J<=N2)
         IF (K(LIDX(I)) .LE. K(RIDX(J))) THEN
            IDX(POS)=LIDX(I); I=I+1
         ELSE
            IDX(POS)=RIDX(J); J=J+1
         END IF
         POS=POS+1
      END DO
      DO T=I,N1
         IDX(POS)=LIDX(T); POS=POS+1
      END DO
      DO T=J,N2
         IDX(POS)=RIDX(T); POS=POS+1
      END DO
      DEALLOCATE(LIDX,RIDX)
      END SUBROUTINE MERGE

      END MODULE SFC_REORDER
