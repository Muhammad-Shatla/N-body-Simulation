CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      PROGRAM LEAP
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      real*4 RERR
      PARAMETER (N= 2048)
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
      DIMENSION Al (N, 3), As (N, 3), PENS (N), PENL(N)
      DIMENSION XS (N, 3), VS (N, 3), XL (N, 3), VL (N, 3)
C
C    
      T = 0.0
      PMASS = 1D0/ DFLOAT(N)
      CALL ICPLUM (PMASS, XS, VS)
      DO I = 1, N
        DO K = 1,3
        XL (I,K) = XS (I, K)
        VL (I,K) = VS (I, K)        
        END DO
      END DO
      CALL ACC (PMASS, XS, As, PE, PENS)
      CALL ACC (PMASS, XL, Al, PE, PENL)
      TEND = 200D0
      DTLARGE = 1.*1D-3
      NLEAP = 2
c
      COUNTER = 0
      DO WHILE (T .LT. TEND)
C
      DT = DTLARGE
      DT2 = DT * DT
      call stepper (DT, DT2, PMASS, XL, VL, PEL, TEL, PENL,al)
C
      DT = DTLARGE/ DFLOAT (NLEAP)
      DT2 = DT * DT 
      DO KL = 1, NLEAP
       call stepper (DT, DT2, PMASS, XS, VS, PES, TES, PENS, as)
      END DO
C      
      T = T + DTLARGE
C      
      KOUNTER = KOUNTER + 1
      IF (MOD(KOUNTER, 10) .EQ. 0) THEN 
        CALL ANALYSER (RERR, XS, VS, XL, VL, PENS, PENL, pmass)
      SSV =0.
      SSC = 0.
      do iter =  1 , n
        do k = 1, 3
        SSC = ssc + xL(Iter,K)**2
        SSV = ssv + VL(Itrer,K)**2
        end do 
      end do
c
c        WRITE(*,*) T, TEL/PEL, tel+Pel, rerr
       WRITE(*,*) T,  sqrt(ssc), sqrt(ssv), tel/pel !rerr, tel+pel, tel/pel, 
      END IF
C     6 PENS (7) + 0.5*PMASS * (VS(7,1)**2+VS(7,2)**2+VS(7,3)**2),
C     6 PENS(7), 0.5*PMASS * (VS(7,1)**2+VS(7,2)**2+VS(7,3)**2)  
C
      END DO
C
C
      END
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE STEPPER (DT, DT2, PMASS, X, V, PE, TE, PEN, a)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (N = 2048)
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
      DIMENSION X (N, 3), V (N, 3), A (N, 3), AOLD (N, 3), PEN (N)
C
      TE = 0D0
C 
      DO II = 1, N
       DO KK = 1,3
        X (II, KK) = X(II, KK) + DT * V(II, KK) + 0.5D0 * DT2* A(II, KK)
c        TE = TE + 0.5 * PMASS * V (II, KK)* V (II, KK)
        AOLD(II, KK) = A (II, KK)
c        CALL ACC (PMASS, X, A, PE, PEN)
        END DO
      END DO
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx       CALL ACC (PMASS, X, A, PE, PEN)
c
      DO II = 1, N
       DO KK = 1,3
        A(II, KK) = 0D0
       END DO
      END DO
      PE = 0D0
      DO NN = 1, N
         PEN(NN) = 0D0
      END DO
C
      DO II = 1,  N - 1
        DO JJ = 1, n ! II +1, N
        dfac= 1.0
        if (jj .lt. ii + 1) dfac =0.0 
      DX = X (II, 1) - X (JJ, 1)
      DY = X (II, 2) - X (JJ, 2)
      DZ = X (II, 3) - X (JJ, 3)   
C
      DIST2 = DX* DX + DY* DY + DZ * DZ
      DIST = DSQRT (DIST2+0.01)
      COFFI = 1d0/(DIST * DIST * DIST)
      DISTI = 1D0/DIST
      PE = PE - PMASS * PMASS * DISTI * dfac
C
      A(JJ, 1) = A(JJ, 1) + (PMASS*COFFI)* DX *dfac
      A(JJ, 2) = A(JJ, 2) + (PMASS*COFFI)* DY *dfac
      A(JJ, 3) = A(JJ, 3) + (PMASS*COFFI)* DZ *dfac
C
      A(II, 1) = A(II, 1) - (PMASS*COFFI)* DX * dfac
      A(II, 2) = A(II, 2) - (PMASS*COFFI)* DY * dfac
      A(II, 3) = A(II, 3) - (PMASS*COFFI)* DZ * dfac
C
      PEN (II) = PEN (II) - PMASS * PMASS * DISTI * dfac
      PEN (JJ) = PEN (JJ) - PMASS * PMASS * DISTI * dfac
c
        END DO
      END DO

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      DO II = 1, N
       DO KK = 1,3
        V(II, KK) = V(II, KK)+ 0.5D0*DT* (AOLD(II, KK) + A(II, KK))        
        TE = TE + 0.5 * PMASS * V (II, KK)* V (II, KK)      
       END DO
      END DO
C      
c     WRITE(*,*) TE/PE, TE+PE

      RETURN
      END
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE ACC (PMASS, X, A, PE, PEN)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (N= 2048)
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
      DIMENSION X (N, 3), A (N, 3), PEN (N)
C
      DO II = 1, N
       DO KK = 1,3
        A(II, KK) = 0D0
       END DO
      END DO
      PE = 0D0
      DO NN = 1, N
         PEN(NN) = 0D0
      END DO
C
      DO II = 1, N -1
        DO JJ = II+1, N
      DX = X (II, 1) - X (JJ, 1)
      DY = X (II, 2) - X (JJ, 2)
      DZ = X (II, 3) - X (JJ, 3)   
C
      DIST2 = DX* DX + DY* DY + DZ * DZ
      DIST = DSQRT (DIST2+0.00001)
      COFFI = 1d0/(DIST * DIST * DIST)
      DISTI = 1D0/DIST
      PE = PE - PMASS * PMASS * DISTI
C
      A(JJ, 1) = A(JJ, 1) + (PMASS*COFFI)* DX
      A(JJ, 2) = A(JJ, 2) + (PMASS*COFFI)* DY
      A(JJ, 3) = A(JJ, 3) + (PMASS*COFFI)* DZ
C
      A(II, 1) = A(II, 1) - (PMASS*COFFI)* DX
      A(II, 2) = A(II, 2) - (PMASS*COFFI)* DY
      A(II, 3) = A(II, 3) - (PMASS*COFFI)* DZ
C
      PEN (II) = PEN (II) - PMASS * PMASS * DISTI
      PEN (JJ) = PEN (JJ) - PMASS * PMASS * DISTI
C
        END DO
      END DO
C
      RETURN 
      END
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE ANALYSER (RERR_P, XS, VS, XL, VL, PENS, PENL, pmass)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (N = 2048)
      real*4 RERR_P
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      DIMENSION XS (N, 3), VS(N, 3), XL(N, 3), VL (N, 3),
     6  PENS(N),PENL(N)     
C
      ERR_C = 0D0  
      ENORM_C =0d0 
      DO I = 1, N
        DO K = 1, 3
        ERR_C = ERR_C + 
     6 (XS(I, K) - XL (I, K)) * (XS(I, K) - XL (I, K))
        ENORM_C = ENORM_C + XS(I, K) * XS(I, K) 
        END DO
      END DO
      RERR_C = DSQRT(ERR_C / ENORM_C)
C
      ERR_V = 0D0  
      ENORM_V = 0D0
      DO I = 1, N
        DO K = 1, 3
        ERR_V = ERR_V + 
     6 (VS(I, K) - VL (I, K)) * (VS(I, K) - VL (I, K))
        ENORM_V = ENORM_V + VS(I, K) * VS(I, K) 
        END DO
      END DO
      RERR_V = sqrt(ERR_V / ENORM_V)
C
      ERR_P = 0D0  
      ENORM_P = 0D0
      DO I = 1, N
        PXS = PMASS*(XS (I, 2) * VS (I, 3) - XS(I, 3) * VS(I, 2)) 
        PYS = PMASS*(XS (I, 3) * VS (I, 1) - XS(I, 1) * VS(I, 3))        
        PZS = PMASS*(XS (I, 1) * VS (I, 2) - XS(I, 2)*  VS(I, 1)) 
        PXL = PMASS*(XL (I, 2) * VL (I, 3) - XL(I, 3) * VL(I, 2)) 
        PYL = PMASS*(XL (I, 3) * VL (I, 1) - XL(I, 1) * VL(I, 3))        
        PZL = PMASS*(XL (I, 1) * VL (I, 2) - XL(I, 2)*  VL(I, 1)) 
C       
       
        ERR_P = ERR_P + 
     6  (PXS - PXL)*(PXS - PXL)+
     6  (PYS - PYL)*(PYS - PYL)+(PZS - PZL) *  (PZS - PZL)  
        ENORM_P = ENORM_P + PXS*PXS+ PYS*PYS + PZS*PZS
      END DO
      RERR_P = sqrt(ERR_P / ENORM_P)
c     

      ERR_E = 0D0  
      ENORM_E = 0D0
      DO I = 1, N
        ES = PENS(I) + 0.5D0 * PMASS * 
     6 (VS(I, 1)*VS(I,1) + VS(I, 2)*VS(I,2) + VS(I, 3)*VS(I,3))
        EL = PENL(I) + 0.5D0 * PMASS * 
     6 (VL(I, 1)*VL(I,1) + VL(I, 2)*VL(I,2) + VL(I, 3)*VL(I,3))

        ERR_E = ERR_E + (ES - EL)* (ES - EL)
        ENORM_E = ENORM_E + ES * ES 
      END DO
      RERR_E = SQRT(ERR_E / ENORM_E)
C
C
      RETURN
      END
C
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE ICPLUM (PMASS, X, V)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      real*4 ran2
      PARAMETER (N = 2048)
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C Generate initial conditions from Plummer model (A & A 37, 183). FROM AARSETH CODE data.f
C
      DIMENSION A(8), X (N, 3), V(N , 3), CMR(3), CMRDOT(3)

C 
      DO  K = 1,3
              CMR(K) = 0D0
              CMRDOT(K) = 0D0
      END DO
      ZMASS = 1d0  ! Total mass...
      kdum = 40000
C
      TWOPI = 8.0*ATAN(1.0D0)
      DO 85 I = 1,N
   82     A(1) = RAN2(KDUM)
          RI = (A(1)**(-0.6666667) - 1.0)**(-0.5)
*       Reject distant particles.
          IF (RI.GT.10.0) GO TO 82
          A(2) = RAN2(KDUM)
          A(3) = RAN2(KDUM)
          X(I,3) = (1.0 - 2.0*A(2))*RI
          X(I,1) = SQRT(RI**2 - X(I, 3)**2)*COS(TWOPI*A(3))
          X(I,2) = SQRT(RI**2 - X(I, 3)**2)*SIN(TWOPI*A(3))
   83     A(4) = RAN2(KDUM)
          A(5) = RAN2(KDUM)
          A(6) = A(4)**2*(1.0 - A(4)**2)**3.5
          IF (0.1*A(5).GT.A(6)) GO TO 83
          A(8) = A(4)*SQRT(2.0)/(1.0 + RI**2)**0.25
          A(6) = RAN2(KDUM)
          A(7) = RAN2(KDUM)
          V(I,3) = (1.0 - 2.0*A(6))*A(8)
          V(I,1) = SQRT(A(8)**2 - V(I, 3)**2)*COS(TWOPI*A(7))
          V(I,2) = SQRT(A(8)**2 - V(I, 3)**2)*SIN(TWOPI*A(7))
*
*       Accumulate the centre of mass terms.
          DO 84 K = 1,3
              CMR(K) = CMR(K) + PMASS* X(I,K)
              CMRDOT(K) = CMRDOT(K) + PMASS *V(I,K)
   84     CONTINUE
   85 CONTINUE
*
*       Scale coordinates & velocities to analytical expectation values.
      SX = 1.5*TWOPI/16.0
      SV = SQRT(ZMASS/SX)
      DO 88 I = 1,N
          DO 86 K = 1,3
              X(I,K) = X(I,K) - CMR(K)/ZMASS
              V(I,K) = V(I,K) - CMRDOT(K)/ZMASS
              X(I,K) = SX*X(I,K)
              V(I,K) = SV*V(I,K)
   86     CONTINUE
   88 CONTINUE
*
   90 RETURN
*
      END
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      FUNCTION RAN2(IDUM)
*
*
*       Random number generator (Press p. 195).
*       ---------------------------------------
*
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1./M)
      COMMON/RAND2/  IY,IFF,IR(97) 
*     DATA  IFF /0/
*
*
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
          IFF = 1
          IDUM = MOD(IC-IDUM,M)
          DO 11 J = 1,97
              IDUM = MOD(IA*IDUM+IC,M)
              IR(J) = IDUM
   11     CONTINUE
          IDUM = MOD(IA*IDUM+IC,M)
          IY = IDUM
      END IF
      J = 1 + (97*IY)/M
      IF (J.GT.97.OR.J.LT.1) WRITE (6,12)  J, IDUM
   12 FORMAT (/,'  TROUBLES IN RAN2   J IDUM ',2I12)
      IY = IR(J)
      RAN2 = IY*RM
      IDUM = MOD(IA*IDUM+IC,M)
      IR(J) = IDUM
*
      RETURN
*
      END 

















