!       COMPUTATION OF SPECIAL FUNCTIONS
!
!          Shanjie Zhang and Jianming Jin
!
!       Copyrighted but permission granted to use code in programs.
!       Buy their book "Computation of Special Functions", 1996, John Wiley & Sons, Inc.
!
!       Scipy changes:
!       - Compiled into a single source file and changed REAL To DBLE throughout.
!       - Changed according to ERRATA.
!       - Changed GAMMA to GAMMA2 and PSI to PSI_SPEC to avoid potential conflicts.
!       - Made functions return sf_error codes in ISFER variables instead
!         of printing warnings. The codes are
!         - SF_ERROR_OK        = 0: no error
!         - SF_ERROR_SINGULAR  = 1: singularity encountered
!         - SF_ERROR_UNDERFLOW = 2: floating point underflow
!         - SF_ERROR_OVERFLOW  = 3: floating point overflow
!         - SF_ERROR_SLOW      = 4: too many iterations required
!         - SF_ERROR_LOSS      = 5: loss of precision
!         - SF_ERROR_NO_RESULT = 6: no result obtained
!         - SF_ERROR_DOMAIN    = 7: out of domain
!         - SF_ERROR_ARG       = 8: invalid input parameter
!         - SF_ERROR_OTHER     = 9: unclassified error
!
      FUNCTION DNAN()
      IMPLICIT NONE
!*--DNAN28
      DOUBLE PRECISION DNAN
      DNAN = 0.0D0
      DNAN = 0.0D0/DNAN
      END
 
      FUNCTION DINF()
      IMPLICIT NONE
!*--DINF37
      DOUBLE PRECISION DINF
      DINF = 1.0D300
      DINF = DINF*DINF
      END
 
      SUBROUTINE CPDSA(N,Z,Cdn)
!
!       ===========================================================
!       Purpose: Compute complex parabolic cylinder function Dn(z)
!                for small argument
!       Input:   z   --- complex argument of D(z)
!                n   --- Order of D(z) (n = 0,-1,-2,...)
!       Output:  CDN --- Dn(z)
!       Routine called: GAIH for computing Г(x), x=n/2 (n=1,2,...)
!       ===========================================================
!
      IMPLICIT NONE
!*--CPDSA56
      COMPLEX*16 ca0 , cb0 , Cdn , cdw , cr , Z
      DOUBLE PRECISION eps , g0 , g1 , ga0 , gm , pd , pi , sq2 , va0 , &
                     & vm , vt , xn
      INTEGER m , N
      eps = 1.0D-15
      pi = 3.141592653589793D0
      sq2 = DSQRT(2.0D0)
      ca0 = EXP(-.25D0*Z*Z)
      va0 = 0.5D0*(1.0D0-N)
      IF ( N==0.0 ) THEN
         Cdn = ca0
      ELSEIF ( ABS(Z)==0.0 ) THEN
         IF ( va0<=0.0 .AND. va0==INT(va0) ) THEN
            Cdn = 0.0D0
         ELSE
            CALL GAIH(va0,ga0)
            pd = DSQRT(pi)/(2.0D0**(-.5D0*N)*ga0)
            Cdn = DCMPLX(pd,0.0D0)
         ENDIF
      ELSE
         xn = -N
         CALL GAIH(xn,g1)
         cb0 = 2.0D0**(-0.5D0*N-1.0D0)*ca0/g1
         vt = -.5D0*N
         CALL GAIH(vt,g0)
         Cdn = DCMPLX(g0,0.0D0)
         cr = (1.0D0,0.0D0)
         DO m = 1 , 250
            vm = .5D0*(m-N)
            CALL GAIH(vm,gm)
            cr = -cr*sq2*Z/m
            cdw = gm*cr
            Cdn = Cdn + cdw
            IF ( ABS(cdw)<ABS(Cdn)*eps ) GOTO 50
         ENDDO
 50      Cdn = cb0*Cdn
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE CFS(Z,Zf,Zd)
!
!       =========================================================
!       Purpose: Compute complex Fresnel Integral S(z) and S'(z)
!       Input :  z  --- Argument of S(z)
!       Output:  ZF --- S(z)
!                ZD --- S'(z)
!       =========================================================
!
      IMPLICIT NONE
!*--CFS113
      COMPLEX*16 cf , cf0 , cf1 , cg , cr , d , s , Z , z0 , Zd , Zf ,  &
               & zp , zp2
      DOUBLE PRECISION eps , pi , w0 , wb , wb0
      INTEGER k , m
      eps = 1.0D-14
      pi = 3.141592653589793D0
      w0 = ABS(Z)
      zp = 0.5D0*pi*Z*Z
      zp2 = zp*zp
      z0 = (0.0D0,0.0D0)
      IF ( Z==z0 ) THEN
         s = z0
      ELSEIF ( w0<=2.5 ) THEN
         s = Z*zp/3.0D0
         cr = s
         wb0 = 0.0D0
         DO k = 1 , 80
            cr = -.5D0*cr*(4.0D0*k-1.0D0)/k/(2.0D0*k+1.0D0)             &
               & /(4.0D0*k+3.0D0)*zp2
            s = s + cr
            wb = ABS(s)
            IF ( DABS(wb-wb0)<eps .AND. k>10 ) GOTO 100
            wb0 = wb
         ENDDO
      ELSEIF ( w0>2.5 .AND. w0<4.5 ) THEN
         m = 85
         s = z0
         cf1 = z0
         cf0 = (1.0D-100,0.0D0)
         DO k = m , 0 , -1
            cf = (2.0D0*k+3.0D0)*cf0/zp - cf1
            IF ( k/=INT(k/2)*2 ) s = s + cf
            cf1 = cf0
            cf0 = cf
         ENDDO
         s = 2.0D0/(pi*Z)*SIN(zp)/cf*s
      ELSE
!          Auxiliary functions f(z) and g(z) can be computed using an
!          asymptotic expansion in the right quadrant |arg(z)| <= pi/4, not pi/2
!          as sometimes suggested. Use the symmetry S(z) = -iS(-iz).
!          Interestingly, most of the expansion code is the same across
!          the quadrants. (The forth power in Z is the equalizer here.)
!          Only one constant has to be adapted.
         IF ( DIMAG(Z)>-DBLE(Z) .AND. DIMAG(Z)<=DBLE(Z) ) THEN
!            right quadrant
            d = DCMPLX(.5D0,0.0D0)
         ELSEIF ( DIMAG(Z)>DBLE(Z) .AND. DIMAG(Z)>=-DBLE(Z) ) THEN
!            upper quadrant
            d = DCMPLX(0.0D0,-.5D0)
         ELSEIF ( DIMAG(Z)<-DBLE(Z) .AND. DIMAG(Z)>=DBLE(Z) ) THEN
!            left quadrant
            d = DCMPLX(-.5D0,0.0D0)
         ELSE
!            lower quadrant
            d = DCMPLX(0.0D0,.5D0)
         ENDIF
         cr = (1.0D0,0.0D0)
         cf = (1.0D0,0.0D0)
         DO k = 1 , 20
            cr = -.25D0*cr*(4.0D0*k-1.0D0)*(4.0D0*k-3.0D0)/zp2
            cf = cf + cr
         ENDDO
         cr = (1.0D0,0.0D0)
         cg = (1.0D0,0.0D0)
         DO k = 1 , 12
            cr = -.25D0*cr*(4.0D0*k+1.0D0)*(4.0D0*k-1.0D0)/zp2
            cg = cg + cr
         ENDDO
         cg = cg/(pi*Z*Z)
         s = d - (cf*COS(zp)+cg*SIN(zp))/(pi*Z)
      ENDIF
 100  Zf = s
      Zd = SIN(0.5*pi*Z*Z)
      END
 
!       **********************************
 
      SUBROUTINE LQMN(Mm,M,N,X,Qm,Qd)
!
!       ==========================================================
!       Purpose: Compute the associated Legendre functions of the
!                second kind, Qmn(x) and Qmn'(x)
!       Input :  x  --- Argument of Qmn(x)
!                m  --- Order of Qmn(x)  ( m = 0,1,2,… )
!                n  --- Degree of Qmn(x) ( n = 0,1,2,… )
!                mm --- Physical dimension of QM and QD
!       Output:  QM(m,n) --- Qmn(x)
!                QD(m,n) --- Qmn'(x)
!       ==========================================================
!
      IMPLICIT NONE
!*--LQMN208
      INTEGER i , j , k , km , ls , M , Mm , N
      DOUBLE PRECISION q0 , q1 , q10 , Qd , qf , qf0 , qf1 , qf2 , Qm , &
                     & X , xq , xs
      DIMENSION Qm(0:Mm,0:N) , Qd(0:Mm,0:N)
      IF ( DABS(X)==1.0D0 ) THEN
         DO i = 0 , M
            DO j = 0 , N
               Qm(i,j) = 1.0D+300
               Qd(i,j) = 1.0D+300
            ENDDO
         ENDDO
         RETURN
      ENDIF
      ls = 1
      IF ( DABS(X)>1.0D0 ) ls = -1
      xs = ls*(1.0D0-X*X)
      xq = DSQRT(xs)
      q0 = 0.5D0*DLOG(DABS((X+1.0D0)/(X-1.0D0)))
      IF ( DABS(X)<1.0001D0 ) THEN
         Qm(0,0) = q0
         Qm(0,1) = X*q0 - 1.0D0
         Qm(1,0) = -1.0D0/xq
         Qm(1,1) = -ls*xq*(q0+X/(1.0D0-X*X))
         DO i = 0 , 1
            DO j = 2 , N
               Qm(i,j) = ((2.0D0*j-1.0D0)*X*Qm(i,j-1)-(j+i-1.0D0)       &
                       & *Qm(i,j-2))/(j-i)
            ENDDO
         ENDDO
         DO j = 0 , N
            DO i = 2 , M
               Qm(i,j) = -2.0D0*(i-1.0D0)*X/xq*Qm(i-1,j)                &
                       & - ls*(j+i-1.0D0)*(j-i+2.0D0)*Qm(i-2,j)
            ENDDO
         ENDDO
      ELSE
         IF ( DABS(X)>1.1D0 ) THEN
            km = 40 + M + N
         ELSE
            km = (40+M+N)*INT(-1.0-1.8*LOG(X-1.0))
         ENDIF
         qf2 = 0.0D0
         qf1 = 1.0D0
         qf0 = 0.0D0
         DO k = km , 0 , -1
            qf0 = ((2*k+3.0D0)*X*qf1-(k+2.0D0)*qf2)/(k+1.0D0)
            IF ( k<=N ) Qm(0,k) = qf0
            qf2 = qf1
            qf1 = qf0
         ENDDO
         DO k = 0 , N
            Qm(0,k) = q0*Qm(0,k)/qf0
         ENDDO
         qf2 = 0.0D0
         qf1 = 1.0D0
         DO k = km , 0 , -1
            qf0 = ((2*k+3.0D0)*X*qf1-(k+1.0D0)*qf2)/(k+2.0D0)
            IF ( k<=N ) Qm(1,k) = qf0
            qf2 = qf1
            qf1 = qf0
         ENDDO
         q10 = -1.0D0/xq
         DO k = 0 , N
            Qm(1,k) = q10*Qm(1,k)/qf0
         ENDDO
         DO j = 0 , N
            q0 = Qm(0,j)
            q1 = Qm(1,j)
            DO i = 0 , M - 2
               qf = -2.0D0*(i+1)*X/xq*q1 + (j-i)*(j+i+1.0D0)*q0
               Qm(i+2,j) = qf
               q0 = q1
               q1 = qf
            ENDDO
         ENDDO
      ENDIF
      Qd(0,0) = ls/xs
      DO j = 1 , N
         Qd(0,j) = ls*j*(Qm(0,j-1)-X*Qm(0,j))/xs
      ENDDO
      DO j = 0 , N
         DO i = 1 , M
            Qd(i,j) = ls*i*X/xs*Qm(i,j) + (i+j)*(j-i+1.0D0)/xq*Qm(i-1,j)
         ENDDO
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE CLPMN(Mm,M,N,X,Y,Ntype,Cpm,Cpd)
!
!       =========================================================
!       Purpose: Compute the associated Legendre functions Pmn(z)
!                and their derivatives Pmn'(z) for a complex
!                argument
!       Input :  x     --- Real part of z
!                y     --- Imaginary part of z
!                m     --- Order of Pmn(z),  m = 0,1,2,...,n
!                n     --- Degree of Pmn(z), n = 0,1,2,...,N
!                mm    --- Physical dimension of CPM and CPD
!                ntype --- type of cut, either 2 or 3
!       Output:  CPM(m,n) --- Pmn(z)
!                CPD(m,n) --- Pmn'(z)
!       =========================================================
!
      IMPLICIT NONE
!*--CLPMN318
      COMPLEX*16 Cpd , Cpm , z , zq , zs
      DOUBLE PRECISION DINF , X , Y
      INTEGER i , j , ls , M , Mm , N , Ntype
      DIMENSION Cpm(0:Mm,0:N) , Cpd(0:Mm,0:N)
      z = DCMPLX(X,Y)
      DO i = 0 , N
         DO j = 0 , M
            Cpm(j,i) = (0.0D0,0.0D0)
            Cpd(j,i) = (0.0D0,0.0D0)
         ENDDO
      ENDDO
      Cpm(0,0) = (1.0D0,0.0D0)
      IF ( N==0 ) RETURN
      IF ( DABS(X)==1.0D0 .AND. Y==0.0D0 ) THEN
         DO i = 1 , N
            Cpm(0,i) = X**i
            Cpd(0,i) = 0.5D0*i*(i+1)*X**(i+1)
         ENDDO
         DO j = 1 , N
            DO i = 1 , M
               IF ( i==1 ) THEN
                  Cpd(i,j) = DINF()
               ELSEIF ( i==2 ) THEN
                  Cpd(i,j) = -0.25D0*(j+2)*(j+1)*j*(j-1)*X**(j+1)
               ENDIF
            ENDDO
         ENDDO
         RETURN
      ENDIF
      IF ( Ntype==2 ) THEN
!       sqrt(1 - z^2) with branch cut on |x|>1
         zs = (1.0D0-z*z)
         zq = -SQRT(zs)
         ls = -1
      ELSE
!       sqrt(z^2 - 1) with branch cut between [-1, 1]
         zs = (z*z-1.0D0)
         zq = SQRT(zs)
         IF ( X<0D0 ) zq = -zq
         ls = 1
      ENDIF
      DO i = 1 , M
!       DLMF 14.7.15
         Cpm(i,i) = (2.0D0*i-1.0D0)*zq*Cpm(i-1,i-1)
      ENDDO
      DO i = 0 , MIN(M,N-1)
!       DLMF 14.10.7
         Cpm(i,i+1) = (2.0D0*i+1.0D0)*z*Cpm(i,i)
      ENDDO
      DO i = 0 , M
         DO j = i + 2 , N
!       DLMF 14.10.3
            Cpm(i,j) = ((2.0D0*j-1.0D0)*z*Cpm(i,j-1)-(i+j-1.0D0)        &
                     & *Cpm(i,j-2))/(j-i)
         ENDDO
      ENDDO
      Cpd(0,0) = (0.0D0,0.0D0)
      DO j = 1 , N
!       DLMF 14.10.5
         Cpd(0,j) = ls*j*(z*Cpm(0,j)-Cpm(0,j-1))/zs
      ENDDO
      DO i = 1 , M
         DO j = i , N
!       derivative of DLMF 14.7.11 & DLMF 14.10.6 for type 3
!       derivative of DLMF 14.7.8 & DLMF 14.10.1 for type 2
            Cpd(i,j) = ls*(-i*z*Cpm(i,j)/zs+(j+i)*(j-i+1.0D0)           &
                     & /zq*Cpm(i-1,j))
         ENDDO
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE VVSA(Va,X,Pv)
!
!       ===================================================
!       Purpose: Compute parabolic cylinder function Vv(x)
!                for small argument
!       Input:   x  --- Argument
!                va --- Order
!       Output:  PV --- Vv(x)
!       Routine called : GAMMA2 for computing Г(x)
!       ===================================================
!
      IMPLICIT NONE
!*--VVSA407
      DOUBLE PRECISION a0 , ep , eps , fac , g1 , ga0 , gm , gw , pi ,  &
                     & Pv , r , r1 , sq2 , sv , sv0 , v1 , Va , va0 ,   &
                     & vb0 , vm
      DOUBLE PRECISION X
      INTEGER m
      eps = 1.0D-15
      pi = 3.141592653589793D0
      ep = EXP(-.25D0*X*X)
      va0 = 1.0D0 + 0.5D0*Va
      IF ( X/=0.0 ) THEN
         sq2 = DSQRT(2.0D0)
         a0 = 2.0D0**(-.5D0*Va)*ep/(2.0D0*pi)
         sv = DSIN(-(Va+.5D0)*pi)
         v1 = -.5D0*Va
         CALL GAMMA2(v1,g1)
         Pv = (sv+1.0D0)*g1
         r = 1.0D0
         fac = 1.0D0
         DO m = 1 , 250
            vm = .5D0*(m-Va)
            CALL GAMMA2(vm,gm)
            r = r*sq2*X/m
            fac = -fac
            gw = fac*sv + 1.0D0
            r1 = gw*r*gm
            Pv = Pv + r1
            IF ( DABS(r1/Pv)<eps .AND. gw/=0.0 ) GOTO 50
         ENDDO
 50      Pv = a0*Pv
      ELSEIF ( va0<=0.0 .AND. va0==INT(va0) .OR. Va==0.0 ) THEN
         Pv = 0.0D0
      ELSE
         vb0 = -0.5D0*Va
         sv0 = DSIN(va0*pi)
         CALL GAMMA2(va0,ga0)
         Pv = 2.0D0**vb0*sv0/ga0
      ENDIF
      END
 
 
 
!       **********************************
!       SciPy: Changed P from a character array to an integer array.
      SUBROUTINE JDZO(Nt,N,M,P,Zo)
!
!       ===========================================================
!       Purpose: Compute the zeros of Bessel functions Jn(x) and
!                Jn'(x), and arrange them in the order of their
!                magnitudes
!       Input :  NT    --- Number of total zeros ( NT ≤ 1200 )
!       Output:  ZO(L) --- Value of the L-th zero of Jn(x)
!                          and Jn'(x)
!                N(L)  --- n, order of Jn(x) or Jn'(x) associated
!                          with the L-th zero
!                M(L)  --- m, serial number of the zeros of Jn(x)
!                          or Jn'(x) associated with the L-th zero
!                          ( L is the serial number of all the
!                            zeros of Jn(x) and Jn'(x) )
!                P(L)  --- 0 (TM) or 1 (TE), a code for designating the
!                          zeros of Jn(x)  or Jn'(x).
!                          In the waveguide applications, the zeros
!                          of Jn(x) correspond to TM modes and
!                          those of Jn'(x) correspond to TE modes
!       Routine called:    BJNDD for computing Jn(x), Jn'(x) and
!                          Jn''(x)
!       =============================================================
!
      IMPLICIT NONE
!*--JDZO479
      DOUBLE PRECISION bj , dj , fj , x , x0 , x1 , x2 , xm , Zo , zoc
      INTEGER i , j , k , l , l0 , l1 , l2 , M , m1 , mm , N , n1 , nm ,&
            & Nt
      INTEGER P(1400) , p1(70)
      DIMENSION N(1400) , M(1400) , Zo(0:1400) , n1(70) , m1(70) ,      &
              & zoc(0:70) , bj(101) , dj(101) , fj(101)
      x = 0
      zoc(0) = 0
      IF ( Nt<600 ) THEN
         xm = -1.0 + 2.248485*Nt**0.5 - .0159382*Nt +                   &
            & 3.208775E-4*Nt**1.5
         nm = INT(14.5+.05875*Nt)
         mm = INT(.02*Nt) + 6
      ELSE
         xm = 5.0 + 1.445389*Nt**.5 + .01889876*Nt - 2.147763E-4*Nt**1.5
         nm = INT(27.8+.0327*Nt)
         mm = INT(.01088*Nt) + 10
      ENDIF
      l0 = 0
      DO i = 1 , nm
         x1 = .407658 + .4795504*(i-1)**.5 + .983618*(i-1)
         x2 = 1.99535 + .8333883*(i-1)**.5 + .984584*(i-1)
         l1 = 0
         DO j = 1 , mm
            IF ( i/=1 .OR. j/=1 ) THEN
               x = x1
 10            CALL BJNDD(i,x,bj,dj,fj)
               x0 = x
               x = x - dj(i)/fj(i)
               IF ( x1>xm ) GOTO 20
               IF ( DABS(x-x0)>1.0D-10 ) GOTO 10
            ENDIF
            l1 = l1 + 1
            n1(l1) = i - 1
            m1(l1) = j
            IF ( i==1 ) m1(l1) = j - 1
            p1(l1) = 1
            zoc(l1) = x
            IF ( i<=15 ) THEN
               x1 = x + 3.057 + .0122*(i-1) + (1.555+.41575*(i-1))/(j+1)&
                  & **2
            ELSE
               x1 = x + 2.918 + .01924*(i-1) + (6.26+.13205*(i-1))/(j+1)&
                  & **2
            ENDIF
 20         x = x2
 40         CALL BJNDD(i,x,bj,dj,fj)
            x0 = x
            x = x - bj(i)/dj(i)
            IF ( x<=xm ) THEN
               IF ( DABS(x-x0)>1.0D-10 ) GOTO 40
               l1 = l1 + 1
               n1(l1) = i - 1
               m1(l1) = j
               p1(l1) = 0
               zoc(l1) = x
               IF ( i<=15 ) THEN
                  x2 = x + 3.11 + .0138*(i-1) + (.04832+.2804*(i-1))    &
                     & /(j+1)**2
               ELSE
                  x2 = x + 3.001 + .0105*(i-1) + (11.52+.48525*(i-1))   &
                     & /(j+3)**2
               ENDIF
            ENDIF
         ENDDO
         l = l0 + l1
         l2 = l
 50      IF ( l0==0 ) THEN
            DO k = 1 , l
               Zo(k) = zoc(k)
               N(k) = n1(k)
               M(k) = m1(k)
               P(k) = p1(k)
            ENDDO
            l1 = 0
         ELSEIF ( l0/=0 ) THEN
            IF ( Zo(l0)>=zoc(l1) ) THEN
               Zo(l0+l1) = Zo(l0)
               N(l0+l1) = N(l0)
               M(l0+l1) = M(l0)
               P(l0+l1) = P(l0)
               l0 = l0 - 1
            ELSE
               Zo(l0+l1) = zoc(l1)
               N(l0+l1) = n1(l1)
               M(l0+l1) = m1(l1)
               P(l0+l1) = p1(l1)
               l1 = l1 - 1
            ENDIF
         ENDIF
         IF ( l1/=0 ) GOTO 50
         l0 = l2
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE CBK(M,N,C,Cv,Qt,Ck,Bk)
!
!       =====================================================
!       Purpose: Compute coefficient Bk's for oblate radial
!                functions with a small argument
!       =====================================================
!
      IMPLICIT NONE
!*--CBK590
      DOUBLE PRECISION Bk , C , Ck , Cv , eps , Qt , r1 , s1 , sw , t , &
                     & u , v , w
      INTEGER i , i1 , ip , j , k , M , N , n2 , nm
      DIMENSION Bk(200) , Ck(200) , u(200) , v(200) , w(200)
      eps = 1.0D-14
      ip = 1
      IF ( N-M==2*INT((N-M)/2) ) ip = 0
      nm = 25 + INT(0.5*(N-M)+C)
      u(1) = 0.0D0
      n2 = nm - 2
      DO j = 2 , n2
         u(j) = C*C
      ENDDO
      DO j = 1 , n2
         v(j) = (2.0*j-1.0-ip)*(2.0*(j-M)-ip) + M*(M-1.0) - Cv
      ENDDO
      DO j = 1 , nm - 1
         w(j) = (2.0*j-ip)*(2.0*j+1.0-ip)
      ENDDO
      IF ( ip==0 ) THEN
         sw = 0.0D0
         DO k = 0 , n2 - 1
            s1 = 0.0D0
            i1 = k - M + 1
            DO i = i1 , nm
               IF ( i>=0 ) THEN
                  r1 = 1.0D0
                  DO j = 1 , k
                     r1 = r1*(i+M-j)/j
                  ENDDO
                  s1 = s1 + Ck(i+1)*(2.0*i+M)*r1
                  IF ( DABS(s1-sw)<DABS(s1)*eps ) GOTO 20
                  sw = s1
               ENDIF
            ENDDO
 20         Bk(k+1) = Qt*s1
         ENDDO
      ELSEIF ( ip==1 ) THEN
         sw = 0.0D0
         DO k = 0 , n2 - 1
            s1 = 0.0D0
            i1 = k - M + 1
            DO i = i1 , nm
               IF ( i>=0 ) THEN
                  r1 = 1.0D0
                  DO j = 1 , k
                     r1 = r1*(i+M-j)/j
                  ENDDO
                  IF ( i>0 ) s1 = s1 + Ck(i)*(2.0*i+M-1)*r1
                  s1 = s1 - Ck(i+1)*(2.0*i+M)*r1
                  IF ( DABS(s1-sw)<DABS(s1)*eps ) GOTO 40
                  sw = s1
               ENDIF
            ENDDO
 40         Bk(k+1) = Qt*s1
         ENDDO
      ENDIF
      w(1) = w(1)/v(1)
      Bk(1) = Bk(1)/v(1)
      DO k = 2 , n2
         t = v(k) - w(k-1)*u(k)
         w(k) = w(k)/t
         Bk(k) = (Bk(k)-Bk(k-1)*u(k))/t
      ENDDO
      DO k = n2 - 1 , 1 , -1
         Bk(k) = Bk(k) - w(k)*Bk(k+1)
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE RMN2SP(M,N,C,X,Cv,Df,Kd,R2f,R2d)
!
!       ======================================================
!       Purpose: Compute prolate spheroidal radial function
!                of the second kind with a small argument
!       Routines called:
!            (1) LPMNS for computing the associated Legendre
!                functions of the first kind
!            (2) LQMNS for computing the associated Legendre
!                functions of the second kind
!            (3) KMN for computing expansion coefficients
!                and joining factors
!       ======================================================
!
      IMPLICIT NONE
!*--RMN2SP682
      DOUBLE PRECISION C , ck1 , ck2 , Cv , Df , dn , eps , ga , gb ,   &
                     & gc , pd , pm , qd , qm , r1 , r2 , R2d , R2f ,   &
                     & r3 , r4
      DOUBLE PRECISION sd , sd0 , sd1 , sd2 , sdm , sf , spd1 , spd2 ,  &
                     & spl , su0 , su1 , su2 , sum , sw , X
      INTEGER ip , j , j1 , j2 , k , Kd , ki , l1 , M , N , nm , nm1 ,  &
            & nm2 , nm3
      DIMENSION pm(0:251) , pd(0:251) , qm(0:251) , qd(0:251) , dn(200) &
              & , Df(200)
      IF ( DABS(Df(1))<1.0D-280 ) THEN
         R2f = 1.0D+300
         R2d = 1.0D+300
         RETURN
      ENDIF
      eps = 1.0D-14
      ip = 1
      nm1 = INT((N-M)/2)
      IF ( N-M==2*nm1 ) ip = 0
      nm = 25 + nm1 + INT(C)
      nm2 = 2*nm + M
      CALL KMN(M,N,C,Cv,Kd,Df,dn,ck1,ck2)
      CALL LPMNS(M,nm2,X,pm,pd)
      CALL LQMNS(M,nm2,X,qm,qd)
      su0 = 0.0D0
      sw = 0.0D0
      DO k = 1 , nm
         j = 2*k - 2 + M + ip
         su0 = su0 + Df(k)*qm(j)
         IF ( k>nm1 .AND. DABS(su0-sw)<DABS(su0)*eps ) GOTO 100
         sw = su0
      ENDDO
 100  sd0 = 0.0D0
      DO k = 1 , nm
         j = 2*k - 2 + M + ip
         sd0 = sd0 + Df(k)*qd(j)
         IF ( k>nm1 .AND. DABS(sd0-sw)<DABS(sd0)*eps ) GOTO 200
         sw = sd0
      ENDDO
 200  su1 = 0.0D0
      sd1 = 0.0D0
      DO k = 1 , M
         j = M - 2*k + ip
         IF ( j<0 ) j = -j - 1
         su1 = su1 + dn(k)*qm(j)
         sd1 = sd1 + dn(k)*qd(j)
      ENDDO
      ga = ((X-1.0D0)/(X+1.0D0))**(0.5D0*M)
      DO k = 1 , M
         j = M - 2*k + ip
         IF ( j<0 ) THEN
            IF ( j<0 ) j = -j - 1
            r1 = 1.0D0
            DO j1 = 1 , j
               r1 = (M+j1)*r1
            ENDDO
            r2 = 1.0D0
            DO j2 = 1 , M - j - 2
               r2 = j2*r2
            ENDDO
            r3 = 1.0D0
            sf = 1.0D0
            DO l1 = 1 , j
               r3 = 0.5D0*r3*(-j+l1-1.0)*(j+l1)/((M+l1)*l1)*(1.0-X)
               sf = sf + r3
            ENDDO
            IF ( M-j>=2 ) gb = (M-j-1.0D0)*r2
            IF ( M-j<=1 ) gb = 1.0D0
            spl = r1*ga*gb*sf
            su1 = su1 + (-1)**(j+M)*dn(k)*spl
            spd1 = M/(X*X-1.0D0)*spl
            gc = 0.5D0*j*(j+1.0)/(M+1.0)
            sd = 1.0D0
            r4 = 1.0D0
            DO l1 = 1 , j - 1
               r4 = 0.5D0*r4*(-j+l1)*(j+l1+1.0)/((M+l1+1.0)*l1)*(1.0-X)
               sd = sd + r4
            ENDDO
            spd2 = r1*ga*gb*gc*sd
            sd1 = sd1 + (-1)**(j+M)*dn(k)*(spd1+spd2)
         ENDIF
      ENDDO
      su2 = 0.0D0
      ki = (2*M+1+ip)/2
      nm3 = nm + ki
      DO k = ki , nm3
         j = 2*k - 1 - M - ip
         su2 = su2 + dn(k)*pm(j)
         IF ( j>M .AND. DABS(su2-sw)<DABS(su2)*eps ) GOTO 300
         sw = su2
      ENDDO
 300  sd2 = 0.0D0
      DO k = ki , nm3
         j = 2*k - 1 - M - ip
         sd2 = sd2 + dn(k)*pd(j)
         IF ( j>M .AND. DABS(sd2-sw)<DABS(sd2)*eps ) GOTO 400
         sw = sd2
      ENDDO
 400  sum = su0 + su1 + su2
      sdm = sd0 + sd1 + sd2
      R2f = sum/ck2
      R2d = sdm/ck2
      END
 
 
 
!       **********************************
 
      SUBROUTINE BERNOB(N,Bn)
!
!       ======================================
!       Purpose: Compute Bernoulli number Bn
!       Input :  n --- Serial number
!       Output:  BN(n) --- Bn
!       ======================================
!
      IMPLICIT NONE
!*--BERNOB802
      DOUBLE PRECISION Bn , r1 , r2 , s , tpi
      INTEGER k , m , N
      DIMENSION Bn(0:N)
      tpi = 6.283185307179586D0
      Bn(0) = 1.0D0
      Bn(1) = -0.5D0
      Bn(2) = 1.0D0/6.0D0
      r1 = (2.0D0/tpi)**2
      DO m = 4 , N , 2
         r1 = -r1*(m-1)*m/(tpi*tpi)
         r2 = 1.0D0
         DO k = 2 , 10000
            s = (1.0D0/k)**m
            r2 = r2 + s
            IF ( s<1.0D-15 ) GOTO 50
         ENDDO
 50      Bn(m) = r1*r2
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE BERNOA(N,Bn)
!
!       ======================================
!       Purpose: Compute Bernoulli number Bn
!       Input :  n --- Serial number
!       Output:  BN(n) --- Bn
!       ======================================
!
      IMPLICIT NONE
!*--BERNOA837
      DOUBLE PRECISION Bn , r , s
      INTEGER j , k , m , N
      DIMENSION Bn(0:N)
      Bn(0) = 1.0D0
      Bn(1) = -0.5D0
      DO m = 2 , N
         s = -(1.0D0/(m+1.0D0)-0.5D0)
         DO k = 2 , m - 1
            r = 1.0D0
            DO j = 2 , k
               r = r*(j+m-k)/j
            ENDDO
            s = s - r*Bn(k)
         ENDDO
         Bn(m) = s
      ENDDO
      DO m = 3 , N , 2
         Bn(m) = 0.0D0
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE QSTAR(M,N,C,Ck,Ck1,Qs,Qt)
!
!       =========================================================
!       Purpose: Compute Q*mn(-ic) for oblate radial functions
!                with a small argument
!       =========================================================
!
      IMPLICIT NONE
!*--QSTAR872
      DOUBLE PRECISION ap , C , Ck , Ck1 , Qs , qs0 , Qt , r , s , sk
      INTEGER i , ip , k , l , M , N
      DIMENSION ap(200) , Ck(200)
      ip = 1
      IF ( N-M==2*INT((N-M)/2) ) ip = 0
      r = 1.0D0/Ck(1)**2
      ap(1) = r
      DO i = 1 , M
         s = 0.0D0
         DO l = 1 , i
            sk = 0.0D0
            DO k = 0 , l
               sk = sk + Ck(k+1)*Ck(l-k+1)
            ENDDO
            s = s + sk*ap(i-l+1)
         ENDDO
         ap(i+1) = -r*s
      ENDDO
      qs0 = ap(M+1)
      DO l = 1 , M
         r = 1.0D0
         DO k = 1 , l
            r = r*(2.0D0*k+ip)*(2.0D0*k-1.0D0+ip)/(2.0D0*k)**2
         ENDDO
         qs0 = qs0 + ap(M-l+1)*r
      ENDDO
      Qs = (-1)**ip*Ck1*(Ck1*qs0)/C
      Qt = -2.0D0/Ck1*Qs
      END
 
 
 
!       **********************************
 
      SUBROUTINE CV0(Kd,M,Q,A0)
!
!       =====================================================
!       Purpose: Compute the initial characteristic value of
!                Mathieu functions for m ≤ 12  or q ≤ 300 or
!                q ≥ m*m
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!       Output:  A0 --- Characteristic value
!       Routines called:
!             (1) CVQM for computing initial characteristic
!                 value for q ≤ 3*m
!             (2) CVQL for computing initial characteristic
!                 value for q ≥ m*m
!       ====================================================
!
      IMPLICIT NONE
!*--CV0927
      DOUBLE PRECISION A0 , Q , q2
      INTEGER Kd , M
      q2 = Q*Q
      IF ( M==0 ) THEN
         IF ( Q<=1.0 ) THEN
            A0 = (((.0036392*q2-.0125868)*q2+.0546875)*q2-.5)*q2
         ELSEIF ( Q<=10.0 ) THEN
            A0 = ((3.999267D-3*Q-9.638957D-2)*Q-.88297)*Q + .5542818
         ELSE
            CALL CVQL(Kd,M,Q,A0)
         ENDIF
      ELSEIF ( M==1 ) THEN
         IF ( Q<=1.0 .AND. Kd==2 ) THEN
            A0 = (((-6.51E-4*Q-.015625)*Q-.125)*Q+1.0)*Q + 1.0
         ELSEIF ( Q<=1.0 .AND. Kd==3 ) THEN
            A0 = (((-6.51E-4*Q+.015625)*Q-.125)*Q-1.0)*Q + 1.0
         ELSEIF ( Q<=10.0 .AND. Kd==2 ) THEN
            A0 = (((-4.94603D-4*Q+1.92917D-2)*Q-.3089229)*Q+1.33372)    &
               & *Q + .811752
         ELSEIF ( Q<=10.0 .AND. Kd==3 ) THEN
            A0 = ((1.971096D-3*Q-5.482465D-2)*Q-1.152218)*Q + 1.10427
         ELSE
            CALL CVQL(Kd,M,Q,A0)
         ENDIF
      ELSEIF ( M==2 ) THEN
         IF ( Q<=1.0 .AND. Kd==1 ) THEN
            A0 = (((-.0036391*q2+.0125888)*q2-.0551939)*q2+.416667)     &
               & *q2 + 4.0
         ELSEIF ( Q<=1.0 .AND. Kd==4 ) THEN
            A0 = (.0003617*q2-.0833333)*q2 + 4.0
         ELSEIF ( Q<=15 .AND. Kd==1 ) THEN
            A0 = (((3.200972D-4*Q-8.667445D-3)*Q-1.829032D-4)           &
               & *Q+.9919999)*Q + 3.3290504
         ELSEIF ( Q<=10.0 .AND. Kd==4 ) THEN
            A0 = ((2.38446D-3*Q-.08725329)*Q-4.732542D-3)*Q + 4.00909
         ELSE
            CALL CVQL(Kd,M,Q,A0)
         ENDIF
      ELSEIF ( M==3 ) THEN
         IF ( Q<=1.0 .AND. Kd==2 ) THEN
            A0 = ((6.348E-4*Q+.015625)*Q+.0625)*q2 + 9.0
         ELSEIF ( Q<=1.0 .AND. Kd==3 ) THEN
            A0 = ((6.348E-4*Q-.015625)*Q+.0625)*q2 + 9.0
         ELSEIF ( Q<=20.0 .AND. Kd==2 ) THEN
            A0 = (((3.035731D-4*Q-1.453021D-2)*Q+.19069602)*Q-.1039356) &
               & *Q + 8.9449274
         ELSEIF ( Q<=15.0 .AND. Kd==3 ) THEN
            A0 = ((9.369364D-5*Q-.03569325)*Q+.2689874)*Q + 8.771735
         ELSE
            CALL CVQL(Kd,M,Q,A0)
         ENDIF
      ELSEIF ( M==4 ) THEN
         IF ( Q<=1.0 .AND. Kd==1 ) THEN
            A0 = ((-2.1E-6*q2+5.012E-4)*q2+.0333333)*q2 + 16.0
         ELSEIF ( Q<=1.0 .AND. Kd==4 ) THEN
            A0 = ((3.7E-6*q2-3.669E-4)*q2+.0333333)*q2 + 16.0
         ELSEIF ( Q<=25.0 .AND. Kd==1 ) THEN
            A0 = (((1.076676D-4*Q-7.9684875D-3)*Q+.17344854)*Q-.5924058)&
               & *Q + 16.620847
         ELSEIF ( Q<=20.0 .AND. Kd==4 ) THEN
            A0 = ((-7.08719D-4*Q+3.8216144D-3)*Q+.1907493)*Q + 15.744
         ELSE
            CALL CVQL(Kd,M,Q,A0)
         ENDIF
      ELSEIF ( M==5 ) THEN
         IF ( Q<=1.0 .AND. Kd==2 ) THEN
            A0 = ((6.8E-6*Q+1.42E-5)*q2+.0208333)*q2 + 25.0
         ELSEIF ( Q<=1.0 .AND. Kd==3 ) THEN
            A0 = ((-6.8E-6*Q+1.42E-5)*q2+.0208333)*q2 + 25.0
         ELSEIF ( Q<=35.0 .AND. Kd==2 ) THEN
            A0 = (((2.238231D-5*Q-2.983416D-3)*Q+.10706975)*Q-.600205)  &
               & *Q + 25.93515
         ELSEIF ( Q<=25.0 .AND. Kd==3 ) THEN
            A0 = ((-7.425364D-4*Q+2.18225D-2)*Q+4.16399D-2)*Q + 24.897
         ELSE
            CALL CVQL(Kd,M,Q,A0)
         ENDIF
      ELSEIF ( M==6 ) THEN
         IF ( Q<=1.0 ) THEN
            A0 = (.4D-6*q2+.0142857)*q2 + 36.0
         ELSEIF ( Q<=40.0 .AND. Kd==1 ) THEN
            A0 = (((-1.66846D-5*Q+4.80263D-4)*Q+2.53998D-2)*Q-.181233)  &
               & *Q + 36.423
         ELSEIF ( Q<=35.0 .AND. Kd==4 ) THEN
            A0 = ((-4.57146D-4*Q+2.16609D-2)*Q-2.349616D-2)*Q + 35.99251
         ELSE
            CALL CVQL(Kd,M,Q,A0)
         ENDIF
      ELSEIF ( M==7 ) THEN
         IF ( Q<=10.0 ) THEN
            CALL CVQM(M,Q,A0)
         ELSEIF ( Q<=50.0 .AND. Kd==2 ) THEN
            A0 = (((-1.411114D-5*Q+9.730514D-4)*Q-3.097887D-3)          &
               & *Q+3.533597D-2)*Q + 49.0547
         ELSEIF ( Q<=40.0 .AND. Kd==3 ) THEN
            A0 = ((-3.043872D-4*Q+2.05511D-2)*Q-9.16292D-2)*Q + 49.19035
         ELSE
            CALL CVQL(Kd,M,Q,A0)
         ENDIF
      ELSEIF ( M>=8 ) THEN
         IF ( Q<=3.*M ) THEN
            CALL CVQM(M,Q,A0)
         ELSEIF ( Q>M*M ) THEN
            CALL CVQL(Kd,M,Q,A0)
         ELSEIF ( M==8 .AND. Kd==1 ) THEN
            A0 = (((8.634308D-6*Q-2.100289D-3)*Q+.169072)*Q-4.64336)    &
               & *Q + 109.4211
         ELSEIF ( M==8 .AND. Kd==4 ) THEN
            A0 = ((-6.7842D-5*Q+2.2057D-3)*Q+.48296)*Q + 56.59
         ELSEIF ( M==9 .AND. Kd==2 ) THEN
            A0 = (((2.906435D-6*Q-1.019893D-3)*Q+.1101965)*Q-3.821851)  &
               & *Q + 127.6098
         ELSEIF ( M==9 .AND. Kd==3 ) THEN
            A0 = ((-9.577289D-5*Q+.01043839)*Q+.06588934)*Q + 78.0198
         ELSEIF ( M==10 .AND. Kd==1 ) THEN
            A0 = (((5.44927D-7*Q-3.926119D-4)*Q+.0612099)*Q-2.600805)   &
               & *Q + 138.1923
         ELSEIF ( M==10 .AND. Kd==4 ) THEN
            A0 = ((-7.660143D-5*Q+.01132506)*Q-.09746023)*Q + 99.29494
         ELSEIF ( M==11 .AND. Kd==2 ) THEN
            A0 = (((-5.67615D-7*Q+7.152722D-6)*Q+.01920291)*Q-1.081583) &
               & *Q + 140.88
         ELSEIF ( M==11 .AND. Kd==3 ) THEN
            A0 = ((-6.310551D-5*Q+.0119247)*Q-.2681195)*Q + 123.667
         ELSEIF ( M==12 .AND. Kd==1 ) THEN
            A0 = (((-2.38351D-7*Q-2.90139D-5)*Q+.02023088)*Q-1.289)     &
               & *Q + 171.2723
         ELSEIF ( M==12 .AND. Kd==4 ) THEN
            A0 = (((3.08902D-7*Q-1.577869D-4)*Q+.0247911)*Q-1.05454)    &
               & *Q + 161.471
         ENDIF
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE CVQM(M,Q,A0)
!
!       =====================================================
!       Purpose: Compute the characteristic value of Mathieu
!                functions for q ≤ m*m
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!       Output:  A0 --- Initial characteristic value
!       =====================================================
!
      IMPLICIT NONE
!*--CVQM1080
      DOUBLE PRECISION A0 , hm1 , hm3 , hm5 , Q
      INTEGER M
      hm1 = .5*Q/(M*M-1.0)
      hm3 = .25*hm1**3/(M*M-4.0)
      hm5 = hm1*hm3*Q/((M*M-1.0)*(M*M-9.0))
      A0 = M*M + Q*(hm1+(5.0*M*M+7.0)*hm3+(9.0*M**4+58.0*M*M+29.0)*hm5)
      END
 
!       **********************************
 
      SUBROUTINE CVQL(Kd,M,Q,A0)
!
!       ========================================================
!       Purpose: Compute the characteristic value of Mathieu
!                functions  for q ≥ 3m
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!       Output:  A0 --- Initial characteristic value
!       ========================================================
!
      IMPLICIT NONE
!*--CVQL1105
      DOUBLE PRECISION A0 , c1 , cv1 , cv2 , d1 , d2 , d3 , d4 , p1 ,   &
                     & p2 , Q , w , w2 , w3 , w4 , w6
      INTEGER Kd , M
      w = 0.0D0
      IF ( Kd==1 .OR. Kd==2 ) w = 2.0D0*M + 1.0D0
      IF ( Kd==3 .OR. Kd==4 ) w = 2.0D0*M - 1.0D0
      w2 = w*w
      w3 = w*w2
      w4 = w2*w2
      w6 = w2*w4
      d1 = 5.0 + 34.0/w2 + 9.0/w4
      d2 = (33.0+410.0/w2+405.0/w4)/w
      d3 = (63.0+1260.0/w2+2943.0/w4+486.0/w6)/w2
      d4 = (527.0+15617.0/w2+69001.0/w4+41607.0/w6)/w3
      c1 = 128.0
      p2 = Q/w4
      p1 = DSQRT(p2)
      cv1 = -2.0*Q + 2.0*w*DSQRT(Q) - (w2+1.0)/8.0
      cv2 = (w+3.0/w) + d1/(32.0*p1) + d2/(8.0*c1*p2)
      cv2 = cv2 + d3/(64.0*c1*p1*p2) + d4/(16.0*c1*c1*p2*p2)
      A0 = cv1 - cv2/(c1*p1)
      END
 
 
 
      INTEGER FUNCTION MSTA1(X,Mp)
!
!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that the magnitude of
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point
!       ===================================================
!
      IMPLICIT NONE
!*--MSTA11146
      DOUBLE PRECISION a0 , ENVJ , f , f0 , f1 , X
      INTEGER it , Mp , n0 , n1 , nn
      a0 = DABS(X)
      n0 = INT(1.1D0*a0) + 1
      f0 = ENVJ(n0,a0) - Mp
      n1 = n0 + 5
      f1 = ENVJ(n1,a0) - Mp
      DO it = 1 , 20
         nn = n1 - (n1-n0)/(1.0D0-f0/f1)
         f = ENVJ(nn,a0) - Mp
         IF ( ABS(nn-n1)<1 ) GOTO 100
         n0 = n1
         f0 = f1
         n1 = nn
         f1 = f
      ENDDO
 100  MSTA1 = nn
      END
 
 
      INTEGER FUNCTION MSTA2(X,N,Mp)
!
!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!       ===================================================
!
      IMPLICIT NONE
!*--MSTA21183
      DOUBLE PRECISION a0 , ejn , ENVJ , f , f0 , f1 , hmp , obj , X
      INTEGER it , Mp , N , n0 , n1 , nn
      a0 = DABS(X)
      hmp = 0.5D0*Mp
      ejn = ENVJ(N,a0)
      IF ( ejn<=hmp ) THEN
         obj = Mp
         n0 = INT(1.1*a0) + 1
      ELSE
         obj = hmp + ejn
         n0 = N
      ENDIF
      f0 = ENVJ(n0,a0) - obj
      n1 = n0 + 5
      f1 = ENVJ(n1,a0) - obj
      DO it = 1 , 20
         nn = n1 - (n1-n0)/(1.0D0-f0/f1)
         f = ENVJ(nn,a0) - obj
         IF ( ABS(nn-n1)<1 ) GOTO 100
         n0 = n1
         f0 = f1
         n1 = nn
         f1 = f
      ENDDO
 100  MSTA2 = nn + 10
      END
 
      REAL*8 FUNCTION ENVJ(N,X)
      IMPLICIT NONE
!*--ENVJ1216
      INTEGER N
      DOUBLE PRECISION X
      ENVJ = 0.5D0*DLOG10(6.28D0*N) - N*DLOG10(1.36D0*X/N)
      END
 
!       **********************************
 
      SUBROUTINE ITTJYB(X,Ttj,Tty)
!
!       ==========================================================
!       Purpose: Integrate [1-J0(t)]/t with respect to t from 0
!                to x, and Y0(t)/t with respect to t from x to ∞
!       Input :  x   --- Variable in the limits  ( x ≥ 0 )
!       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
!                TTY --- Integration of Y0(t)/t from x to ∞
!       ==========================================================
!
      IMPLICIT NONE
!*--ITTJYB1238
      DOUBLE PRECISION e0 , el , f0 , g0 , pi , t , t1 , Ttj , Tty , X ,&
                     & x1 , xt
      pi = 3.141592653589793D0
      el = .5772156649015329D0
      IF ( X==0.0D0 ) THEN
         Ttj = 0.0D0
         Tty = -1.0D+300
      ELSEIF ( X<=4.0D0 ) THEN
         x1 = X/4.0D0
         t = x1*x1
         Ttj = ((((((.35817D-4*t-.639765D-3)*t+.7092535D-2)*t-          &
             & .055544803D0)*t+.296292677D0)*t-.999999326D0)            &
             & *t+1.999999936D0)*t
         Tty = (((((((-.3546D-5*t+.76217D-4)*t-.1059499D-2)*t+          &
             & .010787555D0)*t-.07810271D0)*t+.377255736D0)             &
             & *t-1.114084491D0)*t+1.909859297D0)*t
         e0 = el + DLOG(X/2.0D0)
         Tty = pi/6.0D0 + e0/pi*(2.0D0*Ttj-e0) - Tty
      ELSEIF ( X<=8.0D0 ) THEN
         xt = X + .25D0*pi
         t1 = 4.0D0/X
         t = t1*t1
         f0 = (((((.0145369D0*t-.0666297D0)*t+.1341551D0)*t-.1647797D0) &
            & *t+.1608874D0)*t-.2021547D0)*t + .7977506D0
         g0 = ((((((.0160672D0*t-.0759339D0)*t+.1576116D0)*t-.1960154D0)&
            & *t+.1797457D0)*t-.1702778D0)*t+.3235819D0)*t1
         Ttj = (f0*DCOS(xt)+g0*DSIN(xt))/(DSQRT(X)*X)
         Ttj = Ttj + el + DLOG(X/2.0D0)
         Tty = (f0*DSIN(xt)-g0*DCOS(xt))/(DSQRT(X)*X)
      ELSE
         t = 8.0D0/X
         xt = X + .25D0*pi
         f0 = (((((.18118D-2*t-.91909D-2)*t+.017033D0)*t-.9394D-3)      &
            & *t-.051445D0)*t-.11D-5)*t + .7978846D0
         g0 = (((((-.23731D-2*t+.59842D-2)*t+.24437D-2)*t-.0233178D0)   &
            & *t+.595D-4)*t+.1620695D0)*t
         Ttj = (f0*DCOS(xt)+g0*DSIN(xt))/(DSQRT(X)*X)                   &
             & + el + DLOG(X/2.0D0)
         Tty = (f0*DSIN(xt)-g0*DCOS(xt))/(DSQRT(X)*X)
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE ITTJYA(X,Ttj,Tty)
!
!       =========================================================
!       Purpose: Integrate [1-J0(t)]/t with respect to t from 0
!                to x, and Y0(t)/t with respect to t from x to ∞
!       Input :  x   --- Variable in the limits  ( x ≥ 0 )
!       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
!                TTY --- Integration of Y0(t)/t from x to ∞
!       =========================================================
!
      IMPLICIT NONE
!*--ITTJYA1297
      DOUBLE PRECISION a0 , b1 , bj0 , bj1 , by0 , by1 , e0 , el , g0 , &
                     & g1 , pi , px , qx , r , r0 , r1 , r2 , rs , t ,  &
                     & Ttj
      DOUBLE PRECISION Tty , vt , X , xk
      INTEGER k , l
      pi = 3.141592653589793D0
      el = .5772156649015329D0
      IF ( X==0.0D0 ) THEN
         Ttj = 0.0D0
         Tty = -1.0D+300
      ELSEIF ( X<=20.0D0 ) THEN
         Ttj = 1.0D0
         r = 1.0D0
         DO k = 2 , 100
            r = -.25D0*r*(k-1.0D0)/(k*k*k)*X*X
            Ttj = Ttj + r
            IF ( DABS(r)<DABS(Ttj)*1.0D-12 ) GOTO 50
         ENDDO
 50      Ttj = Ttj*.125D0*X*X
         e0 = .5D0*(pi*pi/6.0D0-el*el) - (.5D0*DLOG(X/2.0D0)+el)        &
            & *DLOG(X/2.0D0)
         b1 = el + DLOG(X/2.0D0) - 1.5D0
         rs = 1.0D0
         r = -1.0D0
         DO k = 2 , 100
            r = -.25D0*r*(k-1.0D0)/(k*k*k)*X*X
            rs = rs + 1.0D0/k
            r2 = r*(rs+1.0D0/(2.0D0*k)-(el+DLOG(X/2.0D0)))
            b1 = b1 + r2
            IF ( DABS(r2)<DABS(b1)*1.0D-12 ) GOTO 100
         ENDDO
 100     Tty = 2.0D0/pi*(e0+.125D0*X*X*b1)
      ELSE
         a0 = DSQRT(2.0D0/(pi*X))
         bj0 = 0.0D0
         by0 = 0.0D0
         bj1 = 0.0D0
         DO l = 0 , 1
            vt = 4.0D0*l*l
            px = 1.0D0
            r = 1.0D0
            DO k = 1 , 14
               r = -.0078125D0*r*(vt-(4.0D0*k-3.0D0)**2)/(X*k)          &
                 & *(vt-(4.0D0*k-1.0D0)**2)/((2.0D0*k-1.0D0)*X)
               px = px + r
               IF ( DABS(r)<DABS(px)*1.0D-12 ) GOTO 120
            ENDDO
 120        qx = 1.0D0
            r = 1.0D0
            DO k = 1 , 14
               r = -.0078125D0*r*(vt-(4.0D0*k-1.0D0)**2)/(X*k)          &
                 & *(vt-(4.0D0*k+1.0D0)**2)/(2.0D0*k+1.0D0)/X
               qx = qx + r
               IF ( DABS(r)<DABS(qx)*1.0D-12 ) GOTO 140
            ENDDO
 140        qx = .125D0*(vt-1.0D0)/X*qx
            xk = X - (.25D0+.5D0*l)*pi
            bj1 = a0*(px*DCOS(xk)-qx*DSIN(xk))
            by1 = a0*(px*DSIN(xk)+qx*DCOS(xk))
            IF ( l==0 ) THEN
               bj0 = bj1
               by0 = by1
            ENDIF
         ENDDO
         t = 2.0D0/X
         g0 = 1.0D0
         r0 = 1.0D0
         DO k = 1 , 10
            r0 = -k*k*t*t*r0
            g0 = g0 + r0
         ENDDO
         g1 = 1.0D0
         r1 = 1.0D0
         DO k = 1 , 10
            r1 = -k*(k+1.0D0)*t*t*r1
            g1 = g1 + r1
         ENDDO
         Ttj = 2.0D0*g1*bj0/(X*X) - g0*bj1/X + el + DLOG(X/2.0D0)
         Tty = 2.0D0*g1*by0/(X*X) - g0*by1/X
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE CJYLV(V,Z,Cbjv,Cdjv,Cbyv,Cdyv)
!
!       ===================================================
!       Purpose: Compute Bessel functions Jv(z) and Yv(z)
!                and their derivatives with a complex
!                argument and a large order
!       Input:   v --- Order of Jv(z) and Yv(z)
!                z --- Complex argument
!       Output:  CBJV --- Jv(z)
!                CDJV --- Jv'(z)
!                CBYV --- Yv(z)
!                CDYV --- Yv'(z)
!       Routine called:
!                CJK to compute the expansion coefficients
!       ===================================================
!
      IMPLICIT NONE
!*--CJYLV1402
      DOUBLE PRECISION a , pi , V , v0 , vr
      COMPLEX*16 Cbjv , Cbyv , Cdjv , Cdyv , ceta , cf , cfj , cfy ,    &
               & csj , csy , ct , ct2 , cws , Z
      INTEGER i , k , km , l , l0 , lf
      DIMENSION cf(12) , a(91)
      km = 12
      CALL CJK(km,a)
      pi = 3.141592653589793D0
      DO l = 1 , 0 , -1
         v0 = V - l
         cws = SQRT(1.0D0-(Z/v0)*(Z/v0))
         ceta = cws + LOG(Z/v0/(1.0D0+cws))
         ct = 1.0D0/cws
         ct2 = ct*ct
         DO k = 1 , km
            l0 = k*(k+1)/2 + 1
            lf = l0 + k
            cf(k) = a(lf)
            DO i = lf - 1 , l0 , -1
               cf(k) = cf(k)*ct2 + a(i)
            ENDDO
            cf(k) = cf(k)*ct**k
         ENDDO
         vr = 1.0D0/v0
         csj = (1.0D0,0.0D0)
         DO k = 1 , km
            csj = csj + cf(k)*vr**k
         ENDDO
         Cbjv = SQRT(ct/(2.0D0*pi*v0))*EXP(v0*ceta)*csj
         IF ( l==1 ) cfj = Cbjv
         csy = (1.0D0,0.0D0)
         DO k = 1 , km
            csy = csy + (-1)**k*cf(k)*vr**k
         ENDDO
         Cbyv = -SQRT(2.0D0*ct/(pi*v0))*EXP(-v0*ceta)*csy
         IF ( l==1 ) cfy = Cbyv
      ENDDO
      Cdjv = -V/Z*Cbjv + cfj
      Cdyv = -V/Z*Cbyv + cfy
      END
 
 
 
!       **********************************
 
      SUBROUTINE RMN2L(M,N,C,X,Df,Kd,R2f,R2d,Id)
!
!       ========================================================
!       Purpose: Compute prolate and oblate spheroidal radial
!                functions of the second kind for given m, n,
!                c and a large cx
!       Routine called:
!                SPHY for computing the spherical Bessel
!                functions of the second kind
!       ========================================================
!
      IMPLICIT NONE
!*--RMN2L1463
      DOUBLE PRECISION a0 , b0 , C , cx , Df , dy , eps , eps1 , eps2 , &
                     & r , r0 , R2d , R2f , reg , suc , sud , sw , sy , &
                     & X
      INTEGER Id , id1 , id2 , ip , j , k , Kd , l , lg , M , N , nm ,  &
            & nm1 , nm2 , np
      DIMENSION Df(200) , sy(0:251) , dy(0:251)
      eps = 1.0D-14
      ip = 1
      nm1 = INT((N-M)/2)
      IF ( N-M==2*nm1 ) ip = 0
      nm = 25 + nm1 + INT(C)
      reg = 1.0D0
      IF ( M+nm>80 ) reg = 1.0D-200
      nm2 = 2*nm + M
      cx = C*X
      CALL SPHY(nm2,cx,nm2,sy,dy)
      r0 = reg
      DO j = 1 , 2*M + ip
         r0 = r0*j
      ENDDO
      r = r0
      suc = r*Df(1)
      sw = 0.0D0
      DO k = 2 , nm
         r = r*(M+k-1.0)*(M+k+ip-1.5D0)/(k-1.0D0)/(k+ip-1.5D0)
         suc = suc + r*Df(k)
         IF ( k>nm1 .AND. DABS(suc-sw)<DABS(suc)*eps ) GOTO 100
         sw = suc
      ENDDO
 100  a0 = (1.0D0-Kd/(X*X))**(0.5D0*M)/suc
      R2f = 0.0D0
      eps1 = 0.0D0
      np = 0
      DO k = 1 , nm
         l = 2*k + M - N - 2 + ip
         lg = 1
         IF ( l/=4*INT(l/4) ) lg = -1
         IF ( k==1 ) THEN
            r = r0
         ELSE
            r = r*(M+k-1.0)*(M+k+ip-1.5D0)/(k-1.0D0)/(k+ip-1.5D0)
         ENDIF
         np = M + 2*k - 2 + ip
         R2f = R2f + lg*r*(Df(k)*sy(np))
         eps1 = DABS(R2f-sw)
         IF ( k>nm1 .AND. eps1<DABS(R2f)*eps ) GOTO 200
         sw = R2f
      ENDDO
 200  id1 = INT(LOG10(eps1/DABS(R2f)+eps))
      R2f = R2f*a0
      IF ( np>=nm2 ) THEN
         Id = 10
         RETURN
      ENDIF
      b0 = Kd*M/X**3.0D0/(1.0-Kd/(X*X))*R2f
      sud = 0.0D0
      eps2 = 0.0D0
      DO k = 1 , nm
         l = 2*k + M - N - 2 + ip
         lg = 1
         IF ( l/=4*INT(l/4) ) lg = -1
         IF ( k==1 ) THEN
            r = r0
         ELSE
            r = r*(M+k-1.0)*(M+k+ip-1.5D0)/(k-1.0D0)/(k+ip-1.5D0)
         ENDIF
         np = M + 2*k - 2 + ip
         sud = sud + lg*r*(Df(k)*dy(np))
         eps2 = DABS(sud-sw)
         IF ( k>nm1 .AND. eps2<DABS(sud)*eps ) GOTO 300
         sw = sud
      ENDDO
 300  R2d = b0 + a0*C*sud
      id2 = INT(LOG10(eps2/DABS(sud)+eps))
      Id = MAX(id1,id2)
      END
 
 
 
!       **********************************
 
      SUBROUTINE PSI_SPEC(X,Ps)
!
!       ======================================
!       Purpose: Compute Psi function
!       Input :  x  --- Argument of psi(x)
!       Output:  PS --- psi(x)
!       ======================================
!
      IMPLICIT NONE
!*--PSI_SPEC1557
      DOUBLE PRECISION a1 , a2 , a3 , a4 , a5 , a6 , a7 , a8 , el , pi ,&
                     & Ps , s , X , x2 , xa
      INTEGER k , n
      xa = DABS(X)
      pi = 3.141592653589793D0
      el = .5772156649015329D0
      s = 0.0D0
      IF ( X==INT(X) .AND. X<=0.0 ) THEN
         Ps = 1.0D+300
         RETURN
      ELSEIF ( xa==INT(xa) ) THEN
         n = xa
         DO k = 1 , n - 1
            s = s + 1.0D0/k
         ENDDO
         Ps = -el + s
      ELSEIF ( xa+.5==INT(xa+.5) ) THEN
         n = xa - .5
         DO k = 1 , n
            s = s + 1.0/(2.0D0*k-1.0D0)
         ENDDO
         Ps = -el + 2.0D0*s - 1.386294361119891D0
      ELSE
         IF ( xa<10.0 ) THEN
            n = 10 - INT(xa)
            DO k = 0 , n - 1
               s = s + 1.0D0/(xa+k)
            ENDDO
            xa = xa + n
         ENDIF
         x2 = 1.0D0/(xa*xa)
         a1 = -.8333333333333D-01
         a2 = .83333333333333333D-02
         a3 = -.39682539682539683D-02
         a4 = .41666666666666667D-02
         a5 = -.75757575757575758D-02
         a6 = .21092796092796093D-01
         a7 = -.83333333333333333D-01
         a8 = .4432598039215686D0
         Ps = DLOG(xa) - .5D0/xa +                                      &
            & x2*(((((((a8*x2+a7)*x2+a6)*x2+a5)*x2+a4)*x2+a3)*x2+a2)    &
            & *x2+a1)
         Ps = Ps - s
      ENDIF
      IF ( X<0.0 ) Ps = Ps - pi*DCOS(pi*X)/DSIN(pi*X) - 1.0D0/X
      END
 
!       **********************************
 
      SUBROUTINE CVA2(Kd,M,Q,A)
!
!       ======================================================
!       Purpose: Calculate a specific characteristic value of
!                Mathieu functions
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!                KD --- Case code
!                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
!                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
!                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
!                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
!       Output:  A  --- Characteristic value
!       Routines called:
!             (1) REFINE for finding accurate characteristic
!                 value using an iteration method
!             (2) CV0 for finding initial characteristic
!                 values using polynomial approximation
!             (3) CVQM for computing initial characteristic
!                 values for q ≤ 3*m
!             (3) CVQL for computing initial characteristic
!                 values for q ≥ m*m
!       ======================================================
!
      IMPLICIT NONE
!*--CVA21635
      DOUBLE PRECISION A , a1 , a2 , delta , Q , q1 , q2 , qq
      INTEGER i , iflag , Kd , M , ndiv , nn
      IF ( M<=12 .OR. Q<=3.0*M .OR. Q>M*M ) THEN
         CALL CV0(Kd,M,Q,A)
         IF ( Q/=0.0D0 .AND. M/=2 ) CALL REFINE(Kd,M,Q,A)
         IF ( Q>2.0D-3 .AND. M==2 ) CALL REFINE(Kd,M,Q,A)
      ELSE
         ndiv = 10
         delta = (M-3.0)*M/ndiv
         IF ( (Q-3.0*M)<=(M*M-Q) ) THEN
 20         nn = INT((Q-3.0*M)/delta) + 1
            delta = (Q-3.0*M)/nn
            q1 = 2.0*M
            CALL CVQM(M,q1,a1)
            q2 = 3.0*M
            CALL CVQM(M,q2,a2)
            qq = 3.0*M
            DO i = 1 , nn
               qq = qq + delta
               A = (a1*q2-a2*q1+(a2-a1)*qq)/(q2-q1)
               iflag = 1
               IF ( i==nn ) iflag = -1
               CALL REFINE(Kd,M,qq,A)
               q1 = q2
               q2 = qq
               a1 = a2
               a2 = A
            ENDDO
            IF ( iflag==-10 ) THEN
               ndiv = ndiv*2
               delta = (M-3.0)*M/ndiv
               GOTO 20
            ENDIF
         ELSE
 40         nn = INT((M*M-Q)/delta) + 1
            delta = (M*M-Q)/nn
            q1 = M*(M-1.0)
            CALL CVQL(Kd,M,q1,a1)
            q2 = M*M
            CALL CVQL(Kd,M,q2,a2)
            qq = M*M
            DO i = 1 , nn
               qq = qq - delta
               A = (a1*q2-a2*q1+(a2-a1)*qq)/(q2-q1)
               iflag = 1
               IF ( i==nn ) iflag = -1
               CALL REFINE(Kd,M,qq,A)
               q1 = q2
               q2 = qq
               a1 = a2
               a2 = A
            ENDDO
            IF ( iflag==-10 ) THEN
               ndiv = ndiv*2
               delta = (M-3.0)*M/ndiv
               GOTO 40
            ENDIF
         ENDIF
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE LPMNS(M,N,X,Pm,Pd)
!
!       ========================================================
!       Purpose: Compute associated Legendre functions Pmn(x)
!                and Pmn'(x) for a given order
!       Input :  x --- Argument of Pmn(x)
!                m --- Order of Pmn(x),  m = 0,1,2,...,n
!                n --- Degree of Pmn(x), n = 0,1,2,...,N
!       Output:  PM(n) --- Pmn(x)
!                PD(n) --- Pmn'(x)
!       ========================================================
!
      IMPLICIT NONE
!*--LPMNS1717
      INTEGER k , M , N
      DOUBLE PRECISION Pd , Pm , pm0 , pm1 , pm2 , pmk , X , x0
      DIMENSION Pm(0:N) , Pd(0:N)
      DO k = 0 , N
         Pm(k) = 0.0D0
         Pd(k) = 0.0D0
      ENDDO
      IF ( DABS(X)==1.0D0 ) THEN
         DO k = 0 , N
            IF ( M==0 ) THEN
               Pm(k) = 1.0D0
               Pd(k) = 0.5D0*k*(k+1.0)
               IF ( X<0.0 ) THEN
                  Pm(k) = (-1)**k*Pm(k)
                  Pd(k) = (-1)**(k+1)*Pd(k)
               ENDIF
            ELSEIF ( M==1 ) THEN
               Pd(k) = 1.0D+300
            ELSEIF ( M==2 ) THEN
               Pd(k) = -0.25D0*(k+2.0)*(k+1.0)*k*(k-1.0)
               IF ( X<0.0 ) Pd(k) = (-1)**(k+1)*Pd(k)
            ENDIF
         ENDDO
         RETURN
      ENDIF
      x0 = DABS(1.0D0-X*X)
      pm0 = 1.0D0
      pmk = pm0
      DO k = 1 , M
         pmk = (2.0D0*k-1.0D0)*DSQRT(x0)*pm0
         pm0 = pmk
      ENDDO
      pm1 = (2.0D0*M+1.0D0)*X*pm0
      Pm(M) = pmk
      Pm(M+1) = pm1
      DO k = M + 2 , N
         pm2 = ((2.0D0*k-1.0D0)*X*pm1-(k+M-1.0D0)*pmk)/(k-M)
         Pm(k) = pm2
         pmk = pm1
         pm1 = pm2
      ENDDO
      Pd(0) = ((1.0D0-M)*Pm(1)-X*Pm(0))/(X*X-1.0)
      DO k = 1 , N
         Pd(k) = (k*X*Pm(k)-(k+M)*Pm(k-1))/(X*X-1.0D0)
      ENDDO
      DO k = 1 , N
         Pm(k) = (-1)**M*Pm(k)
         Pd(k) = (-1)**M*Pd(k)
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE CERF(Z,Cer,Cder)
!
!       ==========================================================
!       Purpose: Compute complex Error function erf(z) & erf'(z)
!       Input:   z   --- Complex argument of erf(z)
!                x   --- Real part of z
!                y   --- Imaginary part of z
!       Output:  CER --- erf(z)
!                CDER --- erf'(z)
!       ==========================================================
      IMPLICIT NONE
!*--CERF1785
      DOUBLE PRECISION c0 , cs , ei1 , ei2 , eps , er , er0 , er1 ,     &
                     & er2 , eri , err , pi , r , ss , w , w1 , w2 , x ,&
                     & x2 , y
      INTEGER k , n
      COMPLEX*16 Z , Cer , Cder
      eps = 1.0D-12
      pi = 3.141592653589793D0
      x = DBLE(Z)
      y = DIMAG(Z)
      x2 = x*x
      IF ( x<=3.5D0 ) THEN
         er = 1.0D0
         r = 1.0D0
         w = 0.0D0
         DO k = 1 , 100
            r = r*x2/(k+0.5D0)
            er = er + r
            IF ( DABS(er-w)<=eps*DABS(er) ) GOTO 50
            w = er
         ENDDO
 50      c0 = 2.0D0/DSQRT(pi)*x*EXP(-x2)
         er0 = c0*er
      ELSE
         er = 1.0D0
         r = 1.0D0
         DO k = 1 , 12
            r = -r*(k-0.5D0)/x2
            er = er + r
         ENDDO
         c0 = EXP(-x2)/(x*DSQRT(pi))
         er0 = 1.0D0 - c0*er
      ENDIF
      IF ( y==0.0D0 ) THEN
         err = er0
         eri = 0.0D0
      ELSE
         cs = DCOS(2.0D0*x*y)
         ss = DSIN(2.0D0*x*y)
         er1 = EXP(-x2)*(1.0D0-cs)/(2.0D0*pi*x)
         ei1 = EXP(-x2)*ss/(2.0D0*pi*x)
         er2 = 0.0D0
         w1 = 0.0D0
         DO n = 1 , 100
            er2 = er2 + EXP(-.25D0*n*n)/(n*n+4.0D0*x2)                  &
                & *(2.0D0*x-2.0D0*x*DCOSH(n*y)*cs+n*DSINH(n*y)*ss)
            IF ( DABS((er2-w1)/er2)<eps ) GOTO 100
            w1 = er2
         ENDDO
 100     c0 = 2.0D0*EXP(-x2)/pi
         err = er0 + er1 + c0*er2
         ei2 = 0.0D0
         w2 = 0.0D0
         DO n = 1 , 100
            ei2 = ei2 + EXP(-.25D0*n*n)/(n*n+4.0D0*x2)                  &
                & *(2.0D0*x*DCOSH(n*y)*ss+n*DSINH(n*y)*cs)
            IF ( DABS((ei2-w2)/ei2)<eps ) GOTO 150
            w2 = ei2
         ENDDO
 150     eri = ei1 + c0*ei2
      ENDIF
      Cer = DCMPLX(err,eri)
      Cder = 2.0D0/DSQRT(pi)*EXP(-Z*Z)
      END
 
!       **********************************
 
      SUBROUTINE RSWFP(M,N,C,X,Cv,Kf,R1f,R1d,R2f,R2d)
!
!       ==============================================================
!       Purpose: Compute prolate spheriodal radial functions of the
!                first and second kinds, and their derivatives
!       Input :  m  --- Mode parameter, m = 0,1,2,...
!                n  --- Mode parameter, n = m,m+1,m+2,...
!                c  --- Spheroidal parameter
!                x  --- Argument of radial function ( x > 1.0 )
!                cv --- Characteristic value
!                KF --- Function code
!                       KF=1 for the first kind
!                       KF=2 for the second kind
!                       KF=3 for both the first and second kinds
!       Output:  R1F --- Radial function of the first kind
!                R1D --- Derivative of the radial function of
!                        the first kind
!                R2F --- Radial function of the second kind
!                R2D --- Derivative of the radial function of
!                        the second kind
!       Routines called:
!            (1) SDMN for computing expansion coefficients dk
!            (2) RMN1 for computing prolate and oblate radial
!                functions of the first kind
!            (3) RMN2L for computing prolate and oblate radial
!                functions of the second kind for a large argument
!            (4) RMN2SP for computing the prolate radial function
!                of the second kind for a small argument
!       ==============================================================
!
      IMPLICIT NONE
!*--RSWFP1886
      DOUBLE PRECISION C , Cv , df , R1d , R1f , R2d , R2f , X
      INTEGER id , kd , Kf , M , N
      DIMENSION df(200)
      kd = 1
      CALL SDMN(M,N,C,Cv,kd,df)
      IF ( Kf/=2 ) CALL RMN1(M,N,C,X,df,kd,R1f,R1d)
      IF ( Kf>1 ) THEN
         CALL RMN2L(M,N,C,X,df,kd,R2f,R2d,id)
         IF ( id>-8 ) CALL RMN2SP(M,N,C,X,Cv,df,kd,R2f,R2d)
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE JYNDD(N,X,Bjn,Djn,Fjn,Byn,Dyn,Fyn)
!
!       ===========================================================
!       Purpose: Compute Bessel functions Jn(x) and Yn(x), and
!                their first and second derivatives
!       Input:   x   ---  Argument of Jn(x) and Yn(x) ( x > 0 )
!                n   ---  Order of Jn(x) and Yn(x)
!       Output:  BJN ---  Jn(x)
!                DJN ---  Jn'(x)
!                FJN ---  Jn"(x)
!                BYN ---  Yn(x)
!                DYN ---  Yn'(x)
!                FYN ---  Yn"(x)
!       Routines called:
!                JYNBH to compute Jn and Yn
!       ===========================================================
!
      IMPLICIT NONE
!*--JYNDD1924
      DOUBLE PRECISION bj , Bjn , by , Byn , Djn , Dyn , Fjn , Fyn , X
      INTEGER N , nm
      DIMENSION bj(2) , by(2)
      CALL JYNBH(N+1,N,X,nm,bj,by)
!       Compute derivatives by differentiation formulas
      Bjn = bj(1)
      Byn = by(1)
      Djn = -bj(2) + N*bj(1)/X
      Dyn = -by(2) + N*by(1)/X
      Fjn = (N*N/(X*X)-1.0D0)*Bjn - Djn/X
      Fyn = (N*N/(X*X)-1.0D0)*Byn - Dyn/X
      END
 
 
!       **********************************
 
      SUBROUTINE GAM0(X,Ga)
!
!       ================================================
!       Purpose: Compute gamma function Г(x)
!       Input :  x  --- Argument of Г(x)  ( |x| ≤ 1 )
!       Output:  GA --- Г(x)
!       ================================================
!
      IMPLICIT NONE
!*--GAM01953
      DOUBLE PRECISION g , Ga , gr , X
      INTEGER k
      DIMENSION g(25)
      DATA g/1.0D0 , 0.5772156649015329D0 , -0.6558780715202538D0 ,     &
         & -0.420026350340952D-1 , 0.1665386113822915D0 ,               &
         & -.421977345555443D-1 , -.96219715278770D-2 ,                 &
         & .72189432466630D-2 , -.11651675918591D-2 ,                   &
         & -.2152416741149D-3 , .1280502823882D-3 , -.201348547807D-4 , &
         & -.12504934821D-5 , .11330272320D-5 , -.2056338417D-6 ,       &
         & .61160950D-8 , .50020075D-8 , -.11812746D-8 , .1043427D-9 ,  &
         & .77823D-11 , -.36968D-11 , .51D-12 , -.206D-13 , -.54D-14 ,  &
         & .14D-14/
      gr = (25)
      DO k = 24 , 1 , -1
         gr = gr*X + g(k)
      ENDDO
      Ga = 1.0D0/(gr*X)
      END
 
 
!       **********************************
 
      SUBROUTINE CISIB(X,Ci,Si)
!
!       =============================================
!       Purpose: Compute cosine and sine integrals
!                Si(x) and Ci(x) ( x ≥ 0 )
!       Input :  x  --- Argument of Ci(x) and Si(x)
!       Output:  CI --- Ci(x)
!                SI --- Si(x)
!       =============================================
!
      IMPLICIT NONE
!*--CISIB1990
      DOUBLE PRECISION Ci , fx , gx , Si , X , x2
      x2 = X*X
      IF ( X==0.0 ) THEN
         Ci = -1.0D+300
         Si = 0.0D0
      ELSEIF ( X<=1.0D0 ) THEN
         Ci = ((((-3.0D-8*x2+3.10D-6)*x2-2.3148D-4)*x2+1.041667D-2)     &
            & *x2-0.25)*x2 + 0.577215665D0 + LOG(X)
         Si = ((((3.1D-7*x2-2.834D-5)*x2+1.66667D-003)*x2-5.555556D-002)&
            & *x2+1.0)*X
      ELSE
         fx = ((((x2+38.027264D0)*x2+265.187033D0)*x2+335.67732D0)      &
            & *x2+38.102495D0)                                          &
            & /((((x2+40.021433D0)*x2+322.624911D0)*x2+570.23628D0)     &
            & *x2+157.105423D0)
         gx = ((((x2+42.242855D0)*x2+302.757865D0)*x2+352.018498D0)     &
            & *x2+21.821899D0)                                          &
            & /((((x2+48.196927D0)*x2+482.485984D0)*x2+1114.978885D0)   &
            & *x2+449.690326D0)/X
         Ci = fx*SIN(X)/X - gx*COS(X)/X
         Si = 1.570796327D0 - fx*COS(X)/X - gx*SIN(X)/X
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE EULERA(N,En)
!
!       ======================================
!       Purpose: Compute Euler number En
!       Input :  n --- Serial number
!       Output:  EN(n) --- En
!       ======================================
!
      IMPLICIT NONE
!*--EULERA2029
      DOUBLE PRECISION En , r , s
      INTEGER j , k , m , N
      DIMENSION En(0:N)
      En(0) = 1.0D0
      DO m = 1 , N/2
         s = 1.0D0
         DO k = 1 , m - 1
            r = 1.0D0
            DO j = 1 , 2*k
               r = r*(2.0D0*m-2.0D0*k+j)/j
            ENDDO
            s = s + r*En(2*k)
         ENDDO
         En(2*m) = -s
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE REFINE(Kd,M,Q,A)
!
!       =====================================================
!       Purpose: calculate the accurate characteristic value
!                by the secant method
!       Input :  m --- Order of Mathieu functions
!                q --- Parameter of Mathieu functions
!                A --- Initial characteristic value
!       Output:  A --- Refineed characteristic value
!       Routine called:  CVF for computing the value of F for
!                        characteristic equation
!       ========================================================
!
      IMPLICIT NONE
!*--REFINE2066
      DOUBLE PRECISION A , ca , delta , eps , f , f0 , f1 , Q , x , x0 ,&
                     & x1
      INTEGER it , Kd , M , mj
      eps = 1.0D-14
      mj = 10 + M
      ca = A
      delta = 0.0D0
      x0 = A
      CALL CVF(Kd,M,Q,x0,mj,f0)
      x1 = 1.002*A
      CALL CVF(Kd,M,Q,x1,mj,f1)
      DO it = 1 , 100
         mj = mj + 1
         x = x1 - (x1-x0)/(1.0D0-f0/f1)
         CALL CVF(Kd,M,Q,x,mj,f)
         IF ( ABS(1.0-x1/x)<eps .OR. f==0.0 ) GOTO 100
         x0 = x1
         f0 = f1
         x1 = x
         f1 = f
      ENDDO
 100  A = x
      END
 
 
 
!       **********************************
 
      SUBROUTINE CISIA(X,Ci,Si)
!
!       =============================================
!       Purpose: Compute cosine and sine integrals
!                Si(x) and Ci(x)  ( x ≥ 0 )
!       Input :  x  --- Argument of Ci(x) and Si(x)
!       Output:  CI --- Ci(x)
!                SI --- Si(x)
!       =============================================
!
      IMPLICIT NONE
!*--CISIA2109
      DOUBLE PRECISION bj , Ci , el , eps , p2 , Si , X , x2 , xa ,     &
                     & xa0 , xa1 , xcs , xf , xg , xg1 , xg2 , xr , xs ,&
                     & xss
      INTEGER k , m
      DIMENSION bj(101)
      p2 = 1.570796326794897D0
      el = .5772156649015329D0
      eps = 1.0D-15
      x2 = X*X
      IF ( X==0.0D0 ) THEN
         Ci = -1.0D+300
         Si = 0.0D0
      ELSEIF ( X<=16.0D0 ) THEN
         xr = -.25D0*x2
         Ci = el + DLOG(X) + xr
         DO k = 2 , 40
            xr = -.5D0*xr*(k-1)/(k*k*(2*k-1))*x2
            Ci = Ci + xr
            IF ( DABS(xr)<DABS(Ci)*eps ) GOTO 50
         ENDDO
 50      xr = X
         Si = X
         DO k = 1 , 40
            xr = -.5D0*xr*(2*k-1)/k/(4*k*k+4*k+1)*x2
            Si = Si + xr
            IF ( DABS(xr)<DABS(Si)*eps ) RETURN
         ENDDO
      ELSEIF ( X<=32.0D0 ) THEN
         m = INT(47.2+.82*X)
         xa1 = 0.0D0
         xa0 = 1.0D-100
         DO k = m , 1 , -1
            xa = 4.0D0*k*xa0/X - xa1
            bj(k) = xa
            xa1 = xa0
            xa0 = xa
         ENDDO
         xs = bj(1)
         DO k = 3 , m , 2
            xs = xs + 2.0D0*bj(k)
         ENDDO
         bj(1) = bj(1)/xs
         DO k = 2 , m
            bj(k) = bj(k)/xs
         ENDDO
         xr = 1.0D0
         xg1 = bj(1)
         DO k = 2 , m
            xr = .25D0*xr*(2.0*k-3.0)**2/((k-1.0)*(2.0*k-1.0)**2)*X
            xg1 = xg1 + bj(k)*xr
         ENDDO
         xr = 1.0D0
         xg2 = bj(1)
         DO k = 2 , m
            xr = .25D0*xr*(2.0*k-5.0)**2/((k-1.0)*(2.0*k-3.0)**2)*X
            xg2 = xg2 + bj(k)*xr
         ENDDO
         xcs = DCOS(X/2.0D0)
         xss = DSIN(X/2.0D0)
         Ci = el + DLOG(X) - X*xss*xg1 + 2*xcs*xg2 - 2*xcs*xcs
         Si = X*xcs*xg1 + 2*xss*xg2 - DSIN(X)
      ELSE
         xr = 1.0D0
         xf = 1.0D0
         DO k = 1 , 9
            xr = -2.0D0*xr*k*(2*k-1)/x2
            xf = xf + xr
         ENDDO
         xr = 1.0D0/X
         xg = xr
         DO k = 1 , 8
            xr = -2.0D0*xr*(2*k+1)*k/x2
            xg = xg + xr
         ENDDO
         Ci = xf*DSIN(X)/X - xg*DCOS(X)/X
         Si = p2 - xf*DCOS(X)/X - xg*DSIN(X)/X
      ENDIF
      END
 
 
!       **********************************
 
      SUBROUTINE ITSL0(X,Tl0)
!
!       ===========================================================
!       Purpose: Evaluate the integral of modified Struve function
!                L0(t) with respect to t from 0 to x
!       Input :  x   --- Upper limit  ( x ≥ 0 )
!       Output:  TL0 --- Integration of L0(t) from 0 to x
!       ===========================================================
!
      IMPLICIT NONE
!*--ITSL02205
      DOUBLE PRECISION a , a0 , a1 , af , el , pi , r , rd , s , s0 ,   &
                     & ti , Tl0 , X
      INTEGER k
      DIMENSION a(18)
      pi = 3.141592653589793D0
      r = 1.0D0
      IF ( X<=20.0 ) THEN
         s = 0.5D0
         DO k = 1 , 100
            rd = 1.0D0
            IF ( k==1 ) rd = 0.5D0
            r = r*rd*k/(k+1.0D0)*(X/(2.0D0*k+1.0D0))**2
            s = s + r
            IF ( DABS(r/s)<1.0D-12 ) GOTO 50
         ENDDO
 50      Tl0 = 2.0D0/pi*X*X*s
      ELSE
         s = 1.0D0
         DO k = 1 , 10
            r = r*k/(k+1.0D0)*((2.0D0*k+1.0D0)/X)**2
            s = s + r
            IF ( DABS(r/s)<1.0D-12 ) GOTO 100
         ENDDO
 100     el = .57721566490153D0
         s0 = -s/(pi*X*X) + 2.0D0/pi*(DLOG(2.0D0*X)+el)
         a0 = 1.0D0
         a1 = 5.0D0/8.0D0
         a(1) = a1
         DO k = 1 , 10
            af = ((1.5D0*(k+.50D0)*(k+5.0D0/6.0D0)*a1-.5D0*(k+.5D0)     &
               & **2*(k-.5D0)*a0))/(k+1.0D0)
            a(k+1) = af
            a0 = a1
            a1 = af
         ENDDO
         ti = 1.0D0
         r = 1.0D0
         DO k = 1 , 11
            r = r/X
            ti = ti + a(k)*r
         ENDDO
         Tl0 = ti/DSQRT(2*pi*X)*EXP(X) + s0
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE CLQN(N,X,Y,Cqn,Cqd)
!
!       ==================================================
!       Purpose: Compute the Legendre functions Qn(z) and
!                their derivatives Qn'(z) for a complex
!                argument
!       Input :  x --- Real part of z
!                y --- Imaginary part of z
!                n --- Degree of Qn(z), n = 0,1,2,...
!       Output:  CQN(n) --- Qn(z)
!                CQD(n) --- Qn'(z)
!       ==================================================
!
      IMPLICIT NONE
!*--CLQN2270
      COMPLEX*16 cq0 , cq1 , Cqd , cqf0 , cqf1 , cqf2 , Cqn , z
      INTEGER k , km , ls , N
      DOUBLE PRECISION X , Y
      DIMENSION Cqn(0:N) , Cqd(0:N)
      z = DCMPLX(X,Y)
      IF ( z==1.0D0 ) THEN
         DO k = 0 , N
            Cqn(k) = (1.0D+300,0.0D0)
            Cqd(k) = (1.0D+300,0.0D0)
         ENDDO
         RETURN
      ENDIF
      ls = 1
      IF ( ABS(z)>1.0D0 ) ls = -1
      cq0 = 0.5D0*LOG(ls*(1.0D0+z)/(1.0D0-z))
      cq1 = z*cq0 - 1.0D0
      Cqn(0) = cq0
      Cqn(1) = cq1
      IF ( ABS(z)<1.0001D0 ) THEN
         cqf0 = cq0
         cqf1 = cq1
         DO k = 2 , N
            cqf2 = ((2.0D0*k-1.0D0)*z*cqf1-(k-1.0D0)*cqf0)/k
            Cqn(k) = cqf2
            cqf0 = cqf1
            cqf1 = cqf2
         ENDDO
      ELSE
         IF ( ABS(z)>1.1D0 ) THEN
            km = 40 + N
         ELSE
            km = (40+N)*INT(-1.0-1.8*LOG(ABS(z-1.0)))
         ENDIF
         cqf2 = 0.0D0
         cqf1 = 1.0D0
         DO k = km , 0 , -1
            cqf0 = ((2*k+3.0D0)*z*cqf1-(k+2.0D0)*cqf2)/(k+1.0D0)
            IF ( k<=N ) Cqn(k) = cqf0
            cqf2 = cqf1
            cqf1 = cqf0
         ENDDO
         DO k = 0 , N
            Cqn(k) = Cqn(k)*cq0/cqf0
         ENDDO
      ENDIF
      Cqd(0) = (Cqn(1)-z*Cqn(0))/(z*z-1.0D0)
      DO k = 1 , N
         Cqd(k) = (k*z*Cqn(k)-k*Cqn(k-1))/(z*z-1.0D0)
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE AIRYZO(Nt,Kf,Xa,Xb,Xc,Xd)
!
!       ========================================================
!       Purpose: Compute the first NT zeros of Airy functions
!                Ai(x) and Ai'(x), a and a', and the associated
!                values of Ai(a') and Ai'(a); and the first NT
!                zeros of Airy functions Bi(x) and Bi'(x), b and
!                b', and the associated values of Bi(b') and
!                Bi'(b)
!       Input :  NT    --- Total number of zeros
!                KF    --- Function code
!                          KF=1 for Ai(x) and Ai'(x)
!                          KF=2 for Bi(x) and Bi'(x)
!       Output:  XA(m) --- a, the m-th zero of Ai(x) or
!                          b, the m-th zero of Bi(x)
!                XB(m) --- a', the m-th zero of Ai'(x) or
!                          b', the m-th zero of Bi'(x)
!                XC(m) --- Ai(a') or Bi(b')
!                XD(m) --- Ai'(a) or Bi'(b)
!                          ( m --- Serial number of zeros )
!       Routine called: AIRYB for computing Airy functions and
!                       their derivatives
!       =======================================================
!
      IMPLICIT NONE
!*--AIRYZO2352
      DOUBLE PRECISION ad , ai , bd , bi , err , pi , rt , rt0 , u ,    &
                     & u1 , x , Xa , Xb , Xc , Xd
      INTEGER i , Kf , Nt
      DIMENSION Xa(Nt) , Xb(Nt) , Xc(Nt) , Xd(Nt)
      pi = 3.141592653589793D0
      rt = 0.0D0
      DO i = 1 , Nt
         rt0 = 0D0
         IF ( Kf==1 ) THEN
            u = 3.0D0*pi*(4.0D0*i-1)/8.0D0
            u1 = 1/(u*u)
         ELSEIF ( Kf==2 ) THEN
            IF ( i==1 ) THEN
               rt0 = -1.17371D0
            ELSE
               u = 3.0D0*pi*(4.0D0*i-3.0D0)/8.0D0
               u1 = 1/(u*u)
            ENDIF
         ENDIF
!             DLMF 9.9.18
         IF ( rt0==0 ) rt0 = -(u*u)**(1.0D0/3.0D0)                      &
                           & *(+1D0+u1*(5D0/48D0+u1*                    &
                           & (-5D0/36D0+u1*(77125D0/82944D0+            &
                           & u1*(-108056875D0/6967296D0)))))
 50      x = rt0
         CALL AIRYB(x,ai,bi,ad,bd)
         IF ( Kf==1 ) rt = rt0 - ai/ad
         IF ( Kf==2 ) rt = rt0 - bi/bd
         err = DABS((rt-rt0)/rt)
         IF ( err>1.D-12 ) THEN
            rt0 = rt
            GOTO 50
         ELSE
            Xa(i) = rt
            IF ( err>1D-14 ) CALL AIRYB(rt,ai,bi,ad,bd)
            IF ( Kf==1 ) Xd(i) = ad
            IF ( Kf==2 ) Xd(i) = bd
         ENDIF
      ENDDO
      DO i = 1 , Nt
         rt0 = 0D0
         IF ( Kf==1 ) THEN
            IF ( i==1 ) THEN
               rt0 = -1.01879D0
            ELSE
               u = 3.0D0*pi*(4.0D0*i-3.0D0)/8.0D0
               u1 = 1/(u*u)
            ENDIF
         ELSEIF ( Kf==2 ) THEN
            IF ( i==1 ) THEN
               rt0 = -2.29444D0
            ELSE
               u = 3.0D0*pi*(4.0D0*i-1.0D0)/8.0D0
               u1 = 1/(u*u)
            ENDIF
         ENDIF
!             DLMF 9.9.19
         IF ( rt0==0 ) rt0 = -(u*u)**(1.0D0/3.0D0)                      &
                           & *(+1D0+u1*(-7D0/48D0+u1*                   &
                           & (+35D0/288D0+u1*(-181223D0/207360D0+       &
                           & u1*(18683371D0/1244160D0)))))
 100     x = rt0
         CALL AIRYB(x,ai,bi,ad,bd)
         IF ( Kf==1 ) rt = rt0 - ad/(ai*x)
         IF ( Kf==2 ) rt = rt0 - bd/(bi*x)
         err = DABS((rt-rt0)/rt)
         IF ( err>1.0D-12 ) THEN
            rt0 = rt
            GOTO 100
         ELSE
            Xb(i) = rt
            IF ( err>1D-14 ) CALL AIRYB(rt,ai,bi,ad,bd)
            IF ( Kf==1 ) Xc(i) = ai
            IF ( Kf==2 ) Xc(i) = bi
         ENDIF
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE ERROR(X,Err)
!
!       =========================================
!       Purpose: Compute error function erf(x)
!       Input:   x   --- Argument of erf(x)
!       Output:  ERR --- erf(x)
!       =========================================
!
      IMPLICIT NONE
!*--ERROR2447
      DOUBLE PRECISION c0 , eps , er , Err , pi , r , X , x2
      INTEGER k
      eps = 1.0D-15
      pi = 3.141592653589793D0
      x2 = X*X
      IF ( DABS(X)<3.5D0 ) THEN
         er = 1.0D0
         r = 1.0D0
         DO k = 1 , 50
            r = r*x2/(k+0.5D0)
            er = er + r
            IF ( DABS(r)<=DABS(er)*eps ) GOTO 50
         ENDDO
 50      c0 = 2.0D0/DSQRT(pi)*X*EXP(-x2)
         Err = c0*er
      ELSE
         er = 1.0D0
         r = 1.0D0
         DO k = 1 , 12
            r = -r*(k-0.5D0)/x2
            er = er + r
         ENDDO
         c0 = EXP(-x2)/(DABS(X)*DSQRT(pi))
         Err = 1.0D0 - c0*er
         IF ( X<0.0 ) Err = -Err
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE CERROR(Z,Cer)
!
!       ====================================================
!       Purpose: Compute error function erf(z) for a complex
!                argument (z=x+iy)
!       Input :  z   --- Complex argument
!       Output:  CER --- erf(z)
!       ====================================================
!
      IMPLICIT NONE
!*--CERROR2491
      COMPLEX*16 c0 , Cer , cl , cr , cs , Z , z1
      INTEGER k
      DOUBLE PRECISION a0 , pi
      a0 = ABS(Z)
      c0 = EXP(-Z*Z)
      pi = 3.141592653589793D0
      z1 = Z
      IF ( DBLE(Z)<0.0 ) z1 = -Z
!
!       Cutoff radius R = 4.36; determined by balancing rounding error
!       and asymptotic expansion error, see below.
!
!       The resulting maximum global accuracy expected is around 1e-8
!
      IF ( a0<=4.36D0 ) THEN
!
!          Rounding error in the Taylor expansion is roughly
!
!          ~ R*R * EPSILON * R**(2 R**2) / (2 R**2 Gamma(R**2 + 1/2))
!
         cs = z1
         cr = z1
         DO k = 1 , 120
            cr = cr*z1*z1/(k+0.5D0)
            cs = cs + cr
            IF ( ABS(cr/cs)<1.0D-15 ) GOTO 50
         ENDDO
 50      Cer = 2.0D0*c0*cs/DSQRT(pi)
      ELSE
         cl = 1.0D0/z1
         cr = cl
!
!          Asymptotic series; maximum K must be at most ~ R^2.
!
!          The maximum accuracy obtainable from this expansion is roughly
!
!          ~ Gamma(2R**2 + 2) / (
!                   (2 R**2)**(R**2 + 1/2) Gamma(R**2 + 3/2) 2**(R**2 + 1/2))
!
         DO k = 1 , 20
            cr = -cr*(k-0.5D0)/(z1*z1)
            cl = cl + cr
            IF ( ABS(cr/cl)<1.0D-15 ) GOTO 100
         ENDDO
 100     Cer = 1.0D0 - c0*cl/DSQRT(pi)
      ENDIF
      IF ( DBLE(Z)<0.0 ) Cer = -Cer
      END
 
 
 
!       **********************************
 
      SUBROUTINE EULERB(N,En)
!
!       ======================================
!       Purpose: Compute Euler number En
!       Input :  n --- Serial number
!       Output:  EN(n) --- En
!       ======================================
!
      IMPLICIT NONE
!*--EULERB2557
      DOUBLE PRECISION En , hpi , r1 , r2 , s
      INTEGER isgn , k , m , N
      DIMENSION En(0:N)
      hpi = 2.0D0/3.141592653589793D0
      En(0) = 1.0D0
      En(2) = -1.0D0
      r1 = -4.0D0*hpi**3
      DO m = 4 , N , 2
         r1 = -r1*(m-1)*m*hpi*hpi
         r2 = 1.0D0
         isgn = 1.0D0
         DO k = 3 , 1000 , 2
            isgn = -isgn
            s = (1.0D0/k)**(m+1)
            r2 = r2 + isgn*s
            IF ( s<1.0D-15 ) GOTO 50
         ENDDO
 50      En(m) = r1*r2
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE CVA1(Kd,M,Q,Cv)
!
!       ============================================================
!       Purpose: Compute a sequence of characteristic values of
!                Mathieu functions
!       Input :  M  --- Maximum order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!                KD --- Case code
!                       KD=1 for cem(x,q)  ( m = 0,2,4,… )
!                       KD=2 for cem(x,q)  ( m = 1,3,5,… )
!                       KD=3 for sem(x,q)  ( m = 1,3,5,… )
!                       KD=4 for sem(x,q)  ( m = 2,4,6,… )
!       Output:  CV(I) --- Characteristic values; I = 1,2,3,...
!                For KD=1, CV(1), CV(2), CV(3),..., correspond to
!                the characteristic values of cem for m = 0,2,4,...
!                For KD=2, CV(1), CV(2), CV(3),..., correspond to
!                the characteristic values of cem for m = 1,3,5,...
!                For KD=3, CV(1), CV(2), CV(3),..., correspond to
!                the characteristic values of sem for m = 1,3,5,...
!                For KD=4, CV(1), CV(2), CV(3),..., correspond to
!                the characteristic values of sem for m = 0,2,4,...
!       ============================================================
!
      IMPLICIT NONE
!*--CVA12608
      DOUBLE PRECISION Cv , d , e , eps , f , g , h , Q , s , t , t1 ,  &
                     & x1 , xa , xb
      INTEGER i , ic , icm , j , k , k1 , Kd , M , nm , nm1
      DIMENSION g(200) , h(200) , d(500) , e(500) , f(500) , Cv(200)
      eps = 1.0D-14
      icm = INT(M/2) + 1
      IF ( Kd==4 ) icm = M/2
      IF ( Q/=0.0D0 ) THEN
         nm = INT(10+1.5*M+0.5*Q)
         e(1) = 0.0D0
         f(1) = 0.0D0
         IF ( Kd==1 ) THEN
            d(1) = 0.0D0
            DO i = 2 , nm
               d(i) = 4.0D0*(i-1.0D0)**2
               e(i) = Q
               f(i) = Q*Q
            ENDDO
            e(2) = DSQRT(2.0D0)*Q
            f(2) = 2.0D0*Q*Q
         ELSEIF ( Kd/=4 ) THEN
            d(1) = 1.0D0 + (-1)**Kd*Q
            DO i = 2 , nm
               d(i) = (2.0D0*i-1.0D0)**2
               e(i) = Q
               f(i) = Q*Q
            ENDDO
         ELSE
            d(1) = 4.0D0
            DO i = 2 , nm
               d(i) = 4.0D0*i*i
               e(i) = Q
               f(i) = Q*Q
            ENDDO
         ENDIF
         xa = d(nm) + DABS(e(nm))
         xb = d(nm) - DABS(e(nm))
         nm1 = nm - 1
         DO i = 1 , nm1
            t = DABS(e(i)) + DABS(e(i+1))
            t1 = d(i) + t
            IF ( xa<t1 ) xa = t1
            t1 = d(i) - t
            IF ( t1<xb ) xb = t1
         ENDDO
         DO i = 1 , icm
            g(i) = xa
            h(i) = xb
         ENDDO
         DO k = 1 , icm
            DO k1 = k , icm
               IF ( g(k1)<g(k) ) THEN
                  g(k) = g(k1)
                  GOTO 20
               ENDIF
            ENDDO
 20         IF ( k/=1 .AND. h(k)<h(k-1) ) h(k) = h(k-1)
 40         x1 = (g(k)+h(k))/2.0D0
            Cv(k) = x1
            IF ( DABS((g(k)-h(k))/x1)<eps ) THEN
               Cv(k) = x1
            ELSE
               j = 0
               s = 1.0D0
               DO i = 1 , nm
                  IF ( s==0.0D0 ) s = s + 1.0D-30
                  t = f(i)/s
                  s = d(i) - t - x1
                  IF ( s<0.0 ) j = j + 1
               ENDDO
               IF ( j<k ) THEN
                  h(k) = x1
               ELSE
                  g(k) = x1
                  IF ( j>=icm ) THEN
                     g(icm) = x1
                  ELSE
                     IF ( h(j+1)<x1 ) h(j+1) = x1
                     IF ( x1<g(j) ) g(j) = x1
                  ENDIF
               ENDIF
               GOTO 40
            ENDIF
         ENDDO
      ELSEIF ( Kd==1 ) THEN
         DO ic = 1 , icm
            Cv(ic) = 4.0D0*(ic-1.0D0)**2
         ENDDO
      ELSEIF ( Kd/=4 ) THEN
         DO ic = 1 , icm
            Cv(ic) = (2.0D0*ic-1.0D0)**2
         ENDDO
      ELSE
         DO ic = 1 , icm
            Cv(ic) = 4.0D0*ic*ic
         ENDDO
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE ITTIKB(X,Tti,Ttk)
!
!       =========================================================
!       Purpose: Integrate [I0(t)-1]/t with respect to t from 0
!                to x, and K0(t)/t with respect to t from x to ∞
!       Input :  x   --- Variable in the limits  ( x ≥ 0 )
!       Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
!                TTK --- Integration of K0(t)/t from x to ∞
!       =========================================================
!
      IMPLICIT NONE
!*--ITTIKB2724
      DOUBLE PRECISION e0 , el , pi , t , t1 , Tti , Ttk , X , x1
      pi = 3.141592653589793D0
      el = .5772156649015329D0
      IF ( X==0.0D0 ) THEN
         Tti = 0.0D0
      ELSEIF ( X<=5.0D0 ) THEN
         x1 = X/5.0D0
         t = x1*x1
         Tti = (((((((.1263D-3*t+.96442D-3)*t+.968217D-2)*t+.06615507D0)&
             & *t+.33116853D0)*t+1.13027241D0)*t+2.44140746D0)          &
             & *t+3.12499991D0)*t
      ELSE
         t = 5.0D0/X
         Tti = (((((((((2.1945464D0*t-3.5195009D0)*t-11.9094395D0)*t+   &
             & 40.394734D0)*t-48.0524115D0)*t+28.1221478D0)             &
             & *t-8.6556013D0)*t+1.4780044D0)*t-.0493843D0)             &
             & *t+.1332055D0)*t + .3989314D0
         Tti = Tti*EXP(X)/(DSQRT(X)*X)
      ENDIF
      IF ( X==0.0D0 ) THEN
         Ttk = 1.0D+300
      ELSEIF ( X<=2.0D0 ) THEN
         t1 = X/2.0D0
         t = t1*t1
         Ttk = (((((.77D-6*t+.1544D-4)*t+.48077D-3)*t+.925821D-2)       &
             & *t+.10937537D0)*t+.74999993D0)*t
         e0 = el + DLOG(X/2.0D0)
         Ttk = pi*pi/24.0D0 + e0*(.5D0*e0+Tti) - Ttk
      ELSEIF ( X<=4.0D0 ) THEN
         t = 2.0D0/X
         Ttk = (((.06084D0*t-.280367D0)*t+.590944D0)*t-.850013D0)       &
             & *t + 1.234684D0
         Ttk = Ttk*EXP(-X)/(DSQRT(X)*X)
      ELSE
         t = 4.0D0/X
         Ttk = (((((.02724D0*t-.1110396D0)*t+.2060126D0)*t-.2621446D0)  &
             & *t+.3219184D0)*t-.5091339D0)*t + 1.2533141D0
         Ttk = Ttk*EXP(-X)/(DSQRT(X)*X)
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE LQNB(N,X,Qn,Qd)
!
!       ====================================================
!       Purpose: Compute Legendre functions Qn(x) & Qn'(x)
!       Input :  x  --- Argument of Qn(x)
!                n  --- Degree of Qn(x)  ( n = 0,1,2,…)
!       Output:  QN(n) --- Qn(x)
!                QD(n) --- Qn'(x)
!       ====================================================
!
      IMPLICIT NONE
!*--LQNB2782
      DOUBLE PRECISION eps , q0 , q1 , qc1 , qc2 , Qd , qf , qf0 , qf1 ,&
                     & qf2 , Qn , qr , X , x2
      INTEGER j , k , l , N , nl
      DIMENSION Qn(0:N) , Qd(0:N)
      eps = 1.0D-14
      IF ( DABS(X)==1.0D0 ) THEN
         DO k = 0 , N
            Qn(k) = 1.0D+300
            Qd(k) = 1.0D+300
         ENDDO
         RETURN
      ENDIF
      IF ( X<=1.021D0 ) THEN
         x2 = DABS((1.0D0+X)/(1.0D0-X))
         q0 = 0.5D0*DLOG(x2)
         q1 = X*q0 - 1.0D0
         Qn(0) = q0
         Qn(1) = q1
         Qd(0) = 1.0D0/(1.0D0-X*X)
         Qd(1) = Qn(0) + X*Qd(0)
         DO k = 2 , N
            qf = ((2.0D0*k-1.0D0)*X*q1-(k-1.0D0)*q0)/k
            Qn(k) = qf
            Qd(k) = (Qn(k-1)-X*qf)*k/(1.0D0-X*X)
            q0 = q1
            q1 = qf
         ENDDO
      ELSE
         qc1 = 0.0D0
         qc2 = 1.0D0/X
         DO j = 1 , N
            qc2 = qc2*j/((2.0*j+1.0D0)*X)
            IF ( j==N-1 ) qc1 = qc2
         ENDDO
         DO l = 0 , 1
            nl = N + l
            qf = 1.0D0
            qr = 1.0D0
            DO k = 1 , 500
               qr = qr*(0.5D0*nl+k-1.0D0)*(0.5D0*(nl-1)+k)              &
                  & /((nl+k-0.5D0)*k*X*X)
               qf = qf + qr
               IF ( DABS(qr/qf)<eps ) GOTO 20
            ENDDO
 20         IF ( l==0 ) THEN
               Qn(N-1) = qf*qc1
            ELSE
               Qn(N) = qf*qc2
            ENDIF
         ENDDO
         qf2 = Qn(N)
         qf1 = Qn(N-1)
         DO k = N , 2 , -1
            qf0 = ((2*k-1.0D0)*X*qf1-k*qf2)/(k-1.0D0)
            Qn(k-2) = qf0
            qf2 = qf1
            qf1 = qf0
         ENDDO
         Qd(0) = 1.0D0/(1.0D0-X*X)
         DO k = 1 , N
            Qd(k) = k*(Qn(k-1)-X*Qn(k))/(1.0D0-X*X)
         ENDDO
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE CJK(Km,A)
!
!       ========================================================
!       Purpose: Compute the expansion coefficients for the
!                asymptotic expansion of Bessel functions
!                with large orders
!       Input :  Km   --- Maximum k
!       Output:  A(L) --- Cj(k) where j and k are related to L
!                         by L=j+1+[k*(k+1)]/2; j,k=0,1,...,Km
!       ========================================================
!
      IMPLICIT NONE
!*--CJK2865
      DOUBLE PRECISION A , f , f0 , g , g0
      INTEGER j , k , Km , l1 , l2 , l3 , l4
      DIMENSION A(*)
      A(1) = 1.0D0
      f0 = 1.0D0
      g0 = 1.0D0
      DO k = 0 , Km - 1
         l1 = (k+1)*(k+2)/2 + 1
         l2 = (k+1)*(k+2)/2 + k + 2
         f = (0.5D0*k+0.125D0/(k+1))*f0
         g = -(1.5D0*k+0.625D0/(3.0*(k+1.0D0)))*g0
         A(l1) = f
         A(l2) = g
         f0 = f
         g0 = g
      ENDDO
      DO k = 1 , Km - 1
         DO j = 1 , k
            l3 = k*(k+1)/2 + j + 1
            l4 = (k+1)*(k+2)/2 + j + 1
            A(l4) = (j+0.5D0*k+0.125D0/(2.0*j+k+1.0))*A(l3)             &
                  & - (j+0.5D0*k-1.0+0.625D0/(2.0*j+k+1.0))*A(l3-1)
         ENDDO
      ENDDO
      END
 
 
!       **********************************
 
      SUBROUTINE ITTIKA(X,Tti,Ttk)
!
!       =========================================================
!       Purpose: Integrate [I0(t)-1]/t with respect to t from 0
!                to x, and K0(t)/t with respect to t from x to ∞
!       Input :  x   --- Variable in the limits  ( x ≥ 0 )
!       Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
!                TTK --- Integration of K0(t)/t from x to ∞
!       =========================================================
!
      IMPLICIT NONE
!*--ITTIKA2909
      DOUBLE PRECISION b1 , c , e0 , el , pi , r , r2 , rc , rs , Tti , &
                     & Ttk , X
      INTEGER k
      DIMENSION c(8)
      pi = 3.141592653589793D0
      el = .5772156649015329D0
      DATA c/1.625D0 , 4.1328125D0 , 1.45380859375D+1 ,                 &
         & 6.553353881835D+1 , 3.6066157150269D+2 , 2.3448727161884D+3 ,&
         & 1.7588273098916D+4 , 1.4950639538279D+5/
      IF ( X==0.0D0 ) THEN
         Tti = 0.0D0
         Ttk = 1.0D+300
         RETURN
      ENDIF
      IF ( X<40.0D0 ) THEN
         Tti = 1.0D0
         r = 1.0D0
         DO k = 2 , 50
            r = .25D0*r*(k-1.0D0)/(k*k*k)*X*X
            Tti = Tti + r
            IF ( DABS(r/Tti)<1.0D-12 ) GOTO 50
         ENDDO
 50      Tti = Tti*.125D0*X*X
      ELSE
         Tti = 1.0D0
         r = 1.0D0
         DO k = 1 , 8
            r = r/X
            Tti = Tti + c(k)*r
         ENDDO
         rc = X*DSQRT(2.0D0*pi*X)
         Tti = Tti*EXP(X)/rc
      ENDIF
      IF ( X<=12.0D0 ) THEN
         e0 = (.5D0*DLOG(X/2.0D0)+el)*DLOG(X/2.0D0) + pi*pi/24.0D0 +    &
            & .5D0*el*el
         b1 = 1.5D0 - (el+DLOG(X/2.0D0))
         rs = 1.0D0
         r = 1.0D0
         DO k = 2 , 50
            r = .25D0*r*(k-1.0D0)/(k*k*k)*X*X
            rs = rs + 1.0D0/k
            r2 = r*(rs+1.0D0/(2.0D0*k)-(el+DLOG(X/2.0D0)))
            b1 = b1 + r2
            IF ( DABS(r2/b1)<1.0D-12 ) GOTO 100
         ENDDO
 100     Ttk = e0 - .125D0*X*X*b1
      ELSE
         Ttk = 1.0D0
         r = 1.0D0
         DO k = 1 , 8
            r = -r/X
            Ttk = Ttk + c(k)*r
         ENDDO
         rc = X*DSQRT(2.0D0/pi*X)
         Ttk = Ttk*EXP(-X)/rc
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE LAMV(V,X,Vm,Vl,Dl)
!
!       =========================================================
!       Purpose: Compute lambda function with arbitrary order v,
!                and their derivative
!       Input :  x --- Argument of lambda function
!                v --- Order of lambda function
!       Output:  VL(n) --- Lambda function of order n+v0
!                DL(n) --- Derivative of lambda function
!                VM --- Highest order computed
!       Routines called:
!            (1) MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!            (2) GAM0 for computing gamma function (|x| ≤ 1)
!       =========================================================
!
      IMPLICIT NONE
!*--LAMV2991
      DOUBLE PRECISION a0 , bjv0 , bjv1 , bk , ck , cs , Dl , f , f0 ,  &
                     & f1 , f2 , fac , ga , pi , px , qx , r , r0 , rc ,&
                     & rp
      DOUBLE PRECISION rp2 , rq , sk , uk , V , v0 , vk , Vl , Vm , vv ,&
                     & X , x2 , xk
      INTEGER i , j , k , k0 , m , MSTA1 , MSTA2 , n
      DIMENSION Vl(0:*) , Dl(0:*)
      pi = 3.141592653589793D0
      rp2 = 0.63661977236758D0
      X = DABS(X)
      x2 = X*X
      n = INT(V)
      v0 = V - n
      Vm = V
      IF ( X<=12.0D0 ) THEN
         DO k = 0 , n
            vk = v0 + k
            bk = 1.0D0
            r = 1.0D0
            DO i = 1 , 50
               r = -0.25D0*r*x2/(i*(i+vk))
               bk = bk + r
               IF ( DABS(r)<DABS(bk)*1.0D-15 ) GOTO 20
            ENDDO
 20         Vl(k) = bk
            uk = 1.0D0
            r = 1.0D0
            DO i = 1 , 50
               r = -0.25D0*r*x2/(i*(i+vk+1.0D0))
               uk = uk + r
               IF ( DABS(r)<DABS(uk)*1.0D-15 ) GOTO 40
            ENDDO
 40         Dl(k) = -0.5D0*X/(vk+1.0D0)*uk
         ENDDO
         RETURN
      ENDIF
      k0 = 11
      IF ( X>=35.0D0 ) k0 = 10
      IF ( X>=50.0D0 ) k0 = 8
      bjv0 = 0.0D0
      bjv1 = 0.0D0
      DO j = 0 , 1
         vv = 4.0D0*(j+v0)*(j+v0)
         px = 1.0D0
         rp = 1.0D0
         DO k = 1 , k0
            rp = -0.78125D-2*rp*(vv-(4.0*k-3.0)**2.0)                   &
               & *(vv-(4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*x2)
            px = px + rp
         ENDDO
         qx = 1.0D0
         rq = 1.0D0
         DO k = 1 , k0
            rq = -0.78125D-2*rq*(vv-(4.0*k-1.0)**2.0)                   &
               & *(vv-(4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*x2)
            qx = qx + rq
         ENDDO
         qx = 0.125D0*(vv-1.0D0)*qx/X
         xk = X - (0.5D0*(j+v0)+0.25D0)*pi
         a0 = DSQRT(rp2/X)
         ck = DCOS(xk)
         sk = DSIN(xk)
         IF ( j==0 ) bjv0 = a0*(px*ck-qx*sk)
         IF ( j==1 ) bjv1 = a0*(px*ck-qx*sk)
      ENDDO
      IF ( v0==0.0D0 ) THEN
         ga = 1.0D0
      ELSE
         CALL GAM0(v0,ga)
         ga = v0*ga
      ENDIF
      fac = (2.0D0/X)**v0*ga
      Vl(0) = bjv0
      Dl(0) = -bjv1 + v0/X*bjv0
      Vl(1) = bjv1
      Dl(1) = bjv0 - (1.0D0+v0)/X*bjv1
      r0 = 2.0D0*(1.0D0+v0)/X
      IF ( n<=1 ) THEN
         Vl(0) = fac*Vl(0)
         Dl(0) = fac*Dl(0) - v0/X*Vl(0)
         Vl(1) = fac*r0*Vl(1)
         Dl(1) = fac*r0*Dl(1) - (1.0D0+v0)/X*Vl(1)
         RETURN
      ENDIF
      IF ( n>=2 .AND. n<=INT(0.9*X) ) THEN
         f0 = bjv0
         f1 = bjv1
         DO k = 2 , n
            f = 2.0D0*(k+v0-1.0D0)/X*f1 - f0
            f0 = f1
            f1 = f
            Vl(k) = f
         ENDDO
      ELSEIF ( n>=2 ) THEN
         m = MSTA1(X,200)
         IF ( m<n ) THEN
            n = m
         ELSE
            m = MSTA2(X,n,15)
         ENDIF
         f = 0.0D0
         f2 = 0.0D0
         f1 = 1.0D-100
         DO k = m , 0 , -1
            f = 2.0D0*(v0+k+1.0D0)/X*f1 - f2
            IF ( k<=n ) Vl(k) = f
            f2 = f1
            f1 = f
         ENDDO
         cs = 0.0D0
         IF ( DABS(bjv0)>DABS(bjv1) ) THEN
            cs = bjv0/f
         ELSE
            cs = bjv1/f2
         ENDIF
         DO k = 0 , n
            Vl(k) = cs*Vl(k)
         ENDDO
      ENDIF
      Vl(0) = fac*Vl(0)
      DO j = 1 , n
         rc = fac*r0
         Vl(j) = rc*Vl(j)
         Dl(j-1) = -0.5D0*X/(j+v0)*Vl(j)
         r0 = 2.0D0*(j+v0+1)/X*r0
      ENDDO
      Dl(n) = 2.0D0*(v0+n)*(Vl(n-1)-Vl(n))/X
      Vm = n + v0
      END
 
 
 
!       **********************************
 
      SUBROUTINE CHGUIT(A,B,X,Hu,Id)
!
!       ======================================================
!       Purpose: Compute hypergeometric function U(a,b,x) by
!                using Gaussian-Legendre integration (n=60)
!       Input  : a  --- Parameter ( a > 0 )
!                b  --- Parameter
!                x  --- Argument ( x > 0 )
!       Output:  HU --- U(a,b,z)
!                ID --- Estimated number of significant digits
!       Routine called: GAMMA2 for computing Г(x)
!       ======================================================
!
      IMPLICIT NONE
!*--CHGUIT3143
      DOUBLE PRECISION A , a1 , B , b1 , c , d , f1 , f2 , g , ga , Hu ,&
                     & hu0 , hu1 , hu2 , s , t , t1 , t2 , t3 , t4
      DOUBLE PRECISION w , X
      INTEGER Id , j , k , m
      DIMENSION t(30) , w(30)
      DATA t/.259597723012478D-01 , .778093339495366D-01 ,              &
         & .129449135396945D+00 , .180739964873425D+00 ,                &
         & .231543551376029D+00 , .281722937423262D+00 ,                &
         & .331142848268448D+00 , .379670056576798D+00 ,                &
         & .427173741583078D+00 , .473525841761707D+00 ,                &
         & .518601400058570D+00 , .562278900753945D+00 ,                &
         & .604440597048510D+00 , .644972828489477D+00 ,                &
         & .683766327381356D+00 , .720716513355730D+00 ,                &
         & .755723775306586D+00 , .788693739932264D+00 ,                &
         & .819537526162146D+00 , .848171984785930D+00 ,                &
         & .874519922646898D+00 , .898510310810046D+00 ,                &
         & .920078476177628D+00 , .939166276116423D+00 ,                &
         & .955722255839996D+00 , .969701788765053D+00 ,                &
         & .981067201752598D+00 , .989787895222222D+00 ,                &
         & .995840525118838D+00 , .999210123227436D+00/
      DATA w/.519078776312206D-01 , .517679431749102D-01 ,              &
         & .514884515009810D-01 , .510701560698557D-01 ,                &
         & .505141845325094D-01 , .498220356905502D-01 ,                &
         & .489955754557568D-01 , .480370318199712D-01 ,                &
         & .469489888489122D-01 , .457343797161145D-01 ,                &
         & .443964787957872D-01 , .429388928359356D-01 ,                &
         & .413655512355848D-01 , .396806954523808D-01 ,                &
         & .378888675692434D-01 , .359948980510845D-01 ,                &
         & .340038927249464D-01 , .319212190192963D-01 ,                &
         & .297524915007890D-01 , .275035567499248D-01 ,                &
         & .251804776215213D-01 , .227895169439978D-01 ,                &
         & .203371207294572D-01 , .178299010142074D-01 ,                &
         & .152746185967848D-01 , .126781664768159D-01 ,                &
         & .100475571822880D-01 , .738993116334531D-02 ,                &
         & .471272992695363D-02 , .202681196887362D-02/
      Id = 9
!       DLMF 13.4.4, integration up to C=12/X
      a1 = A - 1.0D0
      b1 = B - A - 1.0D0
      c = 12.0D0/X
      hu0 = 0.0D0
      DO m = 10 , 100 , 5
         hu1 = 0.0D0
         g = 0.5D0*c/m
         d = g
         DO j = 1 , m
            s = 0.0D0
            DO k = 1 , 30
               t1 = d + g*t(k)
               t2 = d - g*t(k)
               f1 = EXP(-X*t1)*t1**a1*(1.0D0+t1)**b1
               f2 = EXP(-X*t2)*t2**a1*(1.0D0+t2)**b1
               s = s + w(k)*(f1+f2)
            ENDDO
            hu1 = hu1 + s*g
            d = d + 2.0D0*g
         ENDDO
         IF ( DABS(1.0D0-hu0/hu1)<1.0D-9 ) GOTO 100
         hu0 = hu1
      ENDDO
 100  CALL GAMMA2(A,ga)
      hu1 = hu1/ga
!       DLMF 13.4.4 with substitution t=C/(1-u)
!       integration u from 0 to 1, i.e. t from C=12/X to infinity
      DO m = 2 , 10 , 2
         hu2 = 0.0D0
         g = 0.5D0/m
         d = g
         DO j = 1 , m
            s = 0.0D0
            DO k = 1 , 30
               t1 = d + g*t(k)
               t2 = d - g*t(k)
               t3 = c/(1.0D0-t1)
               t4 = c/(1.0D0-t2)
               f1 = t3*t3/c*EXP(-X*t3)*t3**a1*(1.0D0+t3)**b1
               f2 = t4*t4/c*EXP(-X*t4)*t4**a1*(1.0D0+t4)**b1
               s = s + w(k)*(f1+f2)
            ENDDO
            hu2 = hu2 + s*g
            d = d + 2.0D0*g
         ENDDO
         IF ( DABS(1.0D0-hu0/hu2)<1.0D-9 ) GOTO 200
         hu0 = hu2
      ENDDO
 200  CALL GAMMA2(A,ga)
      hu2 = hu2/ga
      Hu = hu1 + hu2
      END
 
 
 
!       **********************************
 
      SUBROUTINE KMN(M,N,C,Cv,Kd,Df,Dn,Ck1,Ck2)
!
!       ===================================================
!       Purpose: Compute the expansion coefficients of the
!                prolate and oblate spheroidal functions
!                and joining factors
!       ===================================================
!
      IMPLICIT NONE
!*--KMN3250
      DOUBLE PRECISION C , Ck1 , Ck2 , cs , Cv , Df , Dn , dnp , g0 ,   &
                     & gk0 , gk1 , gk2 , gk3 , r , r1 , r2 , r3 , r4 ,  &
                     & r5 , rk
      DOUBLE PRECISION sa0 , sb0 , su0 , sw , t , tp , u , v , w
      INTEGER i , ip , j , k , Kd , l , M , N , nm , nm1 , nn
      DIMENSION u(200) , v(200) , w(200) , Df(200) , Dn(200) , tp(200) ,&
              & rk(200)
      nm = 25 + INT(0.5*(N-M)+C)
      nn = nm + M
      cs = C*C*Kd
      ip = 1
      IF ( N-M==2*INT((N-M)/2) ) ip = 0
      k = 0
      DO i = 1 , nn + 3
         IF ( ip==0 ) k = -2*(i-1)
         IF ( ip==1 ) k = -(2*i-3)
         gk0 = 2.0D0*M + k
         gk1 = (M+k)*(M+k+1.0D0)
         gk2 = 2.0D0*(M+k) - 1.0D0
         gk3 = 2.0D0*(M+k) + 3.0D0
         u(i) = gk0*(gk0-1.0D0)*cs/(gk2*(gk2+2.0D0))
         v(i) = gk1 - Cv + (2.0D0*(gk1-M*M)-1.0D0)*cs/(gk2*gk3)
         w(i) = (k+1.0D0)*(k+2.0D0)*cs/((gk2+2.0D0)*gk3)
      ENDDO
      DO k = 1 , M
         t = v(M+1)
         DO l = 0 , M - k - 1
            t = v(M-l) - w(M-l+1)*u(M-l)/t
         ENDDO
         rk(k) = -u(k)/t
      ENDDO
      r = 1.0D0
      DO k = 1 , M
         r = r*rk(k)
         Dn(k) = Df(1)*r
      ENDDO
      tp(nn) = v(nn+1)
      DO k = nn - 1 , M + 1 , -1
         tp(k) = v(k+1) - w(k+2)*u(k+1)/tp(k+1)
         IF ( k>M+1 ) rk(k) = -u(k)/tp(k)
      ENDDO
      IF ( M==0 ) dnp = Df(1)
      IF ( M/=0 ) dnp = Dn(M)
      Dn(M+1) = (-1)**ip*dnp*cs/((2.0*M-1.0)*(2.0*M+1.0-4.0*ip)*tp(M+1))
      DO k = M + 2 , nn
         Dn(k) = rk(k)*Dn(k-1)
      ENDDO
      r1 = 1.0D0
      DO j = 1 , (N+M+ip)/2
         r1 = r1*(j+0.5D0*(N+M+ip))
      ENDDO
      nm1 = (N-M)/2
      r = 1.0D0
      DO j = 1 , 2*M + ip
         r = r*j
      ENDDO
      su0 = r*Df(1)
      sw = 0.0D0
      DO k = 2 , nm
         r = r*(M+k-1.0)*(M+k+ip-1.5D0)/(k-1.0D0)/(k+ip-1.5D0)
         su0 = su0 + r*Df(k)
         IF ( k>nm1 .AND. DABS((su0-sw)/su0)<1.0D-14 ) GOTO 100
         sw = su0
      ENDDO
 100  IF ( Kd/=1 ) THEN
         r2 = 1.0D0
         DO j = 1 , M
            r2 = 2.0D0*C*r2*j
         ENDDO
         r3 = 1.0D0
         DO j = 1 , (N-M-ip)/2
            r3 = r3*j
         ENDDO
         sa0 = (2.0*(M+ip)+1.0)*r1/(2.0**N*C**ip*r2*r3*Df(1))
         Ck1 = sa0*su0
         IF ( Kd==-1 ) RETURN
      ENDIF
      r4 = 1.0D0
      DO j = 1 , (N-M-ip)/2
         r4 = 4.0D0*r4*j
      ENDDO
      r5 = 1.0D0
      DO j = 1 , M
         r5 = r5*(j+M)/C
      ENDDO
      g0 = Dn(M)
      IF ( M==0 ) g0 = Df(1)
      sb0 = (ip+1.0)*C**(ip+1)/(2.0*ip*(M-2.0)+1.0)/(2.0*M-1.0)
      Ck2 = (-1)**ip*sb0*r4*r5*g0/r1*su0
      END
 
 
 
!       **********************************
 
      SUBROUTINE LAGZO(N,X,W)
!
!       =========================================================
!       Purpose : Compute the zeros of Laguerre polynomial Ln(x)
!                 in the interval [0,∞], and the corresponding
!                 weighting coefficients for Gauss-Laguerre
!                 integration
!       Input :   n    --- Order of the Laguerre polynomial
!                 X(n) --- Zeros of the Laguerre polynomial
!                 W(n) --- Corresponding weighting coefficients
!       =========================================================
!
      IMPLICIT NONE
!*--LAGZO3362
      DOUBLE PRECISION f0 , f1 , fd , gd , hn , p , pd , pf , q , W ,   &
                     & wp , X , z , z0
      INTEGER i , it , j , k , N , nr
      DIMENSION X(N) , W(N)
      hn = 1.0D0/N
      pf = 0.0D0
      pd = 0.0D0
      DO nr = 1 , N
         z = hn
         IF ( nr>1 ) z = X(nr-1) + hn*nr**1.27
         it = 0
 50      it = it + 1
         z0 = z
         p = 1.0D0
         DO i = 1 , nr - 1
            p = p*(z-X(i))
         ENDDO
         f0 = 1.0D0
         f1 = 1.0D0 - z
         DO k = 2 , N
            pf = ((2.0D0*k-1.0D0-z)*f1-(k-1.0D0)*f0)/k
            pd = k/z*(pf-f1)
            f0 = f1
            f1 = pf
         ENDDO
         fd = pf/p
         q = 0.0D0
         DO i = 1 , nr - 1
            wp = 1.0D0
            DO j = 1 , nr - 1
               IF ( j/=i ) wp = wp*(z-X(j))
            ENDDO
            q = q + wp
         ENDDO
         gd = (pd-q*fd)/p
         z = z - fd/gd
         IF ( it<=40 .AND. DABS((z-z0)/z)>1.0D-15 ) GOTO 50
         X(nr) = z
         W(nr) = 1.0D0/(z*pd*pd)
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE VVLA(Va,X,Pv)
!
!       ===================================================
!       Purpose: Compute parabolic cylinder function Vv(x)
!                for large argument
!       Input:   x  --- Argument
!                va --- Order
!       Output:  PV --- Vv(x)
!       Routines called:
!             (1) DVLA for computing Dv(x) for large |x|
!             (2) GAMMA2 for computing Г(x)
!       ===================================================
!
      IMPLICIT NONE
!*--VVLA3424
      DOUBLE PRECISION a0 , dsl , eps , gl , pdl , pi , Pv , qe , r ,   &
                     & Va , X , x1
      INTEGER k
      pi = 3.141592653589793D0
      eps = 1.0D-12
      qe = EXP(0.25*X*X)
      a0 = DABS(X)**(-Va-1.0D0)*DSQRT(2.0D0/pi)*qe
      r = 1.0D0
      Pv = 1.0D0
      DO k = 1 , 18
         r = 0.5D0*r*(2.0*k+Va-1.0)*(2.0*k+Va)/(k*X*X)
         Pv = Pv + r
         IF ( DABS(r/Pv)<eps ) GOTO 100
      ENDDO
 100  Pv = a0*Pv
      IF ( X<0.0D0 ) THEN
         x1 = -X
         CALL DVLA(Va,x1,pdl)
         CALL GAMMA2(-Va,gl)
         dsl = DSIN(pi*Va)*DSIN(pi*Va)
         Pv = dsl*gl/pi*pdl - DCOS(pi*Va)*Pv
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE CJYVA(V,Z,Vm,Cbj,Cdj,Cby,Cdy)
!
!       ===========================================================
!       Purpose: Compute Bessel functions Jv(z), Yv(z) and their
!                derivatives for a complex argument
!       Input :  z --- Complex argument
!                v --- Order of Jv(z) and Yv(z)
!                      ( v = n+v0, n = 0,1,2,..., 0 ≤ v0 < 1 )
!       Output:  CBJ(n) --- Jn+v0(z)
!                CDJ(n) --- Jn+v0'(z)
!                CBY(n) --- Yn+v0(z)
!                CDY(n) --- Yn+v0'(z)
!                VM --- Highest order computed
!       Routines called:
!            (1) GAMMA2 for computing the gamma function
!            (2) MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ===========================================================
!
      IMPLICIT NONE
!*--CJYVA3476
      DOUBLE PRECISION a0 , ga , gb , pi , pv0 , pv1 , rp2 , V , v0 ,   &
                     & vg , vl , Vm , vv , w0 , w1 , wa , ya0 , ya1 ,   &
                     & yak
      COMPLEX*16 ca , ca0 , cb , Cbj , Cby , cck , Cdj , Cdy , cec ,    &
               & cf , cf0 , cf1 , cf2 , cfac0 , cfac1 , cg0 , cg1 ,     &
               & ch0 , ch1 , ch2
      COMPLEX*16 ci , cju0 , cju1 , cjv0 , cjv1 , cjvl , cp11 , cp12 ,  &
               & cp21 , cp22 , cpz , cqz , cr , cr0 , cr1 , crp , crq , &
               & cs , cs0 , cs1
      COMPLEX*16 csk , cyk , cyl1 , cyl2 , cylk , cyv0 , cyv1 , Z , z1 ,&
               & z2 , zk
      INTEGER j , k , k0 , l , lb , lb0 , m , MSTA1 , MSTA2 , n
      DIMENSION Cbj(0:*) , Cdj(0:*) , Cby(0:*) , Cdy(0:*)
      pi = 3.141592653589793D0
      rp2 = .63661977236758D0
      ci = (0.0D0,1.0D0)
      a0 = ABS(Z)
      z1 = Z
      z2 = Z*Z
      n = INT(V)
      v0 = V - n
      pv0 = pi*v0
      pv1 = pi*(1.0D0+v0)
      IF ( a0<1.0D-100 ) THEN
         DO k = 0 , n
            Cbj(k) = (0.0D0,0.0D0)
            Cdj(k) = (0.0D0,0.0D0)
            Cby(k) = -(1.0D+300,0.0D0)
            Cdy(k) = (1.0D+300,0.0D0)
         ENDDO
         IF ( v0==0.0 ) THEN
            Cbj(0) = (1.0D0,0.0D0)
            Cdj(1) = (0.5D0,0.0D0)
         ELSE
            Cdj(0) = (1.0D+300,0.0D0)
         ENDIF
         Vm = V
         RETURN
      ENDIF
      lb0 = 0.0D0
      IF ( DBLE(Z)<0.0 ) z1 = -Z
      IF ( a0<=12.0 ) THEN
         DO l = 0 , 1
            vl = v0 + l
            cjvl = (1.0D0,0.0D0)
            cr = (1.0D0,0.0D0)
            DO k = 1 , 40
               cr = -0.25D0*cr*z2/(k*(k+vl))
               cjvl = cjvl + cr
               IF ( ABS(cr)<ABS(cjvl)*1.0D-15 ) GOTO 20
            ENDDO
 20         vg = 1.0D0 + vl
            CALL GAMMA2(vg,ga)
            ca = (0.5D0*z1)**vl/ga
            IF ( l==0 ) cjv0 = cjvl*ca
            IF ( l==1 ) cjv1 = cjvl*ca
         ENDDO
      ELSE
         k0 = 11
         IF ( a0>=35.0 ) k0 = 10
         IF ( a0>=50.0 ) k0 = 8
         DO j = 0 , 1
            vv = 4.0D0*(j+v0)*(j+v0)
            cpz = (1.0D0,0.0D0)
            crp = (1.0D0,0.0D0)
            DO k = 1 , k0
               crp = -0.78125D-2*crp*(vv-(4.0*k-3.0)**2.0)              &
                   & *(vv-(4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*z2)
               cpz = cpz + crp
            ENDDO
            cqz = (1.0D0,0.0D0)
            crq = (1.0D0,0.0D0)
            DO k = 1 , k0
               crq = -0.78125D-2*crq*(vv-(4.0*k-1.0)**2.0)              &
                   & *(vv-(4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*z2)
               cqz = cqz + crq
            ENDDO
            cqz = 0.125D0*(vv-1.0)*cqz/z1
            zk = z1 - (0.5D0*(j+v0)+0.25D0)*pi
            ca0 = SQRT(rp2/z1)
            cck = COS(zk)
            csk = SIN(zk)
            IF ( j==0 ) THEN
               cjv0 = ca0*(cpz*cck-cqz*csk)
               cyv0 = ca0*(cpz*csk+cqz*cck)
            ELSEIF ( j==1 ) THEN
               cjv1 = ca0*(cpz*cck-cqz*csk)
               cyv1 = ca0*(cpz*csk+cqz*cck)
            ENDIF
         ENDDO
      ENDIF
      IF ( a0<=12.0 ) THEN
         IF ( v0/=0.0 ) THEN
            DO l = 0 , 1
               vl = v0 + l
               cjvl = (1.0D0,0.0D0)
               cr = (1.0D0,0.0D0)
               DO k = 1 , 40
                  cr = -0.25D0*cr*z2/(k*(k-vl))
                  cjvl = cjvl + cr
                  IF ( ABS(cr)<ABS(cjvl)*1.0D-15 ) GOTO 30
               ENDDO
 30            vg = 1.0D0 - vl
               CALL GAMMA2(vg,gb)
               cb = (2.0D0/z1)**vl/gb
               IF ( l==0 ) cju0 = cjvl*cb
               IF ( l==1 ) cju1 = cjvl*cb
            ENDDO
            cyv0 = (cjv0*DCOS(pv0)-cju0)/DSIN(pv0)
            cyv1 = (cjv1*DCOS(pv1)-cju1)/DSIN(pv1)
         ELSE
            cec = LOG(z1/2.0D0) + .5772156649015329D0
            cs0 = (0.0D0,0.0D0)
            w0 = 0.0D0
            cr0 = (1.0D0,0.0D0)
            DO k = 1 , 30
               w0 = w0 + 1.0D0/k
               cr0 = -0.25D0*cr0/(k*k)*z2
               cs0 = cs0 + cr0*w0
            ENDDO
            cyv0 = rp2*(cec*cjv0-cs0)
            cs1 = (1.0D0,0.0D0)
            w1 = 0.0D0
            cr1 = (1.0D0,0.0D0)
            DO k = 1 , 30
               w1 = w1 + 1.0D0/k
               cr1 = -0.25D0*cr1/(k*(k+1))*z2
               cs1 = cs1 + cr1*(2.0D0*w1+1.0D0/(k+1.0D0))
            ENDDO
            cyv1 = rp2*(cec*cjv1-1.0D0/z1-0.25D0*z1*cs1)
         ENDIF
      ENDIF
      IF ( DBLE(Z)<0.0D0 ) THEN
         cfac0 = EXP(pv0*ci)
         cfac1 = EXP(pv1*ci)
         IF ( DIMAG(Z)<0.0D0 ) THEN
            cyv0 = cfac0*cyv0 - 2.0D0*ci*DCOS(pv0)*cjv0
            cyv1 = cfac1*cyv1 - 2.0D0*ci*DCOS(pv1)*cjv1
            cjv0 = cjv0/cfac0
            cjv1 = cjv1/cfac1
         ELSEIF ( DIMAG(Z)>0.0D0 ) THEN
            cyv0 = cyv0/cfac0 + 2.0D0*ci*DCOS(pv0)*cjv0
            cyv1 = cyv1/cfac1 + 2.0D0*ci*DCOS(pv1)*cjv1
            cjv0 = cfac0*cjv0
            cjv1 = cfac1*cjv1
         ENDIF
      ENDIF
      Cbj(0) = cjv0
      Cbj(1) = cjv1
      IF ( n>=2 .AND. n<=INT(0.25*a0) ) THEN
         cf0 = cjv0
         cf1 = cjv1
         DO k = 2 , n
            cf = 2.0D0*(k+v0-1.0D0)/Z*cf1 - cf0
            Cbj(k) = cf
            cf0 = cf1
            cf1 = cf
         ENDDO
      ELSEIF ( n>=2 ) THEN
         m = MSTA1(a0,200)
         IF ( m<n ) THEN
            n = m
         ELSE
            m = MSTA2(a0,n,15)
         ENDIF
         cf2 = (0.0D0,0.0D0)
         cf1 = (1.0D-100,0.0D0)
         DO k = m , 0 , -1
            cf = 2.0D0*(v0+k+1.0D0)/Z*cf1 - cf2
            IF ( k<=n ) Cbj(k) = cf
            cf2 = cf1
            cf1 = cf
         ENDDO
         IF ( ABS(cjv0)>ABS(cjv1) ) cs = cjv0/cf
         IF ( ABS(cjv0)<=ABS(cjv1) ) cs = cjv1/cf2
         DO k = 0 , n
            Cbj(k) = cs*Cbj(k)
         ENDDO
      ENDIF
      Cdj(0) = v0/Z*Cbj(0) - Cbj(1)
      DO k = 1 , n
         Cdj(k) = -(k+v0)/Z*Cbj(k) + Cbj(k-1)
      ENDDO
      Cby(0) = cyv0
      Cby(1) = cyv1
      ya0 = ABS(cyv0)
      lb = 0
      cg0 = cyv0
      cg1 = cyv1
      DO k = 2 , n
         cyk = 2.0D0*(v0+k-1.0D0)/Z*cg1 - cg0
         IF ( ABS(cyk)<=1.0D+290 ) THEN
            yak = ABS(cyk)
            ya1 = ABS(cg0)
            IF ( yak<ya0 .AND. yak<ya1 ) lb = k
            Cby(k) = cyk
            cg0 = cg1
            cg1 = cyk
         ENDIF
      ENDDO
      IF ( lb>4 .AND. DIMAG(Z)/=0.0D0 ) THEN
 50      IF ( lb/=lb0 ) THEN
            ch2 = (1.0D0,0.0D0)
            ch1 = (0.0D0,0.0D0)
            lb0 = lb
            DO k = lb , 1 , -1
               ch0 = 2.0D0*(k+v0)/Z*ch1 - ch2
               ch2 = ch1
               ch1 = ch0
            ENDDO
            cp12 = ch0
            cp22 = ch2
            ch2 = (0.0D0,0.0D0)
            ch1 = (1.0D0,0.0D0)
            DO k = lb , 1 , -1
               ch0 = 2.0D0*(k+v0)/Z*ch1 - ch2
               ch2 = ch1
               ch1 = ch0
            ENDDO
            cp11 = ch0
            cp21 = ch2
            IF ( lb==n ) Cbj(lb+1) = 2.0D0*(lb+v0)/Z*Cbj(lb) - Cbj(lb-1)
            IF ( ABS(Cbj(0))>ABS(Cbj(1)) ) THEN
               Cby(lb+1) = (Cbj(lb+1)*cyv0-2.0D0*cp11/(pi*Z))/Cbj(0)
               Cby(lb) = (Cbj(lb)*cyv0+2.0D0*cp12/(pi*Z))/Cbj(0)
            ELSE
               Cby(lb+1) = (Cbj(lb+1)*cyv1-2.0D0*cp21/(pi*Z))/Cbj(1)
               Cby(lb) = (Cbj(lb)*cyv1+2.0D0*cp22/(pi*Z))/Cbj(1)
            ENDIF
            cyl2 = Cby(lb+1)
            cyl1 = Cby(lb)
            DO k = lb - 1 , 0 , -1
               cylk = 2.0D0*(k+v0+1.0D0)/Z*cyl1 - cyl2
               Cby(k) = cylk
               cyl2 = cyl1
               cyl1 = cylk
            ENDDO
            cyl1 = Cby(lb)
            cyl2 = Cby(lb+1)
            DO k = lb + 1 , n - 1
               cylk = 2.0D0*(k+v0)/Z*cyl2 - cyl1
               Cby(k+1) = cylk
               cyl1 = cyl2
               cyl2 = cylk
            ENDDO
            DO k = 2 , n
               wa = ABS(Cby(k))
               IF ( wa<ABS(Cby(k-1)) ) lb = k
            ENDDO
            GOTO 50
         ENDIF
      ENDIF
      Cdy(0) = v0/Z*Cby(0) - Cby(1)
      DO k = 1 , n
         Cdy(k) = Cby(k-1) - (k+v0)/Z*Cby(k)
      ENDDO
      Vm = n + v0
      END
 
 
 
!       **********************************
 
      SUBROUTINE CJYVB(V,Z,Vm,Cbj,Cdj,Cby,Cdy)
!
!       ===========================================================
!       Purpose: Compute Bessel functions Jv(z), Yv(z) and their
!                derivatives for a complex argument
!       Input :  z --- Complex argument
!                v --- Order of Jv(z) and Yv(z)
!                      ( v = n+v0, n = 0,1,2,..., 0 ≤ v0 < 1 )
!       Output:  CBJ(n) --- Jn+v0(z)
!                CDJ(n) --- Jn+v0'(z)
!                CBY(n) --- Yn+v0(z)
!                CDY(n) --- Yn+v0'(z)
!                VM --- Highest order computed
!       Routines called:
!            (1) GAMMA2 for computing the gamma function
!            (2) MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ===========================================================
!
      IMPLICIT NONE
!*--CJYVB3763
      DOUBLE PRECISION a0 , ga , gb , pi , pv0 , rp2 , V , v0 , vg ,    &
                     & Vm , vv , w0
      COMPLEX*16 ca , ca0 , cb , Cbj , Cby , cck , Cdj , Cdy , cec ,    &
               & cf , cf1 , cf2 , cfac0 , ci , cju0 , cjv0 , cjvn ,     &
               & cpz , cqz , cr
      COMPLEX*16 cr0 , crp , crq , cs , cs0 , csk , cyv0 , cyy , Z ,    &
               & z1 , z2 , zk
      INTEGER k , k0 , m , MSTA1 , MSTA2 , n
      DIMENSION Cbj(0:*) , Cdj(0:*) , Cby(0:*) , Cdy(0:*)
      pi = 3.141592653589793D0
      rp2 = .63661977236758D0
      ci = (0.0D0,1.0D0)
      a0 = ABS(Z)
      z1 = Z
      z2 = Z*Z
      n = INT(V)
      v0 = V - n
      pv0 = pi*v0
      IF ( a0<1.0D-100 ) THEN
         DO k = 0 , n
            Cbj(k) = (0.0D0,0.0D0)
            Cdj(k) = (0.0D0,0.0D0)
            Cby(k) = -(1.0D+300,0.0D0)
            Cdy(k) = (1.0D+300,0.0D0)
         ENDDO
         IF ( v0==0.0 ) THEN
            Cbj(0) = (1.0D0,0.0D0)
            Cdj(1) = (0.5D0,0.0D0)
         ELSE
            Cdj(0) = (1.0D+300,0.0D0)
         ENDIF
         Vm = V
         RETURN
      ENDIF
      IF ( DBLE(Z)<0.0D0 ) z1 = -Z
      IF ( a0<=12.0 ) THEN
         cjv0 = (1.0D0,0.0D0)
         cr = (1.0D0,0.0D0)
         DO k = 1 , 40
            cr = -0.25D0*cr*z2/(k*(k+v0))
            cjv0 = cjv0 + cr
            IF ( ABS(cr)<ABS(cjv0)*1.0D-15 ) GOTO 50
         ENDDO
 50      vg = 1.0D0 + v0
         CALL GAMMA2(vg,ga)
         ca = (0.5D0*z1)**v0/ga
         cjv0 = cjv0*ca
      ELSE
         k0 = 11
         IF ( a0>=35.0 ) k0 = 10
         IF ( a0>=50.0 ) k0 = 8
         vv = 4.0D0*v0*v0
         cpz = (1.0D0,0.0D0)
         crp = (1.0D0,0.0D0)
         DO k = 1 , k0
            crp = -0.78125D-2*crp*(vv-(4.0*k-3.0)**2.0)                 &
                & *(vv-(4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*z2)
            cpz = cpz + crp
         ENDDO
         cqz = (1.0D0,0.0D0)
         crq = (1.0D0,0.0D0)
         DO k = 1 , k0
            crq = -0.78125D-2*crq*(vv-(4.0*k-1.0)**2.0)                 &
                & *(vv-(4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*z2)
            cqz = cqz + crq
         ENDDO
         cqz = 0.125D0*(vv-1.0)*cqz/z1
         zk = z1 - (0.5D0*v0+0.25D0)*pi
         ca0 = SQRT(rp2/z1)
         cck = COS(zk)
         csk = SIN(zk)
         cjv0 = ca0*(cpz*cck-cqz*csk)
         cyv0 = ca0*(cpz*csk+cqz*cck)
      ENDIF
      IF ( a0<=12.0 ) THEN
         IF ( v0/=0.0 ) THEN
            cjvn = (1.0D0,0.0D0)
            cr = (1.0D0,0.0D0)
            DO k = 1 , 40
               cr = -0.25D0*cr*z2/(k*(k-v0))
               cjvn = cjvn + cr
               IF ( ABS(cr)<ABS(cjvn)*1.0D-15 ) GOTO 60
            ENDDO
 60         vg = 1.0D0 - v0
            CALL GAMMA2(vg,gb)
            cb = (2.0D0/z1)**v0/gb
            cju0 = cjvn*cb
            cyv0 = (cjv0*DCOS(pv0)-cju0)/DSIN(pv0)
         ELSE
            cec = LOG(z1/2.0D0) + .5772156649015329D0
            cs0 = (0.0D0,0.0D0)
            w0 = 0.0D0
            cr0 = (1.0D0,0.0D0)
            DO k = 1 , 30
               w0 = w0 + 1.0D0/k
               cr0 = -0.25D0*cr0/(k*k)*z2
               cs0 = cs0 + cr0*w0
            ENDDO
            cyv0 = rp2*(cec*cjv0-cs0)
         ENDIF
      ENDIF
      IF ( n==0 ) n = 1
      m = MSTA1(a0,200)
      IF ( m<n ) THEN
         n = m
      ELSE
         m = MSTA2(a0,n,15)
      ENDIF
      cf2 = (0.0D0,0.0D0)
      cf1 = (1.0D-100,0.0D0)
      DO k = m , 0 , -1
         cf = 2.0D0*(v0+k+1.0D0)/z1*cf1 - cf2
         IF ( k<=n ) Cbj(k) = cf
         cf2 = cf1
         cf1 = cf
      ENDDO
      cs = cjv0/cf
      DO k = 0 , n
         Cbj(k) = cs*Cbj(k)
      ENDDO
      IF ( DBLE(Z)<0.0D0 ) THEN
         cfac0 = EXP(pv0*ci)
         IF ( DIMAG(Z)<0.0D0 ) THEN
            cyv0 = cfac0*cyv0 - 2.0D0*ci*DCOS(pv0)*cjv0
         ELSEIF ( DIMAG(Z)>0.0D0 ) THEN
            cyv0 = cyv0/cfac0 + 2.0D0*ci*DCOS(pv0)*cjv0
         ENDIF
         DO k = 0 , n
            IF ( DIMAG(Z)<0.0D0 ) THEN
               Cbj(k) = EXP(-pi*(k+v0)*ci)*Cbj(k)
            ELSEIF ( DIMAG(Z)>0.0D0 ) THEN
               Cbj(k) = EXP(pi*(k+v0)*ci)*Cbj(k)
            ENDIF
         ENDDO
         z1 = z1
      ENDIF
      Cby(0) = cyv0
      DO k = 1 , n
         cyy = (Cbj(k)*Cby(k-1)-2.0D0/(pi*Z))/Cbj(k-1)
         Cby(k) = cyy
      ENDDO
      Cdj(0) = v0/Z*Cbj(0) - Cbj(1)
      DO k = 1 , n
         Cdj(k) = -(k+v0)/Z*Cbj(k) + Cbj(k-1)
      ENDDO
      Cdy(0) = v0/Z*Cby(0) - Cby(1)
      DO k = 1 , n
         Cdy(k) = Cby(k-1) - (k+v0)/Z*Cby(k)
      ENDDO
      Vm = n + v0
      END
 
 
 
!       **********************************
 
      SUBROUTINE JY01A(X,Bj0,Dj0,Bj1,Dj1,By0,Dy0,By1,Dy1)
!
!       =======================================================
!       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
!                Y1(x), and their derivatives
!       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ≥ 0 )
!       Output:  BJ0 --- J0(x)
!                DJ0 --- J0'(x)
!                BJ1 --- J1(x)
!                DJ1 --- J1'(x)
!                BY0 --- Y0(x)
!                DY0 --- Y0'(x)
!                BY1 --- Y1(x)
!                DY1 --- Y1'(x)
!       =======================================================
!
      IMPLICIT NONE
!*--JY01A3940
      DOUBLE PRECISION a , a1 , b , b1 , Bj0 , Bj1 , By0 , By1 , cs0 ,  &
                     & cs1 , cu , Dj0 , Dj1 , Dy0 , Dy1 , ec , p0 , p1 ,&
                     & pi , q0
      DOUBLE PRECISION q1 , r , r0 , r1 , rp2 , t1 , t2 , w0 , w1 , X , &
                     & x2
      INTEGER k , k0
      DIMENSION a(12) , b(12) , a1(12) , b1(12)
      pi = 3.141592653589793D0
      rp2 = 0.63661977236758D0
      x2 = X*X
      IF ( X==0.0D0 ) THEN
         Bj0 = 1.0D0
         Bj1 = 0.0D0
         Dj0 = 0.0D0
         Dj1 = 0.5D0
         By0 = -1.0D+300
         By1 = -1.0D+300
         Dy0 = 1.0D+300
         Dy1 = 1.0D+300
         RETURN
      ENDIF
      IF ( X<=12.0D0 ) THEN
         Bj0 = 1.0D0
         r = 1.0D0
         DO k = 1 , 30
            r = -0.25D0*r*x2/(k*k)
            Bj0 = Bj0 + r
            IF ( DABS(r)<DABS(Bj0)*1.0D-15 ) GOTO 50
         ENDDO
 50      Bj1 = 1.0D0
         r = 1.0D0
         DO k = 1 , 30
            r = -0.25D0*r*x2/(k*(k+1.0D0))
            Bj1 = Bj1 + r
            IF ( DABS(r)<DABS(Bj1)*1.0D-15 ) GOTO 100
         ENDDO
 100     Bj1 = 0.5D0*X*Bj1
         ec = DLOG(X/2.0D0) + 0.5772156649015329D0
         cs0 = 0.0D0
         w0 = 0.0D0
         r0 = 1.0D0
         DO k = 1 , 30
            w0 = w0 + 1.0D0/k
            r0 = -0.25D0*r0/(k*k)*x2
            r = r0*w0
            cs0 = cs0 + r
            IF ( DABS(r)<DABS(cs0)*1.0D-15 ) GOTO 150
         ENDDO
 150     By0 = rp2*(ec*Bj0-cs0)
         cs1 = 1.0D0
         w1 = 0.0D0
         r1 = 1.0D0
         DO k = 1 , 30
            w1 = w1 + 1.0D0/k
            r1 = -0.25D0*r1/(k*(k+1))*x2
            r = r1*(2.0D0*w1+1.0D0/(k+1.0D0))
            cs1 = cs1 + r
            IF ( DABS(r)<DABS(cs1)*1.0D-15 ) GOTO 200
         ENDDO
 200     By1 = rp2*(ec*Bj1-1.0D0/X-0.25D0*X*cs1)
      ELSE
         DATA a/ - .7031250000000000D-01 , .1121520996093750D+00 ,      &
            & -.5725014209747314D+00 , .6074042001273483D+01 ,          &
            & -.1100171402692467D+03 , .3038090510922384D+04 ,          &
            & -.1188384262567832D+06 , .6252951493434797D+07 ,          &
            & -.4259392165047669D+09 , .3646840080706556D+11 ,          &
            & -.3833534661393944D+13 , .4854014686852901D+15/
         DATA b/.7324218750000000D-01 , -.2271080017089844D+00 ,        &
            & .1727727502584457D+01 , -.2438052969955606D+02 ,          &
            & .5513358961220206D+03 , -.1825775547429318D+05 ,          &
            & .8328593040162893D+06 , -.5006958953198893D+08 ,          &
            & .3836255180230433D+10 , -.3649010818849833D+12 ,          &
            & .4218971570284096D+14 , -.5827244631566907D+16/
         DATA a1/.1171875000000000D+00 , -.1441955566406250D+00 ,       &
            & .6765925884246826D+00 , -.6883914268109947D+01 ,          &
            & .1215978918765359D+03 , -.3302272294480852D+04 ,          &
            & .1276412726461746D+06 , -.6656367718817688D+07 ,          &
            & .4502786003050393D+09 , -.3833857520742790D+11 ,          &
            & .4011838599133198D+13 , -.5060568503314727D+15/
         DATA b1/ - .1025390625000000D+00 , .2775764465332031D+00 ,     &
            & -.1993531733751297D+01 , .2724882731126854D+02 ,          &
            & -.6038440767050702D+03 , .1971837591223663D+05 ,          &
            & -.8902978767070678D+06 , .5310411010968522D+08 ,          &
            & -.4043620325107754D+10 , .3827011346598605D+12 ,          &
            & -.4406481417852278D+14 , .6065091351222699D+16/
         k0 = 12
         IF ( X>=35.0 ) k0 = 10
         IF ( X>=50.0 ) k0 = 8
         t1 = X - 0.25D0*pi
         p0 = 1.0D0
         q0 = -0.125D0/X
         DO k = 1 , k0
            p0 = p0 + a(k)*X**(-2*k)
            q0 = q0 + b(k)*X**(-2*k-1)
         ENDDO
         cu = DSQRT(rp2/X)
         Bj0 = cu*(p0*DCOS(t1)-q0*DSIN(t1))
         By0 = cu*(p0*DSIN(t1)+q0*DCOS(t1))
         t2 = X - 0.75D0*pi
         p1 = 1.0D0
         q1 = 0.375D0/X
         DO k = 1 , k0
            p1 = p1 + a1(k)*X**(-2*k)
            q1 = q1 + b1(k)*X**(-2*k-1)
         ENDDO
         cu = DSQRT(rp2/X)
         Bj1 = cu*(p1*DCOS(t2)-q1*DSIN(t2))
         By1 = cu*(p1*DSIN(t2)+q1*DCOS(t2))
      ENDIF
      Dj0 = -Bj1
      Dj1 = Bj0 - Bj1/X
      Dy0 = -By1
      Dy1 = By0 - By1/X
      END
 
!       **********************************
 
      SUBROUTINE INCOG(A,X,Gin,Gim,Gip,Isfer)
!
!       ===================================================
!       Purpose: Compute the incomplete gamma function
!                r(a,x), Г(a,x) and P(a,x)
!       Input :  a   --- Parameter ( a ≤ 170 )
!                x   --- Argument
!       Output:  GIN --- r(a,x)
!                GIM --- Г(a,x)
!                GIP --- P(a,x)
!                ISFER --- Error flag
!       Routine called: GAMMA2 for computing Г(x)
!       ===================================================
!
      IMPLICIT NONE
!*--INCOG4076
      DOUBLE PRECISION A , ga , Gim , Gin , Gip , r , s , t0 , X , xam
      INTEGER Isfer , k
      Isfer = 0
      xam = -X + A*DLOG(X)
      IF ( xam>700.0 .OR. A>170.0 ) THEN
         Isfer = 6
         RETURN
      ENDIF
      IF ( X==0.0 ) THEN
         Gin = 0.0
         CALL GAMMA2(A,ga)
         Gim = ga
         Gip = 0.0
      ELSEIF ( X<=1.0+A ) THEN
         s = 1.0D0/A
         r = s
         DO k = 1 , 60
            r = r*X/(A+k)
            s = s + r
            IF ( DABS(r/s)<1.0D-15 ) GOTO 50
         ENDDO
 50      Gin = EXP(xam)*s
         CALL GAMMA2(A,ga)
         Gip = Gin/ga
         Gim = ga - Gin
      ELSEIF ( X>1.0+A ) THEN
         t0 = 0.0D0
         DO k = 60 , 1 , -1
            t0 = (k-A)/(1.0D0+k/(X+t0))
         ENDDO
         Gim = EXP(xam)/(X+t0)
         CALL GAMMA2(A,ga)
         Gin = ga - Gim
         Gip = 1.0D0 - Gim/ga
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE ITIKB(X,Ti,Tk)
!
!       =======================================================
!       Purpose: Integrate Bessel functions I0(t) and K0(t)
!                with respect to t from 0 to x
!       Input :  x  --- Upper limit of the integral ( x ≥ 0 )
!       Output:  TI --- Integration of I0(t) from 0 to x
!                TK --- Integration of K0(t) from 0 to x
!       =======================================================
!
      IMPLICIT NONE
!*--ITIKB4132
      DOUBLE PRECISION pi , t , t1 , Ti , Tk , X
      pi = 3.141592653589793D0
      IF ( X==0.0D0 ) THEN
         Ti = 0.0D0
      ELSEIF ( X<5.0D0 ) THEN
         t1 = X/5.0D0
         t = t1*t1
         Ti = ((((((((.59434D-3*t+.4500642D-2)*t+.044686921D0)*t+       &
            & .300704878D0)*t+1.471860153D0)*t+4.844024624D0)           &
            & *t+9.765629849D0)*t+10.416666367D0)*t+5.0D0)*t1
      ELSEIF ( X>=5.0 .AND. X<=8.0D0 ) THEN
         t = 5.0D0/X
         Ti = (((-.015166D0*t-.0202292D0)*t+.1294122D0)*t-.0302912D0)   &
            & *t + .4161224D0
         Ti = Ti*EXP(X)/DSQRT(X)
      ELSE
         t = 8.0D0/X
         Ti = (((((-.0073995D0*t+.017744D0)*t-.0114858D0)*t+.55956D-2)  &
            & *t+.59191D-2)*t+.0311734D0)*t + .3989423D0
         Ti = Ti*EXP(X)/DSQRT(X)
      ENDIF
      IF ( X==0.0D0 ) THEN
         Tk = 0.0D0
      ELSEIF ( X<=2.0D0 ) THEN
         t1 = X/2.0D0
         t = t1*t1
         Tk = ((((((.116D-5*t+.2069D-4)*t+.62664D-3)*t+.01110118D0)*t+  &
            & .11227902D0)*t+.50407836D0)*t+.84556868D0)*t1
         Tk = Tk - DLOG(X/2.0D0)*Ti
      ELSEIF ( X>2.0 .AND. X<=4.0D0 ) THEN
         t = 2.0D0/X
         Tk = (((.0160395D0*t-.0781715D0)*t+.185984D0)*t-.3584641D0)    &
            & *t + 1.2494934D0
         Tk = pi/2.0D0 - Tk*EXP(-X)/DSQRT(X)
      ELSEIF ( X>4.0 .AND. X<=7.0D0 ) THEN
         t = 4.0D0/X
         Tk = (((((.37128D-2*t-.0158449D0)*t+.0320504D0)*t-.0481455D0)  &
            & *t+.0787284D0)*t-.1958273D0)*t + 1.2533141D0
         Tk = pi/2.0D0 - Tk*EXP(-X)/DSQRT(X)
      ELSE
         t = 7.0D0/X
         Tk = (((((.33934D-3*t-.163271D-2)*t+.417454D-2)*t-.933944D-2)  &
            & *t+.02576646D0)*t-.11190289D0)*t + 1.25331414D0
         Tk = pi/2.0D0 - Tk*EXP(-X)/DSQRT(X)
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE ITIKA(X,Ti,Tk)
!
!       =======================================================
!       Purpose: Integrate modified Bessel functions I0(t) and
!                K0(t) with respect to t from 0 to x
!       Input :  x  --- Upper limit of the integral  ( x ≥ 0 )
!       Output:  TI --- Integration of I0(t) from 0 to x
!                TK --- Integration of K0(t) from 0 to x
!       =======================================================
!
      IMPLICIT NONE
!*--ITIKA4196
      DOUBLE PRECISION a , b1 , b2 , e0 , el , pi , r , rc1 , rc2 , rs ,&
                     & Ti , Tk , tw , X , x2
      INTEGER k
      DIMENSION a(10)
      pi = 3.141592653589793D0
      el = .5772156649015329D0
      DATA a/.625D0 , 1.0078125D0 , 2.5927734375D0 , 9.1868591308594D0 ,&
         & 4.1567974090576D+1 , 2.2919635891914D+2 , 1.491504060477D+3 ,&
         & 1.1192354495579D+4 , 9.515939374212D+4 , 9.0412425769041D+5/
      IF ( X==0.0D0 ) THEN
         Ti = 0.0D0
         Tk = 0.0D0
         RETURN
      ELSEIF ( X<20.0D0 ) THEN
         x2 = X*X
         Ti = 1.0D0
         r = 1.0D0
         DO k = 1 , 50
            r = .25D0*r*(2*k-1.0D0)/(2*k+1.0D0)/(k*k)*x2
            Ti = Ti + r
            IF ( DABS(r/Ti)<1.0D-12 ) GOTO 50
         ENDDO
 50      Ti = Ti*X
      ELSE
         x2 = 0.0D0
         Ti = 1.0D0
         r = 1.0D0
         DO k = 1 , 10
            r = r/X
            Ti = Ti + a(k)*r
         ENDDO
         rc1 = 1.0D0/DSQRT(2.0D0*pi*X)
         Ti = rc1*EXP(X)*Ti
      ENDIF
      IF ( X<12.0D0 ) THEN
         e0 = el + DLOG(X/2.0D0)
         b1 = 1.0D0 - e0
         b2 = 0.0D0
         rs = 0.0D0
         r = 1.0D0
         tw = 0.0D0
         DO k = 1 , 50
            r = .25D0*r*(2*k-1.0D0)/(2*k+1.0D0)/(k*k)*x2
            b1 = b1 + r*(1.0D0/(2*k+1)-e0)
            rs = rs + 1.0D0/k
            b2 = b2 + r*rs
            Tk = b1 + b2
            IF ( DABS((Tk-tw)/Tk)<1.0D-12 ) GOTO 100
            tw = Tk
         ENDDO
 100     Tk = Tk*X
      ELSE
         Tk = 1.0D0
         r = 1.0D0
         DO k = 1 , 10
            r = -r/X
            Tk = Tk + a(k)*r
         ENDDO
         rc2 = DSQRT(pi/(2.0D0*X))
         Tk = pi/2.0D0 - rc2*Tk*EXP(-X)
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE JYV(V,X,Vm,Bj,Dj,By,Dy)
!
!       =======================================================
!       Purpose: Compute Bessel functions Jv(x) and Yv(x)
!                and their derivatives
!       Input :  x --- Argument of Jv(x) and Yv(x)
!                v --- Order of Jv(x) and Yv(x)
!                      ( v = n+v0, 0 ≤ v0 < 1, n = 0,1,2,... )
!       Output:  BJ(n) --- Jn+v0(x)
!                DJ(n) --- Jn+v0'(x)
!                BY(n) --- Yn+v0(x)
!                DY(n) --- Yn+v0'(x)
!                VM --- Highest order computed
!       Routines called:
!            (1) GAMMA2 for computing gamma function
!            (2) MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       =======================================================
!
      IMPLICIT NONE
!*--JYV4285
      DOUBLE PRECISION a , a0 , b , Bj , bju0 , bju1 , bjv0 , bjv1 ,    &
                     & bjvl , By , byv0 , byv1 , byvk , ck , cs , cs0 , &
                     & cs1 , Dj , Dy , ec
      DOUBLE PRECISION el , f , f0 , f1 , f2 , ga , gb , pi , pv0 ,     &
                     & pv1 , px , qx , r , r0 , r1 , rp , rp2 , rq ,    &
                     & sk , V
      DOUBLE PRECISION v0 , vg , vl , Vm , vv , w0 , w1 , X , x2 , xk
      INTEGER j , k , k0 , l , m , MSTA1 , MSTA2 , n
      DIMENSION Bj(0:*) , Dj(0:*) , By(0:*) , Dy(0:*)
      el = .5772156649015329D0
      pi = 3.141592653589793D0
      rp2 = .63661977236758D0
      x2 = X*X
      n = INT(V)
      v0 = V - n
      IF ( X<1.0D-100 ) THEN
         DO k = 0 , n
            Bj(k) = 0.0D0
            Dj(k) = 0.0D0
            By(k) = -1.0D+300
            Dy(k) = 1.0D+300
         ENDDO
         IF ( v0==0.0 ) THEN
            Bj(0) = 1.0D0
            Dj(1) = 0.5D0
         ELSE
            Dj(0) = 1.0D+300
         ENDIF
         Vm = V
         RETURN
      ENDIF
      bjv0 = 0.0D0
      bjv1 = 0.0D0
      byv0 = 0.0D0
      byv1 = 0.0D0
      IF ( X<=12.0 ) THEN
         DO l = 0 , 1
            vl = v0 + l
            bjvl = 1.0D0
            r = 1.0D0
            DO k = 1 , 40
               r = -0.25D0*r*x2/(k*(k+vl))
               bjvl = bjvl + r
               IF ( DABS(r)<DABS(bjvl)*1.0D-15 ) GOTO 20
            ENDDO
 20         vg = 1.0D0 + vl
            CALL GAMMA2(vg,ga)
            a = (0.5D0*X)**vl/ga
            IF ( l==0 ) bjv0 = bjvl*a
            IF ( l==1 ) bjv1 = bjvl*a
         ENDDO
      ELSE
         k0 = 11
         IF ( X>=35.0 ) k0 = 10
         IF ( X>=50.0 ) k0 = 8
         DO j = 0 , 1
            vv = 4.0D0*(j+v0)*(j+v0)
            px = 1.0D0
            rp = 1.0D0
            DO k = 1 , k0
               rp = -0.78125D-2*rp*(vv-(4.0*k-3.0)**2.0)                &
                  & *(vv-(4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*x2)
               px = px + rp
            ENDDO
            qx = 1.0D0
            rq = 1.0D0
            DO k = 1 , k0
               rq = -0.78125D-2*rq*(vv-(4.0*k-1.0)**2.0)                &
                  & *(vv-(4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*x2)
               qx = qx + rq
            ENDDO
            qx = 0.125D0*(vv-1.0)*qx/X
            xk = X - (0.5D0*(j+v0)+0.25D0)*pi
            a0 = DSQRT(rp2/X)
            ck = DCOS(xk)
            sk = DSIN(xk)
            IF ( j==0 ) THEN
               bjv0 = a0*(px*ck-qx*sk)
               byv0 = a0*(px*sk+qx*ck)
            ELSEIF ( j==1 ) THEN
               bjv1 = a0*(px*ck-qx*sk)
               byv1 = a0*(px*sk+qx*ck)
            ENDIF
         ENDDO
      ENDIF
      Bj(0) = bjv0
      Bj(1) = bjv1
      Dj(0) = v0/X*Bj(0) - Bj(1)
      Dj(1) = -(1.0D0+v0)/X*Bj(1) + Bj(0)
      IF ( n>=2 .AND. n<=INT(0.9*X) ) THEN
         f0 = bjv0
         f1 = bjv1
         DO k = 2 , n
            f = 2.0D0*(k+v0-1.0D0)/X*f1 - f0
            Bj(k) = f
            f0 = f1
            f1 = f
         ENDDO
      ELSEIF ( n>=2 ) THEN
         m = MSTA1(X,200)
         IF ( m<n ) THEN
            n = m
         ELSE
            m = MSTA2(X,n,15)
         ENDIF
         f = 0.0D0
         f2 = 0.0D0
         f1 = 1.0D-100
         DO k = m , 0 , -1
            f = 2.0D0*(v0+k+1.0D0)/X*f1 - f2
            IF ( k<=n ) Bj(k) = f
            f2 = f1
            f1 = f
         ENDDO
         IF ( DABS(bjv0)>DABS(bjv1) ) THEN
            cs = bjv0/f
         ELSE
            cs = bjv1/f2
         ENDIF
         DO k = 0 , n
            Bj(k) = cs*Bj(k)
         ENDDO
      ENDIF
      DO k = 2 , n
         Dj(k) = -(k+v0)/X*Bj(k) + Bj(k-1)
      ENDDO
      IF ( X<=12.0D0 ) THEN
         IF ( v0/=0.0 ) THEN
            bju0 = 0.0D0
            bju1 = 0.0D0
            DO l = 0 , 1
               vl = v0 + l
               bjvl = 1.0D0
               r = 1.0D0
               DO k = 1 , 40
                  r = -0.25D0*r*x2/(k*(k-vl))
                  bjvl = bjvl + r
                  IF ( DABS(r)<DABS(bjvl)*1.0D-15 ) GOTO 30
               ENDDO
 30            vg = 1.0D0 - vl
               CALL GAMMA2(vg,gb)
               b = (2.0D0/X)**vl/gb
               IF ( l==0 ) bju0 = bjvl*b
               IF ( l==1 ) bju1 = bjvl*b
            ENDDO
            pv0 = pi*v0
            pv1 = pi*(1.0D0+v0)
            byv0 = (bjv0*DCOS(pv0)-bju0)/DSIN(pv0)
            byv1 = (bjv1*DCOS(pv1)-bju1)/DSIN(pv1)
         ELSE
            ec = DLOG(X/2.0D0) + el
            cs0 = 0.0D0
            w0 = 0.0D0
            r0 = 1.0D0
            DO k = 1 , 30
               w0 = w0 + 1.0D0/k
               r0 = -0.25D0*r0/(k*k)*x2
               cs0 = cs0 + r0*w0
            ENDDO
            byv0 = rp2*(ec*bjv0-cs0)
            cs1 = 1.0D0
            w1 = 0.0D0
            r1 = 1.0D0
            DO k = 1 , 30
               w1 = w1 + 1.0D0/k
               r1 = -0.25D0*r1/(k*(k+1))*x2
               cs1 = cs1 + r1*(2.0D0*w1+1.0D0/(k+1.0D0))
            ENDDO
            byv1 = rp2*(ec*bjv1-1.0D0/X-0.25D0*X*cs1)
         ENDIF
      ENDIF
      By(0) = byv0
      By(1) = byv1
      DO k = 2 , n
         byvk = 2.0D0*(v0+k-1.0D0)/X*byv1 - byv0
         By(k) = byvk
         byv0 = byv1
         byv1 = byvk
      ENDDO
      Dy(0) = v0/X*By(0) - By(1)
      DO k = 1 , n
         Dy(k) = -(k+v0)/X*By(k) + By(k-1)
      ENDDO
      Vm = n + v0
      END
 
 
 
!       **********************************
 
      SUBROUTINE JYNB(N,X,Nm,Bj,Dj,By,Dy)
!
!       =====================================================
!       Purpose: Compute Bessel functions Jn(x), Yn(x) and
!                their derivatives
!       Input :  x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 )
!                n --- Order of Jn(x) and Yn(x)
!       Output:  BJ(n) --- Jn(x)
!                DJ(n) --- Jn'(x)
!                BY(n) --- Yn(x)
!                DY(n) --- Yn'(x)
!                NM --- Highest order computed
!       Routines called:
!                JYNBH to calculate the Jn and Yn
!       =====================================================
!
      IMPLICIT NONE
!*--JYNB4496
      DOUBLE PRECISION Bj , By , Dj , Dy , X
      INTEGER k , N , Nm
      DIMENSION Bj(0:N) , Dj(0:N) , By(0:N) , Dy(0:N)
      CALL JYNBH(N,0,X,Nm,Bj,By)
!       Compute derivatives by differentiation formulas
      IF ( X<1.0D-100 ) THEN
         DO k = 0 , N
            Dj(k) = 0.0D0
            Dy(k) = 1.0D+300
         ENDDO
         Dj(1) = 0.5D0
      ELSE
         Dj(0) = -Bj(1)
         DO k = 1 , Nm
            Dj(k) = Bj(k-1) - k/X*Bj(k)
         ENDDO
         Dy(0) = -By(1)
         DO k = 1 , Nm
            Dy(k) = By(k-1) - k*By(k)/X
         ENDDO
      ENDIF
      END
 
 
!       **********************************
 
      SUBROUTINE JYNBH(N,Nmin,X,Nm,Bj,By)
!
!       =====================================================
!       Purpose: Compute Bessel functions Jn(x), Yn(x)
!       Input :  x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 )
!                n --- Highest order of Jn(x) and Yn(x) computed  ( n ≥ 0 )
!                nmin -- Lowest order computed  ( nmin ≥ 0 )
!       Output:  BJ(n-NMIN) --- Jn(x)   ; if indexing starts at 0
!                BY(n-NMIN) --- Yn(x)   ; if indexing starts at 0
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 to calculate the starting
!                point for backward recurrence
!       =====================================================
!
      IMPLICIT NONE
!*--JYNBH4542
      DOUBLE PRECISION a , a1 , b , b1 , Bj , bj0 , bj1 , bjk , bs ,    &
                     & By , by0 , by1 , byk , cu , ec , f , f1 , f2 ,   &
                     & p0 , p1
      DOUBLE PRECISION pi , q0 , q1 , r2p , s0 , su , sv , t1 , t2 , X
      INTEGER k , ky , m , MSTA1 , MSTA2 , N , Nm , Nmin
      DIMENSION Bj(0:N-Nmin) , By(0:N-Nmin) , a(4) , b(4) , a1(4) ,     &
              & b1(4)
      pi = 3.141592653589793D0
      r2p = .63661977236758D0
      Nm = N
      IF ( X<1.0D-100 ) THEN
         DO k = Nmin , N
            Bj(k-Nmin) = 0.0D0
            By(k-Nmin) = -1.0D+300
         ENDDO
         IF ( Nmin==0 ) Bj(0) = 1.0D0
         RETURN
      ENDIF
      IF ( X<=300.0 .OR. N>INT(0.9*X) ) THEN
!          Backward recurrence for Jn
         IF ( N==0 ) Nm = 1
         m = MSTA1(X,200)
         IF ( m<Nm ) THEN
            Nm = m
         ELSE
            m = MSTA2(X,Nm,15)
         ENDIF
         bs = 0.0D0
         su = 0.0D0
         sv = 0.0D0
         f2 = 0.0D0
         f1 = 1.0D-100
         f = 0.0D0
         DO k = m , 0 , -1
            f = 2.0D0*(k+1.0D0)/X*f1 - f2
            IF ( k<=Nm .AND. k>=Nmin ) Bj(k-Nmin) = f
            IF ( k==2*INT(k/2) .AND. k/=0 ) THEN
               bs = bs + 2.0D0*f
               su = su + (-1)**(k/2)*f/k
            ELSEIF ( k>1 ) THEN
               sv = sv + (-1)**(k/2)*k/(k*k-1.0D0)*f
            ENDIF
            f2 = f1
            f1 = f
         ENDDO
         s0 = bs + f
         DO k = Nmin , Nm
            Bj(k-Nmin) = Bj(k-Nmin)/s0
         ENDDO
!          Estimates for Yn at start of recurrence
         bj0 = f1/s0
         bj1 = f2/s0
         ec = DLOG(X/2.0D0) + 0.5772156649015329D0
         by0 = r2p*(ec*bj0-4.0D0*su/s0)
         by1 = r2p*((ec-1.0D0)*bj1-bj0/X-4.0D0*sv/s0)
         IF ( 0>=Nmin ) By(0-Nmin) = by0
         IF ( 1>=Nmin ) By(1-Nmin) = by1
         ky = 2
      ELSE
!          Hankel expansion
         DATA a/ - .7031250000000000D-01 , .1121520996093750D+00 ,      &
            & -.5725014209747314D+00 , .6074042001273483D+01/
         DATA b/.7324218750000000D-01 , -.2271080017089844D+00 ,        &
            & .1727727502584457D+01 , -.2438052969955606D+02/
         DATA a1/.1171875000000000D+00 , -.1441955566406250D+00 ,       &
            & .6765925884246826D+00 , -.6883914268109947D+01/
         DATA b1/ - .1025390625000000D+00 , .2775764465332031D+00 ,     &
            & -.1993531733751297D+01 , .2724882731126854D+02/
         t1 = X - 0.25D0*pi
         p0 = 1.0D0
         q0 = -0.125D0/X
         DO k = 1 , 4
            p0 = p0 + a(k)*X**(-2*k)
            q0 = q0 + b(k)*X**(-2*k-1)
         ENDDO
         cu = DSQRT(r2p/X)
         bj0 = cu*(p0*DCOS(t1)-q0*DSIN(t1))
         by0 = cu*(p0*DSIN(t1)+q0*DCOS(t1))
         IF ( 0>=Nmin ) Bj(0-Nmin) = bj0
         IF ( 0>=Nmin ) By(0-Nmin) = by0
         t2 = X - 0.75D0*pi
         p1 = 1.0D0
         q1 = 0.375D0/X
         DO k = 1 , 4
            p1 = p1 + a1(k)*X**(-2*k)
            q1 = q1 + b1(k)*X**(-2*k-1)
         ENDDO
         bj1 = cu*(p1*DCOS(t2)-q1*DSIN(t2))
         by1 = cu*(p1*DSIN(t2)+q1*DCOS(t2))
         IF ( 1>=Nmin ) Bj(1-Nmin) = bj1
         IF ( 1>=Nmin ) By(1-Nmin) = by1
         DO k = 2 , Nm
            bjk = 2.0D0*(k-1.0D0)/X*bj1 - bj0
            IF ( k>=Nmin ) Bj(k-Nmin) = bjk
            bj0 = bj1
            bj1 = bjk
         ENDDO
         ky = 2
      ENDIF
!       Forward recurrence for Yn
      DO k = ky , Nm
         byk = 2.0D0*(k-1.0D0)*by1/X - by0
         IF ( k>=Nmin ) By(k-Nmin) = byk
         by0 = by1
         by1 = byk
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE LEGZO(N,X,W)
!
!       =========================================================
!       Purpose : Compute the zeros of Legendre polynomial Pn(x)
!                 in the interval [-1,1], and the corresponding
!                 weighting coefficients for Gauss-Legendre
!                 integration
!       Input :   n    --- Order of the Legendre polynomial
!       Output:   X(n) --- Zeros of the Legendre polynomial
!                 W(n) --- Corresponding weighting coefficients
!       =========================================================
!
      IMPLICIT NONE
!*--LEGZO4669
      DOUBLE PRECISION f0 , f1 , fd , gd , p , pd , pf , q , W , wp ,   &
                     & X , z , z0
      INTEGER i , j , k , N , n0 , nr
      DIMENSION X(N) , W(N)
      n0 = (N+1)/2
      pf = 0.0D0
      pd = 0.0D0
      DO nr = 1 , n0
         z = DCOS(3.1415926D0*(nr-0.25D0)/N)
 50      z0 = z
         p = 1.0D0
         DO i = 1 , nr - 1
            p = p*(z-X(i))
         ENDDO
         f0 = 1.0D0
         IF ( nr==n0 .AND. N/=2*INT(N/2) ) z = 0.0D0
         f1 = z
         DO k = 2 , N
            pf = (2.0D0-1.0D0/k)*z*f1 - (1.0D0-1.0D0/k)*f0
            pd = k*(f1-z*pf)/(1.0D0-z*z)
            f0 = f1
            f1 = pf
         ENDDO
         IF ( z/=0.0 ) THEN
            fd = pf/p
            q = 0.0D0
            DO i = 1 , nr
               wp = 1.0D0
               DO j = 1 , nr
                  IF ( j/=i ) wp = wp*(z-X(j))
               ENDDO
               q = q + wp
            ENDDO
            gd = (pd-q*fd)/p
            z = z - fd/gd
            IF ( DABS(z-z0)>DABS(z)*1.0D-15 ) GOTO 50
         ENDIF
         X(nr) = z
         X(N+1-nr) = -z
         W(nr) = 2.0D0/((1.0D0-z*z)*pd*pd)
         W(N+1-nr) = W(nr)
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE ASWFA(M,N,C,X,Kd,Cv,S1f,S1d)
!
!       ===========================================================
!       Purpose: Compute the prolate and oblate spheroidal angular
!                functions of the first kind and their derivatives
!       Input :  m  --- Mode parameter,  m = 0,1,2,...
!                n  --- Mode parameter,  n = m,m+1,...
!                c  --- Spheroidal parameter
!                x  --- Argument of angular function, |x| < 1.0
!                KD --- Function code
!                       KD=1 for prolate;  KD=-1 for oblate
!                cv --- Characteristic value
!       Output:  S1F --- Angular function of the first kind
!                S1D --- Derivative of the angular function of
!                        the first kind
!       Routine called:
!                SCKB for computing expansion coefficients ck
!       ===========================================================
!
      IMPLICIT NONE
!*--ASWFA4739
      DOUBLE PRECISION a0 , C , ck , Cv , d0 , d1 , df , eps , r , S1d ,&
                     & S1f , su1 , su2 , X , x0 , x1
      INTEGER ip , k , Kd , M , N , nm , nm2
      DIMENSION ck(200) , df(200)
      eps = 1.0D-14
      x0 = X
      X = DABS(X)
      ip = 1
      IF ( N-M==2*INT((N-M)/2) ) ip = 0
      nm = 40 + INT((N-M)/2+C)
      nm2 = nm/2 - 2
      CALL SDMN(M,N,C,Cv,Kd,df)
      CALL SCKB(M,N,C,df,ck)
      x1 = 1.0D0 - X*X
      IF ( M==0 .AND. x1==0.0D0 ) THEN
         a0 = 1.0D0
      ELSE
         a0 = x1**(0.5D0*M)
      ENDIF
      su1 = ck(1)
      DO k = 1 , nm2
         r = ck(k+1)*x1**k
         su1 = su1 + r
         IF ( k>=10 .AND. DABS(r/su1)<eps ) GOTO 100
      ENDDO
 100  S1f = a0*X**ip*su1
      IF ( X==1.0D0 ) THEN
         IF ( M==0 ) S1d = ip*ck(1) - 2.0D0*ck(2)
         IF ( M==1 ) S1d = -1.0D+100
         IF ( M==2 ) S1d = -2.0D0*ck(1)
         IF ( M>=3 ) S1d = 0.0D0
      ELSE
         d0 = ip - M/x1*X**(ip+1.0D0)
         d1 = -2.0D0*a0*X**(ip+1.0D0)
         su2 = ck(2)
         DO k = 2 , nm2
            r = k*ck(k+1)*x1**(k-1.0D0)
            su2 = su2 + r
            IF ( k>=10 .AND. DABS(r/su2)<eps ) GOTO 150
         ENDDO
 150     S1d = d0*a0*su1 + d1*su2
      ENDIF
      IF ( x0<0.0D0 .AND. ip==0 ) S1d = -S1d
      IF ( x0<0.0D0 .AND. ip==1 ) S1f = -S1f
      X = x0
      END
 
 
 
!       **********************************
 
      SUBROUTINE JYNA(N,X,Nm,Bj,Dj,By,Dy)
!
!       ==========================================================
!       Purpose: Compute Bessel functions Jn(x) & Yn(x) and
!                their derivatives
!       Input :  x --- Argument of Jn(x) & Yn(x)  ( x ≥ 0 )
!                n --- Order of Jn(x) & Yn(x)
!       Output:  BJ(n) --- Jn(x)
!                DJ(n) --- Jn'(x)
!                BY(n) --- Yn(x)
!                DY(n) --- Yn'(x)
!                NM --- Highest order computed
!       Routines called:
!            (1) JY01B to calculate J0(x), J1(x), Y0(x) & Y1(x)
!            (2) MSTA1 and MSTA2 to calculate the starting
!                point for backward recurrence
!       =========================================================
!
      IMPLICIT NONE
!*--JYNA4813
      DOUBLE PRECISION Bj , bj0 , bj1 , bjk , By , by0 , by1 , cs , Dj ,&
                     & dj0 , dj1 , Dy , dy0 , dy1 , f , f0 , f1 , f2 , X
      INTEGER k , m , MSTA1 , MSTA2 , N , Nm
      DIMENSION Bj(0:N) , By(0:N) , Dj(0:N) , Dy(0:N)
      Nm = N
      IF ( X<1.0D-100 ) THEN
         DO k = 0 , N
            Bj(k) = 0.0D0
            Dj(k) = 0.0D0
            By(k) = -1.0D+300
            Dy(k) = 1.0D+300
         ENDDO
         Bj(0) = 1.0D0
         Dj(1) = 0.5D0
         RETURN
      ENDIF
      CALL JY01B(X,bj0,dj0,bj1,dj1,by0,dy0,by1,dy1)
      Bj(0) = bj0
      Bj(1) = bj1
      By(0) = by0
      By(1) = by1
      Dj(0) = dj0
      Dj(1) = dj1
      Dy(0) = dy0
      Dy(1) = dy1
      IF ( N<=1 ) RETURN
      IF ( N<INT(0.9*X) ) THEN
         DO k = 2 , N
            bjk = 2.0D0*(k-1.0D0)/X*bj1 - bj0
            Bj(k) = bjk
            bj0 = bj1
            bj1 = bjk
         ENDDO
      ELSE
         m = MSTA1(X,200)
         IF ( m<N ) THEN
            Nm = m
         ELSE
            m = MSTA2(X,N,15)
         ENDIF
         f2 = 0.0D0
         f1 = 1.0D-100
         f = 0.0D0
         DO k = m , 0 , -1
            f = 2.0D0*(k+1.0D0)/X*f1 - f2
            IF ( k<=Nm ) Bj(k) = f
            f2 = f1
            f1 = f
         ENDDO
         IF ( DABS(bj0)>DABS(bj1) ) THEN
            cs = bj0/f
         ELSE
            cs = bj1/f2
         ENDIF
         DO k = 0 , Nm
            Bj(k) = cs*Bj(k)
         ENDDO
      ENDIF
      DO k = 2 , Nm
         Dj(k) = Bj(k-1) - k/X*Bj(k)
      ENDDO
      f0 = By(0)
      f1 = By(1)
      DO k = 2 , Nm
         f = 2.0D0*(k-1.0D0)/X*f1 - f0
         By(k) = f
         f0 = f1
         f1 = f
      ENDDO
      DO k = 2 , Nm
         Dy(k) = By(k-1) - k*By(k)/X
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE PBDV(V,X,Dv,Dp,Pdf,Pdd)
!
!       ====================================================
!       Purpose: Compute parabolic cylinder functions Dv(x)
!                and their derivatives
!       Input:   x --- Argument of Dv(x)
!                v --- Order of Dv(x)
!       Output:  DV(na) --- Dn+v0(x)
!                DP(na) --- Dn+v0'(x)
!                ( na = |n|, v0 = v-n, |v0| < 1,
!                  n = 0,±1,±2,… )
!                PDF --- Dv(x)
!                PDD --- Dv'(x)
!       Routines called:
!             (1) DVSA for computing Dv(x) for small |x|
!             (2) DVLA for computing Dv(x) for large |x|
!       ====================================================
!
      IMPLICIT NONE
!*--PBDV4914
      DOUBLE PRECISION Dp , Dv , ep , f , f0 , f1 , pd , pd0 , pd1 ,    &
                     & Pdd , Pdf , s0 , V , v0 , v1 , v2 , vh , X , xa
      INTEGER ja , k , l , m , na , nk , nv
      DIMENSION Dv(0:*) , Dp(0:*)
      xa = DABS(X)
      vh = V
      V = V + DSIGN(1.0D0,V)
      nv = INT(V)
      v0 = V - nv
      na = ABS(nv)
      ep = EXP(-.25D0*X*X)
      ja = 0
      IF ( na>=1 ) ja = 1
      IF ( V>=0.0 ) THEN
         IF ( v0==0.0 ) THEN
            pd0 = ep
            pd1 = X*ep
         ELSE
            DO l = 0 , ja
               v1 = v0 + l
               IF ( xa<=5.8 ) CALL DVSA(v1,X,pd1)
               IF ( xa>5.8 ) CALL DVLA(v1,X,pd1)
               IF ( l==0 ) pd0 = pd1
            ENDDO
         ENDIF
         Dv(0) = pd0
         Dv(1) = pd1
         DO k = 2 , na
            Pdf = X*pd1 - (k+v0-1.0D0)*pd0
            Dv(k) = Pdf
            pd0 = pd1
            pd1 = Pdf
         ENDDO
      ELSEIF ( X<=0.0 ) THEN
         IF ( xa<=5.8D0 ) THEN
            CALL DVSA(v0,X,pd0)
            v1 = v0 - 1.0D0
            CALL DVSA(v1,X,pd1)
         ELSE
            CALL DVLA(v0,X,pd0)
            v1 = v0 - 1.0D0
            CALL DVLA(v1,X,pd1)
         ENDIF
         Dv(0) = pd0
         Dv(1) = pd1
         DO k = 2 , na
            pd = (-X*pd1+pd0)/(k-1.0D0-v0)
            Dv(k) = pd
            pd0 = pd1
            pd1 = pd
         ENDDO
      ELSEIF ( X<=2.0 ) THEN
         v2 = nv + v0
         IF ( nv==0 ) v2 = v2 - 1.0D0
         nk = INT(-v2)
         CALL DVSA(v2,X,f1)
         v1 = v2 + 1.0D0
         CALL DVSA(v1,X,f0)
         Dv(nk) = f1
         Dv(nk-1) = f0
         DO k = nk - 2 , 0 , -1
            f = X*f0 + (k-v0+1.0D0)*f1
            Dv(k) = f
            f1 = f0
            f0 = f
         ENDDO
      ELSE
         IF ( xa<=5.8 ) CALL DVSA(v0,X,pd0)
         IF ( xa>5.8 ) CALL DVLA(v0,X,pd0)
         Dv(0) = pd0
         m = 100 + na
         f1 = 0.0D0
         f0 = 1.0D-30
         f = 0.0D0
         DO k = m , 0 , -1
            f = X*f0 + (k-v0+1.0D0)*f1
            IF ( k<=na ) Dv(k) = f
            f1 = f0
            f0 = f
         ENDDO
         s0 = pd0/f
         DO k = 0 , na
            Dv(k) = s0*Dv(k)
         ENDDO
      ENDIF
      DO k = 0 , na - 1
         v1 = ABS(v0) + k
         IF ( V>=0.0D0 ) THEN
            Dp(k) = 0.5D0*X*Dv(k) - Dv(k+1)
         ELSE
            Dp(k) = -0.5D0*X*Dv(k) - v1*Dv(k+1)
         ENDIF
      ENDDO
      Pdf = Dv(na-1)
      Pdd = Dp(na-1)
      V = vh
      END
 
 
 
!       **********************************
 
      SUBROUTINE ITSH0(X,Th0)
!
!       ===================================================
!       Purpose: Evaluate the integral of Struve function
!                H0(t) with respect to t from 0 and x
!       Input :  x   --- Upper limit  ( x ≥ 0 )
!       Output:  TH0 --- Integration of H0(t) from 0 and x
!       ===================================================
!
      IMPLICIT NONE
!*--ITSH05030
      DOUBLE PRECISION a , a0 , a1 , af , bf , bg , el , pi , r , rd ,  &
                     & s , s0 , Th0 , ty , X , xp
      INTEGER k
      DIMENSION a(25)
      pi = 3.141592653589793D0
      r = 1.0D0
      IF ( X<=30.0 ) THEN
         s = 0.5D0
         DO k = 1 , 100
            rd = 1.0D0
            IF ( k==1 ) rd = 0.5D0
            r = -r*rd*k/(k+1.0D0)*(X/(2.0D0*k+1.0D0))**2
            s = s + r
            IF ( DABS(r)<DABS(s)*1.0D-12 ) GOTO 50
         ENDDO
 50      Th0 = 2.0D0/pi*X*X*s
      ELSE
         s = 1.0D0
         DO k = 1 , 12
            r = -r*k/(k+1.0D0)*((2.0D0*k+1.0D0)/X)**2
            s = s + r
            IF ( DABS(r)<DABS(s)*1.0D-12 ) GOTO 100
         ENDDO
 100     el = .57721566490153D0
         s0 = s/(pi*X*X) + 2.0D0/pi*(DLOG(2.0D0*X)+el)
         a0 = 1.0D0
         a1 = 5.0D0/8.0D0
         a(1) = a1
         DO k = 1 , 20
            af = ((1.5D0*(k+.5D0)*(k+5.0D0/6.0D0)*a1-.5D0*(k+.5D0)      &
               & *(k+.5D0)*(k-.5D0)*a0))/(k+1.0D0)
            a(k+1) = af
            a0 = a1
            a1 = af
         ENDDO
         bf = 1.0D0
         r = 1.0D0
         DO k = 1 , 10
            r = -r/(X*X)
            bf = bf + a(2*k)*r
         ENDDO
         bg = a(1)/X
         r = 1.0D0/X
         DO k = 1 , 10
            r = -r/(X*X)
            bg = bg + a(2*k+1)*r
         ENDDO
         xp = X + .25D0*pi
         ty = DSQRT(2.0D0/(pi*X))*(bg*DCOS(xp)-bf*DSIN(xp))
         Th0 = ty + s0
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE CERZO(Nt,Zo)
!
!       ===============================================================
!       Purpose : Evaluate the complex zeros of error function erf(z)
!                 using the modified Newton's iteration method
!       Input :   NT --- Total number of zeros
!       Output:   ZO(L) --- L-th zero of erf(z), L=1,2,...,NT
!       Routine called: CERF for computing erf(z) and erf'(z)
!       ===============================================================
!
      IMPLICIT NONE
!*--CERZO5100
      INTEGER i , it , j , nr , Nt
      DOUBLE PRECISION pi , pu , pv , px , py , w , w0
      COMPLEX*16 z , zd , zf , zfd , zgd , Zo , zp , zq , zw
      DIMENSION Zo(Nt)
      pi = 3.141592653589793D0
      w = 0.0D0
      DO nr = 1 , Nt
         pu = DSQRT(pi*(4.0D0*nr-0.5D0))
         pv = pi*DSQRT(2.0D0*nr-0.25D0)
         px = 0.5*pu - 0.5*DLOG(pv)/pu
         py = 0.5*pu + 0.5*DLOG(pv)/pu
         z = DCMPLX(px,py)
         it = 0
 50      it = it + 1
         CALL CERF(z,zf,zd)
         zp = (1.0D0,0.0D0)
         DO i = 1 , nr - 1
            zp = zp*(z-Zo(i))
         ENDDO
         zfd = zf/zp
         zq = (0.0D0,0.0D0)
         DO i = 1 , nr - 1
            zw = (1.0D0,0.0D0)
            DO j = 1 , nr - 1
               IF ( j/=i ) zw = zw*(z-Zo(j))
            ENDDO
            zq = zq + zw
         ENDDO
         zgd = (zd-zq*zfd)/zp
         z = z - zfd/zgd
         w0 = w
         w = ABS(z)
         IF ( it<=50 .AND. DABS((w-w0)/w)>1.0D-11 ) GOTO 50
         Zo(nr) = z
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE GAMMA2(X,Ga)
!
!       ==================================================
!       Purpose: Compute gamma function Г(x)
!       Input :  x  --- Argument of Г(x)
!                       ( x is not equal to 0,-1,-2,…)
!       Output:  GA --- Г(x)
!       ==================================================
!
      IMPLICIT NONE
!*--GAMMA25155
      DOUBLE PRECISION g , Ga , gr , pi , r , X , z
      INTEGER k , m , m1
      DIMENSION g(26)
      pi = 3.141592653589793D0
      IF ( X/=INT(X) ) THEN
         r = 1.0D0
         IF ( DABS(X)>1.0D0 ) THEN
            z = DABS(X)
            m = INT(z)
            DO k = 1 , m
               r = r*(z-k)
            ENDDO
            z = z - m
         ELSE
            z = X
         ENDIF
         DATA g/1.0D0 , 0.5772156649015329D0 , -0.6558780715202538D0 ,  &
            & -0.420026350340952D-1 , 0.1665386113822915D0 ,            &
            & -.421977345555443D-1 , -.96219715278770D-2 ,              &
            & .72189432466630D-2 , -.11651675918591D-2 ,                &
            & -.2152416741149D-3 , .1280502823882D-3 ,                  &
            & -.201348547807D-4 , -.12504934821D-5 , .11330272320D-5 ,  &
            & -.2056338417D-6 , .61160950D-8 , .50020075D-8 ,           &
            & -.11812746D-8 , .1043427D-9 , .77823D-11 , -.36968D-11 ,  &
            & .51D-12 , -.206D-13 , -.54D-14 , .14D-14 , .1D-15/
         gr = g(26)
         DO k = 25 , 1 , -1
            gr = gr*z + g(k)
         ENDDO
         Ga = 1.0D0/(gr*z)
         IF ( DABS(X)>1.0D0 ) THEN
            Ga = Ga*r
            IF ( X<0.0D0 ) Ga = -pi/(X*Ga*DSIN(pi*X))
         ENDIF
      ELSEIF ( X>0.0D0 ) THEN
         Ga = 1.0D0
         m1 = X - 1
         DO k = 2 , m1
            Ga = Ga*k
         ENDDO
      ELSE
         Ga = 1.0D+300
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE CHGU(A,B,X,Hu,Md,Isfer)
!
!       =======================================================
!       Purpose: Compute the confluent hypergeometric function
!                U(a,b,x)
!       Input  : a  --- Parameter
!                b  --- Parameter
!                x  --- Argument  ( x > 0 )
!       Output:  HU --- U(a,b,x)
!                MD --- Method code
!                ISFER --- Error flag
!       Routines called:
!            (1) CHGUS for small x ( MD=1 )
!            (2) CHGUL for large x ( MD=2 )
!            (3) CHGUBI for integer b ( MD=3 )
!            (4) CHGUIT for numerical integration ( MD=4 )
!       =======================================================
!
      IMPLICIT NONE
!*--CHGU5225
      DOUBLE PRECISION A , a00 , aa , B , b00 , Hu , hu1 , X
      INTEGER id , id1 , Isfer , Md
      LOGICAL il1 , il2 , il3 , bl1 , bl2 , bl3 , bn
      aa = A - B + 1.0D0
      Isfer = 0
      il1 = A==INT(A) .AND. A<=0.0
      il2 = aa==INT(aa) .AND. aa<=0.0
      il3 = ABS(A*(A-B+1.0))/X<=2.0
      bl1 = X<=5.0 .OR. (X<=10.0 .AND. A<=2.0)
      bl2 = (X>5.0 .AND. X<=12.5) .AND. (A>=1.0 .AND. B>=A+4.0)
      bl3 = X>12.5 .AND. A>=5.0 .AND. B>=A + 5.0
      bn = B==INT(B) .AND. B/=0.0
      id1 = -100
      hu1 = 0.0D0
      IF ( B/=INT(B) ) THEN
         CALL CHGUS(A,B,X,Hu,id1)
         Md = 1
         IF ( id1>=9 ) RETURN
         hu1 = Hu
      ENDIF
      IF ( il1 .OR. il2 .OR. il3 ) THEN
         CALL CHGUL(A,B,X,Hu,id)
         Md = 2
         IF ( id>=9 ) RETURN
         IF ( id1>id ) THEN
            Md = 1
            id = id1
            Hu = hu1
         ENDIF
      ENDIF
      IF ( A>=1.0 ) THEN
         IF ( bn .AND. (bl1 .OR. bl2 .OR. bl3) ) THEN
            CALL CHGUBI(A,B,X,Hu,id)
            Md = 3
         ELSE
            CALL CHGUIT(A,B,X,Hu,id)
            Md = 4
         ENDIF
      ELSEIF ( B<=A ) THEN
         a00 = A
         b00 = B
         A = A - B + 1.0D0
         B = 2.0D0 - B
         CALL CHGUIT(A,B,X,Hu,id)
         Hu = X**(1.0D0-b00)*Hu
         A = a00
         B = b00
         Md = 4
      ELSEIF ( bn .AND. (.NOT.il1) ) THEN
         CALL CHGUBI(A,B,X,Hu,id)
         Md = 3
      ENDIF
      IF ( id<6 ) Isfer = 6
      END
 
 
 
!       **********************************
 
      SUBROUTINE LAMN(N,X,Nm,Bl,Dl)
!
!       =========================================================
!       Purpose: Compute lambda functions and their derivatives
!       Input:   x --- Argument of lambda function
!                n --- Order of lambda function
!       Output:  BL(n) --- Lambda function of order n
!                DL(n) --- Derivative of lambda function
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 for computing the start
!                point for backward recurrence
!       =========================================================
!
      IMPLICIT NONE
!*--LAMN5303
      DOUBLE PRECISION bg , bk , Bl , bs , Dl , f , f0 , f1 , r , r0 ,  &
                     & uk , X , x2
      INTEGER i , k , m , MSTA1 , MSTA2 , N , Nm
      DIMENSION Bl(0:N) , Dl(0:N)
      Nm = N
      IF ( DABS(X)<1.0D-100 ) THEN
         DO k = 0 , N
            Bl(k) = 0.0D0
            Dl(k) = 0.0D0
         ENDDO
         Bl(0) = 1.0D0
         Dl(1) = 0.5D0
         RETURN
      ENDIF
      IF ( X<=12.0D0 ) THEN
         x2 = X*X
         DO k = 0 , N
            bk = 1.0D0
            r = 1.0D0
            DO i = 1 , 50
               r = -0.25D0*r*x2/(i*(i+k))
               bk = bk + r
               IF ( DABS(r)<DABS(bk)*1.0D-15 ) GOTO 20
            ENDDO
 20         Bl(k) = bk
            IF ( k>=1 ) Dl(k-1) = -0.5D0*X/k*bk
         ENDDO
         uk = 1.0D0
         r = 1.0D0
         DO i = 1 , 50
            r = -0.25D0*r*x2/(i*(i+N+1.0D0))
            uk = uk + r
            IF ( DABS(r)<DABS(uk)*1.0D-15 ) GOTO 50
         ENDDO
 50      Dl(N) = -0.5D0*X/(N+1.0D0)*uk
         RETURN
      ENDIF
      IF ( N==0 ) Nm = 1
      m = MSTA1(X,200)
      IF ( m<Nm ) THEN
         Nm = m
      ELSE
         m = MSTA2(X,Nm,15)
      ENDIF
      bs = 0.0D0
      f = 0.0D0
      f0 = 0.0D0
      f1 = 1.0D-100
      DO k = m , 0 , -1
         f = 2.0D0*(k+1.0D0)*f1/X - f0
         IF ( k<=Nm ) Bl(k) = f
         IF ( k==2*INT(k/2) ) bs = bs + 2.0D0*f
         f0 = f1
         f1 = f
      ENDDO
      bg = bs - f
      DO k = 0 , Nm
         Bl(k) = Bl(k)/bg
      ENDDO
      r0 = 1.0D0
      DO k = 1 , Nm
         r0 = 2.0D0*r0*k/X
         Bl(k) = r0*Bl(k)
      ENDDO
      Dl(0) = -0.5D0*X*Bl(1)
      DO k = 1 , Nm
         Dl(k) = 2.0D0*k/X*(Bl(k-1)-Bl(k))
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE COMELP(Hk,Ck,Ce)
!
!       ==================================================
!       Purpose: Compute complete elliptic integrals K(k)
!                and E(k)
!       Input  : K  --- Modulus k ( 0 ≤ k ≤ 1 )
!       Output : CK --- K(k)
!                CE --- E(k)
!       ==================================================
!
      IMPLICIT NONE
!*--COMELP5392
      DOUBLE PRECISION ae , ak , be , bk , Ce , Ck , Hk , pk
      pk = 1.0D0 - Hk*Hk
      IF ( Hk==1.0 ) THEN
         Ck = 1.0D+300
         Ce = 1.0D0
      ELSE
         ak = (((.01451196212D0*pk+.03742563713D0)*pk+.03590092383D0)   &
            & *pk+.09666344259D0)*pk + 1.38629436112D0
         bk = (((.00441787012D0*pk+.03328355346D0)*pk+.06880248576D0)   &
            & *pk+.12498593597D0)*pk + .5D0
         Ck = ak - bk*DLOG(pk)
         ae = (((.01736506451D0*pk+.04757383546D0)*pk+.0626060122D0)    &
            & *pk+.44325141463D0)*pk + 1.0D0
         be = (((.00526449639D0*pk+.04069697526D0)*pk+.09200180037D0)   &
            & *pk+.2499836831D0)*pk
         Ce = ae - be*DLOG(pk)
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE INCOB(A,B,X,Bix)
!
!       ========================================================
!       Purpose: Compute the incomplete beta function Ix(a,b)
!       Input :  a --- Parameter
!                b --- Parameter
!                x --- Argument ( 0 ≤ x ≤ 1 )
!       Output:  BIX --- Ix(a,b)
!       Routine called: BETA for computing beta function B(p,q)
!       ========================================================
!
      IMPLICIT NONE
!*--INCOB5429
      DOUBLE PRECISION A , B , Bix , bt , dk , fk , s0 , t1 , t2 , ta , &
                     & tb , X
      INTEGER k
      DIMENSION dk(51) , fk(51)
      s0 = (A+1.0D0)/(A+B+2.0D0)
      CALL BETA(A,B,bt)
      IF ( X<=s0 ) THEN
         DO k = 1 , 20
            dk(2*k) = k*(B-k)*X/(A+2.0D0*k-1.0D0)/(A+2.0D0*k)
         ENDDO
         DO k = 0 , 20
            dk(2*k+1) = -(A+k)*(A+B+k)*X/(A+2.D0*k)/(A+2.0*k+1.0)
         ENDDO
         t1 = 0.0D0
         DO k = 20 , 1 , -1
            t1 = dk(k)/(1.0D0+t1)
         ENDDO
         ta = 1.0D0/(1.0D0+t1)
         Bix = X**A*(1.0D0-X)**B/(A*bt)*ta
      ELSE
         DO k = 1 , 20
            fk(2*k) = k*(A-k)*(1.0D0-X)/(B+2.*k-1.0)/(B+2.0*k)
         ENDDO
         DO k = 0 , 20
            fk(2*k+1) = -(B+k)*(A+B+k)*(1.D0-X)/(B+2.D0*k)              &
                      & /(B+2.D0*k+1.D0)
         ENDDO
         t2 = 0.0D0
         DO k = 20 , 1 , -1
            t2 = fk(k)/(1.0D0+t2)
         ENDDO
         tb = 1.0D0/(1.0D0+t2)
         Bix = 1.0D0 - X**A*(1.0D0-X)**B/(B*bt)*tb
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE CVF(Kd,M,Q,A,Mj,F)
!
!       ======================================================
!       Purpose: Compute the value of F for characteristic
!                equation of Mathieu functions
!       Input :  m --- Order of Mathieu functions
!                q --- Parameter of Mathieu functions
!                A --- Characteristic value
!       Output:  F --- Value of F for characteristic equation
!       ======================================================
!
      IMPLICIT NONE
!*--CVF5485
      DOUBLE PRECISION A , b , F , Q , t0 , t1 , t2
      INTEGER ic , j , j0 , jf , Kd , l , l0 , M , Mj
      b = A
      ic = INT(M/2)
      l = 0
      l0 = 0
      j0 = 2
      jf = ic
      IF ( Kd==1 ) l0 = 2
      IF ( Kd==1 ) j0 = 3
      IF ( Kd==2 .OR. Kd==3 ) l = 1
      IF ( Kd==4 ) jf = ic - 1
      t1 = 0.0D0
      DO j = Mj , ic + 1 , -1
         t1 = -Q*Q/((2.0D0*j+l)**2-b+t1)
      ENDDO
      IF ( M<=2 ) THEN
         t2 = 0.0D0
         IF ( Kd==1 .AND. M==0 ) t1 = t1 + t1
         IF ( Kd==1 .AND. M==2 ) t1 = -2.0D0*Q*Q/(4.0D0-b+t1) - 4.0D0
         IF ( Kd==2 .AND. M==1 ) t1 = t1 + Q
         IF ( Kd==3 .AND. M==1 ) t1 = t1 - Q
      ELSE
         t0 = 0.0D0
         IF ( Kd==1 ) t0 = 4.0D0 - b + 2.0D0*Q*Q/b
         IF ( Kd==2 ) t0 = 1.0D0 - b + Q
         IF ( Kd==3 ) t0 = 1.0D0 - b - Q
         IF ( Kd==4 ) t0 = 4.0D0 - b
         t2 = -Q*Q/t0
         DO j = j0 , jf
            t2 = -Q*Q/((2.0D0*j-l-l0)**2-b+t2)
         ENDDO
      ENDIF
      F = (2.0D0*ic+l)**2 + t1 + t2 - b
      END
 
 
 
!       **********************************
 
      SUBROUTINE CLPN(N,X,Y,Cpn,Cpd)
!
!       ==================================================
!       Purpose: Compute Legendre polynomials Pn(z) and
!                their derivatives Pn'(z) for a complex
!                argument
!       Input :  x --- Real part of z
!                y --- Imaginary part of z
!                n --- Degree of Pn(z), n = 0,1,2,...
!       Output:  CPN(n) --- Pn(z)
!                CPD(n) --- Pn'(z)
!       ==================================================
!
      IMPLICIT NONE
!*--CLPN5543
      COMPLEX*16 cp0 , cp1 , Cpd , cpf , Cpn , z
      INTEGER k , N
      DOUBLE PRECISION X , Y
      DIMENSION Cpn(0:N) , Cpd(0:N)
      z = DCMPLX(X,Y)
      Cpn(0) = (1.0D0,0.0D0)
      Cpn(1) = z
      Cpd(0) = (0.0D0,0.0D0)
      Cpd(1) = (1.0D0,0.0D0)
      cp0 = (1.0D0,0.0D0)
      cp1 = z
      DO k = 2 , N
         cpf = (2.0D0*k-1.0D0)/k*z*cp1 - (k-1.0D0)/k*cp0
         Cpn(k) = cpf
         IF ( DABS(X)==1.0D0 .AND. Y==0.0D0 ) THEN
            Cpd(k) = 0.5D0*X**(k+1)*k*(k+1.0D0)
         ELSE
            Cpd(k) = k*(cp1-z*cpf)/(1.0D0-z*z)
         ENDIF
         cp0 = cp1
         cp1 = cpf
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE LQMNS(M,N,X,Qm,Qd)
!
!       ========================================================
!       Purpose: Compute associated Legendre functions Qmn(x)
!                and Qmn'(x) for a given order
!       Input :  x --- Argument of Qmn(x)
!                m --- Order of Qmn(x),  m = 0,1,2,...
!                n --- Degree of Qmn(x), n = 0,1,2,...
!       Output:  QM(n) --- Qmn(x)
!                QD(n) --- Qmn'(x)
!       ========================================================
!
      IMPLICIT NONE
!*--LQMNS5586
      INTEGER k , km , l , ls , M , N
      DOUBLE PRECISION q0 , q00 , q01 , q0l , q10 , q11 , q1l , Qd ,    &
                     & qf0 , qf1 , qf2 , qg0 , qg1 , qh0 , qh1 , qh2 ,  &
                     & Qm , qm0 , qm1 , qmk
      DOUBLE PRECISION X , xq
      DIMENSION Qm(0:N) , Qd(0:N)
      DO k = 0 , N
         Qm(k) = 0.0D0
         Qd(k) = 0.0D0
      ENDDO
      IF ( DABS(X)==1.0D0 ) THEN
         DO k = 0 , N
            Qm(k) = 1.0D+300
            Qd(k) = 1.0D+300
         ENDDO
         RETURN
      ENDIF
      ls = 1
      IF ( DABS(X)>1.0D0 ) ls = -1
      xq = DSQRT(ls*(1.0D0-X*X))
      q0 = 0.5D0*DLOG(DABS((X+1.0)/(X-1.0)))
      q00 = q0
      q10 = -1.0D0/xq
      q01 = X*q0 - 1.0D0
      q11 = -ls*xq*(q0+X/(1.0D0-X*X))
      qf0 = q00
      qf1 = q10
      qm0 = 0.0D0
      qm1 = 0.0D0
      DO k = 2 , M
         qm0 = -2.0D0*(k-1.0)/xq*X*qf1 - ls*(k-1.0)*(2.0-k)*qf0
         qf0 = qf1
         qf1 = qm0
      ENDDO
      IF ( M==0 ) qm0 = q00
      IF ( M==1 ) qm0 = q10
      Qm(0) = qm0
      IF ( DABS(X)<1.0001D0 ) THEN
         IF ( M==0 .AND. N>0 ) THEN
            qf0 = q00
            qf1 = q01
            DO k = 2 , N
               qf2 = ((2.0*k-1.0D0)*X*qf1-(k-1.0)*qf0)/k
               Qm(k) = qf2
               qf0 = qf1
               qf1 = qf2
            ENDDO
         ENDIF
         qg0 = q01
         qg1 = q11
         DO k = 2 , M
            qm1 = -2.0D0*(k-1.0)/xq*X*qg1 - ls*k*(3.0-k)*qg0
            qg0 = qg1
            qg1 = qm1
         ENDDO
         IF ( M==0 ) qm1 = q01
         IF ( M==1 ) qm1 = q11
         Qm(1) = qm1
         IF ( M==1 .AND. N>1 ) THEN
            qh0 = q10
            qh1 = q11
            DO k = 2 , N
               qh2 = ((2.0*k-1.0D0)*X*qh1-k*qh0)/(k-1.0)
               Qm(k) = qh2
               qh0 = qh1
               qh1 = qh2
            ENDDO
         ELSEIF ( M>=2 ) THEN
            qg0 = q00
            qg1 = q01
            qh0 = q10
            qh1 = q11
            qmk = 0.0D0
            DO l = 2 , N
               q0l = ((2.0D0*l-1.0D0)*X*qg1-(l-1.0D0)*qg0)/l
               q1l = ((2.0*l-1.0D0)*X*qh1-l*qh0)/(l-1.0D0)
               qf0 = q0l
               qf1 = q1l
               DO k = 2 , M
                  qmk = -2.0D0*(k-1.0)/xq*X*qf1 - ls*(k+l-1.0)*(l+2.0-k)&
                      & *qf0
                  qf0 = qf1
                  qf1 = qmk
               ENDDO
               Qm(l) = qmk
               qg0 = qg1
               qg1 = q0l
               qh0 = qh1
               qh1 = q1l
            ENDDO
         ENDIF
      ELSE
         IF ( DABS(X)>1.1 ) THEN
            km = 40 + M + N
         ELSE
            km = (40+M+N)*INT(-1.0-1.8*LOG(X-1.0))
         ENDIF
         qf2 = 0.0D0
         qf1 = 1.0D0
         DO k = km , 0 , -1
            qf0 = ((2.0*k+3.0D0)*X*qf1-(k+2.0-M)*qf2)/(k+M+1.0)
            IF ( k<=N ) Qm(k) = qf0
            qf2 = qf1
            qf1 = qf0
         ENDDO
         DO k = 0 , N
            Qm(k) = Qm(k)*qm0/qf0
         ENDDO
      ENDIF
      IF ( DABS(X)<1.0D0 ) THEN
         DO k = 0 , N
            Qm(k) = (-1)**M*Qm(k)
         ENDDO
      ENDIF
      Qd(0) = ((1.0D0-M)*Qm(1)-X*Qm(0))/(X*X-1.0)
      DO k = 1 , N
         Qd(k) = (k*X*Qm(k)-(k+M)*Qm(k-1))/(X*X-1.0)
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE CIKLV(V,Z,Cbiv,Cdiv,Cbkv,Cdkv)
!
!       =====================================================
!       Purpose: Compute modified Bessel functions Iv(z) and
!                Kv(z) and their derivatives with a complex
!                argument and a large order
!       Input:   v --- Order of Iv(z) and Kv(z)
!                z --- Complex argument
!       Output:  CBIV --- Iv(z)
!                CDIV --- Iv'(z)
!                CBKV --- Kv(z)
!                CDKV --- Kv'(z)
!       Routine called:
!                CJK to compute the expansion coefficients
!       ====================================================
!
      IMPLICIT NONE
!*--CIKLV5729
      DOUBLE PRECISION a , pi , V , v0 , vr
      COMPLEX*16 Cbiv , Cbkv , Cdiv , Cdkv , ceta , cf , cfi , cfk ,    &
               & csi , csk , ct , ct2 , cws , Z
      INTEGER i , k , km , l , l0 , lf
      DIMENSION cf(12) , a(91)
      pi = 3.141592653589793D0
      km = 12
      CALL CJK(km,a)
      DO l = 1 , 0 , -1
         v0 = V - l
         cws = SQRT(1.0D0+(Z/v0)*(Z/v0))
         ceta = cws + LOG(Z/v0/(1.0D0+cws))
         ct = 1.0D0/cws
         ct2 = ct*ct
         DO k = 1 , km
            l0 = k*(k+1)/2 + 1
            lf = l0 + k
            cf(k) = a(lf)
            DO i = lf - 1 , l0 , -1
               cf(k) = cf(k)*ct2 + a(i)
            ENDDO
            cf(k) = cf(k)*ct**k
         ENDDO
         vr = 1.0D0/v0
         csi = (1.0D0,0.0D0)
         DO k = 1 , km
            csi = csi + cf(k)*vr**k
         ENDDO
         Cbiv = SQRT(ct/(2.0D0*pi*v0))*EXP(v0*ceta)*csi
         IF ( l==1 ) cfi = Cbiv
         csk = (1.0D0,0.0D0)
         DO k = 1 , km
            csk = csk + (-1)**k*cf(k)*vr**k
         ENDDO
         Cbkv = SQRT(pi*ct/(2.0D0*v0))*EXP(-v0*ceta)*csk
         IF ( l==1 ) cfk = Cbkv
      ENDDO
      Cdiv = cfi - V/Z*Cbiv
      Cdkv = -cfk - V/Z*Cbkv
      END
 
 
 
!       **********************************
 
      SUBROUTINE ELIT(Hk,Phi,Fe,Ee)
!
!       ==================================================
!       Purpose: Compute complete and incomplete elliptic
!                integrals F(k,phi) and E(k,phi)
!       Input  : HK  --- Modulus k ( 0 ≤ k ≤ 1 )
!                Phi --- Argument ( in degrees )
!       Output : FE  --- F(k,phi)
!                EE  --- E(k,phi)
!       ==================================================
!
      IMPLICIT NONE
!*--ELIT5790
      DOUBLE PRECISION a , a0 , b , b0 , c , ce , ck , d , d0 , Ee ,    &
                     & fac , Fe , g , Hk , Phi , pi , r
      INTEGER n
      g = 0.0D0
      pi = 3.14159265358979D0
      a0 = 1.0D0
      b0 = DSQRT(1.0D0-Hk*Hk)
      d0 = (pi/180.0D0)*Phi
      r = Hk*Hk
      IF ( Hk==1.0D0 .AND. Phi==90.0D0 ) THEN
         Fe = 1.0D+300
         Ee = 1.0D0
      ELSEIF ( Hk==1.0D0 ) THEN
         Fe = DLOG((1.0D0+DSIN(d0))/DCOS(d0))
         Ee = DSIN(d0)
      ELSE
         fac = 1.0D0
         d = 0.0D0
         DO n = 1 , 40
            a = (a0+b0)/2.0D0
            b = DSQRT(a0*b0)
            c = (a0-b0)/2.0D0
            fac = 2.0D0*fac
            r = r + fac*c*c
            IF ( Phi/=90.0D0 ) THEN
               d = d0 + DATAN((b0/a0)*DTAN(d0))
               g = g + c*DSIN(d)
               d0 = d + pi*INT(d/pi+.5D0)
            ENDIF
            a0 = a
            b0 = b
            IF ( c<1.0D-7 ) GOTO 50
         ENDDO
 50      ck = pi/(2.0D0*a)
         ce = pi*(2.0D0-r)/(4.0D0*a)
         IF ( Phi==90.0D0 ) THEN
            Fe = ck
            Ee = ce
         ELSE
            Fe = d/(fac*a)
            Ee = Fe*ce/ck + g
         ENDIF
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE ELIT3(Phi,Hk,C,El3)
!
!       =========================================================
!       Purpose: Compute the elliptic integral of the third kind
!                using Gauss-Legendre quadrature
!       Input :  Phi --- Argument ( in degrees )
!                 k  --- Modulus   ( 0 ≤ k ≤ 1.0 )
!                 c  --- Parameter ( 0 ≤ c ≤ 1.0 )
!       Output:  EL3 --- Value of the elliptic integral of the
!                        third kind
!       =========================================================
!
      IMPLICIT NONE
!*--ELIT35854
      DOUBLE PRECISION C , c0 , c1 , c2 , El3 , f1 , f2 , Hk , Phi , t ,&
                     & t1 , t2 , w
      INTEGER i
      DIMENSION t(10) , w(10)
      LOGICAL lb1 , lb2
      DATA t/.9931285991850949D0 , .9639719272779138D0 ,                &
         & .9122344282513259D0 , .8391169718222188D0 ,                  &
         & .7463319064601508D0 , .6360536807265150D0 ,                  &
         & .5108670019508271D0 , .3737060887154195D0 ,                  &
         & .2277858511416451D0 , .7652652113349734D-1/
      DATA w/.1761400713915212D-1 , .4060142980038694D-1 ,              &
         & .6267204833410907D-1 , .8327674157670475D-1 ,                &
         & .1019301198172404D0 , .1181945319615184D0 ,                  &
         & .1316886384491766D0 , .1420961093183820D0 ,                  &
         & .1491729864726037D0 , .1527533871307258D0/
      lb1 = Hk==1.0D0 .AND. DABS(Phi-90.0)<=1.0D-8
      lb2 = C==1.0D0 .AND. DABS(Phi-90.0)<=1.0D-8
      IF ( lb1 .OR. lb2 ) THEN
         El3 = 1.0D+300
         RETURN
      ENDIF
      c1 = 0.87266462599716D-2*Phi
      c2 = c1
      El3 = 0.0D0
      DO i = 1 , 10
         c0 = c2*t(i)
         t1 = c1 + c0
         t2 = c1 - c0
         f1 = 1.0D0/((1.0D0-C*DSIN(t1)*DSIN(t1))                        &
            & *DSQRT(1.0D0-Hk*Hk*DSIN(t1)*DSIN(t1)))
         f2 = 1.0D0/((1.0D0-C*DSIN(t2)*DSIN(t2))                        &
            & *DSQRT(1.0D0-Hk*Hk*DSIN(t2)*DSIN(t2)))
         El3 = El3 + w(i)*(f1+f2)
      ENDDO
      El3 = c1*El3
      END
 
!       **********************************
 
      SUBROUTINE EIX(X,Ei)
!
!       ============================================
!       Purpose: Compute exponential integral Ei(x)
!       Input :  x  --- Argument of Ei(x)
!       Output:  EI --- Ei(x)
!       ============================================
!
      IMPLICIT NONE
!*--EIX5906
      DOUBLE PRECISION Ei , ga , r , X
      INTEGER k
      IF ( X==0.0 ) THEN
         Ei = -1.0D+300
      ELSEIF ( X<0 ) THEN
         CALL E1XB(-X,Ei)
         Ei = -Ei
      ELSEIF ( DABS(X)<=40.0 ) THEN
!          Power series around x=0
         Ei = 1.0D0
         r = 1.0D0
         DO k = 1 , 100
            r = r*k*X/(k+1.0D0)**2
            Ei = Ei + r
            IF ( DABS(r/Ei)<=1.0D-15 ) GOTO 50
         ENDDO
 50      ga = 0.5772156649015328D0
         Ei = ga + DLOG(X) + X*Ei
      ELSE
!          Asymptotic expansion (the series is not convergent)
         Ei = 1.0D0
         r = 1.0D0
         DO k = 1 , 20
            r = r*k/X
            Ei = Ei + r
         ENDDO
         Ei = EXP(X)/X*Ei
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE EIXZ(Z,Cei)
!
!       ============================================
!       Purpose: Compute exponential integral Ei(x)
!       Input :  x  --- Complex argument of Ei(x)
!       Output:  EI --- Ei(x)
!       ============================================
!
      IMPLICIT NONE
!*--EIXZ5951
      DOUBLE COMPLEX Z , Cei
      DOUBLE PRECISION pi
      pi = 3.141592653589793D0
      CALL E1Z(-Z,Cei)
      Cei = -Cei
      IF ( DIMAG(Z)>0 ) THEN
         Cei = Cei + (0D0,1D0)*pi
      ELSEIF ( DIMAG(Z)<0 ) THEN
         Cei = Cei - (0D0,1D0)*pi
      ELSEIF ( DIMAG(Z)==0 ) THEN
         IF ( DBLE(Z)>0 ) Cei = Cei + (0D0,1D0)*DSIGN(pi,DIMAG(Z))
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE E1XB(X,E1)
!
!       ============================================
!       Purpose: Compute exponential integral E1(x)
!       Input :  x  --- Argument of E1(x)
!       Output:  E1 --- E1(x)  ( x > 0 )
!       ============================================
!
      IMPLICIT NONE
!*--E1XB5978
      DOUBLE PRECISION E1 , ga , r , t , t0 , X
      INTEGER k , m
      IF ( X==0.0 ) THEN
         E1 = 1.0D+300
      ELSEIF ( X<=1.0 ) THEN
         E1 = 1.0D0
         r = 1.0D0
         DO k = 1 , 25
            r = -r*k*X/(k+1.0D0)**2
            E1 = E1 + r
            IF ( DABS(r)<=DABS(E1)*1.0D-15 ) GOTO 50
         ENDDO
 50      ga = 0.5772156649015328D0
         E1 = -ga - DLOG(X) + X*E1
      ELSE
         m = 20 + INT(80.0/X)
         t0 = 0.0D0
         DO k = m , 1 , -1
            t0 = k/(1.0D0+k/(X+t0))
         ENDDO
         t = 1.0D0/(X+t0)
         E1 = EXP(-X)*t
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE CHGM(A,B,X,Hg)
!
!       ===================================================
!       Purpose: Compute confluent hypergeometric function
!                M(a,b,x)
!       Input  : a  --- Parameter
!                b  --- Parameter ( b <> 0,-1,-2,... )
!                x  --- Argument
!       Output:  HG --- M(a,b,x)
!       Routine called: CGAMA for computing complex ln[Г(x)]
!       ===================================================
!
      IMPLICIT NONE
!*--CHGM6022
      DOUBLE PRECISION A , a0 , a1 , B , Hg , hg1 , hg2 , pi , r1 , r2 ,&
                     & rg , sum1 , sum2 , tai , tar , tbai , tbar ,     &
                     & tbi , tbr , X
      DOUBLE PRECISION x0 , xg , y , y0 , y1
      COMPLEX*16 cta , ctb , ctba
      INTEGER i , j , la , n , nl
      pi = 3.141592653589793D0
      a0 = A
      a1 = A
      x0 = X
      Hg = 0.0D0
!       DLMF 13.2.39
      IF ( X<0.0D0 ) THEN
         A = B - A
         a0 = A
         X = DABS(X)
      ENDIF
      nl = 0
      la = 0
      IF ( A>=2.0D0 ) THEN
!       preparing terms for DLMF 13.3.1
         nl = 1
         la = INT(A)
         A = A - la - 1.0D0
      ENDIF
      y0 = 0.0D0
      y1 = 0.0D0
      DO n = 0 , nl
         IF ( a0>=2.0D0 ) A = A + 1.0D0
         IF ( X<=30.0D0+DABS(B) .OR. A<0.0D0 ) THEN
            Hg = 1.0D0
            rg = 1.0D0
            DO j = 1 , 500
               rg = rg*(A+j-1.0D0)/(j*(B+j-1.0D0))*X
               Hg = Hg + rg
               IF ( Hg/=0D0 .AND. DABS(rg/Hg)<1.0D-15 ) THEN
!       DLMF 13.2.39 (cf. above)
                  IF ( x0<0.0D0 ) Hg = Hg*EXP(x0)
                  GOTO 50
               ENDIF
            ENDDO
         ELSE
!       DLMF 13.7.2 & 13.2.4, SUM2 corresponds to first sum
            y = 0.0D0
            CALL CGAMA(A,y,0,tar,tai)
            cta = DCMPLX(tar,tai)
            y = 0.0D0
            CALL CGAMA(B,y,0,tbr,tbi)
            ctb = DCMPLX(tbr,tbi)
            xg = B - A
            y = 0.0D0
            CALL CGAMA(xg,y,0,tbar,tbai)
            ctba = DCMPLX(tbar,tbai)
            sum1 = 1.0D0
            sum2 = 1.0D0
            r1 = 1.0D0
            r2 = 1.0D0
            DO i = 1 , 8
               r1 = -r1*(A+i-1.0D0)*(A-B+i)/(X*i)
               r2 = -r2*(B-A+i-1.0D0)*(A-i)/(X*i)
               sum1 = sum1 + r1
               sum2 = sum2 + r2
            ENDDO
            IF ( x0>=0.0D0 ) THEN
               hg1 = DBLE(EXP(ctb-ctba))*X**(-A)*DCOS(pi*A)*sum1
               hg2 = DBLE(EXP(ctb-cta+X))*X**(A-B)*sum2
            ELSE
!       DLMF 13.2.39 (cf. above)
               hg1 = DBLE(EXP(ctb-ctba+x0))*X**(-A)*DCOS(pi*A)*sum1
               hg2 = DBLE(EXP(ctb-cta))*X**(A-B)*sum2
            ENDIF
            Hg = hg1 + hg2
         ENDIF
 50      IF ( n==0 ) y0 = Hg
         IF ( n==1 ) y1 = Hg
      ENDDO
      IF ( a0>=2.0D0 ) THEN
!       DLMF 13.3.1
         DO i = 1 , la - 1
            Hg = ((2.0D0*A-B+X)*y1+(B-A)*y0)/A
            y0 = y1
            y1 = Hg
            A = A + 1.0D0
         ENDDO
      ENDIF
      A = a1
      X = x0
      END
 
!       **********************************
 
      SUBROUTINE HYGFX(A,B,C,X,Hf,Isfer)
!
!       ====================================================
!       Purpose: Compute hypergeometric function F(a,b,c,x)
!       Input :  a --- Parameter
!                b --- Parameter
!                c --- Parameter, c <> 0,-1,-2,...
!                x --- Argument   ( x < 1 )
!       Output:  HF --- F(a,b,c,x)
!                ISFER --- Error flag
!       Routines called:
!            (1) GAMMA2 for computing gamma function
!            (2) PSI_SPEC for computing psi function
!       ====================================================
!
      IMPLICIT NONE
!*--HYGFX6133
      DOUBLE PRECISION A , a0 , aa , B , bb , C , c0 , c1 , el , eps ,  &
                     & f0 , f1 , g0 , g1 , g2 , g3 , ga , gabc , gam ,  &
                     & gb
      DOUBLE PRECISION gbm , gc , gca , gcab , gcb , gm , Hf , hw , pa ,&
                     & pb , pi , r , r0 , r1 , rm , rp , sm , sp , sp0 ,&
                     & X
      DOUBLE PRECISION x1
      INTEGER Isfer , j , k , m , nm
      LOGICAL l0 , l1 , l2 , l3 , l4 , l5
      pi = 3.141592653589793D0
      el = .5772156649015329D0
      Isfer = 0
      l0 = C==INT(C) .AND. C<0.0
      l1 = 1.0D0 - X<1.0D-15 .AND. C - A - B<=0.0
      l2 = A==INT(A) .AND. A<0.0
      l3 = B==INT(B) .AND. B<0.0
      l4 = C - A==INT(C-A) .AND. C - A<=0.0
      l5 = C - B==INT(C-B) .AND. C - B<=0.0
      IF ( l0 .OR. l1 ) THEN
         Isfer = 3
         RETURN
      ENDIF
      eps = 1.0D-15
      IF ( X>0.95 ) eps = 1.0D-8
      IF ( X==0.0 .OR. A==0.0 .OR. B==0.0 ) THEN
         Hf = 1.0D0
         RETURN
      ELSEIF ( 1.0D0-X==eps .AND. C-A-B>0.0 ) THEN
         CALL GAMMA2(C,gc)
         CALL GAMMA2(C-A-B,gcab)
         CALL GAMMA2(C-A,gca)
         CALL GAMMA2(C-B,gcb)
         Hf = gc*gcab/(gca*gcb)
         RETURN
      ELSEIF ( 1.0D0+X<=eps .AND. DABS(C-A+B-1.0)<=eps ) THEN
         g0 = DSQRT(pi)*2.0D0**(-A)
         CALL GAMMA2(C,g1)
         CALL GAMMA2(1.0D0+A/2.0-B,g2)
         CALL GAMMA2(0.5D0+0.5*A,g3)
         Hf = g0*g1/(g2*g3)
         RETURN
      ELSEIF ( l2 .OR. l3 ) THEN
         IF ( l2 ) nm = INT(ABS(A))
         IF ( l3 ) nm = INT(ABS(B))
         Hf = 1.0D0
         r = 1.0D0
         DO k = 1 , nm
            r = r*(A+k-1.0D0)*(B+k-1.0D0)/(k*(C+k-1.0D0))*X
            Hf = Hf + r
         ENDDO
         RETURN
      ELSEIF ( l4 .OR. l5 ) THEN
         IF ( l4 ) nm = INT(ABS(C-A))
         IF ( l5 ) nm = INT(ABS(C-B))
         Hf = 1.0D0
         r = 1.0D0
         DO k = 1 , nm
            r = r*(C-A+k-1.0D0)*(C-B+k-1.0D0)/(k*(C+k-1.0D0))*X
            Hf = Hf + r
         ENDDO
         Hf = (1.0D0-X)**(C-A-B)*Hf
         RETURN
      ENDIF
      aa = A
      bb = B
      x1 = X
      IF ( X<0.0D0 ) THEN
         X = X/(X-1.0D0)
         IF ( C>A .AND. B<A .AND. B>0.0 ) THEN
            A = bb
            B = aa
         ENDIF
         B = C - B
      ENDIF
      hw = 0.0D0
      IF ( X>=0.75D0 ) THEN
         gm = 0.0D0
         IF ( DABS(C-A-B-INT(C-A-B))<1.0D-15 ) THEN
            m = INT(C-A-B)
            CALL GAMMA2(A,ga)
            CALL GAMMA2(B,gb)
            CALL GAMMA2(C,gc)
            CALL GAMMA2(A+m,gam)
            CALL GAMMA2(B+m,gbm)
            CALL PSI_SPEC(A,pa)
            CALL PSI_SPEC(B,pb)
            IF ( m/=0 ) gm = 1.0D0
            DO j = 1 , ABS(m) - 1
               gm = gm*j
            ENDDO
            rm = 1.0D0
            DO j = 1 , ABS(m)
               rm = rm*j
            ENDDO
            f0 = 1.0D0
            r0 = 1.0D0
            r1 = 1.0D0
            sp0 = 0.D0
            sp = 0.0D0
            IF ( m>=0 ) THEN
               c0 = gm*gc/(gam*gbm)
               c1 = -gc*(X-1.0D0)**m/(ga*gb*rm)
               DO k = 1 , m - 1
                  r0 = r0*(A+k-1.0D0)*(B+k-1.0)/(k*(k-m))*(1.0-X)
                  f0 = f0 + r0
               ENDDO
               DO k = 1 , m
                  sp0 = sp0 + 1.0D0/(A+k-1.0) + 1.0/(B+k-1.0) - 1.0/k
               ENDDO
               f1 = pa + pb + sp0 + 2.0D0*el + DLOG(1.0D0-X)
               DO k = 1 , 250
                  sp = sp + (1.0D0-A)/(k*(A+k-1.0)) + (1.0-B)           &
                     & /(k*(B+k-1.0))
                  sm = 0.0D0
                  DO j = 1 , m
                     sm = sm + (1.0D0-A)/((j+k)*(A+j+k-1.0))            &
                        & + 1.0/(B+j+k-1.0)
                  ENDDO
                  rp = pa + pb + 2.0D0*el + sp + sm + DLOG(1.0D0-X)
                  r1 = r1*(A+m+k-1.0D0)*(B+m+k-1.0)/(k*(m+k))*(1.0-X)
                  f1 = f1 + r1*rp
                  IF ( DABS(f1-hw)<DABS(f1)*eps ) GOTO 10
                  hw = f1
               ENDDO
 10            Hf = f0*c0 + f1*c1
            ELSEIF ( m<0 ) THEN
               m = -m
               c0 = gm*gc/(ga*gb*(1.0D0-X)**m)
               c1 = -(-1)**m*gc/(gam*gbm*rm)
               DO k = 1 , m - 1
                  r0 = r0*(A-m+k-1.0D0)*(B-m+k-1.0)/(k*(k-m))*(1.0-X)
                  f0 = f0 + r0
               ENDDO
               DO k = 1 , m
                  sp0 = sp0 + 1.0D0/k
               ENDDO
               f1 = pa + pb - sp0 + 2.0D0*el + DLOG(1.0D0-X)
               DO k = 1 , 250
                  sp = sp + (1.0D0-A)/(k*(A+k-1.0)) + (1.0-B)           &
                     & /(k*(B+k-1.0))
                  sm = 0.0D0
                  DO j = 1 , m
                     sm = sm + 1.0D0/(j+k)
                  ENDDO
                  rp = pa + pb + 2.0D0*el + sp - sm + DLOG(1.0D0-X)
                  r1 = r1*(A+k-1.0D0)*(B+k-1.0)/(k*(m+k))*(1.0-X)
                  f1 = f1 + r1*rp
                  IF ( DABS(f1-hw)<DABS(f1)*eps ) GOTO 20
                  hw = f1
               ENDDO
 20            Hf = f0*c0 + f1*c1
            ENDIF
         ELSE
            CALL GAMMA2(A,ga)
            CALL GAMMA2(B,gb)
            CALL GAMMA2(C,gc)
            CALL GAMMA2(C-A,gca)
            CALL GAMMA2(C-B,gcb)
            CALL GAMMA2(C-A-B,gcab)
            CALL GAMMA2(A+B-C,gabc)
            c0 = gc*gcab/(gca*gcb)
            c1 = gc*gabc/(ga*gb)*(1.0D0-X)**(C-A-B)
            Hf = 0.0D0
            r0 = c0
            r1 = c1
            DO k = 1 , 250
               r0 = r0*(A+k-1.0D0)*(B+k-1.0)/(k*(A+B-C+k))*(1.0-X)
               r1 = r1*(C-A+k-1.0D0)*(C-B+k-1.0)/(k*(C-A-B+k))*(1.0-X)
               Hf = Hf + r0 + r1
               IF ( DABS(Hf-hw)<DABS(Hf)*eps ) GOTO 40
               hw = Hf
            ENDDO
 40         Hf = Hf + c0 + c1
         ENDIF
      ELSE
         a0 = 1.0D0
         IF ( C>A .AND. C<2.0D0*A .AND. C>B .AND. C<2.0D0*B ) THEN
            a0 = (1.0D0-X)**(C-A-B)
            A = C - A
            B = C - B
         ENDIF
         Hf = 1.0D0
         r = 1.0D0
         DO k = 1 , 250
            r = r*(A+k-1.0D0)*(B+k-1.0D0)/(k*(C+k-1.0D0))*X
            Hf = Hf + r
            IF ( DABS(Hf-hw)<=DABS(Hf)*eps ) GOTO 50
            hw = Hf
         ENDDO
 50      Hf = a0*Hf
      ENDIF
      IF ( x1<0.0D0 ) THEN
         X = x1
         c0 = 1.0D0/(1.0D0-X)**aa
         Hf = c0*Hf
      ENDIF
      A = aa
      B = bb
      IF ( k>120 ) Isfer = 5
      END
 
 
 
!       **********************************
 
      SUBROUTINE CCHG(A,B,Z,Chg)
!
!       ===================================================
!       Purpose: Compute confluent hypergeometric function
!                M(a,b,z) with real parameters a, b and a
!                complex argument z
!       Input :  a --- Parameter
!                b --- Parameter
!                z --- Complex argument
!       Output:  CHG --- M(a,b,z)
!       Routine called: CGAMA for computing complex ln[Г(x)]
!       ===================================================
!
      IMPLICIT NONE
!*--CCHG6356
      DOUBLE PRECISION A , a0 , a1 , B , ba , g1i , g1r , g2i , g2r ,   &
                     & g3i , g3r , phi , pi , x , x0 , y
      COMPLEX*16 cfac , cg1 , cg2 , cg3 , Chg , chg1 , chg2 , chw , ci ,&
               & cr , cr1 , cr2 , crg , cs1 , cs2 , cy0 , cy1 , Z , z0
      INTEGER i , j , k , la , m , n , nl , ns
      pi = 3.141592653589793D0
      ci = (0.0D0,1.0D0)
      a0 = A
      a1 = A
      z0 = Z
      IF ( B==0.0 .OR. B==-INT(ABS(B)) ) THEN
         Chg = (1.0D+300,0.0D0)
      ELSEIF ( A==0.0D0 .OR. Z==0.0D0 ) THEN
         Chg = (1.0D0,0.0D0)
      ELSEIF ( A==-1.0D0 ) THEN
         Chg = 1.0D0 - Z/B
      ELSEIF ( A==B ) THEN
         Chg = EXP(Z)
      ELSEIF ( A-B==1.0D0 ) THEN
         Chg = (1.0D0+Z/B)*EXP(Z)
      ELSEIF ( A==1.0D0 .AND. B==2.0D0 ) THEN
         Chg = (EXP(Z)-1.0D0)/Z
      ELSEIF ( A==INT(A) .AND. A<0.0D0 ) THEN
         m = INT(-A)
         cr = (1.0D0,0.0D0)
         Chg = (1.0D0,0.0D0)
         DO k = 1 , m
            cr = cr*(A+k-1.0D0)/k/(B+k-1.0D0)*Z
            Chg = Chg + cr
         ENDDO
      ELSE
         x0 = DBLE(Z)
         IF ( x0<0.0D0 ) THEN
            A = B - A
            a0 = A
            Z = -Z
         ENDIF
         nl = 0
         la = 0
         IF ( A>=2.0D0 ) THEN
            nl = 1
            la = INT(A)
            A = A - la - 1.0D0
         ENDIF
         ns = 0
         DO n = 0 , nl
            IF ( a0>=2.0D0 ) A = A + 1.0D0
            IF ( ABS(Z)<20.0D0+ABS(B) .OR. A<0.0D0 ) THEN
               Chg = (1.0D0,0.0D0)
               crg = (1.0D0,0.0D0)
               DO j = 1 , 500
                  crg = crg*(A+j-1.0D0)/(j*(B+j-1.0D0))*Z
                  Chg = Chg + crg
                  IF ( ABS((Chg-chw)/Chg)<1.D-15 ) GOTO 20
                  chw = Chg
               ENDDO
            ELSE
               y = 0.0D0
               CALL CGAMA(A,y,0,g1r,g1i)
               cg1 = DCMPLX(g1r,g1i)
               y = 0.0D0
               CALL CGAMA(B,y,0,g2r,g2i)
               cg2 = DCMPLX(g2r,g2i)
               ba = B - A
               y = 0.0D0
               CALL CGAMA(ba,y,0,g3r,g3i)
               cg3 = DCMPLX(g3r,g3i)
               cs1 = (1.0D0,0.0D0)
               cs2 = (1.0D0,0.0D0)
               cr1 = (1.0D0,0.0D0)
               cr2 = (1.0D0,0.0D0)
               DO i = 1 , 8
                  cr1 = -cr1*(A+i-1.0D0)*(A-B+i)/(Z*i)
                  cr2 = cr2*(B-A+i-1.0D0)*(i-A)/(Z*i)
                  cs1 = cs1 + cr1
                  cs2 = cs2 + cr2
               ENDDO
               x = DBLE(Z)
               y = DIMAG(Z)
               IF ( x==0.0 .AND. y>=0.0 ) THEN
                  phi = 0.5D0*pi
               ELSEIF ( x==0.0 .AND. y<=0.0 ) THEN
                  phi = -0.5D0*pi
               ELSE
                  phi = DATAN(y/x)
               ENDIF
               IF ( phi>-0.5*pi .AND. phi<1.5*pi ) ns = 1
               IF ( phi>-1.5*pi .AND. phi<=-0.5*pi ) ns = -1
               cfac = EXP(ns*ci*pi*A)
               IF ( y==0.0D0 ) cfac = DCOS(pi*A)
               chg1 = EXP(cg2-cg3)*Z**(-A)*cfac*cs1
               chg2 = EXP(cg2-cg1+Z)*Z**(A-B)*cs2
               Chg = chg1 + chg2
            ENDIF
 20         IF ( n==0 ) cy0 = Chg
            IF ( n==1 ) cy1 = Chg
         ENDDO
         IF ( a0>=2.0D0 ) THEN
            DO i = 1 , la - 1
               Chg = ((2.0D0*A-B+Z)*cy1+(B-A)*cy0)/A
               cy0 = cy1
               cy1 = Chg
               A = A + 1.0D0
            ENDDO
         ENDIF
         IF ( x0<0.0D0 ) Chg = Chg*EXP(-Z)
      ENDIF
      A = a1
      Z = z0
      END
 
 
 
!       **********************************
 
      SUBROUTINE HYGFZ(A,B,C,Z,Zhf,Isfer)
!
!       ======================================================
!       Purpose: Compute the hypergeometric function for a
!                complex argument, F(a,b,c,z)
!       Input :  a --- Parameter
!                b --- Parameter
!                c --- Parameter,  c <> 0,-1,-2,...
!                z --- Complex argument
!       Output:  ZHF --- F(a,b,c,z)
!                ISFER --- Error flag
!       Routines called:
!            (1) GAMMA2 for computing gamma function
!            (2) PSI_SPEC for computing psi function
!       ======================================================
!
      IMPLICIT NONE
!*--HYGFZ6492
      DOUBLE PRECISION A , a0 , aa , B , bb , C , ca , cb , el , eps ,  &
                     & g0 , g1 , g2 , g3 , ga , gab , gabc , gam , gb , &
                     & gba
      DOUBLE PRECISION gbm , gc , gca , gcab , gcb , gcbk , gm , pa ,   &
                     & pac , pb , pca , pi , rk1 , rk2 , rm , sj1 ,     &
                     & sj2 , sm , sp , sp0
      DOUBLE PRECISION sq , t0 , w0 , ws , x , y
      INTEGER Isfer , j , k , m , mab , mcab , nca , ncb , nm
      COMPLEX*16 Z , z00 , z1 , zc0 , zc1 , zf0 , zf1 , Zhf , zp , zp0 ,&
               & zr , zr0 , zr1 , zw
      LOGICAL l0 , l1 , l2 , l3 , l4 , l5 , l6
      x = DBLE(Z)
      y = DIMAG(Z)
      eps = 1.0D-15
      Isfer = 0
      l0 = C==INT(C) .AND. C<0.0D0
      l1 = DABS(1.0D0-x)<eps .AND. y==0.0D0 .AND. C - A - B<=0.0D0
      l2 = ABS(Z+1.0D0)<eps .AND. DABS(C-A+B-1.0D0)<eps
      l3 = A==INT(A) .AND. A<0.0D0
      l4 = B==INT(B) .AND. B<0.0D0
      l5 = C - A==INT(C-A) .AND. C - A<=0.0D0
      l6 = C - B==INT(C-B) .AND. C - B<=0.0D0
      aa = A
      bb = B
      a0 = ABS(Z)
      IF ( a0>0.95D0 ) eps = 1.0D-8
      pi = 3.141592653589793D0
      el = .5772156649015329D0
      IF ( l0 .OR. l1 ) THEN
         Isfer = 3
         RETURN
      ENDIF
      nm = 0
      IF ( a0==0.0D0 .OR. A==0.0D0 .OR. B==0.0D0 ) THEN
         Zhf = (1.0D0,0.0D0)
      ELSEIF ( Z==1.0D0 .AND. C-A-B>0.0D0 ) THEN
         CALL GAMMA2(C,gc)
         CALL GAMMA2(C-A-B,gcab)
         CALL GAMMA2(C-A,gca)
         CALL GAMMA2(C-B,gcb)
         Zhf = gc*gcab/(gca*gcb)
      ELSEIF ( l2 ) THEN
         g0 = DSQRT(pi)*2.0D0**(-A)
         CALL GAMMA2(C,g1)
         CALL GAMMA2(1.0D0+A/2.0D0-B,g2)
         CALL GAMMA2(0.5D0+0.5D0*A,g3)
         Zhf = g0*g1/(g2*g3)
      ELSEIF ( l3 .OR. l4 ) THEN
         IF ( l3 ) nm = INT(ABS(A))
         IF ( l4 ) nm = INT(ABS(B))
         Zhf = (1.0D0,0.0D0)
         zr = (1.0D0,0.0D0)
         DO k = 1 , nm
            zr = zr*(A+k-1.0D0)*(B+k-1.0D0)/(k*(C+k-1.0D0))*Z
            Zhf = Zhf + zr
         ENDDO
      ELSEIF ( l5 .OR. l6 ) THEN
         IF ( l5 ) nm = INT(ABS(C-A))
         IF ( l6 ) nm = INT(ABS(C-B))
         Zhf = (1.0D0,0.0D0)
         zr = (1.0D0,0.0D0)
         DO k = 1 , nm
            zr = zr*(C-A+k-1.0D0)*(C-B+k-1.0D0)/(k*(C+k-1.0D0))*Z
            Zhf = Zhf + zr
         ENDDO
         Zhf = (1.0D0-Z)**(C-A-B)*Zhf
      ELSEIF ( a0<=1.0D0 ) THEN
         IF ( x<0.0D0 ) THEN
            z1 = Z/(Z-1.0D0)
            IF ( C>A .AND. B<A .AND. B>0.0 ) THEN
               A = bb
               B = aa
            ENDIF
            zc0 = 1.0D0/((1.0D0-Z)**A)
            Zhf = (1.0D0,0.0D0)
            zr0 = (1.0D0,0.0D0)
            DO k = 1 , 500
               zr0 = zr0*(A+k-1.0D0)*(C-B+k-1.0D0)/(k*(C+k-1.0D0))*z1
               Zhf = Zhf + zr0
               IF ( ABS(Zhf-zw)<ABS(Zhf)*eps ) GOTO 20
               zw = Zhf
            ENDDO
 20         Zhf = zc0*Zhf
         ELSEIF ( a0>=0.90D0 ) THEN
            gm = 0.0D0
            mcab = INT(C-A-B+eps*DSIGN(1.0D0,C-A-B))
            IF ( DABS(C-A-B-mcab)<eps ) THEN
               m = INT(C-A-B)
               CALL GAMMA2(A,ga)
               CALL GAMMA2(B,gb)
               CALL GAMMA2(C,gc)
               CALL GAMMA2(A+m,gam)
               CALL GAMMA2(B+m,gbm)
               CALL PSI_SPEC(A,pa)
               CALL PSI_SPEC(B,pb)
               IF ( m/=0 ) gm = 1.0D0
               DO j = 1 , ABS(m) - 1
                  gm = gm*j
               ENDDO
               rm = 1.0D0
               DO j = 1 , ABS(m)
                  rm = rm*j
               ENDDO
               zf0 = (1.0D0,0.0D0)
               zr0 = (1.0D0,0.0D0)
               zr1 = (1.0D0,0.0D0)
               sp0 = 0.D0
               sp = 0.0D0
               IF ( m>=0 ) THEN
                  zc0 = gm*gc/(gam*gbm)
                  zc1 = -gc*(Z-1.0D0)**m/(ga*gb*rm)
                  DO k = 1 , m - 1
                     zr0 = zr0*(A+k-1.D0)*(B+k-1.D0)/(k*(k-m))*(1.D0-Z)
                     zf0 = zf0 + zr0
                  ENDDO
                  DO k = 1 , m
                     sp0 = sp0 + 1.0D0/(A+k-1.0D0) + 1.0/(B+k-1.0D0)    &
                         & - 1.D0/k
                  ENDDO
                  zf1 = pa + pb + sp0 + 2.0D0*el + LOG(1.0D0-Z)
                  DO k = 1 , 500
                     sp = sp + (1.0D0-A)/(k*(A+k-1.0D0)) + (1.0D0-B)    &
                        & /(k*(B+k-1.0D0))
                     sm = 0.0D0
                     DO j = 1 , m
                        sm = sm + (1.0D0-A)/((j+k)*(A+j+k-1.0D0))       &
                           & + 1.0D0/(B+j+k-1.0D0)
                     ENDDO
                     zp = pa + pb + 2.0D0*el + sp + sm + LOG(1.0D0-Z)
                     zr1 = zr1*(A+m+k-1.0D0)*(B+m+k-1.0D0)/(k*(m+k))    &
                         & *(1.0D0-Z)
                     zf1 = zf1 + zr1*zp
                     IF ( ABS(zf1-zw)<ABS(zf1)*eps ) GOTO 25
                     zw = zf1
                  ENDDO
 25               Zhf = zf0*zc0 + zf1*zc1
               ELSEIF ( m<0 ) THEN
                  m = -m
                  zc0 = gm*gc/(ga*gb*(1.0D0-Z)**m)
                  zc1 = -(-1)**m*gc/(gam*gbm*rm)
                  DO k = 1 , m - 1
                     zr0 = zr0*(A-m+k-1.0D0)*(B-m+k-1.0D0)/(k*(k-m))    &
                         & *(1.0D0-Z)
                     zf0 = zf0 + zr0
                  ENDDO
                  DO k = 1 , m
                     sp0 = sp0 + 1.0D0/k
                  ENDDO
                  zf1 = pa + pb - sp0 + 2.0D0*el + LOG(1.0D0-Z)
                  DO k = 1 , 500
                     sp = sp + (1.0D0-A)/(k*(A+k-1.0D0)) + (1.0D0-B)    &
                        & /(k*(B+k-1.0D0))
                     sm = 0.0D0
                     DO j = 1 , m
                        sm = sm + 1.0D0/(j+k)
                     ENDDO
                     zp = pa + pb + 2.0D0*el + sp - sm + LOG(1.0D0-Z)
                     zr1 = zr1*(A+k-1.D0)*(B+k-1.D0)/(k*(m+k))*(1.D0-Z)
                     zf1 = zf1 + zr1*zp
                     IF ( ABS(zf1-zw)<ABS(zf1)*eps ) GOTO 30
                     zw = zf1
                  ENDDO
 30               Zhf = zf0*zc0 + zf1*zc1
               ENDIF
            ELSE
               CALL GAMMA2(A,ga)
               CALL GAMMA2(B,gb)
               CALL GAMMA2(C,gc)
               CALL GAMMA2(C-A,gca)
               CALL GAMMA2(C-B,gcb)
               CALL GAMMA2(C-A-B,gcab)
               CALL GAMMA2(A+B-C,gabc)
               zc0 = gc*gcab/(gca*gcb)
               zc1 = gc*gabc/(ga*gb)*(1.0D0-Z)**(C-A-B)
               Zhf = (0.0D0,0.0D0)
               zr0 = zc0
               zr1 = zc1
               DO k = 1 , 500
                  zr0 = zr0*(A+k-1.D0)*(B+k-1.D0)/(k*(A+B-C+k))*(1.D0-Z)
                  zr1 = zr1*(C-A+k-1.0D0)*(C-B+k-1.0D0)/(k*(C-A-B+k))   &
                      & *(1.0D0-Z)
                  Zhf = Zhf + zr0 + zr1
                  IF ( ABS(Zhf-zw)<ABS(Zhf)*eps ) GOTO 40
                  zw = Zhf
               ENDDO
 40            Zhf = Zhf + zc0 + zc1
            ENDIF
         ELSE
            z00 = (1.0D0,0.0D0)
            IF ( C-A<A .AND. C-B<B ) THEN
               z00 = (1.0D0-Z)**(C-A-B)
               A = C - A
               B = C - B
            ENDIF
            Zhf = (1.0D0,0.D0)
            zr = (1.0D0,0.0D0)
            DO k = 1 , 1500
               zr = zr*(A+k-1.0D0)*(B+k-1.0D0)/(k*(C+k-1.0D0))*Z
               Zhf = Zhf + zr
               IF ( ABS(Zhf-zw)<=ABS(Zhf)*eps ) GOTO 60
               zw = Zhf
            ENDDO
 60         Zhf = z00*Zhf
         ENDIF
      ELSEIF ( a0>1.0D0 ) THEN
         mab = INT(A-B+eps*DSIGN(1.0D0,A-B))
         IF ( DABS(A-B-mab)<eps .AND. a0<=1.1D0 ) B = B + eps
         IF ( DABS(A-B-mab)>eps ) THEN
            CALL GAMMA2(A,ga)
            CALL GAMMA2(B,gb)
            CALL GAMMA2(C,gc)
            CALL GAMMA2(A-B,gab)
            CALL GAMMA2(B-A,gba)
            CALL GAMMA2(C-A,gca)
            CALL GAMMA2(C-B,gcb)
            zc0 = gc*gba/(gca*gb*(-Z)**A)
            zc1 = gc*gab/(gcb*ga*(-Z)**B)
            zr0 = zc0
            zr1 = zc1
            Zhf = (0.0D0,0.0D0)
            DO k = 1 , 500
               zr0 = zr0*(A+k-1.0D0)*(A-C+k)/((A-B+k)*k*Z)
               zr1 = zr1*(B+k-1.0D0)*(B-C+k)/((B-A+k)*k*Z)
               Zhf = Zhf + zr0 + zr1
               IF ( ABS((Zhf-zw)/Zhf)<=eps ) GOTO 80
               zw = Zhf
            ENDDO
 80         Zhf = Zhf + zc0 + zc1
         ELSE
            IF ( A-B<0.0D0 ) THEN
               A = bb
               B = aa
            ENDIF
            ca = C - A
            cb = C - B
            nca = INT(ca+eps*DSIGN(1.0D0,ca))
            ncb = INT(cb+eps*DSIGN(1.0D0,cb))
            IF ( DABS(ca-nca)<eps .OR. DABS(cb-ncb)<eps ) C = C + eps
            CALL GAMMA2(A,ga)
            CALL GAMMA2(C,gc)
            CALL GAMMA2(C-B,gcb)
            CALL PSI_SPEC(A,pa)
            CALL PSI_SPEC(C-A,pca)
            CALL PSI_SPEC(A-C,pac)
            mab = INT(A-B+eps)
            zc0 = gc/(ga*(-Z)**B)
            CALL GAMMA2(A-B,gm)
            zf0 = gm/gcb*zc0
            zr = zc0
            DO k = 1 , mab - 1
               zr = zr*(B+k-1.0D0)/(k*Z)
               t0 = A - B - k
               CALL GAMMA2(t0,g0)
               CALL GAMMA2(C-B-k,gcbk)
               zf0 = zf0 + zr*g0/gcbk
            ENDDO
            IF ( mab==0 ) zf0 = (0.0D0,0.0D0)
            zc1 = gc/(ga*gcb*(-Z)**A)
            sp = -2.0D0*el - pa - pca
            DO j = 1 , mab
               sp = sp + 1.0D0/j
            ENDDO
            zp0 = sp + LOG(-Z)
            sq = 1.0D0
            DO j = 1 , mab
               sq = sq*(B+j-1.0D0)*(B-C+j)/j
            ENDDO
            zf1 = (sq*zp0)*zc1
            zr = zc1
            rk1 = 1.0D0
            sj1 = 0.0D0
            w0 = 0.0D0
            DO k = 1 , 10000
               zr = zr/Z
               rk1 = rk1*(B+k-1.0D0)*(B-C+k)/(k*k)
               rk2 = rk1
               DO j = k + 1 , k + mab
                  rk2 = rk2*(B+j-1.0D0)*(B-C+j)/j
               ENDDO
               sj1 = sj1 + (A-1.0D0)/(k*(A+k-1.0D0)) + (A-C-1.0D0)      &
                   & /(k*(A-C+k-1.0D0))
               sj2 = sj1
               DO j = k + 1 , k + mab
                  sj2 = sj2 + 1.0D0/j
               ENDDO
               zp = -2.0D0*el - pa - pac + sj2 - 1.0D0/(k+A-C)          &
                  & - pi/DTAN(pi*(k+A-C)) + LOG(-Z)
               zf1 = zf1 + rk2*zr*zp
               ws = ABS(zf1)
               IF ( DABS((ws-w0)/ws)<eps ) GOTO 100
               w0 = ws
            ENDDO
 100        Zhf = zf0 + zf1
         ENDIF
      ENDIF
      A = aa
      B = bb
      IF ( k>150 ) Isfer = 5
      END
 
 
 
!       **********************************
 
      SUBROUTINE ITAIRY(X,Apt,Bpt,Ant,Bnt)
!
!       ======================================================
!       Purpose: Compute the integrals of Airy fnctions with
!                respect to t from 0 and x ( x ≥ 0 )
!       Input  : x   --- Upper limit of the integral
!       Output : APT --- Integration of Ai(t) from 0 and x
!                BPT --- Integration of Bi(t) from 0 and x
!                ANT --- Integration of Ai(-t) from 0 and x
!                BNT --- Integration of Bi(-t) from 0 and x
!       ======================================================
!
      IMPLICIT NONE
!*--ITAIRY6813
      DOUBLE PRECISION a , Ant , Apt , Bnt , Bpt , c1 , c2 , eps , fx , &
                     & gx , pi , q0 , q1 , q2 , r , sr3 , su1 , su2 ,   &
                     & su3 , su4
      DOUBLE PRECISION su5 , su6 , X , xe , xp6 , xr1 , xr2
      INTEGER k , l
      DIMENSION a(16)
      eps = 1.0D-15
      pi = 3.141592653589793D0
      c1 = .355028053887817D0
      c2 = .258819403792807D0
      sr3 = 1.732050807568877D0
      IF ( X==0.0D0 ) THEN
         Apt = 0.0D0
         Bpt = 0.0D0
         Ant = 0.0D0
         Bnt = 0.0D0
      ELSEIF ( DABS(X)<=9.25D0 ) THEN
         DO l = 0 , 1
            X = (-1)**l*X
            fx = X
            r = X
            DO k = 1 , 40
               r = r*(3.0*k-2.0D0)/(3.0*k+1.0D0)*X/(3.0*k)              &
                 & *X/(3.0*k-1.0D0)*X
               fx = fx + r
               IF ( DABS(r)<DABS(fx)*eps ) GOTO 20
            ENDDO
 20         gx = .5D0*X*X
            r = gx
            DO k = 1 , 40
               r = r*(3.0*k-1.0D0)/(3.0*k+2.0D0)*X/(3.0*k)              &
                 & *X/(3.0*k+1.0D0)*X
               gx = gx + r
               IF ( DABS(r)<DABS(gx)*eps ) GOTO 40
            ENDDO
 40         Ant = c1*fx - c2*gx
            Bnt = sr3*(c1*fx+c2*gx)
            IF ( l==0 ) THEN
               Apt = Ant
               Bpt = Bnt
            ELSE
               Ant = -Ant
               Bnt = -Bnt
               X = -X
            ENDIF
         ENDDO
      ELSE
         DATA a/.569444444444444D0 , .891300154320988D0 ,               &
            & .226624344493027D+01 , .798950124766861D+01 ,             &
            & .360688546785343D+02 , .198670292131169D+03 ,             &
            & .129223456582211D+04 , .969483869669600D+04 ,             &
            & .824184704952483D+05 , .783031092490225D+06 ,             &
            & .822210493622814D+07 , .945557399360556D+08 ,             &
            & .118195595640730D+10 , .159564653040121D+11 ,             &
            & .231369166433050D+12 , .358622522796969D+13/
         q2 = 1.414213562373095D0
         q0 = .3333333333333333D0
         q1 = .6666666666666667D0
         xe = X*DSQRT(X)/1.5D0
         xp6 = 1.0D0/DSQRT(6.0D0*pi*xe)
         su1 = 1.0D0
         r = 1.0D0
         xr1 = 1.0D0/xe
         DO k = 1 , 16
            r = -r*xr1
            su1 = su1 + a(k)*r
         ENDDO
         su2 = 1.0D0
         r = 1.0D0
         DO k = 1 , 16
            r = r*xr1
            su2 = su2 + a(k)*r
         ENDDO
         Apt = q0 - EXP(-xe)*xp6*su1
         Bpt = 2.0D0*EXP(xe)*xp6*su2
         su3 = 1.0D0
         r = 1.0D0
         xr2 = 1.0D0/(xe*xe)
         DO k = 1 , 8
            r = -r*xr2
            su3 = su3 + a(2*k)*r
         ENDDO
         su4 = a(1)*xr1
         r = xr1
         DO k = 1 , 7
            r = -r*xr2
            su4 = su4 + a(2*k+1)*r
         ENDDO
         su5 = su3 + su4
         su6 = su3 - su4
         Ant = q1 - q2*xp6*(su5*DCOS(xe)-su6*DSIN(xe))
         Bnt = q2*xp6*(su5*DSIN(xe)+su6*DCOS(xe))
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE IKNA(N,X,Nm,Bi,Di,Bk,Dk)
!
!       ========================================================
!       Purpose: Compute modified Bessel functions In(x) and
!                Kn(x), and their derivatives
!       Input:   x --- Argument of In(x) and Kn(x) ( x ≥ 0 )
!                n --- Order of In(x) and Kn(x)
!       Output:  BI(n) --- In(x)
!                DI(n) --- In'(x)
!                BK(n) --- Kn(x)
!                DK(n) --- Kn'(x)
!                NM --- Highest order computed
!       Routines called:
!            (1) IK01A for computing I0(x),I1(x),K0(x) & K1(x)
!            (2) MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ========================================================
!
      IMPLICIT NONE
!*--IKNA6933
      DOUBLE PRECISION Bi , bi0 , bi1 , Bk , bk0 , bk1 , Di , di0 ,     &
                     & di1 , Dk , dk0 , dk1 , f , f0 , f1 , g , g0 ,    &
                     & g1 , h , h0
      DOUBLE PRECISION h1 , s0 , X
      INTEGER k , m , MSTA1 , MSTA2 , N , Nm
      DIMENSION Bi(0:N) , Di(0:N) , Bk(0:N) , Dk(0:N)
      Nm = N
      IF ( X<=1.0D-100 ) THEN
         DO k = 0 , N
            Bi(k) = 0.0D0
            Di(k) = 0.0D0
            Bk(k) = 1.0D+300
            Dk(k) = -1.0D+300
         ENDDO
         Bi(0) = 1.0D0
         Di(1) = 0.5D0
         RETURN
      ENDIF
      CALL IK01A(X,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1)
      Bi(0) = bi0
      Bi(1) = bi1
      Bk(0) = bk0
      Bk(1) = bk1
      Di(0) = di0
      Di(1) = di1
      Dk(0) = dk0
      Dk(1) = dk1
      IF ( N<=1 ) RETURN
      IF ( X>40.0 .AND. N<INT(0.25*X) ) THEN
         h0 = bi0
         h1 = bi1
         DO k = 2 , N
            h = -2.0D0*(k-1.0D0)/X*h1 + h0
            Bi(k) = h
            h0 = h1
            h1 = h
         ENDDO
      ELSE
         m = MSTA1(X,200)
         IF ( m<N ) THEN
            Nm = m
         ELSE
            m = MSTA2(X,N,15)
         ENDIF
         f0 = 0.0D0
         f1 = 1.0D-100
         f = 0.0D0
         DO k = m , 0 , -1
            f = 2.0D0*(k+1.0D0)*f1/X + f0
            IF ( k<=Nm ) Bi(k) = f
            f0 = f1
            f1 = f
         ENDDO
         s0 = bi0/f
         DO k = 0 , Nm
            Bi(k) = s0*Bi(k)
         ENDDO
      ENDIF
      g0 = bk0
      g1 = bk1
      DO k = 2 , Nm
         g = 2.0D0*(k-1.0D0)/X*g1 + g0
         Bk(k) = g
         g0 = g1
         g1 = g
      ENDDO
      DO k = 2 , Nm
         Di(k) = Bi(k-1) - k/X*Bi(k)
         Dk(k) = -Bk(k-1) - k/X*Bk(k)
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE CJYNB(N,Z,Nm,Cbj,Cdj,Cby,Cdy)
!
!       =======================================================
!       Purpose: Compute Bessel functions Jn(z), Yn(z) and
!                their derivatives for a complex argument
!       Input :  z --- Complex argument of Jn(z) and Yn(z)
!                n --- Order of Jn(z) and Yn(z)
!       Output:  CBJ(n) --- Jn(z)
!                CDJ(n) --- Jn'(z)
!                CBY(n) --- Yn(z)
!                CDY(n) --- Yn'(z)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 to calculate the starting
!                point for backward recurrence
!       =======================================================
!
      IMPLICIT NONE
!*--CJYNB7031
      DOUBLE PRECISION a , a0 , a1 , b , b1 , el , pi , r2p , y0
      COMPLEX*16 Cbj , cbj0 , cbj1 , cbjk , cbs , Cby , cby0 , cby1 ,   &
               & Cdj , Cdy , ce , cf , cf1 , cf2 , cp0 , cp1 , cq0 ,    &
               & cq1 , cs0 , csu
      COMPLEX*16 csv , ct1 , ct2 , cu , cyy , Z
      INTEGER k , m , MSTA1 , MSTA2 , N , Nm
      DIMENSION Cbj(0:N) , Cdj(0:N) , Cby(0:N) , Cdy(0:N) , a(4) ,      &
              & b(4) , a1(4) , b1(4)
      el = 0.5772156649015329D0
      pi = 3.141592653589793D0
      r2p = .63661977236758D0
      y0 = DABS(DIMAG(Z))
      a0 = ABS(Z)
      Nm = N
      IF ( a0<1.0D-100 ) THEN
         DO k = 0 , N
            Cbj(k) = (0.0D0,0.0D0)
            Cdj(k) = (0.0D0,0.0D0)
            Cby(k) = -(1.0D+300,0.0D0)
            Cdy(k) = (1.0D+300,0.0D0)
         ENDDO
         Cbj(0) = (1.0D0,0.0D0)
         Cdj(1) = (0.5D0,0.0D0)
         RETURN
      ENDIF
      IF ( a0<=300.D0 .OR. N>80 ) THEN
         IF ( N==0 ) Nm = 1
         m = MSTA1(a0,200)
         IF ( m<Nm ) THEN
            Nm = m
         ELSE
            m = MSTA2(a0,Nm,15)
         ENDIF
         cbs = (0.0D0,0.0D0)
         csu = (0.0D0,0.0D0)
         csv = (0.0D0,0.0D0)
         cf2 = (0.0D0,0.0D0)
         cf1 = (1.0D-100,0.0D0)
         DO k = m , 0 , -1
            cf = 2.0D0*(k+1.0D0)/Z*cf1 - cf2
            IF ( k<=Nm ) Cbj(k) = cf
            IF ( k==2*INT(k/2) .AND. k/=0 ) THEN
               IF ( y0<=1.0D0 ) THEN
                  cbs = cbs + 2.0D0*cf
               ELSE
                  cbs = cbs + (-1)**(k/2)*2.0D0*cf
               ENDIF
               csu = csu + (-1)**(k/2)*cf/k
            ELSEIF ( k>1 ) THEN
               csv = csv + (-1)**(k/2)*k/(k*k-1.0D0)*cf
            ENDIF
            cf2 = cf1
            cf1 = cf
         ENDDO
         IF ( y0<=1.0D0 ) THEN
            cs0 = cbs + cf
         ELSE
            cs0 = (cbs+cf)/COS(Z)
         ENDIF
         DO k = 0 , Nm
            Cbj(k) = Cbj(k)/cs0
         ENDDO
         ce = LOG(Z/2.0D0) + el
         Cby(0) = r2p*(ce*Cbj(0)-4.0D0*csu/cs0)
         Cby(1) = r2p*(-Cbj(0)/Z+(ce-1.0D0)*Cbj(1)-4.0D0*csv/cs0)
      ELSE
         DATA a/ - .7031250000000000D-01 , .1121520996093750D+00 ,      &
            & -.5725014209747314D+00 , .6074042001273483D+01/
         DATA b/.7324218750000000D-01 , -.2271080017089844D+00 ,        &
            & .1727727502584457D+01 , -.2438052969955606D+02/
         DATA a1/.1171875000000000D+00 , -.1441955566406250D+00 ,       &
            & .6765925884246826D+00 , -.6883914268109947D+01/
         DATA b1/ - .1025390625000000D+00 , .2775764465332031D+00 ,     &
            & -.1993531733751297D+01 , .2724882731126854D+02/
         ct1 = Z - 0.25D0*pi
         cp0 = (1.0D0,0.0D0)
         DO k = 1 , 4
            cp0 = cp0 + a(k)*Z**(-2*k)
         ENDDO
         cq0 = -0.125D0/Z
         DO k = 1 , 4
            cq0 = cq0 + b(k)*Z**(-2*k-1)
         ENDDO
         cu = SQRT(r2p/Z)
         cbj0 = cu*(cp0*COS(ct1)-cq0*SIN(ct1))
         cby0 = cu*(cp0*SIN(ct1)+cq0*COS(ct1))
         Cbj(0) = cbj0
         Cby(0) = cby0
         ct2 = Z - 0.75D0*pi
         cp1 = (1.0D0,0.0D0)
         DO k = 1 , 4
            cp1 = cp1 + a1(k)*Z**(-2*k)
         ENDDO
         cq1 = 0.375D0/Z
         DO k = 1 , 4
            cq1 = cq1 + b1(k)*Z**(-2*k-1)
         ENDDO
         cbj1 = cu*(cp1*COS(ct2)-cq1*SIN(ct2))
         cby1 = cu*(cp1*SIN(ct2)+cq1*COS(ct2))
         Cbj(1) = cbj1
         Cby(1) = cby1
         DO k = 2 , Nm
            cbjk = 2.0D0*(k-1.0D0)/Z*cbj1 - cbj0
            Cbj(k) = cbjk
            cbj0 = cbj1
            cbj1 = cbjk
         ENDDO
      ENDIF
      Cdj(0) = -Cbj(1)
      DO k = 1 , Nm
         Cdj(k) = Cbj(k-1) - k/Z*Cbj(k)
      ENDDO
      IF ( ABS(Cbj(0))>1.0D0 ) Cby(1) = (Cbj(1)*Cby(0)-2.0D0/(pi*Z))    &
                                      & /Cbj(0)
      DO k = 2 , Nm
         IF ( ABS(Cbj(k-1))>=ABS(Cbj(k-2)) ) THEN
            cyy = (Cbj(k)*Cby(k-1)-2.0D0/(pi*Z))/Cbj(k-1)
         ELSE
            cyy = (Cbj(k)*Cby(k-2)-4.0D0*(k-1.0D0)/(pi*Z*Z))/Cbj(k-2)
         ENDIF
         Cby(k) = cyy
      ENDDO
      Cdy(0) = -Cby(1)
      DO k = 1 , Nm
         Cdy(k) = Cby(k-1) - k/Z*Cby(k)
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE IKNB(N,X,Nm,Bi,Di,Bk,Dk)
!
!       ============================================================
!       Purpose: Compute modified Bessel functions In(x) and Kn(x),
!                and their derivatives
!       Input:   x --- Argument of In(x) and Kn(x) ( 0 ≤ x ≤ 700 )
!                n --- Order of In(x) and Kn(x)
!       Output:  BI(n) --- In(x)
!                DI(n) --- In'(x)
!                BK(n) --- Kn(x)
!                DK(n) --- Kn'(x)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 for computing the starting point
!                for backward recurrence
!       ===========================================================
!
      IMPLICIT NONE
!*--IKNB7185
      DOUBLE PRECISION a0 , Bi , Bk , bkl , bs , Di , Dk , el , f , f0 ,&
                     & f1 , g , g0 , g1 , pi , r , s0 , sk0 , vt , X
      INTEGER k , k0 , l , m , MSTA1 , MSTA2 , N , Nm
      DIMENSION Bi(0:N) , Di(0:N) , Bk(0:N) , Dk(0:N)
      pi = 3.141592653589793D0
      el = 0.5772156649015329D0
      Nm = N
      IF ( X<=1.0D-100 ) THEN
         DO k = 0 , N
            Bi(k) = 0.0D0
            Di(k) = 0.0D0
            Bk(k) = 1.0D+300
            Dk(k) = -1.0D+300
         ENDDO
         Bi(0) = 1.0D0
         Di(1) = 0.5D0
         RETURN
      ENDIF
      IF ( N==0 ) Nm = 1
      m = MSTA1(X,200)
      IF ( m<Nm ) THEN
         Nm = m
      ELSE
         m = MSTA2(X,Nm,15)
      ENDIF
      bs = 0.0D0
      sk0 = 0.0D0
      f = 0.0D0
      f0 = 0.0D0
      f1 = 1.0D-100
      DO k = m , 0 , -1
         f = 2.0D0*(k+1.0D0)/X*f1 + f0
         IF ( k<=Nm ) Bi(k) = f
         IF ( k/=0 .AND. k==2*INT(k/2) ) sk0 = sk0 + 4.0D0*f/k
         bs = bs + 2.0D0*f
         f0 = f1
         f1 = f
      ENDDO
      s0 = EXP(X)/(bs-f)
      DO k = 0 , Nm
         Bi(k) = s0*Bi(k)
      ENDDO
      IF ( X<=8.0D0 ) THEN
         Bk(0) = -(DLOG(0.5D0*X)+el)*Bi(0) + s0*sk0
         Bk(1) = (1.0D0/X-Bi(1)*Bk(0))/Bi(0)
      ELSE
         a0 = DSQRT(pi/(2.0D0*X))*EXP(-X)
         k0 = 16
         IF ( X>=25.0 ) k0 = 10
         IF ( X>=80.0 ) k0 = 8
         IF ( X>=200.0 ) k0 = 6
         DO l = 0 , 1
            bkl = 1.0D0
            vt = 4.0D0*l
            r = 1.0D0
            DO k = 1 , k0
               r = 0.125D0*r*(vt-(2.0*k-1.0)**2)/(k*X)
               bkl = bkl + r
            ENDDO
            Bk(l) = a0*bkl
         ENDDO
      ENDIF
      g0 = Bk(0)
      g1 = Bk(1)
      DO k = 2 , Nm
         g = 2.0D0*(k-1.0D0)/X*g1 + g0
         Bk(k) = g
         g0 = g1
         g1 = g
      ENDDO
      Di(0) = Bi(1)
      Dk(0) = -Bk(1)
      DO k = 1 , Nm
         Di(k) = Bi(k-1) - k/X*Bi(k)
         Dk(k) = -Bk(k-1) - k/X*Bk(k)
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE LPMN(Mm,M,N,X,Pm,Pd)
!
!       =====================================================
!       Purpose: Compute the associated Legendre functions
!                Pmn(x) and their derivatives Pmn'(x) for
!                real argument
!       Input :  x  --- Argument of Pmn(x)
!                m  --- Order of Pmn(x),  m = 0,1,2,...,n
!                n  --- Degree of Pmn(x), n = 0,1,2,...,N
!                mm --- Physical dimension of PM and PD
!       Output:  PM(m,n) --- Pmn(x)
!                PD(m,n) --- Pmn'(x)
!       =====================================================
!
      IMPLICIT NONE
!*--LPMN7286
      DOUBLE PRECISION DINF , Pd , Pm , X , xq , xs
      INTEGER i , j , ls , M , Mm , N
      DIMENSION Pm(0:Mm,0:N) , Pd(0:Mm,0:N)
      INTRINSIC MIN
      DO i = 0 , N
         DO j = 0 , M
            Pm(j,i) = 0.0D0
            Pd(j,i) = 0.0D0
         ENDDO
      ENDDO
      Pm(0,0) = 1.0D0
      IF ( N==0 ) RETURN
      IF ( DABS(X)==1.0D0 ) THEN
         DO i = 1 , N
            Pm(0,i) = X**i
            Pd(0,i) = 0.5D0*i*(i+1.0D0)*X**(i+1)
         ENDDO
         DO j = 1 , N
            DO i = 1 , M
               IF ( i==1 ) THEN
                  Pd(i,j) = DINF()
               ELSEIF ( i==2 ) THEN
                  Pd(i,j) = -0.25D0*(j+2)*(j+1)*j*(j-1)*X**(j+1)
               ENDIF
            ENDDO
         ENDDO
         RETURN
      ENDIF
      ls = 1
      IF ( DABS(X)>1.0D0 ) ls = -1
      xq = DSQRT(ls*(1.0D0-X*X))
!       Ensure connection to the complex-valued function for |x| > 1
      IF ( X<-1D0 ) xq = -xq
      xs = ls*(1.0D0-X*X)
      DO i = 1 , M
         Pm(i,i) = -ls*(2.0D0*i-1.0D0)*xq*Pm(i-1,i-1)
      ENDDO
      DO i = 0 , MIN(M,N-1)
         Pm(i,i+1) = (2.0D0*i+1.0D0)*X*Pm(i,i)
      ENDDO
      DO i = 0 , M
         DO j = i + 2 , N
            Pm(i,j) = ((2.0D0*j-1.0D0)*X*Pm(i,j-1)-(i+j-1.0D0)*Pm(i,j-2)&
                    & )/(j-i)
         ENDDO
      ENDDO
      Pd(0,0) = 0.0D0
      DO j = 1 , N
         Pd(0,j) = ls*j*(Pm(0,j-1)-X*Pm(0,j))/xs
      ENDDO
      DO i = 1 , M
         DO j = i , N
            Pd(i,j) = ls*i*X*Pm(i,j)/xs + (j+i)*(j-i+1.0D0)/xq*Pm(i-1,j)
         ENDDO
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE MTU0(Kf,M,Q,X,Csf,Csd)
!
!       ===============================================================
!       Purpose: Compute Mathieu functions cem(x,q) and sem(x,q)
!                and their derivatives ( q ≥ 0 )
!       Input :  KF  --- Function code
!                        KF=1 for computing cem(x,q) and cem'(x,q)
!                        KF=2 for computing sem(x,q) and sem'(x,q)
!                m   --- Order of Mathieu functions
!                q   --- Parameter of Mathieu functions
!                x   --- Argument of Mathieu functions (in degrees)
!       Output:  CSF --- cem(x,q) or sem(x,q)
!                CSD --- cem'x,q) or sem'x,q)
!       Routines called:
!            (1) CVA2 for computing the characteristic values
!            (2) FCOEF for computing the expansion coefficients
!       ===============================================================
!
      IMPLICIT NONE
!*--MTU07368
      DOUBLE PRECISION a , Csd , Csf , DNAN , eps , fg , Q , qm , rd ,  &
                     & X , xr
      INTEGER ic , k , kd , Kf , km , M
      DIMENSION fg(251)
      eps = 1.0D-14
      IF ( Kf==1 .AND. M==2*INT(M/2) ) kd = 1
      IF ( Kf==1 .AND. M/=2*INT(M/2) ) kd = 2
      IF ( Kf==2 .AND. M/=2*INT(M/2) ) kd = 3
      IF ( Kf==2 .AND. M==2*INT(M/2) ) kd = 4
      CALL CVA2(kd,M,Q,a)
      IF ( Q<=1.0D0 ) THEN
         qm = 7.5 + 56.1*SQRT(Q) - 134.7*Q + 90.7*SQRT(Q)*Q
      ELSE
         qm = 17.0 + 3.1*SQRT(Q) - .126*Q + .0037*SQRT(Q)*Q
      ENDIF
      km = INT(qm+0.5*M)
      IF ( km>251 ) THEN
         Csf = DNAN()
         Csd = DNAN()
         RETURN
      ENDIF
      CALL FCOEF(kd,M,Q,a,fg)
      ic = INT(M/2) + 1
      rd = 1.74532925199433D-2
      xr = X*rd
      Csf = 0.0D0
      DO k = 1 , km
         IF ( kd==1 ) THEN
            Csf = Csf + fg(k)*DCOS((2*k-2)*xr)
         ELSEIF ( kd==2 ) THEN
            Csf = Csf + fg(k)*DCOS((2*k-1)*xr)
         ELSEIF ( kd==3 ) THEN
            Csf = Csf + fg(k)*DSIN((2*k-1)*xr)
         ELSEIF ( kd==4 ) THEN
            Csf = Csf + fg(k)*DSIN(2*k*xr)
         ENDIF
         IF ( k>=ic .AND. DABS(fg(k))<DABS(Csf)*eps ) GOTO 100
      ENDDO
 100  Csd = 0.0D0
      DO k = 1 , km
         IF ( kd==1 ) THEN
            Csd = Csd - (2*k-2)*fg(k)*DSIN((2*k-2)*xr)
         ELSEIF ( kd==2 ) THEN
            Csd = Csd - (2*k-1)*fg(k)*DSIN((2*k-1)*xr)
         ELSEIF ( kd==3 ) THEN
            Csd = Csd + (2*k-1)*fg(k)*DCOS((2*k-1)*xr)
         ELSEIF ( kd==4 ) THEN
            Csd = Csd + 2.0D0*k*fg(k)*DCOS(2*k*xr)
         ENDIF
         IF ( k>=ic .AND. DABS(fg(k))<DABS(Csd)*eps ) GOTO 99999
      ENDDO
99999 END
 
 
 
!       **********************************
 
      SUBROUTINE CY01(Kf,Z,Zf,Zd)
!
!       ===========================================================
!       Purpose: Compute complex Bessel functions Y0(z), Y1(z)
!                and their derivatives
!       Input :  z  --- Complex argument of Yn(z) ( n=0,1 )
!                KF --- Function choice code
!                    KF=0 for ZF=Y0(z) and ZD=Y0'(z)
!                    KF=1 for ZF=Y1(z) and ZD=Y1'(z)
!                    KF=2 for ZF=Y1'(z) and ZD=Y1''(z)
!       Output:  ZF --- Y0(z) or Y1(z) or Y1'(z)
!                ZD --- Y0'(z) or Y1'(z) or Y1''(z)
!       ===========================================================
!
      IMPLICIT NONE
!*--CY017444
      DOUBLE PRECISION a , a0 , a1 , b , b1 , el , pi , rp2 , w0 , w1
      COMPLEX*16 cbj0 , cbj1 , cby0 , cby1 , cdy0 , cdy1 , ci , cp ,    &
               & cp0 , cp1 , cq0 , cq1 , cr , cs , ct1 , ct2 , cu , Z , &
               & z1 , z2
      COMPLEX*16 Zd , Zf
      INTEGER k , k0 , Kf
      DIMENSION a(12) , b(12) , a1(12) , b1(12)
      pi = 3.141592653589793D0
      el = 0.5772156649015329D0
      rp2 = 2.0D0/pi
      ci = (0.0D0,1.0D0)
      a0 = ABS(Z)
      z2 = Z*Z
      z1 = Z
      IF ( a0==0.0D0 ) THEN
         cbj0 = (1.0D0,0.0D0)
         cbj1 = (0.0D0,0.0D0)
         cby0 = -(1.0D300,0.0D0)
         cby1 = -(1.0D300,0.0D0)
         cdy0 = (1.0D300,0.0D0)
         cdy1 = (1.0D300,0.0D0)
         GOTO 300
      ENDIF
      IF ( DBLE(Z)<0.0 ) z1 = -Z
      IF ( a0<=12.0 ) THEN
         cbj0 = (1.0D0,0.0D0)
         cr = (1.0D0,0.0D0)
         DO k = 1 , 40
            cr = -0.25D0*cr*z2/(k*k)
            cbj0 = cbj0 + cr
            IF ( ABS(cr)<ABS(cbj0)*1.0D-15 ) GOTO 50
         ENDDO
 50      cbj1 = (1.0D0,0.0D0)
         cr = (1.0D0,0.0D0)
         DO k = 1 , 40
            cr = -0.25D0*cr*z2/(k*(k+1.0D0))
            cbj1 = cbj1 + cr
            IF ( ABS(cr)<ABS(cbj1)*1.0D-15 ) GOTO 100
         ENDDO
 100     cbj1 = 0.5D0*z1*cbj1
         w0 = 0.0D0
         cr = (1.0D0,0.0D0)
         cs = (0.0D0,0.0D0)
         DO k = 1 , 40
            w0 = w0 + 1.0D0/k
            cr = -0.25D0*cr/(k*k)*z2
            cp = cr*w0
            cs = cs + cp
            IF ( ABS(cp)<ABS(cs)*1.0D-15 ) GOTO 150
         ENDDO
 150     cby0 = rp2*(LOG(z1/2.0D0)+el)*cbj0 - rp2*cs
         w1 = 0.0D0
         cr = (1.0D0,0.0D0)
         cs = (1.0D0,0.0D0)
         DO k = 1 , 40
            w1 = w1 + 1.0D0/k
            cr = -0.25D0*cr/(k*(k+1))*z2
            cp = cr*(2.0D0*w1+1.0D0/(k+1.0D0))
            cs = cs + cp
            IF ( ABS(cp)<ABS(cs)*1.0D-15 ) GOTO 200
         ENDDO
 200     cby1 = rp2*((LOG(z1/2.0D0)+el)*cbj1-1.0D0/z1-.25D0*z1*cs)
      ELSE
         DATA a/ - .703125D-01 , .112152099609375D+00 ,                 &
            & -.5725014209747314D+00 , .6074042001273483D+01 ,          &
            & -.1100171402692467D+03 , .3038090510922384D+04 ,          &
            & -.1188384262567832D+06 , .6252951493434797D+07 ,          &
            & -.4259392165047669D+09 , .3646840080706556D+11 ,          &
            & -.3833534661393944D+13 , .4854014686852901D+15/
         DATA b/.732421875D-01 , -.2271080017089844D+00 ,               &
            & .1727727502584457D+01 , -.2438052969955606D+02 ,          &
            & .5513358961220206D+03 , -.1825775547429318D+05 ,          &
            & .8328593040162893D+06 , -.5006958953198893D+08 ,          &
            & .3836255180230433D+10 , -.3649010818849833D+12 ,          &
            & .4218971570284096D+14 , -.5827244631566907D+16/
         DATA a1/.1171875D+00 , -.144195556640625D+00 ,                 &
            & .6765925884246826D+00 , -.6883914268109947D+01 ,          &
            & .1215978918765359D+03 , -.3302272294480852D+04 ,          &
            & .1276412726461746D+06 , -.6656367718817688D+07 ,          &
            & .4502786003050393D+09 , -.3833857520742790D+11 ,          &
            & .4011838599133198D+13 , -.5060568503314727D+15/
         DATA b1/ - .1025390625D+00 , .2775764465332031D+00 ,           &
            & -.1993531733751297D+01 , .2724882731126854D+02 ,          &
            & -.6038440767050702D+03 , .1971837591223663D+05 ,          &
            & -.8902978767070678D+06 , .5310411010968522D+08 ,          &
            & -.4043620325107754D+10 , .3827011346598605D+12 ,          &
            & -.4406481417852278D+14 , .6065091351222699D+16/
         k0 = 12
         IF ( a0>=35.0 ) k0 = 10
         IF ( a0>=50.0 ) k0 = 8
         ct1 = z1 - .25D0*pi
         cp0 = (1.0D0,0.0D0)
         DO k = 1 , k0
            cp0 = cp0 + a(k)*z1**(-2*k)
         ENDDO
         cq0 = -0.125D0/z1
         DO k = 1 , k0
            cq0 = cq0 + b(k)*z1**(-2*k-1)
         ENDDO
         cu = SQRT(rp2/z1)
         cbj0 = cu*(cp0*COS(ct1)-cq0*SIN(ct1))
         cby0 = cu*(cp0*SIN(ct1)+cq0*COS(ct1))
         ct2 = z1 - .75D0*pi
         cp1 = (1.0D0,0.0D0)
         DO k = 1 , k0
            cp1 = cp1 + a1(k)*z1**(-2*k)
         ENDDO
         cq1 = 0.375D0/z1
         DO k = 1 , k0
            cq1 = cq1 + b1(k)*z1**(-2*k-1)
         ENDDO
         cbj1 = cu*(cp1*COS(ct2)-cq1*SIN(ct2))
         cby1 = cu*(cp1*SIN(ct2)+cq1*COS(ct2))
      ENDIF
      IF ( DBLE(Z)<0.0 ) THEN
         IF ( DIMAG(Z)<0.0 ) cby0 = cby0 - 2.0D0*ci*cbj0
         IF ( DIMAG(Z)>0.0 ) cby0 = cby0 + 2.0D0*ci*cbj0
         IF ( DIMAG(Z)<0.0 ) cby1 = -(cby1-2.0D0*ci*cbj1)
         IF ( DIMAG(Z)>0.0 ) cby1 = -(cby1+2.0D0*ci*cbj1)
         cbj1 = -cbj1
      ENDIF
      cdy0 = -cby1
      cdy1 = cby0 - 1.0D0/Z*cby1
 300  IF ( Kf==0 ) THEN
         Zf = cby0
         Zd = cdy0
      ELSEIF ( Kf==1 ) THEN
         Zf = cby1
         Zd = cdy1
      ELSEIF ( Kf==2 ) THEN
         Zf = cdy1
         Zd = -cdy1/Z - (1.0D0-1.0D0/(Z*Z))*cby1
      ENDIF
      END
 
 
!       **********************************
 
      SUBROUTINE FFK(Ks,X,Fr,Fi,Fm,Fa,Gr,Gi,Gm,Ga)
!
!       =======================================================
!       Purpose: Compute modified Fresnel integrals F±(x)
!                and K±(x)
!       Input :  x   --- Argument of F±(x) and K±(x)
!                KS  --- Sign code
!                        KS=0 for calculating F+(x) and K+(x)
!                        KS=1 for calculating F_(x) and K_(x)
!       Output:  FR  --- Re[F±(x)]
!                FI  --- Im[F±(x)]
!                FM  --- |F±(x)|
!                FA  --- Arg[F±(x)]  (Degs.)
!                GR  --- Re[K±(x)]
!                GI  --- Im[K±(x)]
!                GM  --- |K±(x)|
!                GA  --- Arg[K±(x)]  (Degs.)
!       ======================================================
!
      IMPLICIT NONE
!*--FFK7606
      DOUBLE PRECISION c1 , cs , eps , Fa , Fi , fi0 , Fm , Fr , Ga ,   &
                     & Gi , Gm , Gr , p2p , pi , pp2 , s1 , srd , ss ,  &
                     & X , x2
      DOUBLE PRECISION x4 , xa , xc , xf , xf0 , xf1 , xg , xp , xq ,   &
                     & xq2 , xr , xs , xsu , xw
      INTEGER k , Ks , m
      srd = 57.29577951308233D0
      eps = 1.0D-15
      pi = 3.141592653589793D0
      pp2 = 1.2533141373155D0
      p2p = .7978845608028654D0
      xa = DABS(X)
      x2 = X*X
      x4 = x2*x2
      IF ( X==0.0D0 ) THEN
         Fr = .5D0*DSQRT(0.5D0*pi)
         Fi = (-1)**Ks*Fr
         Fm = DSQRT(0.25D0*pi)
         Fa = (-1)**Ks*45.0D0
         Gr = .5D0
         Gi = 0.0D0
         Gm = .5D0
         Ga = 0.0D0
      ELSE
         IF ( xa<=2.5D0 ) THEN
            xr = p2p*xa
            c1 = xr
            DO k = 1 , 50
               xr = -.5D0*xr*(4.0D0*k-3.0D0)/k/(2.0D0*k-1.0D0)          &
                  & /(4.0D0*k+1.0D0)*x4
               c1 = c1 + xr
               IF ( DABS(xr/c1)<eps ) GOTO 20
            ENDDO
 20         s1 = p2p*xa*xa*xa/3.0D0
            xr = s1
            DO k = 1 , 50
               xr = -.5D0*xr*(4.0D0*k-1.0D0)/k/(2.0D0*k+1.0D0)          &
                  & /(4.0D0*k+3.0D0)*x4
               s1 = s1 + xr
               IF ( DABS(xr/s1)<eps ) GOTO 50
            ENDDO
         ELSEIF ( xa<5.5D0 ) THEN
            m = INT(42+1.75*x2)
            xsu = 0.0D0
            xc = 0.0D0
            xs = 0.0D0
            xf1 = 0.0D0
            xf0 = 1D-100
            DO k = m , 0 , -1
               xf = (2.0D0*k+3.0D0)*xf0/x2 - xf1
               IF ( k==2*INT(k/2) ) THEN
                  xc = xc + xf
               ELSE
                  xs = xs + xf
               ENDIF
               xsu = xsu + (2.0D0*k+1.0D0)*xf*xf
               xf1 = xf0
               xf0 = xf
            ENDDO
            xq = DSQRT(xsu)
            xw = p2p*xa/xq
            c1 = xc*xw
            s1 = xs*xw
         ELSE
            xr = 1.0D0
            xf = 1.0D0
            DO k = 1 , 12
               xr = -.25D0*xr*(4.0D0*k-1.0D0)*(4.0D0*k-3.0D0)/x4
               xf = xf + xr
            ENDDO
            xr = 1.0D0/(2.0D0*xa*xa)
            xg = xr
            DO k = 1 , 12
               xr = -.25D0*xr*(4.0D0*k+1.0D0)*(4.0D0*k-1.0D0)/x4
               xg = xg + xr
            ENDDO
            c1 = .5D0 + (xf*DSIN(x2)-xg*DCOS(x2))/DSQRT(2.0D0*pi)/xa
            s1 = .5D0 - (xf*DCOS(x2)+xg*DSIN(x2))/DSQRT(2.0D0*pi)/xa
         ENDIF
 50      Fr = pp2*(.5D0-c1)
         fi0 = pp2*(.5D0-s1)
         Fi = (-1)**Ks*fi0
         Fm = DSQRT(Fr*Fr+Fi*Fi)
         IF ( Fr>=0.0 ) THEN
            Fa = srd*DATAN(Fi/Fr)
         ELSEIF ( Fi>0.0 ) THEN
            Fa = srd*(DATAN(Fi/Fr)+pi)
         ELSEIF ( Fi<0.0 ) THEN
            Fa = srd*(DATAN(Fi/Fr)-pi)
         ENDIF
         xp = X*X + pi/4.0D0
         cs = DCOS(xp)
         ss = DSIN(xp)
         xq2 = 1.0D0/DSQRT(pi)
         Gr = xq2*(Fr*cs+fi0*ss)
         Gi = (-1)**Ks*xq2*(fi0*cs-Fr*ss)
         Gm = DSQRT(Gr*Gr+Gi*Gi)
         IF ( Gr>=0.0 ) THEN
            Ga = srd*DATAN(Gi/Gr)
         ELSEIF ( Gi>0.0 ) THEN
            Ga = srd*(DATAN(Gi/Gr)+pi)
         ELSEIF ( Gi<0.0 ) THEN
            Ga = srd*(DATAN(Gi/Gr)-pi)
         ENDIF
         IF ( X<0.0D0 ) THEN
            Fr = pp2 - Fr
            Fi = (-1)**Ks*pp2 - Fi
            Fm = DSQRT(Fr*Fr+Fi*Fi)
            Fa = srd*DATAN(Fi/Fr)
            Gr = DCOS(X*X) - Gr
            Gi = -(-1)**Ks*DSIN(X*X) - Gi
            Gm = DSQRT(Gr*Gr+Gi*Gi)
            Ga = srd*DATAN(Gi/Gr)
         ENDIF
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE AIRYA(X,Ai,Bi,Ad,Bd)
!
!       ======================================================
!       Purpose: Compute Airy functions and their derivatives
!       Input:   x  --- Argument of Airy function
!       Output:  AI --- Ai(x)
!                BI --- Bi(x)
!                AD --- Ai'(x)
!                BD --- Bi'(x)
!       Routine called:
!                AJYIK for computing Jv(x), Yv(x), Iv(x) and
!                Kv(x) with v=1/3 and 2/3
!       ======================================================
!
      IMPLICIT NONE
!*--AIRYA7744
      DOUBLE PRECISION Ad , Ai , Bd , Bi , c1 , c2 , pir , sr3 , vi1 ,  &
                     & vi2 , vj1 , vj2 , vk1 , vk2 , vy1 , vy2 , X ,    &
                     & xa , xq , z
      xa = DABS(X)
      pir = 0.318309886183891D0
      c1 = 0.355028053887817D0
      c2 = 0.258819403792807D0
      sr3 = 1.732050807568877D0
      z = xa**1.5/1.5D0
      xq = DSQRT(xa)
      CALL AJYIK(z,vj1,vj2,vy1,vy2,vi1,vi2,vk1,vk2)
      IF ( X==0.0D0 ) THEN
         Ai = c1
         Bi = sr3*c1
         Ad = -c2
         Bd = sr3*c2
      ELSEIF ( X>0.0D0 ) THEN
         Ai = pir*xq/sr3*vk1
         Bi = xq*(pir*vk1+2.0D0/sr3*vi1)
         Ad = -xa/sr3*pir*vk2
         Bd = xa*(pir*vk2+2.0D0/sr3*vi2)
      ELSE
         Ai = 0.5D0*xq*(vj1-vy1/sr3)
         Bi = -0.5D0*xq*(vj1/sr3+vy1)
         Ad = 0.5D0*xa*(vj2+vy2/sr3)
         Bd = 0.5D0*xa*(vj2/sr3-vy2)
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE AIRYB(X,Ai,Bi,Ad,Bd)
!
!       =======================================================
!       Purpose: Compute Airy functions and their derivatives
!       Input:   x  --- Argument of Airy function
!       Output:  AI --- Ai(x)
!                BI --- Bi(x)
!                AD --- Ai'(x)
!                BD --- Bi'(x)
!       =======================================================
!
      IMPLICIT NONE
!*--AIRYB7793
      DOUBLE PRECISION Ad , Ai , Bd , Bi , c1 , c2 , ck , df , dg , dk ,&
                     & eps , fx , gx , pi , r , rp , sad , sai , sbd ,  &
                     & sbi
      DOUBLE PRECISION sda , sdb , sr3 , ssa , ssb , X , xa , xar ,     &
                     & xcs , xe , xf , xm , xp1 , xq , xr1 , xr2 , xss
      INTEGER k , km , km2 , kmax
      DIMENSION ck(51) , dk(51)
      eps = 1.0D-15
      pi = 3.141592653589793D0
      c1 = 0.355028053887817D0
      c2 = 0.258819403792807D0
      sr3 = 1.732050807568877D0
      xa = DABS(X)
      xq = DSQRT(xa)
      xm = 8.0D0
      IF ( X>0.0D0 ) xm = 5.0D0
      IF ( X==0.0D0 ) THEN
         Ai = c1
         Bi = sr3*c1
         Ad = -c2
         Bd = sr3*c2
         RETURN
      ENDIF
      IF ( xa<=xm ) THEN
         fx = 1.0D0
         r = 1.0D0
         DO k = 1 , 40
            r = r*X/(3.0D0*k)*X/(3.0D0*k-1.0D0)*X
            fx = fx + r
            IF ( DABS(r)<DABS(fx)*eps ) GOTO 50
         ENDDO
 50      gx = X
         r = X
         DO k = 1 , 40
            r = r*X/(3.0D0*k)*X/(3.0D0*k+1.0D0)*X
            gx = gx + r
            IF ( DABS(r)<DABS(gx)*eps ) GOTO 100
         ENDDO
 100     Ai = c1*fx - c2*gx
         Bi = sr3*(c1*fx+c2*gx)
         df = 0.5D0*X*X
         r = df
         DO k = 1 , 40
            r = r*X/(3.0D0*k)*X/(3.0D0*k+2.0D0)*X
            df = df + r
            IF ( DABS(r)<DABS(df)*eps ) GOTO 150
         ENDDO
 150     dg = 1.0D0
         r = 1.0D0
         DO k = 1 , 40
            r = r*X/(3.0D0*k)*X/(3.0D0*k-2.0D0)*X
            dg = dg + r
            IF ( DABS(r)<DABS(dg)*eps ) GOTO 200
         ENDDO
 200     Ad = c1*df - c2*dg
         Bd = sr3*(c1*df+c2*dg)
      ELSE
         km = INT(24.5-xa)
         IF ( xa<6.0 ) km = 14
         IF ( xa>15.0 ) km = 10
         IF ( X>0.0D0 ) THEN
            kmax = km
         ELSE
!             Choose cutoffs so that the remainder term in asymptotic
!             expansion is epsilon size. The X<0 branch needs to be fast
!             in order to make AIRYZO efficient
            IF ( xa>70.0 ) km = 3
            IF ( xa>500.0 ) km = 2
            IF ( xa>1000.0 ) km = 1
            km2 = km
            IF ( xa>150.0 ) km2 = 1
            IF ( xa>3000.0 ) km2 = 0
            kmax = 2*km + 1
         ENDIF
         xe = xa*xq/1.5D0
         xr1 = 1.0D0/xe
         xar = 1.0D0/xq
         xf = DSQRT(xar)
         rp = 0.5641895835477563D0
         r = 1.0D0
         DO k = 1 , kmax
            r = r*(6.0D0*k-1.0D0)/216.0D0*(6.0D0*k-3.0D0)               &
              & /k*(6.0D0*k-5.0D0)/(2.0D0*k-1.0D0)
            ck(k) = r
            dk(k) = -(6.0D0*k+1.0D0)/(6.0D0*k-1.0D0)*ck(k)
         ENDDO
         IF ( X>0.0D0 ) THEN
            sai = 1.0D0
            sad = 1.0D0
            r = 1.0D0
            DO k = 1 , km
               r = -r*xr1
               sai = sai + ck(k)*r
               sad = sad + dk(k)*r
            ENDDO
            sbi = 1.0D0
            sbd = 1.0D0
            r = 1.0D0
            DO k = 1 , km
               r = r*xr1
               sbi = sbi + ck(k)*r
               sbd = sbd + dk(k)*r
            ENDDO
            xp1 = EXP(-xe)
            Ai = 0.5D0*rp*xf*xp1*sai
            Bi = rp*xf/xp1*sbi
            Ad = -.5D0*rp/xf*xp1*sad
            Bd = rp/xf/xp1*sbd
         ELSE
            xcs = DCOS(xe+pi/4.0D0)
            xss = DSIN(xe+pi/4.0D0)
            ssa = 1.0D0
            sda = 1.0D0
            r = 1.0D0
            xr2 = 1.0D0/(xe*xe)
            DO k = 1 , km
               r = -r*xr2
               ssa = ssa + ck(2*k)*r
               sda = sda + dk(2*k)*r
            ENDDO
            ssb = ck(1)*xr1
            sdb = dk(1)*xr1
            r = xr1
            DO k = 1 , km2
               r = -r*xr2
               ssb = ssb + ck(2*k+1)*r
               sdb = sdb + dk(2*k+1)*r
            ENDDO
            Ai = rp*xf*(xss*ssa-xcs*ssb)
            Bi = rp*xf*(xcs*ssa+xss*ssb)
            Ad = -rp/xf*(xcs*sda+xss*sdb)
            Bd = rp/xf*(xss*sda-xcs*sdb)
         ENDIF
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE SCKA(M,N,C,Cv,Kd,Ck)
!
!       ======================================================
!       Purpose: Compute the expansion coefficients of the
!                prolate and oblate spheroidal functions, c2k
!       Input :  m  --- Mode parameter
!                n  --- Mode parameter
!                c  --- Spheroidal parameter
!                cv --- Characteristic value
!                KD --- Function code
!                       KD=1 for prolate; KD=-1 for oblate
!       Output:  CK(k) --- Expansion coefficients ck;
!                          CK(1), CK(2),... correspond to
!                          c0, c2,...
!       ======================================================
!
      IMPLICIT NONE
!*--SCKA7952
      DOUBLE PRECISION C , Ck , cs , Cv , f , f0 , f1 , f2 , fl , fs ,  &
                     & r1 , r2 , s0 , su1 , su2
      INTEGER ip , j , k , k1 , kb , Kd , M , N , nm
      DIMENSION Ck(200)
      IF ( C<=1.0D-10 ) C = 1.0D-10
      nm = 25 + INT((N-M)/2+C)
      cs = C*C*Kd
      ip = 1
      IF ( N-M==2*INT((N-M)/2) ) ip = 0
      fs = 1.0D0
      f1 = 0.0D0
      f0 = 1.0D-100
      kb = 0
      Ck(nm+1) = 0.0D0
      fl = 0.0D0
      DO k = nm , 1 , -1
         f = (((2.0D0*k+M+ip)*(2.0D0*k+M+1.0D0+ip)-Cv+cs)               &
           & *f0-4.0D0*(k+1.0D0)*(k+M+1.0D0)*f1)/cs
         IF ( DABS(f)>DABS(Ck(k+1)) ) THEN
            Ck(k) = f
            f1 = f0
            f0 = f
            IF ( DABS(f)>1.0D+100 ) THEN
               DO k1 = nm , k , -1
                  Ck(k1) = Ck(k1)*1.0D-100
               ENDDO
               f1 = f1*1.0D-100
               f0 = f0*1.0D-100
            ENDIF
         ELSE
            kb = k
            fl = Ck(k+1)
            f1 = 1.0D0
            f2 = 0.25D0*((M+ip)*(M+ip+1.0)-Cv+cs)/(M+1.0)*f1
            Ck(1) = f1
            IF ( kb==1 ) THEN
               fs = f2
            ELSEIF ( kb==2 ) THEN
               Ck(2) = f2
               fs = 0.125D0*(((M+ip+2.0)*(M+ip+3.0)-Cv+cs)*f2-cs*f1)    &
                  & /(M+2.0)
            ELSE
               Ck(2) = f2
               DO j = 3 , kb + 1
                  f = 0.25D0*(((2.0*j+M+ip-4.0)*(2.0*j+M+ip-3.0)-Cv+cs) &
                    & *f2-cs*f1)/((j-1.0)*(j+M-1.0))
                  IF ( j<=kb ) Ck(j) = f
                  f1 = f2
                  f2 = f
               ENDDO
               fs = f
            ENDIF
            GOTO 100
         ENDIF
      ENDDO
 100  su1 = 0.0D0
      DO k = 1 , kb
         su1 = su1 + Ck(k)
      ENDDO
      su2 = 0.0D0
      DO k = kb + 1 , nm
         su2 = su2 + Ck(k)
      ENDDO
      r1 = 1.0D0
      DO j = 1 , (N+M+ip)/2
         r1 = r1*(j+0.5D0*(N+M+ip))
      ENDDO
      r2 = 1.0D0
      DO j = 1 , (N-M-ip)/2
         r2 = -r2*j
      ENDDO
      IF ( kb==0 ) THEN
         s0 = r1/(2.0D0**N*r2*su2)
      ELSE
         s0 = r1/(2.0D0**N*r2*(fl/fs*su1+su2))
      ENDIF
      DO k = 1 , kb
         Ck(k) = fl/fs*s0*Ck(k)
      ENDDO
      DO k = kb + 1 , nm
         Ck(k) = s0*Ck(k)
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE SCKB(M,N,C,Df,Ck)
!
!       ======================================================
!       Purpose: Compute the expansion coefficients of the
!                prolate and oblate spheroidal functions
!       Input :  m  --- Mode parameter
!                n  --- Mode parameter
!                c  --- Spheroidal parameter
!                DF(k) --- Expansion coefficients dk
!       Output:  CK(k) --- Expansion coefficients ck;
!                          CK(1), CK(2), ... correspond to
!                          c0, c2, ...
!       ======================================================
!
      IMPLICIT NONE
!*--SCKB8059
      DOUBLE PRECISION C , Ck , d1 , d2 , d3 , Df , fac , r , r1 , reg ,&
                     & sum , sw
      INTEGER i , i1 , i2 , ip , k , M , N , nm
      DIMENSION Df(200) , Ck(200)
      IF ( C<=1.0D-10 ) C = 1.0D-10
      nm = 25 + INT(0.5*(N-M)+C)
      ip = 1
      IF ( N-M==2*INT((N-M)/2) ) ip = 0
      reg = 1.0D0
      IF ( M+nm>80 ) reg = 1.0D-200
      fac = -0.5D0**M
      sw = 0.0D0
      DO k = 0 , nm - 1
         fac = -fac
         i1 = 2*k + ip + 1
         r = reg
         DO i = i1 , i1 + 2*M - 1
            r = r*i
         ENDDO
         i2 = k + M + ip
         DO i = i2 , i2 + k - 1
            r = r*(i+0.5D0)
         ENDDO
         sum = r*Df(k+1)
         DO i = k + 1 , nm
            d1 = 2.0D0*i + ip
            d2 = 2.0D0*M + d1
            d3 = i + M + ip - 0.5D0
            r = r*d2*(d2-1.0D0)*i*(d3+k)/(d1*(d1-1.0D0)*(i-k)*d3)
            sum = sum + r*Df(i+1)
            IF ( DABS(sw-sum)<DABS(sum)*1.0D-14 ) GOTO 50
            sw = sum
         ENDDO
 50      r1 = reg
         DO i = 2 , M + k
            r1 = r1*i
         ENDDO
         Ck(k+1) = fac*sum/r1
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE CPDLA(N,Z,Cdn)
!
!       ===========================================================
!       Purpose: Compute complex parabolic cylinder function Dn(z)
!                for large argument
!       Input:   z   --- Complex argument of Dn(z)
!                n   --- Order of Dn(z) (n = 0,±1,±2,…)
!       Output:  CDN --- Dn(z)
!       ===========================================================
!
      IMPLICIT NONE
!*--CPDLA8119
      COMPLEX*16 cb0 , Cdn , cr , Z
      INTEGER k , N
      cb0 = Z**N*EXP(-.25D0*Z*Z)
      cr = (1.0D0,0.0D0)
      Cdn = (1.0D0,0.0D0)
      DO k = 1 , 16
         cr = -0.5D0*cr*(2.0*k-N-1.0)*(2.0*k-N-2.0)/(k*Z*Z)
         Cdn = Cdn + cr
         IF ( ABS(cr)<ABS(Cdn)*1.0D-12 ) GOTO 100
      ENDDO
 100  Cdn = cb0*Cdn
      END
 
 
 
!       **********************************
 
      SUBROUTINE FCSZO(Kf,Nt,Zo)
!
!       ===============================================================
!       Purpose: Compute the complex zeros of Fresnel integral C(z)
!                or S(z) using modified Newton's iteration method
!       Input :  KF  --- Function code
!                        KF=1 for C(z) or KF=2 for S(z)
!                NT  --- Total number of zeros
!       Output:  ZO(L) --- L-th zero of C(z) or S(z)
!       Routines called:
!            (1) CFC for computing Fresnel integral C(z)
!            (2) CFS for computing Fresnel integral S(z)
!       ==============================================================
!
      IMPLICIT NONE
!*--FCSZO8155
      INTEGER i , it , j , Kf , nr , Nt
      DOUBLE PRECISION pi , psq , px , py , w , w0
      COMPLEX*16 z , zd , zf , zfd , zgd , Zo , zp , zq , zw
      DIMENSION Zo(Nt)
      pi = 3.141592653589793D0
      psq = 0.0D0
      w = 0.0D0
      DO nr = 1 , Nt
         IF ( Kf==1 ) psq = DSQRT(4.0D0*nr-1.0D0)
         IF ( Kf==2 ) psq = 2.0D0*nr**(0.5)
         px = psq - DLOG(pi*psq)/(pi*pi*psq**3.0)
         py = DLOG(pi*psq)/(pi*psq)
         z = DCMPLX(px,py)
         IF ( Kf==2 ) THEN
            IF ( nr==2 ) z = (2.8334,0.2443)
            IF ( nr==3 ) z = (3.4674,0.2185)
            IF ( nr==4 ) z = (4.0025,0.2008)
         ENDIF
         it = 0
 50      it = it + 1
         IF ( Kf==1 ) CALL CFC(z,zf,zd)
         IF ( Kf==2 ) CALL CFS(z,zf,zd)
         zp = (1.0D0,0.0D0)
         DO i = 1 , nr - 1
            zp = zp*(z-Zo(i))
         ENDDO
         zfd = zf/zp
         zq = (0.0D0,0.0D0)
         DO i = 1 , nr - 1
            zw = (1.0D0,0.0D0)
            DO j = 1 , nr - 1
               IF ( j/=i ) zw = zw*(z-Zo(j))
            ENDDO
            zq = zq + zw
         ENDDO
         zgd = (zd-zq*zfd)/zp
         z = z - zfd/zgd
         w0 = w
         w = ABS(z)
         IF ( it<=50 .AND. DABS((w-w0)/w)>1.0D-12 ) GOTO 50
         Zo(nr) = z
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE E1XA(X,E1)
!
!       ============================================
!       Purpose: Compute exponential integral E1(x)
!       Input :  x  --- Argument of E1(x)
!       Output:  E1 --- E1(x) ( x > 0 )
!       ============================================
!
      IMPLICIT NONE
!*--E1XA8216
      DOUBLE PRECISION E1 , es1 , es2 , X
      IF ( X==0.0 ) THEN
         E1 = 1.0D+300
      ELSEIF ( X<=1.0 ) THEN
         E1 = -DLOG(X) + ((((1.07857D-3*X-9.76004D-3)*X+5.519968D-2)*X- &
            & 0.24991055D0)*X+0.99999193D0)*X - 0.57721566D0
      ELSE
         es1 = (((X+8.5733287401D0)*X+18.059016973D0)*X+8.6347608925D0) &
             & *X + 0.2677737343D0
         es2 = (((X+9.5733223454D0)*X+25.6329561486D0)                  &
             & *X+21.0996530827D0)*X + 3.9584969228D0
         E1 = EXP(-X)/X*es1/es2
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE LPMV0(V,M,X,Pmv)
!
!       =======================================================
!       Purpose: Compute the associated Legendre function
!                Pmv(x) with an integer order and an arbitrary
!                nonnegative degree v
!       Input :  x   --- Argument of Pm(x)  ( -1 ≤ x ≤ 1 )
!                m   --- Order of Pmv(x)
!                v   --- Degree of Pmv(x)
!       Output:  PMV --- Pmv(x)
!       Routine called:  PSI_SPEC for computing Psi function
!       =======================================================
!
      IMPLICIT NONE
!*--LPMV08251
      DOUBLE PRECISION c0 , el , eps , pa , pi , Pmv , pss , psv , pv0 ,&
                     & qr , r , r0 , r1 , r2 , rg , s , s0 , s1 , s2 , V
      DOUBLE PRECISION v0 , vs , X , xq
      INTEGER j , k , M , nv
      pi = 3.141592653589793D0
      el = .5772156649015329D0
      eps = 1.0D-14
      nv = INT(V)
      v0 = V - nv
      IF ( X==-1.0D0 .AND. V/=nv ) THEN
         IF ( M==0 ) Pmv = -1.0D+300
         IF ( M/=0 ) Pmv = 1.0D+300
         RETURN
      ENDIF
      c0 = 1.0D0
      IF ( M/=0 ) THEN
         rg = V*(V+M)
         DO j = 1 , M - 1
            rg = rg*(V*V-j*j)
         ENDDO
         xq = DSQRT(1.0D0-X*X)
         r0 = 1.0D0
         DO j = 1 , M
            r0 = .5D0*r0*xq/j
         ENDDO
         c0 = r0*rg
      ENDIF
      IF ( v0==0.0D0 ) THEN
!          DLMF 14.3.4, 14.7.17, 15.2.4
         Pmv = 1.0D0
         r = 1.0D0
         DO k = 1 , nv - M
            r = 0.5D0*r*(-nv+M+k-1.0D0)*(nv+M+k)/(k*(k+M))*(1.0D0+X)
            Pmv = Pmv + r
         ENDDO
         Pmv = (-1)**nv*c0*Pmv
      ELSEIF ( X>=-0.35D0 ) THEN
!             DLMF 14.3.4, 15.2.1
         Pmv = 1.0D0
         r = 1.0D0
         DO k = 1 , 100
            r = 0.5D0*r*(-V+M+k-1.0D0)*(V+M+k)/(k*(M+k))*(1.0D0-X)
            Pmv = Pmv + r
            IF ( k>12 .AND. DABS(r/Pmv)<eps ) GOTO 50
         ENDDO
 50      Pmv = (-1)**M*c0*Pmv
      ELSE
!             DLMF 14.3.5, 15.8.10
         vs = DSIN(V*pi)/pi
         pv0 = 0.0D0
         IF ( M/=0 ) THEN
            qr = DSQRT((1.0D0-X)/(1.0D0+X))
            r2 = 1.0D0
            DO j = 1 , M
               r2 = r2*qr*j
            ENDDO
            s0 = 1.0D0
            r1 = 1.0D0
            DO k = 1 , M - 1
               r1 = 0.5D0*r1*(-V+k-1)*(V+k)/(k*(k-M))*(1.0D0+X)
               s0 = s0 + r1
            ENDDO
            pv0 = -vs*r2/M*s0
         ENDIF
         CALL PSI_SPEC(V,psv)
         pa = 2.0D0*(psv+el) + pi/DTAN(pi*V) + 1.0D0/V
         s1 = 0.0D0
         DO j = 1 , M
            s1 = s1 + (j*j+V*V)/(j*(j*j-V*V))
         ENDDO
         Pmv = pa + s1 - 1.0D0/(M-V) + DLOG(0.5D0*(1.0D0+X))
         r = 1.0D0
         DO k = 1 , 100
            r = 0.5D0*r*(-V+M+k-1.0D0)*(V+M+k)/(k*(k+M))*(1.0D0+X)
            s = 0.0D0
            DO j = 1 , M
               s = s + ((k+j)**2+V*V)/((k+j)*((k+j)**2-V*V))
            ENDDO
            s2 = 0.0D0
            DO j = 1 , k
               s2 = s2 + 1.0D0/(j*(j*j-V*V))
            ENDDO
            pss = pa + s + 2.0D0*V*V*s2 - 1.0D0/(M+k-V)                 &
                & + DLOG(0.5D0*(1.0D0+X))
            r2 = pss*r
            Pmv = Pmv + r2
            IF ( DABS(r2/Pmv)<eps ) GOTO 100
         ENDDO
 100     Pmv = pv0 + Pmv*vs*c0
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE LPMV(V,M,X,Pmv)
!
!       =======================================================
!       Purpose: Compute the associated Legendre function
!                Pmv(x) with an integer order and an arbitrary
!                degree v, using recursion for large degrees
!       Input :  x   --- Argument of Pm(x)  ( -1 ≤ x ≤ 1 )
!                m   --- Order of Pmv(x)
!                v   --- Degree of Pmv(x)
!       Output:  PMV --- Pmv(x)
!       Routine called:  LPMV0
!       =======================================================
!
      IMPLICIT NONE
!*--LPMV8363
      DOUBLE PRECISION DINF , DNAN , g1 , g2 , p0 , p1 , Pmv , V , v0 , &
                     & vx , X
      INTEGER j , M , mx , neg_m , nv
      IF ( X==-1.0D0 .AND. V/=INT(V) ) THEN
         IF ( M==0 ) Pmv = -DINF()
         IF ( M/=0 ) Pmv = DINF()
         RETURN
      ENDIF
      vx = V
      mx = M
!       DLMF 14.9.5
      IF ( V<0 ) vx = -vx - 1
      neg_m = 0
      IF ( M<0 ) THEN
         IF ( (vx+M+1)>0D0 .OR. vx/=INT(vx) ) THEN
            neg_m = 1
            mx = -M
         ELSE
!             We don't handle cases where DLMF 14.9.3 doesn't help
            Pmv = DNAN()
            RETURN
         ENDIF
      ENDIF
      nv = INT(vx)
      v0 = vx - nv
      IF ( nv>2 .AND. nv>mx ) THEN
!          Up-recursion on degree, AMS 8.5.3 / DLMF 14.10.3
         CALL LPMV0(v0+mx,mx,X,p0)
         CALL LPMV0(v0+mx+1,mx,X,p1)
         Pmv = p1
         DO j = mx + 2 , nv
            Pmv = ((2*(v0+j)-1)*X*p1-(v0+j-1+mx)*p0)/(v0+j-mx)
            p0 = p1
            p1 = Pmv
         ENDDO
      ELSE
         CALL LPMV0(vx,mx,X,Pmv)
      ENDIF
      IF ( neg_m/=0 .AND. ABS(Pmv)<1.0D+300 ) THEN
!          DLMF 14.9.3
         CALL GAMMA2(vx-mx+1,g1)
         CALL GAMMA2(vx+mx+1,g2)
         Pmv = Pmv*g1/g2*(-1)**mx
      ENDIF
      END
 
 
!       **********************************
 
      SUBROUTINE CGAMA(X,Y,Kf,Gr,Gi)
!
!       =========================================================
!       Purpose: Compute the gamma function Г(z) or ln[Г(z)]
!                for a complex argument
!       Input :  x  --- Real part of z
!                y  --- Imaginary part of z
!                KF --- Function code
!                       KF=0 for ln[Г(z)]
!                       KF=1 for Г(z)
!       Output:  GR --- Real part of ln[Г(z)] or Г(z)
!                GI --- Imaginary part of ln[Г(z)] or Г(z)
!       ========================================================
!
      IMPLICIT NONE
!*--CGAMA8431
      DOUBLE PRECISION a , g0 , Gi , gi1 , Gr , gr1 , pi , si , sr , t ,&
                     & th , th1 , th2 , X , x0 , x1 , Y , y1 , z1 , z2
      INTEGER j , k , Kf , na
      DIMENSION a(10)
      pi = 3.141592653589793D0
      DATA a/8.333333333333333D-02 , -2.777777777777778D-03 ,           &
         & 7.936507936507937D-04 , -5.952380952380952D-04 ,             &
         & 8.417508417508418D-04 , -1.917526917526918D-03 ,             &
         & 6.410256410256410D-03 , -2.955065359477124D-02 ,             &
         & 1.796443723688307D-01 , -1.39243221690590D+00/
      IF ( Y==0.0D0 .AND. X==INT(X) .AND. X<=0.0D0 ) THEN
         Gr = 1.0D+300
         Gi = 0.0D0
         RETURN
      ELSEIF ( X<0.0D0 ) THEN
         x1 = X
         y1 = Y
         X = -X
         Y = -Y
      ELSE
         y1 = 0.0D0
         x1 = X
      ENDIF
      x0 = X
      na = 0
      IF ( X<=7.0 ) THEN
         na = INT(7-X)
         x0 = X + na
      ENDIF
      z1 = DSQRT(x0*x0+Y*Y)
      th = DATAN(Y/x0)
      Gr = (x0-.5D0)*DLOG(z1) - th*Y - x0 + 0.5D0*DLOG(2.0D0*pi)
      Gi = th*(x0-0.5D0) + Y*DLOG(z1) - Y
      DO k = 1 , 10
         t = z1**(1-2*k)
         Gr = Gr + a(k)*t*DCOS((2.0D0*k-1.0D0)*th)
         Gi = Gi - a(k)*t*DSIN((2.0D0*k-1.0D0)*th)
      ENDDO
      IF ( X<=7.0 ) THEN
         gr1 = 0.0D0
         gi1 = 0.0D0
         DO j = 0 , na - 1
            gr1 = gr1 + .5D0*DLOG((X+j)**2+Y*Y)
            gi1 = gi1 + DATAN(Y/(X+j))
         ENDDO
         Gr = Gr - gr1
         Gi = Gi - gi1
      ENDIF
      IF ( x1<0.0D0 ) THEN
         z1 = DSQRT(X*X+Y*Y)
         th1 = DATAN(Y/X)
         sr = -DSIN(pi*X)*DCOSH(pi*Y)
         si = -DCOS(pi*X)*DSINH(pi*Y)
         z2 = DSQRT(sr*sr+si*si)
         th2 = DATAN(si/sr)
         IF ( sr<0.0D0 ) th2 = pi + th2
         Gr = DLOG(pi/(z1*z2)) - Gr
         Gi = -th1 - th2 - Gi
         X = x1
         Y = y1
      ENDIF
      IF ( Kf==1 ) THEN
         g0 = EXP(Gr)
         Gr = g0*DCOS(Gi)
         Gi = g0*DSIN(Gi)
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE ASWFB(M,N,C,X,Kd,Cv,S1f,S1d)
!
!       ===========================================================
!       Purpose: Compute the prolate and oblate spheroidal angular
!                functions of the first kind and their derivatives
!       Input :  m  --- Mode parameter,  m = 0,1,2,...
!                n  --- Mode parameter,  n = m,m+1,...
!                c  --- Spheroidal parameter
!                x  --- Argument of angular function, |x| < 1.0
!                KD --- Function code
!                       KD=1 for prolate;  KD=-1 for oblate
!                cv --- Characteristic value
!       Output:  S1F --- Angular function of the first kind
!                S1D --- Derivative of the angular function of
!                        the first kind
!       Routines called:
!            (1) SDMN for computing expansion coefficients dk
!            (2) LPMNS for computing associated Legendre function
!                of the first kind Pmn(x)
!       ===========================================================
!
      IMPLICIT NONE
!*--ASWFB8527
      DOUBLE PRECISION C , Cv , df , eps , pd , pm , S1d , S1f , su1 ,  &
                     & sw , X
      INTEGER ip , k , Kd , M , mk , N , nm , nm2
      DIMENSION df(200) , pm(0:251) , pd(0:251)
      eps = 1.0D-14
      ip = 1
      IF ( N-M==2*INT((N-M)/2) ) ip = 0
      nm = 25 + INT((N-M)/2+C)
      nm2 = 2*nm + M
      CALL SDMN(M,N,C,Cv,Kd,df)
      CALL LPMNS(M,nm2,X,pm,pd)
      sw = 0.0D0
      su1 = 0.0D0
      DO k = 1 , nm
         mk = M + 2*(k-1) + ip
         su1 = su1 + df(k)*pm(mk)
         IF ( DABS(sw-su1)<DABS(su1)*eps ) GOTO 100
         sw = su1
      ENDDO
 100  S1f = (-1)**M*su1
      su1 = 0.0D0
      DO k = 1 , nm
         mk = M + 2*(k-1) + ip
         su1 = su1 + df(k)*pd(mk)
         IF ( DABS(sw-su1)<DABS(su1)*eps ) GOTO 200
         sw = su1
      ENDDO
 200  S1d = (-1)**M*su1
      END
 
 
 
!       **********************************
 
      SUBROUTINE CHGUS(A,B,X,Hu,Id)
!
!       ======================================================
!       Purpose: Compute confluent hypergeometric function
!                U(a,b,x) for small argument x
!       Input  : a  --- Parameter
!                b  --- Parameter ( b <> 0,-1,-2,...)
!                x  --- Argument
!       Output:  HU --- U(a,b,x)
!                ID --- Estimated number of significant digits
!       Routine called: GAMMA2 for computing gamma function
!       ======================================================
!
!       DLMF 13.2.42 with prefactors rewritten according to
!       DLMF 5.5.3, M(a, b, x) with DLMF 13.2.2
!
      IMPLICIT NONE
!*--CHGUS8582
      DOUBLE PRECISION A , B , d1 , d2 , ga , gab , gb , gb2 , h0 ,     &
                     & hmax , hmin , Hu , hu0 , hua , pi , r1 , r2 , X ,&
                     & xg1 , xg2
      INTEGER Id , j
      Id = -100
      pi = 3.141592653589793D0
      CALL GAMMA2(A,ga)
      CALL GAMMA2(B,gb)
      xg1 = 1.0D0 + A - B
      CALL GAMMA2(xg1,gab)
      xg2 = 2.0D0 - B
      CALL GAMMA2(xg2,gb2)
      hu0 = pi/DSIN(pi*B)
      r1 = hu0/(gab*gb)
      r2 = hu0*X**(1.0D0-B)/(ga*gb2)
      Hu = r1 - r2
      hmax = 0.0D0
      hmin = 1.0D+300
      h0 = 0.0D0
      DO j = 1 , 150
         r1 = r1*(A+j-1.0D0)/(j*(B+j-1.0D0))*X
         r2 = r2*(A-B+j)/(j*(1.0D0-B+j))*X
         Hu = Hu + r1 - r2
         hua = DABS(Hu)
         IF ( hua>hmax ) hmax = hua
         IF ( hua<hmin ) hmin = hua
         IF ( DABS(Hu-h0)<DABS(Hu)*1.0D-15 ) GOTO 100
         h0 = Hu
      ENDDO
 100  d1 = LOG10(hmax)
      d2 = 0.0D0
      IF ( hmin/=0.0 ) d2 = LOG10(hmin)
      Id = 15 - ABS(d1-d2)
      END
 
 
 
!       **********************************
 
      SUBROUTINE ITTH0(X,Tth)
!
!       ===========================================================
!       Purpose: Evaluate the integral H0(t)/t with respect to t
!                from x to infinity
!       Input :  x   --- Lower limit  ( x ≥ 0 )
!       Output:  TTH --- Integration of H0(t)/t from x to infinity
!       ===========================================================
!
      IMPLICIT NONE
!*--ITTH08635
      DOUBLE PRECISION f0 , g0 , pi , r , s , t , Tth , tty , X , xt
      INTEGER k
      pi = 3.141592653589793D0
      s = 1.0D0
      r = 1.0D0
      IF ( X<24.5D0 ) THEN
         DO k = 1 , 60
            r = -r*X*X*(2.0*k-1.0D0)/(2.0*k+1.0D0)**3
            s = s + r
            IF ( DABS(r)<DABS(s)*1.0D-12 ) GOTO 50
         ENDDO
 50      Tth = pi/2.0D0 - 2.0D0/pi*X*s
      ELSE
         DO k = 1 , 10
            r = -r*(2.0*k-1.0D0)**3/((2.0*k+1.0D0)*X*X)
            s = s + r
            IF ( DABS(r)<DABS(s)*1.0D-12 ) GOTO 100
         ENDDO
 100     Tth = 2.0D0/(pi*X)*s
         t = 8.0D0/X
         xt = X + .25D0*pi
         f0 = (((((.18118D-2*t-.91909D-2)*t+.017033D0)*t-.9394D-3)      &
            & *t-.051445D0)*t-.11D-5)*t + .7978846D0
         g0 = (((((-.23731D-2*t+.59842D-2)*t+.24437D-2)*t-.0233178D0)   &
            & *t+.595D-4)*t+.1620695D0)*t
         tty = (f0*DSIN(xt)-g0*DCOS(xt))/(DSQRT(X)*X)
         Tth = Tth + tty
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE LGAMA(Kf,X,Gl)
!
!       ==================================================
!       Purpose: Compute gamma function Г(x) or ln[Г(x)]
!       Input:   x  --- Argument of Г(x) ( x > 0 )
!                KF --- Function code
!                       KF=1 for Г(x); KF=0 for ln[Г(x)]
!       Output:  GL --- Г(x) or ln[Г(x)]
!       ==================================================
!
      IMPLICIT NONE
!*--LGAMA8682
      DOUBLE PRECISION a , Gl , gl0 , X , x0 , x2 , xp
      INTEGER k , Kf , n
      DIMENSION a(10)
      DATA a/8.333333333333333D-02 , -2.777777777777778D-03 ,           &
         & 7.936507936507937D-04 , -5.952380952380952D-04 ,             &
         & 8.417508417508418D-04 , -1.917526917526918D-03 ,             &
         & 6.410256410256410D-03 , -2.955065359477124D-02 ,             &
         & 1.796443723688307D-01 , -1.39243221690590D+00/
      x0 = X
      n = 0
      IF ( X==1.0 .OR. X==2.0 ) THEN
         Gl = 0.0D0
         GOTO 100
      ELSEIF ( X<=7.0 ) THEN
         n = INT(7-X)
         x0 = X + n
      ENDIF
      x2 = 1.0D0/(x0*x0)
      xp = 6.283185307179586477D0
      gl0 = a(10)
      DO k = 9 , 1 , -1
         gl0 = gl0*x2 + a(k)
      ENDDO
      Gl = gl0/x0 + 0.5D0*DLOG(xp) + (x0-.5D0)*DLOG(x0) - x0
      IF ( X<=7.0 ) THEN
         DO k = 1 , n
            Gl = Gl - DLOG(x0-1.0D0)
            x0 = x0 - 1.0D0
         ENDDO
      ENDIF
 100  IF ( Kf==1 ) Gl = EXP(Gl)
      END
 
!       **********************************
 
      SUBROUTINE LQNA(N,X,Qn,Qd)
!
!       =====================================================
!       Purpose: Compute Legendre functions Qn(x) and Qn'(x)
!       Input :  x  --- Argument of Qn(x) ( -1 ≤ x ≤ 1 )
!                n  --- Degree of Qn(x) ( n = 0,1,2,… )
!       Output:  QN(n) --- Qn(x)
!                QD(n) --- Qn'(x)
!                ( 1.0D+300 stands for infinity )
!       =====================================================
!
      IMPLICIT NONE
!*--LQNA8733
      INTEGER k , N
      DOUBLE PRECISION q0 , q1 , Qd , qf , Qn , X
      DIMENSION Qn(0:N) , Qd(0:N)
      IF ( DABS(X)==1.0D0 ) THEN
         DO k = 0 , N
            Qn(k) = 1.0D+300
            Qd(k) = -1.0D+300
         ENDDO
      ELSEIF ( DABS(X)<1.0D0 ) THEN
         q0 = 0.5D0*DLOG((1.0D0+X)/(1.0D0-X))
         q1 = X*q0 - 1.0D0
         Qn(0) = q0
         Qn(1) = q1
         Qd(0) = 1.0D0/(1.0D0-X*X)
         Qd(1) = Qn(0) + X*Qd(0)
         DO k = 2 , N
            qf = ((2*k-1)*X*q1-(k-1)*q0)/k
            Qn(k) = qf
            Qd(k) = (Qn(k-1)-X*qf)*k/(1.0D0-X*X)
            q0 = q1
            q1 = qf
         ENDDO
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE DVLA(Va,X,Pd)
!
!       ====================================================
!       Purpose: Compute parabolic cylinder functions Dv(x)
!                for large argument
!       Input:   x  --- Argument
!                va --- Order
!       Output:  PD --- Dv(x)
!       Routines called:
!             (1) VVLA for computing Vv(x) for large |x|
!             (2) GAMMA2 for computing Г(x)
!       ====================================================
!
      IMPLICIT NONE
!*--DVLA8778
      DOUBLE PRECISION a0 , ep , eps , gl , Pd , pi , r , Va , vl , X , &
                     & x1
      INTEGER k
      pi = 3.141592653589793D0
      eps = 1.0D-12
      ep = EXP(-.25*X*X)
      a0 = DABS(X)**Va*ep
      r = 1.0D0
      Pd = 1.0D0
      DO k = 1 , 16
         r = -0.5D0*r*(2.0*k-Va-1.0)*(2.0*k-Va-2.0)/(k*X*X)
         Pd = Pd + r
         IF ( DABS(r/Pd)<eps ) GOTO 100
      ENDDO
 100  Pd = a0*Pd
      IF ( X<0.0D0 ) THEN
         x1 = -X
         CALL VVLA(Va,x1,vl)
         CALL GAMMA2(-Va,gl)
         Pd = pi*vl/gl + DCOS(pi*Va)*Pd
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE IK01A(X,Bi0,Di0,Bi1,Di1,Bk0,Dk0,Bk1,Dk1)
!
!       =========================================================
!       Purpose: Compute modified Bessel functions I0(x), I1(1),
!                K0(x) and K1(x), and their derivatives
!       Input :  x   --- Argument ( x ≥ 0 )
!       Output:  BI0 --- I0(x)
!                DI0 --- I0'(x)
!                BI1 --- I1(x)
!                DI1 --- I1'(x)
!                BK0 --- K0(x)
!                DK0 --- K0'(x)
!                BK1 --- K1(x)
!                DK1 --- K1'(x)
!       =========================================================
!
      IMPLICIT NONE
!*--IK01A8826
      DOUBLE PRECISION a , a1 , b , Bi0 , Bi1 , Bk0 , Bk1 , ca , cb ,   &
                     & ct , Di0 , Di1 , Dk0 , Dk1 , el , pi , r , w0 ,  &
                     & ww , X
      DOUBLE PRECISION x2 , xr , xr2
      INTEGER k , k0
      DIMENSION a(12) , b(12) , a1(8)
      pi = 3.141592653589793D0
      el = 0.5772156649015329D0
      x2 = X*X
      IF ( X==0.0D0 ) THEN
         Bi0 = 1.0D0
         Bi1 = 0.0D0
         Bk0 = 1.0D+300
         Bk1 = 1.0D+300
         Di0 = 0.0D0
         Di1 = 0.5D0
         Dk0 = -1.0D+300
         Dk1 = -1.0D+300
         RETURN
      ELSEIF ( X<=18.0D0 ) THEN
         Bi0 = 1.0D0
         r = 1.0D0
         DO k = 1 , 50
            r = 0.25D0*r*x2/(k*k)
            Bi0 = Bi0 + r
            IF ( DABS(r/Bi0)<1.0D-15 ) GOTO 50
         ENDDO
 50      Bi1 = 1.0D0
         r = 1.0D0
         DO k = 1 , 50
            r = 0.25D0*r*x2/(k*(k+1))
            Bi1 = Bi1 + r
            IF ( DABS(r/Bi1)<1.0D-15 ) GOTO 100
         ENDDO
 100     Bi1 = 0.5D0*X*Bi1
      ELSE
         DATA a/0.125D0 , 7.03125D-2 , 7.32421875D-2 ,                  &
            & 1.1215209960938D-1 , 2.2710800170898D-1 ,                 &
            & 5.7250142097473D-1 , 1.7277275025845D0 ,                  &
            & 6.0740420012735D0 , 2.4380529699556D01 ,                  &
            & 1.1001714026925D02 , 5.5133589612202D02 ,                 &
            & 3.0380905109224D03/
         DATA b/ - 0.375D0 , -1.171875D-1 , -1.025390625D-1 ,           &
            & -1.4419555664063D-1 , -2.7757644653320D-1 ,               &
            & -6.7659258842468D-1 , -1.9935317337513D0 ,                &
            & -6.8839142681099D0 , -2.7248827311269D01 ,                &
            & -1.2159789187654D02 , -6.0384407670507D02 ,               &
            & -3.3022722944809D03/
         k0 = 12
         IF ( X>=35.0 ) k0 = 9
         IF ( X>=50.0 ) k0 = 7
         ca = EXP(X)/DSQRT(2.0D0*pi*X)
         Bi0 = 1.0D0
         xr = 1.0D0/X
         DO k = 1 , k0
            Bi0 = Bi0 + a(k)*xr**k
         ENDDO
         Bi0 = ca*Bi0
         Bi1 = 1.0D0
         DO k = 1 , k0
            Bi1 = Bi1 + b(k)*xr**k
         ENDDO
         Bi1 = ca*Bi1
      ENDIF
      ww = 0.0D0
      IF ( X<=9.0D0 ) THEN
         ct = -(DLOG(X/2.0D0)+el)
         Bk0 = 0.0D0
         w0 = 0.0D0
         r = 1.0D0
         DO k = 1 , 50
            w0 = w0 + 1.0D0/k
            r = 0.25D0*r/(k*k)*x2
            Bk0 = Bk0 + r*(w0+ct)
            IF ( DABS((Bk0-ww)/Bk0)<1.0D-15 ) GOTO 150
            ww = Bk0
         ENDDO
 150     Bk0 = Bk0 + ct
      ELSE
         DATA a1/0.125D0 , 0.2109375D0 , 1.0986328125D0 ,               &
            & 1.1775970458984D01 , 2.1461706161499D02 ,                 &
            & 5.9511522710323D03 , 2.3347645606175D05 ,                 &
            & 1.2312234987631D07/
         cb = 0.5D0/X
         xr2 = 1.0D0/x2
         Bk0 = 1.0D0
         DO k = 1 , 8
            Bk0 = Bk0 + a1(k)*xr2**k
         ENDDO
         Bk0 = cb*Bk0/Bi0
      ENDIF
      Bk1 = (1.0D0/X-Bi1*Bk0)/Bi0
      Di0 = Bi1
      Di1 = Bi0 - Bi1/X
      Dk0 = -Bk1
      Dk1 = -Bk0 - Bk1/X
      END
 
!       **********************************
 
      SUBROUTINE CPBDN(N,Z,Cpb,Cpd)
!
!       ==================================================
!       Purpose: Compute the parabolic cylinder functions
!                 Dn(z) and Dn'(z) for a complex argument
!       Input:   z --- Complex argument of Dn(z)
!                n --- Order of Dn(z)  ( n=0,±1,±2,… )
!       Output:  CPB(|n|) --- Dn(z)
!                CPD(|n|) --- Dn'(z)
!       Routines called:
!            (1) CPDSA for computing Dn(z) for a small |z|
!            (2) CPDLA for computing Dn(z) for a large |z|
!       ==================================================
!
      IMPLICIT NONE
!*--CPBDN8945
      DOUBLE PRECISION a0 , pi , x
      COMPLEX*16 c0 , ca0 , cf , cf0 , cf1 , cfa , cfb , Cpb , Cpd ,    &
               & cs0 , Z , z1
      INTEGER k , m , N , n0 , n1 , nm1
      DIMENSION Cpb(0:*) , Cpd(0:*)
      pi = 3.141592653589793D0
      x = DBLE(Z)
      a0 = ABS(Z)
      c0 = (0.0D0,0.0D0)
      ca0 = EXP(-0.25D0*Z*Z)
      n0 = 0
      IF ( N>=0 ) THEN
         cf0 = ca0
         cf1 = Z*ca0
         Cpb(0) = cf0
         Cpb(1) = cf1
         DO k = 2 , N
            cf = Z*cf1 - (k-1.0D0)*cf0
            Cpb(k) = cf
            cf0 = cf1
            cf1 = cf
         ENDDO
      ELSE
         n0 = -N
         IF ( x<=0.0 .OR. ABS(Z)==0.0 ) THEN
            cf0 = ca0
            Cpb(0) = cf0
            z1 = -Z
            IF ( a0<=7.0 ) THEN
               CALL CPDSA(-1,z1,cf1)
            ELSE
               CALL CPDLA(-1,z1,cf1)
            ENDIF
            cf1 = DSQRT(2.0D0*pi)/ca0 - cf1
            Cpb(1) = cf1
            DO k = 2 , n0
               cf = (-Z*cf1+cf0)/(k-1.0D0)
               Cpb(k) = cf
               cf0 = cf1
               cf1 = cf
            ENDDO
         ELSEIF ( a0<=3.0 ) THEN
            CALL CPDSA(-n0,Z,cfa)
            Cpb(n0) = cfa
            n1 = n0 + 1
            CALL CPDSA(-n1,Z,cfb)
            Cpb(n1) = cfb
            nm1 = n0 - 1
            DO k = nm1 , 0 , -1
               cf = Z*cfa + (k+1.0D0)*cfb
               Cpb(k) = cf
               cfb = cfa
               cfa = cf
            ENDDO
         ELSE
            m = 100 + ABS(N)
            cfa = c0
            cfb = (1.0D-30,0.0D0)
            DO k = m , 0 , -1
               cf = Z*cfb + (k+1.0D0)*cfa
               IF ( k<=n0 ) Cpb(k) = cf
               cfa = cfb
               cfb = cf
            ENDDO
            cs0 = ca0/cf
            DO k = 0 , n0
               Cpb(k) = cs0*Cpb(k)
            ENDDO
         ENDIF
      ENDIF
      Cpd(0) = -0.5D0*Z*Cpb(0)
      IF ( N>=0 ) THEN
         DO k = 1 , N
            Cpd(k) = -0.5D0*Z*Cpb(k) + k*Cpb(k-1)
         ENDDO
      ELSE
         DO k = 1 , n0
            Cpd(k) = 0.5D0*Z*Cpb(k) - Cpb(k-1)
         ENDDO
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE IK01B(X,Bi0,Di0,Bi1,Di1,Bk0,Dk0,Bk1,Dk1)
!
!       =========================================================
!       Purpose: Compute modified Bessel functions I0(x), I1(1),
!                K0(x) and K1(x), and their derivatives
!       Input :  x   --- Argument ( x ≥ 0 )
!       Output:  BI0 --- I0(x)
!                DI0 --- I0'(x)
!                BI1 --- I1(x)
!                DI1 --- I1'(x)
!                BK0 --- K0(x)
!                DK0 --- K0'(x)
!                BK1 --- K1(x)
!                DK1 --- K1'(x)
!       =========================================================
!
      IMPLICIT NONE
!*--IK01B9052
      DOUBLE PRECISION Bi0 , Bi1 , Bk0 , Bk1 , Di0 , Di1 , Dk0 , Dk1 ,  &
                     & t , t2 , X
      IF ( X==0.0D0 ) THEN
         Bi0 = 1.0D0
         Bi1 = 0.0D0
         Bk0 = 1.0D+300
         Bk1 = 1.0D+300
         Di0 = 0.0D0
         Di1 = 0.5D0
         Dk0 = -1.0D+300
         Dk1 = -1.0D+300
         RETURN
      ELSEIF ( X<=3.75D0 ) THEN
         t = X/3.75D0
         t2 = t*t
         Bi0 = (((((.0045813D0*t2+.0360768D0)*t2+.2659732D0)*t2+        &
             & 1.2067492D0)*t2+3.0899424D0)*t2+3.5156229D0)*t2 + 1.0D0
         Bi1 = X*((((((.00032411D0*t2+.00301532D0)*t2+.02658733D0)*t2+  &
             & .15084934D0)*t2+.51498869D0)*t2+.87890594D0)*t2+.5D0)
      ELSE
         t = 3.75D0/X
         Bi0 = ((((((((.00392377D0*t-.01647633D0)*t+.02635537D0)*t-     &
             & .02057706D0)*t+.916281D-2)*t-.157565D-2)*t+.225319D-2)   &
             & *t+.01328592D0)*t+.39894228D0)*EXP(X)/DSQRT(X)
         Bi1 = ((((((((-.420059D-2*t+.01787654D0)*t-.02895312D0)*t+     &
             & .02282967D0)*t-.01031555D0)*t+.163801D-2)*t-.00362018D0) &
             & *t-.03988024D0)*t+.39894228D0)*EXP(X)/DSQRT(X)
      ENDIF
      IF ( X<=2.0D0 ) THEN
         t = X/2.0D0
         t2 = t*t
         Bk0 = (((((.0000074D0*t2+.0001075D0)*t2+.00262698D0)*t2+       &
             & .0348859D0)*t2+.23069756D0)*t2+.4227842D0)               &
             & *t2 - .57721566D0 - Bi0*DLOG(t)
         Bk1 = ((((((-.00004686D0*t2-.00110404D0)*t2-.01919402D0)*t2-   &
             & .18156897D0)*t2-.67278579D0)*t2+.15443144D0)*t2+1.0D0)   &
             & /X + Bi1*DLOG(t)
      ELSE
         t = 2.0D0/X
         t2 = t*t
         Bk0 = ((((((.00053208D0*t-.0025154D0)*t+.00587872D0)*t-        &
             & .01062446D0)*t+.02189568D0)*t-.07832358D0)               &
             & *t+1.25331414D0)*EXP(-X)/DSQRT(X)
         Bk1 = ((((((-.00068245D0*t+.00325614D0)*t-.00780353D0)*t+      &
             & .01504268D0)*t-.0365562D0)*t+.23498619D0)*t+1.25331414D0)&
             & *EXP(-X)/DSQRT(X)
      ENDIF
      Di0 = Bi1
      Di1 = Bi0 - Bi1/X
      Dk0 = -Bk1
      Dk1 = -Bk0 - Bk1/X
      END
 
!       **********************************
 
      SUBROUTINE BETA(P,Q,Bt)
!
!       ==========================================
!       Purpose: Compute the beta function B(p,q)
!       Input :  p  --- Parameter  ( p > 0 )
!                q  --- Parameter  ( q > 0 )
!       Output:  BT --- B(p,q)
!       Routine called: GAMMA2 for computing Г(x)
!       ==========================================
!
      IMPLICIT NONE
!*--BETA9122
      DOUBLE PRECISION Bt , gp , gpq , gq , P , ppq , Q
      CALL GAMMA2(P,gp)
      CALL GAMMA2(Q,gq)
      ppq = P + Q
      CALL GAMMA2(ppq,gpq)
      Bt = gp*gq/gpq
      END
 
 
 
!       **********************************
 
      SUBROUTINE LPN(N,X,Pn,Pd)
!
!       ===============================================
!       Purpose: Compute Legendre polynomials Pn(x)
!                and their derivatives Pn'(x)
!       Input :  x --- Argument of Pn(x)
!                n --- Degree of Pn(x) ( n = 0,1,...)
!       Output:  PN(n) --- Pn(x)
!                PD(n) --- Pn'(x)
!       ===============================================
!
      IMPLICIT NONE
!*--LPN9150
      INTEGER k , N
      DOUBLE PRECISION p0 , p1 , Pd , pf , Pn , X
      DIMENSION Pn(0:N) , Pd(0:N)
      Pn(0) = 1.0D0
      Pn(1) = X
      Pd(0) = 0.0D0
      Pd(1) = 1.0D0
      p0 = 1.0D0
      p1 = X
      DO k = 2 , N
         pf = (2.0D0*k-1.0D0)/k*X*p1 - (k-1.0D0)/k*p0
         Pn(k) = pf
         IF ( DABS(X)==1.0D0 ) THEN
            Pd(k) = 0.5D0*X**(k+1)*k*(k+1.0D0)
         ELSE
            Pd(k) = k*(p1-X*pf)/(1.0D0-X*X)
         ENDIF
         p0 = p1
         p1 = pf
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE FCOEF(Kd,M,Q,A,Fc)
!
!       =====================================================
!       Purpose: Compute expansion coefficients for Mathieu
!                functions and modified Mathieu functions
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!                KD --- Case code
!                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
!                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
!                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
!                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
!                A  --- Characteristic value of Mathieu
!                       functions for given m and q
!       Output:  FC(k) --- Expansion coefficients of Mathieu
!                       functions ( k= 1,2,...,KM )
!                       FC(1),FC(2),FC(3),... correspond to
!                       A0,A2,A4,... for KD=1 case, A1,A3,
!                       A5,... for KD=2 case, B1,B3,B5,...
!                       for KD=3 case and B2,B4,B6,... for
!                       KD=4 case
!       =====================================================
!
      IMPLICIT NONE
!*--FCOEF9202
      DOUBLE PRECISION A , DNAN , f , f1 , f2 , f3 , Fc , fnan , Q ,    &
                     & qm , s , s0 , sp , ss , u , v
      INTEGER i , j , jm , k , kb , Kd , km , M
      DIMENSION Fc(251)
      DO i = 1 , 251
         Fc(i) = 0.0D0
      ENDDO
      IF ( DABS(Q)<=1.0D-7 ) THEN
!          Expansion up to order Q^1 (Abramowitz & Stegun 20.2.27-28)
         IF ( Kd==1 ) THEN
            jm = M/2 + 1
         ELSEIF ( Kd==2 .OR. Kd==3 ) THEN
            jm = (M-1)/2 + 1
         ELSEIF ( Kd==4 ) THEN
            jm = M/2
         ENDIF
!          Check for overflow
         IF ( jm+1>251 ) THEN
            fnan = DNAN()
            DO i = 1 , 251
               Fc(i) = fnan
            ENDDO
            RETURN
         ENDIF
!          Proceed using the simplest expansion
         IF ( Kd==1 .OR. Kd==2 ) THEN
            IF ( M==0 ) THEN
               Fc(1) = 1/SQRT(2.0D0)
               Fc(2) = -Q/2.0D0/SQRT(2.0D0)
            ELSEIF ( M==1 ) THEN
               Fc(1) = 1.0D0
               Fc(2) = -Q/8.0D0
            ELSEIF ( M==2 ) THEN
               Fc(1) = Q/4.0D0
               Fc(2) = 1.0D0
               Fc(3) = -Q/12.0D0
            ELSE
               Fc(jm) = 1.0D0
               Fc(jm+1) = -Q/(4.0D0*(M+1))
               Fc(jm-1) = Q/(4.0D0*(M-1))
            ENDIF
         ELSEIF ( Kd==3 .OR. Kd==4 ) THEN
            IF ( M==1 ) THEN
               Fc(1) = 1.0D0
               Fc(2) = -Q/8.0D0
            ELSEIF ( M==2 ) THEN
               Fc(1) = 1.0D0
               Fc(2) = -Q/12.0D0
            ELSE
               Fc(jm) = 1.0D0
               Fc(jm+1) = -Q/(4.0D0*(M+1))
               Fc(jm-1) = Q/(4.0D0*(M-1))
            ENDIF
         ENDIF
         RETURN
      ELSEIF ( Q<=1.0D0 ) THEN
         qm = 7.5 + 56.1*SQRT(Q) - 134.7*Q + 90.7*SQRT(Q)*Q
      ELSE
         qm = 17.0 + 3.1*SQRT(Q) - .126*Q + .0037*SQRT(Q)*Q
      ENDIF
      km = INT(qm+0.5*M)
      IF ( km>251 ) THEN
!          Overflow, generate NaNs
         fnan = DNAN()
         DO i = 1 , 251
            Fc(i) = fnan
         ENDDO
         RETURN
      ENDIF
      kb = 0
      s = 0.0D0
      f = 1.0D-100
      u = 0.0D0
      Fc(km) = 0.0D0
      f2 = 0.0D0
      IF ( Kd==1 ) THEN
         DO k = km , 3 , -1
            v = u
            u = f
            f = (A-4.0D0*k*k)*u/Q - v
            IF ( DABS(f)<DABS(Fc(k+1)) ) THEN
               kb = k
               Fc(1) = 1.0D-100
               sp = 0.0D0
               f3 = Fc(k+1)
               Fc(2) = A/Q*Fc(1)
               Fc(3) = (A-4.0D0)*Fc(2)/Q - 2.0D0*Fc(1)
               u = Fc(2)
               f1 = Fc(3)
               DO i = 3 , kb
                  v = u
                  u = f1
                  f1 = (A-4.0D0*(i-1.0D0)**2)*u/Q - v
                  Fc(i+1) = f1
                  IF ( i==kb ) f2 = f1
                  IF ( i/=kb ) sp = sp + f1*f1
               ENDDO
               sp = sp + 2.0D0*Fc(1)**2 + Fc(2)**2 + Fc(3)**2
               ss = s + sp*(f3/f2)**2
               s0 = DSQRT(1.0D0/ss)
               DO j = 1 , km
                  IF ( j<=kb+1 ) THEN
                     Fc(j) = s0*Fc(j)*f3/f2
                  ELSE
                     Fc(j) = s0*Fc(j)
                  ENDIF
               ENDDO
               GOTO 200
            ELSE
               Fc(k) = f
               s = s + f*f
            ENDIF
         ENDDO
         Fc(2) = Q*Fc(3)/(A-4.0D0-2.0D0*Q*Q/A)
         Fc(1) = Q/A*Fc(2)
         s = s + 2.0D0*Fc(1)**2 + Fc(2)**2
         s0 = DSQRT(1.0D0/s)
         DO k = 1 , km
            Fc(k) = s0*Fc(k)
         ENDDO
      ELSEIF ( Kd==2 .OR. Kd==3 ) THEN
         DO k = km , 3 , -1
            v = u
            u = f
            f = (A-(2.0D0*k-1)**2)*u/Q - v
            IF ( DABS(f)>=DABS(Fc(k)) ) THEN
               Fc(k-1) = f
               s = s + f*f
            ELSE
               kb = k
               f3 = Fc(k)
               GOTO 50
            ENDIF
         ENDDO
         Fc(1) = Q/(A-1.0D0-(-1)**Kd*Q)*Fc(2)
         s = s + Fc(1)*Fc(1)
         s0 = DSQRT(1.0D0/s)
         DO k = 1 , km
            Fc(k) = s0*Fc(k)
         ENDDO
         GOTO 200
 50      Fc(1) = 1.0D-100
         Fc(2) = (A-1.0D0-(-1)**Kd*Q)/Q*Fc(1)
         sp = 0.0D0
         u = Fc(1)
         f1 = Fc(2)
         DO i = 2 , kb - 1
            v = u
            u = f1
            f1 = (A-(2.0D0*i-1.0D0)**2)*u/Q - v
            IF ( i/=kb-1 ) THEN
               Fc(i+1) = f1
               sp = sp + f1*f1
            ELSE
               f2 = f1
            ENDIF
         ENDDO
         sp = sp + Fc(1)**2 + Fc(2)**2
         ss = s + sp*(f3/f2)**2
         s0 = 1.0D0/DSQRT(ss)
         DO j = 1 , km
            IF ( j<kb ) Fc(j) = s0*Fc(j)*f3/f2
            IF ( j>=kb ) Fc(j) = s0*Fc(j)
         ENDDO
      ELSEIF ( Kd==4 ) THEN
         DO k = km , 3 , -1
            v = u
            u = f
            f = (A-4.0D0*k*k)*u/Q - v
            IF ( DABS(f)>=DABS(Fc(k)) ) THEN
               Fc(k-1) = f
               s = s + f*f
            ELSE
               kb = k
               f3 = Fc(k)
               GOTO 100
            ENDIF
         ENDDO
         Fc(1) = Q/(A-4.0D0)*Fc(2)
         s = s + Fc(1)*Fc(1)
         s0 = DSQRT(1.0D0/s)
         DO k = 1 , km
            Fc(k) = s0*Fc(k)
         ENDDO
         GOTO 200
 100     Fc(1) = 1.0D-100
         Fc(2) = (A-4.0D0)/Q*Fc(1)
         sp = 0.0D0
         u = Fc(1)
         f1 = Fc(2)
         DO i = 2 , kb - 1
            v = u
            u = f1
            f1 = (A-4.0D0*i*i)*u/Q - v
            IF ( i/=kb-1 ) THEN
               Fc(i+1) = f1
               sp = sp + f1*f1
            ELSE
               f2 = f1
            ENDIF
         ENDDO
         sp = sp + Fc(1)**2 + Fc(2)**2
         ss = s + sp*(f3/f2)**2
         s0 = 1.0D0/DSQRT(ss)
         DO j = 1 , km
            IF ( j<kb ) Fc(j) = s0*Fc(j)*f3/f2
            IF ( j>=kb ) Fc(j) = s0*Fc(j)
         ENDDO
      ENDIF
 200  IF ( Fc(1)<0.0D0 ) THEN
         DO j = 1 , km
            Fc(j) = -Fc(j)
         ENDDO
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE SPHI(N,X,Nm,Si,Di)
!
!       ========================================================
!       Purpose: Compute modified spherical Bessel functions
!                of the first kind, in(x) and in'(x)
!       Input :  x --- Argument of in(x)
!                n --- Order of in(x) ( n = 0,1,2,... )
!       Output:  SI(n) --- in(x)
!                DI(n) --- in'(x)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ========================================================
!
      IMPLICIT NONE
!*--SPHI9442
      DOUBLE PRECISION cs , Di , f , f0 , f1 , Si , si0 , X
      INTEGER k , m , MSTA1 , MSTA2 , N , Nm
      DIMENSION Si(0:N) , Di(0:N)
      Nm = N
      IF ( DABS(X)<1.0D-100 ) THEN
         DO k = 0 , N
            Si(k) = 0.0D0
            Di(k) = 0.0D0
         ENDDO
         Si(0) = 1.0D0
         Di(1) = 0.333333333333333D0
         RETURN
      ENDIF
      Si(0) = DSINH(X)/X
      Si(1) = -(DSINH(X)/X-DCOSH(X))/X
      si0 = Si(0)
      IF ( N>=2 ) THEN
         m = MSTA1(X,200)
         IF ( m<N ) THEN
            Nm = m
         ELSE
            m = MSTA2(X,N,15)
         ENDIF
         f = 0.0D0
         f0 = 0.0D0
         f1 = 1.0D0 - 100
         DO k = m , 0 , -1
            f = (2.0D0*k+3.0D0)*f1/X + f0
            IF ( k<=Nm ) Si(k) = f
            f0 = f1
            f1 = f
         ENDDO
         cs = si0/f
         DO k = 0 , Nm
            Si(k) = cs*Si(k)
         ENDDO
      ENDIF
      Di(0) = Si(1)
      DO k = 1 , Nm
         Di(k) = Si(k-1) - (k+1.0D0)/X*Si(k)
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE PBWA(A,X,W1f,W1d,W2f,W2d)
!
!       ======================================================
!       Purpose: Compute parabolic cylinder functions W(a,±x)
!                and their derivatives
!       Input  : a --- Parameter  ( 0 ≤ |a| ≤ 5 )
!                x --- Argument of W(a,±x)  ( 0 ≤ |x| ≤ 5 )
!       Output : W1F --- W(a,x)
!                W1D --- W'(a,x)
!                W2F --- W(a,-x)
!                W2D --- W'(a,-x)
!       Routine called:
!               CGAMA for computing complex gamma function
!       ======================================================
!
      IMPLICIT NONE
!*--PBWA9509
      DOUBLE PRECISION A , d , d1 , d2 , dl , eps , f1 , f2 , g1 , g2 , &
                     & h , h0 , h1 , hl , p0 , r , r1 , ugi , ugr , vgi
      DOUBLE PRECISION vgr , W1d , W1f , W2d , W2f , X , x1 , x2 , y1 , &
                     & y1d , y1f , y2d , y2f
      INTEGER k , l1 , l2 , m
      DIMENSION h(100) , d(80)
      eps = 1.0D-15
      p0 = 0.59460355750136D0
      IF ( A==0.0D0 ) THEN
         g1 = 3.625609908222D0
         g2 = 1.225416702465D0
      ELSE
         x1 = 0.25D0
         y1 = 0.5D0*A
         CALL CGAMA(x1,y1,1,ugr,ugi)
         g1 = DSQRT(ugr*ugr+ugi*ugi)
         x2 = 0.75D0
         CALL CGAMA(x2,y1,1,vgr,vgi)
         g2 = DSQRT(vgr*vgr+vgi*vgi)
      ENDIF
      f1 = DSQRT(g1/g2)
      f2 = DSQRT(2.0D0*g2/g1)
      h0 = 1.0D0
      h1 = A
      h(1) = A
      DO l1 = 4 , 200 , 2
         m = l1/2
         hl = A*h1 - 0.25D0*(l1-2.0D0)*(l1-3.0D0)*h0
         h(m) = hl
         h0 = h1
         h1 = hl
      ENDDO
      y1f = 1.0D0
      r = 1.0D0
      DO k = 1 , 100
         r = 0.5D0*r*X*X/(k*(2.0D0*k-1.0D0))
         r1 = h(k)*r
         y1f = y1f + r1
         IF ( DABS(r1)<=eps*DABS(y1f) .AND. k>30 ) GOTO 100
      ENDDO
 100  y1d = A
      r = 1.0D0
      DO k = 1 , 99
         r = 0.5D0*r*X*X/(k*(2.0D0*k+1.0D0))
         r1 = h(k+1)*r
         y1d = y1d + r1
         IF ( DABS(r1)<=eps*DABS(y1d) .AND. k>30 ) GOTO 200
      ENDDO
 200  y1d = X*y1d
      d1 = 1.0D0
      d2 = A
      d(1) = 1.0D0
      d(2) = A
      DO l2 = 5 , 160 , 2
         m = (l2+1)/2
         dl = A*d2 - 0.25D0*(l2-2.0D0)*(l2-3.0D0)*d1
         d(m) = dl
         d1 = d2
         d2 = dl
      ENDDO
      y2f = 1.0D0
      r = 1.0D0
      DO k = 1 , 79
         r = 0.5D0*r*X*X/(k*(2.0D0*k+1.0D0))
         r1 = d(k+1)*r
         y2f = y2f + r1
         IF ( DABS(r1)<=eps*DABS(y2f) .AND. k>30 ) GOTO 300
      ENDDO
 300  y2f = X*y2f
      y2d = 1.0D0
      r = 1.0D0
      DO k = 1 , 79
         r = 0.5D0*r*X*X/(k*(2.0D0*k-1.0D0))
         r1 = d(k+1)*r
         y2d = y2d + r1
         IF ( DABS(r1)<=eps*DABS(y2f) .AND. k>30 ) GOTO 400
      ENDDO
 400  W1f = p0*(f1*y1f-f2*y2f)
      W2f = p0*(f1*y1f+f2*y2f)
      W1d = p0*(f1*y1d-f2*y2d)
      W2d = p0*(f1*y1d+f2*y2d)
      END
 
 
 
!       **********************************
 
      SUBROUTINE RMN1(M,N,C,X,Df,Kd,R1f,R1d)
!
!       =======================================================
!       Purpose: Compute prolate and oblate spheroidal radial
!                functions of the first kind for given m, n,
!                c and x
!       Routines called:
!            (1) SCKB for computing expansion coefficients c2k
!            (2) SPHJ for computing the spherical Bessel
!                functions of the first kind
!       =======================================================
!
      IMPLICIT NONE
!*--RMN19613
      DOUBLE PRECISION a0 , b0 , C , ck , cx , Df , dj , eps , r , r0 , &
                     & r1 , R1d , R1f , r2 , r3 , reg , sa0 , sj , suc ,&
                     & sud
      DOUBLE PRECISION sum , sw , sw1 , X
      INTEGER ip , j , k , Kd , l , lg , M , N , nm , nm1 , nm2 , np
      DIMENSION ck(200) , Df(200) , sj(0:251) , dj(0:251)
      eps = 1.0D-14
      ip = 1
      nm1 = INT((N-M)/2)
      IF ( N-M==2*nm1 ) ip = 0
      nm = 25 + nm1 + INT(C)
      reg = 1.0D0
      IF ( M+nm>80 ) reg = 1.0D-200
      r0 = reg
      DO j = 1 , 2*M + ip
         r0 = r0*j
      ENDDO
      r = r0
      suc = r*Df(1)
      sw = 0.0D0
      DO k = 2 , nm
         r = r*(M+k-1.0)*(M+k+ip-1.5D0)/(k-1.0D0)/(k+ip-1.5D0)
         suc = suc + r*Df(k)
         IF ( k>nm1 .AND. DABS(suc-sw)<DABS(suc)*eps ) GOTO 100
         sw = suc
      ENDDO
 100  IF ( X==0.0 ) THEN
         CALL SCKB(M,N,C,Df,ck)
         sum = 0.0D0
         sw1 = 0.0D0
         DO j = 1 , nm
            sum = sum + ck(j)
            IF ( DABS(sum-sw1)<DABS(sum)*eps ) GOTO 150
            sw1 = sum
         ENDDO
 150     r1 = 1.0D0
         DO j = 1 , (N+M+ip)/2
            r1 = r1*(j+0.5D0*(N+M+ip))
         ENDDO
         r2 = 1.0D0
         DO j = 1 , M
            r2 = 2.0D0*C*r2*j
         ENDDO
         r3 = 1.0D0
         DO j = 1 , (N-M-ip)/2
            r3 = r3*j
         ENDDO
         sa0 = (2.0*(M+ip)+1.0)*r1/(2.0**N*C**ip*r2*r3)
         IF ( ip==0 ) THEN
            R1f = sum/(sa0*suc)*Df(1)*reg
            R1d = 0.0D0
         ELSEIF ( ip==1 ) THEN
            R1f = 0.0D0
            R1d = sum/(sa0*suc)*Df(1)*reg
         ENDIF
         RETURN
      ENDIF
      cx = C*X
      nm2 = 2*nm + M
      CALL SPHJ(nm2,cx,nm2,sj,dj)
      a0 = (1.0D0-Kd/(X*X))**(0.5D0*M)/suc
      R1f = 0.0D0
      sw = 0.0D0
      lg = 0
      DO k = 1 , nm
         l = 2*k + M - N - 2 + ip
         IF ( l==4*INT(l/4) ) lg = 1
         IF ( l/=4*INT(l/4) ) lg = -1
         IF ( k==1 ) THEN
            r = r0
         ELSE
            r = r*(M+k-1.0)*(M+k+ip-1.5D0)/(k-1.0D0)/(k+ip-1.5D0)
         ENDIF
         np = M + 2*k - 2 + ip
         R1f = R1f + lg*r*Df(k)*sj(np)
         IF ( k>nm1 .AND. DABS(R1f-sw)<DABS(R1f)*eps ) GOTO 200
         sw = R1f
      ENDDO
 200  R1f = R1f*a0
      b0 = Kd*M/X**3.0D0/(1.0-Kd/(X*X))*R1f
      sud = 0.0D0
      sw = 0.0D0
      DO k = 1 , nm
         l = 2*k + M - N - 2 + ip
         IF ( l==4*INT(l/4) ) lg = 1
         IF ( l/=4*INT(l/4) ) lg = -1
         IF ( k==1 ) THEN
            r = r0
         ELSE
            r = r*(M+k-1.0)*(M+k+ip-1.5D0)/(k-1.0D0)/(k+ip-1.5D0)
         ENDIF
         np = M + 2*k - 2 + ip
         sud = sud + lg*r*Df(k)*dj(np)
         IF ( k>nm1 .AND. DABS(sud-sw)<DABS(sud)*eps ) GOTO 300
         sw = sud
      ENDDO
 300  R1d = b0 + a0*C*sud
      END
 
 
 
!       **********************************
 
      SUBROUTINE DVSA(Va,X,Pd)
!
!       ===================================================
!       Purpose: Compute parabolic cylinder function Dv(x)
!                for small argument
!       Input:   x  --- Argument
!                va --- Order
!       Output:  PD --- Dv(x)
!       Routine called: GAMMA2 for computing Г(x)
!       ===================================================
!
      IMPLICIT NONE
!*--DVSA9732
      DOUBLE PRECISION a0 , ep , eps , g0 , g1 , ga0 , gm , Pd , pi ,   &
                     & r , r1 , sq2 , Va , va0 , vm , vt , X
      INTEGER m
      eps = 1.0D-15
      pi = 3.141592653589793D0
      sq2 = DSQRT(2.0D0)
      ep = EXP(-.25D0*X*X)
      va0 = 0.5D0*(1.0D0-Va)
      IF ( Va==0.0 ) THEN
         Pd = ep
      ELSEIF ( X==0.0 ) THEN
         IF ( va0<=0.0 .AND. va0==INT(va0) ) THEN
            Pd = 0.0D0
         ELSE
            CALL GAMMA2(va0,ga0)
            Pd = DSQRT(pi)/(2.0D0**(-.5D0*Va)*ga0)
         ENDIF
      ELSE
         CALL GAMMA2(-Va,g1)
         a0 = 2.0D0**(-0.5D0*Va-1.0D0)*ep/g1
         vt = -.5D0*Va
         CALL GAMMA2(vt,g0)
         Pd = g0
         r = 1.0D0
         DO m = 1 , 250
            vm = .5D0*(m-Va)
            CALL GAMMA2(vm,gm)
            r = -r*sq2*X/m
            r1 = gm*r
            Pd = Pd + r1
            IF ( DABS(r1)<DABS(Pd)*eps ) GOTO 50
         ENDDO
 50      Pd = a0*Pd
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE E1Z(Z,Ce1)
!
!       ====================================================
!       Purpose: Compute complex exponential integral E1(z)
!       Input :  z   --- Argument of E1(z)
!       Output:  CE1 --- E1(z)
!       ====================================================
!
      IMPLICIT NONE
!*--E1Z9785
      DOUBLE PRECISION a0 , el , pi , x , xt
      COMPLEX*16 Ce1 , cr , Z , zc , zd , zdc
      INTEGER k
      pi = 3.141592653589793D0
      el = 0.5772156649015328D0
      x = DBLE(Z)
      a0 = ABS(Z)
!       Continued fraction converges slowly near negative real axis,
!       so use power series in a wedge around it until radius 40.0
      xt = -2*DABS(DIMAG(Z))
      IF ( a0==0.0D0 ) THEN
         Ce1 = (1.0D+300,0.0D0)
      ELSEIF ( a0<=5.0 .OR. x<xt .AND. a0<40.0 ) THEN
!          Power series
         Ce1 = (1.0D0,0.0D0)
         cr = (1.0D0,0.0D0)
         DO k = 1 , 500
            cr = -cr*k*Z/(k+1.0D0)**2
            Ce1 = Ce1 + cr
            IF ( ABS(cr)<=ABS(Ce1)*1.0D-15 ) GOTO 50
         ENDDO
 50      IF ( x<=0.0 .AND. DIMAG(Z)==0.0 ) THEN
!     Careful on the branch cut -- use the sign of the imaginary part
!     to get the right sign on the factor if pi.
            Ce1 = -el - LOG(-Z) + Z*Ce1 - DSIGN(pi,DIMAG(Z))            &
                & *(0.0D0,1.0D0)
         ELSE
            Ce1 = -el - LOG(Z) + Z*Ce1
         ENDIF
      ELSE
!          Continued fraction https://dlmf.nist.gov/6.9
!
!                           1     1     1     2     2     3     3
!          E1 = exp(-z) * ----- ----- ----- ----- ----- ----- ----- ...
!                         Z +   1 +   Z +   1 +   Z +   1 +   Z +
         zc = 0D0
         zd = 1/Z
         zdc = 1*zd
         zc = zc + zdc
         DO k = 1 , 500
            zd = 1/(zd*k+1)
            zdc = (1*zd-1)*zdc
            zc = zc + zdc
 
            zd = 1/(zd*k+Z)
            zdc = (Z*zd-1)*zdc
            zc = zc + zdc
 
            IF ( ABS(zdc)<=ABS(zc)*1.0D-15 .AND. k>20 ) GOTO 100
         ENDDO
 100     Ce1 = EXP(-Z)*zc
         IF ( x<=0.0 .AND. DIMAG(Z)==0.0 ) Ce1 = Ce1 - pi*(0.0D0,1.0D0)
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE ITJYB(X,Tj,Ty)
!
!       =======================================================
!       Purpose: Integrate Bessel functions J0(t) and Y0(t)
!                with respect to t from 0 to x ( x ≥ 0 )
!       Input :  x  --- Upper limit of the integral
!       Output:  TJ --- Integration of J0(t) from 0 to x
!                TY --- Integration of Y0(t) from 0 to x
!       =======================================================
!
      IMPLICIT NONE
!*--ITJYB9857
      DOUBLE PRECISION f0 , g0 , pi , t , Tj , Ty , X , x1 , xt
      pi = 3.141592653589793D0
      IF ( X==0.0D0 ) THEN
         Tj = 0.0D0
         Ty = 0.0D0
      ELSEIF ( X<=4.0D0 ) THEN
         x1 = X/4.0D0
         t = x1*x1
         Tj = (((((((-.133718D-3*t+.2362211D-2)*t-.025791036D0)*t+      &
            & .197492634D0)*t-1.015860606D0)*t+3.199997842D0)           &
            & *t-5.333333161D0)*t+4.0D0)*x1
         Ty = ((((((((.13351D-4*t-.235002D-3)*t+.3034322D-2)*t-         &
            & .029600855D0)*t+.203380298D0)*t-.904755062D0)             &
            & *t+2.287317974D0)*t-2.567250468D0)*t+1.076611469D0)*x1
         Ty = 2.0D0/pi*DLOG(X/2.0D0)*Tj - Ty
      ELSEIF ( X<=8.0D0 ) THEN
         xt = X - .25D0*pi
         t = 16.0D0/(X*X)
         f0 = ((((((.1496119D-2*t-.739083D-2)*t+.016236617D0)*t-        &
            & .022007499D0)*t+.023644978D0)*t-.031280848D0)             &
            & *t+.124611058D0)*4.0D0/X
         g0 = (((((.1076103D-2*t-.5434851D-2)*t+.01242264D0)*t-         &
            & .018255209)*t+.023664841D0)*t-.049635633D0)               &
            & *t + .79784879D0
         Tj = 1.0D0 - (f0*DCOS(xt)-g0*DSIN(xt))/DSQRT(X)
         Ty = -(f0*DSIN(xt)+g0*DCOS(xt))/DSQRT(X)
      ELSE
         t = 64.0D0/(X*X)
         xt = X - .25D0*pi
         f0 = (((((((-.268482D-4*t+.1270039D-3)*t-.2755037D-3)*t+       &
            & .3992825D-3)*t-.5366169D-3)*t+.10089872D-2)               &
            & *t-.40403539D-2)*t+.0623347304D0)*8.0D0/X
         g0 = ((((((-.226238D-4*t+.1107299D-3)*t-.2543955D-3)*t+        &
            & .4100676D-3)*t-.6740148D-3)*t+.17870944D-2)               &
            & *t-.01256424405D0)*t + .79788456D0
         Tj = 1.0D0 - (f0*DCOS(xt)-g0*DSIN(xt))/DSQRT(X)
         Ty = -(f0*DSIN(xt)+g0*DCOS(xt))/DSQRT(X)
      ENDIF
      END
 
 
!       **********************************
 
      SUBROUTINE CHGUL(A,B,X,Hu,Id)
!
!       =======================================================
!       Purpose: Compute the confluent hypergeometric function
!                U(a,b,x) for large argument x
!       Input  : a  --- Parameter
!                b  --- Parameter
!                x  --- Argument
!       Output:  HU --- U(a,b,x)
!                ID --- Estimated number of significant digits
!       =======================================================
!
      IMPLICIT NONE
!*--CHGUL9917
      DOUBLE PRECISION A , aa , B , Hu , r , r0 , ra , X
      INTEGER Id , k , nm
      LOGICAL il1 , il2
      Id = -100
      aa = A - B + 1.0D0
      il1 = A==INT(A) .AND. A<=0.0
      il2 = aa==INT(aa) .AND. aa<=0.0
      nm = 0
      IF ( il1 ) nm = ABS(A)
      IF ( il2 ) nm = ABS(aa)
!       IL1: DLMF 13.2.7 with k=-s-a
!       IL2: DLMF 13.2.8
      IF ( il1 .OR. il2 ) THEN
         Hu = 1.0D0
         r = 1.0D0
         DO k = 1 , nm
            r = -r*(A+k-1.0D0)*(A-B+k)/(k*X)
            Hu = Hu + r
         ENDDO
         Hu = X**(-A)*Hu
         Id = 10
      ELSE
!       DLMF 13.7.3
         Hu = 1.0D0
         r = 1.0D0
         DO k = 1 , 25
            r = -r*(A+k-1.0D0)*(A-B+k)/(k*X)
            ra = DABS(r)
            IF ( k>5 .AND. ra>=r0 .OR. ra<1.0D-15 ) GOTO 50
            r0 = ra
            Hu = Hu + r
         ENDDO
 50      Id = ABS(LOG10(ra))
         Hu = X**(-A)*Hu
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE GMN(M,N,C,X,Bk,Gf,Gd)
!
!       ===========================================================
!       Purpose: Compute gmn(-ic,ix) and its derivative for oblate
!                radial functions with a small argument
!       ===========================================================
!
      IMPLICIT NONE
!*--GMN9970
      DOUBLE PRECISION Bk , C , eps , Gd , gd0 , gd1 , Gf , gf0 , gw ,  &
                     & X , xm
      INTEGER ip , k , M , N , nm
      DIMENSION Bk(200)
      eps = 1.0D-14
      ip = 1
      IF ( N-M==2*INT((N-M)/2) ) ip = 0
      nm = 25 + INT(0.5*(N-M)+C)
      xm = (1.0D0+X*X)**(-0.5D0*M)
      gf0 = 0.0D0
      gw = 0.0D0
      DO k = 1 , nm
         gf0 = gf0 + Bk(k)*X**(2.0*k-2.0)
         IF ( DABS((gf0-gw)/gf0)<eps .AND. k>=10 ) GOTO 100
         gw = gf0
      ENDDO
 100  Gf = xm*gf0*X**(1-ip)
      gd1 = -M*X/(1.0D0+X*X)*Gf
      gd0 = 0.0D0
      DO k = 1 , nm
         IF ( ip==0 ) THEN
            gd0 = gd0 + (2.0D0*k-1.0)*Bk(k)*X**(2.0*k-2.0)
         ELSE
            gd0 = gd0 + 2.0D0*k*Bk(k+1)*X**(2.0*k-1.0)
         ENDIF
         IF ( DABS((gd0-gw)/gd0)<eps .AND. k>=10 ) GOTO 200
         gw = gd0
      ENDDO
 200  Gd = gd1 + xm*gd0
      END
 
 
 
!       **********************************
 
      SUBROUTINE ITJYA(X,Tj,Ty)
!
!       ==========================================================
!       Purpose: Integrate Bessel functions J0(t) & Y0(t) with
!                respect to t from 0 to x
!       Input :  x  --- Upper limit of the integral ( x >= 0 )
!       Output:  TJ --- Integration of J0(t) from 0 to x
!                TY --- Integration of Y0(t) from 0 to x
!       =======================================================
!
      IMPLICIT NONE
!*--ITJYA10020
      DOUBLE PRECISION a , a0 , a1 , af , bf , bg , el , eps , pi , r , &
                     & r2 , rc , rs , Tj , Ty , ty1 , ty2 , X , x2 , xp
      INTEGER k
      DIMENSION a(18)
      pi = 3.141592653589793D0
      el = .5772156649015329D0
      eps = 1.0D-12
      IF ( X==0.0D0 ) THEN
         Tj = 0.0D0
         Ty = 0.0D0
      ELSEIF ( X<=20.0D0 ) THEN
         x2 = X*X
         Tj = X
         r = X
         DO k = 1 , 60
            r = -.25D0*r*(2*k-1.0D0)/(2*k+1.0D0)/(k*k)*x2
            Tj = Tj + r
            IF ( DABS(r)<DABS(Tj)*eps ) GOTO 50
         ENDDO
 50      ty1 = (el+DLOG(X/2.0D0))*Tj
         rs = 0.0D0
         ty2 = 1.0D0
         r = 1.0D0
         DO k = 1 , 60
            r = -.25D0*r*(2*k-1.0D0)/(2*k+1.0D0)/(k*k)*x2
            rs = rs + 1.0D0/k
            r2 = r*(rs+1.0D0/(2.0D0*k+1.0D0))
            ty2 = ty2 + r2
            IF ( DABS(r2)<DABS(ty2)*eps ) GOTO 100
         ENDDO
 100     Ty = (ty1-X*ty2)*2.0D0/pi
      ELSE
         a0 = 1.0D0
         a1 = 5.0D0/8.0D0
         a(1) = a1
         DO k = 1 , 16
            af = ((1.5D0*(k+.5D0)*(k+5.0D0/6.0D0)*a1-.5D0*(k+.5D0)      &
               & *(k+.5D0)*(k-.5D0)*a0))/(k+1.0D0)
            a(k+1) = af
            a0 = a1
            a1 = af
         ENDDO
         bf = 1.0D0
         r = 1.0D0
         DO k = 1 , 8
            r = -r/(X*X)
            bf = bf + a(2*k)*r
         ENDDO
         bg = a(1)/X
         r = 1.0D0/X
         DO k = 1 , 8
            r = -r/(X*X)
            bg = bg + a(2*k+1)*r
         ENDDO
         xp = X + .25D0*pi
         rc = DSQRT(2.0D0/(pi*X))
         Tj = 1.0D0 - rc*(bf*DCOS(xp)+bg*DSIN(xp))
         Ty = rc*(bg*DCOS(xp)-bf*DSIN(xp))
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE RCTY(N,X,Nm,Ry,Dy)
!
!       ========================================================
!       Purpose: Compute Riccati-Bessel functions of the second
!                kind and their derivatives
!       Input:   x --- Argument of Riccati-Bessel function
!                n --- Order of yn(x)
!       Output:  RY(n) --- x·yn(x)
!                DY(n) --- [x·yn(x)]'
!                NM --- Highest order computed
!       ========================================================
!
      IMPLICIT NONE
!*--RCTY10100
      DOUBLE PRECISION Dy , rf0 , rf1 , rf2 , Ry , X
      INTEGER k , N , Nm
      DIMENSION Ry(0:N) , Dy(0:N)
      Nm = N
      IF ( X<1.0D-60 ) THEN
         DO k = 0 , N
            Ry(k) = -1.0D+300
            Dy(k) = 1.0D+300
         ENDDO
         Ry(0) = -1.0D0
         Dy(0) = 0.0D0
         RETURN
      ENDIF
      Ry(0) = -DCOS(X)
      Ry(1) = Ry(0)/X - DSIN(X)
      rf0 = Ry(0)
      rf1 = Ry(1)
      DO k = 2 , N
         rf2 = (2.0D0*k-1.0D0)*rf1/X - rf0
         IF ( DABS(rf2)>1.0D+300 ) GOTO 100
         Ry(k) = rf2
         rf0 = rf1
         rf1 = rf2
      ENDDO
 100  Nm = k - 1
      Dy(0) = DSIN(X)
      DO k = 1 , Nm
         Dy(k) = -k*Ry(k)/X + Ry(k-1)
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE LPNI(N,X,Pn,Pd,Pl)
!
!       =====================================================
!       Purpose: Compute Legendre polynomials Pn(x), Pn'(x)
!                and the integral of Pn(t) from 0 to x
!       Input :  x --- Argument of Pn(x)
!                n --- Degree of Pn(x) ( n = 0,1,... )
!       Output:  PN(n) --- Pn(x)
!                PD(n) --- Pn'(x)
!                PL(n) --- Integral of Pn(t) from 0 to x
!       =====================================================
!
      IMPLICIT NONE
!*--LPNI10150
      INTEGER j , k , N , n1
      DOUBLE PRECISION p0 , p1 , Pd , pf , Pl , Pn , r , X
      DIMENSION Pn(0:N) , Pd(0:N) , Pl(0:N)
      Pn(0) = 1.0D0
      Pn(1) = X
      Pd(0) = 0.0D0
      Pd(1) = 1.0D0
      Pl(0) = X
      Pl(1) = 0.5D0*X*X
      p0 = 1.0D0
      p1 = X
      DO k = 2 , N
         pf = (2.0D0*k-1.0D0)/k*X*p1 - (k-1.0D0)/k*p0
         Pn(k) = pf
         IF ( DABS(X)==1.0D0 ) THEN
            Pd(k) = 0.5D0*X**(k+1)*k*(k+1.0D0)
         ELSE
            Pd(k) = k*(p1-X*pf)/(1.0D0-X*X)
         ENDIF
         Pl(k) = (X*Pn(k)-Pn(k-1))/(k+1.0D0)
         p0 = p1
         p1 = pf
         IF ( k/=2*INT(k/2) ) THEN
            r = 1.0D0/(k+1.0D0)
            n1 = (k-1)/2
            DO j = 1 , n1
               r = (0.5D0/j-1.0D0)*r
            ENDDO
            Pl(k) = Pl(k) + r
         ENDIF
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE KLVNA(X,Ber,Bei,Ger,Gei,Der,Dei,Her,Hei)
!
!       ======================================================
!       Purpose: Compute Kelvin functions ber x, bei x, ker x
!                and kei x, and their derivatives  ( x > 0 )
!       Input :  x   --- Argument of Kelvin functions
!       Output:  BER --- ber x
!                BEI --- bei x
!                GER --- ker x
!                GEI --- kei x
!                DER --- ber'x
!                DEI --- bei'x
!                HER --- ker'x
!                HEI --- kei'x
!       ================================================
!
      IMPLICIT NONE
!*--KLVNA10206
      DOUBLE PRECISION Bei , Ber , cn0 , cp0 , cs , Dei , Der , el ,    &
                     & eps , fac , Gei , Ger , gs , Hei , Her , pi ,    &
                     & pn0 , pn1 , pp0 , pp1
      DOUBLE PRECISION qn0 , qn1 , qp0 , qp1 , r , r0 , r1 , rc , rs ,  &
                     & sn0 , sp0 , ss , X , x2 , x4 , xc1 , xc2 , xd ,  &
                     & xe1 , xe2
      DOUBLE PRECISION xt
      INTEGER k , km , m
      pi = 3.141592653589793D0
      el = .5772156649015329D0
      eps = 1.0D-15
      IF ( X==0.0D0 ) THEN
         Ber = 1.0D0
         Bei = 0.0D0
         Ger = 1.0D+300
         Gei = -0.25D0*pi
         Der = 0.0D0
         Dei = 0.0D0
         Her = -1.0D+300
         Hei = 0.0D0
         RETURN
      ENDIF
      x2 = 0.25D0*X*X
      x4 = x2*x2
      IF ( DABS(X)<10.0D0 ) THEN
         Ber = 1.0D0
         r = 1.0D0
         DO m = 1 , 60
            r = -0.25D0*r/(m*m)/(2.0D0*m-1.0D0)**2*x4
            Ber = Ber + r
            IF ( DABS(r)<DABS(Ber)*eps ) GOTO 50
         ENDDO
 50      Bei = x2
         r = x2
         DO m = 1 , 60
            r = -0.25D0*r/(m*m)/(2.0D0*m+1.0D0)**2*x4
            Bei = Bei + r
            IF ( DABS(r)<DABS(Bei)*eps ) GOTO 100
         ENDDO
 100     Ger = -(DLOG(X/2.0D0)+el)*Ber + 0.25D0*pi*Bei
         r = 1.0D0
         gs = 0.0D0
         DO m = 1 , 60
            r = -0.25D0*r/(m*m)/(2.0D0*m-1.0D0)**2*x4
            gs = gs + 1.0D0/(2.0D0*m-1.0D0) + 1.0D0/(2.0D0*m)
            Ger = Ger + r*gs
            IF ( DABS(r*gs)<DABS(Ger)*eps ) GOTO 150
         ENDDO
 150     Gei = x2 - (DLOG(X/2.0D0)+el)*Bei - 0.25D0*pi*Ber
         r = x2
         gs = 1.0D0
         DO m = 1 , 60
            r = -0.25D0*r/(m*m)/(2.0D0*m+1.0D0)**2*x4
            gs = gs + 1.0D0/(2.0D0*m) + 1.0D0/(2.0D0*m+1.0D0)
            Gei = Gei + r*gs
            IF ( DABS(r*gs)<DABS(Gei)*eps ) GOTO 200
         ENDDO
 200     Der = -0.25D0*X*x2
         r = Der
         DO m = 1 , 60
            r = -0.25D0*r/m/(m+1.0D0)/(2.0D0*m+1.0D0)**2*x4
            Der = Der + r
            IF ( DABS(r)<DABS(Der)*eps ) GOTO 250
         ENDDO
 250     Dei = 0.5D0*X
         r = Dei
         DO m = 1 , 60
            r = -0.25D0*r/(m*m)/(2.D0*m-1.D0)/(2.D0*m+1.D0)*x4
            Dei = Dei + r
            IF ( DABS(r)<DABS(Dei)*eps ) GOTO 300
         ENDDO
 300     r = -0.25D0*X*x2
         gs = 1.5D0
         Her = 1.5D0*r - Ber/X - (DLOG(X/2.D0)+el)*Der + 0.25*pi*Dei
         DO m = 1 , 60
            r = -0.25D0*r/m/(m+1.0D0)/(2.0D0*m+1.0D0)**2*x4
            gs = gs + 1.0D0/(2*m+1.0D0) + 1.0D0/(2*m+2.0D0)
            Her = Her + r*gs
            IF ( DABS(r*gs)<DABS(Her)*eps ) GOTO 350
         ENDDO
 350     r = 0.5D0*X
         gs = 1.0D0
         Hei = 0.5D0*X - Bei/X - (DLOG(X/2.D0)+el)*Dei - 0.25*pi*Der
         DO m = 1 , 60
            r = -0.25D0*r/(m*m)/(2*m-1.0D0)/(2*m+1.0D0)*x4
            gs = gs + 1.0D0/(2.0D0*m) + 1.0D0/(2*m+1.0D0)
            Hei = Hei + r*gs
            IF ( DABS(r*gs)<DABS(Hei)*eps ) RETURN
         ENDDO
      ELSE
         pp0 = 1.0D0
         pn0 = 1.0D0
         qp0 = 0.0D0
         qn0 = 0.0D0
         r0 = 1.0D0
         km = 18
         IF ( DABS(X)>=40.0 ) km = 10
         fac = 1.0D0
         DO k = 1 , km
            fac = -fac
            xt = 0.25D0*k*pi - INT(0.125D0*k)*2.0D0*pi
            cs = COS(xt)
            ss = SIN(xt)
            r0 = 0.125D0*r0*(2.0D0*k-1.0D0)**2/k/X
            rc = r0*cs
            rs = r0*ss
            pp0 = pp0 + rc
            pn0 = pn0 + fac*rc
            qp0 = qp0 + rs
            qn0 = qn0 + fac*rs
         ENDDO
         xd = X/DSQRT(2.0D0)
         xe1 = EXP(xd)
         xe2 = EXP(-xd)
         xc1 = 1.D0/DSQRT(2.0D0*pi*X)
         xc2 = DSQRT(.5D0*pi/X)
         cp0 = DCOS(xd+0.125D0*pi)
         cn0 = DCOS(xd-0.125D0*pi)
         sp0 = DSIN(xd+0.125D0*pi)
         sn0 = DSIN(xd-0.125D0*pi)
         Ger = xc2*xe2*(pn0*cp0-qn0*sp0)
         Gei = xc2*xe2*(-pn0*sp0-qn0*cp0)
         Ber = xc1*xe1*(pp0*cn0+qp0*sn0) - Gei/pi
         Bei = xc1*xe1*(pp0*sn0-qp0*cn0) + Ger/pi
         pp1 = 1.0D0
         pn1 = 1.0D0
         qp1 = 0.0D0
         qn1 = 0.0D0
         r1 = 1.0D0
         fac = 1.0D0
         DO k = 1 , km
            fac = -fac
            xt = 0.25D0*k*pi - INT(0.125D0*k)*2.0D0*pi
            cs = DCOS(xt)
            ss = DSIN(xt)
            r1 = 0.125D0*r1*(4.D0-(2.0D0*k-1.0D0)**2)/k/X
            rc = r1*cs
            rs = r1*ss
            pp1 = pp1 + fac*rc
            pn1 = pn1 + rc
            qp1 = qp1 + fac*rs
            qn1 = qn1 + rs
         ENDDO
         Her = xc2*xe2*(-pn1*cn0+qn1*sn0)
         Hei = xc2*xe2*(pn1*sn0+qn1*cn0)
         Der = xc1*xe1*(pp1*cp0+qp1*sp0) - Hei/pi
         Dei = xc1*xe1*(pp1*sp0-qp1*cp0) + Her/pi
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE CHGUBI(A,B,X,Hu,Id)
!
!       ======================================================
!       Purpose: Compute confluent hypergeometric function
!                U(a,b,x) with integer b ( b = ±1,±2,... )
!       Input  : a  --- Parameter
!                b  --- Parameter
!                x  --- Argument
!       Output:  HU --- U(a,b,x)
!                ID --- Estimated number of significant digits
!       Routines called:
!            (1) GAMMA2 for computing gamma function Г(x)
!            (2) PSI_SPEC for computing psi function
!       ======================================================
!
      IMPLICIT NONE
!*--CHGUBI10378
      DOUBLE PRECISION A , a0 , a1 , a2 , B , da1 , da2 , db1 , db2 ,   &
                     & el , ga , ga1 , h0 , hm1 , hm2 , hm3 , hmax ,    &
                     & hmin , Hu , hu1
      DOUBLE PRECISION hu2 , hw , ps , r , rn , rn1 , s0 , s1 , s2 ,    &
                     & sa , sb , ua , ub , X
      INTEGER Id , id1 , id2 , j , k , m , n
      Id = -100
      el = 0.5772156649015329D0
      n = ABS(B-1)
      rn1 = 1.0D0
      rn = 1.0D0
      DO j = 1 , n
         rn = rn*j
         IF ( j==n-1 ) rn1 = rn
      ENDDO
      CALL PSI_SPEC(A,ps)
      CALL GAMMA2(A,ga)
      IF ( B>0.0 ) THEN
         a0 = A
         a1 = A - n
         a2 = a1
         CALL GAMMA2(a1,ga1)
         ua = (-1)**(n-1)/(rn*ga1)
         ub = rn1/ga*X**(-n)
      ELSE
         a0 = A + n
         a1 = a0
         a2 = A
         CALL GAMMA2(a1,ga1)
         ua = (-1)**(n-1)/(rn*ga)*X**n
         ub = rn1/ga1
      ENDIF
      hm1 = 1.0D0
      r = 1.0D0
      hmax = 0.0D0
      hmin = 1.0D+300
      h0 = 0D0
      DO k = 1 , 150
         r = r*(a0+k-1.0D0)*X/((n+k)*k)
         hm1 = hm1 + r
         hu1 = DABS(hm1)
         IF ( hu1>hmax ) hmax = hu1
         IF ( hu1<hmin ) hmin = hu1
         IF ( DABS(hm1-h0)<DABS(hm1)*1.0D-15 ) GOTO 100
         h0 = hm1
      ENDDO
 100  da1 = LOG10(hmax)
      da2 = 0.0D0
      IF ( hmin/=0.0 ) da2 = LOG10(hmin)
      Id = 15 - ABS(da1-da2)
      hm1 = hm1*DLOG(X)
      s0 = 0.0D0
      DO m = 1 , n
         IF ( B>=0.0 ) s0 = s0 - 1.0D0/m
         IF ( B<0.0 ) s0 = s0 + (1.0D0-A)/(m*(A+m-1.0D0))
      ENDDO
      hm2 = ps + 2.0D0*el + s0
      r = 1.0D0
      hmax = 0.0D0
      hmin = 1.0D+300
      DO k = 1 , 150
         s1 = 0.0D0
         s2 = 0.0D0
         IF ( B>0.0 ) THEN
            DO m = 1 , k
               s1 = s1 - (m+2.0D0*A-2.0D0)/(m*(m+A-1.0D0))
            ENDDO
            DO m = 1 , n
               s2 = s2 + 1.0D0/(k+m)
            ENDDO
         ELSE
            DO m = 1 , k + n
               s1 = s1 + (1.0D0-A)/(m*(m+A-1.0D0))
            ENDDO
            DO m = 1 , k
               s2 = s2 + 1.0D0/m
            ENDDO
         ENDIF
         hw = 2.0D0*el + ps + s1 - s2
         r = r*(a0+k-1.0D0)*X/((n+k)*k)
         hm2 = hm2 + r*hw
         hu2 = DABS(hm2)
         IF ( hu2>hmax ) hmax = hu2
         IF ( hu2<hmin ) hmin = hu2
         IF ( DABS((hm2-h0)/hm2)<1.0D-15 ) GOTO 200
         h0 = hm2
      ENDDO
 200  db1 = LOG10(hmax)
      db2 = 0.0D0
      IF ( hmin/=0.0 ) db2 = LOG10(hmin)
      id1 = 15 - ABS(db1-db2)
      IF ( id1<Id ) Id = id1
      hm3 = 1.0D0
      IF ( n==0 ) hm3 = 0.0D0
      r = 1.0D0
      DO k = 1 , n - 1
         r = r*(a2+k-1.0D0)/((k-n)*k)*X
         hm3 = hm3 + r
      ENDDO
      sa = ua*(hm1+hm2)
      sb = ub*hm3
      Hu = sa + sb
      id2 = 0.0D0
      IF ( sa/=0.0 ) id1 = INT(LOG10(ABS(sa)))
      IF ( Hu/=0.0 ) id2 = INT(LOG10(ABS(Hu)))
      IF ( sa*sb<0.0 ) Id = Id - ABS(id1-id2)
      END
 
 
 
!       **********************************
 
      SUBROUTINE CYZO(Nt,Kf,Kc,Zo,Zv)
!
!       ===========================================================
!       Purpose : Compute the complex zeros of Y0(z), Y1(z) and
!                 Y1'(z), and their associated values at the zeros
!                 using the modified Newton's iteration method
!       Input:    NT --- Total number of zeros/roots
!                 KF --- Function choice code
!                        KF=0 for  Y0(z) & Y1(z0)
!                        KF=1 for  Y1(z) & Y0(z1)
!                        KF=2 for  Y1'(z) & Y1(z1')
!                 KC --- Choice code
!                        KC=0 for complex roots
!                        KC=1 for real roots
!       Output:   ZO(L) --- L-th zero of Y0(z) or Y1(z) or Y1'(z)
!                 ZV(L) --- Value of Y0'(z) or Y1'(z) or Y1(z)
!                           at the L-th zero
!       Routine called: CY01 for computing Y0(z) and Y1(z), and
!                       their derivatives
!       ===========================================================
      IMPLICIT NONE
!*--CYZO10515
      DOUBLE PRECISION h , w , w0 , x , y
      INTEGER i , it , j , Kc , Kf , nr , Nt
      COMPLEX*16 z , zd , zero , zf , zfd , zgd , Zo , zp , zq , Zv , zw
      DIMENSION Zo(Nt) , Zv(Nt)
      x = 0.0D0
      y = 0.0D0
      h = 0.0D0
      IF ( Kc==0 ) THEN
         x = -2.4D0
         y = 0.54D0
         h = 3.14D0
      ELSEIF ( Kc==1 ) THEN
         x = 0.89
         y = 0.0
         h = -3.14
      ENDIF
      IF ( Kf==1 ) x = -0.503
      IF ( Kf==2 ) x = 0.577
      zero = DCMPLX(x,y)
      z = zero
      w = 0.0D0
      DO nr = 1 , Nt
         IF ( nr/=1 ) z = Zo(nr-1) - h
         it = 0
 50      it = it + 1
         CALL CY01(Kf,z,zf,zd)
         zp = (1.0D0,0.0D0)
         DO i = 1 , nr - 1
            zp = zp*(z-Zo(i))
         ENDDO
         zfd = zf/zp
         zq = (0.0D0,0.0D0)
         DO i = 1 , nr - 1
            zw = (1.0D0,0.0D0)
            DO j = 1 , nr - 1
               IF ( j/=i ) zw = zw*(z-Zo(j))
            ENDDO
            zq = zq + zw
         ENDDO
         zgd = (zd-zq*zfd)/zp
         z = z - zfd/zgd
         w0 = w
         w = ABS(z)
         IF ( it<=50 .AND. DABS((w-w0)/w)>1.0D-12 ) GOTO 50
         Zo(nr) = z
      ENDDO
      DO i = 1 , Nt
         z = Zo(i)
         IF ( Kf==0 .OR. Kf==2 ) THEN
            CALL CY01(1,z,zf,zd)
            Zv(i) = zf
         ELSEIF ( Kf==1 ) THEN
            CALL CY01(0,z,zf,zd)
            Zv(i) = zf
         ENDIF
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE KLVNB(X,Ber,Bei,Ger,Gei,Der,Dei,Her,Hei)
!
!       ======================================================
!       Purpose: Compute Kelvin functions ber x, bei x, ker x
!                and kei x, and their derivatives  ( x > 0 )
!       Input :  x   --- Argument of Kelvin functions
!       Output:  BER --- ber x
!                BEI --- bei x
!                GER --- ker x
!                GEI --- kei x
!                DER --- ber'x
!                DEI --- bei'x
!                HER --- ker'x
!                HEI --- kei'x
!       ================================================
!
      IMPLICIT NONE
!*--KLVNB10598
      DOUBLE PRECISION Bei , Ber , csn , csp , Dei , Der , fxi , fxr ,  &
                     & Gei , Ger , Hei , Her , pi , pni , pnr , ppi ,   &
                     & ppr , ssn , ssp , t
      DOUBLE PRECISION t2 , tni , tnr , tpi , tpr , u , v , X , yc1 ,   &
                     & yc2 , yd , ye1 , ye2
      INTEGER l
      pi = 3.141592653589793D0
      IF ( X==0.0D0 ) THEN
         Ber = 1.0D0
         Bei = 0.0D0
         Ger = 1.0D+300
         Gei = -.25D0*pi
         Der = 0.0D0
         Dei = 0.0D0
         Her = -1.0D+300
         Hei = 0.0D0
      ELSEIF ( X<8.0D0 ) THEN
         t = X/8.0D0
         t2 = t*t
         u = t2*t2
         Ber = ((((((-.901D-5*u+.122552D-2)*u-.08349609D0)*u+           &
             & 2.64191397D0)*u-32.36345652D0)*u+113.77777774D0)         &
             & *u-64.0D0)*u + 1.0D0
         Bei = t*t*((((((.11346D-3*u-.01103667D0)*u+.52185615D0)*u-     &
             & 10.56765779D0)*u+72.81777742D0)*u-113.77777774D0)        &
             & *u+16.0D0)
         Ger = ((((((-.2458D-4*u+.309699D-2)*u-.19636347D0)*u+          &
             & 5.65539121D0)*u-60.60977451D0)*u+171.36272133D0)         &
             & *u-59.05819744D0)*u - .57721566D0
         Ger = Ger - DLOG(.5D0*X)*Ber + .25D0*pi*Bei
         Gei = t2*((((((.29532D-3*u-.02695875D0)*u+1.17509064D0)*u-     &
             & 21.30060904D0)*u+124.2356965D0)*u-142.91827687D0)        &
             & *u+6.76454936D0)
         Gei = Gei - DLOG(.5D0*X)*Bei - .25D0*pi*Ber
         Der = X*t2*                                                    &
             & ((((((-.394D-5*u+.45957D-3)*u-.02609253D0)*u+.66047849D0)&
             & *u-6.0681481D0)*u+14.22222222D0)*u-4.0D0)
         Dei = X*((((((.4609D-4*u-.379386D-2)*u+.14677204D0)*u-         &
             & 2.31167514D0)*u+11.37777772D0)*u-10.66666666D0)*u+.5D0)
         Her = X*t2*((((((-.1075D-4*u+.116137D-2)*u-.06136358D0)*u+     &
             & 1.4138478D0)*u-11.36433272D0)*u+21.42034017D0)           &
             & *u-3.69113734D0)
         Her = Her - DLOG(.5D0*X)*Der - Ber/X + .25D0*pi*Dei
         Hei = X*((((((.11997D-3*u-.926707D-2)*u+.33049424D0)*u-        &
             & 4.65950823D0)*u+19.41182758D0)*u-13.39858846D0)          &
             & *u+.21139217D0)
         Hei = Hei - DLOG(.5D0*X)*Dei - Bei/X - .25D0*pi*Der
      ELSE
         t = 8.0D0/X
         tnr = 0.0D0
         tni = 0.0D0
         DO l = 1 , 2
            v = (-1)**l*t
            tpr = ((((.6D-6*v-.34D-5)*v-.252D-4)*v-.906D-4)             &
                & *v*v+.0110486D0)*v
            tpi = ((((.19D-5*v+.51D-5)*v*v-.901D-4)*v-.9765D-3)         &
                & *v-.0110485D0)*v - .3926991D0
            IF ( l==1 ) THEN
               tnr = tpr
               tni = tpi
            ENDIF
         ENDDO
         yd = X/DSQRT(2.0D0)
         ye1 = EXP(yd+tpr)
         ye2 = EXP(-yd+tnr)
         yc1 = 1.0D0/DSQRT(2.0D0*pi*X)
         yc2 = DSQRT(pi/(2.0D0*X))
         csp = DCOS(yd+tpi)
         ssp = DSIN(yd+tpi)
         csn = DCOS(-yd+tni)
         ssn = DSIN(-yd+tni)
         Ger = yc2*ye2*csn
         Gei = yc2*ye2*ssn
         fxr = yc1*ye1*csp
         fxi = yc1*ye1*ssp
         Ber = fxr - Gei/pi
         Bei = fxi + Ger/pi
         pnr = 0.0D0
         pni = 0.0D0
         DO l = 1 , 2
            v = (-1)**l*t
            ppr = (((((.16D-5*v+.117D-4)*v+.346D-4)*v+.5D-6)*v-.13813D-2&
                & )*v-.0625001D0)*v + .7071068D0
            ppi = (((((-.32D-5*v-.24D-5)*v+.338D-4)*v+.2452D-3)*v+      &
                & .13811D-2)*v-.1D-6)*v + .7071068D0
            IF ( l==1 ) THEN
               pnr = ppr
               pni = ppi
            ENDIF
         ENDDO
         Her = Gei*pni - Ger*pnr
         Hei = -(Gei*pnr+Ger*pni)
         Der = fxr*ppr - fxi*ppi - Hei/pi
         Dei = fxi*ppr + fxr*ppi + Her/pi
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE RMN2SO(M,N,C,X,Cv,Df,Kd,R2f,R2d)
!
!       =============================================================
!       Purpose: Compute oblate radial functions of the second kind
!                with a small argument, Rmn(-ic,ix) & Rmn'(-ic,ix)
!       Routines called:
!            (1) SCKB for computing the expansion coefficients c2k
!            (2) KMN for computing the joining factors
!            (3) QSTAR for computing the factor defined in (15.7.3)
!            (4) CBK for computing the the expansion coefficient
!                defined in (15.7.6)
!            (5) GMN for computing the function defined in (15.7.4)
!            (6) RMN1 for computing the radial function of the first
!                kind
!       =============================================================
!
      IMPLICIT NONE
!*--RMN2SO10718
      DOUBLE PRECISION bk , C , ck , ck1 , ck2 , Cv , Df , dn , eps ,   &
                     & gd , gf , h0 , pi , qs , qt , r1d , r1f , R2d ,  &
                     & R2f , sum
      DOUBLE PRECISION sw , X
      INTEGER ip , j , Kd , M , N , nm
      DIMENSION bk(200) , ck(200) , Df(200) , dn(200)
      IF ( DABS(Df(1))<=1.0D-280 ) THEN
         R2f = 1.0D+300
         R2d = 1.0D+300
         RETURN
      ENDIF
      eps = 1.0D-14
      pi = 3.141592653589793D0
      nm = 25 + INT((N-M)/2+C)
      ip = 1
      IF ( N-M==2*INT((N-M)/2) ) ip = 0
      CALL SCKB(M,N,C,Df,ck)
      CALL KMN(M,N,C,Cv,Kd,Df,dn,ck1,ck2)
      CALL QSTAR(M,N,C,ck,ck1,qs,qt)
      CALL CBK(M,N,C,Cv,qt,ck,bk)
      IF ( X==0.0D0 ) THEN
         sum = 0.0D0
         sw = 0.0D0
         DO j = 1 , nm
            sum = sum + ck(j)
            IF ( DABS(sum-sw)<DABS(sum)*eps ) GOTO 50
            sw = sum
         ENDDO
 50      IF ( ip==0 ) THEN
            r1f = sum/ck1
            R2f = -0.5D0*pi*qs*r1f
            R2d = qs*r1f + bk(1)
         ELSEIF ( ip==1 ) THEN
            r1d = sum/ck1
            R2f = bk(1)
            R2d = -0.5D0*pi*qs*r1d
         ENDIF
         RETURN
      ELSE
         CALL GMN(M,N,C,X,bk,gf,gd)
         CALL RMN1(M,N,C,X,Df,Kd,r1f,r1d)
         h0 = DATAN(X) - 0.5D0*pi
         R2f = qs*r1f*h0 + gf
         R2d = qs*(r1d*h0+r1f/(1.0D0+X*X)) + gd
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE BJNDD(N,X,Bj,Dj,Fj)
!
!       =====================================================
!       Purpose: Compute Bessel functions Jn(x) and their
!                first and second derivatives ( n= 0,1,… )
!       Input:   x ---  Argument of Jn(x)  ( x ≥ 0 )
!                n ---  Order of Jn(x)
!       Output:  BJ(n+1) ---  Jn(x)
!                DJ(n+1) ---  Jn'(x)
!                FJ(n+1) ---  Jn"(x)
!       =====================================================
!
      IMPLICIT NONE
!*--BJNDD10786
      DOUBLE PRECISION Bj , bs , Dj , f , f0 , f1 , Fj , X
      INTEGER k , m , mt , N , nt
      DIMENSION Bj(101) , Dj(101) , Fj(101)
      DO nt = 1 , 900
         mt = INT(0.5*LOG10(6.28*nt)-nt*LOG10(1.36*DABS(X)/nt))
         IF ( mt>20 ) GOTO 100
      ENDDO
 100  m = nt
      bs = 0.0D0
      f = 0.0D0
      f0 = 0.0D0
      f1 = 1.0D-35
      DO k = m , 0 , -1
         f = 2.0D0*(k+1.0D0)*f1/X - f0
         IF ( k<=N ) Bj(k+1) = f
         IF ( k==2*INT(k/2) ) bs = bs + 2.0D0*f
         f0 = f1
         f1 = f
      ENDDO
      DO k = 0 , N
         Bj(k+1) = Bj(k+1)/(bs-f)
      ENDDO
      Dj(1) = -Bj(2)
      Fj(1) = -1.0D0*Bj(1) - Dj(1)/X
      DO k = 1 , N
         Dj(k+1) = Bj(k) - k*Bj(k+1)/X
         Fj(k+1) = (k*k/(X*X)-1.0D0)*Bj(k+1) - Dj(k+1)/X
      ENDDO
      END
 
!       **********************************
 
 
      SUBROUTINE SPHJ(N,X,Nm,Sj,Dj)
!       MODIFIED to ALLOW N=0 CASE (ALSO IN SPHY)
!
!       =======================================================
!       Purpose: Compute spherical Bessel functions jn(x) and
!                their derivatives
!       Input :  x --- Argument of jn(x)
!                n --- Order of jn(x)  ( n = 0,1,… )
!       Output:  SJ(n) --- jn(x)
!                DJ(n) --- jn'(x)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       =======================================================
!
      IMPLICIT NONE
!*--SPHJ10840
      DOUBLE PRECISION cs , Dj , f , f0 , f1 , sa , sb , Sj , X
      INTEGER k , m , MSTA1 , MSTA2 , N , Nm
      DIMENSION Sj(0:N) , Dj(0:N)
      Nm = N
      IF ( DABS(X)<1.0D-100 ) THEN
         DO k = 0 , N
            Sj(k) = 0.0D0
            Dj(k) = 0.0D0
         ENDDO
         Sj(0) = 1.0D0
         IF ( N>0 ) Dj(1) = .3333333333333333D0
         RETURN
      ENDIF
      Sj(0) = DSIN(X)/X
      Dj(0) = (DCOS(X)-DSIN(X)/X)/X
      IF ( N<1 ) RETURN
      Sj(1) = (Sj(0)-DCOS(X))/X
      IF ( N>=2 ) THEN
         sa = Sj(0)
         sb = Sj(1)
         m = MSTA1(X,200)
         IF ( m<N ) THEN
            Nm = m
         ELSE
            m = MSTA2(X,N,15)
         ENDIF
         f = 0.0D0
         f0 = 0.0D0
         f1 = 1.0D0 - 100
         DO k = m , 0 , -1
            f = (2.0D0*k+3.0D0)*f1/X - f0
            IF ( k<=Nm ) Sj(k) = f
            f0 = f1
            f1 = f
         ENDDO
         cs = 0.0D0
         IF ( DABS(sa)>DABS(sb) ) cs = sa/f
         IF ( DABS(sa)<=DABS(sb) ) cs = sb/f0
         DO k = 0 , Nm
            Sj(k) = cs*Sj(k)
         ENDDO
      ENDIF
      DO k = 1 , Nm
         Dj(k) = Sj(k-1) - (k+1.0D0)*Sj(k)/X
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE OTHPL(Kf,N,X,Pl,Dpl)
!
!       ==========================================================
!       Purpose: Compute orthogonal polynomials: Tn(x) or Un(x),
!                or Ln(x) or Hn(x), and their derivatives
!       Input :  KF --- Function code
!                       KF=1 for Chebyshev polynomial Tn(x)
!                       KF=2 for Chebyshev polynomial Un(x)
!                       KF=3 for Laguerre polynomial Ln(x)
!                       KF=4 for Hermite polynomial Hn(x)
!                n ---  Order of orthogonal polynomials
!                x ---  Argument of orthogonal polynomials
!       Output:  PL(n) --- Tn(x) or Un(x) or Ln(x) or Hn(x)
!                DPL(n)--- Tn'(x) or Un'(x) or Ln'(x) or Hn'(x)
!       =========================================================
!
      IMPLICIT NONE
!*--OTHPL10912
      DOUBLE PRECISION a , b , c , Dpl , dy0 , dy1 , dyn , Pl , X , y0 ,&
                     & y1 , yn
      INTEGER k , Kf , N
      DIMENSION Pl(0:N) , Dpl(0:N)
      a = 2.0D0
      b = 0.0D0
      c = 1.0D0
      y0 = 1.0D0
      y1 = 2.0D0*X
      dy0 = 0.0D0
      dy1 = 2.0D0
      Pl(0) = 1.0D0
      Pl(1) = 2.0D0*X
      Dpl(0) = 0.0D0
      Dpl(1) = 2.0D0
      IF ( Kf==1 ) THEN
         y1 = X
         dy1 = 1.0D0
         Pl(1) = X
         Dpl(1) = 1.0D0
      ELSEIF ( Kf==3 ) THEN
         y1 = 1.0D0 - X
         dy1 = -1.0D0
         Pl(1) = 1.0D0 - X
         Dpl(1) = -1.0D0
      ENDIF
      DO k = 2 , N
         IF ( Kf==3 ) THEN
            a = -1.0D0/k
            b = 2.0D0 + a
            c = 1.0D0 + a
         ELSEIF ( Kf==4 ) THEN
            c = 2.0D0*(k-1.0D0)
         ENDIF
         yn = (a*X+b)*y1 - c*y0
         dyn = a*y1 + (a*X+b)*dy1 - c*dy0
         Pl(k) = yn
         Dpl(k) = dyn
         y0 = y1
         y1 = yn
         dy0 = dy1
         dy1 = dyn
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE KLVNZO(Nt,Kd,Zo)
!
!       ====================================================
!       Purpose: Compute the zeros of Kelvin functions
!       Input :  NT  --- Total number of zeros
!                KD  --- Function code
!                KD=1 to 8 for ber x, bei x, ker x, kei x,
!                          ber'x, bei'x, ker'x and kei'x,
!                          respectively.
!       Output:  ZO(M) --- the M-th zero of Kelvin function
!                          for code KD
!       Routine called:
!                KLVNA for computing Kelvin functions and
!                their derivatives
!       ====================================================
!
      IMPLICIT NONE
!*--KLVNZO10980
      DOUBLE PRECISION bei , ber , ddi , ddr , dei , der , gdi , gdr ,  &
                     & gei , ger , hei , her , rt , rt0 , Zo
      INTEGER Kd , m , Nt
      DIMENSION Zo(Nt) , rt0(8)
      rt0(1) = 2.84891
      rt0(2) = 5.02622
      rt0(3) = 1.71854
      rt0(4) = 3.91467
      rt0(5) = 6.03871
      rt0(6) = 3.77268
      rt0(7) = 2.66584
      rt0(8) = 4.93181
      rt = rt0(Kd)
      DO m = 1 , Nt
 50      CALL KLVNA(rt,ber,bei,ger,gei,der,dei,her,hei)
         IF ( Kd==1 ) THEN
            rt = rt - ber/der
         ELSEIF ( Kd==2 ) THEN
            rt = rt - bei/dei
         ELSEIF ( Kd==3 ) THEN
            rt = rt - ger/her
         ELSEIF ( Kd==4 ) THEN
            rt = rt - gei/hei
         ELSEIF ( Kd==5 ) THEN
            ddr = -bei - der/rt
            rt = rt - der/ddr
         ELSEIF ( Kd==6 ) THEN
            ddi = ber - dei/rt
            rt = rt - dei/ddi
         ELSEIF ( Kd==7 ) THEN
            gdr = -gei - her/rt
            rt = rt - her/gdr
         ELSE
            gdi = ger - hei/rt
            rt = rt - hei/gdi
         ENDIF
         IF ( DABS(rt-rt0(Kd))>5.0D-10 ) THEN
            rt0(Kd) = rt
            GOTO 50
         ENDIF
         Zo(m) = rt
         rt = rt + 4.44D0
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE RSWFO(M,N,C,X,Cv,Kf,R1f,R1d,R2f,R2d)
!
!       ==========================================================
!       Purpose: Compute oblate radial functions of the first
!                and second kinds, and their derivatives
!       Input :  m  --- Mode parameter,  m = 0,1,2,...
!                n  --- Mode parameter,  n = m,m+1,m+2,...
!                c  --- Spheroidal parameter
!                x  --- Argument (x ≥ 0)
!                cv --- Characteristic value
!                KF --- Function code
!                       KF=1 for the first kind
!                       KF=2 for the second kind
!                       KF=3 for both the first and second kinds
!       Output:  R1F --- Radial function of the first kind
!                R1D --- Derivative of the radial function of
!                        the first kind
!                R2F --- Radial function of the second kind
!                R2D --- Derivative of the radial function of
!                        the second kind
!       Routines called:
!            (1) SDMN for computing expansion coefficients dk
!            (2) RMN1 for computing prolate or oblate radial
!                function of the first kind
!            (3) RMN2L for computing prolate or oblate radial
!                function of the second kind for a large argument
!            (4) RMN2SO for computing oblate radial functions of
!                the second kind for a small argument
!       ==========================================================
!
      IMPLICIT NONE
!*--RSWFO11064
      DOUBLE PRECISION C , Cv , df , R1d , R1f , R2d , R2f , X
      INTEGER id , kd , Kf , M , N
      DIMENSION df(200)
      kd = -1
      CALL SDMN(M,N,C,Cv,kd,df)
      IF ( Kf/=2 ) CALL RMN1(M,N,C,X,df,kd,R1f,R1d)
      IF ( Kf>1 ) THEN
         id = 10
         IF ( X>1.0D-8 ) CALL RMN2L(M,N,C,X,df,kd,R2f,R2d,id)
         IF ( id>-1 ) CALL RMN2SO(M,N,C,X,Cv,df,kd,R2f,R2d)
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE CH12N(N,Z,Nm,Chf1,Chd1,Chf2,Chd2)
!
!       ====================================================
!       Purpose: Compute Hankel functions of the first and
!                second kinds and their derivatives for a
!                complex argument
!       Input :  z --- Complex argument
!                n --- Order of Hn(1)(z) and Hn(2)(z)
!       Output:  CHF1(n) --- Hn(1)(z)
!                CHD1(n) --- Hn(1)'(z)
!                CHF2(n) --- Hn(2)(z)
!                CHD2(n) --- Hn(2)'(z)
!                NM --- Highest order computed
!       Routines called:
!             (1) CJYNB for computing Jn(z) and Yn(z)
!             (2) CIKNB for computing In(z) and Kn(z)
!       ====================================================
!
      IMPLICIT NONE
!*--CH12N11104
      COMPLEX*16 cbi , cbj , cbk , cby , cdi , cdj , cdk , cdy , cf1 ,  &
               & cfac , Chd1 , Chd2 , Chf1 , Chf2 , ci , Z , zi
      INTEGER k , N , Nm
      DOUBLE PRECISION pi
      DIMENSION cbj(0:250) , cdj(0:250) , cby(0:250) , cdy(0:250) ,     &
              & cbi(0:250) , cdi(0:250) , cbk(0:250) , cdk(0:250)
      DIMENSION Chf1(0:N) , Chd1(0:N) , Chf2(0:N) , Chd2(0:N)
      ci = (0.0D0,1.0D0)
      pi = 3.141592653589793D0
      IF ( DIMAG(Z)<0.0D0 ) THEN
         CALL CJYNB(N,Z,Nm,cbj,cdj,cby,cdy)
         DO k = 0 , Nm
            Chf1(k) = cbj(k) + ci*cby(k)
            Chd1(k) = cdj(k) + ci*cdy(k)
         ENDDO
         zi = ci*Z
         CALL CIKNB(N,zi,Nm,cbi,cdi,cbk,cdk)
         cfac = -2.0D0/(pi*ci)
         DO k = 0 , Nm
            Chf2(k) = cfac*cbk(k)
            Chd2(k) = cfac*ci*cdk(k)
            cfac = cfac*ci
         ENDDO
      ELSEIF ( DIMAG(Z)>0.0D0 ) THEN
         zi = -ci*Z
         CALL CIKNB(N,zi,Nm,cbi,cdi,cbk,cdk)
         cf1 = -ci
         cfac = 2.0D0/(pi*ci)
         DO k = 0 , Nm
            Chf1(k) = cfac*cbk(k)
            Chd1(k) = -cfac*ci*cdk(k)
            cfac = cfac*cf1
         ENDDO
         CALL CJYNB(N,Z,Nm,cbj,cdj,cby,cdy)
         DO k = 0 , Nm
            Chf2(k) = cbj(k) - ci*cby(k)
            Chd2(k) = cdj(k) - ci*cdy(k)
         ENDDO
      ELSE
         CALL CJYNB(N,Z,Nm,cbj,cdj,cby,cdy)
         DO k = 0 , Nm
            Chf1(k) = cbj(k) + ci*cby(k)
            Chd1(k) = cdj(k) + ci*cdy(k)
            Chf2(k) = cbj(k) - ci*cby(k)
            Chd2(k) = cdj(k) - ci*cdy(k)
         ENDDO
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE JYZO(N,Nt,Rj0,Rj1,Ry0,Ry1)
!
!       ======================================================
!       Purpose: Compute the zeros of Bessel functions Jn(x),
!                Yn(x), and their derivatives
!       Input :  n  --- Order of Bessel functions  (n >= 0)
!                NT --- Number of zeros (roots)
!       Output:  RJ0(L) --- L-th zero of Jn(x),  L=1,2,...,NT
!                RJ1(L) --- L-th zero of Jn'(x), L=1,2,...,NT
!                RY0(L) --- L-th zero of Yn(x),  L=1,2,...,NT
!                RY1(L) --- L-th zero of Yn'(x), L=1,2,...,NT
!       Routine called: JYNDD for computing Jn(x), Yn(x), and
!                       their first and second derivatives
!       ======================================================
!
      IMPLICIT NONE
!*--JYZO11177
      DOUBLE PRECISION bjn , byn , djn , dyn , fjn , fyn , pi , Rj0 ,   &
                     & Rj1 , Ry0 , Ry1 , x , x0 , xguess
      INTEGER l , N , Nt
      DIMENSION Rj0(Nt) , Rj1(Nt) , Ry0(Nt) , Ry1(Nt)
      pi = 3.141592653589793D0
!       -- Newton method for j_{N,L}
!       1) initial guess for j_{N,1}
      IF ( N<=20 ) THEN
         x = 2.82141 + 1.15859*N
      ELSE
!          Abr & Stg (9.5.14)
         x = N + 1.85576*N**0.33333 + 1.03315/N**0.33333
      ENDIF
      l = 0
!       2) iterate
      xguess = x
 100  x0 = x
      CALL JYNDD(N,x,bjn,djn,fjn,byn,dyn,fyn)
      x = x - bjn/djn
      IF ( x-x0<-1 ) x = x0 - 1
      IF ( x-x0>1 ) x = x0 + 1
      IF ( DABS(x-x0)>1.0D-11 ) GOTO 100
!       3) initial guess for j_{N,L+1}
      IF ( l>=1 ) THEN
         IF ( x<=Rj0(l)+0.5 ) THEN
            x = xguess + pi
            xguess = x
            GOTO 100
         ENDIF
      ENDIF
      l = l + 1
      Rj0(l) = x
!       XXX: should have a better initial guess for large N ~> 100 here
      x = x + pi + MAX((0.0972D0+0.0679*N-0.000354*N**2)/l,0D0)
      IF ( l<Nt ) GOTO 100
!       -- Newton method for j_{N,L}'
      IF ( N<=20 ) THEN
         x = 0.961587 + 1.07703*N
      ELSE
         x = N + 0.80861*N**0.33333 + 0.07249/N**0.33333
      ENDIF
      IF ( N==0 ) x = 3.8317
      l = 0
      xguess = x
 200  x0 = x
      CALL JYNDD(N,x,bjn,djn,fjn,byn,dyn,fyn)
      x = x - djn/fjn
      IF ( x-x0<-1 ) x = x0 - 1
      IF ( x-x0>1 ) x = x0 + 1
      IF ( DABS(x-x0)>1.0D-11 ) GOTO 200
      IF ( l>=1 ) THEN
         IF ( x<=Rj1(l)+0.5 ) THEN
            x = xguess + pi
            xguess = x
            GOTO 200
         ENDIF
      ENDIF
      l = l + 1
      Rj1(l) = x
!       XXX: should have a better initial guess for large N ~> 100 here
      x = x + pi + MAX((0.4955D0+0.0915*N-0.000435*N**2)/l,0D0)
      IF ( l<Nt ) GOTO 200
!       -- Newton method for y_{N,L}
      IF ( N<=20 ) THEN
         x = 1.19477 + 1.08933*N
      ELSE
         x = N + 0.93158*N**0.33333 + 0.26035/N**0.33333
      ENDIF
      l = 0
      xguess = x
 300  x0 = x
      CALL JYNDD(N,x,bjn,djn,fjn,byn,dyn,fyn)
      x = x - byn/dyn
      IF ( x-x0<-1 ) x = x0 - 1
      IF ( x-x0>1 ) x = x0 + 1
      IF ( DABS(x-x0)>1.0D-11 ) GOTO 300
      IF ( l>=1 ) THEN
         IF ( x<=Ry0(l)+0.5 ) THEN
            x = xguess + pi
            xguess = x
            GOTO 300
         ENDIF
      ENDIF
      l = l + 1
      Ry0(l) = x
!       XXX: should have a better initial guess for large N ~> 100 here
      x = x + pi + MAX((0.312D0+0.0852*N-0.000403*N**2)/l,0D0)
      IF ( l<Nt ) GOTO 300
!       -- Newton method for y_{N,L}'
      IF ( N<=20 ) THEN
         x = 2.67257 + 1.16099*N
      ELSE
         x = N + 1.8211*N**0.33333 + 0.94001/N**0.33333
      ENDIF
      l = 0
      xguess = x
 400  x0 = x
      CALL JYNDD(N,x,bjn,djn,fjn,byn,dyn,fyn)
      x = x - dyn/fyn
      IF ( DABS(x-x0)>1.0D-11 ) GOTO 400
      IF ( l>=1 ) THEN
         IF ( x<=Ry1(l)+0.5 ) THEN
            x = xguess + pi
            xguess = x
            GOTO 400
         ENDIF
      ENDIF
      l = l + 1
      Ry1(l) = x
!       XXX: should have a better initial guess for large N ~> 100 here
      x = x + pi + MAX((0.197D0+0.0643*N-0.000286*N**2)/l,0D0)
      IF ( l<Nt ) GOTO 400
      END
 
 
 
!       **********************************
 
      SUBROUTINE IKV(V,X,Vm,Bi,Di,Bk,Dk)
!
!       =======================================================
!       Purpose: Compute modified Bessel functions Iv(x) and
!                Kv(x), and their derivatives
!       Input :  x --- Argument ( x ≥ 0 )
!                v --- Order of Iv(x) and Kv(x)
!                      ( v = n+v0, n = 0,1,2,..., 0 ≤ v0 < 1 )
!       Output:  BI(n) --- In+v0(x)
!                DI(n) --- In+v0'(x)
!                BK(n) --- Kn+v0(x)
!                DK(n) --- Kn+v0'(x)
!                VM --- Highest order computed
!       Routines called:
!            (1) GAMMA2 for computing the gamma function
!            (2) MSTA1 and MSTA2 to compute the starting
!                point for backward recurrence
!       =======================================================
!
      IMPLICIT NONE
!*--IKV11319
      DOUBLE PRECISION a1 , a2 , Bi , bi0 , Bk , bk0 , bk1 , bk2 , ca , &
                     & cb , cs , ct , Di , Dk , f , f1 , f2 , gan ,     &
                     & gap , pi
      DOUBLE PRECISION piv , r , r1 , r2 , sum , V , v0 , v0n , v0p ,   &
                     & Vm , vt , w0 , wa , ww , X , x2
      INTEGER k , k0 , m , MSTA1 , MSTA2 , n
      DIMENSION Bi(0:*) , Di(0:*) , Bk(0:*) , Dk(0:*)
      pi = 3.141592653589793D0
      x2 = X*X
      n = INT(V)
      v0 = V - n
      IF ( n==0 ) n = 1
      IF ( X<1.0D-100 ) THEN
         DO k = 0 , n
            Bi(k) = 0.0D0
            Di(k) = 0.0D0
            Bk(k) = -1.0D+300
            Dk(k) = 1.0D+300
         ENDDO
         IF ( V==0.0 ) THEN
            Bi(0) = 1.0D0
            Di(1) = 0.5D0
         ENDIF
         Vm = V
         RETURN
      ENDIF
      piv = pi*v0
      vt = 4.0D0*v0*v0
      IF ( v0==0.0D0 ) THEN
         a1 = 1.0D0
      ELSE
         v0p = 1.0D0 + v0
         CALL GAMMA2(v0p,gap)
         a1 = (0.5D0*X)**v0/gap
      ENDIF
      k0 = 14
      IF ( X>=35.0 ) k0 = 10
      IF ( X>=50.0 ) k0 = 8
      IF ( X<=18.0 ) THEN
         bi0 = 1.0D0
         r = 1.0D0
         DO k = 1 , 30
            r = 0.25D0*r*x2/(k*(k+v0))
            bi0 = bi0 + r
            IF ( DABS(r/bi0)<1.0D-15 ) GOTO 50
         ENDDO
 50      bi0 = bi0*a1
      ELSE
         ca = EXP(X)/DSQRT(2.0D0*pi*X)
         sum = 1.0D0
         r = 1.0D0
         DO k = 1 , k0
            r = -0.125D0*r*(vt-(2.0D0*k-1.0D0)**2.0)/(k*X)
            sum = sum + r
         ENDDO
         bi0 = ca*sum
      ENDIF
      m = MSTA1(X,200)
      IF ( m<n ) THEN
         n = m
      ELSE
         m = MSTA2(X,n,15)
      ENDIF
      f = 0.0D0
      f2 = 0.0D0
      f1 = 1.0D-100
      ww = 0.0D0
      DO k = m , 0 , -1
         f = 2.0D0*(v0+k+1.0D0)/X*f1 + f2
         IF ( k<=n ) Bi(k) = f
         f2 = f1
         f1 = f
      ENDDO
      cs = bi0/f
      DO k = 0 , n
         Bi(k) = cs*Bi(k)
      ENDDO
      Di(0) = v0/X*Bi(0) + Bi(1)
      DO k = 1 , n
         Di(k) = -(k+v0)/X*Bi(k) + Bi(k-1)
      ENDDO
      IF ( X>9.0D0 ) THEN
         cb = EXP(-X)*DSQRT(0.5D0*pi/X)
         sum = 1.0D0
         r = 1.0D0
         DO k = 1 , k0
            r = 0.125D0*r*(vt-(2.0*k-1.0)**2.0)/(k*X)
            sum = sum + r
         ENDDO
         bk0 = cb*sum
      ELSEIF ( v0==0.0D0 ) THEN
         ct = -DLOG(0.5D0*X) - 0.5772156649015329D0
         cs = 0.0D0
         w0 = 0.0D0
         r = 1.0D0
         DO k = 1 , 50
            w0 = w0 + 1.0D0/k
            r = 0.25D0*r/(k*k)*x2
            cs = cs + r*(w0+ct)
            wa = DABS(cs)
            IF ( DABS((wa-ww)/wa)<1.0D-15 ) GOTO 100
            ww = wa
         ENDDO
 100     bk0 = ct + cs
      ELSE
         v0n = 1.0D0 - v0
         CALL GAMMA2(v0n,gan)
         a2 = 1.0D0/(gan*(0.5D0*X)**v0)
         a1 = (0.5D0*X)**v0/gap
         sum = a2 - a1
         r1 = 1.0D0
         r2 = 1.0D0
         DO k = 1 , 120
            r1 = 0.25D0*r1*x2/(k*(k-v0))
            r2 = 0.25D0*r2*x2/(k*(k+v0))
            sum = sum + a2*r1 - a1*r2
            wa = DABS(sum)
            IF ( DABS((wa-ww)/wa)<1.0D-15 ) GOTO 150
            ww = wa
         ENDDO
 150     bk0 = 0.5D0*pi*sum/DSIN(piv)
      ENDIF
      bk1 = (1.0D0/X-Bi(1)*bk0)/Bi(0)
      Bk(0) = bk0
      Bk(1) = bk1
      DO k = 2 , n
         bk2 = 2.0D0*(v0+k-1.0D0)/X*bk1 + bk0
         Bk(k) = bk2
         bk0 = bk1
         bk1 = bk2
      ENDDO
      Dk(0) = v0/X*Bk(0) - Bk(1)
      DO k = 1 , n
         Dk(k) = -(k+v0)/X*Bk(k) - Bk(k-1)
      ENDDO
      Vm = n + v0
      END
 
 
 
!       **********************************
 
      SUBROUTINE SDMN(M,N,C,Cv,Kd,Df)
!
!       =====================================================
!       Purpose: Compute the expansion coefficients of the
!                prolate and oblate spheroidal functions, dk
!       Input :  m  --- Mode parameter
!                n  --- Mode parameter
!                c  --- Spheroidal parameter
!                cv --- Characteristic value
!                KD --- Function code
!                       KD=1 for prolate; KD=-1 for oblate
!       Output:  DF(k) --- Expansion coefficients dk;
!                          DF(1), DF(2), ... correspond to
!                          d0, d2, ... for even n-m and d1,
!                          d3, ... for odd n-m
!       =====================================================
!
      IMPLICIT NONE
!*--SDMN11483
      DOUBLE PRECISION a , C , cs , Cv , d , d2k , Df , dk0 , dk1 ,     &
                     & dk2 , f , f0 , f1 , f2 , fl , fs , g , r1 , r3 , &
                     & r4
      DOUBLE PRECISION s0 , su1 , su2 , sw
      INTEGER i , ip , j , k , k1 , kb , Kd , M , N , nm
      DIMENSION a(200) , d(200) , g(200) , Df(200)
      nm = 25 + INT(0.5*(N-M)+C)
      IF ( C<1.0D-10 ) THEN
         DO i = 1 , nm
            Df(i) = 0D0
         ENDDO
         Df((N-M)/2+1) = 1.0D0
         RETURN
      ENDIF
      cs = C*C*Kd
      ip = 1
      k = 0
      IF ( N-M==2*INT((N-M)/2) ) ip = 0
      DO i = 1 , nm + 2
         IF ( ip==0 ) k = 2*(i-1)
         IF ( ip==1 ) k = 2*i - 1
         dk0 = M + k
         dk1 = M + k + 1
         dk2 = 2*(M+k)
         d2k = 2*M + k
         a(i) = (d2k+2.0)*(d2k+1.0)/((dk2+3.0)*(dk2+5.0))*cs
         d(i) = dk0*dk1 + (2.0*dk0*dk1-2.0*M*M-1.0)                     &
              & /((dk2-1.0)*(dk2+3.0))*cs
         g(i) = k*(k-1.0)/((dk2-3.0)*(dk2-1.0))*cs
      ENDDO
      fs = 1.0D0
      f1 = 0.0D0
      f0 = 1.0D-100
      kb = 0
      Df(nm+1) = 0.0D0
      fl = 0.0D0
      DO k = nm , 1 , -1
         f = -((d(k+1)-Cv)*f0+a(k+1)*f1)/g(k+1)
         IF ( DABS(f)>DABS(Df(k+1)) ) THEN
            Df(k) = f
            f1 = f0
            f0 = f
            IF ( DABS(f)>1.0D+100 ) THEN
               DO k1 = k , nm
                  Df(k1) = Df(k1)*1.0D-100
               ENDDO
               f1 = f1*1.0D-100
               f0 = f0*1.0D-100
            ENDIF
         ELSE
            kb = k
            fl = Df(k+1)
            f1 = 1.0D-100
            f2 = -(d(1)-Cv)/a(1)*f1
            Df(1) = f1
            IF ( kb==1 ) THEN
               fs = f2
            ELSEIF ( kb==2 ) THEN
               Df(2) = f2
               fs = -((d(2)-Cv)*f2+g(2)*f1)/a(2)
            ELSE
               Df(2) = f2
               DO j = 3 , kb + 1
                  f = -((d(j-1)-Cv)*f2+g(j-1)*f1)/a(j-1)
                  IF ( j<=kb ) Df(j) = f
                  IF ( DABS(f)>1.0D+100 ) THEN
                     DO k1 = 1 , j
                        Df(k1) = Df(k1)*1.0D-100
                     ENDDO
                     f = f*1.0D-100
                     f2 = f2*1.0D-100
                  ENDIF
                  f1 = f2
                  f2 = f
               ENDDO
               fs = f
            ENDIF
            GOTO 100
         ENDIF
      ENDDO
 100  su1 = 0.0D0
      r1 = 1.0D0
      DO j = M + ip + 1 , 2*(M+ip)
         r1 = r1*j
      ENDDO
      su1 = Df(1)*r1
      DO k = 2 , kb
         r1 = -r1*(k+M+ip-1.5D0)/(k-1.0D0)
         su1 = su1 + r1*Df(k)
      ENDDO
      su2 = 0.0D0
      sw = 0.0D0
      DO k = kb + 1 , nm
         IF ( k/=1 ) r1 = -r1*(k+M+ip-1.5D0)/(k-1.0D0)
         su2 = su2 + r1*Df(k)
         IF ( DABS(sw-su2)<DABS(su2)*1.0D-14 ) GOTO 200
         sw = su2
      ENDDO
 200  r3 = 1.0D0
      DO j = 1 , (M+N+ip)/2
         r3 = r3*(j+0.5D0*(N+M+ip))
      ENDDO
      r4 = 1.0D0
      DO j = 1 , (N-M-ip)/2
         r4 = -4.0D0*r4*j
      ENDDO
      s0 = r3/(fl*(su1/fs)+su2)/r4
      DO k = 1 , kb
         Df(k) = fl/fs*s0*Df(k)
      ENDDO
      DO k = kb + 1 , nm
         Df(k) = s0*Df(k)
      ENDDO
      END
 
 
 
 
!       **********************************
 
      SUBROUTINE AJYIK(X,Vj1,Vj2,Vy1,Vy2,Vi1,Vi2,Vk1,Vk2)
!
!       =======================================================
!       Purpose: Compute Bessel functions Jv(x) and Yv(x),
!                and modified Bessel functions Iv(x) and
!                Kv(x), and their derivatives with v=1/3,2/3
!       Input :  x --- Argument of Jv(x),Yv(x),Iv(x) and
!                      Kv(x) ( x ≥ 0 )
!       Output:  VJ1 --- J1/3(x)
!                VJ2 --- J2/3(x)
!                VY1 --- Y1/3(x)
!                VY2 --- Y2/3(x)
!                VI1 --- I1/3(x)
!                VI2 --- I2/3(x)
!                VK1 --- K1/3(x)
!                VK2 --- K2/3(x)
!       =======================================================
!
      IMPLICIT NONE
!*--AJYIK11626
      DOUBLE PRECISION a0 , b0 , c0 , ck , gn , gn1 , gn2 , gp1 , gp2 , &
                     & pi , pv1 , pv2 , px , qx , r , rp , rp2 , rq ,   &
                     & sk , sum
      DOUBLE PRECISION uj1 , uj2 , uu0 , Vi1 , Vi2 , vil , Vj1 , Vj2 ,  &
                     & vjl , Vk1 , Vk2 , vl , vsl , vv , vv0 , Vy1 ,    &
                     & Vy2 , X , x2 , xk
      INTEGER k , k0 , l
      IF ( X==0.0D0 ) THEN
         Vj1 = 0.0D0
         Vj2 = 0.0D0
         Vy1 = -1.0D+300
         Vy2 = 1.0D+300
         Vi1 = 0.0D0
         Vi2 = 0.0D0
         Vk1 = -1.0D+300
         Vk2 = -1.0D+300
         RETURN
      ENDIF
      pi = 3.141592653589793D0
      rp2 = .63661977236758D0
      gp1 = .892979511569249D0
      gp2 = .902745292950934D0
      gn1 = 1.3541179394264D0
      gn2 = 2.678938534707747D0
      vv0 = 0.444444444444444D0
      uu0 = 1.1547005383793D0
      x2 = X*X
      k0 = 12
      IF ( X>=35.0 ) k0 = 10
      IF ( X>=50.0 ) k0 = 8
      IF ( X<=12.0 ) THEN
         DO l = 1 , 2
            vl = l/3.0D0
            vjl = 1.0D0
            r = 1.0D0
            DO k = 1 , 40
               r = -0.25D0*r*x2/(k*(k+vl))
               vjl = vjl + r
               IF ( DABS(r)<1.0D-15 ) GOTO 20
            ENDDO
 20         a0 = (0.5D0*X)**vl
            IF ( l==1 ) Vj1 = a0/gp1*vjl
            IF ( l==2 ) Vj2 = a0/gp2*vjl
         ENDDO
      ELSE
         DO l = 1 , 2
            vv = vv0*l*l
            px = 1.0D0
            rp = 1.0D0
            DO k = 1 , k0
               rp = -0.78125D-2*rp*(vv-(4.0*k-3.0)**2.0)                &
                  & *(vv-(4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*x2)
               px = px + rp
            ENDDO
            qx = 1.0D0
            rq = 1.0D0
            DO k = 1 , k0
               rq = -0.78125D-2*rq*(vv-(4.0*k-1.0)**2.0)                &
                  & *(vv-(4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*x2)
               qx = qx + rq
            ENDDO
            qx = 0.125D0*(vv-1.0)*qx/X
            xk = X - (0.5D0*l/3.0D0+0.25D0)*pi
            a0 = DSQRT(rp2/X)
            ck = DCOS(xk)
            sk = DSIN(xk)
            IF ( l==1 ) THEN
               Vj1 = a0*(px*ck-qx*sk)
               Vy1 = a0*(px*sk+qx*ck)
            ELSEIF ( l==2 ) THEN
               Vj2 = a0*(px*ck-qx*sk)
               Vy2 = a0*(px*sk+qx*ck)
            ENDIF
         ENDDO
      ENDIF
      IF ( X<=12.0D0 ) THEN
         uj1 = 0.0D0
         uj2 = 0.0D0
         DO l = 1 , 2
            vl = l/3.0D0
            vjl = 1.0D0
            r = 1.0D0
            DO k = 1 , 40
               r = -0.25D0*r*x2/(k*(k-vl))
               vjl = vjl + r
               IF ( DABS(r)<1.0D-15 ) GOTO 40
            ENDDO
 40         b0 = (2.0D0/X)**vl
            IF ( l==1 ) uj1 = b0*vjl/gn1
            IF ( l==2 ) uj2 = b0*vjl/gn2
         ENDDO
         pv1 = pi/3.0D0
         pv2 = pi/1.5D0
         Vy1 = uu0*(Vj1*DCOS(pv1)-uj1)
         Vy2 = uu0*(Vj2*DCOS(pv2)-uj2)
      ENDIF
      IF ( X<=18.0 ) THEN
         DO l = 1 , 2
            vl = l/3.0D0
            vil = 1.0D0
            r = 1.0D0
            DO k = 1 , 40
               r = 0.25D0*r*x2/(k*(k+vl))
               vil = vil + r
               IF ( DABS(r)<1.0D-15 ) GOTO 60
            ENDDO
 60         a0 = (0.5D0*X)**vl
            IF ( l==1 ) Vi1 = a0/gp1*vil
            IF ( l==2 ) Vi2 = a0/gp2*vil
         ENDDO
      ELSE
         c0 = EXP(X)/DSQRT(2.0D0*pi*X)
         DO l = 1 , 2
            vv = vv0*l*l
            vsl = 1.0D0
            r = 1.0D0
            DO k = 1 , k0
               r = -0.125D0*r*(vv-(2.0D0*k-1.0D0)**2.0)/(k*X)
               vsl = vsl + r
            ENDDO
            IF ( l==1 ) Vi1 = c0*vsl
            IF ( l==2 ) Vi2 = c0*vsl
         ENDDO
      ENDIF
      IF ( X<=9.0D0 ) THEN
         gn = 0.0D0
         DO l = 1 , 2
            vl = l/3.0D0
            IF ( l==1 ) gn = gn1
            IF ( l==2 ) gn = gn2
            a0 = (2.0D0/X)**vl/gn
            sum = 1.0D0
            r = 1.0D0
            DO k = 1 , 60
               r = 0.25D0*r*x2/(k*(k-vl))
               sum = sum + r
               IF ( DABS(r)<1.0D-15 ) GOTO 80
            ENDDO
 80         IF ( l==1 ) Vk1 = 0.5D0*uu0*pi*(sum*a0-Vi1)
            IF ( l==2 ) Vk2 = 0.5D0*uu0*pi*(sum*a0-Vi2)
         ENDDO
      ELSE
         c0 = EXP(-X)*DSQRT(0.5D0*pi/X)
         DO l = 1 , 2
            vv = vv0*l*l
            sum = 1.0D0
            r = 1.0D0
            DO k = 1 , k0
               r = 0.125D0*r*(vv-(2.0*k-1.0)**2.0)/(k*X)
               sum = sum + r
            ENDDO
            IF ( l==1 ) Vk1 = c0*sum
            IF ( l==2 ) Vk2 = c0*sum
         ENDDO
      ENDIF
      END
 
 
 
!       **********************************
 
      SUBROUTINE CIKVB(V,Z,Vm,Cbi,Cdi,Cbk,Cdk)
!
!       ===========================================================
!       Purpose: Compute the modified Bessel functions Iv(z), Kv(z)
!                and their derivatives for an arbitrary order and
!                complex argument
!       Input :  z --- Complex argument z
!                v --- Real order of Iv(z) and Kv(z)
!                      ( v =n+v0, n = 0,1,2,..., 0 ≤ v0 < 1 )
!       Output:  CBI(n) --- In+v0(z)
!                CDI(n) --- In+v0'(z)
!                CBK(n) --- Kn+v0(z)
!                CDK(n) --- Kn+v0'(z)
!                VM --- Highest order computed
!       Routines called:
!            (1) GAMMA2 for computing the gamma function
!            (2) MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ===========================================================
!
      IMPLICIT NONE
!*--CIKVB11812
      DOUBLE PRECISION a0 , gan , gap , pi , piv , V , v0 , v0n , v0p , &
                     & Vm , vt , w0
      COMPLEX*16 ca , ca1 , ca2 , cb , Cbi , cbi0 , Cbk , cbk0 , Cdi ,  &
               & Cdk , cf , cf1 , cf2 , ci , ci0 , ckk , cp , cr , cr1 ,&
               & cr2
      COMPLEX*16 cs , csu , ct , cvk , Z , z1 , z2
      INTEGER k , k0 , m , MSTA1 , MSTA2 , n
      DIMENSION Cbi(0:*) , Cdi(0:*) , Cbk(0:*) , Cdk(0:*)
      z1 = Z
      z2 = Z*Z
      a0 = ABS(Z)
      pi = 3.141592653589793D0
      ci = (0.0D0,1.0D0)
      n = INT(V)
      v0 = V - n
      piv = pi*v0
      vt = 4.0D0*v0*v0
      IF ( n==0 ) n = 1
      IF ( a0<1.0D-100 ) THEN
         DO k = 0 , n
            Cbi(k) = 0.0D0
            Cdi(k) = 0.0D0
            Cbk(k) = -1.0D+300
            Cdk(k) = 1.0D+300
         ENDDO
         IF ( v0==0.0 ) THEN
            Cbi(0) = (1.0D0,0.0D0)
            Cdi(1) = (0.5D0,0.0D0)
         ENDIF
         Vm = V
         RETURN
      ENDIF
      k0 = 14
      IF ( a0>=35.0 ) k0 = 10
      IF ( a0>=50.0 ) k0 = 8
      IF ( DBLE(Z)<0.0 ) z1 = -Z
      IF ( a0<18.0 ) THEN
         IF ( v0==0.0 ) THEN
            ca1 = (1.0D0,0.0D0)
         ELSE
            v0p = 1.0D0 + v0
            CALL GAMMA2(v0p,gap)
            ca1 = (0.5D0*z1)**v0/gap
         ENDIF
         ci0 = (1.0D0,0.0D0)
         cr = (1.0D0,0.0D0)
         DO k = 1 , 50
            cr = 0.25D0*cr*z2/(k*(k+v0))
            ci0 = ci0 + cr
            IF ( ABS(cr/ci0)<1.0D-15 ) GOTO 50
         ENDDO
 50      cbi0 = ci0*ca1
      ELSE
         ca = EXP(z1)/SQRT(2.0D0*pi*z1)
         cs = (1.0D0,0.0D0)
         cr = (1.0D0,0.0D0)
         DO k = 1 , k0
            cr = -0.125D0*cr*(vt-(2.0D0*k-1.0D0)**2.0)/(k*z1)
            cs = cs + cr
         ENDDO
         cbi0 = ca*cs
      ENDIF
      m = MSTA1(a0,200)
      IF ( m<n ) THEN
         n = m
      ELSE
         m = MSTA2(a0,n,15)
      ENDIF
      cf2 = (0.0D0,0.0D0)
      cf1 = (1.0D-100,0.0D0)
      DO k = m , 0 , -1
         cf = 2.0D0*(v0+k+1.0D0)/z1*cf1 + cf2
         IF ( k<=n ) Cbi(k) = cf
         cf2 = cf1
         cf1 = cf
      ENDDO
      cs = cbi0/cf
      DO k = 0 , n
         Cbi(k) = cs*Cbi(k)
      ENDDO
      IF ( a0>9.0 ) THEN
         cb = EXP(-z1)*SQRT(0.5D0*pi/z1)
         cs = (1.0D0,0.0D0)
         cr = (1.0D0,0.0D0)
         DO k = 1 , k0
            cr = 0.125D0*cr*(vt-(2.0D0*k-1.0D0)**2.0)/(k*z1)
            cs = cs + cr
         ENDDO
         cbk0 = cb*cs
      ELSEIF ( v0==0.0 ) THEN
         ct = -LOG(0.5D0*z1) - 0.5772156649015329D0
         cs = (0.0D0,0.0D0)
         w0 = 0.0D0
         cr = (1.0D0,0.0D0)
         DO k = 1 , 50
            w0 = w0 + 1.0D0/k
            cr = 0.25D0*cr/(k*k)*z2
            cp = cr*(w0+ct)
            cs = cs + cp
            IF ( k>=10 .AND. ABS(cp/cs)<1.0D-15 ) GOTO 100
         ENDDO
 100     cbk0 = ct + cs
      ELSE
         v0n = 1.0D0 - v0
         CALL GAMMA2(v0n,gan)
         ca2 = 1.0D0/(gan*(0.5D0*z1)**v0)
         ca1 = (0.5D0*z1)**v0/gap
         csu = ca2 - ca1
         cr1 = (1.0D0,0.0D0)
         cr2 = (1.0D0,0.0D0)
         DO k = 1 , 50
            cr1 = 0.25D0*cr1*z2/(k*(k-v0))
            cr2 = 0.25D0*cr2*z2/(k*(k+v0))
            cp = ca2*cr1 - ca1*cr2
            csu = csu + cp
            IF ( k>=10 .AND. ABS(cp/csu)<1.0D-15 ) GOTO 150
         ENDDO
 150     cbk0 = 0.5D0*pi*csu/DSIN(piv)
      ENDIF
      Cbk(0) = cbk0
      IF ( DBLE(Z)<0.0 ) THEN
         DO k = 0 , n
            cvk = EXP((k+v0)*pi*ci)
            IF ( DIMAG(Z)<0.0D0 ) THEN
               Cbk(k) = cvk*Cbk(k) + pi*ci*Cbi(k)
               Cbi(k) = Cbi(k)/cvk
            ELSEIF ( DIMAG(Z)>0.0 ) THEN
               Cbk(k) = Cbk(k)/cvk - pi*ci*Cbi(k)
               Cbi(k) = cvk*Cbi(k)
            ENDIF
         ENDDO
      ENDIF
      DO k = 1 , n
         ckk = (1.0D0/Z-Cbi(k)*Cbk(k-1))/Cbi(k-1)
         Cbk(k) = ckk
      ENDDO
      Cdi(0) = v0/Z*Cbi(0) + Cbi(1)
      Cdk(0) = v0/Z*Cbk(0) - Cbk(1)
      DO k = 1 , n
         Cdi(k) = -(k+v0)/Z*Cbi(k) + Cbi(k-1)
         Cdk(k) = -(k+v0)/Z*Cbk(k) - Cbk(k-1)
      ENDDO
      Vm = n + v0
      END
 
 
 
!       **********************************
 
      SUBROUTINE CIKVA(V,Z,Vm,Cbi,Cdi,Cbk,Cdk)
!
!       ============================================================
!       Purpose: Compute the modified Bessel functions Iv(z), Kv(z)
!                and their derivatives for an arbitrary order and
!                complex argument
!       Input :  z --- Complex argument
!                v --- Real order of Iv(z) and Kv(z)
!                      ( v = n+v0, n = 0,1,2,…, 0 ≤ v0 < 1 )
!       Output:  CBI(n) --- In+v0(z)
!                CDI(n) --- In+v0'(z)
!                CBK(n) --- Kn+v0(z)
!                CDK(n) --- Kn+v0'(z)
!                VM --- Highest order computed
!       Routines called:
!            (1) GAMMA2 for computing the gamma function
!            (2) MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ============================================================
!
      IMPLICIT NONE
!*--CIKVA11986
      DOUBLE PRECISION a0 , gan , gap , pi , piv , V , v0 , v0n , v0p , &
                     & Vm , vt , w0 , ws , ws0
      COMPLEX*16 ca , ca1 , ca2 , cb , Cbi , cbi0 , Cbk , cbk0 , cbk1 , &
               & Cdi , Cdk , cf , cf1 , cf2 , cg0 , cg1 , cgk , ci ,    &
               & ci0 , cp
      COMPLEX*16 cr , cr1 , cr2 , cs , csu , ct , cvk , Z , z1 , z2
      INTEGER k , k0 , m , MSTA1 , MSTA2 , n
      DIMENSION Cbi(0:*) , Cdi(0:*) , Cbk(0:*) , Cdk(0:*)
      pi = 3.141592653589793D0
      ci = (0.0D0,1.0D0)
      a0 = ABS(Z)
      z1 = Z
      z2 = Z*Z
      n = INT(V)
      v0 = V - n
      piv = pi*v0
      vt = 4.0D0*v0*v0
      IF ( n==0 ) n = 1
      IF ( a0<1.0D-100 ) THEN
         DO k = 0 , n
            Cbi(k) = 0.0D0
            Cdi(k) = 0.0D0
            Cbk(k) = -1.0D+300
            Cdk(k) = 1.0D+300
         ENDDO
         IF ( v0==0.0 ) THEN
            Cbi(0) = (1.0D0,0.0D0)
            Cdi(1) = (0.5D0,0.0D0)
         ENDIF
         Vm = V
         RETURN
      ENDIF
      k0 = 14
      IF ( a0>=35.0 ) k0 = 10
      IF ( a0>=50.0 ) k0 = 8
      IF ( DBLE(Z)<0.0 ) z1 = -Z
      IF ( a0<18.0 ) THEN
         IF ( v0==0.0 ) THEN
            ca1 = (1.0D0,0.0D0)
         ELSE
            v0p = 1.0D0 + v0
            CALL GAMMA2(v0p,gap)
            ca1 = (0.5D0*z1)**v0/gap
         ENDIF
         ci0 = (1.0D0,0.0D0)
         cr = (1.0D0,0.0D0)
         DO k = 1 , 50
            cr = 0.25D0*cr*z2/(k*(k+v0))
            ci0 = ci0 + cr
            IF ( ABS(cr)<ABS(ci0)*1.0D-15 ) GOTO 50
         ENDDO
 50      cbi0 = ci0*ca1
      ELSE
         ca = EXP(z1)/SQRT(2.0D0*pi*z1)
         cs = (1.0D0,0.0D0)
         cr = (1.0D0,0.0D0)
         DO k = 1 , k0
            cr = -0.125D0*cr*(vt-(2.0D0*k-1.0D0)**2.0)/(k*z1)
            cs = cs + cr
         ENDDO
         cbi0 = ca*cs
      ENDIF
      m = MSTA1(a0,200)
      IF ( m<n ) THEN
         n = m
      ELSE
         m = MSTA2(a0,n,15)
      ENDIF
      cf2 = (0.0D0,0.0D0)
      cf1 = (1.0D-100,0.0D0)
      DO k = m , 0 , -1
         cf = 2.0D0*(v0+k+1.0D0)/z1*cf1 + cf2
         IF ( k<=n ) Cbi(k) = cf
         cf2 = cf1
         cf1 = cf
      ENDDO
      cs = cbi0/cf
      DO k = 0 , n
         Cbi(k) = cs*Cbi(k)
      ENDDO
      IF ( a0>9.0 ) THEN
         cb = EXP(-z1)*SQRT(0.5D0*pi/z1)
         cs = (1.0D0,0.0D0)
         cr = (1.0D0,0.0D0)
         DO k = 1 , k0
            cr = 0.125D0*cr*(vt-(2.0D0*k-1.0D0)**2.0)/(k*z1)
            cs = cs + cr
         ENDDO
         cbk0 = cb*cs
      ELSEIF ( v0==0.0 ) THEN
         ct = -LOG(0.5D0*z1) - 0.5772156649015329D0
         cs = (0.0D0,0.0D0)
         w0 = 0.0D0
         cr = (1.0D0,0.0D0)
         DO k = 1 , 50
            w0 = w0 + 1.0D0/k
            cr = 0.25D0*cr/(k*k)*z2
            cp = cr*(w0+ct)
            cs = cs + cp
            IF ( k>=10 .AND. ABS(cp/cs)<1.0D-15 ) GOTO 100
         ENDDO
 100     cbk0 = ct + cs
      ELSE
         v0n = 1.0D0 - v0
         CALL GAMMA2(v0n,gan)
         ca2 = 1.0D0/(gan*(0.5D0*z1)**v0)
         ca1 = (0.5D0*z1)**v0/gap
         csu = ca2 - ca1
         cr1 = (1.0D0,0.0D0)
         cr2 = (1.0D0,0.0D0)
         ws0 = 0.0D0
         DO k = 1 , 50
            cr1 = 0.25D0*cr1*z2/(k*(k-v0))
            cr2 = 0.25D0*cr2*z2/(k*(k+v0))
            csu = csu + ca2*cr1 - ca1*cr2
            ws = ABS(csu)
            IF ( k>=10 .AND. DABS(ws-ws0)/ws<1.0D-15 ) GOTO 150
            ws0 = ws
         ENDDO
 150     cbk0 = 0.5D0*pi*csu/DSIN(piv)
      ENDIF
      cbk1 = (1.0D0/z1-Cbi(1)*cbk0)/Cbi(0)
      Cbk(0) = cbk0
      Cbk(1) = cbk1
      cg0 = cbk0
      cg1 = cbk1
      DO k = 2 , n
         cgk = 2.0D0*(v0+k-1.0D0)/z1*cg1 + cg0
         Cbk(k) = cgk
         cg0 = cg1
         cg1 = cgk
      ENDDO
      IF ( DBLE(Z)<0.0 ) THEN
         DO k = 0 , n
            cvk = EXP((k+v0)*pi*ci)
            IF ( DIMAG(Z)<0.0D0 ) THEN
               Cbk(k) = cvk*Cbk(k) + pi*ci*Cbi(k)
               Cbi(k) = Cbi(k)/cvk
            ELSEIF ( DIMAG(Z)>0.0 ) THEN
               Cbk(k) = Cbk(k)/cvk - pi*ci*Cbi(k)
               Cbi(k) = cvk*Cbi(k)
            ENDIF
         ENDDO
      ENDIF
      Cdi(0) = v0/Z*Cbi(0) + Cbi(1)
      Cdk(0) = v0/Z*Cbk(0) - Cbk(1)
      DO k = 1 , n
         Cdi(k) = -(k+v0)/Z*Cbi(k) + Cbi(k-1)
         Cdk(k) = -(k+v0)/Z*Cbk(k) - Cbk(k-1)
      ENDDO
      Vm = n + v0
      END
 
 
 
!       **********************************
 
      SUBROUTINE CFC(Z,Zf,Zd)
!
!       =========================================================
!       Purpose: Compute complex Fresnel integral C(z) and C'(z)
!       Input :  z --- Argument of C(z)
!       Output:  ZF --- C(z)
!                ZD --- C'(z)
!       =========================================================
!
      IMPLICIT NONE
!*--CFC12157
      COMPLEX*16 c , cf , cf0 , cf1 , cg , cr , d , Z , z0 , Zd , Zf ,  &
               & zp , zp2
      DOUBLE PRECISION eps , pi , w0 , wa , wa0
      INTEGER k , m
      eps = 1.0D-14
      pi = 3.141592653589793D0
      w0 = ABS(Z)
      zp = 0.5D0*pi*Z*Z
      zp2 = zp*zp
      z0 = (0.0D0,0.0D0)
      IF ( Z==z0 ) THEN
         c = z0
      ELSEIF ( w0<=2.5 ) THEN
         cr = Z
         c = cr
         wa0 = 0.0D0
         DO k = 1 , 80
            cr = -.5D0*cr*(4.0D0*k-3.0D0)/k/(2.0D0*k-1.0D0)             &
               & /(4.0D0*k+1.0D0)*zp2
            c = c + cr
            wa = ABS(c)
            IF ( DABS((wa-wa0)/wa)<eps .AND. k>10 ) GOTO 100
            wa0 = wa
         ENDDO
      ELSEIF ( w0>2.5 .AND. w0<4.5 ) THEN
         m = 85
         c = z0
         cf1 = z0
         cf0 = (1.0D-100,0.0D0)
         DO k = m , 0 , -1
            cf = (2.0D0*k+3.0D0)*cf0/zp - cf1
            IF ( k==INT(k/2)*2 ) c = c + cf
            cf1 = cf0
            cf0 = cf
         ENDDO
         c = 2.0D0/(pi*Z)*SIN(zp)/cf*c
      ELSE
!          See comment at CFS(), use C(z) = iC(-iz)
         IF ( DIMAG(Z)>-DBLE(Z) .AND. DIMAG(Z)<=DBLE(Z) ) THEN
!            right quadrant
            d = DCMPLX(.5D0,0.0D0)
         ELSEIF ( DIMAG(Z)>DBLE(Z) .AND. DIMAG(Z)>=-DBLE(Z) ) THEN
!            upper quadrant
            d = DCMPLX(0.0D0,.5D0)
         ELSEIF ( DIMAG(Z)<-DBLE(Z) .AND. DIMAG(Z)>=DBLE(Z) ) THEN
!            left quadrant
            d = DCMPLX(-.5D0,0.0D0)
         ELSE
!            lower quadrant
            d = DCMPLX(0.0D0,-.5D0)
         ENDIF
         cr = (1.0D0,0.0D0)
         cf = (1.0D0,0.0D0)
         DO k = 1 , 20
            cr = -.25D0*cr*(4.0D0*k-1.0D0)*(4.0D0*k-3.0D0)/zp2
            cf = cf + cr
         ENDDO
         cr = 1.0D0/(pi*Z*Z)
         cg = cr
         DO k = 1 , 12
            cr = -.25D0*cr*(4.0D0*k+1.0D0)*(4.0D0*k-1.0D0)/zp2
            cg = cg + cr
         ENDDO
         c = d + (cf*SIN(zp)-cg*COS(zp))/(pi*Z)
      ENDIF
 100  Zf = c
      Zd = COS(0.5*pi*Z*Z)
      END
 
 
 
!       **********************************
 
      SUBROUTINE FCS(X,C,S)
!
!       =================================================
!       Purpose: Compute Fresnel integrals C(x) and S(x)
!       Input :  x --- Argument of C(x) and S(x)
!       Output:  C --- C(x)
!                S --- S(x)
!       =================================================
!
      IMPLICIT NONE
!*--FCS12244
      DOUBLE PRECISION C , eps , f , f0 , f1 , g , pi , px , q , r , S ,&
                     & su , t , t0 , t2 , X , xa
      INTEGER k , m
      eps = 1.0D-15
      pi = 3.141592653589793D0
      xa = DABS(X)
      px = pi*xa
      t = .5D0*px*xa
      t2 = t*t
      IF ( xa==0.0 ) THEN
         C = 0.0D0
         S = 0.0D0
      ELSEIF ( xa<2.5D0 ) THEN
         r = xa
         C = r
         DO k = 1 , 50
            r = -.5D0*r*(4.0D0*k-3.0D0)/k/(2.0D0*k-1.0D0)               &
              & /(4.0D0*k+1.0D0)*t2
            C = C + r
            IF ( DABS(r)<DABS(C)*eps ) GOTO 50
         ENDDO
 50      S = xa*t/3.0D0
         r = S
         DO k = 1 , 50
            r = -.5D0*r*(4.0D0*k-1.0D0)/k/(2.0D0*k+1.0D0)               &
              & /(4.0D0*k+3.0D0)*t2
            S = S + r
            IF ( DABS(r)<DABS(S)*eps ) GOTO 100
         ENDDO
      ELSEIF ( xa<4.5D0 ) THEN
         m = INT(42.0+1.75*t)
         su = 0.0D0
         C = 0.0D0
         S = 0.0D0
         f1 = 0.0D0
         f0 = 1.0D-100
         DO k = m , 0 , -1
            f = (2.0D0*k+3.0D0)*f0/t - f1
            IF ( k==INT(k/2)*2 ) THEN
               C = C + f
            ELSE
               S = S + f
            ENDIF
            su = su + (2.0D0*k+1.0D0)*f*f
            f1 = f0
            f0 = f
         ENDDO
         q = DSQRT(su)
         C = C*xa/q
         S = S*xa/q
      ELSE
         r = 1.0D0
         f = 1.0D0
         DO k = 1 , 20
            r = -.25D0*r*(4.0D0*k-1.0D0)*(4.0D0*k-3.0D0)/t2
            f = f + r
         ENDDO
         r = 1.0D0/(px*xa)
         g = r
         DO k = 1 , 12
            r = -.25D0*r*(4.0D0*k+1.0D0)*(4.0D0*k-1.0D0)/t2
            g = g + r
         ENDDO
         t0 = t - INT(t/(2.0D0*pi))*2.0D0*pi
         C = .5D0 + (f*DSIN(t0)-g*DCOS(t0))/px
         S = .5D0 - (f*DCOS(t0)+g*DSIN(t0))/px
      ENDIF
 100  IF ( X<0.0D0 ) THEN
         C = -C
         S = -S
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE RCTJ(N,X,Nm,Rj,Dj)
!
!       ========================================================
!       Purpose: Compute Riccati-Bessel functions of the first
!                kind and their derivatives
!       Input:   x --- Argument of Riccati-Bessel function
!                n --- Order of jn(x)  ( n = 0,1,2,... )
!       Output:  RJ(n) --- x·jn(x)
!                DJ(n) --- [x·jn(x)]'
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ========================================================
!
      IMPLICIT NONE
!*--RCTJ12339
      DOUBLE PRECISION cs , Dj , f , f0 , f1 , Rj , rj0 , rj1 , X
      INTEGER k , m , MSTA1 , MSTA2 , N , Nm
      DIMENSION Rj(0:N) , Dj(0:N)
      Nm = N
      IF ( DABS(X)<1.0D-100 ) THEN
         DO k = 0 , N
            Rj(k) = 0.0D0
            Dj(k) = 0.0D0
         ENDDO
         Dj(0) = 1.0D0
         RETURN
      ENDIF
      Rj(0) = DSIN(X)
      Rj(1) = Rj(0)/X - DCOS(X)
      rj0 = Rj(0)
      rj1 = Rj(1)
      cs = 0.0D0
      f = 0.0D0
      IF ( N>=2 ) THEN
         m = MSTA1(X,200)
         IF ( m<N ) THEN
            Nm = m
         ELSE
            m = MSTA2(X,N,15)
         ENDIF
         f0 = 0.0D0
         f1 = 1.0D-100
         DO k = m , 0 , -1
            f = (2.0D0*k+3.0D0)*f1/X - f0
            IF ( k<=Nm ) Rj(k) = f
            f0 = f1
            f1 = f
         ENDDO
         IF ( DABS(rj0)>DABS(rj1) ) cs = rj0/f
         IF ( DABS(rj0)<=DABS(rj1) ) cs = rj1/f0
         DO k = 0 , Nm
            Rj(k) = cs*Rj(k)
         ENDDO
      ENDIF
      Dj(0) = DCOS(X)
      DO k = 1 , Nm
         Dj(k) = -k*Rj(k)/X + Rj(k-1)
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE HERZO(N,X,W)
!
!       ========================================================
!       Purpose : Compute the zeros of Hermite polynomial Ln(x)
!                 in the interval [-∞,∞], and the corresponding
!                 weighting coefficients for Gauss-Hermite
!                 integration
!       Input :   n    --- Order of the Hermite polynomial
!                 X(n) --- Zeros of the Hermite polynomial
!                 W(n) --- Corresponding weighting coefficients
!       ========================================================
!
      IMPLICIT NONE
!*--HERZO12405
      DOUBLE PRECISION f0 , f1 , fd , gd , hd , hf , hn , p , q , r ,   &
                     & r1 , r2 , W , wp , X , z , z0 , zl
      INTEGER i , it , j , k , N , nr
      DIMENSION X(N) , W(N)
      hn = 1.0D0/N
      zl = -1.1611D0 + 1.46D0*N**0.5
      z = 0.0D0
      hf = 0.0D0
      hd = 0.0D0
      DO nr = 1 , N/2
         IF ( nr==1 ) z = zl
         IF ( nr/=1 ) z = z - hn*(N/2+1-nr)
         it = 0
 50      it = it + 1
         z0 = z
         f0 = 1.0D0
         f1 = 2.0D0*z
         DO k = 2 , N
            hf = 2.0D0*z*f1 - 2.0D0*(k-1.0D0)*f0
            hd = 2.0D0*k*f1
            f0 = f1
            f1 = hf
         ENDDO
         p = 1.0D0
         DO i = 1 , nr - 1
            p = p*(z-X(i))
         ENDDO
         fd = hf/p
         q = 0.0D0
         DO i = 1 , nr - 1
            wp = 1.0D0
            DO j = 1 , nr - 1
               IF ( j/=i ) wp = wp*(z-X(j))
            ENDDO
            q = q + wp
         ENDDO
         gd = (hd-q*fd)/p
         z = z - fd/gd
         IF ( it<=40 .AND. DABS((z-z0)/z)>1.0D-15 ) GOTO 50
         X(nr) = z
         X(N+1-nr) = -z
         r = 1.0D0
         DO k = 1 , N
            r = 2.0D0*r*k
         ENDDO
         W(nr) = 3.544907701811D0*r/(hd*hd)
         W(N+1-nr) = W(nr)
      ENDDO
      IF ( N/=2*INT(N/2) ) THEN
         r1 = 1.0D0
         r2 = 1.0D0
         DO j = 1 , N
            r1 = 2.0D0*r1*j
            IF ( j>=(N+1)/2 ) r2 = r2*j
         ENDDO
         W(N/2+1) = 0.88622692545276D0*r1/(r2*r2)
         X(N/2+1) = 0.0D0
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE JY01B(X,Bj0,Dj0,Bj1,Dj1,By0,Dy0,By1,Dy1)
!
!       =======================================================
!       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
!                Y1(x), and their derivatives
!       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ≥ 0 )
!       Output:  BJ0 --- J0(x)
!                DJ0 --- J0'(x)
!                BJ1 --- J1(x)
!                DJ1 --- J1'(x)
!                BY0 --- Y0(x)
!                DY0 --- Y0'(x)
!                BY1 --- Y1(x)
!                DY1 --- Y1'(x)
!       =======================================================
!
      IMPLICIT NONE
!*--JY01B12488
      DOUBLE PRECISION a0 , Bj0 , Bj1 , By0 , By1 , Dj0 , Dj1 , Dy0 ,   &
                     & Dy1 , p0 , p1 , pi , q0 , q1 , t , t2 , ta0 ,    &
                     & ta1 , X
      pi = 3.141592653589793D0
      IF ( X==0.0D0 ) THEN
         Bj0 = 1.0D0
         Bj1 = 0.0D0
         Dj0 = 0.0D0
         Dj1 = 0.5D0
         By0 = -1.0D+300
         By1 = -1.0D+300
         Dy0 = 1.0D+300
         Dy1 = 1.0D+300
         RETURN
      ELSEIF ( X<=4.0D0 ) THEN
         t = X/4.0D0
         t2 = t*t
         Bj0 = ((((((-.5014415D-3*t2+.76771853D-2)*t2-.0709253492D0)*t2+&
             & .4443584263D0)*t2-1.7777560599D0)*t2+3.9999973021D0)     &
             & *t2-3.9999998721D0)*t2 + 1.0D0
         Bj1 = t*                                                       &
             & (((((((-.1289769D-3*t2+.22069155D-2)*t2-.0236616773D0)*t2&
             & +.1777582922D0)*t2-.8888839649D0)*t2+2.6666660544D0)     &
             & *t2-3.9999999710D0)*t2+1.9999999998D0)
         By0 = (((((((-.567433D-4*t2+.859977D-3)*t2-.94855882D-2)*t2+   &
             & .0772975809D0)*t2-.4261737419D0)*t2+1.4216421221D0)      &
             & *t2-2.3498519931D0)*t2+1.0766115157D0)*t2 + .3674669052D0
         By0 = 2.0D0/pi*DLOG(X/2.0D0)*Bj0 + By0
         By1 = ((((((((.6535773D-3*t2-.0108175626D0)*t2+.107657606D0)*t2&
             & -.7268945577D0)*t2+3.1261399273D0)*t2-7.3980241381D0)    &
             & *t2+6.8529236342D0)*t2+.3932562018D0)*t2-.6366197726D0)/X
         By1 = 2.0D0/pi*DLOG(X/2.0D0)*Bj1 + By1
      ELSE
         t = 4.0D0/X
         t2 = t*t
         a0 = DSQRT(2.0D0/(pi*X))
         p0 = ((((-.9285D-5*t2+.43506D-4)*t2-.122226D-3)*t2+.434725D-3) &
            & *t2-.4394275D-2)*t2 + .999999997D0
         q0 = t*(((((.8099D-5*t2-.35614D-4)*t2+.85844D-4)*t2-.218024D-3)&
            & *t2+.1144106D-2)*t2-.031249995D0)
         ta0 = X - .25D0*pi
         Bj0 = a0*(p0*DCOS(ta0)-q0*DSIN(ta0))
         By0 = a0*(p0*DSIN(ta0)+q0*DCOS(ta0))
         p1 = ((((.10632D-4*t2-.50363D-4)*t2+.145575D-3)*t2-.559487D-3) &
            & *t2+.7323931D-2)*t2 + 1.000000004D0
         q1 = t*                                                        &
            & (((((-.9173D-5*t2+.40658D-4)*t2-.99941D-4)*t2+.266891D-3) &
            & *t2-.1601836D-2)*t2+.093749994D0)
         ta1 = X - .75D0*pi
         Bj1 = a0*(p1*DCOS(ta1)-q1*DSIN(ta1))
         By1 = a0*(p1*DSIN(ta1)+q1*DCOS(ta1))
      ENDIF
      Dj0 = -Bj1
      Dj1 = Bj0 - Bj1/X
      Dy0 = -By1
      Dy1 = By0 - By1/X
      END
 
!       **********************************
 
      SUBROUTINE ENXB(N,X,En)
!
!       ===============================================
!       Purpose: Compute exponential integral En(x)
!       Input :  x --- Argument of En(x)
!                n --- Order of En(x)  (n = 0,1,2,...)
!       Output:  EN(n) --- En(x)
!       ===============================================
!
      IMPLICIT NONE
!*--ENXB12562
      DOUBLE PRECISION En , ens , ps , r , rp , s , s0 , t , t0 , X
      INTEGER j , k , l , m , N
      DIMENSION En(0:N)
      IF ( X==0.0 ) THEN
         En(0) = 1.0D+300
         En(1) = 1.0D+300
         DO k = 2 , N
            En(k) = 1.0D0/(k-1.0)
         ENDDO
         RETURN
      ELSEIF ( X<=1.0 ) THEN
         En(0) = EXP(-X)/X
         s0 = 0.0D0
         DO l = 1 , N
            rp = 1.0D0
            DO j = 1 , l - 1
               rp = -rp*X/j
            ENDDO
            ps = -0.5772156649015328D0
            DO m = 1 , l - 1
               ps = ps + 1.0D0/m
            ENDDO
            ens = rp*(-DLOG(X)+ps)
            s = 0.0D0
            DO m = 0 , 20
               IF ( m/=l-1 ) THEN
                  r = 1.0D0
                  DO j = 1 , m
                     r = -r*X/j
                  ENDDO
                  s = s + r/(m-l+1.0D0)
                  IF ( DABS(s-s0)<DABS(s)*1.0D-15 ) GOTO 20
                  s0 = s
               ENDIF
            ENDDO
 20         En(l) = ens - s
         ENDDO
      ELSE
         En(0) = EXP(-X)/X
         m = 15 + INT(100.0/X)
         DO l = 1 , N
            t0 = 0.0D0
            DO k = m , 1 , -1
               t0 = (l+k-1.0D0)/(1.0D0+k/(X+t0))
            ENDDO
            t = 1.0D0/(X+t0)
            En(l) = EXP(-X)*t
         ENDDO
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE SPHK(N,X,Nm,Sk,Dk)
!
!       =====================================================
!       Purpose: Compute modified spherical Bessel functions
!                of the second kind, kn(x) and kn'(x)
!       Input :  x --- Argument of kn(x)  ( x ≥ 0 )
!                n --- Order of kn(x) ( n = 0,1,2,... )
!       Output:  SK(n) --- kn(x)
!                DK(n) --- kn'(x)
!                NM --- Highest order computed
!       =====================================================
!
      IMPLICIT NONE
!*--SPHK12632
      DOUBLE PRECISION Dk , f , f0 , f1 , pi , Sk , X
      INTEGER k , N , Nm
      DIMENSION Sk(0:N) , Dk(0:N)
      pi = 3.141592653589793D0
      Nm = N
      IF ( X<1.0D-60 ) THEN
         DO k = 0 , N
            Sk(k) = 1.0D+300
            Dk(k) = -1.0D+300
         ENDDO
         RETURN
      ENDIF
      Sk(0) = 0.5D0*pi/X*EXP(-X)
      Sk(1) = Sk(0)*(1.0D0+1.0D0/X)
      f0 = Sk(0)
      f1 = Sk(1)
      DO k = 2 , N
         f = (2.0D0*k-1.0D0)*f1/X + f0
         Sk(k) = f
         IF ( DABS(f)>1.0D+300 ) GOTO 100
         f0 = f1
         f1 = f
      ENDDO
 100  Nm = k - 1
      Dk(0) = -Sk(1)
      DO k = 1 , Nm
         Dk(k) = -Sk(k-1) - (k+1.0D0)/X*Sk(k)
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE ENXA(N,X,En)
!
!       ============================================
!       Purpose: Compute exponential integral En(x)
!       Input :  x --- Argument of En(x) ( x ≤ 20 )
!                n --- Order of En(x)
!       Output:  EN(n) --- En(x)
!       Routine called: E1XB for computing E1(x)
!       ============================================
!
      IMPLICIT NONE
!*--ENXA12679
      DOUBLE PRECISION e1 , ek , En , X
      INTEGER k , N
      DIMENSION En(0:N)
      En(0) = EXP(-X)/X
      CALL E1XB(X,e1)
      En(1) = e1
      DO k = 2 , N
         ek = (EXP(-X)-X*e1)/(k-1.0D0)
         En(k) = ek
         e1 = ek
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE GAIH(X,Ga)
!
!       =====================================================
!       Purpose: Compute gamma function Г(x)
!       Input :  x  --- Argument of Г(x), x = n/2, n=1,2,…
!       Output:  GA --- Г(x)
!       =====================================================
!
      IMPLICIT NONE
!*--GAIH12709
      DOUBLE PRECISION Ga , pi , X
      INTEGER k , m , m1
      pi = 3.141592653589793D0
      IF ( X==INT(X) .AND. X>0.0 ) THEN
         Ga = 1.0D0
         m1 = INT(X-1.0)
         DO k = 2 , m1
            Ga = Ga*k
         ENDDO
      ELSEIF ( X+.5D0==INT(X+.5D0) .AND. X>0.0 ) THEN
         m = INT(X)
         Ga = DSQRT(pi)
         DO k = 1 , m
            Ga = 0.5D0*Ga*(2.0D0*k-1.0D0)
         ENDDO
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE PBVV(V,X,Vv,Vp,Pvf,Pvd)
!
!       ===================================================
!       Purpose: Compute parabolic cylinder functions Vv(x)
!                and their derivatives
!       Input:   x --- Argument of Vv(x)
!                v --- Order of Vv(x)
!       Output:  VV(na) --- Vv(x)
!                VP(na) --- Vv'(x)
!                ( na = |n|, v = n+v0, |v0| < 1
!                  n = 0,±1,±2,… )
!                PVF --- Vv(x)
!                PVD --- Vv'(x)
!       Routines called:
!             (1) VVSA for computing Vv(x) for small |x|
!             (2) VVLA for computing Vv(x) for large |x|
!       ===================================================
!
      IMPLICIT NONE
!*--PBVV12752
      DOUBLE PRECISION f , f0 , f1 , pi , pv0 , Pvd , Pvf , q2p , qe ,  &
                     & s0 , V , v0 , v1 , v2 , vh , Vp , Vv , X , xa
      INTEGER ja , k , kv , l , m , na , nv
      DIMENSION Vv(0:*) , Vp(0:*)
      pi = 3.141592653589793D0
      xa = DABS(X)
      vh = V
      V = V + DSIGN(1.0D0,V)
      nv = INT(V)
      v0 = V - nv
      na = ABS(nv)
      qe = EXP(0.25D0*X*X)
      q2p = DSQRT(2.0D0/pi)
      ja = 0
      IF ( na>=1 ) ja = 1
      f = 0.0D0
      IF ( V<=0.0 ) THEN
         IF ( v0==0.0 ) THEN
            IF ( xa<=7.5 ) CALL VVSA(v0,X,pv0)
            IF ( xa>7.5 ) CALL VVLA(v0,X,pv0)
            f0 = q2p*qe
            f1 = X*f0
            Vv(0) = pv0
            Vv(1) = f0
            Vv(2) = f1
         ELSE
            DO l = 0 , ja
               v1 = v0 - l
               IF ( xa<=7.5 ) CALL VVSA(v1,X,f1)
               IF ( xa>7.5 ) CALL VVLA(v1,X,f1)
               IF ( l==0 ) f0 = f1
            ENDDO
            Vv(0) = f0
            Vv(1) = f1
         ENDIF
         kv = 2
         IF ( v0==0.0 ) kv = 3
         DO k = kv , na
            f = X*f1 + (k-v0-2.0D0)*f0
            Vv(k) = f
            f0 = f1
            f1 = f
         ENDDO
      ELSEIF ( X>=0.0 .AND. X<=7.5D0 ) THEN
         v2 = V
         IF ( v2<1.0 ) v2 = v2 + 1.0D0
         CALL VVSA(v2,X,f1)
         v1 = v2 - 1.0D0
         kv = INT(v2)
         CALL VVSA(v1,X,f0)
         Vv(kv) = f1
         Vv(kv-1) = f0
         DO k = kv - 2 , 0 , -1
            f = X*f0 - (k+v0+2.0D0)*f1
            IF ( k<=na ) Vv(k) = f
            f1 = f0
            f0 = f
         ENDDO
      ELSEIF ( X>7.5D0 ) THEN
         CALL VVLA(v0,X,pv0)
         m = 100 + ABS(na)
         Vv(1) = pv0
         f1 = 0.0D0
         f0 = 1.0D-40
         DO k = m , 0 , -1
            f = X*f0 - (k+v0+2.0D0)*f1
            IF ( k<=na ) Vv(k) = f
            f1 = f0
            f0 = f
         ENDDO
         s0 = pv0/f
         DO k = 0 , na
            Vv(k) = s0*Vv(k)
         ENDDO
      ELSE
         IF ( xa<=7.5D0 ) THEN
            CALL VVSA(v0,X,f0)
            v1 = v0 + 1.0
            CALL VVSA(v1,X,f1)
         ELSE
            CALL VVLA(v0,X,f0)
            v1 = v0 + 1.0D0
            CALL VVLA(v1,X,f1)
         ENDIF
         Vv(0) = f0
         Vv(1) = f1
         DO k = 2 , na
            f = (X*f1-f0)/(k+v0)
            Vv(k) = f
            f0 = f1
            f1 = f
         ENDDO
      ENDIF
      DO k = 0 , na - 1
         v1 = v0 + k
         IF ( V>=0.0D0 ) THEN
            Vp(k) = 0.5D0*X*Vv(k) - (v1+1.0D0)*Vv(k+1)
         ELSE
            Vp(k) = -0.5D0*X*Vv(k) + Vv(k+1)
         ENDIF
      ENDDO
      Pvf = Vv(na-1)
      Pvd = Vp(na-1)
      V = vh
      END
 
 
 
!       **********************************
 
      SUBROUTINE CLQMN(Mm,M,N,X,Y,Cqm,Cqd)
!
!       =======================================================
!       Purpose: Compute the associated Legendre functions of
!                the second kind, Qmn(z) and Qmn'(z), for a
!                complex argument
!       Input :  x  --- Real part of z
!                y  --- Imaginary part of z
!                m  --- Order of Qmn(z)  ( m = 0,1,2,… )
!                n  --- Degree of Qmn(z) ( n = 0,1,2,… )
!                mm --- Physical dimension of CQM and CQD
!       Output:  CQM(m,n) --- Qmn(z)
!                CQD(m,n) --- Qmn'(z)
!       =======================================================
!
      IMPLICIT NONE
!*--CLQMN12882
      COMPLEX*16 cq0 , cq1 , cq10 , Cqd , cqf , cqf0 , cqf1 , cqf2 ,    &
               & Cqm , z , zq , zs
      INTEGER i , j , k , km , ls , M , Mm , N
      DOUBLE PRECISION X , xc , Y
      DIMENSION Cqm(0:Mm,0:N) , Cqd(0:Mm,0:N)
      z = DCMPLX(X,Y)
      IF ( DABS(X)==1.0D0 .AND. Y==0.0D0 ) THEN
         DO i = 0 , M
            DO j = 0 , N
               Cqm(i,j) = (1.0D+300,0.0D0)
               Cqd(i,j) = (1.0D+300,0.0D0)
            ENDDO
         ENDDO
         RETURN
      ENDIF
      xc = ABS(z)
      ls = 0
      IF ( DIMAG(z)==0.0D0 .OR. xc<1.0D0 ) ls = 1
      IF ( xc>1.0D0 ) ls = -1
      zq = SQRT(ls*(1.0D0-z*z))
      zs = ls*(1.0D0-z*z)
      cq0 = 0.5D0*LOG(ls*(1.0D0+z)/(1.0D0-z))
      IF ( xc<1.0001D0 ) THEN
         Cqm(0,0) = cq0
         Cqm(0,1) = z*cq0 - 1.0D0
         Cqm(1,0) = -1.0D0/zq
         Cqm(1,1) = -zq*(cq0+z/(1.0D0-z*z))
         DO i = 0 , 1
            DO j = 2 , N
               Cqm(i,j) = ((2.0D0*j-1.0D0)*z*Cqm(i,j-1)-(j+i-1.0D0)     &
                        & *Cqm(i,j-2))/(j-i)
            ENDDO
         ENDDO
         DO j = 0 , N
            DO i = 2 , M
               Cqm(i,j) = -2.0D0*(i-1.0D0)*z/zq*Cqm(i-1,j)              &
                        & - ls*(j+i-1.0D0)*(j-i+2.0D0)*Cqm(i-2,j)
            ENDDO
         ENDDO
      ELSE
         IF ( xc>1.1 ) THEN
            km = 40 + M + N
         ELSE
            km = (40+M+N)*INT(-1.0-1.8*LOG(xc-1.0))
         ENDIF
         cqf2 = (0.0D0,0.0D0)
         cqf1 = (1.0D0,0.0D0)
         DO k = km , 0 , -1
            cqf0 = ((2*k+3.0D0)*z*cqf1-(k+2.0D0)*cqf2)/(k+1.0D0)
            IF ( k<=N ) Cqm(0,k) = cqf0
            cqf2 = cqf1
            cqf1 = cqf0
         ENDDO
         DO k = 0 , N
            Cqm(0,k) = cq0*Cqm(0,k)/cqf0
         ENDDO
         cqf2 = 0.0D0
         cqf1 = 1.0D0
         DO k = km , 0 , -1
            cqf0 = ((2*k+3.0D0)*z*cqf1-(k+1.0D0)*cqf2)/(k+2.0D0)
            IF ( k<=N ) Cqm(1,k) = cqf0
            cqf2 = cqf1
            cqf1 = cqf0
         ENDDO
         cq10 = -1.0D0/zq
         DO k = 0 , N
            Cqm(1,k) = cq10*Cqm(1,k)/cqf0
         ENDDO
         DO j = 0 , N
            cq0 = Cqm(0,j)
            cq1 = Cqm(1,j)
            DO i = 0 , M - 2
               cqf = -2.0D0*(i+1)*z/zq*cq1 + (j-i)*(j+i+1.0D0)*cq0
               Cqm(i+2,j) = cqf
               cq0 = cq1
               cq1 = cqf
            ENDDO
         ENDDO
      ENDIF
      Cqd(0,0) = ls/zs
      DO j = 1 , N
         Cqd(0,j) = ls*j*(Cqm(0,j-1)-z*Cqm(0,j))/zs
      ENDDO
      DO j = 0 , N
         DO i = 1 , M
            Cqd(i,j) = ls*i*z/zs*Cqm(i,j) + (i+j)*(j-i+1.0D0)           &
                     & /zq*Cqm(i-1,j)
         ENDDO
      ENDDO
      END
 
 
!       **********************************
 
      SUBROUTINE SEGV(M,N,C,Kd,Cv,Eg)
!
!       =========================================================
!       Purpose: Compute the characteristic values of spheroidal
!                wave functions
!       Input :  m  --- Mode parameter
!                n  --- Mode parameter
!                c  --- Spheroidal parameter
!                KD --- Function code
!                       KD=1 for Prolate; KD=-1 for Oblate
!       Output:  CV --- Characteristic value for given m, n and c
!                EG(L) --- Characteristic value for mode m and n'
!                          ( L = n' - m + 1 )
!       =========================================================
!
      IMPLICIT NONE
!*--SEGV12996
      DOUBLE PRECISION a , b , C , cs , Cv , cv0 , d , d2k , dk0 , dk1 ,&
                     & dk2 , e , Eg , f , g , h , s , t , t1 , x1
      DOUBLE PRECISION xa , xb
      INTEGER i , icm , j , k , k1 , Kd , l , M , N , nm , nm1
      DIMENSION b(100) , h(100) , d(300) , e(300) , f(300) , cv0(100) , &
              & a(300) , g(300) , Eg(200)
      IF ( C<1.0D-10 ) THEN
         DO i = 1 , N - M + 1
            Eg(i) = (i+M)*(i+M-1.0D0)
         ENDDO
         GOTO 100
      ENDIF
      icm = (N-M+2)/2
      nm = 10 + INT(0.5*(N-M)+C)
      cs = C*C*Kd
      k = 0
      DO l = 0 , 1
         DO i = 1 , nm
            IF ( l==0 ) k = 2*(i-1)
            IF ( l==1 ) k = 2*i - 1
            dk0 = M + k
            dk1 = M + k + 1
            dk2 = 2*(M+k)
            d2k = 2*M + k
            a(i) = (d2k+2.0)*(d2k+1.0)/((dk2+3.0)*(dk2+5.0))*cs
            d(i) = dk0*dk1 + (2.0*dk0*dk1-2.0*M*M-1.0)                  &
                 & /((dk2-1.0)*(dk2+3.0))*cs
            g(i) = k*(k-1.0)/((dk2-3.0)*(dk2-1.0))*cs
         ENDDO
         DO k = 2 , nm
            e(k) = DSQRT(a(k-1)*g(k))
            f(k) = e(k)*e(k)
         ENDDO
         f(1) = 0.0D0
         e(1) = 0.0D0
         xa = d(nm) + DABS(e(nm))
         xb = d(nm) - DABS(e(nm))
         nm1 = nm - 1
         DO i = 1 , nm1
            t = DABS(e(i)) + DABS(e(i+1))
            t1 = d(i) + t
            IF ( xa<t1 ) xa = t1
            t1 = d(i) - t
            IF ( t1<xb ) xb = t1
         ENDDO
         DO i = 1 , icm
            b(i) = xa
            h(i) = xb
         ENDDO
         DO k = 1 , icm
            DO k1 = k , icm
               IF ( b(k1)<b(k) ) THEN
                  b(k) = b(k1)
                  GOTO 20
               ENDIF
            ENDDO
 20         IF ( k/=1 ) THEN
               IF ( h(k)<h(k-1) ) h(k) = h(k-1)
            ENDIF
 40         x1 = (b(k)+h(k))/2.0D0
            cv0(k) = x1
            IF ( DABS((b(k)-h(k))/x1)<1.0D-14 ) THEN
               cv0(k) = x1
               IF ( l==0 ) Eg(2*k-1) = cv0(k)
               IF ( l==1 ) Eg(2*k) = cv0(k)
            ELSE
               j = 0
               s = 1.0D0
               DO i = 1 , nm
                  IF ( s==0.0D0 ) s = s + 1.0D-30
                  t = f(i)/s
                  s = d(i) - t - x1
                  IF ( s<0.0D0 ) j = j + 1
               ENDDO
               IF ( j<k ) THEN
                  h(k) = x1
               ELSE
                  b(k) = x1
                  IF ( j>=icm ) THEN
                     b(icm) = x1
                  ELSE
                     IF ( h(j+1)<x1 ) h(j+1) = x1
                     IF ( x1<b(j) ) b(j) = x1
                  ENDIF
               ENDIF
               GOTO 40
            ENDIF
         ENDDO
      ENDDO
 100  Cv = Eg(N-M+1)
      END
 
 
!       **********************************
 
      SUBROUTINE CIKNB(N,Z,Nm,Cbi,Cdi,Cbk,Cdk)
!
!       ============================================================
!       Purpose: Compute modified Bessel functions In(z) and Kn(z),
!                and their derivatives for a complex argument
!       Input:   z --- Complex argument
!                n --- Order of In(z) and Kn(z)
!       Output:  CBI(n) --- In(z)
!                CDI(n) --- In'(z)
!                CBK(n) --- Kn(z)
!                CDK(n) --- Kn'(z)
!                NM --- Highest order computed
!       Routones called:
!                MSTA1 and MSTA2 to compute the starting point for
!                backward recurrence
!       ===========================================================
!
      IMPLICIT NONE
!*--CIKNB13113
      DOUBLE PRECISION a0 , el , fac , pi , vt
      COMPLEX*16 ca0 , Cbi , Cbk , cbkl , cbs , Cdi , Cdk , cf , cf0 ,  &
               & cf1 , cg , cg0 , cg1 , ci , cr , cs0 , csk0 , Z , z1
      INTEGER k , k0 , l , m , MSTA1 , MSTA2 , N , Nm
      DIMENSION Cbi(0:N) , Cdi(0:N) , Cbk(0:N) , Cdk(0:N)
      pi = 3.141592653589793D0
      el = 0.57721566490153D0
      a0 = ABS(Z)
      Nm = N
      IF ( a0<1.0D-100 ) THEN
         DO k = 0 , N
            Cbi(k) = (0.0D0,0.0D0)
            Cbk(k) = (1.0D+300,0.0D0)
            Cdi(k) = (0.0D0,0.0D0)
            Cdk(k) = -(1.0D+300,0.0D0)
         ENDDO
         Cbi(0) = (1.0D0,0.0D0)
         Cdi(1) = (0.5D0,0.0D0)
         RETURN
      ENDIF
      z1 = Z
      ci = (0.0D0,1.0D0)
      IF ( DBLE(Z)<0.0 ) z1 = -Z
      IF ( N==0 ) Nm = 1
      m = MSTA1(a0,200)
      IF ( m<Nm ) THEN
         Nm = m
      ELSE
         m = MSTA2(a0,Nm,15)
      ENDIF
      cbs = 0.0D0
      csk0 = 0.0D0
      cf0 = 0.0D0
      cf1 = 1.0D-100
      DO k = m , 0 , -1
         cf = 2.0D0*(k+1.0D0)*cf1/z1 + cf0
         IF ( k<=Nm ) Cbi(k) = cf
         IF ( k/=0 .AND. k==2*INT(k/2) ) csk0 = csk0 + 4.0D0*cf/k
         cbs = cbs + 2.0D0*cf
         cf0 = cf1
         cf1 = cf
      ENDDO
      cs0 = EXP(z1)/(cbs-cf)
      DO k = 0 , Nm
         Cbi(k) = cs0*Cbi(k)
      ENDDO
      IF ( a0<=9.0 ) THEN
         Cbk(0) = -(LOG(0.5D0*z1)+el)*Cbi(0) + cs0*csk0
         Cbk(1) = (1.0D0/z1-Cbi(1)*Cbk(0))/Cbi(0)
      ELSE
         ca0 = SQRT(pi/(2.0D0*z1))*EXP(-z1)
         k0 = 16
         IF ( a0>=25.0 ) k0 = 10
         IF ( a0>=80.0 ) k0 = 8
         IF ( a0>=200.0 ) k0 = 6
         DO l = 0 , 1
            cbkl = 1.0D0
            vt = 4.0D0*l
            cr = (1.0D0,0.0D0)
            DO k = 1 , k0
               cr = 0.125D0*cr*(vt-(2.0*k-1.0)**2)/(k*z1)
               cbkl = cbkl + cr
            ENDDO
            Cbk(l) = ca0*cbkl
         ENDDO
      ENDIF
      cg0 = Cbk(0)
      cg1 = Cbk(1)
      DO k = 2 , Nm
         cg = 2.0D0*(k-1.0D0)/z1*cg1 + cg0
         Cbk(k) = cg
         cg0 = cg1
         cg1 = cg
      ENDDO
      IF ( DBLE(Z)<0.0 ) THEN
         fac = 1.0D0
         DO k = 0 , Nm
            IF ( DIMAG(Z)<0.0 ) THEN
               Cbk(k) = fac*Cbk(k) + ci*pi*Cbi(k)
            ELSE
               Cbk(k) = fac*Cbk(k) - ci*pi*Cbi(k)
            ENDIF
            Cbi(k) = fac*Cbi(k)
            fac = -fac
         ENDDO
      ENDIF
      Cdi(0) = Cbi(1)
      Cdk(0) = -Cbk(1)
      DO k = 1 , Nm
         Cdi(k) = Cbi(k-1) - k/Z*Cbi(k)
         Cdk(k) = -Cbk(k-1) - k/Z*Cbk(k)
      ENDDO
      END
 
 
!       **********************************
 
      SUBROUTINE CIKNA(N,Z,Nm,Cbi,Cdi,Cbk,Cdk)
!
!       ========================================================
!       Purpose: Compute modified Bessel functions In(z), Kn(x)
!                and their derivatives for a complex argument
!       Input :  z --- Complex argument of In(z) and Kn(z)
!                n --- Order of In(z) and Kn(z)
!       Output:  CBI(n) --- In(z)
!                CDI(n) --- In'(z)
!                CBK(n) --- Kn(z)
!                CDK(n) --- Kn'(z)
!                NM --- Highest order computed
!       Routines called:
!             (1) CIK01 to compute I0(z), I1(z) K0(z) & K1(z)
!             (2) MSTA1 and MSTA2 to compute the starting
!                 point for backward recurrence
!       ========================================================
!
      IMPLICIT NONE
!*--CIKNA13233
      DOUBLE PRECISION a0
      COMPLEX*16 Cbi , cbi0 , cbi1 , Cbk , cbk0 , cbk1 , Cdi , cdi0 ,   &
               & cdi1 , Cdk , cdk0 , cdk1 , cf , cf1 , cf2 , ckk , cs , &
               & Z
      INTEGER k , m , MSTA1 , MSTA2 , N , Nm
      DIMENSION Cbi(0:N) , Cdi(0:N) , Cbk(0:N) , Cdk(0:N)
      a0 = ABS(Z)
      Nm = N
      IF ( a0<1.0D-100 ) THEN
         DO k = 0 , N
            Cbi(k) = (0.0D0,0.0D0)
            Cdi(k) = (0.0D0,0.0D0)
            Cbk(k) = -(1.0D+300,0.0D0)
            Cdk(k) = (1.0D+300,0.0D0)
         ENDDO
         Cbi(0) = (1.0D0,0.0D0)
         Cdi(1) = (0.5D0,0.0D0)
         RETURN
      ENDIF
      CALL CIK01(Z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1)
      Cbi(0) = cbi0
      Cbi(1) = cbi1
      Cbk(0) = cbk0
      Cbk(1) = cbk1
      Cdi(0) = cdi0
      Cdi(1) = cdi1
      Cdk(0) = cdk0
      Cdk(1) = cdk1
      IF ( N<=1 ) RETURN
      m = MSTA1(a0,200)
      IF ( m<N ) THEN
         Nm = m
      ELSE
         m = MSTA2(a0,N,15)
      ENDIF
      cf2 = (0.0D0,0.0D0)
      cf1 = (1.0D-100,0.0D0)
      DO k = m , 0 , -1
         cf = 2.0D0*(k+1.0D0)/Z*cf1 + cf2
         IF ( k<=Nm ) Cbi(k) = cf
         cf2 = cf1
         cf1 = cf
      ENDDO
      cs = cbi0/cf
      DO k = 0 , Nm
         Cbi(k) = cs*Cbi(k)
      ENDDO
      DO k = 2 , Nm
         IF ( ABS(Cbi(k-1))>ABS(Cbi(k-2)) ) THEN
            ckk = (1.0D0/Z-Cbi(k)*Cbk(k-1))/Cbi(k-1)
         ELSE
            ckk = (Cbi(k)*Cbk(k-2)+2.0D0*(k-1.0D0)/(Z*Z))/Cbi(k-2)
         ENDIF
         Cbk(k) = ckk
      ENDDO
      DO k = 2 , Nm
         Cdi(k) = Cbi(k-1) - k/Z*Cbi(k)
         Cdk(k) = -Cbk(k-1) - k/Z*Cbk(k)
      ENDDO
      END
 
 
 
!       **********************************
 
      SUBROUTINE MTU12(Kf,Kc,M,Q,X,F1r,D1r,F2r,D2r)
!
!       ==============================================================
!       Purpose: Compute modified Mathieu functions of the first and
!                second kinds, Mcm(1)(2)(x,q) and Msm(1)(2)(x,q),
!                and their derivatives
!       Input:   KF --- Function code
!                       KF=1 for computing Mcm(x,q)
!                       KF=2 for computing Msm(x,q)
!                KC --- Function Code
!                       KC=1 for computing the first kind
!                       KC=2 for computing the second kind
!                            or Msm(2)(x,q) and Msm(2)'(x,q)
!                       KC=3 for computing both the first
!                            and second kinds
!                m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions ( q ≥ 0 )
!                x  --- Argument of Mathieu functions
!       Output:  F1R --- Mcm(1)(x,q) or Msm(1)(x,q)
!                D1R --- Derivative of Mcm(1)(x,q) or Msm(1)(x,q)
!                F2R --- Mcm(2)(x,q) or Msm(2)(x,q)
!                D2R --- Derivative of Mcm(2)(x,q) or Msm(2)(x,q)
!       Routines called:
!            (1) CVA2 for computing the characteristic values
!            (2) FCOEF for computing expansion coefficients
!            (3) JYNB for computing Jn(x), Yn(x) and their
!                derivatives
!       ==============================================================
!
      IMPLICIT NONE
!*--MTU1213332
      DOUBLE PRECISION a , bj1 , bj2 , by1 , by2 , c1 , c2 , D1r , D2r ,&
                     & dj1 , dj2 , DNAN , dy1 , dy2 , eps , F1r , F2r , &
                     & fg , Q , qm
      DOUBLE PRECISION u1 , u2 , w1 , w2 , X
      INTEGER ic , k , Kc , kd , Kf , km , M , nm
      DIMENSION fg(251) , bj1(0:251) , dj1(0:251) , bj2(0:251) ,        &
              & dj2(0:251) , by1(0:251) , dy1(0:251) , by2(0:251) ,     &
              & dy2(0:251)
      eps = 1.0D-14
      IF ( Kf==1 .AND. M==2*INT(M/2) ) kd = 1
      IF ( Kf==1 .AND. M/=2*INT(M/2) ) kd = 2
      IF ( Kf==2 .AND. M/=2*INT(M/2) ) kd = 3
      IF ( Kf==2 .AND. M==2*INT(M/2) ) kd = 4
      CALL CVA2(kd,M,Q,a)
      IF ( Q<=1.0D0 ) THEN
         qm = 7.5 + 56.1*SQRT(Q) - 134.7*Q + 90.7*SQRT(Q)*Q
      ELSE
         qm = 17.0 + 3.1*SQRT(Q) - .126*Q + .0037*SQRT(Q)*Q
      ENDIF
      km = INT(qm+0.5*M)
      IF ( km>=251 ) THEN
         F1r = DNAN()
         D1r = DNAN()
         F2r = DNAN()
         D2r = DNAN()
         RETURN
      ENDIF
      CALL FCOEF(kd,M,Q,a,fg)
      ic = INT(M/2) + 1
      IF ( kd==4 ) ic = M/2
      c1 = EXP(-X)
      c2 = EXP(X)
      u1 = DSQRT(Q)*c1
      u2 = DSQRT(Q)*c2
      CALL JYNB(km+1,u1,nm,bj1,dj1,by1,dy1)
      CALL JYNB(km+1,u2,nm,bj2,dj2,by2,dy2)
      w1 = 0.0D0
      w2 = 0.0D0
      IF ( Kc/=2 ) THEN
         F1r = 0.0D0
         DO k = 1 , km
            IF ( kd==1 ) THEN
               F1r = F1r + (-1)**(ic+k)*fg(k)*bj1(k-1)*bj2(k-1)
            ELSEIF ( kd==2 .OR. kd==3 ) THEN
               F1r = F1r + (-1)**(ic+k)*fg(k)                           &
                   & *(bj1(k-1)*bj2(k)+(-1)**kd*bj1(k)*bj2(k-1))
            ELSE
               F1r = F1r + (-1)**(ic+k)*fg(k)                           &
                   & *(bj1(k-1)*bj2(k+1)-bj1(k+1)*bj2(k-1))
            ENDIF
            IF ( k>=5 .AND. DABS(F1r-w1)<DABS(F1r)*eps ) GOTO 50
            w1 = F1r
         ENDDO
 50      F1r = F1r/fg(1)
         D1r = 0.0D0
         DO k = 1 , km
            IF ( kd==1 ) THEN
               D1r = D1r + (-1)**(ic+k)*fg(k)                           &
                   & *(c2*bj1(k-1)*dj2(k-1)-c1*dj1(k-1)*bj2(k-1))
            ELSEIF ( kd==2 .OR. kd==3 ) THEN
               D1r = D1r + (-1)**(ic+k)*fg(k)                           &
                   & *(c2*(bj1(k-1)*dj2(k)+(-1)**kd*bj1(k)*dj2(k-1))    &
                   & -c1*(dj1(k-1)*bj2(k)+(-1)**kd*dj1(k)*bj2(k-1)))
            ELSE
               D1r = D1r + (-1)**(ic+k)*fg(k)                           &
                   & *(c2*(bj1(k-1)*dj2(k+1)-bj1(k+1)*dj2(k-1))         &
                   & -c1*(dj1(k-1)*bj2(k+1)-dj1(k+1)*bj2(k-1)))
            ENDIF
            IF ( k>=5 .AND. DABS(D1r-w2)<DABS(D1r)*eps ) GOTO 100
            w2 = D1r
         ENDDO
 100     D1r = D1r*DSQRT(Q)/fg(1)
         IF ( Kc==1 ) RETURN
      ENDIF
      F2r = 0.0D0
      DO k = 1 , km
         IF ( kd==1 ) THEN
            F2r = F2r + (-1)**(ic+k)*fg(k)*bj1(k-1)*by2(k-1)
         ELSEIF ( kd==2 .OR. kd==3 ) THEN
            F2r = F2r + (-1)**(ic+k)*fg(k)                              &
                & *(bj1(k-1)*by2(k)+(-1)**kd*bj1(k)*by2(k-1))
         ELSE
            F2r = F2r + (-1)**(ic+k)*fg(k)                              &
                & *(bj1(k-1)*by2(k+1)-bj1(k+1)*by2(k-1))
         ENDIF
         IF ( k>=5 .AND. DABS(F2r-w1)<DABS(F2r)*eps ) GOTO 200
         w1 = F2r
      ENDDO
 200  F2r = F2r/fg(1)
      D2r = 0.0D0
      DO k = 1 , km
         IF ( kd==1 ) THEN
            D2r = D2r + (-1)**(ic+k)*fg(k)                              &
                & *(c2*bj1(k-1)*dy2(k-1)-c1*dj1(k-1)*by2(k-1))
         ELSEIF ( kd==2 .OR. kd==3 ) THEN
            D2r = D2r + (-1)**(ic+k)*fg(k)                              &
                & *(c2*(bj1(k-1)*dy2(k)+(-1)**kd*bj1(k)*dy2(k-1))       &
                & -c1*(dj1(k-1)*by2(k)+(-1)**kd*dj1(k)*by2(k-1)))
         ELSE
            D2r = D2r + (-1)**(ic+k)*fg(k)                              &
                & *(c2*(bj1(k-1)*dy2(k+1)-bj1(k+1)*dy2(k-1))            &
                & -c1*(dj1(k-1)*by2(k+1)-dj1(k+1)*by2(k-1)))
         ENDIF
         IF ( k>=5 .AND. DABS(D2r-w2)<DABS(D2r)*eps ) GOTO 300
         w2 = D2r
      ENDDO
 300  D2r = D2r*DSQRT(Q)/fg(1)
      END
 
 
 
!       **********************************
 
      SUBROUTINE CIK01(Z,Cbi0,Cdi0,Cbi1,Cdi1,Cbk0,Cdk0,Cbk1,Cdk1)
!
!       ==========================================================
!       Purpose: Compute modified Bessel functions I0(z), I1(z),
!                K0(z), K1(z), and their derivatives for a
!                complex argument
!       Input :  z --- Complex argument
!       Output:  CBI0 --- I0(z)
!                CDI0 --- I0'(z)
!                CBI1 --- I1(z)
!                CDI1 --- I1'(z)
!                CBK0 --- K0(z)
!                CDK0 --- K0'(z)
!                CBK1 --- K1(z)
!                CDK1 --- K1'(z)
!       ==========================================================
!
      IMPLICIT NONE
!*--CIK0113467
      DOUBLE PRECISION a , a0 , a1 , b , pi , w0
      COMPLEX*16 ca , cb , Cbi0 , Cbi1 , Cbk0 , Cbk1 , Cdi0 , Cdi1 ,    &
               & Cdk0 , Cdk1 , ci , cr , cs , ct , cw , Z , z1 , z2 ,   &
               & zr , zr2
      INTEGER k , k0
      DIMENSION a(12) , b(12) , a1(10)
      pi = 3.141592653589793D0
      ci = (0.0D0,1.0D0)
      a0 = ABS(Z)
      z2 = Z*Z
      z1 = Z
      IF ( a0==0.0D0 ) THEN
         Cbi0 = (1.0D0,0.0D0)
         Cbi1 = (0.0D0,0.0D0)
         Cdi0 = (0.0D0,0.0D0)
         Cdi1 = (0.5D0,0.0D0)
         Cbk0 = (1.0D+300,0.0D0)
         Cbk1 = (1.0D+300,0.0D0)
         Cdk0 = -(1.0D+300,0.0D0)
         Cdk1 = -(1.0D+300,0.0D0)
         RETURN
      ENDIF
      IF ( DBLE(Z)<0.0 ) z1 = -Z
      IF ( a0<=18.0 ) THEN
         Cbi0 = (1.0D0,0.0D0)
         cr = (1.0D0,0.0D0)
         DO k = 1 , 50
            cr = 0.25D0*cr*z2/(k*k)
            Cbi0 = Cbi0 + cr
            IF ( ABS(cr/Cbi0)<1.0D-15 ) GOTO 50
         ENDDO
 50      Cbi1 = (1.0D0,0.0D0)
         cr = (1.0D0,0.0D0)
         DO k = 1 , 50
            cr = 0.25D0*cr*z2/(k*(k+1))
            Cbi1 = Cbi1 + cr
            IF ( ABS(cr/Cbi1)<1.0D-15 ) GOTO 100
         ENDDO
 100     Cbi1 = 0.5D0*z1*Cbi1
      ELSE
         DATA a/0.125D0 , 7.03125D-2 , 7.32421875D-2 ,                  &
            & 1.1215209960938D-1 , 2.2710800170898D-1 ,                 &
            & 5.7250142097473D-1 , 1.7277275025845D0 ,                  &
            & 6.0740420012735D0 , 2.4380529699556D01 ,                  &
            & 1.1001714026925D02 , 5.5133589612202D02 ,                 &
            & 3.0380905109224D03/
         DATA b/ - 0.375D0 , -1.171875D-1 , -1.025390625D-1 ,           &
            & -1.4419555664063D-1 , -2.7757644653320D-1 ,               &
            & -6.7659258842468D-1 , -1.9935317337513D0 ,                &
            & -6.8839142681099D0 , -2.7248827311269D01 ,                &
            & -1.2159789187654D02 , -6.0384407670507D02 ,               &
            & -3.3022722944809D03/
         k0 = 12
         IF ( a0>=35.0 ) k0 = 9
         IF ( a0>=50.0 ) k0 = 7
         ca = EXP(z1)/SQRT(2.0D0*pi*z1)
         Cbi0 = (1.0D0,0.0D0)
         zr = 1.0D0/z1
         DO k = 1 , k0
            Cbi0 = Cbi0 + a(k)*zr**k
         ENDDO
         Cbi0 = ca*Cbi0
         Cbi1 = (1.0D0,0.0D0)
         DO k = 1 , k0
            Cbi1 = Cbi1 + b(k)*zr**k
         ENDDO
         Cbi1 = ca*Cbi1
      ENDIF
      IF ( a0<=9.0 ) THEN
         cs = (0.0D0,0.0D0)
         ct = -LOG(0.5D0*z1) - 0.5772156649015329D0
         w0 = 0.0D0
         cr = (1.0D0,0.0D0)
         DO k = 1 , 50
            w0 = w0 + 1.0D0/k
            cr = 0.25D0*cr/(k*k)*z2
            cs = cs + cr*(w0+ct)
            IF ( ABS((cs-cw)/cs)<1.0D-15 ) GOTO 150
            cw = cs
         ENDDO
 150     Cbk0 = ct + cs
      ELSE
         DATA a1/0.125D0 , 0.2109375D0 , 1.0986328125D0 ,               &
            & 1.1775970458984D01 , 2.1461706161499D02 ,                 &
            & 5.9511522710323D03 , 2.3347645606175D05 ,                 &
            & 1.2312234987631D07 , 8.401390346421D08 ,                  &
            & 7.2031420482627D10/
         cb = 0.5D0/z1
         zr2 = 1.0D0/z2
         Cbk0 = (1.0D0,0.0D0)
         DO k = 1 , 10
            Cbk0 = Cbk0 + a1(k)*zr2**k
         ENDDO
         Cbk0 = cb*Cbk0/Cbi0
      ENDIF
      Cbk1 = (1.0D0/z1-Cbi1*Cbk0)/Cbi0
      IF ( DBLE(Z)<0.0 ) THEN
         IF ( DIMAG(Z)<0.0 ) Cbk0 = Cbk0 + ci*pi*Cbi0
         IF ( DIMAG(Z)>0.0 ) Cbk0 = Cbk0 - ci*pi*Cbi0
         IF ( DIMAG(Z)<0.0 ) Cbk1 = -Cbk1 + ci*pi*Cbi1
         IF ( DIMAG(Z)>0.0 ) Cbk1 = -Cbk1 - ci*pi*Cbi1
         Cbi1 = -Cbi1
      ENDIF
      Cdi0 = Cbi1
      Cdi1 = Cbi0 - 1.0D0/Z*Cbi1
      Cdk0 = -Cbk1
      Cdk1 = -Cbk0 - 1.0D0/Z*Cbk1
      END
 
!       **********************************
 
      SUBROUTINE CPSI(X,Y,Psr,Psi)
!
!       =============================================
!       Purpose: Compute the psi function for a
!                complex argument
!       Input :  x   --- Real part of z
!                y   --- Imaginary part of z
!       Output:  PSR --- Real part of psi(z)
!                PSI --- Imaginary part of psi(z)
!       =============================================
!
      IMPLICIT NONE
!*--CPSI13594
      DOUBLE PRECISION a , ct2 , pi , Psi , Psr , ri , rr , th , tm ,   &
                     & tn , X , x0 , x1 , Y , y1 , z0 , z2
      INTEGER k , n
      DIMENSION a(8)
      DATA a/ - .8333333333333D-01 , .83333333333333333D-02 ,           &
         & -.39682539682539683D-02 , .41666666666666667D-02 ,           &
         & -.75757575757575758D-02 , .21092796092796093D-01 ,           &
         & -.83333333333333333D-01 , .4432598039215686D0/
      pi = 3.141592653589793D0
      IF ( Y==0.0D0 .AND. X==INT(X) .AND. X<=0.0D0 ) THEN
         Psr = 1.0D+300
         Psi = 0.0D0
      ELSE
         x1 = X
         y1 = Y
         IF ( X<0.0D0 ) THEN
            X = -X
            Y = -Y
         ENDIF
         x0 = X
         n = 0
         IF ( X<8.0D0 ) THEN
            n = 8 - INT(X)
            x0 = X + n
         ENDIF
         th = 0.0D0
         IF ( x0==0.0D0 .AND. Y/=0.0D0 ) th = 0.5D0*pi
         IF ( x0/=0.0D0 ) th = DATAN(Y/x0)
         z2 = x0*x0 + Y*Y
         z0 = DSQRT(z2)
         Psr = DLOG(z0) - 0.5D0*x0/z2
         Psi = th + 0.5D0*Y/z2
         DO k = 1 , 8
            Psr = Psr + a(k)*z2**(-k)*DCOS(2.0D0*k*th)
            Psi = Psi - a(k)*z2**(-k)*DSIN(2.0D0*k*th)
         ENDDO
         IF ( X<8.0D0 ) THEN
            rr = 0.0D0
            ri = 0.0D0
            DO k = 1 , n
               rr = rr + (x0-k)/((x0-k)**2.0D0+Y*Y)
               ri = ri + Y/((x0-k)**2.0D0+Y*Y)
            ENDDO
            Psr = Psr - rr
            Psi = Psi + ri
         ENDIF
         IF ( x1<0.0D0 ) THEN
            tn = DTAN(pi*X)
            tm = DTANH(pi*Y)
            ct2 = tn*tn + tm*tm
            Psr = Psr + X/(X*X+Y*Y) + pi*(tn-tn*tm*tm)/ct2
            Psi = Psi - Y/(X*X+Y*Y) - pi*tm*(1.0D0+tn*tn)/ct2
            X = x1
            Y = y1
         ENDIF
      ENDIF
      END
 
!       **********************************
 
      SUBROUTINE SPHY(N,X,Nm,Sy,Dy)
!
!       ======================================================
!       Purpose: Compute spherical Bessel functions yn(x) and
!                their derivatives
!       Input :  x --- Argument of yn(x) ( x ≥ 0 )
!                n --- Order of yn(x) ( n = 0,1,… )
!       Output:  SY(n) --- yn(x)
!                DY(n) --- yn'(x)
!                NM --- Highest order computed
!       ======================================================
!
      IMPLICIT NONE
!*--SPHY13671
      DOUBLE PRECISION Dy , f , f0 , f1 , Sy , X
      INTEGER k , N , Nm
      DIMENSION Sy(0:N) , Dy(0:N)
      Nm = N
      IF ( X<1.0D-60 ) THEN
         DO k = 0 , N
            Sy(k) = -1.0D+300
            Dy(k) = 1.0D+300
         ENDDO
         RETURN
      ENDIF
      Sy(0) = -DCOS(X)/X
      f0 = Sy(0)
      Dy(0) = (DSIN(X)+DCOS(X)/X)/X
      IF ( N<1 ) RETURN
      Sy(1) = (Sy(0)-DSIN(X))/X
      f1 = Sy(1)
      DO k = 2 , N
         f = (2.0D0*k-1.0D0)*f1/X - f0
         Sy(k) = f
         IF ( DABS(f)>=1.0D+300 ) GOTO 100
         f0 = f1
         f1 = f
      ENDDO
 100  Nm = k - 1
      DO k = 1 , Nm
         Dy(k) = Sy(k-1) - (k+1.0D0)*Sy(k)/X
      ENDDO
      END
 
!       **********************************
 
      SUBROUTINE JELP(U,Hk,Esn,Ecn,Edn,Eph)
!
!       ========================================================
!       Purpose: Compute Jacobian elliptic functions sn u, cn u
!                and dn u
!       Input  : u   --- Argument of Jacobian elliptic functions
!                Hk  --- Modulus k ( 0 ≤ k ≤ 1 )
!       Output : ESN --- sn u
!                ECN --- cn u
!                EDN --- dn u
!                EPH --- phi ( in degrees )
!       ========================================================
!
      IMPLICIT NONE
!*--JELP13721
      DOUBLE PRECISION a , a0 , b , b0 , c , d , dn , Ecn , Edn , Eph , &
                     & Esn , Hk , pi , r , sa , t , U
      INTEGER j , n
      DIMENSION r(40)
      pi = 3.14159265358979D0
      a0 = 1.0D0
      b0 = DSQRT(1.0D0-Hk*Hk)
      DO n = 1 , 40
         a = (a0+b0)/2.0D0
         b = DSQRT(a0*b0)
         c = (a0-b0)/2.0D0
         r(n) = c/a
         IF ( c<1.0D-7 ) GOTO 100
         a0 = a
         b0 = b
      ENDDO
 100  dn = 2.0D0**n*a*U
      d = 0.0D0
      DO j = n , 1 , -1
         t = r(j)*DSIN(dn)
         sa = DATAN(t/DSQRT(DABS(1.0D0-t*t)))
         d = .5D0*(dn+sa)
         dn = d
      ENDDO
      Eph = d*180.0D0/pi
      Esn = DSIN(d)
      Ecn = DCOS(d)
      Edn = DSQRT(1.0D0-Hk*Hk*Esn*Esn)
      END
