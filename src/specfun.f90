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

    module specfun_module

      use iso_fortran_env, only: wp => real64

      implicit none

      contains

      function dnan()
      implicit none
      real(wp) dnan
      dnan = 0.0d0
      dnan = 0.0d0/dnan
      end

      function dinf()
      implicit none
      real(wp) dinf
      dinf = 1.0d300
      dinf = dinf*dinf
      end

      subroutine cpdsa(n,z,Cdn)
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
      implicit none
      complex(wp) ca0 , cb0 , Cdn , cdw , cr , z
      real(wp) eps , g0 , g1 , ga0 , gm , pd , pi , sq2 , va0 , &
                     & vm , vt , xn
      integer m , n
      eps = 1.0d-15
      pi = 3.141592653589793d0
      sq2 = sqrt(2.0d0)
      ca0 = exp(-.25d0*z*z)
      va0 = 0.5d0*(1.0d0-n)
      if ( n==0.0 ) then
         Cdn = ca0
      elseif ( abs(z)==0.0 ) then
         if ( va0<=0.0 .and. va0==int(va0) ) then
            Cdn = 0.0d0
         else
            call gaih(va0,ga0)
            pd = sqrt(pi)/(2.0d0**(-.5d0*n)*ga0)
            Cdn = dcmplx(pd,0.0d0)
         endif
      else
         xn = -n
         call gaih(xn,g1)
         cb0 = 2.0d0**(-0.5d0*n-1.0d0)*ca0/g1
         vt = -.5d0*n
         call gaih(vt,g0)
         Cdn = dcmplx(g0,0.0d0)
         cr = (1.0d0,0.0d0)
         do m = 1 , 250
            vm = .5d0*(m-n)
            call gaih(vm,gm)
            cr = -cr*sq2*z/m
            cdw = gm*cr
            Cdn = Cdn + cdw
            if ( abs(cdw)<abs(Cdn)*eps ) exit
         enddo
         Cdn = cb0*Cdn
      endif
      end



!       **********************************

      subroutine cfs(z,Zf,Zd)
!
!       =========================================================
!       Purpose: Compute complex Fresnel Integral S(z) and S'(z)
!       Input :  z  --- Argument of S(z)
!       Output:  ZF --- S(z)
!                ZD --- S'(z)
!       =========================================================
!
      implicit none
      complex(wp) cf , cf0 , cf1 , cg , cr , d , s , z , z0 , Zd , Zf ,  &
               & zp , zp2
      real(wp) eps , pi , w0 , wb , wb0
      integer k , m
      eps = 1.0d-14
      pi = 3.141592653589793d0
      w0 = abs(z)
      zp = 0.5d0*pi*z*z
      zp2 = zp*zp
      z0 = (0.0d0,0.0d0)
      if ( z==z0 ) then
         s = z0
      elseif ( w0<=2.5 ) then
         s = z*zp/3.0d0
         cr = s
         wb0 = 0.0d0
         do k = 1 , 80
            cr = -.5d0*cr*(4.0d0*k-1.0d0)/k/(2.0d0*k+1.0d0)             &
               & /(4.0d0*k+3.0d0)*zp2
            s = s + cr
            wb = abs(s)
            if ( abs(wb-wb0)<eps .and. k>10 ) exit
            wb0 = wb
         enddo
      elseif ( w0>2.5 .and. w0<4.5 ) then
         m = 85
         s = z0
         cf1 = z0
         cf0 = (1.0d-100,0.0d0)
         do k = m , 0 , -1
            cf = (2.0d0*k+3.0d0)*cf0/zp - cf1
            if ( k/=int(k/2)*2 ) s = s + cf
            cf1 = cf0
            cf0 = cf
         enddo
         s = 2.0d0/(pi*z)*sin(zp)/cf*s
      else
!          Auxiliary functions f(z) and g(z) can be computed using an
!          asymptotic expansion in the right quadrant |arg(z)| <= pi/4, not pi/2
!          as sometimes suggested. Use the symmetry S(z) = -iS(-iz).
!          Interestingly, most of the expansion code is the same across
!          the quadrants. (The forth power in Z is the equalizer here.)
!          Only one constant has to be adapted.
         if ( dimag(z)>-dble(z) .and. dimag(z)<=dble(z) ) then
!            right quadrant
            d = dcmplx(.5d0,0.0d0)
         elseif ( dimag(z)>dble(z) .and. dimag(z)>=-dble(z) ) then
!            upper quadrant
            d = dcmplx(0.0d0,-.5d0)
         elseif ( dimag(z)<-dble(z) .and. dimag(z)>=dble(z) ) then
!            left quadrant
            d = dcmplx(-.5d0,0.0d0)
         else
!            lower quadrant
            d = dcmplx(0.0d0,.5d0)
         endif
         cr = (1.0d0,0.0d0)
         cf = (1.0d0,0.0d0)
         do k = 1 , 20
            cr = -.25d0*cr*(4.0d0*k-1.0d0)*(4.0d0*k-3.0d0)/zp2
            cf = cf + cr
         enddo
         cr = (1.0d0,0.0d0)
         cg = (1.0d0,0.0d0)
         do k = 1 , 12
            cr = -.25d0*cr*(4.0d0*k+1.0d0)*(4.0d0*k-1.0d0)/zp2
            cg = cg + cr
         enddo
         cg = cg/(pi*z*z)
         s = d - (cf*cos(zp)+cg*sin(zp))/(pi*z)
      endif
      Zf = s
      Zd = sin(0.5*pi*z*z)
      end

!       **********************************

      subroutine lqmn(Mm,m,n,x,Qm,Qd)
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
      implicit none
      integer i , j , k , km , ls , m , Mm , n
      real(wp) q0 , q1 , q10 , Qd , qf , qf0 , qf1 , qf2 , Qm , &
                     & x , xq , xs
      dimension Qm(0:Mm,0:n) , Qd(0:Mm,0:n)
      if ( abs(x)==1.0d0 ) then
         do i = 0 , m
            do j = 0 , n
               Qm(i,j) = 1.0d+300
               Qd(i,j) = 1.0d+300
            enddo
         enddo
         return
      endif
      ls = 1
      if ( abs(x)>1.0d0 ) ls = -1
      xs = ls*(1.0d0-x*x)
      xq = sqrt(xs)
      q0 = 0.5d0*log(abs((x+1.0d0)/(x-1.0d0)))
      if ( abs(x)<1.0001d0 ) then
         Qm(0,0) = q0
         Qm(0,1) = x*q0 - 1.0d0
         Qm(1,0) = -1.0d0/xq
         Qm(1,1) = -ls*xq*(q0+x/(1.0d0-x*x))
         do i = 0 , 1
            do j = 2 , n
               Qm(i,j) = ((2.0d0*j-1.0d0)*x*Qm(i,j-1)-(j+i-1.0d0)       &
                       & *Qm(i,j-2))/(j-i)
            enddo
         enddo
         do j = 0 , n
            do i = 2 , m
               Qm(i,j) = -2.0d0*(i-1.0d0)*x/xq*Qm(i-1,j)                &
                       & - ls*(j+i-1.0d0)*(j-i+2.0d0)*Qm(i-2,j)
            enddo
         enddo
      else
         if ( abs(x)>1.1d0 ) then
            km = 40 + m + n
         else
            km = (40+m+n)*int(-1.0-1.8*log(x-1.0))
         endif
         qf2 = 0.0d0
         qf1 = 1.0d0
         qf0 = 0.0d0
         do k = km , 0 , -1
            qf0 = ((2*k+3.0d0)*x*qf1-(k+2.0d0)*qf2)/(k+1.0d0)
            if ( k<=n ) Qm(0,k) = qf0
            qf2 = qf1
            qf1 = qf0
         enddo
         do k = 0 , n
            Qm(0,k) = q0*Qm(0,k)/qf0
         enddo
         qf2 = 0.0d0
         qf1 = 1.0d0
         do k = km , 0 , -1
            qf0 = ((2*k+3.0d0)*x*qf1-(k+1.0d0)*qf2)/(k+2.0d0)
            if ( k<=n ) Qm(1,k) = qf0
            qf2 = qf1
            qf1 = qf0
         enddo
         q10 = -1.0d0/xq
         do k = 0 , n
            Qm(1,k) = q10*Qm(1,k)/qf0
         enddo
         do j = 0 , n
            q0 = Qm(0,j)
            q1 = Qm(1,j)
            do i = 0 , m - 2
               qf = -2.0d0*(i+1)*x/xq*q1 + (j-i)*(j+i+1.0d0)*q0
               Qm(i+2,j) = qf
               q0 = q1
               q1 = qf
            enddo
         enddo
      endif
      Qd(0,0) = ls/xs
      do j = 1 , n
         Qd(0,j) = ls*j*(Qm(0,j-1)-x*Qm(0,j))/xs
      enddo
      do j = 0 , n
         do i = 1 , m
            Qd(i,j) = ls*i*x/xs*Qm(i,j) + (i+j)*(j-i+1.0d0)/xq*Qm(i-1,j)
         enddo
      enddo
      end

!       **********************************

      subroutine clpmn(Mm,m,n,x,y,Ntype,Cpm,Cpd)
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
      implicit none
      complex(wp) Cpd , Cpm , z , zq , zs
      real(wp)  x , y
      integer i , j , ls , m , Mm , n , Ntype
      dimension Cpm(0:Mm,0:n) , Cpd(0:Mm,0:n)
      z = dcmplx(x,y)
      do i = 0 , n
         do j = 0 , m
            Cpm(j,i) = (0.0d0,0.0d0)
            Cpd(j,i) = (0.0d0,0.0d0)
         enddo
      enddo
      Cpm(0,0) = (1.0d0,0.0d0)
      if ( n==0 ) return
      if ( abs(x)==1.0d0 .and. y==0.0d0 ) then
         do i = 1 , n
            Cpm(0,i) = x**i
            Cpd(0,i) = 0.5d0*i*(i+1)*x**(i+1)
         enddo
         do j = 1 , n
            do i = 1 , m
               if ( i==1 ) then
                  Cpd(i,j) = dinf()
               elseif ( i==2 ) then
                  Cpd(i,j) = -0.25d0*(j+2)*(j+1)*j*(j-1)*x**(j+1)
               endif
            enddo
         enddo
         return
      endif
      if ( Ntype==2 ) then
!       sqrt(1 - z^2) with branch cut on |x|>1
         zs = (1.0d0-z*z)
         zq = -sqrt(zs)
         ls = -1
      else
!       sqrt(z^2 - 1) with branch cut between [-1, 1]
         zs = (z*z-1.0d0)
         zq = sqrt(zs)
         if ( x<0d0 ) zq = -zq
         ls = 1
      endif
      do i = 1 , m
!       DLMF 14.7.15
         Cpm(i,i) = (2.0d0*i-1.0d0)*zq*Cpm(i-1,i-1)
      enddo
      do i = 0 , min(m,n-1)
!       DLMF 14.10.7
         Cpm(i,i+1) = (2.0d0*i+1.0d0)*z*Cpm(i,i)
      enddo
      do i = 0 , m
         do j = i + 2 , n
!       DLMF 14.10.3
            Cpm(i,j) = ((2.0d0*j-1.0d0)*z*Cpm(i,j-1)-(i+j-1.0d0)        &
                     & *Cpm(i,j-2))/(j-i)
         enddo
      enddo
      Cpd(0,0) = (0.0d0,0.0d0)
      do j = 1 , n
!       DLMF 14.10.5
         Cpd(0,j) = ls*j*(z*Cpm(0,j)-Cpm(0,j-1))/zs
      enddo
      do i = 1 , m
         do j = i , n
!       derivative of DLMF 14.7.11 & DLMF 14.10.6 for type 3
!       derivative of DLMF 14.7.8 & DLMF 14.10.1 for type 2
            Cpd(i,j) = ls*(-i*z*Cpm(i,j)/zs+(j+i)*(j-i+1.0d0)           &
                     & /zq*Cpm(i-1,j))
         enddo
      enddo
      end

!       **********************************

      subroutine vvsa(Va,x,Pv)
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
      implicit none
      real(wp) a0 , ep , eps , fac , g1 , ga0 , gm , gw , pi ,  &
                     & Pv , r , r1 , sq2 , sv , sv0 , v1 , Va , va0 ,   &
                     & vb0 , vm
      real(wp) x
      integer m
      eps = 1.0d-15
      pi = 3.141592653589793d0
      ep = exp(-.25d0*x*x)
      va0 = 1.0d0 + 0.5d0*Va
      if ( x/=0.0 ) then
         sq2 = sqrt(2.0d0)
         a0 = 2.0d0**(-.5d0*Va)*ep/(2.0d0*pi)
         sv = sin(-(Va+.5d0)*pi)
         v1 = -.5d0*Va
         call gamma2(v1,g1)
         Pv = (sv+1.0d0)*g1
         r = 1.0d0
         fac = 1.0d0
         do m = 1 , 250
            vm = .5d0*(m-Va)
            call gamma2(vm,gm)
            r = r*sq2*x/m
            fac = -fac
            gw = fac*sv + 1.0d0
            r1 = gw*r*gm
            Pv = Pv + r1
            if ( abs(r1/Pv)<eps .and. gw/=0.0 ) exit
         enddo
         Pv = a0*Pv
      elseif ( va0<=0.0 .and. va0==int(va0) .or. Va==0.0 ) then
         Pv = 0.0d0
      else
         vb0 = -0.5d0*Va
         sv0 = sin(va0*pi)
         call gamma2(va0,ga0)
         Pv = 2.0d0**vb0*sv0/ga0
      endif
      end



!       **********************************
!       SciPy: Changed P from a character array to an integer array.
      subroutine jdzo(Nt,n,m,p,Zo)
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
      implicit none
      real(wp) bj , dj , fj , x , x0 , x1 , x2 , xm , Zo , zoc
      integer i , j , k , l , l0 , l1 , l2 , m , m1 , mm , n , n1 , nm ,&
            & Nt
      integer p(1400) , p1(70)
      dimension n(1400) , m(1400) , Zo(0:1400) , n1(70) , m1(70) ,      &
              & zoc(0:70) , bj(101) , dj(101) , fj(101)
      x = 0
      zoc(0) = 0
      if ( Nt<600 ) then
         xm = -1.0 + 2.248485*Nt**0.5 - .0159382*Nt +                   &
            & 3.208775e-4*Nt**1.5
         nm = int(14.5+.05875*Nt)
         mm = int(.02*Nt) + 6
      else
         xm = 5.0 + 1.445389*Nt**.5 + .01889876*Nt - 2.147763e-4*Nt**1.5
         nm = int(27.8+.0327*Nt)
         mm = int(.01088*Nt) + 10
      endif
      l0 = 0
      do i = 1 , nm
         x1 = .407658 + .4795504*(i-1)**.5 + .983618*(i-1)
         x2 = 1.99535 + .8333883*(i-1)**.5 + .984584*(i-1)
         l1 = 0
         do j = 1 , mm
            if ( i/=1 .or. j/=1 ) then
               x = x1
 10            call bjndd(i,x,bj,dj,fj)
               x0 = x
               x = x - dj(i)/fj(i)
               if ( x1>xm ) goto 20
               if ( abs(x-x0)>1.0d-10 ) goto 10
            endif
            l1 = l1 + 1
            n1(l1) = i - 1
            m1(l1) = j
            if ( i==1 ) m1(l1) = j - 1
            p1(l1) = 1
            zoc(l1) = x
            if ( i<=15 ) then
               x1 = x + 3.057 + .0122*(i-1) + (1.555+.41575*(i-1))/(j+1)&
                  & **2
            else
               x1 = x + 2.918 + .01924*(i-1) + (6.26+.13205*(i-1))/(j+1)&
                  & **2
            endif
 20         x = x2
 40         call bjndd(i,x,bj,dj,fj)
            x0 = x
            x = x - bj(i)/dj(i)
            if ( x<=xm ) then
               if ( abs(x-x0)>1.0d-10 ) goto 40
               l1 = l1 + 1
               n1(l1) = i - 1
               m1(l1) = j
               p1(l1) = 0
               zoc(l1) = x
               if ( i<=15 ) then
                  x2 = x + 3.11 + .0138*(i-1) + (.04832+.2804*(i-1))    &
                     & /(j+1)**2
               else
                  x2 = x + 3.001 + .0105*(i-1) + (11.52+.48525*(i-1))   &
                     & /(j+3)**2
               endif
            endif
         enddo
         l = l0 + l1
         l2 = l
 50      if ( l0==0 ) then
            do k = 1 , l
               Zo(k) = zoc(k)
               n(k) = n1(k)
               m(k) = m1(k)
               p(k) = p1(k)
            enddo
            l1 = 0
         elseif ( l0/=0 ) then
            if ( Zo(l0)>=zoc(l1) ) then
               Zo(l0+l1) = Zo(l0)
               n(l0+l1) = n(l0)
               m(l0+l1) = m(l0)
               p(l0+l1) = p(l0)
               l0 = l0 - 1
            else
               Zo(l0+l1) = zoc(l1)
               n(l0+l1) = n1(l1)
               m(l0+l1) = m1(l1)
               p(l0+l1) = p1(l1)
               l1 = l1 - 1
            endif
         endif
         if ( l1/=0 ) goto 50
         l0 = l2
      enddo
      end



!       **********************************

      subroutine cbk(m,n,c,Cv,Qt,Ck,Bk)
!
!       =====================================================
!       Purpose: Compute coefficient Bk's for oblate radial
!                functions with a small argument
!       =====================================================
!
      implicit none
      real(wp) Bk , c , Ck , Cv , eps , Qt , r1 , s1 , sw , t , &
                     & u , v , w
      integer i , i1 , ip , j , k , m , n , n2 , nm
      dimension Bk(200) , Ck(200) , u(200) , v(200) , w(200)
      eps = 1.0d-14
      ip = 1
      if ( n-m==2*int((n-m)/2) ) ip = 0
      nm = 25 + int(0.5*(n-m)+c)
      u(1) = 0.0d0
      n2 = nm - 2
      do j = 2 , n2
         u(j) = c*c
      enddo
      do j = 1 , n2
         v(j) = (2.0*j-1.0-ip)*(2.0*(j-m)-ip) + m*(m-1.0) - Cv
      enddo
      do j = 1 , nm - 1
         w(j) = (2.0*j-ip)*(2.0*j+1.0-ip)
      enddo
      if ( ip==0 ) then
         sw = 0.0d0
         do k = 0 , n2 - 1
            s1 = 0.0d0
            i1 = k - m + 1
            do i = i1 , nm
               if ( i>=0 ) then
                  r1 = 1.0d0
                  do j = 1 , k
                     r1 = r1*(i+m-j)/j
                  enddo
                  s1 = s1 + Ck(i+1)*(2.0*i+m)*r1
                  if ( abs(s1-sw)<abs(s1)*eps ) exit
                  sw = s1
               endif
            enddo
            Bk(k+1) = Qt*s1
         enddo
      elseif ( ip==1 ) then
         sw = 0.0d0
         do k = 0 , n2 - 1
            s1 = 0.0d0
            i1 = k - m + 1
            do i = i1 , nm
               if ( i>=0 ) then
                  r1 = 1.0d0
                  do j = 1 , k
                     r1 = r1*(i+m-j)/j
                  enddo
                  if ( i>0 ) s1 = s1 + Ck(i)*(2.0*i+m-1)*r1
                  s1 = s1 - Ck(i+1)*(2.0*i+m)*r1
                  if ( abs(s1-sw)<abs(s1)*eps ) exit
                  sw = s1
               endif
            enddo
            Bk(k+1) = Qt*s1
         enddo
      endif
      w(1) = w(1)/v(1)
      Bk(1) = Bk(1)/v(1)
      do k = 2 , n2
         t = v(k) - w(k-1)*u(k)
         w(k) = w(k)/t
         Bk(k) = (Bk(k)-Bk(k-1)*u(k))/t
      enddo
      do k = n2 - 1 , 1 , -1
         Bk(k) = Bk(k) - w(k)*Bk(k+1)
      enddo
      end



!       **********************************

      subroutine rmn2sp(m,n,c,x,Cv,Df,Kd,R2f,R2d)
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
      implicit none
      real(wp) c , ck1 , ck2 , Cv , Df , dn , eps , ga , gb ,   &
                     & gc , pd , pm , qd , qm , r1 , r2 , R2d , R2f ,   &
                     & r3 , r4
      real(wp) sd , sd0 , sd1 , sd2 , sdm , sf , spd1 , spd2 ,  &
                     & spl , su0 , su1 , su2 , sum , sw , x
      integer ip , j , j1 , j2 , k , Kd , ki , l1 , m , n , nm , nm1 ,  &
            & nm2 , nm3
      dimension pm(0:251) , pd(0:251) , qm(0:251) , qd(0:251) , dn(200) &
              & , Df(200)
      if ( abs(Df(1))<1.0d-280 ) then
         R2f = 1.0d+300
         R2d = 1.0d+300
         return
      endif
      eps = 1.0d-14
      ip = 1
      nm1 = int((n-m)/2)
      if ( n-m==2*nm1 ) ip = 0
      nm = 25 + nm1 + int(c)
      nm2 = 2*nm + m
      call kmn(m,n,c,Cv,Kd,Df,dn,ck1,ck2)
      call lpmns(m,nm2,x,pm,pd)
      call lqmns(m,nm2,x,qm,qd)
      su0 = 0.0d0
      sw = 0.0d0
      do k = 1 , nm
         j = 2*k - 2 + m + ip
         su0 = su0 + Df(k)*qm(j)
         if ( k>nm1 .and. abs(su0-sw)<abs(su0)*eps ) exit
         sw = su0
      enddo
      sd0 = 0.0d0
      do k = 1 , nm
         j = 2*k - 2 + m + ip
         sd0 = sd0 + Df(k)*qd(j)
         if ( k>nm1 .and. abs(sd0-sw)<abs(sd0)*eps ) exit
         sw = sd0
      enddo
      su1 = 0.0d0
      sd1 = 0.0d0
      do k = 1 , m
         j = m - 2*k + ip
         if ( j<0 ) j = -j - 1
         su1 = su1 + dn(k)*qm(j)
         sd1 = sd1 + dn(k)*qd(j)
      enddo
      ga = ((x-1.0d0)/(x+1.0d0))**(0.5d0*m)
      do k = 1 , m
         j = m - 2*k + ip
         if ( j<0 ) then
            if ( j<0 ) j = -j - 1
            r1 = 1.0d0
            do j1 = 1 , j
               r1 = (m+j1)*r1
            enddo
            r2 = 1.0d0
            do j2 = 1 , m - j - 2
               r2 = j2*r2
            enddo
            r3 = 1.0d0
            sf = 1.0d0
            do l1 = 1 , j
               r3 = 0.5d0*r3*(-j+l1-1.0)*(j+l1)/((m+l1)*l1)*(1.0-x)
               sf = sf + r3
            enddo
            if ( m-j>=2 ) gb = (m-j-1.0d0)*r2
            if ( m-j<=1 ) gb = 1.0d0
            spl = r1*ga*gb*sf
            su1 = su1 + (-1)**(j+m)*dn(k)*spl
            spd1 = m/(x*x-1.0d0)*spl
            gc = 0.5d0*j*(j+1.0)/(m+1.0)
            sd = 1.0d0
            r4 = 1.0d0
            do l1 = 1 , j - 1
               r4 = 0.5d0*r4*(-j+l1)*(j+l1+1.0)/((m+l1+1.0)*l1)*(1.0-x)
               sd = sd + r4
            enddo
            spd2 = r1*ga*gb*gc*sd
            sd1 = sd1 + (-1)**(j+m)*dn(k)*(spd1+spd2)
         endif
      enddo
      su2 = 0.0d0
      ki = (2*m+1+ip)/2
      nm3 = nm + ki
      do k = ki , nm3
         j = 2*k - 1 - m - ip
         su2 = su2 + dn(k)*pm(j)
         if ( j>m .and. abs(su2-sw)<abs(su2)*eps ) exit
         sw = su2
      enddo
      sd2 = 0.0d0
      do k = ki , nm3
         j = 2*k - 1 - m - ip
         sd2 = sd2 + dn(k)*pd(j)
         if ( j>m .and. abs(sd2-sw)<abs(sd2)*eps ) exit
         sw = sd2
      enddo
      sum = su0 + su1 + su2
      sdm = sd0 + sd1 + sd2
      R2f = sum/ck2
      R2d = sdm/ck2
      end



!       **********************************

      subroutine bernob(n,Bn)
!
!       ======================================
!       Purpose: Compute Bernoulli number Bn
!       Input :  n --- Serial number
!       Output:  BN(n) --- Bn
!       ======================================
!
      implicit none
      real(wp) Bn , r1 , r2 , s , tpi
      integer k , m , n
      dimension Bn(0:n)
      tpi = 6.283185307179586d0
      Bn(0) = 1.0d0
      Bn(1) = -0.5d0
      Bn(2) = 1.0d0/6.0d0
      r1 = (2.0d0/tpi)**2
      do m = 4 , n , 2
         r1 = -r1*(m-1)*m/(tpi*tpi)
         r2 = 1.0d0
         do k = 2 , 10000
            s = (1.0d0/k)**m
            r2 = r2 + s
            if ( s<1.0d-15 ) exit
         enddo
         Bn(m) = r1*r2
      enddo
      end

!       **********************************

      subroutine bernoa(n,Bn)
!
!       ======================================
!       Purpose: Compute Bernoulli number Bn
!       Input :  n --- Serial number
!       Output:  BN(n) --- Bn
!       ======================================
!
      implicit none
      real(wp) Bn , r , s
      integer j , k , m , n
      dimension Bn(0:n)
      Bn(0) = 1.0d0
      Bn(1) = -0.5d0
      do m = 2 , n
         s = -(1.0d0/(m+1.0d0)-0.5d0)
         do k = 2 , m - 1
            r = 1.0d0
            do j = 2 , k
               r = r*(j+m-k)/j
            enddo
            s = s - r*Bn(k)
         enddo
         Bn(m) = s
      enddo
      do m = 3 , n , 2
         Bn(m) = 0.0d0
      enddo
      end

!       **********************************

      subroutine qstar(m,n,c,Ck,Ck1,Qs,Qt)
!
!       =========================================================
!       Purpose: Compute Q*mn(-ic) for oblate radial functions
!                with a small argument
!       =========================================================
!
      implicit none
      real(wp) ap , c , Ck , Ck1 , Qs , qs0 , Qt , r , s , sk
      integer i , ip , k , l , m , n
      dimension ap(200) , Ck(200)
      ip = 1
      if ( n-m==2*int((n-m)/2) ) ip = 0
      r = 1.0d0/Ck(1)**2
      ap(1) = r
      do i = 1 , m
         s = 0.0d0
         do l = 1 , i
            sk = 0.0d0
            do k = 0 , l
               sk = sk + Ck(k+1)*Ck(l-k+1)
            enddo
            s = s + sk*ap(i-l+1)
         enddo
         ap(i+1) = -r*s
      enddo
      qs0 = ap(m+1)
      do l = 1 , m
         r = 1.0d0
         do k = 1 , l
            r = r*(2.0d0*k+ip)*(2.0d0*k-1.0d0+ip)/(2.0d0*k)**2
         enddo
         qs0 = qs0 + ap(m-l+1)*r
      enddo
      Qs = (-1)**ip*Ck1*(Ck1*qs0)/c
      Qt = -2.0d0/Ck1*Qs
      end



!       **********************************

      subroutine cv0(Kd,m,q,a0)
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
      implicit none
      real(wp) a0 , q , q2
      integer Kd , m
      q2 = q*q
      if ( m==0 ) then
         if ( q<=1.0 ) then
            a0 = (((.0036392*q2-.0125868)*q2+.0546875)*q2-.5)*q2
         elseif ( q<=10.0 ) then
            a0 = ((3.999267d-3*q-9.638957d-2)*q-.88297)*q + .5542818
         else
            call cvql(Kd,m,q,a0)
         endif
      elseif ( m==1 ) then
         if ( q<=1.0 .and. Kd==2 ) then
            a0 = (((-6.51e-4*q-.015625)*q-.125)*q+1.0)*q + 1.0
         elseif ( q<=1.0 .and. Kd==3 ) then
            a0 = (((-6.51e-4*q+.015625)*q-.125)*q-1.0)*q + 1.0
         elseif ( q<=10.0 .and. Kd==2 ) then
            a0 = (((-4.94603d-4*q+1.92917d-2)*q-.3089229)*q+1.33372)    &
               & *q + .811752
         elseif ( q<=10.0 .and. Kd==3 ) then
            a0 = ((1.971096d-3*q-5.482465d-2)*q-1.152218)*q + 1.10427
         else
            call cvql(Kd,m,q,a0)
         endif
      elseif ( m==2 ) then
         if ( q<=1.0 .and. Kd==1 ) then
            a0 = (((-.0036391*q2+.0125888)*q2-.0551939)*q2+.416667)     &
               & *q2 + 4.0
         elseif ( q<=1.0 .and. Kd==4 ) then
            a0 = (.0003617*q2-.0833333)*q2 + 4.0
         elseif ( q<=15 .and. Kd==1 ) then
            a0 = (((3.200972d-4*q-8.667445d-3)*q-1.829032d-4)           &
               & *q+.9919999)*q + 3.3290504
         elseif ( q<=10.0 .and. Kd==4 ) then
            a0 = ((2.38446d-3*q-.08725329)*q-4.732542d-3)*q + 4.00909
         else
            call cvql(Kd,m,q,a0)
         endif
      elseif ( m==3 ) then
         if ( q<=1.0 .and. Kd==2 ) then
            a0 = ((6.348e-4*q+.015625)*q+.0625)*q2 + 9.0
         elseif ( q<=1.0 .and. Kd==3 ) then
            a0 = ((6.348e-4*q-.015625)*q+.0625)*q2 + 9.0
         elseif ( q<=20.0 .and. Kd==2 ) then
            a0 = (((3.035731d-4*q-1.453021d-2)*q+.19069602)*q-.1039356) &
               & *q + 8.9449274
         elseif ( q<=15.0 .and. Kd==3 ) then
            a0 = ((9.369364d-5*q-.03569325)*q+.2689874)*q + 8.771735
         else
            call cvql(Kd,m,q,a0)
         endif
      elseif ( m==4 ) then
         if ( q<=1.0 .and. Kd==1 ) then
            a0 = ((-2.1e-6*q2+5.012e-4)*q2+.0333333)*q2 + 16.0
         elseif ( q<=1.0 .and. Kd==4 ) then
            a0 = ((3.7e-6*q2-3.669e-4)*q2+.0333333)*q2 + 16.0
         elseif ( q<=25.0 .and. Kd==1 ) then
            a0 = (((1.076676d-4*q-7.9684875d-3)*q+.17344854)*q-.5924058)&
               & *q + 16.620847
         elseif ( q<=20.0 .and. Kd==4 ) then
            a0 = ((-7.08719d-4*q+3.8216144d-3)*q+.1907493)*q + 15.744
         else
            call cvql(Kd,m,q,a0)
         endif
      elseif ( m==5 ) then
         if ( q<=1.0 .and. Kd==2 ) then
            a0 = ((6.8e-6*q+1.42e-5)*q2+.0208333)*q2 + 25.0
         elseif ( q<=1.0 .and. Kd==3 ) then
            a0 = ((-6.8e-6*q+1.42e-5)*q2+.0208333)*q2 + 25.0
         elseif ( q<=35.0 .and. Kd==2 ) then
            a0 = (((2.238231d-5*q-2.983416d-3)*q+.10706975)*q-.600205)  &
               & *q + 25.93515
         elseif ( q<=25.0 .and. Kd==3 ) then
            a0 = ((-7.425364d-4*q+2.18225d-2)*q+4.16399d-2)*q + 24.897
         else
            call cvql(Kd,m,q,a0)
         endif
      elseif ( m==6 ) then
         if ( q<=1.0 ) then
            a0 = (.4d-6*q2+.0142857)*q2 + 36.0
         elseif ( q<=40.0 .and. Kd==1 ) then
            a0 = (((-1.66846d-5*q+4.80263d-4)*q+2.53998d-2)*q-.181233)  &
               & *q + 36.423
         elseif ( q<=35.0 .and. Kd==4 ) then
            a0 = ((-4.57146d-4*q+2.16609d-2)*q-2.349616d-2)*q + 35.99251
         else
            call cvql(Kd,m,q,a0)
         endif
      elseif ( m==7 ) then
         if ( q<=10.0 ) then
            call cvqm(m,q,a0)
         elseif ( q<=50.0 .and. Kd==2 ) then
            a0 = (((-1.411114d-5*q+9.730514d-4)*q-3.097887d-3)          &
               & *q+3.533597d-2)*q + 49.0547
         elseif ( q<=40.0 .and. Kd==3 ) then
            a0 = ((-3.043872d-4*q+2.05511d-2)*q-9.16292d-2)*q + 49.19035
         else
            call cvql(Kd,m,q,a0)
         endif
      elseif ( m>=8 ) then
         if ( q<=3.*m ) then
            call cvqm(m,q,a0)
         elseif ( q>m*m ) then
            call cvql(Kd,m,q,a0)
         elseif ( m==8 .and. Kd==1 ) then
            a0 = (((8.634308d-6*q-2.100289d-3)*q+.169072)*q-4.64336)    &
               & *q + 109.4211
         elseif ( m==8 .and. Kd==4 ) then
            a0 = ((-6.7842d-5*q+2.2057d-3)*q+.48296)*q + 56.59
         elseif ( m==9 .and. Kd==2 ) then
            a0 = (((2.906435d-6*q-1.019893d-3)*q+.1101965)*q-3.821851)  &
               & *q + 127.6098
         elseif ( m==9 .and. Kd==3 ) then
            a0 = ((-9.577289d-5*q+.01043839)*q+.06588934)*q + 78.0198
         elseif ( m==10 .and. Kd==1 ) then
            a0 = (((5.44927d-7*q-3.926119d-4)*q+.0612099)*q-2.600805)   &
               & *q + 138.1923
         elseif ( m==10 .and. Kd==4 ) then
            a0 = ((-7.660143d-5*q+.01132506)*q-.09746023)*q + 99.29494
         elseif ( m==11 .and. Kd==2 ) then
            a0 = (((-5.67615d-7*q+7.152722d-6)*q+.01920291)*q-1.081583) &
               & *q + 140.88
         elseif ( m==11 .and. Kd==3 ) then
            a0 = ((-6.310551d-5*q+.0119247)*q-.2681195)*q + 123.667
         elseif ( m==12 .and. Kd==1 ) then
            a0 = (((-2.38351d-7*q-2.90139d-5)*q+.02023088)*q-1.289)     &
               & *q + 171.2723
         elseif ( m==12 .and. Kd==4 ) then
            a0 = (((3.08902d-7*q-1.577869d-4)*q+.0247911)*q-1.05454)    &
               & *q + 161.471
         endif
      endif
      end



!       **********************************

      subroutine cvqm(m,q,a0)
!
!       =====================================================
!       Purpose: Compute the characteristic value of Mathieu
!                functions for q ≤ m*m
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!       Output:  A0 --- Initial characteristic value
!       =====================================================
!
      implicit none
      real(wp) a0 , hm1 , hm3 , hm5 , q
      integer m
      hm1 = .5*q/(m*m-1.0)
      hm3 = .25*hm1**3/(m*m-4.0)
      hm5 = hm1*hm3*q/((m*m-1.0)*(m*m-9.0))
      a0 = m*m + q*(hm1+(5.0*m*m+7.0)*hm3+(9.0*m**4+58.0*m*m+29.0)*hm5)
      end

!       **********************************

      subroutine cvql(Kd,m,q,a0)
!
!       ========================================================
!       Purpose: Compute the characteristic value of Mathieu
!                functions  for q ≥ 3m
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!       Output:  A0 --- Initial characteristic value
!       ========================================================
!
      implicit none
      real(wp) a0 , c1 , cv1 , cv2 , d1 , d2 , d3 , d4 , p1 ,   &
                     & p2 , q , w , w2 , w3 , w4 , w6
      integer Kd , m
      w = 0.0d0
      if ( Kd==1 .or. Kd==2 ) w = 2.0d0*m + 1.0d0
      if ( Kd==3 .or. Kd==4 ) w = 2.0d0*m - 1.0d0
      w2 = w*w
      w3 = w*w2
      w4 = w2*w2
      w6 = w2*w4
      d1 = 5.0 + 34.0/w2 + 9.0/w4
      d2 = (33.0+410.0/w2+405.0/w4)/w
      d3 = (63.0+1260.0/w2+2943.0/w4+486.0/w6)/w2
      d4 = (527.0+15617.0/w2+69001.0/w4+41607.0/w6)/w3
      c1 = 128.0
      p2 = q/w4
      p1 = sqrt(p2)
      cv1 = -2.0*q + 2.0*w*sqrt(q) - (w2+1.0)/8.0
      cv2 = (w+3.0/w) + d1/(32.0*p1) + d2/(8.0*c1*p2)
      cv2 = cv2 + d3/(64.0*c1*p1*p2) + d4/(16.0*c1*c1*p2*p2)
      a0 = cv1 - cv2/(c1*p1)
      end



      integer function msta1(x,Mp)
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
      implicit none
      real(wp) a0 , f , f0 , f1 , x
      integer it , Mp , n0 , n1 , nn
      a0 = abs(x)
      n0 = int(1.1d0*a0) + 1
      f0 = envj(n0,a0) - Mp
      n1 = n0 + 5
      f1 = envj(n1,a0) - Mp
      do it = 1 , 20
         nn = n1 - (n1-n0)/(1.0d0-f0/f1)
         f = envj(nn,a0) - Mp
         if ( abs(nn-n1)<1 ) exit
         n0 = n1
         f0 = f1
         n1 = nn
         f1 = f
      enddo
      msta1 = nn
      end


      integer function msta2(x,n,Mp)
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
      implicit none
      real(wp) a0 , ejn , f , f0 , f1 , hmp , obj , x
      integer it , Mp , n , n0 , n1 , nn
      a0 = abs(x)
      hmp = 0.5d0*Mp
      ejn = envj(n,a0)
      if ( ejn<=hmp ) then
         obj = Mp
         n0 = int(1.1*a0) + 1
      else
         obj = hmp + ejn
         n0 = n
      endif
      f0 = envj(n0,a0) - obj
      n1 = n0 + 5
      f1 = envj(n1,a0) - obj
      do it = 1 , 20
         nn = n1 - (n1-n0)/(1.0d0-f0/f1)
         f = envj(nn,a0) - obj
         if ( abs(nn-n1)<1 ) exit
         n0 = n1
         f0 = f1
         n1 = nn
         f1 = f
      enddo
      msta2 = nn + 10
      end

      real(wp) function envj(n,x)
      implicit none
      integer n
      real(wp) x
      envj = 0.5d0*log10(6.28d0*n) - n*log10(1.36d0*x/n)
      end

!       **********************************

      subroutine ittjyb(x,Ttj,Tty)
!
!       ==========================================================
!       Purpose: Integrate [1-J0(t)]/t with respect to t from 0
!                to x, and Y0(t)/t with respect to t from x to ∞
!       Input :  x   --- Variable in the limits  ( x ≥ 0 )
!       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
!                TTY --- Integration of Y0(t)/t from x to ∞
!       ==========================================================
!
      implicit none
      real(wp) e0 , el , f0 , g0 , pi , t , t1 , Ttj , Tty , x ,&
                     & x1 , xt
      pi = 3.141592653589793d0
      el = .5772156649015329d0
      if ( x==0.0d0 ) then
         Ttj = 0.0d0
         Tty = -1.0d+300
      elseif ( x<=4.0d0 ) then
         x1 = x/4.0d0
         t = x1*x1
         Ttj = ((((((.35817d-4*t-.639765d-3)*t+.7092535d-2)*t-          &
             & .055544803d0)*t+.296292677d0)*t-.999999326d0)            &
             & *t+1.999999936d0)*t
         Tty = (((((((-.3546d-5*t+.76217d-4)*t-.1059499d-2)*t+          &
             & .010787555d0)*t-.07810271d0)*t+.377255736d0)             &
             & *t-1.114084491d0)*t+1.909859297d0)*t
         e0 = el + log(x/2.0d0)
         Tty = pi/6.0d0 + e0/pi*(2.0d0*Ttj-e0) - Tty
      elseif ( x<=8.0d0 ) then
         xt = x + .25d0*pi
         t1 = 4.0d0/x
         t = t1*t1
         f0 = (((((.0145369d0*t-.0666297d0)*t+.1341551d0)*t-.1647797d0) &
            & *t+.1608874d0)*t-.2021547d0)*t + .7977506d0
         g0 = ((((((.0160672d0*t-.0759339d0)*t+.1576116d0)*t-.1960154d0)&
            & *t+.1797457d0)*t-.1702778d0)*t+.3235819d0)*t1
         Ttj = (f0*cos(xt)+g0*sin(xt))/(sqrt(x)*x)
         Ttj = Ttj + el + log(x/2.0d0)
         Tty = (f0*sin(xt)-g0*cos(xt))/(sqrt(x)*x)
      else
         t = 8.0d0/x
         xt = x + .25d0*pi
         f0 = (((((.18118d-2*t-.91909d-2)*t+.017033d0)*t-.9394d-3)      &
            & *t-.051445d0)*t-.11d-5)*t + .7978846d0
         g0 = (((((-.23731d-2*t+.59842d-2)*t+.24437d-2)*t-.0233178d0)   &
            & *t+.595d-4)*t+.1620695d0)*t
         Ttj = (f0*cos(xt)+g0*sin(xt))/(sqrt(x)*x)                   &
             & + el + log(x/2.0d0)
         Tty = (f0*sin(xt)-g0*cos(xt))/(sqrt(x)*x)
      endif
      end

!       **********************************

      subroutine ittjya(x,Ttj,Tty)
!
!       =========================================================
!       Purpose: Integrate [1-J0(t)]/t with respect to t from 0
!                to x, and Y0(t)/t with respect to t from x to ∞
!       Input :  x   --- Variable in the limits  ( x ≥ 0 )
!       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
!                TTY --- Integration of Y0(t)/t from x to ∞
!       =========================================================
!
      implicit none
      real(wp) a0 , b1 , bj0 , bj1 , by0 , by1 , e0 , el , g0 , &
                     & g1 , pi , px , qx , r , r0 , r1 , r2 , rs , t ,  &
                     & Ttj
      real(wp) Tty , vt , x , xk
      integer k , l
      pi = 3.141592653589793d0
      el = .5772156649015329d0
      if ( x==0.0d0 ) then
         Ttj = 0.0d0
         Tty = -1.0d+300
      elseif ( x<=20.0d0 ) then
         Ttj = 1.0d0
         r = 1.0d0
         do k = 2 , 100
            r = -.25d0*r*(k-1.0d0)/(k*k*k)*x*x
            Ttj = Ttj + r
            if ( abs(r)<abs(Ttj)*1.0d-12 ) exit
         enddo
         Ttj = Ttj*.125d0*x*x
         e0 = .5d0*(pi*pi/6.0d0-el*el) - (.5d0*log(x/2.0d0)+el)        &
            & *log(x/2.0d0)
         b1 = el + log(x/2.0d0) - 1.5d0
         rs = 1.0d0
         r = -1.0d0
         do k = 2 , 100
            r = -.25d0*r*(k-1.0d0)/(k*k*k)*x*x
            rs = rs + 1.0d0/k
            r2 = r*(rs+1.0d0/(2.0d0*k)-(el+log(x/2.0d0)))
            b1 = b1 + r2
            if ( abs(r2)<abs(b1)*1.0d-12 ) exit
         enddo
         Tty = 2.0d0/pi*(e0+.125d0*x*x*b1)
      else
         a0 = sqrt(2.0d0/(pi*x))
         bj0 = 0.0d0
         by0 = 0.0d0
         bj1 = 0.0d0
         do l = 0 , 1
            vt = 4.0d0*l*l
            px = 1.0d0
            r = 1.0d0
            do k = 1 , 14
               r = -.0078125d0*r*(vt-(4.0d0*k-3.0d0)**2)/(x*k)          &
                 & *(vt-(4.0d0*k-1.0d0)**2)/((2.0d0*k-1.0d0)*x)
               px = px + r
               if ( abs(r)<abs(px)*1.0d-12 ) exit
            enddo
            qx = 1.0d0
            r = 1.0d0
            do k = 1 , 14
               r = -.0078125d0*r*(vt-(4.0d0*k-1.0d0)**2)/(x*k)          &
                 & *(vt-(4.0d0*k+1.0d0)**2)/(2.0d0*k+1.0d0)/x
               qx = qx + r
               if ( abs(r)<abs(qx)*1.0d-12 ) exit
            enddo
            qx = .125d0*(vt-1.0d0)/x*qx
            xk = x - (.25d0+.5d0*l)*pi
            bj1 = a0*(px*cos(xk)-qx*sin(xk))
            by1 = a0*(px*sin(xk)+qx*cos(xk))
            if ( l==0 ) then
               bj0 = bj1
               by0 = by1
            endif
         enddo
         t = 2.0d0/x
         g0 = 1.0d0
         r0 = 1.0d0
         do k = 1 , 10
            r0 = -k*k*t*t*r0
            g0 = g0 + r0
         enddo
         g1 = 1.0d0
         r1 = 1.0d0
         do k = 1 , 10
            r1 = -k*(k+1.0d0)*t*t*r1
            g1 = g1 + r1
         enddo
         Ttj = 2.0d0*g1*bj0/(x*x) - g0*bj1/x + el + log(x/2.0d0)
         Tty = 2.0d0*g1*by0/(x*x) - g0*by1/x
      endif
      end

!       **********************************

      subroutine cjylv(v,z,Cbjv,Cdjv,Cbyv,Cdyv)
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
      implicit none
      real(wp) a , pi , v , v0 , vr
      complex(wp) Cbjv , Cbyv , Cdjv , Cdyv , ceta , cf , cfj , cfy ,    &
               & csj , csy , ct , ct2 , cws , z
      integer i , k , km , l , l0 , lf
      dimension cf(12) , a(91)
      km = 12
      call cjk(km,a)
      pi = 3.141592653589793d0
      do l = 1 , 0 , -1
         v0 = v - l
         cws = sqrt(1.0d0-(z/v0)*(z/v0))
         ceta = cws + log(z/v0/(1.0d0+cws))
         ct = 1.0d0/cws
         ct2 = ct*ct
         do k = 1 , km
            l0 = k*(k+1)/2 + 1
            lf = l0 + k
            cf(k) = a(lf)
            do i = lf - 1 , l0 , -1
               cf(k) = cf(k)*ct2 + a(i)
            enddo
            cf(k) = cf(k)*ct**k
         enddo
         vr = 1.0d0/v0
         csj = (1.0d0,0.0d0)
         do k = 1 , km
            csj = csj + cf(k)*vr**k
         enddo
         Cbjv = sqrt(ct/(2.0d0*pi*v0))*exp(v0*ceta)*csj
         if ( l==1 ) cfj = Cbjv
         csy = (1.0d0,0.0d0)
         do k = 1 , km
            csy = csy + (-1)**k*cf(k)*vr**k
         enddo
         Cbyv = -sqrt(2.0d0*ct/(pi*v0))*exp(-v0*ceta)*csy
         if ( l==1 ) cfy = Cbyv
      enddo
      Cdjv = -v/z*Cbjv + cfj
      Cdyv = -v/z*Cbyv + cfy
      end



!       **********************************

      subroutine rmn2l(m,n,c,x,Df,Kd,R2f,R2d,Id)
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
      implicit none
      real(wp) a0 , b0 , c , cx , Df , dy , eps , eps1 , eps2 , &
                     & r , r0 , R2d , R2f , reg , suc , sud , sw , sy , &
                     & x
      integer Id , id1 , id2 , ip , j , k , Kd , l , lg , m , n , nm ,  &
            & nm1 , nm2 , np
      dimension Df(200) , sy(0:251) , dy(0:251)
      eps = 1.0d-14
      ip = 1
      nm1 = int((n-m)/2)
      if ( n-m==2*nm1 ) ip = 0
      nm = 25 + nm1 + int(c)
      reg = 1.0d0
      if ( m+nm>80 ) reg = 1.0d-200
      nm2 = 2*nm + m
      cx = c*x
      call sphy(nm2,cx,nm2,sy,dy)
      r0 = reg
      do j = 1 , 2*m + ip
         r0 = r0*j
      enddo
      r = r0
      suc = r*Df(1)
      sw = 0.0d0
      do k = 2 , nm
         r = r*(m+k-1.0)*(m+k+ip-1.5d0)/(k-1.0d0)/(k+ip-1.5d0)
         suc = suc + r*Df(k)
         if ( k>nm1 .and. abs(suc-sw)<abs(suc)*eps ) exit
         sw = suc
      enddo
      a0 = (1.0d0-Kd/(x*x))**(0.5d0*m)/suc
      R2f = 0.0d0
      eps1 = 0.0d0
      np = 0
      do k = 1 , nm
         l = 2*k + m - n - 2 + ip
         lg = 1
         if ( l/=4*int(l/4) ) lg = -1
         if ( k==1 ) then
            r = r0
         else
            r = r*(m+k-1.0)*(m+k+ip-1.5d0)/(k-1.0d0)/(k+ip-1.5d0)
         endif
         np = m + 2*k - 2 + ip
         R2f = R2f + lg*r*(Df(k)*sy(np))
         eps1 = abs(R2f-sw)
         if ( k>nm1 .and. eps1<abs(R2f)*eps ) exit
         sw = R2f
      enddo
      id1 = int(log10(eps1/abs(R2f)+eps))
      R2f = R2f*a0
      if ( np>=nm2 ) then
         Id = 10
         return
      endif
      b0 = Kd*m/x**3.0d0/(1.0-Kd/(x*x))*R2f
      sud = 0.0d0
      eps2 = 0.0d0
      do k = 1 , nm
         l = 2*k + m - n - 2 + ip
         lg = 1
         if ( l/=4*int(l/4) ) lg = -1
         if ( k==1 ) then
            r = r0
         else
            r = r*(m+k-1.0)*(m+k+ip-1.5d0)/(k-1.0d0)/(k+ip-1.5d0)
         endif
         np = m + 2*k - 2 + ip
         sud = sud + lg*r*(Df(k)*dy(np))
         eps2 = abs(sud-sw)
         if ( k>nm1 .and. eps2<abs(sud)*eps ) exit
         sw = sud
      enddo
      R2d = b0 + a0*c*sud
      id2 = int(log10(eps2/abs(sud)+eps))
      Id = max(id1,id2)
      end



!       **********************************

      subroutine psi_spec(x,Ps)
!
!       ======================================
!       Purpose: Compute Psi function
!       Input :  x  --- Argument of psi(x)
!       Output:  PS --- psi(x)
!       ======================================
!
      implicit none
      real(wp) a1 , a2 , a3 , a4 , a5 , a6 , a7 , a8 , el , pi ,&
                     & Ps , s , x , x2 , xa
      integer k , n
      xa = abs(x)
      pi = 3.141592653589793d0
      el = .5772156649015329d0
      s = 0.0d0
      if ( x==int(x) .and. x<=0.0 ) then
         Ps = 1.0d+300
         return
      elseif ( xa==int(xa) ) then
         n = xa
         do k = 1 , n - 1
            s = s + 1.0d0/k
         enddo
         Ps = -el + s
      elseif ( xa+.5==int(xa+.5) ) then
         n = xa - .5
         do k = 1 , n
            s = s + 1.0/(2.0d0*k-1.0d0)
         enddo
         Ps = -el + 2.0d0*s - 1.386294361119891d0
      else
         if ( xa<10.0 ) then
            n = 10 - int(xa)
            do k = 0 , n - 1
               s = s + 1.0d0/(xa+k)
            enddo
            xa = xa + n
         endif
         x2 = 1.0d0/(xa*xa)
         a1 = -.8333333333333d-01
         a2 = .83333333333333333d-02
         a3 = -.39682539682539683d-02
         a4 = .41666666666666667d-02
         a5 = -.75757575757575758d-02
         a6 = .21092796092796093d-01
         a7 = -.83333333333333333d-01
         a8 = .4432598039215686d0
         Ps = log(xa) - .5d0/xa +                                      &
            & x2*(((((((a8*x2+a7)*x2+a6)*x2+a5)*x2+a4)*x2+a3)*x2+a2)    &
            & *x2+a1)
         Ps = Ps - s
      endif
      if ( x<0.0 ) Ps = Ps - pi*cos(pi*x)/sin(pi*x) - 1.0d0/x
      end

!       **********************************

      subroutine cva2(Kd,m,q,a)
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
      implicit none
      real(wp) a , a1 , a2 , delta , q , q1 , q2 , qq
      integer i , iflag , Kd , m , ndiv , nn
      if ( m<=12 .or. q<=3.0*m .or. q>m*m ) then
         call cv0(Kd,m,q,a)
         if ( q/=0.0d0 .and. m/=2 ) call refine(Kd,m,q,a)
         if ( q>2.0d-3 .and. m==2 ) call refine(Kd,m,q,a)
      else
         ndiv = 10
         delta = (m-3.0)*m/ndiv
         if ( (q-3.0*m)<=(m*m-q) ) then
 20         nn = int((q-3.0*m)/delta) + 1
            delta = (q-3.0*m)/nn
            q1 = 2.0*m
            call cvqm(m,q1,a1)
            q2 = 3.0*m
            call cvqm(m,q2,a2)
            qq = 3.0*m
            do i = 1 , nn
               qq = qq + delta
               a = (a1*q2-a2*q1+(a2-a1)*qq)/(q2-q1)
               iflag = 1
               if ( i==nn ) iflag = -1
               call refine(Kd,m,qq,a)
               q1 = q2
               q2 = qq
               a1 = a2
               a2 = a
            enddo
            if ( iflag==-10 ) then
               ndiv = ndiv*2
               delta = (m-3.0)*m/ndiv
               goto 20
            endif
         else
 40         nn = int((m*m-q)/delta) + 1
            delta = (m*m-q)/nn
            q1 = m*(m-1.0)
            call cvql(Kd,m,q1,a1)
            q2 = m*m
            call cvql(Kd,m,q2,a2)
            qq = m*m
            do i = 1 , nn
               qq = qq - delta
               a = (a1*q2-a2*q1+(a2-a1)*qq)/(q2-q1)
               iflag = 1
               if ( i==nn ) iflag = -1
               call refine(Kd,m,qq,a)
               q1 = q2
               q2 = qq
               a1 = a2
               a2 = a
            enddo
            if ( iflag==-10 ) then
               ndiv = ndiv*2
               delta = (m-3.0)*m/ndiv
               goto 40
            endif
         endif
      endif
      end



!       **********************************

      subroutine lpmns(m,n,x,Pm,Pd)
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
      implicit none
      integer k , m , n
      real(wp) Pd , Pm , pm0 , pm1 , pm2 , pmk , x , x0
      dimension Pm(0:n) , Pd(0:n)
      do k = 0 , n
         Pm(k) = 0.0d0
         Pd(k) = 0.0d0
      enddo
      if ( abs(x)==1.0d0 ) then
         do k = 0 , n
            if ( m==0 ) then
               Pm(k) = 1.0d0
               Pd(k) = 0.5d0*k*(k+1.0)
               if ( x<0.0 ) then
                  Pm(k) = (-1)**k*Pm(k)
                  Pd(k) = (-1)**(k+1)*Pd(k)
               endif
            elseif ( m==1 ) then
               Pd(k) = 1.0d+300
            elseif ( m==2 ) then
               Pd(k) = -0.25d0*(k+2.0)*(k+1.0)*k*(k-1.0)
               if ( x<0.0 ) Pd(k) = (-1)**(k+1)*Pd(k)
            endif
         enddo
         return
      endif
      x0 = abs(1.0d0-x*x)
      pm0 = 1.0d0
      pmk = pm0
      do k = 1 , m
         pmk = (2.0d0*k-1.0d0)*sqrt(x0)*pm0
         pm0 = pmk
      enddo
      pm1 = (2.0d0*m+1.0d0)*x*pm0
      Pm(m) = pmk
      Pm(m+1) = pm1
      do k = m + 2 , n
         pm2 = ((2.0d0*k-1.0d0)*x*pm1-(k+m-1.0d0)*pmk)/(k-m)
         Pm(k) = pm2
         pmk = pm1
         pm1 = pm2
      enddo
      Pd(0) = ((1.0d0-m)*Pm(1)-x*Pm(0))/(x*x-1.0)
      do k = 1 , n
         Pd(k) = (k*x*Pm(k)-(k+m)*Pm(k-1))/(x*x-1.0d0)
      enddo
      do k = 1 , n
         Pm(k) = (-1)**m*Pm(k)
         Pd(k) = (-1)**m*Pd(k)
      enddo
      end

!       **********************************

      subroutine cerf(z,Cer,Cder)
!
!       ==========================================================
!       Purpose: Compute complex Error function erf(z) & erf'(z)
!       Input:   z   --- Complex argument of erf(z)
!                x   --- Real part of z
!                y   --- Imaginary part of z
!       Output:  CER --- erf(z)
!                CDER --- erf'(z)
!       ==========================================================
      implicit none
      real(wp) c0 , cs , ei1 , ei2 , eps , er , er0 , er1 ,     &
                     & er2 , eri , err , pi , r , ss , w , w1 , w2 , x ,&
                     & x2 , y
      integer k , n
      complex(wp) z , Cer , Cder
      eps = 1.0d-12
      pi = 3.141592653589793d0
      x = dble(z)
      y = dimag(z)
      x2 = x*x
      if ( x<=3.5d0 ) then
         er = 1.0d0
         r = 1.0d0
         w = 0.0d0
         do k = 1 , 100
            r = r*x2/(k+0.5d0)
            er = er + r
            if ( abs(er-w)<=eps*abs(er) ) exit
            w = er
         enddo
         c0 = 2.0d0/sqrt(pi)*x*exp(-x2)
         er0 = c0*er
      else
         er = 1.0d0
         r = 1.0d0
         do k = 1 , 12
            r = -r*(k-0.5d0)/x2
            er = er + r
         enddo
         c0 = exp(-x2)/(x*sqrt(pi))
         er0 = 1.0d0 - c0*er
      endif
      if ( y==0.0d0 ) then
         err = er0
         eri = 0.0d0
      else
         cs = cos(2.0d0*x*y)
         ss = sin(2.0d0*x*y)
         er1 = exp(-x2)*(1.0d0-cs)/(2.0d0*pi*x)
         ei1 = exp(-x2)*ss/(2.0d0*pi*x)
         er2 = 0.0d0
         w1 = 0.0d0
         do n = 1 , 100
            er2 = er2 + exp(-.25d0*n*n)/(n*n+4.0d0*x2)                  &
                & *(2.0d0*x-2.0d0*x*dcosh(n*y)*cs+n*dsinh(n*y)*ss)
            if ( abs((er2-w1)/er2)<eps ) exit
            w1 = er2
         enddo
         c0 = 2.0d0*exp(-x2)/pi
         err = er0 + er1 + c0*er2
         ei2 = 0.0d0
         w2 = 0.0d0
         do n = 1 , 100
            ei2 = ei2 + exp(-.25d0*n*n)/(n*n+4.0d0*x2)                  &
                & *(2.0d0*x*dcosh(n*y)*ss+n*dsinh(n*y)*cs)
            if ( abs((ei2-w2)/ei2)<eps ) exit
            w2 = ei2
         enddo
         eri = ei1 + c0*ei2
      endif
      Cer = dcmplx(err,eri)
      Cder = 2.0d0/sqrt(pi)*exp(-z*z)
      end

!       **********************************

      subroutine rswfp(m,n,c,x,Cv,Kf,R1f,R1d,R2f,R2d)
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
      implicit none
      real(wp) c , Cv , df , R1d , R1f , R2d , R2f , x
      integer id , kd , Kf , m , n
      dimension df(200)
      kd = 1
      call sdmn(m,n,c,Cv,kd,df)
      if ( Kf/=2 ) call rmn1(m,n,c,x,df,kd,R1f,R1d)
      if ( Kf>1 ) then
         call rmn2l(m,n,c,x,df,kd,R2f,R2d,id)
         if ( id>-8 ) call rmn2sp(m,n,c,x,Cv,df,kd,R2f,R2d)
      endif
      end



!       **********************************

      subroutine jyndd(n,x,Bjn,Djn,Fjn,Byn,Dyn,Fyn)
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
      implicit none
      real(wp) bj , Bjn , by , Byn , Djn , Dyn , Fjn , Fyn , x
      integer n , nm
      dimension bj(2) , by(2)
      call jynbh(n+1,n,x,nm,bj,by)
!       Compute derivatives by differentiation formulas
      Bjn = bj(1)
      Byn = by(1)
      Djn = -bj(2) + n*bj(1)/x
      Dyn = -by(2) + n*by(1)/x
      Fjn = (n*n/(x*x)-1.0d0)*Bjn - Djn/x
      Fyn = (n*n/(x*x)-1.0d0)*Byn - Dyn/x
      end


!       **********************************

      subroutine gam0(x,Ga)
!
!       ================================================
!       Purpose: Compute gamma function Г(x)
!       Input :  x  --- Argument of Г(x)  ( |x| ≤ 1 )
!       Output:  GA --- Г(x)
!       ================================================
!
      implicit none
      real(wp) g , Ga , gr , x
      integer k
      dimension g(25)
      data g/1.0d0 , 0.5772156649015329d0 , -0.6558780715202538d0 ,     &
         & -0.420026350340952d-1 , 0.1665386113822915d0 ,               &
         & -.421977345555443d-1 , -.96219715278770d-2 ,                 &
         & .72189432466630d-2 , -.11651675918591d-2 ,                   &
         & -.2152416741149d-3 , .1280502823882d-3 , -.201348547807d-4 , &
         & -.12504934821d-5 , .11330272320d-5 , -.2056338417d-6 ,       &
         & .61160950d-8 , .50020075d-8 , -.11812746d-8 , .1043427d-9 ,  &
         & .77823d-11 , -.36968d-11 , .51d-12 , -.206d-13 , -.54d-14 ,  &
         & .14d-14/
      gr = (25)
      do k = 24 , 1 , -1
         gr = gr*x + g(k)
      enddo
      Ga = 1.0d0/(gr*x)
      end


!       **********************************

      subroutine cisib(x,Ci,Si)
!
!       =============================================
!       Purpose: Compute cosine and sine integrals
!                Si(x) and Ci(x) ( x ≥ 0 )
!       Input :  x  --- Argument of Ci(x) and Si(x)
!       Output:  CI --- Ci(x)
!                SI --- Si(x)
!       =============================================
!
      implicit none
      real(wp) Ci , fx , gx , Si , x , x2
      x2 = x*x
      if ( x==0.0 ) then
         Ci = -1.0d+300
         Si = 0.0d0
      elseif ( x<=1.0d0 ) then
         Ci = ((((-3.0d-8*x2+3.10d-6)*x2-2.3148d-4)*x2+1.041667d-2)     &
            & *x2-0.25)*x2 + 0.577215665d0 + log(x)
         Si = ((((3.1d-7*x2-2.834d-5)*x2+1.66667d-003)*x2-5.555556d-002)&
            & *x2+1.0)*x
      else
         fx = ((((x2+38.027264d0)*x2+265.187033d0)*x2+335.67732d0)      &
            & *x2+38.102495d0)                                          &
            & /((((x2+40.021433d0)*x2+322.624911d0)*x2+570.23628d0)     &
            & *x2+157.105423d0)
         gx = ((((x2+42.242855d0)*x2+302.757865d0)*x2+352.018498d0)     &
            & *x2+21.821899d0)                                          &
            & /((((x2+48.196927d0)*x2+482.485984d0)*x2+1114.978885d0)   &
            & *x2+449.690326d0)/x
         Ci = fx*sin(x)/x - gx*cos(x)/x
         Si = 1.570796327d0 - fx*cos(x)/x - gx*sin(x)/x
      endif
      end

!       **********************************

      subroutine eulera(n,En)
!
!       ======================================
!       Purpose: Compute Euler number En
!       Input :  n --- Serial number
!       Output:  EN(n) --- En
!       ======================================
!
      implicit none
      real(wp) En , r , s
      integer j , k , m , n
      dimension En(0:n)
      En(0) = 1.0d0
      do m = 1 , n/2
         s = 1.0d0
         do k = 1 , m - 1
            r = 1.0d0
            do j = 1 , 2*k
               r = r*(2.0d0*m-2.0d0*k+j)/j
            enddo
            s = s + r*En(2*k)
         enddo
         En(2*m) = -s
      enddo
      end

!       **********************************

      subroutine refine(Kd,m,q,a)
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
      implicit none
      real(wp) a , ca , delta , eps , f , f0 , f1 , q , x , x0 ,&
                     & x1
      integer it , Kd , m , mj
      eps = 1.0d-14
      mj = 10 + m
      ca = a
      delta = 0.0d0
      x0 = a
      call cvf(Kd,m,q,x0,mj,f0)
      x1 = 1.002*a
      call cvf(Kd,m,q,x1,mj,f1)
      do it = 1 , 100
         mj = mj + 1
         x = x1 - (x1-x0)/(1.0d0-f0/f1)
         call cvf(Kd,m,q,x,mj,f)
         if ( abs(1.0-x1/x)<eps .or. f==0.0 ) exit
         x0 = x1
         f0 = f1
         x1 = x
         f1 = f
      enddo
      a = x
      end



!       **********************************

      subroutine cisia(x,Ci,Si)
!
!       =============================================
!       Purpose: Compute cosine and sine integrals
!                Si(x) and Ci(x)  ( x ≥ 0 )
!       Input :  x  --- Argument of Ci(x) and Si(x)
!       Output:  CI --- Ci(x)
!                SI --- Si(x)
!       =============================================
!
      implicit none
      real(wp) bj , Ci , el , eps , p2 , Si , x , x2 , xa ,     &
                     & xa0 , xa1 , xcs , xf , xg , xg1 , xg2 , xr , xs ,&
                     & xss
      integer k , m
      dimension bj(101)
      p2 = 1.570796326794897d0
      el = .5772156649015329d0
      eps = 1.0d-15
      x2 = x*x
      if ( x==0.0d0 ) then
         Ci = -1.0d+300
         Si = 0.0d0
      elseif ( x<=16.0d0 ) then
         xr = -.25d0*x2
         Ci = el + log(x) + xr
         do k = 2 , 40
            xr = -.5d0*xr*(k-1)/(k*k*(2*k-1))*x2
            Ci = Ci + xr
            if ( abs(xr)<abs(Ci)*eps ) exit
         enddo
         xr = x
         Si = x
         do k = 1 , 40
            xr = -.5d0*xr*(2*k-1)/k/(4*k*k+4*k+1)*x2
            Si = Si + xr
            if ( abs(xr)<abs(Si)*eps ) return
         enddo
      elseif ( x<=32.0d0 ) then
         m = int(47.2+.82*x)
         xa1 = 0.0d0
         xa0 = 1.0d-100
         do k = m , 1 , -1
            xa = 4.0d0*k*xa0/x - xa1
            bj(k) = xa
            xa1 = xa0
            xa0 = xa
         enddo
         xs = bj(1)
         do k = 3 , m , 2
            xs = xs + 2.0d0*bj(k)
         enddo
         bj(1) = bj(1)/xs
         do k = 2 , m
            bj(k) = bj(k)/xs
         enddo
         xr = 1.0d0
         xg1 = bj(1)
         do k = 2 , m
            xr = .25d0*xr*(2.0*k-3.0)**2/((k-1.0)*(2.0*k-1.0)**2)*x
            xg1 = xg1 + bj(k)*xr
         enddo
         xr = 1.0d0
         xg2 = bj(1)
         do k = 2 , m
            xr = .25d0*xr*(2.0*k-5.0)**2/((k-1.0)*(2.0*k-3.0)**2)*x
            xg2 = xg2 + bj(k)*xr
         enddo
         xcs = cos(x/2.0d0)
         xss = sin(x/2.0d0)
         Ci = el + log(x) - x*xss*xg1 + 2*xcs*xg2 - 2*xcs*xcs
         Si = x*xcs*xg1 + 2*xss*xg2 - sin(x)
      else
         xr = 1.0d0
         xf = 1.0d0
         do k = 1 , 9
            xr = -2.0d0*xr*k*(2*k-1)/x2
            xf = xf + xr
         enddo
         xr = 1.0d0/x
         xg = xr
         do k = 1 , 8
            xr = -2.0d0*xr*(2*k+1)*k/x2
            xg = xg + xr
         enddo
         Ci = xf*sin(x)/x - xg*cos(x)/x
         Si = p2 - xf*cos(x)/x - xg*sin(x)/x
      endif
      end


!       **********************************

      subroutine itsl0(x,Tl0)
!
!       ===========================================================
!       Purpose: Evaluate the integral of modified Struve function
!                L0(t) with respect to t from 0 to x
!       Input :  x   --- Upper limit  ( x ≥ 0 )
!       Output:  TL0 --- Integration of L0(t) from 0 to x
!       ===========================================================
!
      implicit none
      real(wp) a , a0 , a1 , af , el , pi , r , rd , s , s0 ,   &
                     & ti , Tl0 , x
      integer k
      dimension a(18)
      pi = 3.141592653589793d0
      r = 1.0d0
      if ( x<=20.0 ) then
         s = 0.5d0
         do k = 1 , 100
            rd = 1.0d0
            if ( k==1 ) rd = 0.5d0
            r = r*rd*k/(k+1.0d0)*(x/(2.0d0*k+1.0d0))**2
            s = s + r
            if ( abs(r/s)<1.0d-12 ) exit
         enddo
         Tl0 = 2.0d0/pi*x*x*s
      else
         s = 1.0d0
         do k = 1 , 10
            r = r*k/(k+1.0d0)*((2.0d0*k+1.0d0)/x)**2
            s = s + r
            if ( abs(r/s)<1.0d-12 ) exit
         enddo
         el = .57721566490153d0
         s0 = -s/(pi*x*x) + 2.0d0/pi*(log(2.0d0*x)+el)
         a0 = 1.0d0
         a1 = 5.0d0/8.0d0
         a(1) = a1
         do k = 1 , 10
            af = ((1.5d0*(k+.50d0)*(k+5.0d0/6.0d0)*a1-.5d0*(k+.5d0)     &
               & **2*(k-.5d0)*a0))/(k+1.0d0)
            a(k+1) = af
            a0 = a1
            a1 = af
         enddo
         ti = 1.0d0
         r = 1.0d0
         do k = 1 , 11
            r = r/x
            ti = ti + a(k)*r
         enddo
         Tl0 = ti/sqrt(2*pi*x)*exp(x) + s0
      endif
      end

!       **********************************

      subroutine clqn(n,x,y,Cqn,Cqd)
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
      implicit none
      complex(wp) cq0 , cq1 , Cqd , cqf0 , cqf1 , cqf2 , Cqn , z
      integer k , km , ls , n
      real(wp) x , y
      dimension Cqn(0:n) , Cqd(0:n)
      z = dcmplx(x,y)
      if ( z==1.0d0 ) then
         do k = 0 , n
            Cqn(k) = (1.0d+300,0.0d0)
            Cqd(k) = (1.0d+300,0.0d0)
         enddo
         return
      endif
      ls = 1
      if ( abs(z)>1.0d0 ) ls = -1
      cq0 = 0.5d0*log(ls*(1.0d0+z)/(1.0d0-z))
      cq1 = z*cq0 - 1.0d0
      Cqn(0) = cq0
      Cqn(1) = cq1
      if ( abs(z)<1.0001d0 ) then
         cqf0 = cq0
         cqf1 = cq1
         do k = 2 , n
            cqf2 = ((2.0d0*k-1.0d0)*z*cqf1-(k-1.0d0)*cqf0)/k
            Cqn(k) = cqf2
            cqf0 = cqf1
            cqf1 = cqf2
         enddo
      else
         if ( abs(z)>1.1d0 ) then
            km = 40 + n
         else
            km = (40+n)*int(-1.0-1.8*log(abs(z-1.0)))
         endif
         cqf2 = 0.0d0
         cqf1 = 1.0d0
         do k = km , 0 , -1
            cqf0 = ((2*k+3.0d0)*z*cqf1-(k+2.0d0)*cqf2)/(k+1.0d0)
            if ( k<=n ) Cqn(k) = cqf0
            cqf2 = cqf1
            cqf1 = cqf0
         enddo
         do k = 0 , n
            Cqn(k) = Cqn(k)*cq0/cqf0
         enddo
      endif
      Cqd(0) = (Cqn(1)-z*Cqn(0))/(z*z-1.0d0)
      do k = 1 , n
         Cqd(k) = (k*z*Cqn(k)-k*Cqn(k-1))/(z*z-1.0d0)
      enddo
      end

!       **********************************

      subroutine airyzo(Nt,Kf,Xa,Xb,Xc,Xd)
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
      implicit none
      real(wp) ad , ai , bd , bi , err , pi , rt , rt0 , u ,    &
                     & u1 , x , Xa , Xb , Xc , Xd
      integer i , Kf , Nt
      dimension Xa(Nt) , Xb(Nt) , Xc(Nt) , Xd(Nt)
      pi = 3.141592653589793d0
      rt = 0.0d0
      do i = 1 , Nt
         rt0 = 0d0
         if ( Kf==1 ) then
            u = 3.0d0*pi*(4.0d0*i-1)/8.0d0
            u1 = 1/(u*u)
         elseif ( Kf==2 ) then
            if ( i==1 ) then
               rt0 = -1.17371d0
            else
               u = 3.0d0*pi*(4.0d0*i-3.0d0)/8.0d0
               u1 = 1/(u*u)
            endif
         endif
!             DLMF 9.9.18
         if ( rt0==0 ) rt0 = -(u*u)**(1.0d0/3.0d0)                      &
                           & *(+1d0+u1*(5d0/48d0+u1*                    &
                           & (-5d0/36d0+u1*(77125d0/82944d0+            &
                           & u1*(-108056875d0/6967296d0)))))
 50      x = rt0
         call airyb(x,ai,bi,ad,bd)
         if ( Kf==1 ) rt = rt0 - ai/ad
         if ( Kf==2 ) rt = rt0 - bi/bd
         err = abs((rt-rt0)/rt)
         if ( err>1.d-12 ) then
            rt0 = rt
            goto 50
         else
            Xa(i) = rt
            if ( err>1d-14 ) call airyb(rt,ai,bi,ad,bd)
            if ( Kf==1 ) Xd(i) = ad
            if ( Kf==2 ) Xd(i) = bd
         endif
      enddo
      do i = 1 , Nt
         rt0 = 0d0
         if ( Kf==1 ) then
            if ( i==1 ) then
               rt0 = -1.01879d0
            else
               u = 3.0d0*pi*(4.0d0*i-3.0d0)/8.0d0
               u1 = 1/(u*u)
            endif
         elseif ( Kf==2 ) then
            if ( i==1 ) then
               rt0 = -2.29444d0
            else
               u = 3.0d0*pi*(4.0d0*i-1.0d0)/8.0d0
               u1 = 1/(u*u)
            endif
         endif
!             DLMF 9.9.19
         if ( rt0==0 ) rt0 = -(u*u)**(1.0d0/3.0d0)                      &
                           & *(+1d0+u1*(-7d0/48d0+u1*                   &
                           & (+35d0/288d0+u1*(-181223d0/207360d0+       &
                           & u1*(18683371d0/1244160d0)))))
 100     x = rt0
         call airyb(x,ai,bi,ad,bd)
         if ( Kf==1 ) rt = rt0 - ad/(ai*x)
         if ( Kf==2 ) rt = rt0 - bd/(bi*x)
         err = abs((rt-rt0)/rt)
         if ( err>1.0d-12 ) then
            rt0 = rt
            goto 100
         else
            Xb(i) = rt
            if ( err>1d-14 ) call airyb(rt,ai,bi,ad,bd)
            if ( Kf==1 ) Xc(i) = ai
            if ( Kf==2 ) Xc(i) = bi
         endif
      enddo
      end



!       **********************************

      subroutine error(x,Err)
!
!       =========================================
!       Purpose: Compute error function erf(x)
!       Input:   x   --- Argument of erf(x)
!       Output:  ERR --- erf(x)
!       =========================================
!
      implicit none
      real(wp) c0 , eps , er , Err , pi , r , x , x2
      integer k
      eps = 1.0d-15
      pi = 3.141592653589793d0
      x2 = x*x
      if ( abs(x)<3.5d0 ) then
         er = 1.0d0
         r = 1.0d0
         do k = 1 , 50
            r = r*x2/(k+0.5d0)
            er = er + r
            if ( abs(r)<=abs(er)*eps ) exit
         enddo
         c0 = 2.0d0/sqrt(pi)*x*exp(-x2)
         Err = c0*er
      else
         er = 1.0d0
         r = 1.0d0
         do k = 1 , 12
            r = -r*(k-0.5d0)/x2
            er = er + r
         enddo
         c0 = exp(-x2)/(abs(x)*sqrt(pi))
         Err = 1.0d0 - c0*er
         if ( x<0.0 ) Err = -Err
      endif
      end

!       **********************************

      subroutine cerror(z,Cer)
!
!       ====================================================
!       Purpose: Compute error function erf(z) for a complex
!                argument (z=x+iy)
!       Input :  z   --- Complex argument
!       Output:  CER --- erf(z)
!       ====================================================
!
      implicit none
      complex(wp) c0 , Cer , cl , cr , cs , z , z1
      integer k
      real(wp) a0 , pi
      a0 = abs(z)
      c0 = exp(-z*z)
      pi = 3.141592653589793d0
      z1 = z
      if ( dble(z)<0.0 ) z1 = -z
!
!       Cutoff radius R = 4.36; determined by balancing rounding error
!       and asymptotic expansion error, see below.
!
!       The resulting maximum global accuracy expected is around 1e-8
!
      if ( a0<=4.36d0 ) then
!
!          Rounding error in the Taylor expansion is roughly
!
!          ~ R*R * EPSILON * R**(2 R**2) / (2 R**2 Gamma(R**2 + 1/2))
!
         cs = z1
         cr = z1
         do k = 1 , 120
            cr = cr*z1*z1/(k+0.5d0)
            cs = cs + cr
            if ( abs(cr/cs)<1.0d-15 ) exit
         enddo
         Cer = 2.0d0*c0*cs/sqrt(pi)
      else
         cl = 1.0d0/z1
         cr = cl
!
!          Asymptotic series; maximum K must be at most ~ R^2.
!
!          The maximum accuracy obtainable from this expansion is roughly
!
!          ~ Gamma(2R**2 + 2) / (
!                   (2 R**2)**(R**2 + 1/2) Gamma(R**2 + 3/2) 2**(R**2 + 1/2))
!
         do k = 1 , 20
            cr = -cr*(k-0.5d0)/(z1*z1)
            cl = cl + cr
            if ( abs(cr/cl)<1.0d-15 ) exit
         enddo
         Cer = 1.0d0 - c0*cl/sqrt(pi)
      endif
      if ( dble(z)<0.0 ) Cer = -Cer
      end



!       **********************************

      subroutine eulerb(n,En)
!
!       ======================================
!       Purpose: Compute Euler number En
!       Input :  n --- Serial number
!       Output:  EN(n) --- En
!       ======================================
!
      implicit none
      real(wp) En , hpi , r1 , r2 , s
      integer isgn , k , m , n
      dimension En(0:n)
      hpi = 2.0d0/3.141592653589793d0
      En(0) = 1.0d0
      En(2) = -1.0d0
      r1 = -4.0d0*hpi**3
      do m = 4 , n , 2
         r1 = -r1*(m-1)*m*hpi*hpi
         r2 = 1.0d0
         isgn = 1.0d0
         do k = 3 , 1000 , 2
            isgn = -isgn
            s = (1.0d0/k)**(m+1)
            r2 = r2 + isgn*s
            if ( s<1.0d-15 ) exit
         enddo
         En(m) = r1*r2
      enddo
      end

!       **********************************

      subroutine cva1(Kd,m,q,Cv)
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
      implicit none
      real(wp) Cv , d , e , eps , f , g , h , q , s , t , t1 ,  &
                     & x1 , xa , xb
      integer i , ic , icm , j , k , k1 , Kd , m , nm , nm1
      dimension g(200) , h(200) , d(500) , e(500) , f(500) , Cv(200)
      eps = 1.0d-14
      icm = int(m/2) + 1
      if ( Kd==4 ) icm = m/2
      if ( q/=0.0d0 ) then
         nm = int(10+1.5*m+0.5*q)
         e(1) = 0.0d0
         f(1) = 0.0d0
         if ( Kd==1 ) then
            d(1) = 0.0d0
            do i = 2 , nm
               d(i) = 4.0d0*(i-1.0d0)**2
               e(i) = q
               f(i) = q*q
            enddo
            e(2) = sqrt(2.0d0)*q
            f(2) = 2.0d0*q*q
         elseif ( Kd/=4 ) then
            d(1) = 1.0d0 + (-1)**Kd*q
            do i = 2 , nm
               d(i) = (2.0d0*i-1.0d0)**2
               e(i) = q
               f(i) = q*q
            enddo
         else
            d(1) = 4.0d0
            do i = 2 , nm
               d(i) = 4.0d0*i*i
               e(i) = q
               f(i) = q*q
            enddo
         endif
         xa = d(nm) + abs(e(nm))
         xb = d(nm) - abs(e(nm))
         nm1 = nm - 1
         do i = 1 , nm1
            t = abs(e(i)) + abs(e(i+1))
            t1 = d(i) + t
            if ( xa<t1 ) xa = t1
            t1 = d(i) - t
            if ( t1<xb ) xb = t1
         enddo
         do i = 1 , icm
            g(i) = xa
            h(i) = xb
         enddo
         do k = 1 , icm
            do k1 = k , icm
               if ( g(k1)<g(k) ) then
                  g(k) = g(k1)
                  exit
               endif
            enddo
            if ( k/=1 .and. h(k)<h(k-1) ) h(k) = h(k-1)
 40         x1 = (g(k)+h(k))/2.0d0
            Cv(k) = x1
            if ( abs((g(k)-h(k))/x1)<eps ) then
               Cv(k) = x1
            else
               j = 0
               s = 1.0d0
               do i = 1 , nm
                  if ( s==0.0d0 ) s = s + 1.0d-30
                  t = f(i)/s
                  s = d(i) - t - x1
                  if ( s<0.0 ) j = j + 1
               enddo
               if ( j<k ) then
                  h(k) = x1
               else
                  g(k) = x1
                  if ( j>=icm ) then
                     g(icm) = x1
                  else
                     if ( h(j+1)<x1 ) h(j+1) = x1
                     if ( x1<g(j) ) g(j) = x1
                  endif
               endif
               goto 40
            endif
         enddo
      elseif ( Kd==1 ) then
         do ic = 1 , icm
            Cv(ic) = 4.0d0*(ic-1.0d0)**2
         enddo
      elseif ( Kd/=4 ) then
         do ic = 1 , icm
            Cv(ic) = (2.0d0*ic-1.0d0)**2
         enddo
      else
         do ic = 1 , icm
            Cv(ic) = 4.0d0*ic*ic
         enddo
      endif
      end

!       **********************************

      subroutine ittikb(x,Tti,Ttk)
!
!       =========================================================
!       Purpose: Integrate [I0(t)-1]/t with respect to t from 0
!                to x, and K0(t)/t with respect to t from x to ∞
!       Input :  x   --- Variable in the limits  ( x ≥ 0 )
!       Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
!                TTK --- Integration of K0(t)/t from x to ∞
!       =========================================================
!
      implicit none
      real(wp) e0 , el , pi , t , t1 , Tti , Ttk , x , x1
      pi = 3.141592653589793d0
      el = .5772156649015329d0
      if ( x==0.0d0 ) then
         Tti = 0.0d0
      elseif ( x<=5.0d0 ) then
         x1 = x/5.0d0
         t = x1*x1
         Tti = (((((((.1263d-3*t+.96442d-3)*t+.968217d-2)*t+.06615507d0)&
             & *t+.33116853d0)*t+1.13027241d0)*t+2.44140746d0)          &
             & *t+3.12499991d0)*t
      else
         t = 5.0d0/x
         Tti = (((((((((2.1945464d0*t-3.5195009d0)*t-11.9094395d0)*t+   &
             & 40.394734d0)*t-48.0524115d0)*t+28.1221478d0)             &
             & *t-8.6556013d0)*t+1.4780044d0)*t-.0493843d0)             &
             & *t+.1332055d0)*t + .3989314d0
         Tti = Tti*exp(x)/(sqrt(x)*x)
      endif
      if ( x==0.0d0 ) then
         Ttk = 1.0d+300
      elseif ( x<=2.0d0 ) then
         t1 = x/2.0d0
         t = t1*t1
         Ttk = (((((.77d-6*t+.1544d-4)*t+.48077d-3)*t+.925821d-2)       &
             & *t+.10937537d0)*t+.74999993d0)*t
         e0 = el + log(x/2.0d0)
         Ttk = pi*pi/24.0d0 + e0*(.5d0*e0+Tti) - Ttk
      elseif ( x<=4.0d0 ) then
         t = 2.0d0/x
         Ttk = (((.06084d0*t-.280367d0)*t+.590944d0)*t-.850013d0)       &
             & *t + 1.234684d0
         Ttk = Ttk*exp(-x)/(sqrt(x)*x)
      else
         t = 4.0d0/x
         Ttk = (((((.02724d0*t-.1110396d0)*t+.2060126d0)*t-.2621446d0)  &
             & *t+.3219184d0)*t-.5091339d0)*t + 1.2533141d0
         Ttk = Ttk*exp(-x)/(sqrt(x)*x)
      endif
      end

!       **********************************

      subroutine lqnb(n,x,Qn,Qd)
!
!       ====================================================
!       Purpose: Compute Legendre functions Qn(x) & Qn'(x)
!       Input :  x  --- Argument of Qn(x)
!                n  --- Degree of Qn(x)  ( n = 0,1,2,…)
!       Output:  QN(n) --- Qn(x)
!                QD(n) --- Qn'(x)
!       ====================================================
!
      implicit none
      real(wp) eps , q0 , q1 , qc1 , qc2 , Qd , qf , qf0 , qf1 ,&
                     & qf2 , Qn , qr , x , x2
      integer j , k , l , n , nl
      dimension Qn(0:n) , Qd(0:n)
      eps = 1.0d-14
      if ( abs(x)==1.0d0 ) then
         do k = 0 , n
            Qn(k) = 1.0d+300
            Qd(k) = 1.0d+300
         enddo
         return
      endif
      if ( x<=1.021d0 ) then
         x2 = abs((1.0d0+x)/(1.0d0-x))
         q0 = 0.5d0*log(x2)
         q1 = x*q0 - 1.0d0
         Qn(0) = q0
         Qn(1) = q1
         Qd(0) = 1.0d0/(1.0d0-x*x)
         Qd(1) = Qn(0) + x*Qd(0)
         do k = 2 , n
            qf = ((2.0d0*k-1.0d0)*x*q1-(k-1.0d0)*q0)/k
            Qn(k) = qf
            Qd(k) = (Qn(k-1)-x*qf)*k/(1.0d0-x*x)
            q0 = q1
            q1 = qf
         enddo
      else
         qc1 = 0.0d0
         qc2 = 1.0d0/x
         do j = 1 , n
            qc2 = qc2*j/((2.0*j+1.0d0)*x)
            if ( j==n-1 ) qc1 = qc2
         enddo
         do l = 0 , 1
            nl = n + l
            qf = 1.0d0
            qr = 1.0d0
            do k = 1 , 500
               qr = qr*(0.5d0*nl+k-1.0d0)*(0.5d0*(nl-1)+k)              &
                  & /((nl+k-0.5d0)*k*x*x)
               qf = qf + qr
               if ( abs(qr/qf)<eps ) exit
            enddo
            if ( l==0 ) then
               Qn(n-1) = qf*qc1
            else
               Qn(n) = qf*qc2
            endif
         enddo
         qf2 = Qn(n)
         qf1 = Qn(n-1)
         do k = n , 2 , -1
            qf0 = ((2*k-1.0d0)*x*qf1-k*qf2)/(k-1.0d0)
            Qn(k-2) = qf0
            qf2 = qf1
            qf1 = qf0
         enddo
         Qd(0) = 1.0d0/(1.0d0-x*x)
         do k = 1 , n
            Qd(k) = k*(Qn(k-1)-x*Qn(k))/(1.0d0-x*x)
         enddo
      endif
      end

!       **********************************

      subroutine cjk(Km,a)
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
      implicit none
      real(wp) a , f , f0 , g , g0
      integer j , k , Km , l1 , l2 , l3 , l4
      dimension a(*)
      a(1) = 1.0d0
      f0 = 1.0d0
      g0 = 1.0d0
      do k = 0 , Km - 1
         l1 = (k+1)*(k+2)/2 + 1
         l2 = (k+1)*(k+2)/2 + k + 2
         f = (0.5d0*k+0.125d0/(k+1))*f0
         g = -(1.5d0*k+0.625d0/(3.0*(k+1.0d0)))*g0
         a(l1) = f
         a(l2) = g
         f0 = f
         g0 = g
      enddo
      do k = 1 , Km - 1
         do j = 1 , k
            l3 = k*(k+1)/2 + j + 1
            l4 = (k+1)*(k+2)/2 + j + 1
            a(l4) = (j+0.5d0*k+0.125d0/(2.0*j+k+1.0))*a(l3)             &
                  & - (j+0.5d0*k-1.0+0.625d0/(2.0*j+k+1.0))*a(l3-1)
         enddo
      enddo
      end


!       **********************************

      subroutine ittika(x,Tti,Ttk)
!
!       =========================================================
!       Purpose: Integrate [I0(t)-1]/t with respect to t from 0
!                to x, and K0(t)/t with respect to t from x to ∞
!       Input :  x   --- Variable in the limits  ( x ≥ 0 )
!       Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
!                TTK --- Integration of K0(t)/t from x to ∞
!       =========================================================
!
      implicit none
      real(wp) b1 , c , e0 , el , pi , r , r2 , rc , rs , Tti , &
                     & Ttk , x
      integer k
      dimension c(8)
      pi = 3.141592653589793d0
      el = .5772156649015329d0
      data c/1.625d0 , 4.1328125d0 , 1.45380859375d+1 ,                 &
         & 6.553353881835d+1 , 3.6066157150269d+2 , 2.3448727161884d+3 ,&
         & 1.7588273098916d+4 , 1.4950639538279d+5/
      if ( x==0.0d0 ) then
         Tti = 0.0d0
         Ttk = 1.0d+300
         return
      endif
      if ( x<40.0d0 ) then
         Tti = 1.0d0
         r = 1.0d0
         do k = 2 , 50
            r = .25d0*r*(k-1.0d0)/(k*k*k)*x*x
            Tti = Tti + r
            if ( abs(r/Tti)<1.0d-12 ) exit
         enddo
         Tti = Tti*.125d0*x*x
      else
         Tti = 1.0d0
         r = 1.0d0
         do k = 1 , 8
            r = r/x
            Tti = Tti + c(k)*r
         enddo
         rc = x*sqrt(2.0d0*pi*x)
         Tti = Tti*exp(x)/rc
      endif
      if ( x<=12.0d0 ) then
         e0 = (.5d0*log(x/2.0d0)+el)*log(x/2.0d0) + pi*pi/24.0d0 +    &
            & .5d0*el*el
         b1 = 1.5d0 - (el+log(x/2.0d0))
         rs = 1.0d0
         r = 1.0d0
         do k = 2 , 50
            r = .25d0*r*(k-1.0d0)/(k*k*k)*x*x
            rs = rs + 1.0d0/k
            r2 = r*(rs+1.0d0/(2.0d0*k)-(el+log(x/2.0d0)))
            b1 = b1 + r2
            if ( abs(r2/b1)<1.0d-12 ) exit
         enddo
         Ttk = e0 - .125d0*x*x*b1
      else
         Ttk = 1.0d0
         r = 1.0d0
         do k = 1 , 8
            r = -r/x
            Ttk = Ttk + c(k)*r
         enddo
         rc = x*sqrt(2.0d0/pi*x)
         Ttk = Ttk*exp(-x)/rc
      endif
      end

!       **********************************

      subroutine lamv(v,x,Vm,Vl,Dl)
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
      implicit none
      real(wp) a0 , bjv0 , bjv1 , bk , ck , cs , Dl , f , f0 ,  &
                     & f1 , f2 , fac , ga , pi , px , qx , r , r0 , rc ,&
                     & rp
      real(wp) rp2 , rq , sk , uk , v , v0 , vk , Vl , Vm , vv ,&
                     & x , x2 , xk
      integer i , j , k , k0 , m , n
      dimension Vl(0:*) , Dl(0:*)
      pi = 3.141592653589793d0
      rp2 = 0.63661977236758d0
      x = abs(x)
      x2 = x*x
      n = int(v)
      v0 = v - n
      Vm = v
      if ( x<=12.0d0 ) then
         do k = 0 , n
            vk = v0 + k
            bk = 1.0d0
            r = 1.0d0
            do i = 1 , 50
               r = -0.25d0*r*x2/(i*(i+vk))
               bk = bk + r
               if ( abs(r)<abs(bk)*1.0d-15 ) exit
            enddo
            Vl(k) = bk
            uk = 1.0d0
            r = 1.0d0
            do i = 1 , 50
               r = -0.25d0*r*x2/(i*(i+vk+1.0d0))
               uk = uk + r
               if ( abs(r)<abs(uk)*1.0d-15 ) exit
            enddo
            Dl(k) = -0.5d0*x/(vk+1.0d0)*uk
         enddo
         return
      endif
      k0 = 11
      if ( x>=35.0d0 ) k0 = 10
      if ( x>=50.0d0 ) k0 = 8
      bjv0 = 0.0d0
      bjv1 = 0.0d0
      do j = 0 , 1
         vv = 4.0d0*(j+v0)*(j+v0)
         px = 1.0d0
         rp = 1.0d0
         do k = 1 , k0
            rp = -0.78125d-2*rp*(vv-(4.0*k-3.0)**2.0)                   &
               & *(vv-(4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*x2)
            px = px + rp
         enddo
         qx = 1.0d0
         rq = 1.0d0
         do k = 1 , k0
            rq = -0.78125d-2*rq*(vv-(4.0*k-1.0)**2.0)                   &
               & *(vv-(4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*x2)
            qx = qx + rq
         enddo
         qx = 0.125d0*(vv-1.0d0)*qx/x
         xk = x - (0.5d0*(j+v0)+0.25d0)*pi
         a0 = sqrt(rp2/x)
         ck = cos(xk)
         sk = sin(xk)
         if ( j==0 ) bjv0 = a0*(px*ck-qx*sk)
         if ( j==1 ) bjv1 = a0*(px*ck-qx*sk)
      enddo
      if ( v0==0.0d0 ) then
         ga = 1.0d0
      else
         call gam0(v0,ga)
         ga = v0*ga
      endif
      fac = (2.0d0/x)**v0*ga
      Vl(0) = bjv0
      Dl(0) = -bjv1 + v0/x*bjv0
      Vl(1) = bjv1
      Dl(1) = bjv0 - (1.0d0+v0)/x*bjv1
      r0 = 2.0d0*(1.0d0+v0)/x
      if ( n<=1 ) then
         Vl(0) = fac*Vl(0)
         Dl(0) = fac*Dl(0) - v0/x*Vl(0)
         Vl(1) = fac*r0*Vl(1)
         Dl(1) = fac*r0*Dl(1) - (1.0d0+v0)/x*Vl(1)
         return
      endif
      if ( n>=2 .and. n<=int(0.9*x) ) then
         f0 = bjv0
         f1 = bjv1
         do k = 2 , n
            f = 2.0d0*(k+v0-1.0d0)/x*f1 - f0
            f0 = f1
            f1 = f
            Vl(k) = f
         enddo
      elseif ( n>=2 ) then
         m = msta1(x,200)
         if ( m<n ) then
            n = m
         else
            m = msta2(x,n,15)
         endif
         f = 0.0d0
         f2 = 0.0d0
         f1 = 1.0d-100
         do k = m , 0 , -1
            f = 2.0d0*(v0+k+1.0d0)/x*f1 - f2
            if ( k<=n ) Vl(k) = f
            f2 = f1
            f1 = f
         enddo
         cs = 0.0d0
         if ( abs(bjv0)>abs(bjv1) ) then
            cs = bjv0/f
         else
            cs = bjv1/f2
         endif
         do k = 0 , n
            Vl(k) = cs*Vl(k)
         enddo
      endif
      Vl(0) = fac*Vl(0)
      do j = 1 , n
         rc = fac*r0
         Vl(j) = rc*Vl(j)
         Dl(j-1) = -0.5d0*x/(j+v0)*Vl(j)
         r0 = 2.0d0*(j+v0+1)/x*r0
      enddo
      Dl(n) = 2.0d0*(v0+n)*(Vl(n-1)-Vl(n))/x
      Vm = n + v0
      end



!       **********************************

      subroutine chguit(a,b,x,Hu,Id)
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
      implicit none
      real(wp) a , a1 , b , b1 , c , d , f1 , f2 , g , ga , Hu ,&
                     & hu0 , hu1 , hu2 , s , t , t1 , t2 , t3 , t4
      real(wp) w , x
      integer Id , j , k , m
      dimension t(30) , w(30)
      data t/.259597723012478d-01 , .778093339495366d-01 ,              &
         & .129449135396945d+00 , .180739964873425d+00 ,                &
         & .231543551376029d+00 , .281722937423262d+00 ,                &
         & .331142848268448d+00 , .379670056576798d+00 ,                &
         & .427173741583078d+00 , .473525841761707d+00 ,                &
         & .518601400058570d+00 , .562278900753945d+00 ,                &
         & .604440597048510d+00 , .644972828489477d+00 ,                &
         & .683766327381356d+00 , .720716513355730d+00 ,                &
         & .755723775306586d+00 , .788693739932264d+00 ,                &
         & .819537526162146d+00 , .848171984785930d+00 ,                &
         & .874519922646898d+00 , .898510310810046d+00 ,                &
         & .920078476177628d+00 , .939166276116423d+00 ,                &
         & .955722255839996d+00 , .969701788765053d+00 ,                &
         & .981067201752598d+00 , .989787895222222d+00 ,                &
         & .995840525118838d+00 , .999210123227436d+00/
      data w/.519078776312206d-01 , .517679431749102d-01 ,              &
         & .514884515009810d-01 , .510701560698557d-01 ,                &
         & .505141845325094d-01 , .498220356905502d-01 ,                &
         & .489955754557568d-01 , .480370318199712d-01 ,                &
         & .469489888489122d-01 , .457343797161145d-01 ,                &
         & .443964787957872d-01 , .429388928359356d-01 ,                &
         & .413655512355848d-01 , .396806954523808d-01 ,                &
         & .378888675692434d-01 , .359948980510845d-01 ,                &
         & .340038927249464d-01 , .319212190192963d-01 ,                &
         & .297524915007890d-01 , .275035567499248d-01 ,                &
         & .251804776215213d-01 , .227895169439978d-01 ,                &
         & .203371207294572d-01 , .178299010142074d-01 ,                &
         & .152746185967848d-01 , .126781664768159d-01 ,                &
         & .100475571822880d-01 , .738993116334531d-02 ,                &
         & .471272992695363d-02 , .202681196887362d-02/
      Id = 9
!       DLMF 13.4.4, integration up to C=12/X
      a1 = a - 1.0d0
      b1 = b - a - 1.0d0
      c = 12.0d0/x
      hu0 = 0.0d0
      do m = 10 , 100 , 5
         hu1 = 0.0d0
         g = 0.5d0*c/m
         d = g
         do j = 1 , m
            s = 0.0d0
            do k = 1 , 30
               t1 = d + g*t(k)
               t2 = d - g*t(k)
               f1 = exp(-x*t1)*t1**a1*(1.0d0+t1)**b1
               f2 = exp(-x*t2)*t2**a1*(1.0d0+t2)**b1
               s = s + w(k)*(f1+f2)
            enddo
            hu1 = hu1 + s*g
            d = d + 2.0d0*g
         enddo
         if ( abs(1.0d0-hu0/hu1)<1.0d-9 ) exit
         hu0 = hu1
      enddo
      call gamma2(a,ga)
      hu1 = hu1/ga
!       DLMF 13.4.4 with substitution t=C/(1-u)
!       integration u from 0 to 1, i.e. t from C=12/X to infinity
      do m = 2 , 10 , 2
         hu2 = 0.0d0
         g = 0.5d0/m
         d = g
         do j = 1 , m
            s = 0.0d0
            do k = 1 , 30
               t1 = d + g*t(k)
               t2 = d - g*t(k)
               t3 = c/(1.0d0-t1)
               t4 = c/(1.0d0-t2)
               f1 = t3*t3/c*exp(-x*t3)*t3**a1*(1.0d0+t3)**b1
               f2 = t4*t4/c*exp(-x*t4)*t4**a1*(1.0d0+t4)**b1
               s = s + w(k)*(f1+f2)
            enddo
            hu2 = hu2 + s*g
            d = d + 2.0d0*g
         enddo
         if ( abs(1.0d0-hu0/hu2)<1.0d-9 ) exit
         hu0 = hu2
      enddo
      call gamma2(a,ga)
      hu2 = hu2/ga
      Hu = hu1 + hu2
      end



!       **********************************

      subroutine kmn(m,n,c,Cv,Kd,Df,Dn,Ck1,Ck2)
!
!       ===================================================
!       Purpose: Compute the expansion coefficients of the
!                prolate and oblate spheroidal functions
!                and joining factors
!       ===================================================
!
      implicit none
      real(wp) c , Ck1 , Ck2 , cs , Cv , Df , Dn , dnp , g0 ,   &
                     & gk0 , gk1 , gk2 , gk3 , r , r1 , r2 , r3 , r4 ,  &
                     & r5 , rk
      real(wp) sa0 , sb0 , su0 , sw , t , tp , u , v , w
      integer i , ip , j , k , Kd , l , m , n , nm , nm1 , nn
      dimension u(200) , v(200) , w(200) , Df(200) , Dn(200) , tp(200) ,&
              & rk(200)
      nm = 25 + int(0.5*(n-m)+c)
      nn = nm + m
      cs = c*c*Kd
      ip = 1
      if ( n-m==2*int((n-m)/2) ) ip = 0
      k = 0
      do i = 1 , nn + 3
         if ( ip==0 ) k = -2*(i-1)
         if ( ip==1 ) k = -(2*i-3)
         gk0 = 2.0d0*m + k
         gk1 = (m+k)*(m+k+1.0d0)
         gk2 = 2.0d0*(m+k) - 1.0d0
         gk3 = 2.0d0*(m+k) + 3.0d0
         u(i) = gk0*(gk0-1.0d0)*cs/(gk2*(gk2+2.0d0))
         v(i) = gk1 - Cv + (2.0d0*(gk1-m*m)-1.0d0)*cs/(gk2*gk3)
         w(i) = (k+1.0d0)*(k+2.0d0)*cs/((gk2+2.0d0)*gk3)
      enddo
      do k = 1 , m
         t = v(m+1)
         do l = 0 , m - k - 1
            t = v(m-l) - w(m-l+1)*u(m-l)/t
         enddo
         rk(k) = -u(k)/t
      enddo
      r = 1.0d0
      do k = 1 , m
         r = r*rk(k)
         Dn(k) = Df(1)*r
      enddo
      tp(nn) = v(nn+1)
      do k = nn - 1 , m + 1 , -1
         tp(k) = v(k+1) - w(k+2)*u(k+1)/tp(k+1)
         if ( k>m+1 ) rk(k) = -u(k)/tp(k)
      enddo
      if ( m==0 ) dnp = Df(1)
      if ( m/=0 ) dnp = Dn(m)
      Dn(m+1) = (-1)**ip*dnp*cs/((2.0*m-1.0)*(2.0*m+1.0-4.0*ip)*tp(m+1))
      do k = m + 2 , nn
         Dn(k) = rk(k)*Dn(k-1)
      enddo
      r1 = 1.0d0
      do j = 1 , (n+m+ip)/2
         r1 = r1*(j+0.5d0*(n+m+ip))
      enddo
      nm1 = (n-m)/2
      r = 1.0d0
      do j = 1 , 2*m + ip
         r = r*j
      enddo
      su0 = r*Df(1)
      sw = 0.0d0
      do k = 2 , nm
         r = r*(m+k-1.0)*(m+k+ip-1.5d0)/(k-1.0d0)/(k+ip-1.5d0)
         su0 = su0 + r*Df(k)
         if ( k>nm1 .and. abs((su0-sw)/su0)<1.0d-14 ) exit
         sw = su0
      enddo
      if ( Kd/=1 ) then
         r2 = 1.0d0
         do j = 1 , m
            r2 = 2.0d0*c*r2*j
         enddo
         r3 = 1.0d0
         do j = 1 , (n-m-ip)/2
            r3 = r3*j
         enddo
         sa0 = (2.0*(m+ip)+1.0)*r1/(2.0**n*c**ip*r2*r3*Df(1))
         Ck1 = sa0*su0
         if ( Kd==-1 ) return
      endif
      r4 = 1.0d0
      do j = 1 , (n-m-ip)/2
         r4 = 4.0d0*r4*j
      enddo
      r5 = 1.0d0
      do j = 1 , m
         r5 = r5*(j+m)/c
      enddo
      g0 = Dn(m)
      if ( m==0 ) g0 = Df(1)
      sb0 = (ip+1.0)*c**(ip+1)/(2.0*ip*(m-2.0)+1.0)/(2.0*m-1.0)
      Ck2 = (-1)**ip*sb0*r4*r5*g0/r1*su0
      end



!       **********************************

      subroutine lagzo(n,x,w)
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
      implicit none
      real(wp) f0 , f1 , fd , gd , hn , p , pd , pf , q , w ,   &
                     & wp , x , z , z0
      integer i , it , j , k , n , nr
      dimension x(n) , w(n)
      hn = 1.0d0/n
      pf = 0.0d0
      pd = 0.0d0
      do nr = 1 , n
         z = hn
         if ( nr>1 ) z = x(nr-1) + hn*nr**1.27
         it = 0
 50      it = it + 1
         z0 = z
         p = 1.0d0
         do i = 1 , nr - 1
            p = p*(z-x(i))
         enddo
         f0 = 1.0d0
         f1 = 1.0d0 - z
         do k = 2 , n
            pf = ((2.0d0*k-1.0d0-z)*f1-(k-1.0d0)*f0)/k
            pd = k/z*(pf-f1)
            f0 = f1
            f1 = pf
         enddo
         fd = pf/p
         q = 0.0d0
         do i = 1 , nr - 1
            wp = 1.0d0
            do j = 1 , nr - 1
               if ( j/=i ) wp = wp*(z-x(j))
            enddo
            q = q + wp
         enddo
         gd = (pd-q*fd)/p
         z = z - fd/gd
         if ( it<=40 .and. abs((z-z0)/z)>1.0d-15 ) goto 50
         x(nr) = z
         w(nr) = 1.0d0/(z*pd*pd)
      enddo
      end

!       **********************************

      subroutine vvla(Va,x,Pv)
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
      implicit none
      real(wp) a0 , dsl , eps , gl , pdl , pi , Pv , qe , r ,   &
                     & Va , x , x1
      integer k
      pi = 3.141592653589793d0
      eps = 1.0d-12
      qe = exp(0.25*x*x)
      a0 = abs(x)**(-Va-1.0d0)*sqrt(2.0d0/pi)*qe
      r = 1.0d0
      Pv = 1.0d0
      do k = 1 , 18
         r = 0.5d0*r*(2.0*k+Va-1.0)*(2.0*k+Va)/(k*x*x)
         Pv = Pv + r
         if ( abs(r/Pv)<eps ) exit
      enddo
      Pv = a0*Pv
      if ( x<0.0d0 ) then
         x1 = -x
         call dvla(Va,x1,pdl)
         call gamma2(-Va,gl)
         dsl = sin(pi*Va)*sin(pi*Va)
         Pv = dsl*gl/pi*pdl - cos(pi*Va)*Pv
      endif
      end



!       **********************************

      subroutine cjyva(v,z,Vm,Cbj,Cdj,Cby,Cdy)
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
      implicit none
      real(wp) a0 , ga , gb , pi , pv0 , pv1 , rp2 , v , v0 ,   &
                     & vg , vl , Vm , vv , w0 , w1 , wa , ya0 , ya1 ,   &
                     & yak
      complex(wp) ca , ca0 , cb , Cbj , Cby , cck , Cdj , Cdy , cec ,    &
               & cf , cf0 , cf1 , cf2 , cfac0 , cfac1 , cg0 , cg1 ,     &
               & ch0 , ch1 , ch2
      complex(wp) ci , cju0 , cju1 , cjv0 , cjv1 , cjvl , cp11 , cp12 ,  &
               & cp21 , cp22 , cpz , cqz , cr , cr0 , cr1 , crp , crq , &
               & cs , cs0 , cs1
      complex(wp) csk , cyk , cyl1 , cyl2 , cylk , cyv0 , cyv1 , z , z1 ,&
               & z2 , zk
      integer j , k , k0 , l , lb , lb0 , m , n
      dimension Cbj(0:*) , Cdj(0:*) , Cby(0:*) , Cdy(0:*)
      pi = 3.141592653589793d0
      rp2 = .63661977236758d0
      ci = (0.0d0,1.0d0)
      a0 = abs(z)
      z1 = z
      z2 = z*z
      n = int(v)
      v0 = v - n
      pv0 = pi*v0
      pv1 = pi*(1.0d0+v0)
      if ( a0<1.0d-100 ) then
         do k = 0 , n
            Cbj(k) = (0.0d0,0.0d0)
            Cdj(k) = (0.0d0,0.0d0)
            Cby(k) = -(1.0d+300,0.0d0)
            Cdy(k) = (1.0d+300,0.0d0)
         enddo
         if ( v0==0.0 ) then
            Cbj(0) = (1.0d0,0.0d0)
            Cdj(1) = (0.5d0,0.0d0)
         else
            Cdj(0) = (1.0d+300,0.0d0)
         endif
         Vm = v
         return
      endif
      lb0 = 0.0d0
      if ( dble(z)<0.0 ) z1 = -z
      if ( a0<=12.0 ) then
         do l = 0 , 1
            vl = v0 + l
            cjvl = (1.0d0,0.0d0)
            cr = (1.0d0,0.0d0)
            do k = 1 , 40
               cr = -0.25d0*cr*z2/(k*(k+vl))
               cjvl = cjvl + cr
               if ( abs(cr)<abs(cjvl)*1.0d-15 ) exit
            enddo
            vg = 1.0d0 + vl
            call gamma2(vg,ga)
            ca = (0.5d0*z1)**vl/ga
            if ( l==0 ) cjv0 = cjvl*ca
            if ( l==1 ) cjv1 = cjvl*ca
         enddo
      else
         k0 = 11
         if ( a0>=35.0 ) k0 = 10
         if ( a0>=50.0 ) k0 = 8
         do j = 0 , 1
            vv = 4.0d0*(j+v0)*(j+v0)
            cpz = (1.0d0,0.0d0)
            crp = (1.0d0,0.0d0)
            do k = 1 , k0
               crp = -0.78125d-2*crp*(vv-(4.0*k-3.0)**2.0)              &
                   & *(vv-(4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*z2)
               cpz = cpz + crp
            enddo
            cqz = (1.0d0,0.0d0)
            crq = (1.0d0,0.0d0)
            do k = 1 , k0
               crq = -0.78125d-2*crq*(vv-(4.0*k-1.0)**2.0)              &
                   & *(vv-(4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*z2)
               cqz = cqz + crq
            enddo
            cqz = 0.125d0*(vv-1.0)*cqz/z1
            zk = z1 - (0.5d0*(j+v0)+0.25d0)*pi
            ca0 = sqrt(rp2/z1)
            cck = cos(zk)
            csk = sin(zk)
            if ( j==0 ) then
               cjv0 = ca0*(cpz*cck-cqz*csk)
               cyv0 = ca0*(cpz*csk+cqz*cck)
            elseif ( j==1 ) then
               cjv1 = ca0*(cpz*cck-cqz*csk)
               cyv1 = ca0*(cpz*csk+cqz*cck)
            endif
         enddo
      endif
      if ( a0<=12.0 ) then
         if ( v0/=0.0 ) then
            do l = 0 , 1
               vl = v0 + l
               cjvl = (1.0d0,0.0d0)
               cr = (1.0d0,0.0d0)
               do k = 1 , 40
                  cr = -0.25d0*cr*z2/(k*(k-vl))
                  cjvl = cjvl + cr
                  if ( abs(cr)<abs(cjvl)*1.0d-15 ) exit
               enddo
               vg = 1.0d0 - vl
               call gamma2(vg,gb)
               cb = (2.0d0/z1)**vl/gb
               if ( l==0 ) cju0 = cjvl*cb
               if ( l==1 ) cju1 = cjvl*cb
            enddo
            cyv0 = (cjv0*cos(pv0)-cju0)/sin(pv0)
            cyv1 = (cjv1*cos(pv1)-cju1)/sin(pv1)
         else
            cec = log(z1/2.0d0) + .5772156649015329d0
            cs0 = (0.0d0,0.0d0)
            w0 = 0.0d0
            cr0 = (1.0d0,0.0d0)
            do k = 1 , 30
               w0 = w0 + 1.0d0/k
               cr0 = -0.25d0*cr0/(k*k)*z2
               cs0 = cs0 + cr0*w0
            enddo
            cyv0 = rp2*(cec*cjv0-cs0)
            cs1 = (1.0d0,0.0d0)
            w1 = 0.0d0
            cr1 = (1.0d0,0.0d0)
            do k = 1 , 30
               w1 = w1 + 1.0d0/k
               cr1 = -0.25d0*cr1/(k*(k+1))*z2
               cs1 = cs1 + cr1*(2.0d0*w1+1.0d0/(k+1.0d0))
            enddo
            cyv1 = rp2*(cec*cjv1-1.0d0/z1-0.25d0*z1*cs1)
         endif
      endif
      if ( dble(z)<0.0d0 ) then
         cfac0 = exp(pv0*ci)
         cfac1 = exp(pv1*ci)
         if ( dimag(z)<0.0d0 ) then
            cyv0 = cfac0*cyv0 - 2.0d0*ci*cos(pv0)*cjv0
            cyv1 = cfac1*cyv1 - 2.0d0*ci*cos(pv1)*cjv1
            cjv0 = cjv0/cfac0
            cjv1 = cjv1/cfac1
         elseif ( dimag(z)>0.0d0 ) then
            cyv0 = cyv0/cfac0 + 2.0d0*ci*cos(pv0)*cjv0
            cyv1 = cyv1/cfac1 + 2.0d0*ci*cos(pv1)*cjv1
            cjv0 = cfac0*cjv0
            cjv1 = cfac1*cjv1
         endif
      endif
      Cbj(0) = cjv0
      Cbj(1) = cjv1
      if ( n>=2 .and. n<=int(0.25*a0) ) then
         cf0 = cjv0
         cf1 = cjv1
         do k = 2 , n
            cf = 2.0d0*(k+v0-1.0d0)/z*cf1 - cf0
            Cbj(k) = cf
            cf0 = cf1
            cf1 = cf
         enddo
      elseif ( n>=2 ) then
         m = msta1(a0,200)
         if ( m<n ) then
            n = m
         else
            m = msta2(a0,n,15)
         endif
         cf2 = (0.0d0,0.0d0)
         cf1 = (1.0d-100,0.0d0)
         do k = m , 0 , -1
            cf = 2.0d0*(v0+k+1.0d0)/z*cf1 - cf2
            if ( k<=n ) Cbj(k) = cf
            cf2 = cf1
            cf1 = cf
         enddo
         if ( abs(cjv0)>abs(cjv1) ) cs = cjv0/cf
         if ( abs(cjv0)<=abs(cjv1) ) cs = cjv1/cf2
         do k = 0 , n
            Cbj(k) = cs*Cbj(k)
         enddo
      endif
      Cdj(0) = v0/z*Cbj(0) - Cbj(1)
      do k = 1 , n
         Cdj(k) = -(k+v0)/z*Cbj(k) + Cbj(k-1)
      enddo
      Cby(0) = cyv0
      Cby(1) = cyv1
      ya0 = abs(cyv0)
      lb = 0
      cg0 = cyv0
      cg1 = cyv1
      do k = 2 , n
         cyk = 2.0d0*(v0+k-1.0d0)/z*cg1 - cg0
         if ( abs(cyk)<=1.0d+290 ) then
            yak = abs(cyk)
            ya1 = abs(cg0)
            if ( yak<ya0 .and. yak<ya1 ) lb = k
            Cby(k) = cyk
            cg0 = cg1
            cg1 = cyk
         endif
      enddo
      if ( lb>4 .and. dimag(z)/=0.0d0 ) then
 50      if ( lb/=lb0 ) then
            ch2 = (1.0d0,0.0d0)
            ch1 = (0.0d0,0.0d0)
            lb0 = lb
            do k = lb , 1 , -1
               ch0 = 2.0d0*(k+v0)/z*ch1 - ch2
               ch2 = ch1
               ch1 = ch0
            enddo
            cp12 = ch0
            cp22 = ch2
            ch2 = (0.0d0,0.0d0)
            ch1 = (1.0d0,0.0d0)
            do k = lb , 1 , -1
               ch0 = 2.0d0*(k+v0)/z*ch1 - ch2
               ch2 = ch1
               ch1 = ch0
            enddo
            cp11 = ch0
            cp21 = ch2
            if ( lb==n ) Cbj(lb+1) = 2.0d0*(lb+v0)/z*Cbj(lb) - Cbj(lb-1)
            if ( abs(Cbj(0))>abs(Cbj(1)) ) then
               Cby(lb+1) = (Cbj(lb+1)*cyv0-2.0d0*cp11/(pi*z))/Cbj(0)
               Cby(lb) = (Cbj(lb)*cyv0+2.0d0*cp12/(pi*z))/Cbj(0)
            else
               Cby(lb+1) = (Cbj(lb+1)*cyv1-2.0d0*cp21/(pi*z))/Cbj(1)
               Cby(lb) = (Cbj(lb)*cyv1+2.0d0*cp22/(pi*z))/Cbj(1)
            endif
            cyl2 = Cby(lb+1)
            cyl1 = Cby(lb)
            do k = lb - 1 , 0 , -1
               cylk = 2.0d0*(k+v0+1.0d0)/z*cyl1 - cyl2
               Cby(k) = cylk
               cyl2 = cyl1
               cyl1 = cylk
            enddo
            cyl1 = Cby(lb)
            cyl2 = Cby(lb+1)
            do k = lb + 1 , n - 1
               cylk = 2.0d0*(k+v0)/z*cyl2 - cyl1
               Cby(k+1) = cylk
               cyl1 = cyl2
               cyl2 = cylk
            enddo
            do k = 2 , n
               wa = abs(Cby(k))
               if ( wa<abs(Cby(k-1)) ) lb = k
            enddo
            goto 50
         endif
      endif
      Cdy(0) = v0/z*Cby(0) - Cby(1)
      do k = 1 , n
         Cdy(k) = Cby(k-1) - (k+v0)/z*Cby(k)
      enddo
      Vm = n + v0
      end



!       **********************************

      subroutine cjyvb(v,z,Vm,Cbj,Cdj,Cby,Cdy)
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
      implicit none
      real(wp) a0 , ga , gb , pi , pv0 , rp2 , v , v0 , vg ,    &
                     & Vm , vv , w0
      complex(wp) ca , ca0 , cb , Cbj , Cby , cck , Cdj , Cdy , cec ,    &
               & cf , cf1 , cf2 , cfac0 , ci , cju0 , cjv0 , cjvn ,     &
               & cpz , cqz , cr
      complex(wp) cr0 , crp , crq , cs , cs0 , csk , cyv0 , cyy , z ,    &
               & z1 , z2 , zk
      integer k , k0 , m , n
      dimension Cbj(0:*) , Cdj(0:*) , Cby(0:*) , Cdy(0:*)
      pi = 3.141592653589793d0
      rp2 = .63661977236758d0
      ci = (0.0d0,1.0d0)
      a0 = abs(z)
      z1 = z
      z2 = z*z
      n = int(v)
      v0 = v - n
      pv0 = pi*v0
      if ( a0<1.0d-100 ) then
         do k = 0 , n
            Cbj(k) = (0.0d0,0.0d0)
            Cdj(k) = (0.0d0,0.0d0)
            Cby(k) = -(1.0d+300,0.0d0)
            Cdy(k) = (1.0d+300,0.0d0)
         enddo
         if ( v0==0.0 ) then
            Cbj(0) = (1.0d0,0.0d0)
            Cdj(1) = (0.5d0,0.0d0)
         else
            Cdj(0) = (1.0d+300,0.0d0)
         endif
         Vm = v
         return
      endif
      if ( dble(z)<0.0d0 ) z1 = -z
      if ( a0<=12.0 ) then
         cjv0 = (1.0d0,0.0d0)
         cr = (1.0d0,0.0d0)
         do k = 1 , 40
            cr = -0.25d0*cr*z2/(k*(k+v0))
            cjv0 = cjv0 + cr
            if ( abs(cr)<abs(cjv0)*1.0d-15 ) exit
         enddo
         vg = 1.0d0 + v0
         call gamma2(vg,ga)
         ca = (0.5d0*z1)**v0/ga
         cjv0 = cjv0*ca
      else
         k0 = 11
         if ( a0>=35.0 ) k0 = 10
         if ( a0>=50.0 ) k0 = 8
         vv = 4.0d0*v0*v0
         cpz = (1.0d0,0.0d0)
         crp = (1.0d0,0.0d0)
         do k = 1 , k0
            crp = -0.78125d-2*crp*(vv-(4.0*k-3.0)**2.0)                 &
                & *(vv-(4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*z2)
            cpz = cpz + crp
         enddo
         cqz = (1.0d0,0.0d0)
         crq = (1.0d0,0.0d0)
         do k = 1 , k0
            crq = -0.78125d-2*crq*(vv-(4.0*k-1.0)**2.0)                 &
                & *(vv-(4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*z2)
            cqz = cqz + crq
         enddo
         cqz = 0.125d0*(vv-1.0)*cqz/z1
         zk = z1 - (0.5d0*v0+0.25d0)*pi
         ca0 = sqrt(rp2/z1)
         cck = cos(zk)
         csk = sin(zk)
         cjv0 = ca0*(cpz*cck-cqz*csk)
         cyv0 = ca0*(cpz*csk+cqz*cck)
      endif
      if ( a0<=12.0 ) then
         if ( v0/=0.0 ) then
            cjvn = (1.0d0,0.0d0)
            cr = (1.0d0,0.0d0)
            do k = 1 , 40
               cr = -0.25d0*cr*z2/(k*(k-v0))
               cjvn = cjvn + cr
               if ( abs(cr)<abs(cjvn)*1.0d-15 ) exit
            enddo
            vg = 1.0d0 - v0
            call gamma2(vg,gb)
            cb = (2.0d0/z1)**v0/gb
            cju0 = cjvn*cb
            cyv0 = (cjv0*cos(pv0)-cju0)/sin(pv0)
         else
            cec = log(z1/2.0d0) + .5772156649015329d0
            cs0 = (0.0d0,0.0d0)
            w0 = 0.0d0
            cr0 = (1.0d0,0.0d0)
            do k = 1 , 30
               w0 = w0 + 1.0d0/k
               cr0 = -0.25d0*cr0/(k*k)*z2
               cs0 = cs0 + cr0*w0
            enddo
            cyv0 = rp2*(cec*cjv0-cs0)
         endif
      endif
      if ( n==0 ) n = 1
      m = msta1(a0,200)
      if ( m<n ) then
         n = m
      else
         m = msta2(a0,n,15)
      endif
      cf2 = (0.0d0,0.0d0)
      cf1 = (1.0d-100,0.0d0)
      do k = m , 0 , -1
         cf = 2.0d0*(v0+k+1.0d0)/z1*cf1 - cf2
         if ( k<=n ) Cbj(k) = cf
         cf2 = cf1
         cf1 = cf
      enddo
      cs = cjv0/cf
      do k = 0 , n
         Cbj(k) = cs*Cbj(k)
      enddo
      if ( dble(z)<0.0d0 ) then
         cfac0 = exp(pv0*ci)
         if ( dimag(z)<0.0d0 ) then
            cyv0 = cfac0*cyv0 - 2.0d0*ci*cos(pv0)*cjv0
         elseif ( dimag(z)>0.0d0 ) then
            cyv0 = cyv0/cfac0 + 2.0d0*ci*cos(pv0)*cjv0
         endif
         do k = 0 , n
            if ( dimag(z)<0.0d0 ) then
               Cbj(k) = exp(-pi*(k+v0)*ci)*Cbj(k)
            elseif ( dimag(z)>0.0d0 ) then
               Cbj(k) = exp(pi*(k+v0)*ci)*Cbj(k)
            endif
         enddo
         z1 = z1
      endif
      Cby(0) = cyv0
      do k = 1 , n
         cyy = (Cbj(k)*Cby(k-1)-2.0d0/(pi*z))/Cbj(k-1)
         Cby(k) = cyy
      enddo
      Cdj(0) = v0/z*Cbj(0) - Cbj(1)
      do k = 1 , n
         Cdj(k) = -(k+v0)/z*Cbj(k) + Cbj(k-1)
      enddo
      Cdy(0) = v0/z*Cby(0) - Cby(1)
      do k = 1 , n
         Cdy(k) = Cby(k-1) - (k+v0)/z*Cby(k)
      enddo
      Vm = n + v0
      end



!       **********************************

      subroutine jy01a(x,Bj0,Dj0,Bj1,Dj1,By0,Dy0,By1,Dy1)
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
      implicit none
      real(wp) a , a1 , b , b1 , Bj0 , Bj1 , By0 , By1 , cs0 ,  &
                     & cs1 , cu , Dj0 , Dj1 , Dy0 , Dy1 , ec , p0 , p1 ,&
                     & pi , q0
      real(wp) q1 , r , r0 , r1 , rp2 , t1 , t2 , w0 , w1 , x , &
                     & x2
      integer k , k0
      dimension a(12) , b(12) , a1(12) , b1(12)
      pi = 3.141592653589793d0
      rp2 = 0.63661977236758d0
      x2 = x*x
      if ( x==0.0d0 ) then
         Bj0 = 1.0d0
         Bj1 = 0.0d0
         Dj0 = 0.0d0
         Dj1 = 0.5d0
         By0 = -1.0d+300
         By1 = -1.0d+300
         Dy0 = 1.0d+300
         Dy1 = 1.0d+300
         return
      endif
      if ( x<=12.0d0 ) then
         Bj0 = 1.0d0
         r = 1.0d0
         do k = 1 , 30
            r = -0.25d0*r*x2/(k*k)
            Bj0 = Bj0 + r
            if ( abs(r)<abs(Bj0)*1.0d-15 ) exit
         enddo
         Bj1 = 1.0d0
         r = 1.0d0
         do k = 1 , 30
            r = -0.25d0*r*x2/(k*(k+1.0d0))
            Bj1 = Bj1 + r
            if ( abs(r)<abs(Bj1)*1.0d-15 ) exit
         enddo
         Bj1 = 0.5d0*x*Bj1
         ec = log(x/2.0d0) + 0.5772156649015329d0
         cs0 = 0.0d0
         w0 = 0.0d0
         r0 = 1.0d0
         do k = 1 , 30
            w0 = w0 + 1.0d0/k
            r0 = -0.25d0*r0/(k*k)*x2
            r = r0*w0
            cs0 = cs0 + r
            if ( abs(r)<abs(cs0)*1.0d-15 ) exit
         enddo
         By0 = rp2*(ec*Bj0-cs0)
         cs1 = 1.0d0
         w1 = 0.0d0
         r1 = 1.0d0
         do k = 1 , 30
            w1 = w1 + 1.0d0/k
            r1 = -0.25d0*r1/(k*(k+1))*x2
            r = r1*(2.0d0*w1+1.0d0/(k+1.0d0))
            cs1 = cs1 + r
            if ( abs(r)<abs(cs1)*1.0d-15 ) exit
         enddo
         By1 = rp2*(ec*Bj1-1.0d0/x-0.25d0*x*cs1)
      else
         data a/ - .7031250000000000d-01 , .1121520996093750d+00 ,      &
            & -.5725014209747314d+00 , .6074042001273483d+01 ,          &
            & -.1100171402692467d+03 , .3038090510922384d+04 ,          &
            & -.1188384262567832d+06 , .6252951493434797d+07 ,          &
            & -.4259392165047669d+09 , .3646840080706556d+11 ,          &
            & -.3833534661393944d+13 , .4854014686852901d+15/
         data b/.7324218750000000d-01 , -.2271080017089844d+00 ,        &
            & .1727727502584457d+01 , -.2438052969955606d+02 ,          &
            & .5513358961220206d+03 , -.1825775547429318d+05 ,          &
            & .8328593040162893d+06 , -.5006958953198893d+08 ,          &
            & .3836255180230433d+10 , -.3649010818849833d+12 ,          &
            & .4218971570284096d+14 , -.5827244631566907d+16/
         data a1/.1171875000000000d+00 , -.1441955566406250d+00 ,       &
            & .6765925884246826d+00 , -.6883914268109947d+01 ,          &
            & .1215978918765359d+03 , -.3302272294480852d+04 ,          &
            & .1276412726461746d+06 , -.6656367718817688d+07 ,          &
            & .4502786003050393d+09 , -.3833857520742790d+11 ,          &
            & .4011838599133198d+13 , -.5060568503314727d+15/
         data b1/ - .1025390625000000d+00 , .2775764465332031d+00 ,     &
            & -.1993531733751297d+01 , .2724882731126854d+02 ,          &
            & -.6038440767050702d+03 , .1971837591223663d+05 ,          &
            & -.8902978767070678d+06 , .5310411010968522d+08 ,          &
            & -.4043620325107754d+10 , .3827011346598605d+12 ,          &
            & -.4406481417852278d+14 , .6065091351222699d+16/
         k0 = 12
         if ( x>=35.0 ) k0 = 10
         if ( x>=50.0 ) k0 = 8
         t1 = x - 0.25d0*pi
         p0 = 1.0d0
         q0 = -0.125d0/x
         do k = 1 , k0
            p0 = p0 + a(k)*x**(-2*k)
            q0 = q0 + b(k)*x**(-2*k-1)
         enddo
         cu = sqrt(rp2/x)
         Bj0 = cu*(p0*cos(t1)-q0*sin(t1))
         By0 = cu*(p0*sin(t1)+q0*cos(t1))
         t2 = x - 0.75d0*pi
         p1 = 1.0d0
         q1 = 0.375d0/x
         do k = 1 , k0
            p1 = p1 + a1(k)*x**(-2*k)
            q1 = q1 + b1(k)*x**(-2*k-1)
         enddo
         cu = sqrt(rp2/x)
         Bj1 = cu*(p1*cos(t2)-q1*sin(t2))
         By1 = cu*(p1*sin(t2)+q1*cos(t2))
      endif
      Dj0 = -Bj1
      Dj1 = Bj0 - Bj1/x
      Dy0 = -By1
      Dy1 = By0 - By1/x
      end

!       **********************************

      subroutine incog(a,x,Gin,Gim,Gip,Isfer)
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
      implicit none
      real(wp) a , ga , Gim , Gin , Gip , r , s , t0 , x , xam
      integer Isfer , k
      Isfer = 0
      xam = -x + a*log(x)
      if ( xam>700.0 .or. a>170.0 ) then
         Isfer = 6
         return
      endif
      if ( x==0.0 ) then
         Gin = 0.0
         call gamma2(a,ga)
         Gim = ga
         Gip = 0.0
      elseif ( x<=1.0+a ) then
         s = 1.0d0/a
         r = s
         do k = 1 , 60
            r = r*x/(a+k)
            s = s + r
            if ( abs(r/s)<1.0d-15 ) exit
         enddo
         Gin = exp(xam)*s
         call gamma2(a,ga)
         Gip = Gin/ga
         Gim = ga - Gin
      elseif ( x>1.0+a ) then
         t0 = 0.0d0
         do k = 60 , 1 , -1
            t0 = (k-a)/(1.0d0+k/(x+t0))
         enddo
         Gim = exp(xam)/(x+t0)
         call gamma2(a,ga)
         Gin = ga - Gim
         Gip = 1.0d0 - Gim/ga
      endif
      end



!       **********************************

      subroutine itikb(x,Ti,Tk)
!
!       =======================================================
!       Purpose: Integrate Bessel functions I0(t) and K0(t)
!                with respect to t from 0 to x
!       Input :  x  --- Upper limit of the integral ( x ≥ 0 )
!       Output:  TI --- Integration of I0(t) from 0 to x
!                TK --- Integration of K0(t) from 0 to x
!       =======================================================
!
      implicit none
      real(wp) pi , t , t1 , Ti , Tk , x
      pi = 3.141592653589793d0
      if ( x==0.0d0 ) then
         Ti = 0.0d0
      elseif ( x<5.0d0 ) then
         t1 = x/5.0d0
         t = t1*t1
         Ti = ((((((((.59434d-3*t+.4500642d-2)*t+.044686921d0)*t+       &
            & .300704878d0)*t+1.471860153d0)*t+4.844024624d0)           &
            & *t+9.765629849d0)*t+10.416666367d0)*t+5.0d0)*t1
      elseif ( x>=5.0 .and. x<=8.0d0 ) then
         t = 5.0d0/x
         Ti = (((-.015166d0*t-.0202292d0)*t+.1294122d0)*t-.0302912d0)   &
            & *t + .4161224d0
         Ti = Ti*exp(x)/sqrt(x)
      else
         t = 8.0d0/x
         Ti = (((((-.0073995d0*t+.017744d0)*t-.0114858d0)*t+.55956d-2)  &
            & *t+.59191d-2)*t+.0311734d0)*t + .3989423d0
         Ti = Ti*exp(x)/sqrt(x)
      endif
      if ( x==0.0d0 ) then
         Tk = 0.0d0
      elseif ( x<=2.0d0 ) then
         t1 = x/2.0d0
         t = t1*t1
         Tk = ((((((.116d-5*t+.2069d-4)*t+.62664d-3)*t+.01110118d0)*t+  &
            & .11227902d0)*t+.50407836d0)*t+.84556868d0)*t1
         Tk = Tk - log(x/2.0d0)*Ti
      elseif ( x>2.0 .and. x<=4.0d0 ) then
         t = 2.0d0/x
         Tk = (((.0160395d0*t-.0781715d0)*t+.185984d0)*t-.3584641d0)    &
            & *t + 1.2494934d0
         Tk = pi/2.0d0 - Tk*exp(-x)/sqrt(x)
      elseif ( x>4.0 .and. x<=7.0d0 ) then
         t = 4.0d0/x
         Tk = (((((.37128d-2*t-.0158449d0)*t+.0320504d0)*t-.0481455d0)  &
            & *t+.0787284d0)*t-.1958273d0)*t + 1.2533141d0
         Tk = pi/2.0d0 - Tk*exp(-x)/sqrt(x)
      else
         t = 7.0d0/x
         Tk = (((((.33934d-3*t-.163271d-2)*t+.417454d-2)*t-.933944d-2)  &
            & *t+.02576646d0)*t-.11190289d0)*t + 1.25331414d0
         Tk = pi/2.0d0 - Tk*exp(-x)/sqrt(x)
      endif
      end

!       **********************************

      subroutine itika(x,Ti,Tk)
!
!       =======================================================
!       Purpose: Integrate modified Bessel functions I0(t) and
!                K0(t) with respect to t from 0 to x
!       Input :  x  --- Upper limit of the integral  ( x ≥ 0 )
!       Output:  TI --- Integration of I0(t) from 0 to x
!                TK --- Integration of K0(t) from 0 to x
!       =======================================================
!
      implicit none
      real(wp) a , b1 , b2 , e0 , el , pi , r , rc1 , rc2 , rs ,&
                     & Ti , Tk , tw , x , x2
      integer k
      dimension a(10)
      pi = 3.141592653589793d0
      el = .5772156649015329d0
      data a/.625d0 , 1.0078125d0 , 2.5927734375d0 , 9.1868591308594d0 ,&
         & 4.1567974090576d+1 , 2.2919635891914d+2 , 1.491504060477d+3 ,&
         & 1.1192354495579d+4 , 9.515939374212d+4 , 9.0412425769041d+5/
      if ( x==0.0d0 ) then
         Ti = 0.0d0
         Tk = 0.0d0
         return
      elseif ( x<20.0d0 ) then
         x2 = x*x
         Ti = 1.0d0
         r = 1.0d0
         do k = 1 , 50
            r = .25d0*r*(2*k-1.0d0)/(2*k+1.0d0)/(k*k)*x2
            Ti = Ti + r
            if ( abs(r/Ti)<1.0d-12 ) exit
         enddo
         Ti = Ti*x
      else
         x2 = 0.0d0
         Ti = 1.0d0
         r = 1.0d0
         do k = 1 , 10
            r = r/x
            Ti = Ti + a(k)*r
         enddo
         rc1 = 1.0d0/sqrt(2.0d0*pi*x)
         Ti = rc1*exp(x)*Ti
      endif
      if ( x<12.0d0 ) then
         e0 = el + log(x/2.0d0)
         b1 = 1.0d0 - e0
         b2 = 0.0d0
         rs = 0.0d0
         r = 1.0d0
         tw = 0.0d0
         do k = 1 , 50
            r = .25d0*r*(2*k-1.0d0)/(2*k+1.0d0)/(k*k)*x2
            b1 = b1 + r*(1.0d0/(2*k+1)-e0)
            rs = rs + 1.0d0/k
            b2 = b2 + r*rs
            Tk = b1 + b2
            if ( abs((Tk-tw)/Tk)<1.0d-12 ) exit
            tw = Tk
         enddo
         Tk = Tk*x
      else
         Tk = 1.0d0
         r = 1.0d0
         do k = 1 , 10
            r = -r/x
            Tk = Tk + a(k)*r
         enddo
         rc2 = sqrt(pi/(2.0d0*x))
         Tk = pi/2.0d0 - rc2*Tk*exp(-x)
      endif
      end

!       **********************************

      subroutine jyv(v,x,Vm,Bj,Dj,By,Dy)
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
      implicit none
      real(wp) a , a0 , b , Bj , bju0 , bju1 , bjv0 , bjv1 ,    &
                     & bjvl , By , byv0 , byv1 , byvk , ck , cs , cs0 , &
                     & cs1 , Dj , Dy , ec
      real(wp) el , f , f0 , f1 , f2 , ga , gb , pi , pv0 ,     &
                     & pv1 , px , qx , r , r0 , r1 , rp , rp2 , rq ,    &
                     & sk , v
      real(wp) v0 , vg , vl , Vm , vv , w0 , w1 , x , x2 , xk
      integer j , k , k0 , l , m , n
      dimension Bj(0:*) , Dj(0:*) , By(0:*) , Dy(0:*)
      el = .5772156649015329d0
      pi = 3.141592653589793d0
      rp2 = .63661977236758d0
      x2 = x*x
      n = int(v)
      v0 = v - n
      if ( x<1.0d-100 ) then
         do k = 0 , n
            Bj(k) = 0.0d0
            Dj(k) = 0.0d0
            By(k) = -1.0d+300
            Dy(k) = 1.0d+300
         enddo
         if ( v0==0.0 ) then
            Bj(0) = 1.0d0
            Dj(1) = 0.5d0
         else
            Dj(0) = 1.0d+300
         endif
         Vm = v
         return
      endif
      bjv0 = 0.0d0
      bjv1 = 0.0d0
      byv0 = 0.0d0
      byv1 = 0.0d0
      if ( x<=12.0 ) then
         do l = 0 , 1
            vl = v0 + l
            bjvl = 1.0d0
            r = 1.0d0
            do k = 1 , 40
               r = -0.25d0*r*x2/(k*(k+vl))
               bjvl = bjvl + r
               if ( abs(r)<abs(bjvl)*1.0d-15 ) exit
            enddo
            vg = 1.0d0 + vl
            call gamma2(vg,ga)
            a = (0.5d0*x)**vl/ga
            if ( l==0 ) bjv0 = bjvl*a
            if ( l==1 ) bjv1 = bjvl*a
         enddo
      else
         k0 = 11
         if ( x>=35.0 ) k0 = 10
         if ( x>=50.0 ) k0 = 8
         do j = 0 , 1
            vv = 4.0d0*(j+v0)*(j+v0)
            px = 1.0d0
            rp = 1.0d0
            do k = 1 , k0
               rp = -0.78125d-2*rp*(vv-(4.0*k-3.0)**2.0)                &
                  & *(vv-(4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*x2)
               px = px + rp
            enddo
            qx = 1.0d0
            rq = 1.0d0
            do k = 1 , k0
               rq = -0.78125d-2*rq*(vv-(4.0*k-1.0)**2.0)                &
                  & *(vv-(4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*x2)
               qx = qx + rq
            enddo
            qx = 0.125d0*(vv-1.0)*qx/x
            xk = x - (0.5d0*(j+v0)+0.25d0)*pi
            a0 = sqrt(rp2/x)
            ck = cos(xk)
            sk = sin(xk)
            if ( j==0 ) then
               bjv0 = a0*(px*ck-qx*sk)
               byv0 = a0*(px*sk+qx*ck)
            elseif ( j==1 ) then
               bjv1 = a0*(px*ck-qx*sk)
               byv1 = a0*(px*sk+qx*ck)
            endif
         enddo
      endif
      Bj(0) = bjv0
      Bj(1) = bjv1
      Dj(0) = v0/x*Bj(0) - Bj(1)
      Dj(1) = -(1.0d0+v0)/x*Bj(1) + Bj(0)
      if ( n>=2 .and. n<=int(0.9*x) ) then
         f0 = bjv0
         f1 = bjv1
         do k = 2 , n
            f = 2.0d0*(k+v0-1.0d0)/x*f1 - f0
            Bj(k) = f
            f0 = f1
            f1 = f
         enddo
      elseif ( n>=2 ) then
         m = msta1(x,200)
         if ( m<n ) then
            n = m
         else
            m = msta2(x,n,15)
         endif
         f = 0.0d0
         f2 = 0.0d0
         f1 = 1.0d-100
         do k = m , 0 , -1
            f = 2.0d0*(v0+k+1.0d0)/x*f1 - f2
            if ( k<=n ) Bj(k) = f
            f2 = f1
            f1 = f
         enddo
         if ( abs(bjv0)>abs(bjv1) ) then
            cs = bjv0/f
         else
            cs = bjv1/f2
         endif
         do k = 0 , n
            Bj(k) = cs*Bj(k)
         enddo
      endif
      do k = 2 , n
         Dj(k) = -(k+v0)/x*Bj(k) + Bj(k-1)
      enddo
      if ( x<=12.0d0 ) then
         if ( v0/=0.0 ) then
            bju0 = 0.0d0
            bju1 = 0.0d0
            do l = 0 , 1
               vl = v0 + l
               bjvl = 1.0d0
               r = 1.0d0
               do k = 1 , 40
                  r = -0.25d0*r*x2/(k*(k-vl))
                  bjvl = bjvl + r
                  if ( abs(r)<abs(bjvl)*1.0d-15 ) exit
               enddo
               vg = 1.0d0 - vl
               call gamma2(vg,gb)
               b = (2.0d0/x)**vl/gb
               if ( l==0 ) bju0 = bjvl*b
               if ( l==1 ) bju1 = bjvl*b
            enddo
            pv0 = pi*v0
            pv1 = pi*(1.0d0+v0)
            byv0 = (bjv0*cos(pv0)-bju0)/sin(pv0)
            byv1 = (bjv1*cos(pv1)-bju1)/sin(pv1)
         else
            ec = log(x/2.0d0) + el
            cs0 = 0.0d0
            w0 = 0.0d0
            r0 = 1.0d0
            do k = 1 , 30
               w0 = w0 + 1.0d0/k
               r0 = -0.25d0*r0/(k*k)*x2
               cs0 = cs0 + r0*w0
            enddo
            byv0 = rp2*(ec*bjv0-cs0)
            cs1 = 1.0d0
            w1 = 0.0d0
            r1 = 1.0d0
            do k = 1 , 30
               w1 = w1 + 1.0d0/k
               r1 = -0.25d0*r1/(k*(k+1))*x2
               cs1 = cs1 + r1*(2.0d0*w1+1.0d0/(k+1.0d0))
            enddo
            byv1 = rp2*(ec*bjv1-1.0d0/x-0.25d0*x*cs1)
         endif
      endif
      By(0) = byv0
      By(1) = byv1
      do k = 2 , n
         byvk = 2.0d0*(v0+k-1.0d0)/x*byv1 - byv0
         By(k) = byvk
         byv0 = byv1
         byv1 = byvk
      enddo
      Dy(0) = v0/x*By(0) - By(1)
      do k = 1 , n
         Dy(k) = -(k+v0)/x*By(k) + By(k-1)
      enddo
      Vm = n + v0
      end



!       **********************************

      subroutine jynb(n,x,Nm,Bj,Dj,By,Dy)
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
      implicit none
      real(wp) Bj , By , Dj , Dy , x
      integer k , n , Nm
      dimension Bj(0:n) , Dj(0:n) , By(0:n) , Dy(0:n)
      call jynbh(n,0,x,Nm,Bj,By)
!       Compute derivatives by differentiation formulas
      if ( x<1.0d-100 ) then
         do k = 0 , n
            Dj(k) = 0.0d0
            Dy(k) = 1.0d+300
         enddo
         Dj(1) = 0.5d0
      else
         Dj(0) = -Bj(1)
         do k = 1 , Nm
            Dj(k) = Bj(k-1) - k/x*Bj(k)
         enddo
         Dy(0) = -By(1)
         do k = 1 , Nm
            Dy(k) = By(k-1) - k*By(k)/x
         enddo
      endif
      end


!       **********************************

      subroutine jynbh(n,Nmin,x,Nm,Bj,By)
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
      implicit none
      real(wp) a , a1 , b , b1 , Bj , bj0 , bj1 , bjk , bs ,    &
                     & By , by0 , by1 , byk , cu , ec , f , f1 , f2 ,   &
                     & p0 , p1
      real(wp) pi , q0 , q1 , r2p , s0 , su , sv , t1 , t2 , x
      integer k , ky , m , n , Nm , Nmin
      dimension Bj(0:n-Nmin) , By(0:n-Nmin) , a(4) , b(4) , a1(4) ,     &
              & b1(4)
      pi = 3.141592653589793d0
      r2p = .63661977236758d0
      Nm = n
      if ( x<1.0d-100 ) then
         do k = Nmin , n
            Bj(k-Nmin) = 0.0d0
            By(k-Nmin) = -1.0d+300
         enddo
         if ( Nmin==0 ) Bj(0) = 1.0d0
         return
      endif
      if ( x<=300.0 .or. n>int(0.9*x) ) then
!          Backward recurrence for Jn
         if ( n==0 ) Nm = 1
         m = msta1(x,200)
         if ( m<Nm ) then
            Nm = m
         else
            m = msta2(x,Nm,15)
         endif
         bs = 0.0d0
         su = 0.0d0
         sv = 0.0d0
         f2 = 0.0d0
         f1 = 1.0d-100
         f = 0.0d0
         do k = m , 0 , -1
            f = 2.0d0*(k+1.0d0)/x*f1 - f2
            if ( k<=Nm .and. k>=Nmin ) Bj(k-Nmin) = f
            if ( k==2*int(k/2) .and. k/=0 ) then
               bs = bs + 2.0d0*f
               su = su + (-1)**(k/2)*f/k
            elseif ( k>1 ) then
               sv = sv + (-1)**(k/2)*k/(k*k-1.0d0)*f
            endif
            f2 = f1
            f1 = f
         enddo
         s0 = bs + f
         do k = Nmin , Nm
            Bj(k-Nmin) = Bj(k-Nmin)/s0
         enddo
!          Estimates for Yn at start of recurrence
         bj0 = f1/s0
         bj1 = f2/s0
         ec = log(x/2.0d0) + 0.5772156649015329d0
         by0 = r2p*(ec*bj0-4.0d0*su/s0)
         by1 = r2p*((ec-1.0d0)*bj1-bj0/x-4.0d0*sv/s0)
         if ( 0>=Nmin ) By(0-Nmin) = by0
         if ( 1>=Nmin ) By(1-Nmin) = by1
         ky = 2
      else
!          Hankel expansion
         data a/ - .7031250000000000d-01 , .1121520996093750d+00 ,      &
            & -.5725014209747314d+00 , .6074042001273483d+01/
         data b/.7324218750000000d-01 , -.2271080017089844d+00 ,        &
            & .1727727502584457d+01 , -.2438052969955606d+02/
         data a1/.1171875000000000d+00 , -.1441955566406250d+00 ,       &
            & .6765925884246826d+00 , -.6883914268109947d+01/
         data b1/ - .1025390625000000d+00 , .2775764465332031d+00 ,     &
            & -.1993531733751297d+01 , .2724882731126854d+02/
         t1 = x - 0.25d0*pi
         p0 = 1.0d0
         q0 = -0.125d0/x
         do k = 1 , 4
            p0 = p0 + a(k)*x**(-2*k)
            q0 = q0 + b(k)*x**(-2*k-1)
         enddo
         cu = sqrt(r2p/x)
         bj0 = cu*(p0*cos(t1)-q0*sin(t1))
         by0 = cu*(p0*sin(t1)+q0*cos(t1))
         if ( 0>=Nmin ) Bj(0-Nmin) = bj0
         if ( 0>=Nmin ) By(0-Nmin) = by0
         t2 = x - 0.75d0*pi
         p1 = 1.0d0
         q1 = 0.375d0/x
         do k = 1 , 4
            p1 = p1 + a1(k)*x**(-2*k)
            q1 = q1 + b1(k)*x**(-2*k-1)
         enddo
         bj1 = cu*(p1*cos(t2)-q1*sin(t2))
         by1 = cu*(p1*sin(t2)+q1*cos(t2))
         if ( 1>=Nmin ) Bj(1-Nmin) = bj1
         if ( 1>=Nmin ) By(1-Nmin) = by1
         do k = 2 , Nm
            bjk = 2.0d0*(k-1.0d0)/x*bj1 - bj0
            if ( k>=Nmin ) Bj(k-Nmin) = bjk
            bj0 = bj1
            bj1 = bjk
         enddo
         ky = 2
      endif
!       Forward recurrence for Yn
      do k = ky , Nm
         byk = 2.0d0*(k-1.0d0)*by1/x - by0
         if ( k>=Nmin ) By(k-Nmin) = byk
         by0 = by1
         by1 = byk
      enddo
      end

!       **********************************

      subroutine legzo(n,x,w)
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
      implicit none
      real(wp) f0 , f1 , fd , gd , p , pd , pf , q , w , wp ,   &
                     & x , z , z0
      integer i , j , k , n , n0 , nr
      dimension x(n) , w(n)
      n0 = (n+1)/2
      pf = 0.0d0
      pd = 0.0d0
      do nr = 1 , n0
         z = cos(3.1415926d0*(nr-0.25d0)/n)
 50      z0 = z
         p = 1.0d0
         do i = 1 , nr - 1
            p = p*(z-x(i))
         enddo
         f0 = 1.0d0
         if ( nr==n0 .and. n/=2*int(n/2) ) z = 0.0d0
         f1 = z
         do k = 2 , n
            pf = (2.0d0-1.0d0/k)*z*f1 - (1.0d0-1.0d0/k)*f0
            pd = k*(f1-z*pf)/(1.0d0-z*z)
            f0 = f1
            f1 = pf
         enddo
         if ( z/=0.0 ) then
            fd = pf/p
            q = 0.0d0
            do i = 1 , nr
               wp = 1.0d0
               do j = 1 , nr
                  if ( j/=i ) wp = wp*(z-x(j))
               enddo
               q = q + wp
            enddo
            gd = (pd-q*fd)/p
            z = z - fd/gd
            if ( abs(z-z0)>abs(z)*1.0d-15 ) goto 50
         endif
         x(nr) = z
         x(n+1-nr) = -z
         w(nr) = 2.0d0/((1.0d0-z*z)*pd*pd)
         w(n+1-nr) = w(nr)
      enddo
      end

!       **********************************

      subroutine aswfa(m,n,c,x,Kd,Cv,S1f,S1d)
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
      implicit none
      real(wp) a0 , c , ck , Cv , d0 , d1 , df , eps , r , S1d ,&
                     & S1f , su1 , su2 , x , x0 , x1
      integer ip , k , Kd , m , n , nm , nm2
      dimension ck(200) , df(200)
      eps = 1.0d-14
      x0 = x
      x = abs(x)
      ip = 1
      if ( n-m==2*int((n-m)/2) ) ip = 0
      nm = 40 + int((n-m)/2+c)
      nm2 = nm/2 - 2
      call sdmn(m,n,c,Cv,Kd,df)
      call sckb(m,n,c,df,ck)
      x1 = 1.0d0 - x*x
      if ( m==0 .and. x1==0.0d0 ) then
         a0 = 1.0d0
      else
         a0 = x1**(0.5d0*m)
      endif
      su1 = ck(1)
      do k = 1 , nm2
         r = ck(k+1)*x1**k
         su1 = su1 + r
         if ( k>=10 .and. abs(r/su1)<eps ) exit
      enddo
      S1f = a0*x**ip*su1
      if ( x==1.0d0 ) then
         if ( m==0 ) S1d = ip*ck(1) - 2.0d0*ck(2)
         if ( m==1 ) S1d = -1.0d+100
         if ( m==2 ) S1d = -2.0d0*ck(1)
         if ( m>=3 ) S1d = 0.0d0
      else
         d0 = ip - m/x1*x**(ip+1.0d0)
         d1 = -2.0d0*a0*x**(ip+1.0d0)
         su2 = ck(2)
         do k = 2 , nm2
            r = k*ck(k+1)*x1**(k-1.0d0)
            su2 = su2 + r
            if ( k>=10 .and. abs(r/su2)<eps ) exit
         enddo
         S1d = d0*a0*su1 + d1*su2
      endif
      if ( x0<0.0d0 .and. ip==0 ) S1d = -S1d
      if ( x0<0.0d0 .and. ip==1 ) S1f = -S1f
      x = x0
      end



!       **********************************

      subroutine jyna(n,x,Nm,Bj,Dj,By,Dy)
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
      implicit none
      real(wp) Bj , bj0 , bj1 , bjk , By , by0 , by1 , cs , Dj ,&
                     & dj0 , dj1 , Dy , dy0 , dy1 , f , f0 , f1 , f2 , x
      integer k , m , n , Nm
      dimension Bj(0:n) , By(0:n) , Dj(0:n) , Dy(0:n)
      Nm = n
      if ( x<1.0d-100 ) then
         do k = 0 , n
            Bj(k) = 0.0d0
            Dj(k) = 0.0d0
            By(k) = -1.0d+300
            Dy(k) = 1.0d+300
         enddo
         Bj(0) = 1.0d0
         Dj(1) = 0.5d0
         return
      endif
      call jy01b(x,bj0,dj0,bj1,dj1,by0,dy0,by1,dy1)
      Bj(0) = bj0
      Bj(1) = bj1
      By(0) = by0
      By(1) = by1
      Dj(0) = dj0
      Dj(1) = dj1
      Dy(0) = dy0
      Dy(1) = dy1
      if ( n<=1 ) return
      if ( n<int(0.9*x) ) then
         do k = 2 , n
            bjk = 2.0d0*(k-1.0d0)/x*bj1 - bj0
            Bj(k) = bjk
            bj0 = bj1
            bj1 = bjk
         enddo
      else
         m = msta1(x,200)
         if ( m<n ) then
            Nm = m
         else
            m = msta2(x,n,15)
         endif
         f2 = 0.0d0
         f1 = 1.0d-100
         f = 0.0d0
         do k = m , 0 , -1
            f = 2.0d0*(k+1.0d0)/x*f1 - f2
            if ( k<=Nm ) Bj(k) = f
            f2 = f1
            f1 = f
         enddo
         if ( abs(bj0)>abs(bj1) ) then
            cs = bj0/f
         else
            cs = bj1/f2
         endif
         do k = 0 , Nm
            Bj(k) = cs*Bj(k)
         enddo
      endif
      do k = 2 , Nm
         Dj(k) = Bj(k-1) - k/x*Bj(k)
      enddo
      f0 = By(0)
      f1 = By(1)
      do k = 2 , Nm
         f = 2.0d0*(k-1.0d0)/x*f1 - f0
         By(k) = f
         f0 = f1
         f1 = f
      enddo
      do k = 2 , Nm
         Dy(k) = By(k-1) - k*By(k)/x
      enddo
      end



!       **********************************

      subroutine pbdv(v,x,Dv,Dp,Pdf,Pdd)
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
      implicit none
      real(wp) Dp , Dv , ep , f , f0 , f1 , pd , pd0 , pd1 ,    &
                     & Pdd , Pdf , s0 , v , v0 , v1 , v2 , vh , x , xa
      integer ja , k , l , m , na , nk , nv
      dimension Dv(0:*) , Dp(0:*)
      xa = abs(x)
      vh = v
      v = v + sign(1.0d0,v)
      nv = int(v)
      v0 = v - nv
      na = abs(nv)
      ep = exp(-.25d0*x*x)
      ja = 0
      if ( na>=1 ) ja = 1
      if ( v>=0.0 ) then
         if ( v0==0.0 ) then
            pd0 = ep
            pd1 = x*ep
         else
            do l = 0 , ja
               v1 = v0 + l
               if ( xa<=5.8 ) call dvsa(v1,x,pd1)
               if ( xa>5.8 ) call dvla(v1,x,pd1)
               if ( l==0 ) pd0 = pd1
            enddo
         endif
         Dv(0) = pd0
         Dv(1) = pd1
         do k = 2 , na
            Pdf = x*pd1 - (k+v0-1.0d0)*pd0
            Dv(k) = Pdf
            pd0 = pd1
            pd1 = Pdf
         enddo
      elseif ( x<=0.0 ) then
         if ( xa<=5.8d0 ) then
            call dvsa(v0,x,pd0)
            v1 = v0 - 1.0d0
            call dvsa(v1,x,pd1)
         else
            call dvla(v0,x,pd0)
            v1 = v0 - 1.0d0
            call dvla(v1,x,pd1)
         endif
         Dv(0) = pd0
         Dv(1) = pd1
         do k = 2 , na
            pd = (-x*pd1+pd0)/(k-1.0d0-v0)
            Dv(k) = pd
            pd0 = pd1
            pd1 = pd
         enddo
      elseif ( x<=2.0 ) then
         v2 = nv + v0
         if ( nv==0 ) v2 = v2 - 1.0d0
         nk = int(-v2)
         call dvsa(v2,x,f1)
         v1 = v2 + 1.0d0
         call dvsa(v1,x,f0)
         Dv(nk) = f1
         Dv(nk-1) = f0
         do k = nk - 2 , 0 , -1
            f = x*f0 + (k-v0+1.0d0)*f1
            Dv(k) = f
            f1 = f0
            f0 = f
         enddo
      else
         if ( xa<=5.8 ) call dvsa(v0,x,pd0)
         if ( xa>5.8 ) call dvla(v0,x,pd0)
         Dv(0) = pd0
         m = 100 + na
         f1 = 0.0d0
         f0 = 1.0d-30
         f = 0.0d0
         do k = m , 0 , -1
            f = x*f0 + (k-v0+1.0d0)*f1
            if ( k<=na ) Dv(k) = f
            f1 = f0
            f0 = f
         enddo
         s0 = pd0/f
         do k = 0 , na
            Dv(k) = s0*Dv(k)
         enddo
      endif
      do k = 0 , na - 1
         v1 = abs(v0) + k
         if ( v>=0.0d0 ) then
            Dp(k) = 0.5d0*x*Dv(k) - Dv(k+1)
         else
            Dp(k) = -0.5d0*x*Dv(k) - v1*Dv(k+1)
         endif
      enddo
      Pdf = Dv(na-1)
      Pdd = Dp(na-1)
      v = vh
      end



!       **********************************

      subroutine itsh0(x,Th0)
!
!       ===================================================
!       Purpose: Evaluate the integral of Struve function
!                H0(t) with respect to t from 0 and x
!       Input :  x   --- Upper limit  ( x ≥ 0 )
!       Output:  TH0 --- Integration of H0(t) from 0 and x
!       ===================================================
!
      implicit none
      real(wp) a , a0 , a1 , af , bf , bg , el , pi , r , rd ,  &
                     & s , s0 , Th0 , ty , x , xp
      integer k
      dimension a(25)
      pi = 3.141592653589793d0
      r = 1.0d0
      if ( x<=30.0 ) then
         s = 0.5d0
         do k = 1 , 100
            rd = 1.0d0
            if ( k==1 ) rd = 0.5d0
            r = -r*rd*k/(k+1.0d0)*(x/(2.0d0*k+1.0d0))**2
            s = s + r
            if ( abs(r)<abs(s)*1.0d-12 ) exit
         enddo
         Th0 = 2.0d0/pi*x*x*s
      else
         s = 1.0d0
         do k = 1 , 12
            r = -r*k/(k+1.0d0)*((2.0d0*k+1.0d0)/x)**2
            s = s + r
            if ( abs(r)<abs(s)*1.0d-12 ) exit
         enddo
         el = .57721566490153d0
         s0 = s/(pi*x*x) + 2.0d0/pi*(log(2.0d0*x)+el)
         a0 = 1.0d0
         a1 = 5.0d0/8.0d0
         a(1) = a1
         do k = 1 , 20
            af = ((1.5d0*(k+.5d0)*(k+5.0d0/6.0d0)*a1-.5d0*(k+.5d0)      &
               & *(k+.5d0)*(k-.5d0)*a0))/(k+1.0d0)
            a(k+1) = af
            a0 = a1
            a1 = af
         enddo
         bf = 1.0d0
         r = 1.0d0
         do k = 1 , 10
            r = -r/(x*x)
            bf = bf + a(2*k)*r
         enddo
         bg = a(1)/x
         r = 1.0d0/x
         do k = 1 , 10
            r = -r/(x*x)
            bg = bg + a(2*k+1)*r
         enddo
         xp = x + .25d0*pi
         ty = sqrt(2.0d0/(pi*x))*(bg*cos(xp)-bf*sin(xp))
         Th0 = ty + s0
      endif
      end

!       **********************************

      subroutine cerzo(Nt,Zo)
!
!       ===============================================================
!       Purpose : Evaluate the complex zeros of error function erf(z)
!                 using the modified Newton's iteration method
!       Input :   NT --- Total number of zeros
!       Output:   ZO(L) --- L-th zero of erf(z), L=1,2,...,NT
!       Routine called: CERF for computing erf(z) and erf'(z)
!       ===============================================================
!
      implicit none
      integer i , it , j , nr , Nt
      real(wp) pi , pu , pv , px , py , w , w0
      complex(wp) z , zd , zf , zfd , zgd , Zo , zp , zq , zw
      dimension Zo(Nt)
      pi = 3.141592653589793d0
      w = 0.0d0
      do nr = 1 , Nt
         pu = sqrt(pi*(4.0d0*nr-0.5d0))
         pv = pi*sqrt(2.0d0*nr-0.25d0)
         px = 0.5*pu - 0.5*log(pv)/pu
         py = 0.5*pu + 0.5*log(pv)/pu
         z = dcmplx(px,py)
         it = 0
 50      it = it + 1
         call cerf(z,zf,zd)
         zp = (1.0d0,0.0d0)
         do i = 1 , nr - 1
            zp = zp*(z-Zo(i))
         enddo
         zfd = zf/zp
         zq = (0.0d0,0.0d0)
         do i = 1 , nr - 1
            zw = (1.0d0,0.0d0)
            do j = 1 , nr - 1
               if ( j/=i ) zw = zw*(z-Zo(j))
            enddo
            zq = zq + zw
         enddo
         zgd = (zd-zq*zfd)/zp
         z = z - zfd/zgd
         w0 = w
         w = abs(z)
         if ( it<=50 .and. abs((w-w0)/w)>1.0d-11 ) goto 50
         Zo(nr) = z
      enddo
      end



!       **********************************

      subroutine gamma2(x,Ga)
!
!       ==================================================
!       Purpose: Compute gamma function Г(x)
!       Input :  x  --- Argument of Г(x)
!                       ( x is not equal to 0,-1,-2,…)
!       Output:  GA --- Г(x)
!       ==================================================
!
      implicit none
      real(wp) g , Ga , gr , pi , r , x , z
      integer k , m , m1
      dimension g(26)
      pi = 3.141592653589793d0
      if ( x/=int(x) ) then
         r = 1.0d0
         if ( abs(x)>1.0d0 ) then
            z = abs(x)
            m = int(z)
            do k = 1 , m
               r = r*(z-k)
            enddo
            z = z - m
         else
            z = x
         endif
         data g/1.0d0 , 0.5772156649015329d0 , -0.6558780715202538d0 ,  &
            & -0.420026350340952d-1 , 0.1665386113822915d0 ,            &
            & -.421977345555443d-1 , -.96219715278770d-2 ,              &
            & .72189432466630d-2 , -.11651675918591d-2 ,                &
            & -.2152416741149d-3 , .1280502823882d-3 ,                  &
            & -.201348547807d-4 , -.12504934821d-5 , .11330272320d-5 ,  &
            & -.2056338417d-6 , .61160950d-8 , .50020075d-8 ,           &
            & -.11812746d-8 , .1043427d-9 , .77823d-11 , -.36968d-11 ,  &
            & .51d-12 , -.206d-13 , -.54d-14 , .14d-14 , .1d-15/
         gr = g(26)
         do k = 25 , 1 , -1
            gr = gr*z + g(k)
         enddo
         Ga = 1.0d0/(gr*z)
         if ( abs(x)>1.0d0 ) then
            Ga = Ga*r
            if ( x<0.0d0 ) Ga = -pi/(x*Ga*sin(pi*x))
         endif
      elseif ( x>0.0d0 ) then
         Ga = 1.0d0
         m1 = x - 1
         do k = 2 , m1
            Ga = Ga*k
         enddo
      else
         Ga = 1.0d+300
      endif
      end

!       **********************************

      subroutine chgu(a,b,x,Hu,Md,Isfer)
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
      implicit none
      real(wp) a , a00 , aa , b , b00 , Hu , hu1 , x
      integer id , id1 , Isfer , Md
      logical il1 , il2 , il3 , bl1 , bl2 , bl3 , bn
      aa = a - b + 1.0d0
      Isfer = 0
      il1 = a==int(a) .and. a<=0.0
      il2 = aa==int(aa) .and. aa<=0.0
      il3 = abs(a*(a-b+1.0))/x<=2.0
      bl1 = x<=5.0 .or. (x<=10.0 .and. a<=2.0)
      bl2 = (x>5.0 .and. x<=12.5) .and. (a>=1.0 .and. b>=a+4.0)
      bl3 = x>12.5 .and. a>=5.0 .and. b>=a + 5.0
      bn = b==int(b) .and. b/=0.0
      id1 = -100
      hu1 = 0.0d0
      if ( b/=int(b) ) then
         call chgus(a,b,x,Hu,id1)
         Md = 1
         if ( id1>=9 ) return
         hu1 = Hu
      endif
      if ( il1 .or. il2 .or. il3 ) then
         call chgul(a,b,x,Hu,id)
         Md = 2
         if ( id>=9 ) return
         if ( id1>id ) then
            Md = 1
            id = id1
            Hu = hu1
         endif
      endif
      if ( a>=1.0 ) then
         if ( bn .and. (bl1 .or. bl2 .or. bl3) ) then
            call chgubi(a,b,x,Hu,id)
            Md = 3
         else
            call chguit(a,b,x,Hu,id)
            Md = 4
         endif
      elseif ( b<=a ) then
         a00 = a
         b00 = b
         a = a - b + 1.0d0
         b = 2.0d0 - b
         call chguit(a,b,x,Hu,id)
         Hu = x**(1.0d0-b00)*Hu
         a = a00
         b = b00
         Md = 4
      elseif ( bn .and. (.not.il1) ) then
         call chgubi(a,b,x,Hu,id)
         Md = 3
      endif
      if ( id<6 ) Isfer = 6
      end



!       **********************************

      subroutine lamn(n,x,Nm,Bl,Dl)
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
      implicit none
      real(wp) bg , bk , Bl , bs , Dl , f , f0 , f1 , r , r0 ,  &
                     & uk , x , x2
      integer i , k , m , n , Nm
      dimension Bl(0:n) , Dl(0:n)
      Nm = n
      if ( abs(x)<1.0d-100 ) then
         do k = 0 , n
            Bl(k) = 0.0d0
            Dl(k) = 0.0d0
         enddo
         Bl(0) = 1.0d0
         Dl(1) = 0.5d0
         return
      endif
      if ( x<=12.0d0 ) then
         x2 = x*x
         do k = 0 , n
            bk = 1.0d0
            r = 1.0d0
            do i = 1 , 50
               r = -0.25d0*r*x2/(i*(i+k))
               bk = bk + r
               if ( abs(r)<abs(bk)*1.0d-15 ) exit
            enddo
            Bl(k) = bk
            if ( k>=1 ) Dl(k-1) = -0.5d0*x/k*bk
         enddo
         uk = 1.0d0
         r = 1.0d0
         do i = 1 , 50
            r = -0.25d0*r*x2/(i*(i+n+1.0d0))
            uk = uk + r
            if ( abs(r)<abs(uk)*1.0d-15 ) exit
         enddo
         Dl(n) = -0.5d0*x/(n+1.0d0)*uk
         return
      endif
      if ( n==0 ) Nm = 1
      m = msta1(x,200)
      if ( m<Nm ) then
         Nm = m
      else
         m = msta2(x,Nm,15)
      endif
      bs = 0.0d0
      f = 0.0d0
      f0 = 0.0d0
      f1 = 1.0d-100
      do k = m , 0 , -1
         f = 2.0d0*(k+1.0d0)*f1/x - f0
         if ( k<=Nm ) Bl(k) = f
         if ( k==2*int(k/2) ) bs = bs + 2.0d0*f
         f0 = f1
         f1 = f
      enddo
      bg = bs - f
      do k = 0 , Nm
         Bl(k) = Bl(k)/bg
      enddo
      r0 = 1.0d0
      do k = 1 , Nm
         r0 = 2.0d0*r0*k/x
         Bl(k) = r0*Bl(k)
      enddo
      Dl(0) = -0.5d0*x*Bl(1)
      do k = 1 , Nm
         Dl(k) = 2.0d0*k/x*(Bl(k-1)-Bl(k))
      enddo
      end



!       **********************************

      subroutine comelp(Hk,Ck,Ce)
!
!       ==================================================
!       Purpose: Compute complete elliptic integrals K(k)
!                and E(k)
!       Input  : K  --- Modulus k ( 0 ≤ k ≤ 1 )
!       Output : CK --- K(k)
!                CE --- E(k)
!       ==================================================
!
      implicit none
      real(wp) ae , ak , be , bk , Ce , Ck , Hk , pk
      pk = 1.0d0 - Hk*Hk
      if ( Hk==1.0 ) then
         Ck = 1.0d+300
         Ce = 1.0d0
      else
         ak = (((.01451196212d0*pk+.03742563713d0)*pk+.03590092383d0)   &
            & *pk+.09666344259d0)*pk + 1.38629436112d0
         bk = (((.00441787012d0*pk+.03328355346d0)*pk+.06880248576d0)   &
            & *pk+.12498593597d0)*pk + .5d0
         Ck = ak - bk*log(pk)
         ae = (((.01736506451d0*pk+.04757383546d0)*pk+.0626060122d0)    &
            & *pk+.44325141463d0)*pk + 1.0d0
         be = (((.00526449639d0*pk+.04069697526d0)*pk+.09200180037d0)   &
            & *pk+.2499836831d0)*pk
         Ce = ae - be*log(pk)
      endif
      end

!       **********************************

      subroutine incob(a,b,x,Bix)
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
      implicit none
      real(wp) a , b , Bix , bt , dk , fk , s0 , t1 , t2 , ta , &
                     & tb , x
      integer k
      dimension dk(51) , fk(51)
      s0 = (a+1.0d0)/(a+b+2.0d0)
      call beta(a,b,bt)
      if ( x<=s0 ) then
         do k = 1 , 20
            dk(2*k) = k*(b-k)*x/(a+2.0d0*k-1.0d0)/(a+2.0d0*k)
         enddo
         do k = 0 , 20
            dk(2*k+1) = -(a+k)*(a+b+k)*x/(a+2.d0*k)/(a+2.0*k+1.0)
         enddo
         t1 = 0.0d0
         do k = 20 , 1 , -1
            t1 = dk(k)/(1.0d0+t1)
         enddo
         ta = 1.0d0/(1.0d0+t1)
         Bix = x**a*(1.0d0-x)**b/(a*bt)*ta
      else
         do k = 1 , 20
            fk(2*k) = k*(a-k)*(1.0d0-x)/(b+2.*k-1.0)/(b+2.0*k)
         enddo
         do k = 0 , 20
            fk(2*k+1) = -(b+k)*(a+b+k)*(1.d0-x)/(b+2.d0*k)              &
                      & /(b+2.d0*k+1.d0)
         enddo
         t2 = 0.0d0
         do k = 20 , 1 , -1
            t2 = fk(k)/(1.0d0+t2)
         enddo
         tb = 1.0d0/(1.0d0+t2)
         Bix = 1.0d0 - x**a*(1.0d0-x)**b/(b*bt)*tb
      endif
      end



!       **********************************

      subroutine cvf(Kd,m,q,a,Mj,f)
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
      implicit none
      real(wp) a , b , f , q , t0 , t1 , t2
      integer ic , j , j0 , jf , Kd , l , l0 , m , Mj
      b = a
      ic = int(m/2)
      l = 0
      l0 = 0
      j0 = 2
      jf = ic
      if ( Kd==1 ) l0 = 2
      if ( Kd==1 ) j0 = 3
      if ( Kd==2 .or. Kd==3 ) l = 1
      if ( Kd==4 ) jf = ic - 1
      t1 = 0.0d0
      do j = Mj , ic + 1 , -1
         t1 = -q*q/((2.0d0*j+l)**2-b+t1)
      enddo
      if ( m<=2 ) then
         t2 = 0.0d0
         if ( Kd==1 .and. m==0 ) t1 = t1 + t1
         if ( Kd==1 .and. m==2 ) t1 = -2.0d0*q*q/(4.0d0-b+t1) - 4.0d0
         if ( Kd==2 .and. m==1 ) t1 = t1 + q
         if ( Kd==3 .and. m==1 ) t1 = t1 - q
      else
         t0 = 0.0d0
         if ( Kd==1 ) t0 = 4.0d0 - b + 2.0d0*q*q/b
         if ( Kd==2 ) t0 = 1.0d0 - b + q
         if ( Kd==3 ) t0 = 1.0d0 - b - q
         if ( Kd==4 ) t0 = 4.0d0 - b
         t2 = -q*q/t0
         do j = j0 , jf
            t2 = -q*q/((2.0d0*j-l-l0)**2-b+t2)
         enddo
      endif
      f = (2.0d0*ic+l)**2 + t1 + t2 - b
      end



!       **********************************

      subroutine clpn(n,x,y,Cpn,Cpd)
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
      implicit none
      complex(wp) cp0 , cp1 , Cpd , cpf , Cpn , z
      integer k , n
      real(wp) x , y
      dimension Cpn(0:n) , Cpd(0:n)
      z = dcmplx(x,y)
      Cpn(0) = (1.0d0,0.0d0)
      Cpn(1) = z
      Cpd(0) = (0.0d0,0.0d0)
      Cpd(1) = (1.0d0,0.0d0)
      cp0 = (1.0d0,0.0d0)
      cp1 = z
      do k = 2 , n
         cpf = (2.0d0*k-1.0d0)/k*z*cp1 - (k-1.0d0)/k*cp0
         Cpn(k) = cpf
         if ( abs(x)==1.0d0 .and. y==0.0d0 ) then
            Cpd(k) = 0.5d0*x**(k+1)*k*(k+1.0d0)
         else
            Cpd(k) = k*(cp1-z*cpf)/(1.0d0-z*z)
         endif
         cp0 = cp1
         cp1 = cpf
      enddo
      end

!       **********************************

      subroutine lqmns(m,n,x,Qm,Qd)
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
      implicit none
      integer k , km , l , ls , m , n
      real(wp) q0 , q00 , q01 , q0l , q10 , q11 , q1l , Qd ,    &
                     & qf0 , qf1 , qf2 , qg0 , qg1 , qh0 , qh1 , qh2 ,  &
                     & Qm , qm0 , qm1 , qmk
      real(wp) x , xq
      dimension Qm(0:n) , Qd(0:n)
      do k = 0 , n
         Qm(k) = 0.0d0
         Qd(k) = 0.0d0
      enddo
      if ( abs(x)==1.0d0 ) then
         do k = 0 , n
            Qm(k) = 1.0d+300
            Qd(k) = 1.0d+300
         enddo
         return
      endif
      ls = 1
      if ( abs(x)>1.0d0 ) ls = -1
      xq = sqrt(ls*(1.0d0-x*x))
      q0 = 0.5d0*log(abs((x+1.0)/(x-1.0)))
      q00 = q0
      q10 = -1.0d0/xq
      q01 = x*q0 - 1.0d0
      q11 = -ls*xq*(q0+x/(1.0d0-x*x))
      qf0 = q00
      qf1 = q10
      qm0 = 0.0d0
      qm1 = 0.0d0
      do k = 2 , m
         qm0 = -2.0d0*(k-1.0)/xq*x*qf1 - ls*(k-1.0)*(2.0-k)*qf0
         qf0 = qf1
         qf1 = qm0
      enddo
      if ( m==0 ) qm0 = q00
      if ( m==1 ) qm0 = q10
      Qm(0) = qm0
      if ( abs(x)<1.0001d0 ) then
         if ( m==0 .and. n>0 ) then
            qf0 = q00
            qf1 = q01
            do k = 2 , n
               qf2 = ((2.0*k-1.0d0)*x*qf1-(k-1.0)*qf0)/k
               Qm(k) = qf2
               qf0 = qf1
               qf1 = qf2
            enddo
         endif
         qg0 = q01
         qg1 = q11
         do k = 2 , m
            qm1 = -2.0d0*(k-1.0)/xq*x*qg1 - ls*k*(3.0-k)*qg0
            qg0 = qg1
            qg1 = qm1
         enddo
         if ( m==0 ) qm1 = q01
         if ( m==1 ) qm1 = q11
         Qm(1) = qm1
         if ( m==1 .and. n>1 ) then
            qh0 = q10
            qh1 = q11
            do k = 2 , n
               qh2 = ((2.0*k-1.0d0)*x*qh1-k*qh0)/(k-1.0)
               Qm(k) = qh2
               qh0 = qh1
               qh1 = qh2
            enddo
         elseif ( m>=2 ) then
            qg0 = q00
            qg1 = q01
            qh0 = q10
            qh1 = q11
            qmk = 0.0d0
            do l = 2 , n
               q0l = ((2.0d0*l-1.0d0)*x*qg1-(l-1.0d0)*qg0)/l
               q1l = ((2.0*l-1.0d0)*x*qh1-l*qh0)/(l-1.0d0)
               qf0 = q0l
               qf1 = q1l
               do k = 2 , m
                  qmk = -2.0d0*(k-1.0)/xq*x*qf1 - ls*(k+l-1.0)*(l+2.0-k)&
                      & *qf0
                  qf0 = qf1
                  qf1 = qmk
               enddo
               Qm(l) = qmk
               qg0 = qg1
               qg1 = q0l
               qh0 = qh1
               qh1 = q1l
            enddo
         endif
      else
         if ( abs(x)>1.1 ) then
            km = 40 + m + n
         else
            km = (40+m+n)*int(-1.0-1.8*log(x-1.0))
         endif
         qf2 = 0.0d0
         qf1 = 1.0d0
         do k = km , 0 , -1
            qf0 = ((2.0*k+3.0d0)*x*qf1-(k+2.0-m)*qf2)/(k+m+1.0)
            if ( k<=n ) Qm(k) = qf0
            qf2 = qf1
            qf1 = qf0
         enddo
         do k = 0 , n
            Qm(k) = Qm(k)*qm0/qf0
         enddo
      endif
      if ( abs(x)<1.0d0 ) then
         do k = 0 , n
            Qm(k) = (-1)**m*Qm(k)
         enddo
      endif
      Qd(0) = ((1.0d0-m)*Qm(1)-x*Qm(0))/(x*x-1.0)
      do k = 1 , n
         Qd(k) = (k*x*Qm(k)-(k+m)*Qm(k-1))/(x*x-1.0)
      enddo
      end

!       **********************************

      subroutine ciklv(v,z,Cbiv,Cdiv,Cbkv,Cdkv)
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
      implicit none
      real(wp) a , pi , v , v0 , vr
      complex(wp) Cbiv , Cbkv , Cdiv , Cdkv , ceta , cf , cfi , cfk ,    &
               & csi , csk , ct , ct2 , cws , z
      integer i , k , km , l , l0 , lf
      dimension cf(12) , a(91)
      pi = 3.141592653589793d0
      km = 12
      call cjk(km,a)
      do l = 1 , 0 , -1
         v0 = v - l
         cws = sqrt(1.0d0+(z/v0)*(z/v0))
         ceta = cws + log(z/v0/(1.0d0+cws))
         ct = 1.0d0/cws
         ct2 = ct*ct
         do k = 1 , km
            l0 = k*(k+1)/2 + 1
            lf = l0 + k
            cf(k) = a(lf)
            do i = lf - 1 , l0 , -1
               cf(k) = cf(k)*ct2 + a(i)
            enddo
            cf(k) = cf(k)*ct**k
         enddo
         vr = 1.0d0/v0
         csi = (1.0d0,0.0d0)
         do k = 1 , km
            csi = csi + cf(k)*vr**k
         enddo
         Cbiv = sqrt(ct/(2.0d0*pi*v0))*exp(v0*ceta)*csi
         if ( l==1 ) cfi = Cbiv
         csk = (1.0d0,0.0d0)
         do k = 1 , km
            csk = csk + (-1)**k*cf(k)*vr**k
         enddo
         Cbkv = sqrt(pi*ct/(2.0d0*v0))*exp(-v0*ceta)*csk
         if ( l==1 ) cfk = Cbkv
      enddo
      Cdiv = cfi - v/z*Cbiv
      Cdkv = -cfk - v/z*Cbkv
      end



!       **********************************

      subroutine elit(Hk,Phi,Fe,Ee)
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
      implicit none
      real(wp) a , a0 , b , b0 , c , ce , ck , d , d0 , Ee ,    &
                     & fac , Fe , g , Hk , Phi , pi , r
      integer n
      g = 0.0d0
      pi = 3.14159265358979d0
      a0 = 1.0d0
      b0 = sqrt(1.0d0-Hk*Hk)
      d0 = (pi/180.0d0)*Phi
      r = Hk*Hk
      if ( Hk==1.0d0 .and. Phi==90.0d0 ) then
         Fe = 1.0d+300
         Ee = 1.0d0
      elseif ( Hk==1.0d0 ) then
         Fe = log((1.0d0+sin(d0))/cos(d0))
         Ee = sin(d0)
      else
         fac = 1.0d0
         d = 0.0d0
         do n = 1 , 40
            a = (a0+b0)/2.0d0
            b = sqrt(a0*b0)
            c = (a0-b0)/2.0d0
            fac = 2.0d0*fac
            r = r + fac*c*c
            if ( Phi/=90.0d0 ) then
               d = d0 + atan((b0/a0)*tan(d0))
               g = g + c*sin(d)
               d0 = d + pi*int(d/pi+.5d0)
            endif
            a0 = a
            b0 = b
            if ( c<1.0d-7 ) exit
         enddo
         ck = pi/(2.0d0*a)
         ce = pi*(2.0d0-r)/(4.0d0*a)
         if ( Phi==90.0d0 ) then
            Fe = ck
            Ee = ce
         else
            Fe = d/(fac*a)
            Ee = Fe*ce/ck + g
         endif
      endif
      end

!       **********************************

      subroutine elit3(Phi,Hk,c,El3)
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
      implicit none
      real(wp) c , c0 , c1 , c2 , El3 , f1 , f2 , Hk , Phi , t ,&
                     & t1 , t2 , w
      integer i
      dimension t(10) , w(10)
      logical lb1 , lb2
      data t/.9931285991850949d0 , .9639719272779138d0 ,                &
         & .9122344282513259d0 , .8391169718222188d0 ,                  &
         & .7463319064601508d0 , .6360536807265150d0 ,                  &
         & .5108670019508271d0 , .3737060887154195d0 ,                  &
         & .2277858511416451d0 , .7652652113349734d-1/
      data w/.1761400713915212d-1 , .4060142980038694d-1 ,              &
         & .6267204833410907d-1 , .8327674157670475d-1 ,                &
         & .1019301198172404d0 , .1181945319615184d0 ,                  &
         & .1316886384491766d0 , .1420961093183820d0 ,                  &
         & .1491729864726037d0 , .1527533871307258d0/
      lb1 = Hk==1.0d0 .and. abs(Phi-90.0)<=1.0d-8
      lb2 = c==1.0d0 .and. abs(Phi-90.0)<=1.0d-8
      if ( lb1 .or. lb2 ) then
         El3 = 1.0d+300
         return
      endif
      c1 = 0.87266462599716d-2*Phi
      c2 = c1
      El3 = 0.0d0
      do i = 1 , 10
         c0 = c2*t(i)
         t1 = c1 + c0
         t2 = c1 - c0
         f1 = 1.0d0/((1.0d0-c*sin(t1)*sin(t1))                        &
            & *sqrt(1.0d0-Hk*Hk*sin(t1)*sin(t1)))
         f2 = 1.0d0/((1.0d0-c*sin(t2)*sin(t2))                        &
            & *sqrt(1.0d0-Hk*Hk*sin(t2)*sin(t2)))
         El3 = El3 + w(i)*(f1+f2)
      enddo
      El3 = c1*El3
      end

!       **********************************

      subroutine eix(x,Ei)
!
!       ============================================
!       Purpose: Compute exponential integral Ei(x)
!       Input :  x  --- Argument of Ei(x)
!       Output:  EI --- Ei(x)
!       ============================================
!
      implicit none
      real(wp) Ei , ga , r , x
      integer k
      if ( x==0.0 ) then
         Ei = -1.0d+300
      elseif ( x<0 ) then
         call e1xb(-x,Ei)
         Ei = -Ei
      elseif ( abs(x)<=40.0 ) then
!          Power series around x=0
         Ei = 1.0d0
         r = 1.0d0
         do k = 1 , 100
            r = r*k*x/(k+1.0d0)**2
            Ei = Ei + r
            if ( abs(r/Ei)<=1.0d-15 ) exit
         enddo
         ga = 0.5772156649015328d0
         Ei = ga + log(x) + x*Ei
      else
!          Asymptotic expansion (the series is not convergent)
         Ei = 1.0d0
         r = 1.0d0
         do k = 1 , 20
            r = r*k/x
            Ei = Ei + r
         enddo
         Ei = exp(x)/x*Ei
      endif
      end

!       **********************************

      subroutine eixz(z,Cei)
!
!       ============================================
!       Purpose: Compute exponential integral Ei(x)
!       Input :  x  --- Complex argument of Ei(x)
!       Output:  EI --- Ei(x)
!       ============================================
!
      implicit none
      complex(wp) z , Cei
      real(wp) pi
      pi = 3.141592653589793d0
      call e1z(-z,Cei)
      Cei = -Cei
      if ( dimag(z)>0 ) then
         Cei = Cei + (0d0,1d0)*pi
      elseif ( dimag(z)<0 ) then
         Cei = Cei - (0d0,1d0)*pi
      elseif ( dimag(z)==0 ) then
         if ( dble(z)>0 ) Cei = Cei + (0d0,1d0)*sign(pi,dimag(z))
      endif
      end

!       **********************************

      subroutine e1xb(x,e1)
!
!       ============================================
!       Purpose: Compute exponential integral E1(x)
!       Input :  x  --- Argument of E1(x)
!       Output:  E1 --- E1(x)  ( x > 0 )
!       ============================================
!
      implicit none
      real(wp) e1 , ga , r , t , t0 , x
      integer k , m
      if ( x==0.0 ) then
         e1 = 1.0d+300
      elseif ( x<=1.0 ) then
         e1 = 1.0d0
         r = 1.0d0
         do k = 1 , 25
            r = -r*k*x/(k+1.0d0)**2
            e1 = e1 + r
            if ( abs(r)<=abs(e1)*1.0d-15 ) exit
         enddo
         ga = 0.5772156649015328d0
         e1 = -ga - log(x) + x*e1
      else
         m = 20 + int(80.0/x)
         t0 = 0.0d0
         do k = m , 1 , -1
            t0 = k/(1.0d0+k/(x+t0))
         enddo
         t = 1.0d0/(x+t0)
         e1 = exp(-x)*t
      endif
      end

!       **********************************

      subroutine chgm(a,b,x,Hg)
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
      implicit none
      real(wp) a , a0 , a1 , b , Hg , hg1 , hg2 , pi , r1 , r2 ,&
                     & rg , sum1 , sum2 , tai , tar , tbai , tbar ,     &
                     & tbi , tbr , x
      real(wp) x0 , xg , y , y0 , y1
      complex(wp) cta , ctb , ctba
      integer i , j , la , n , nl
      pi = 3.141592653589793d0
      a0 = a
      a1 = a
      x0 = x
      Hg = 0.0d0
!       DLMF 13.2.39
      if ( x<0.0d0 ) then
         a = b - a
         a0 = a
         x = abs(x)
      endif
      nl = 0
      la = 0
      if ( a>=2.0d0 ) then
!       preparing terms for DLMF 13.3.1
         nl = 1
         la = int(a)
         a = a - la - 1.0d0
      endif
      y0 = 0.0d0
      y1 = 0.0d0
      do n = 0 , nl
         if ( a0>=2.0d0 ) a = a + 1.0d0
         if ( x<=30.0d0+abs(b) .or. a<0.0d0 ) then
            Hg = 1.0d0
            rg = 1.0d0
            do j = 1 , 500
               rg = rg*(a+j-1.0d0)/(j*(b+j-1.0d0))*x
               Hg = Hg + rg
               if ( Hg/=0d0 .and. abs(rg/Hg)<1.0d-15 ) then
!       DLMF 13.2.39 (cf. above)
                  if ( x0<0.0d0 ) Hg = Hg*exp(x0)
                  goto 50
               endif
            enddo
         else
!       DLMF 13.7.2 & 13.2.4, SUM2 corresponds to first sum
            y = 0.0d0
            call cgama(a,y,0,tar,tai)
            cta = dcmplx(tar,tai)
            y = 0.0d0
            call cgama(b,y,0,tbr,tbi)
            ctb = dcmplx(tbr,tbi)
            xg = b - a
            y = 0.0d0
            call cgama(xg,y,0,tbar,tbai)
            ctba = dcmplx(tbar,tbai)
            sum1 = 1.0d0
            sum2 = 1.0d0
            r1 = 1.0d0
            r2 = 1.0d0
            do i = 1 , 8
               r1 = -r1*(a+i-1.0d0)*(a-b+i)/(x*i)
               r2 = -r2*(b-a+i-1.0d0)*(a-i)/(x*i)
               sum1 = sum1 + r1
               sum2 = sum2 + r2
            enddo
            if ( x0>=0.0d0 ) then
               hg1 = dble(exp(ctb-ctba))*x**(-a)*cos(pi*a)*sum1
               hg2 = dble(exp(ctb-cta+x))*x**(a-b)*sum2
            else
!       DLMF 13.2.39 (cf. above)
               hg1 = dble(exp(ctb-ctba+x0))*x**(-a)*cos(pi*a)*sum1
               hg2 = dble(exp(ctb-cta))*x**(a-b)*sum2
            endif
            Hg = hg1 + hg2
         endif
 50      if ( n==0 ) y0 = Hg
         if ( n==1 ) y1 = Hg
      enddo
      if ( a0>=2.0d0 ) then
!       DLMF 13.3.1
         do i = 1 , la - 1
            Hg = ((2.0d0*a-b+x)*y1+(b-a)*y0)/a
            y0 = y1
            y1 = Hg
            a = a + 1.0d0
         enddo
      endif
      a = a1
      x = x0
      end

!       **********************************

      subroutine hygfx(a,b,c,x,Hf,Isfer)
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
      implicit none
      real(wp) a , a0 , aa , b , bb , c , c0 , c1 , el , eps ,  &
                     & f0 , f1 , g0 , g1 , g2 , g3 , ga , gabc , gam ,  &
                     & gb
      real(wp) gbm , gc , gca , gcab , gcb , gm , Hf , hw , pa ,&
                     & pb , pi , r , r0 , r1 , rm , rp , sm , sp , sp0 ,&
                     & x
      real(wp) x1
      integer Isfer , j , k , m , nm
      logical l0 , l1 , l2 , l3 , l4 , l5
      pi = 3.141592653589793d0
      el = .5772156649015329d0
      Isfer = 0
      l0 = c==int(c) .and. c<0.0
      l1 = 1.0d0 - x<1.0d-15 .and. c - a - b<=0.0
      l2 = a==int(a) .and. a<0.0
      l3 = b==int(b) .and. b<0.0
      l4 = c - a==int(c-a) .and. c - a<=0.0
      l5 = c - b==int(c-b) .and. c - b<=0.0
      if ( l0 .or. l1 ) then
         Isfer = 3
         return
      endif
      eps = 1.0d-15
      if ( x>0.95 ) eps = 1.0d-8
      if ( x==0.0 .or. a==0.0 .or. b==0.0 ) then
         Hf = 1.0d0
         return
      elseif ( 1.0d0-x==eps .and. c-a-b>0.0 ) then
         call gamma2(c,gc)
         call gamma2(c-a-b,gcab)
         call gamma2(c-a,gca)
         call gamma2(c-b,gcb)
         Hf = gc*gcab/(gca*gcb)
         return
      elseif ( 1.0d0+x<=eps .and. abs(c-a+b-1.0)<=eps ) then
         g0 = sqrt(pi)*2.0d0**(-a)
         call gamma2(c,g1)
         call gamma2(1.0d0+a/2.0-b,g2)
         call gamma2(0.5d0+0.5*a,g3)
         Hf = g0*g1/(g2*g3)
         return
      elseif ( l2 .or. l3 ) then
         if ( l2 ) nm = int(abs(a))
         if ( l3 ) nm = int(abs(b))
         Hf = 1.0d0
         r = 1.0d0
         do k = 1 , nm
            r = r*(a+k-1.0d0)*(b+k-1.0d0)/(k*(c+k-1.0d0))*x
            Hf = Hf + r
         enddo
         return
      elseif ( l4 .or. l5 ) then
         if ( l4 ) nm = int(abs(c-a))
         if ( l5 ) nm = int(abs(c-b))
         Hf = 1.0d0
         r = 1.0d0
         do k = 1 , nm
            r = r*(c-a+k-1.0d0)*(c-b+k-1.0d0)/(k*(c+k-1.0d0))*x
            Hf = Hf + r
         enddo
         Hf = (1.0d0-x)**(c-a-b)*Hf
         return
      endif
      aa = a
      bb = b
      x1 = x
      if ( x<0.0d0 ) then
         x = x/(x-1.0d0)
         if ( c>a .and. b<a .and. b>0.0 ) then
            a = bb
            b = aa
         endif
         b = c - b
      endif
      hw = 0.0d0
      if ( x>=0.75d0 ) then
         gm = 0.0d0
         if ( abs(c-a-b-int(c-a-b))<1.0d-15 ) then
            m = int(c-a-b)
            call gamma2(a,ga)
            call gamma2(b,gb)
            call gamma2(c,gc)
            call gamma2(a+m,gam)
            call gamma2(b+m,gbm)
            call psi_spec(a,pa)
            call psi_spec(b,pb)
            if ( m/=0 ) gm = 1.0d0
            do j = 1 , abs(m) - 1
               gm = gm*j
            enddo
            rm = 1.0d0
            do j = 1 , abs(m)
               rm = rm*j
            enddo
            f0 = 1.0d0
            r0 = 1.0d0
            r1 = 1.0d0
            sp0 = 0.d0
            sp = 0.0d0
            if ( m>=0 ) then
               c0 = gm*gc/(gam*gbm)
               c1 = -gc*(x-1.0d0)**m/(ga*gb*rm)
               do k = 1 , m - 1
                  r0 = r0*(a+k-1.0d0)*(b+k-1.0)/(k*(k-m))*(1.0-x)
                  f0 = f0 + r0
               enddo
               do k = 1 , m
                  sp0 = sp0 + 1.0d0/(a+k-1.0) + 1.0/(b+k-1.0) - 1.0/k
               enddo
               f1 = pa + pb + sp0 + 2.0d0*el + log(1.0d0-x)
               do k = 1 , 250
                  sp = sp + (1.0d0-a)/(k*(a+k-1.0)) + (1.0-b)           &
                     & /(k*(b+k-1.0))
                  sm = 0.0d0
                  do j = 1 , m
                     sm = sm + (1.0d0-a)/((j+k)*(a+j+k-1.0))            &
                        & + 1.0/(b+j+k-1.0)
                  enddo
                  rp = pa + pb + 2.0d0*el + sp + sm + log(1.0d0-x)
                  r1 = r1*(a+m+k-1.0d0)*(b+m+k-1.0)/(k*(m+k))*(1.0-x)
                  f1 = f1 + r1*rp
                  if ( abs(f1-hw)<abs(f1)*eps ) exit
                  hw = f1
               enddo
               Hf = f0*c0 + f1*c1
            elseif ( m<0 ) then
               m = -m
               c0 = gm*gc/(ga*gb*(1.0d0-x)**m)
               c1 = -(-1)**m*gc/(gam*gbm*rm)
               do k = 1 , m - 1
                  r0 = r0*(a-m+k-1.0d0)*(b-m+k-1.0)/(k*(k-m))*(1.0-x)
                  f0 = f0 + r0
               enddo
               do k = 1 , m
                  sp0 = sp0 + 1.0d0/k
               enddo
               f1 = pa + pb - sp0 + 2.0d0*el + log(1.0d0-x)
               do k = 1 , 250
                  sp = sp + (1.0d0-a)/(k*(a+k-1.0)) + (1.0-b)           &
                     & /(k*(b+k-1.0))
                  sm = 0.0d0
                  do j = 1 , m
                     sm = sm + 1.0d0/(j+k)
                  enddo
                  rp = pa + pb + 2.0d0*el + sp - sm + log(1.0d0-x)
                  r1 = r1*(a+k-1.0d0)*(b+k-1.0)/(k*(m+k))*(1.0-x)
                  f1 = f1 + r1*rp
                  if ( abs(f1-hw)<abs(f1)*eps ) exit
                  hw = f1
               enddo
               Hf = f0*c0 + f1*c1
            endif
         else
            call gamma2(a,ga)
            call gamma2(b,gb)
            call gamma2(c,gc)
            call gamma2(c-a,gca)
            call gamma2(c-b,gcb)
            call gamma2(c-a-b,gcab)
            call gamma2(a+b-c,gabc)
            c0 = gc*gcab/(gca*gcb)
            c1 = gc*gabc/(ga*gb)*(1.0d0-x)**(c-a-b)
            Hf = 0.0d0
            r0 = c0
            r1 = c1
            do k = 1 , 250
               r0 = r0*(a+k-1.0d0)*(b+k-1.0)/(k*(a+b-c+k))*(1.0-x)
               r1 = r1*(c-a+k-1.0d0)*(c-b+k-1.0)/(k*(c-a-b+k))*(1.0-x)
               Hf = Hf + r0 + r1
               if ( abs(Hf-hw)<abs(Hf)*eps ) exit
               hw = Hf
            enddo
            Hf = Hf + c0 + c1
         endif
      else
         a0 = 1.0d0
         if ( c>a .and. c<2.0d0*a .and. c>b .and. c<2.0d0*b ) then
            a0 = (1.0d0-x)**(c-a-b)
            a = c - a
            b = c - b
         endif
         Hf = 1.0d0
         r = 1.0d0
         do k = 1 , 250
            r = r*(a+k-1.0d0)*(b+k-1.0d0)/(k*(c+k-1.0d0))*x
            Hf = Hf + r
            if ( abs(Hf-hw)<=abs(Hf)*eps ) exit
            hw = Hf
         enddo
         Hf = a0*Hf
      endif
      if ( x1<0.0d0 ) then
         x = x1
         c0 = 1.0d0/(1.0d0-x)**aa
         Hf = c0*Hf
      endif
      a = aa
      b = bb
      if ( k>120 ) Isfer = 5
      end



!       **********************************

      subroutine cchg(a,b,z,Chg)
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
      implicit none
      real(wp) a , a0 , a1 , b , ba , g1i , g1r , g2i , g2r ,   &
                     & g3i , g3r , phi , pi , x , x0 , y
      complex(wp) cfac , cg1 , cg2 , cg3 , Chg , chg1 , chg2 , chw , ci ,&
               & cr , cr1 , cr2 , crg , cs1 , cs2 , cy0 , cy1 , z , z0
      integer i , j , k , la , m , n , nl , ns
      pi = 3.141592653589793d0
      ci = (0.0d0,1.0d0)
      a0 = a
      a1 = a
      z0 = z
      if ( b==0.0 .or. b==-int(abs(b)) ) then
         Chg = (1.0d+300,0.0d0)
      elseif ( a==0.0d0 .or. z==0.0d0 ) then
         Chg = (1.0d0,0.0d0)
      elseif ( a==-1.0d0 ) then
         Chg = 1.0d0 - z/b
      elseif ( a==b ) then
         Chg = exp(z)
      elseif ( a-b==1.0d0 ) then
         Chg = (1.0d0+z/b)*exp(z)
      elseif ( a==1.0d0 .and. b==2.0d0 ) then
         Chg = (exp(z)-1.0d0)/z
      elseif ( a==int(a) .and. a<0.0d0 ) then
         m = int(-a)
         cr = (1.0d0,0.0d0)
         Chg = (1.0d0,0.0d0)
         do k = 1 , m
            cr = cr*(a+k-1.0d0)/k/(b+k-1.0d0)*z
            Chg = Chg + cr
         enddo
      else
         x0 = dble(z)
         if ( x0<0.0d0 ) then
            a = b - a
            a0 = a
            z = -z
         endif
         nl = 0
         la = 0
         if ( a>=2.0d0 ) then
            nl = 1
            la = int(a)
            a = a - la - 1.0d0
         endif
         ns = 0
         do n = 0 , nl
            if ( a0>=2.0d0 ) a = a + 1.0d0
            if ( abs(z)<20.0d0+abs(b) .or. a<0.0d0 ) then
               Chg = (1.0d0,0.0d0)
               crg = (1.0d0,0.0d0)
               do j = 1 , 500
                  crg = crg*(a+j-1.0d0)/(j*(b+j-1.0d0))*z
                  Chg = Chg + crg
                  if ( abs((Chg-chw)/Chg)<1.d-15 ) goto 20
                  chw = Chg
               enddo
            else
               y = 0.0d0
               call cgama(a,y,0,g1r,g1i)
               cg1 = dcmplx(g1r,g1i)
               y = 0.0d0
               call cgama(b,y,0,g2r,g2i)
               cg2 = dcmplx(g2r,g2i)
               ba = b - a
               y = 0.0d0
               call cgama(ba,y,0,g3r,g3i)
               cg3 = dcmplx(g3r,g3i)
               cs1 = (1.0d0,0.0d0)
               cs2 = (1.0d0,0.0d0)
               cr1 = (1.0d0,0.0d0)
               cr2 = (1.0d0,0.0d0)
               do i = 1 , 8
                  cr1 = -cr1*(a+i-1.0d0)*(a-b+i)/(z*i)
                  cr2 = cr2*(b-a+i-1.0d0)*(i-a)/(z*i)
                  cs1 = cs1 + cr1
                  cs2 = cs2 + cr2
               enddo
               x = dble(z)
               y = dimag(z)
               if ( x==0.0 .and. y>=0.0 ) then
                  phi = 0.5d0*pi
               elseif ( x==0.0 .and. y<=0.0 ) then
                  phi = -0.5d0*pi
               else
                  phi = atan(y/x)
               endif
               if ( phi>-0.5*pi .and. phi<1.5*pi ) ns = 1
               if ( phi>-1.5*pi .and. phi<=-0.5*pi ) ns = -1
               cfac = exp(ns*ci*pi*a)
               if ( y==0.0d0 ) cfac = cos(pi*a)
               chg1 = exp(cg2-cg3)*z**(-a)*cfac*cs1
               chg2 = exp(cg2-cg1+z)*z**(a-b)*cs2
               Chg = chg1 + chg2
            endif
 20         if ( n==0 ) cy0 = Chg
            if ( n==1 ) cy1 = Chg
         enddo
         if ( a0>=2.0d0 ) then
            do i = 1 , la - 1
               Chg = ((2.0d0*a-b+z)*cy1+(b-a)*cy0)/a
               cy0 = cy1
               cy1 = Chg
               a = a + 1.0d0
            enddo
         endif
         if ( x0<0.0d0 ) Chg = Chg*exp(-z)
      endif
      a = a1
      z = z0
      end



!       **********************************

      subroutine hygfz(a,b,c,z,Zhf,Isfer)
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
      implicit none
      real(wp) a , a0 , aa , b , bb , c , ca , cb , el , eps ,  &
                     & g0 , g1 , g2 , g3 , ga , gab , gabc , gam , gb , &
                     & gba
      real(wp) gbm , gc , gca , gcab , gcb , gcbk , gm , pa ,   &
                     & pac , pb , pca , pi , rk1 , rk2 , rm , sj1 ,     &
                     & sj2 , sm , sp , sp0
      real(wp) sq , t0 , w0 , ws , x , y
      integer Isfer , j , k , m , mab , mcab , nca , ncb , nm
      complex(wp) z , z00 , z1 , zc0 , zc1 , zf0 , zf1 , Zhf , zp , zp0 ,&
               & zr , zr0 , zr1 , zw
      logical l0 , l1 , l2 , l3 , l4 , l5 , l6
      x = dble(z)
      y = dimag(z)
      eps = 1.0d-15
      Isfer = 0
      l0 = c==int(c) .and. c<0.0d0
      l1 = abs(1.0d0-x)<eps .and. y==0.0d0 .and. c - a - b<=0.0d0
      l2 = abs(z+1.0d0)<eps .and. abs(c-a+b-1.0d0)<eps
      l3 = a==int(a) .and. a<0.0d0
      l4 = b==int(b) .and. b<0.0d0
      l5 = c - a==int(c-a) .and. c - a<=0.0d0
      l6 = c - b==int(c-b) .and. c - b<=0.0d0
      aa = a
      bb = b
      a0 = abs(z)
      if ( a0>0.95d0 ) eps = 1.0d-8
      pi = 3.141592653589793d0
      el = .5772156649015329d0
      if ( l0 .or. l1 ) then
         Isfer = 3
         return
      endif
      nm = 0
      if ( a0==0.0d0 .or. a==0.0d0 .or. b==0.0d0 ) then
         Zhf = (1.0d0,0.0d0)
      elseif ( z==1.0d0 .and. c-a-b>0.0d0 ) then
         call gamma2(c,gc)
         call gamma2(c-a-b,gcab)
         call gamma2(c-a,gca)
         call gamma2(c-b,gcb)
         Zhf = gc*gcab/(gca*gcb)
      elseif ( l2 ) then
         g0 = sqrt(pi)*2.0d0**(-a)
         call gamma2(c,g1)
         call gamma2(1.0d0+a/2.0d0-b,g2)
         call gamma2(0.5d0+0.5d0*a,g3)
         Zhf = g0*g1/(g2*g3)
      elseif ( l3 .or. l4 ) then
         if ( l3 ) nm = int(abs(a))
         if ( l4 ) nm = int(abs(b))
         Zhf = (1.0d0,0.0d0)
         zr = (1.0d0,0.0d0)
         do k = 1 , nm
            zr = zr*(a+k-1.0d0)*(b+k-1.0d0)/(k*(c+k-1.0d0))*z
            Zhf = Zhf + zr
         enddo
      elseif ( l5 .or. l6 ) then
         if ( l5 ) nm = int(abs(c-a))
         if ( l6 ) nm = int(abs(c-b))
         Zhf = (1.0d0,0.0d0)
         zr = (1.0d0,0.0d0)
         do k = 1 , nm
            zr = zr*(c-a+k-1.0d0)*(c-b+k-1.0d0)/(k*(c+k-1.0d0))*z
            Zhf = Zhf + zr
         enddo
         Zhf = (1.0d0-z)**(c-a-b)*Zhf
      elseif ( a0<=1.0d0 ) then
         if ( x<0.0d0 ) then
            z1 = z/(z-1.0d0)
            if ( c>a .and. b<a .and. b>0.0 ) then
               a = bb
               b = aa
            endif
            zc0 = 1.0d0/((1.0d0-z)**a)
            Zhf = (1.0d0,0.0d0)
            zr0 = (1.0d0,0.0d0)
            do k = 1 , 500
               zr0 = zr0*(a+k-1.0d0)*(c-b+k-1.0d0)/(k*(c+k-1.0d0))*z1
               Zhf = Zhf + zr0
               if ( abs(Zhf-zw)<abs(Zhf)*eps ) exit
               zw = Zhf
            enddo
            Zhf = zc0*Zhf
         elseif ( a0>=0.90d0 ) then
            gm = 0.0d0
            mcab = int(c-a-b+eps*sign(1.0d0,c-a-b))
            if ( abs(c-a-b-mcab)<eps ) then
               m = int(c-a-b)
               call gamma2(a,ga)
               call gamma2(b,gb)
               call gamma2(c,gc)
               call gamma2(a+m,gam)
               call gamma2(b+m,gbm)
               call psi_spec(a,pa)
               call psi_spec(b,pb)
               if ( m/=0 ) gm = 1.0d0
               do j = 1 , abs(m) - 1
                  gm = gm*j
               enddo
               rm = 1.0d0
               do j = 1 , abs(m)
                  rm = rm*j
               enddo
               zf0 = (1.0d0,0.0d0)
               zr0 = (1.0d0,0.0d0)
               zr1 = (1.0d0,0.0d0)
               sp0 = 0.d0
               sp = 0.0d0
               if ( m>=0 ) then
                  zc0 = gm*gc/(gam*gbm)
                  zc1 = -gc*(z-1.0d0)**m/(ga*gb*rm)
                  do k = 1 , m - 1
                     zr0 = zr0*(a+k-1.d0)*(b+k-1.d0)/(k*(k-m))*(1.d0-z)
                     zf0 = zf0 + zr0
                  enddo
                  do k = 1 , m
                     sp0 = sp0 + 1.0d0/(a+k-1.0d0) + 1.0/(b+k-1.0d0)    &
                         & - 1.d0/k
                  enddo
                  zf1 = pa + pb + sp0 + 2.0d0*el + log(1.0d0-z)
                  do k = 1 , 500
                     sp = sp + (1.0d0-a)/(k*(a+k-1.0d0)) + (1.0d0-b)    &
                        & /(k*(b+k-1.0d0))
                     sm = 0.0d0
                     do j = 1 , m
                        sm = sm + (1.0d0-a)/((j+k)*(a+j+k-1.0d0))       &
                           & + 1.0d0/(b+j+k-1.0d0)
                     enddo
                     zp = pa + pb + 2.0d0*el + sp + sm + log(1.0d0-z)
                     zr1 = zr1*(a+m+k-1.0d0)*(b+m+k-1.0d0)/(k*(m+k))    &
                         & *(1.0d0-z)
                     zf1 = zf1 + zr1*zp
                     if ( abs(zf1-zw)<abs(zf1)*eps ) exit
                     zw = zf1
                  enddo
                  Zhf = zf0*zc0 + zf1*zc1
               elseif ( m<0 ) then
                  m = -m
                  zc0 = gm*gc/(ga*gb*(1.0d0-z)**m)
                  zc1 = -(-1)**m*gc/(gam*gbm*rm)
                  do k = 1 , m - 1
                     zr0 = zr0*(a-m+k-1.0d0)*(b-m+k-1.0d0)/(k*(k-m))    &
                         & *(1.0d0-z)
                     zf0 = zf0 + zr0
                  enddo
                  do k = 1 , m
                     sp0 = sp0 + 1.0d0/k
                  enddo
                  zf1 = pa + pb - sp0 + 2.0d0*el + log(1.0d0-z)
                  do k = 1 , 500
                     sp = sp + (1.0d0-a)/(k*(a+k-1.0d0)) + (1.0d0-b)    &
                        & /(k*(b+k-1.0d0))
                     sm = 0.0d0
                     do j = 1 , m
                        sm = sm + 1.0d0/(j+k)
                     enddo
                     zp = pa + pb + 2.0d0*el + sp - sm + log(1.0d0-z)
                     zr1 = zr1*(a+k-1.d0)*(b+k-1.d0)/(k*(m+k))*(1.d0-z)
                     zf1 = zf1 + zr1*zp
                     if ( abs(zf1-zw)<abs(zf1)*eps ) exit
                     zw = zf1
                  enddo
                  Zhf = zf0*zc0 + zf1*zc1
               endif
            else
               call gamma2(a,ga)
               call gamma2(b,gb)
               call gamma2(c,gc)
               call gamma2(c-a,gca)
               call gamma2(c-b,gcb)
               call gamma2(c-a-b,gcab)
               call gamma2(a+b-c,gabc)
               zc0 = gc*gcab/(gca*gcb)
               zc1 = gc*gabc/(ga*gb)*(1.0d0-z)**(c-a-b)
               Zhf = (0.0d0,0.0d0)
               zr0 = zc0
               zr1 = zc1
               do k = 1 , 500
                  zr0 = zr0*(a+k-1.d0)*(b+k-1.d0)/(k*(a+b-c+k))*(1.d0-z)
                  zr1 = zr1*(c-a+k-1.0d0)*(c-b+k-1.0d0)/(k*(c-a-b+k))   &
                      & *(1.0d0-z)
                  Zhf = Zhf + zr0 + zr1
                  if ( abs(Zhf-zw)<abs(Zhf)*eps ) exit
                  zw = Zhf
               enddo
               Zhf = Zhf + zc0 + zc1
            endif
         else
            z00 = (1.0d0,0.0d0)
            if ( c-a<a .and. c-b<b ) then
               z00 = (1.0d0-z)**(c-a-b)
               a = c - a
               b = c - b
            endif
            Zhf = (1.0d0,0.d0)
            zr = (1.0d0,0.0d0)
            do k = 1 , 1500
               zr = zr*(a+k-1.0d0)*(b+k-1.0d0)/(k*(c+k-1.0d0))*z
               Zhf = Zhf + zr
               if ( abs(Zhf-zw)<=abs(Zhf)*eps ) exit
               zw = Zhf
            enddo
            Zhf = z00*Zhf
         endif
      elseif ( a0>1.0d0 ) then
         mab = int(a-b+eps*sign(1.0d0,a-b))
         if ( abs(a-b-mab)<eps .and. a0<=1.1d0 ) b = b + eps
         if ( abs(a-b-mab)>eps ) then
            call gamma2(a,ga)
            call gamma2(b,gb)
            call gamma2(c,gc)
            call gamma2(a-b,gab)
            call gamma2(b-a,gba)
            call gamma2(c-a,gca)
            call gamma2(c-b,gcb)
            zc0 = gc*gba/(gca*gb*(-z)**a)
            zc1 = gc*gab/(gcb*ga*(-z)**b)
            zr0 = zc0
            zr1 = zc1
            Zhf = (0.0d0,0.0d0)
            do k = 1 , 500
               zr0 = zr0*(a+k-1.0d0)*(a-c+k)/((a-b+k)*k*z)
               zr1 = zr1*(b+k-1.0d0)*(b-c+k)/((b-a+k)*k*z)
               Zhf = Zhf + zr0 + zr1
               if ( abs((Zhf-zw)/Zhf)<=eps ) exit
               zw = Zhf
            enddo
            Zhf = Zhf + zc0 + zc1
         else
            if ( a-b<0.0d0 ) then
               a = bb
               b = aa
            endif
            ca = c - a
            cb = c - b
            nca = int(ca+eps*sign(1.0d0,ca))
            ncb = int(cb+eps*sign(1.0d0,cb))
            if ( abs(ca-nca)<eps .or. abs(cb-ncb)<eps ) c = c + eps
            call gamma2(a,ga)
            call gamma2(c,gc)
            call gamma2(c-b,gcb)
            call psi_spec(a,pa)
            call psi_spec(c-a,pca)
            call psi_spec(a-c,pac)
            mab = int(a-b+eps)
            zc0 = gc/(ga*(-z)**b)
            call gamma2(a-b,gm)
            zf0 = gm/gcb*zc0
            zr = zc0
            do k = 1 , mab - 1
               zr = zr*(b+k-1.0d0)/(k*z)
               t0 = a - b - k
               call gamma2(t0,g0)
               call gamma2(c-b-k,gcbk)
               zf0 = zf0 + zr*g0/gcbk
            enddo
            if ( mab==0 ) zf0 = (0.0d0,0.0d0)
            zc1 = gc/(ga*gcb*(-z)**a)
            sp = -2.0d0*el - pa - pca
            do j = 1 , mab
               sp = sp + 1.0d0/j
            enddo
            zp0 = sp + log(-z)
            sq = 1.0d0
            do j = 1 , mab
               sq = sq*(b+j-1.0d0)*(b-c+j)/j
            enddo
            zf1 = (sq*zp0)*zc1
            zr = zc1
            rk1 = 1.0d0
            sj1 = 0.0d0
            w0 = 0.0d0
            do k = 1 , 10000
               zr = zr/z
               rk1 = rk1*(b+k-1.0d0)*(b-c+k)/(k*k)
               rk2 = rk1
               do j = k + 1 , k + mab
                  rk2 = rk2*(b+j-1.0d0)*(b-c+j)/j
               enddo
               sj1 = sj1 + (a-1.0d0)/(k*(a+k-1.0d0)) + (a-c-1.0d0)      &
                   & /(k*(a-c+k-1.0d0))
               sj2 = sj1
               do j = k + 1 , k + mab
                  sj2 = sj2 + 1.0d0/j
               enddo
               zp = -2.0d0*el - pa - pac + sj2 - 1.0d0/(k+a-c)          &
                  & - pi/tan(pi*(k+a-c)) + log(-z)
               zf1 = zf1 + rk2*zr*zp
               ws = abs(zf1)
               if ( abs((ws-w0)/ws)<eps ) exit
               w0 = ws
            enddo
            Zhf = zf0 + zf1
         endif
      endif
      a = aa
      b = bb
      if ( k>150 ) Isfer = 5
      end



!       **********************************

      subroutine itairy(x,Apt,Bpt,Ant,Bnt)
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
      implicit none
      real(wp) a , Ant , Apt , Bnt , Bpt , c1 , c2 , eps , fx , &
                     & gx , pi , q0 , q1 , q2 , r , sr3 , su1 , su2 ,   &
                     & su3 , su4
      real(wp) su5 , su6 , x , xe , xp6 , xr1 , xr2
      integer k , l
      dimension a(16)
      eps = 1.0d-15
      pi = 3.141592653589793d0
      c1 = .355028053887817d0
      c2 = .258819403792807d0
      sr3 = 1.732050807568877d0
      if ( x==0.0d0 ) then
         Apt = 0.0d0
         Bpt = 0.0d0
         Ant = 0.0d0
         Bnt = 0.0d0
      elseif ( abs(x)<=9.25d0 ) then
         do l = 0 , 1
            x = (-1)**l*x
            fx = x
            r = x
            do k = 1 , 40
               r = r*(3.0*k-2.0d0)/(3.0*k+1.0d0)*x/(3.0*k)              &
                 & *x/(3.0*k-1.0d0)*x
               fx = fx + r
               if ( abs(r)<abs(fx)*eps ) exit
            enddo
            gx = .5d0*x*x
            r = gx
            do k = 1 , 40
               r = r*(3.0*k-1.0d0)/(3.0*k+2.0d0)*x/(3.0*k)              &
                 & *x/(3.0*k+1.0d0)*x
               gx = gx + r
               if ( abs(r)<abs(gx)*eps ) exit
            enddo
            Ant = c1*fx - c2*gx
            Bnt = sr3*(c1*fx+c2*gx)
            if ( l==0 ) then
               Apt = Ant
               Bpt = Bnt
            else
               Ant = -Ant
               Bnt = -Bnt
               x = -x
            endif
         enddo
      else
         data a/.569444444444444d0 , .891300154320988d0 ,               &
            & .226624344493027d+01 , .798950124766861d+01 ,             &
            & .360688546785343d+02 , .198670292131169d+03 ,             &
            & .129223456582211d+04 , .969483869669600d+04 ,             &
            & .824184704952483d+05 , .783031092490225d+06 ,             &
            & .822210493622814d+07 , .945557399360556d+08 ,             &
            & .118195595640730d+10 , .159564653040121d+11 ,             &
            & .231369166433050d+12 , .358622522796969d+13/
         q2 = 1.414213562373095d0
         q0 = .3333333333333333d0
         q1 = .6666666666666667d0
         xe = x*sqrt(x)/1.5d0
         xp6 = 1.0d0/sqrt(6.0d0*pi*xe)
         su1 = 1.0d0
         r = 1.0d0
         xr1 = 1.0d0/xe
         do k = 1 , 16
            r = -r*xr1
            su1 = su1 + a(k)*r
         enddo
         su2 = 1.0d0
         r = 1.0d0
         do k = 1 , 16
            r = r*xr1
            su2 = su2 + a(k)*r
         enddo
         Apt = q0 - exp(-xe)*xp6*su1
         Bpt = 2.0d0*exp(xe)*xp6*su2
         su3 = 1.0d0
         r = 1.0d0
         xr2 = 1.0d0/(xe*xe)
         do k = 1 , 8
            r = -r*xr2
            su3 = su3 + a(2*k)*r
         enddo
         su4 = a(1)*xr1
         r = xr1
         do k = 1 , 7
            r = -r*xr2
            su4 = su4 + a(2*k+1)*r
         enddo
         su5 = su3 + su4
         su6 = su3 - su4
         Ant = q1 - q2*xp6*(su5*cos(xe)-su6*sin(xe))
         Bnt = q2*xp6*(su5*sin(xe)+su6*cos(xe))
      endif
      end

!       **********************************

      subroutine ikna(n,x,Nm,Bi,Di,Bk,Dk)
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
      implicit none
      real(wp) Bi , bi0 , bi1 , Bk , bk0 , bk1 , Di , di0 ,     &
                     & di1 , Dk , dk0 , dk1 , f , f0 , f1 , g , g0 ,    &
                     & g1 , h , h0
      real(wp) h1 , s0 , x
      integer k , m , n , Nm
      dimension Bi(0:n) , Di(0:n) , Bk(0:n) , Dk(0:n)
      Nm = n
      if ( x<=1.0d-100 ) then
         do k = 0 , n
            Bi(k) = 0.0d0
            Di(k) = 0.0d0
            Bk(k) = 1.0d+300
            Dk(k) = -1.0d+300
         enddo
         Bi(0) = 1.0d0
         Di(1) = 0.5d0
         return
      endif
      call ik01a(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1)
      Bi(0) = bi0
      Bi(1) = bi1
      Bk(0) = bk0
      Bk(1) = bk1
      Di(0) = di0
      Di(1) = di1
      Dk(0) = dk0
      Dk(1) = dk1
      if ( n<=1 ) return
      if ( x>40.0 .and. n<int(0.25*x) ) then
         h0 = bi0
         h1 = bi1
         do k = 2 , n
            h = -2.0d0*(k-1.0d0)/x*h1 + h0
            Bi(k) = h
            h0 = h1
            h1 = h
         enddo
      else
         m = msta1(x,200)
         if ( m<n ) then
            Nm = m
         else
            m = msta2(x,n,15)
         endif
         f0 = 0.0d0
         f1 = 1.0d-100
         f = 0.0d0
         do k = m , 0 , -1
            f = 2.0d0*(k+1.0d0)*f1/x + f0
            if ( k<=Nm ) Bi(k) = f
            f0 = f1
            f1 = f
         enddo
         s0 = bi0/f
         do k = 0 , Nm
            Bi(k) = s0*Bi(k)
         enddo
      endif
      g0 = bk0
      g1 = bk1
      do k = 2 , Nm
         g = 2.0d0*(k-1.0d0)/x*g1 + g0
         Bk(k) = g
         g0 = g1
         g1 = g
      enddo
      do k = 2 , Nm
         Di(k) = Bi(k-1) - k/x*Bi(k)
         Dk(k) = -Bk(k-1) - k/x*Bk(k)
      enddo
      end



!       **********************************

      subroutine cjynb(n,z,Nm,Cbj,Cdj,Cby,Cdy)
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
      implicit none
      real(wp) a , a0 , a1 , b , b1 , el , pi , r2p , y0
      complex(wp) Cbj , cbj0 , cbj1 , cbjk , cbs , Cby , cby0 , cby1 ,   &
               & Cdj , Cdy , ce , cf , cf1 , cf2 , cp0 , cp1 , cq0 ,    &
               & cq1 , cs0 , csu
      complex(wp) csv , ct1 , ct2 , cu , cyy , z
      integer k , m , n , Nm
      dimension Cbj(0:n) , Cdj(0:n) , Cby(0:n) , Cdy(0:n) , a(4) ,      &
              & b(4) , a1(4) , b1(4)
      el = 0.5772156649015329d0
      pi = 3.141592653589793d0
      r2p = .63661977236758d0
      y0 = abs(dimag(z))
      a0 = abs(z)
      Nm = n
      if ( a0<1.0d-100 ) then
         do k = 0 , n
            Cbj(k) = (0.0d0,0.0d0)
            Cdj(k) = (0.0d0,0.0d0)
            Cby(k) = -(1.0d+300,0.0d0)
            Cdy(k) = (1.0d+300,0.0d0)
         enddo
         Cbj(0) = (1.0d0,0.0d0)
         Cdj(1) = (0.5d0,0.0d0)
         return
      endif
      if ( a0<=300.d0 .or. n>80 ) then
         if ( n==0 ) Nm = 1
         m = msta1(a0,200)
         if ( m<Nm ) then
            Nm = m
         else
            m = msta2(a0,Nm,15)
         endif
         cbs = (0.0d0,0.0d0)
         csu = (0.0d0,0.0d0)
         csv = (0.0d0,0.0d0)
         cf2 = (0.0d0,0.0d0)
         cf1 = (1.0d-100,0.0d0)
         do k = m , 0 , -1
            cf = 2.0d0*(k+1.0d0)/z*cf1 - cf2
            if ( k<=Nm ) Cbj(k) = cf
            if ( k==2*int(k/2) .and. k/=0 ) then
               if ( y0<=1.0d0 ) then
                  cbs = cbs + 2.0d0*cf
               else
                  cbs = cbs + (-1)**(k/2)*2.0d0*cf
               endif
               csu = csu + (-1)**(k/2)*cf/k
            elseif ( k>1 ) then
               csv = csv + (-1)**(k/2)*k/(k*k-1.0d0)*cf
            endif
            cf2 = cf1
            cf1 = cf
         enddo
         if ( y0<=1.0d0 ) then
            cs0 = cbs + cf
         else
            cs0 = (cbs+cf)/cos(z)
         endif
         do k = 0 , Nm
            Cbj(k) = Cbj(k)/cs0
         enddo
         ce = log(z/2.0d0) + el
         Cby(0) = r2p*(ce*Cbj(0)-4.0d0*csu/cs0)
         Cby(1) = r2p*(-Cbj(0)/z+(ce-1.0d0)*Cbj(1)-4.0d0*csv/cs0)
      else
         data a/ - .7031250000000000d-01 , .1121520996093750d+00 ,      &
            & -.5725014209747314d+00 , .6074042001273483d+01/
         data b/.7324218750000000d-01 , -.2271080017089844d+00 ,        &
            & .1727727502584457d+01 , -.2438052969955606d+02/
         data a1/.1171875000000000d+00 , -.1441955566406250d+00 ,       &
            & .6765925884246826d+00 , -.6883914268109947d+01/
         data b1/ - .1025390625000000d+00 , .2775764465332031d+00 ,     &
            & -.1993531733751297d+01 , .2724882731126854d+02/
         ct1 = z - 0.25d0*pi
         cp0 = (1.0d0,0.0d0)
         do k = 1 , 4
            cp0 = cp0 + a(k)*z**(-2*k)
         enddo
         cq0 = -0.125d0/z
         do k = 1 , 4
            cq0 = cq0 + b(k)*z**(-2*k-1)
         enddo
         cu = sqrt(r2p/z)
         cbj0 = cu*(cp0*cos(ct1)-cq0*sin(ct1))
         cby0 = cu*(cp0*sin(ct1)+cq0*cos(ct1))
         Cbj(0) = cbj0
         Cby(0) = cby0
         ct2 = z - 0.75d0*pi
         cp1 = (1.0d0,0.0d0)
         do k = 1 , 4
            cp1 = cp1 + a1(k)*z**(-2*k)
         enddo
         cq1 = 0.375d0/z
         do k = 1 , 4
            cq1 = cq1 + b1(k)*z**(-2*k-1)
         enddo
         cbj1 = cu*(cp1*cos(ct2)-cq1*sin(ct2))
         cby1 = cu*(cp1*sin(ct2)+cq1*cos(ct2))
         Cbj(1) = cbj1
         Cby(1) = cby1
         do k = 2 , Nm
            cbjk = 2.0d0*(k-1.0d0)/z*cbj1 - cbj0
            Cbj(k) = cbjk
            cbj0 = cbj1
            cbj1 = cbjk
         enddo
      endif
      Cdj(0) = -Cbj(1)
      do k = 1 , Nm
         Cdj(k) = Cbj(k-1) - k/z*Cbj(k)
      enddo
      if ( abs(Cbj(0))>1.0d0 ) Cby(1) = (Cbj(1)*Cby(0)-2.0d0/(pi*z))    &
                                      & /Cbj(0)
      do k = 2 , Nm
         if ( abs(Cbj(k-1))>=abs(Cbj(k-2)) ) then
            cyy = (Cbj(k)*Cby(k-1)-2.0d0/(pi*z))/Cbj(k-1)
         else
            cyy = (Cbj(k)*Cby(k-2)-4.0d0*(k-1.0d0)/(pi*z*z))/Cbj(k-2)
         endif
         Cby(k) = cyy
      enddo
      Cdy(0) = -Cby(1)
      do k = 1 , Nm
         Cdy(k) = Cby(k-1) - k/z*Cby(k)
      enddo
      end



!       **********************************

      subroutine iknb(n,x,Nm,Bi,Di,Bk,Dk)
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
      implicit none
      real(wp) a0 , Bi , Bk , bkl , bs , Di , Dk , el , f , f0 ,&
                     & f1 , g , g0 , g1 , pi , r , s0 , sk0 , vt , x
      integer k , k0 , l , m , n , Nm
      dimension Bi(0:n) , Di(0:n) , Bk(0:n) , Dk(0:n)
      pi = 3.141592653589793d0
      el = 0.5772156649015329d0
      Nm = n
      if ( x<=1.0d-100 ) then
         do k = 0 , n
            Bi(k) = 0.0d0
            Di(k) = 0.0d0
            Bk(k) = 1.0d+300
            Dk(k) = -1.0d+300
         enddo
         Bi(0) = 1.0d0
         Di(1) = 0.5d0
         return
      endif
      if ( n==0 ) Nm = 1
      m = msta1(x,200)
      if ( m<Nm ) then
         Nm = m
      else
         m = msta2(x,Nm,15)
      endif
      bs = 0.0d0
      sk0 = 0.0d0
      f = 0.0d0
      f0 = 0.0d0
      f1 = 1.0d-100
      do k = m , 0 , -1
         f = 2.0d0*(k+1.0d0)/x*f1 + f0
         if ( k<=Nm ) Bi(k) = f
         if ( k/=0 .and. k==2*int(k/2) ) sk0 = sk0 + 4.0d0*f/k
         bs = bs + 2.0d0*f
         f0 = f1
         f1 = f
      enddo
      s0 = exp(x)/(bs-f)
      do k = 0 , Nm
         Bi(k) = s0*Bi(k)
      enddo
      if ( x<=8.0d0 ) then
         Bk(0) = -(log(0.5d0*x)+el)*Bi(0) + s0*sk0
         Bk(1) = (1.0d0/x-Bi(1)*Bk(0))/Bi(0)
      else
         a0 = sqrt(pi/(2.0d0*x))*exp(-x)
         k0 = 16
         if ( x>=25.0 ) k0 = 10
         if ( x>=80.0 ) k0 = 8
         if ( x>=200.0 ) k0 = 6
         do l = 0 , 1
            bkl = 1.0d0
            vt = 4.0d0*l
            r = 1.0d0
            do k = 1 , k0
               r = 0.125d0*r*(vt-(2.0*k-1.0)**2)/(k*x)
               bkl = bkl + r
            enddo
            Bk(l) = a0*bkl
         enddo
      endif
      g0 = Bk(0)
      g1 = Bk(1)
      do k = 2 , Nm
         g = 2.0d0*(k-1.0d0)/x*g1 + g0
         Bk(k) = g
         g0 = g1
         g1 = g
      enddo
      Di(0) = Bi(1)
      Dk(0) = -Bk(1)
      do k = 1 , Nm
         Di(k) = Bi(k-1) - k/x*Bi(k)
         Dk(k) = -Bk(k-1) - k/x*Bk(k)
      enddo
      end



!       **********************************

      subroutine lpmn(Mm,m,n,x,Pm,Pd)
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
      implicit none
      real(wp) Pd , Pm , x , xq , xs
      integer i , j , ls , m , Mm , n
      dimension Pm(0:Mm,0:n) , Pd(0:Mm,0:n)
      intrinsic min
      do i = 0 , n
         do j = 0 , m
            Pm(j,i) = 0.0d0
            Pd(j,i) = 0.0d0
         enddo
      enddo
      Pm(0,0) = 1.0d0
      if ( n==0 ) return
      if ( abs(x)==1.0d0 ) then
         do i = 1 , n
            Pm(0,i) = x**i
            Pd(0,i) = 0.5d0*i*(i+1.0d0)*x**(i+1)
         enddo
         do j = 1 , n
            do i = 1 , m
               if ( i==1 ) then
                  Pd(i,j) = dinf()
               elseif ( i==2 ) then
                  Pd(i,j) = -0.25d0*(j+2)*(j+1)*j*(j-1)*x**(j+1)
               endif
            enddo
         enddo
         return
      endif
      ls = 1
      if ( abs(x)>1.0d0 ) ls = -1
      xq = sqrt(ls*(1.0d0-x*x))
!       Ensure connection to the complex-valued function for |x| > 1
      if ( x<-1d0 ) xq = -xq
      xs = ls*(1.0d0-x*x)
      do i = 1 , m
         Pm(i,i) = -ls*(2.0d0*i-1.0d0)*xq*Pm(i-1,i-1)
      enddo
      do i = 0 , min(m,n-1)
         Pm(i,i+1) = (2.0d0*i+1.0d0)*x*Pm(i,i)
      enddo
      do i = 0 , m
         do j = i + 2 , n
            Pm(i,j) = ((2.0d0*j-1.0d0)*x*Pm(i,j-1)-(i+j-1.0d0)*Pm(i,j-2)&
                    & )/(j-i)
         enddo
      enddo
      Pd(0,0) = 0.0d0
      do j = 1 , n
         Pd(0,j) = ls*j*(Pm(0,j-1)-x*Pm(0,j))/xs
      enddo
      do i = 1 , m
         do j = i , n
            Pd(i,j) = ls*i*x*Pm(i,j)/xs + (j+i)*(j-i+1.0d0)/xq*Pm(i-1,j)
         enddo
      enddo
      end

!       **********************************

      subroutine mtu0(Kf,m,q,x,Csf,Csd)
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
      implicit none
      real(wp) a , Csd , Csf , eps , fg , q , qm , rd ,  &
                     & x , xr
      integer ic , k , kd , Kf , km , m
      dimension fg(251)
      eps = 1.0d-14
      if ( Kf==1 .and. m==2*int(m/2) ) kd = 1
      if ( Kf==1 .and. m/=2*int(m/2) ) kd = 2
      if ( Kf==2 .and. m/=2*int(m/2) ) kd = 3
      if ( Kf==2 .and. m==2*int(m/2) ) kd = 4
      call cva2(kd,m,q,a)
      if ( q<=1.0d0 ) then
         qm = 7.5 + 56.1*sqrt(q) - 134.7*q + 90.7*sqrt(q)*q
      else
         qm = 17.0 + 3.1*sqrt(q) - .126*q + .0037*sqrt(q)*q
      endif
      km = int(qm+0.5*m)
      if ( km>251 ) then
         Csf = dnan()
         Csd = dnan()
         return
      endif
      call fcoef(kd,m,q,a,fg)
      ic = int(m/2) + 1
      rd = 1.74532925199433d-2
      xr = x*rd
      Csf = 0.0d0
      do k = 1 , km
         if ( kd==1 ) then
            Csf = Csf + fg(k)*cos((2*k-2)*xr)
         elseif ( kd==2 ) then
            Csf = Csf + fg(k)*cos((2*k-1)*xr)
         elseif ( kd==3 ) then
            Csf = Csf + fg(k)*sin((2*k-1)*xr)
         elseif ( kd==4 ) then
            Csf = Csf + fg(k)*sin(2*k*xr)
         endif
         if ( k>=ic .and. abs(fg(k))<abs(Csf)*eps ) exit
      enddo
      Csd = 0.0d0
      do k = 1 , km
         if ( kd==1 ) then
            Csd = Csd - (2*k-2)*fg(k)*sin((2*k-2)*xr)
         elseif ( kd==2 ) then
            Csd = Csd - (2*k-1)*fg(k)*sin((2*k-1)*xr)
         elseif ( kd==3 ) then
            Csd = Csd + (2*k-1)*fg(k)*cos((2*k-1)*xr)
         elseif ( kd==4 ) then
            Csd = Csd + 2.0d0*k*fg(k)*cos(2*k*xr)
         endif
         if ( k>=ic .and. abs(fg(k))<abs(Csd)*eps ) exit
      enddo
      end



!       **********************************

      subroutine cy01(Kf,z,Zf,Zd)
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
      implicit none
      real(wp) a , a0 , a1 , b , b1 , el , pi , rp2 , w0 , w1
      complex(wp) cbj0 , cbj1 , cby0 , cby1 , cdy0 , cdy1 , ci , cp ,    &
               & cp0 , cp1 , cq0 , cq1 , cr , cs , ct1 , ct2 , cu , z , &
               & z1 , z2
      complex(wp) Zd , Zf
      integer k , k0 , Kf
      dimension a(12) , b(12) , a1(12) , b1(12)
      pi = 3.141592653589793d0
      el = 0.5772156649015329d0
      rp2 = 2.0d0/pi
      ci = (0.0d0,1.0d0)
      a0 = abs(z)
      z2 = z*z
      z1 = z
      if ( a0==0.0d0 ) then
         cbj0 = (1.0d0,0.0d0)
         cbj1 = (0.0d0,0.0d0)
         cby0 = -(1.0d300,0.0d0)
         cby1 = -(1.0d300,0.0d0)
         cdy0 = (1.0d300,0.0d0)
         cdy1 = (1.0d300,0.0d0)
         goto 300
      endif
      if ( dble(z)<0.0 ) z1 = -z
      if ( a0<=12.0 ) then
         cbj0 = (1.0d0,0.0d0)
         cr = (1.0d0,0.0d0)
         do k = 1 , 40
            cr = -0.25d0*cr*z2/(k*k)
            cbj0 = cbj0 + cr
            if ( abs(cr)<abs(cbj0)*1.0d-15 ) exit
         enddo
         cbj1 = (1.0d0,0.0d0)
         cr = (1.0d0,0.0d0)
         do k = 1 , 40
            cr = -0.25d0*cr*z2/(k*(k+1.0d0))
            cbj1 = cbj1 + cr
            if ( abs(cr)<abs(cbj1)*1.0d-15 ) exit
         enddo
         cbj1 = 0.5d0*z1*cbj1
         w0 = 0.0d0
         cr = (1.0d0,0.0d0)
         cs = (0.0d0,0.0d0)
         do k = 1 , 40
            w0 = w0 + 1.0d0/k
            cr = -0.25d0*cr/(k*k)*z2
            cp = cr*w0
            cs = cs + cp
            if ( abs(cp)<abs(cs)*1.0d-15 ) exit
         enddo
         cby0 = rp2*(log(z1/2.0d0)+el)*cbj0 - rp2*cs
         w1 = 0.0d0
         cr = (1.0d0,0.0d0)
         cs = (1.0d0,0.0d0)
         do k = 1 , 40
            w1 = w1 + 1.0d0/k
            cr = -0.25d0*cr/(k*(k+1))*z2
            cp = cr*(2.0d0*w1+1.0d0/(k+1.0d0))
            cs = cs + cp
            if ( abs(cp)<abs(cs)*1.0d-15 ) exit
         enddo
         cby1 = rp2*((log(z1/2.0d0)+el)*cbj1-1.0d0/z1-.25d0*z1*cs)
      else
         data a/ - .703125d-01 , .112152099609375d+00 ,                 &
            & -.5725014209747314d+00 , .6074042001273483d+01 ,          &
            & -.1100171402692467d+03 , .3038090510922384d+04 ,          &
            & -.1188384262567832d+06 , .6252951493434797d+07 ,          &
            & -.4259392165047669d+09 , .3646840080706556d+11 ,          &
            & -.3833534661393944d+13 , .4854014686852901d+15/
         data b/.732421875d-01 , -.2271080017089844d+00 ,               &
            & .1727727502584457d+01 , -.2438052969955606d+02 ,          &
            & .5513358961220206d+03 , -.1825775547429318d+05 ,          &
            & .8328593040162893d+06 , -.5006958953198893d+08 ,          &
            & .3836255180230433d+10 , -.3649010818849833d+12 ,          &
            & .4218971570284096d+14 , -.5827244631566907d+16/
         data a1/.1171875d+00 , -.144195556640625d+00 ,                 &
            & .6765925884246826d+00 , -.6883914268109947d+01 ,          &
            & .1215978918765359d+03 , -.3302272294480852d+04 ,          &
            & .1276412726461746d+06 , -.6656367718817688d+07 ,          &
            & .4502786003050393d+09 , -.3833857520742790d+11 ,          &
            & .4011838599133198d+13 , -.5060568503314727d+15/
         data b1/ - .1025390625d+00 , .2775764465332031d+00 ,           &
            & -.1993531733751297d+01 , .2724882731126854d+02 ,          &
            & -.6038440767050702d+03 , .1971837591223663d+05 ,          &
            & -.8902978767070678d+06 , .5310411010968522d+08 ,          &
            & -.4043620325107754d+10 , .3827011346598605d+12 ,          &
            & -.4406481417852278d+14 , .6065091351222699d+16/
         k0 = 12
         if ( a0>=35.0 ) k0 = 10
         if ( a0>=50.0 ) k0 = 8
         ct1 = z1 - .25d0*pi
         cp0 = (1.0d0,0.0d0)
         do k = 1 , k0
            cp0 = cp0 + a(k)*z1**(-2*k)
         enddo
         cq0 = -0.125d0/z1
         do k = 1 , k0
            cq0 = cq0 + b(k)*z1**(-2*k-1)
         enddo
         cu = sqrt(rp2/z1)
         cbj0 = cu*(cp0*cos(ct1)-cq0*sin(ct1))
         cby0 = cu*(cp0*sin(ct1)+cq0*cos(ct1))
         ct2 = z1 - .75d0*pi
         cp1 = (1.0d0,0.0d0)
         do k = 1 , k0
            cp1 = cp1 + a1(k)*z1**(-2*k)
         enddo
         cq1 = 0.375d0/z1
         do k = 1 , k0
            cq1 = cq1 + b1(k)*z1**(-2*k-1)
         enddo
         cbj1 = cu*(cp1*cos(ct2)-cq1*sin(ct2))
         cby1 = cu*(cp1*sin(ct2)+cq1*cos(ct2))
      endif
      if ( dble(z)<0.0 ) then
         if ( dimag(z)<0.0 ) cby0 = cby0 - 2.0d0*ci*cbj0
         if ( dimag(z)>0.0 ) cby0 = cby0 + 2.0d0*ci*cbj0
         if ( dimag(z)<0.0 ) cby1 = -(cby1-2.0d0*ci*cbj1)
         if ( dimag(z)>0.0 ) cby1 = -(cby1+2.0d0*ci*cbj1)
         cbj1 = -cbj1
      endif
      cdy0 = -cby1
      cdy1 = cby0 - 1.0d0/z*cby1
 300  if ( Kf==0 ) then
         Zf = cby0
         Zd = cdy0
      elseif ( Kf==1 ) then
         Zf = cby1
         Zd = cdy1
      elseif ( Kf==2 ) then
         Zf = cdy1
         Zd = -cdy1/z - (1.0d0-1.0d0/(z*z))*cby1
      endif
      end


!       **********************************

      subroutine ffk(Ks,x,Fr,Fi,Fm,Fa,Gr,Gi,Gm,Ga)
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
      implicit none
      real(wp) c1 , cs , eps , Fa , Fi , fi0 , Fm , Fr , Ga ,   &
                     & Gi , Gm , Gr , p2p , pi , pp2 , s1 , srd , ss ,  &
                     & x , x2
      real(wp) x4 , xa , xc , xf , xf0 , xf1 , xg , xp , xq ,   &
                     & xq2 , xr , xs , xsu , xw
      integer k , Ks , m
      srd = 57.29577951308233d0
      eps = 1.0d-15
      pi = 3.141592653589793d0
      pp2 = 1.2533141373155d0
      p2p = .7978845608028654d0
      xa = abs(x)
      x2 = x*x
      x4 = x2*x2
      if ( x==0.0d0 ) then
         Fr = .5d0*sqrt(0.5d0*pi)
         Fi = (-1)**Ks*Fr
         Fm = sqrt(0.25d0*pi)
         Fa = (-1)**Ks*45.0d0
         Gr = .5d0
         Gi = 0.0d0
         Gm = .5d0
         Ga = 0.0d0
      else
         if ( xa<=2.5d0 ) then
            xr = p2p*xa
            c1 = xr
            do k = 1 , 50
               xr = -.5d0*xr*(4.0d0*k-3.0d0)/k/(2.0d0*k-1.0d0)          &
                  & /(4.0d0*k+1.0d0)*x4
               c1 = c1 + xr
               if ( abs(xr/c1)<eps ) exit
            enddo
            s1 = p2p*xa*xa*xa/3.0d0
            xr = s1
            do k = 1 , 50
               xr = -.5d0*xr*(4.0d0*k-1.0d0)/k/(2.0d0*k+1.0d0)          &
                  & /(4.0d0*k+3.0d0)*x4
               s1 = s1 + xr
               if ( abs(xr/s1)<eps ) goto 50
            enddo
         elseif ( xa<5.5d0 ) then
            m = int(42+1.75*x2)
            xsu = 0.0d0
            xc = 0.0d0
            xs = 0.0d0
            xf1 = 0.0d0
            xf0 = 1d-100
            do k = m , 0 , -1
               xf = (2.0d0*k+3.0d0)*xf0/x2 - xf1
               if ( k==2*int(k/2) ) then
                  xc = xc + xf
               else
                  xs = xs + xf
               endif
               xsu = xsu + (2.0d0*k+1.0d0)*xf*xf
               xf1 = xf0
               xf0 = xf
            enddo
            xq = sqrt(xsu)
            xw = p2p*xa/xq
            c1 = xc*xw
            s1 = xs*xw
         else
            xr = 1.0d0
            xf = 1.0d0
            do k = 1 , 12
               xr = -.25d0*xr*(4.0d0*k-1.0d0)*(4.0d0*k-3.0d0)/x4
               xf = xf + xr
            enddo
            xr = 1.0d0/(2.0d0*xa*xa)
            xg = xr
            do k = 1 , 12
               xr = -.25d0*xr*(4.0d0*k+1.0d0)*(4.0d0*k-1.0d0)/x4
               xg = xg + xr
            enddo
            c1 = .5d0 + (xf*sin(x2)-xg*cos(x2))/sqrt(2.0d0*pi)/xa
            s1 = .5d0 - (xf*cos(x2)+xg*sin(x2))/sqrt(2.0d0*pi)/xa
         endif
 50      Fr = pp2*(.5d0-c1)
         fi0 = pp2*(.5d0-s1)
         Fi = (-1)**Ks*fi0
         Fm = sqrt(Fr*Fr+Fi*Fi)
         if ( Fr>=0.0 ) then
            Fa = srd*atan(Fi/Fr)
         elseif ( Fi>0.0 ) then
            Fa = srd*(atan(Fi/Fr)+pi)
         elseif ( Fi<0.0 ) then
            Fa = srd*(atan(Fi/Fr)-pi)
         endif
         xp = x*x + pi/4.0d0
         cs = cos(xp)
         ss = sin(xp)
         xq2 = 1.0d0/sqrt(pi)
         Gr = xq2*(Fr*cs+fi0*ss)
         Gi = (-1)**Ks*xq2*(fi0*cs-Fr*ss)
         Gm = sqrt(Gr*Gr+Gi*Gi)
         if ( Gr>=0.0 ) then
            Ga = srd*atan(Gi/Gr)
         elseif ( Gi>0.0 ) then
            Ga = srd*(atan(Gi/Gr)+pi)
         elseif ( Gi<0.0 ) then
            Ga = srd*(atan(Gi/Gr)-pi)
         endif
         if ( x<0.0d0 ) then
            Fr = pp2 - Fr
            Fi = (-1)**Ks*pp2 - Fi
            Fm = sqrt(Fr*Fr+Fi*Fi)
            Fa = srd*atan(Fi/Fr)
            Gr = cos(x*x) - Gr
            Gi = -(-1)**Ks*sin(x*x) - Gi
            Gm = sqrt(Gr*Gr+Gi*Gi)
            Ga = srd*atan(Gi/Gr)
         endif
      endif
      end

!       **********************************

      subroutine airya(x,Ai,Bi,Ad,Bd)
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
      implicit none
      real(wp) Ad , Ai , Bd , Bi , c1 , c2 , pir , sr3 , vi1 ,  &
                     & vi2 , vj1 , vj2 , vk1 , vk2 , vy1 , vy2 , x ,    &
                     & xa , xq , z
      xa = abs(x)
      pir = 0.318309886183891d0
      c1 = 0.355028053887817d0
      c2 = 0.258819403792807d0
      sr3 = 1.732050807568877d0
      z = xa**1.5/1.5d0
      xq = sqrt(xa)
      call ajyik(z,vj1,vj2,vy1,vy2,vi1,vi2,vk1,vk2)
      if ( x==0.0d0 ) then
         Ai = c1
         Bi = sr3*c1
         Ad = -c2
         Bd = sr3*c2
      elseif ( x>0.0d0 ) then
         Ai = pir*xq/sr3*vk1
         Bi = xq*(pir*vk1+2.0d0/sr3*vi1)
         Ad = -xa/sr3*pir*vk2
         Bd = xa*(pir*vk2+2.0d0/sr3*vi2)
      else
         Ai = 0.5d0*xq*(vj1-vy1/sr3)
         Bi = -0.5d0*xq*(vj1/sr3+vy1)
         Ad = 0.5d0*xa*(vj2+vy2/sr3)
         Bd = 0.5d0*xa*(vj2/sr3-vy2)
      endif
      end



!       **********************************

      subroutine airyb(x,Ai,Bi,Ad,Bd)
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
      implicit none
      real(wp) Ad , Ai , Bd , Bi , c1 , c2 , ck , df , dg , dk ,&
                     & eps , fx , gx , pi , r , rp , sad , sai , sbd ,  &
                     & sbi
      real(wp) sda , sdb , sr3 , ssa , ssb , x , xa , xar ,     &
                     & xcs , xe , xf , xm , xp1 , xq , xr1 , xr2 , xss
      integer k , km , km2 , kmax
      dimension ck(51) , dk(51)
      eps = 1.0d-15
      pi = 3.141592653589793d0
      c1 = 0.355028053887817d0
      c2 = 0.258819403792807d0
      sr3 = 1.732050807568877d0
      xa = abs(x)
      xq = sqrt(xa)
      xm = 8.0d0
      if ( x>0.0d0 ) xm = 5.0d0
      if ( x==0.0d0 ) then
         Ai = c1
         Bi = sr3*c1
         Ad = -c2
         Bd = sr3*c2
         return
      endif
      if ( xa<=xm ) then
         fx = 1.0d0
         r = 1.0d0
         do k = 1 , 40
            r = r*x/(3.0d0*k)*x/(3.0d0*k-1.0d0)*x
            fx = fx + r
            if ( abs(r)<abs(fx)*eps ) exit
         enddo
         gx = x
         r = x
         do k = 1 , 40
            r = r*x/(3.0d0*k)*x/(3.0d0*k+1.0d0)*x
            gx = gx + r
            if ( abs(r)<abs(gx)*eps ) exit
         enddo
         Ai = c1*fx - c2*gx
         Bi = sr3*(c1*fx+c2*gx)
         df = 0.5d0*x*x
         r = df
         do k = 1 , 40
            r = r*x/(3.0d0*k)*x/(3.0d0*k+2.0d0)*x
            df = df + r
            if ( abs(r)<abs(df)*eps ) exit
         enddo
         dg = 1.0d0
         r = 1.0d0
         do k = 1 , 40
            r = r*x/(3.0d0*k)*x/(3.0d0*k-2.0d0)*x
            dg = dg + r
            if ( abs(r)<abs(dg)*eps ) exit
         enddo
         Ad = c1*df - c2*dg
         Bd = sr3*(c1*df+c2*dg)
      else
         km = int(24.5-xa)
         if ( xa<6.0 ) km = 14
         if ( xa>15.0 ) km = 10
         if ( x>0.0d0 ) then
            kmax = km
         else
!             Choose cutoffs so that the remainder term in asymptotic
!             expansion is epsilon size. The X<0 branch needs to be fast
!             in order to make AIRYZO efficient
            if ( xa>70.0 ) km = 3
            if ( xa>500.0 ) km = 2
            if ( xa>1000.0 ) km = 1
            km2 = km
            if ( xa>150.0 ) km2 = 1
            if ( xa>3000.0 ) km2 = 0
            kmax = 2*km + 1
         endif
         xe = xa*xq/1.5d0
         xr1 = 1.0d0/xe
         xar = 1.0d0/xq
         xf = sqrt(xar)
         rp = 0.5641895835477563d0
         r = 1.0d0
         do k = 1 , kmax
            r = r*(6.0d0*k-1.0d0)/216.0d0*(6.0d0*k-3.0d0)               &
              & /k*(6.0d0*k-5.0d0)/(2.0d0*k-1.0d0)
            ck(k) = r
            dk(k) = -(6.0d0*k+1.0d0)/(6.0d0*k-1.0d0)*ck(k)
         enddo
         if ( x>0.0d0 ) then
            sai = 1.0d0
            sad = 1.0d0
            r = 1.0d0
            do k = 1 , km
               r = -r*xr1
               sai = sai + ck(k)*r
               sad = sad + dk(k)*r
            enddo
            sbi = 1.0d0
            sbd = 1.0d0
            r = 1.0d0
            do k = 1 , km
               r = r*xr1
               sbi = sbi + ck(k)*r
               sbd = sbd + dk(k)*r
            enddo
            xp1 = exp(-xe)
            Ai = 0.5d0*rp*xf*xp1*sai
            Bi = rp*xf/xp1*sbi
            Ad = -.5d0*rp/xf*xp1*sad
            Bd = rp/xf/xp1*sbd
         else
            xcs = cos(xe+pi/4.0d0)
            xss = sin(xe+pi/4.0d0)
            ssa = 1.0d0
            sda = 1.0d0
            r = 1.0d0
            xr2 = 1.0d0/(xe*xe)
            do k = 1 , km
               r = -r*xr2
               ssa = ssa + ck(2*k)*r
               sda = sda + dk(2*k)*r
            enddo
            ssb = ck(1)*xr1
            sdb = dk(1)*xr1
            r = xr1
            do k = 1 , km2
               r = -r*xr2
               ssb = ssb + ck(2*k+1)*r
               sdb = sdb + dk(2*k+1)*r
            enddo
            Ai = rp*xf*(xss*ssa-xcs*ssb)
            Bi = rp*xf*(xcs*ssa+xss*ssb)
            Ad = -rp/xf*(xcs*sda+xss*sdb)
            Bd = rp/xf*(xss*sda-xcs*sdb)
         endif
      endif
      end

!       **********************************

      subroutine scka(m,n,c,Cv,Kd,Ck)
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
      implicit none
      real(wp) c , Ck , cs , Cv , f , f0 , f1 , f2 , fl , fs ,  &
                     & r1 , r2 , s0 , su1 , su2
      integer ip , j , k , k1 , kb , Kd , m , n , nm
      dimension Ck(200)
      if ( c<=1.0d-10 ) c = 1.0d-10
      nm = 25 + int((n-m)/2+c)
      cs = c*c*Kd
      ip = 1
      if ( n-m==2*int((n-m)/2) ) ip = 0
      fs = 1.0d0
      f1 = 0.0d0
      f0 = 1.0d-100
      kb = 0
      Ck(nm+1) = 0.0d0
      fl = 0.0d0
      do k = nm , 1 , -1
         f = (((2.0d0*k+m+ip)*(2.0d0*k+m+1.0d0+ip)-Cv+cs)               &
           & *f0-4.0d0*(k+1.0d0)*(k+m+1.0d0)*f1)/cs
         if ( abs(f)>abs(Ck(k+1)) ) then
            Ck(k) = f
            f1 = f0
            f0 = f
            if ( abs(f)>1.0d+100 ) then
               do k1 = nm , k , -1
                  Ck(k1) = Ck(k1)*1.0d-100
               enddo
               f1 = f1*1.0d-100
               f0 = f0*1.0d-100
            endif
         else
            kb = k
            fl = Ck(k+1)
            f1 = 1.0d0
            f2 = 0.25d0*((m+ip)*(m+ip+1.0)-Cv+cs)/(m+1.0)*f1
            Ck(1) = f1
            if ( kb==1 ) then
               fs = f2
            elseif ( kb==2 ) then
               Ck(2) = f2
               fs = 0.125d0*(((m+ip+2.0)*(m+ip+3.0)-Cv+cs)*f2-cs*f1)    &
                  & /(m+2.0)
            else
               Ck(2) = f2
               do j = 3 , kb + 1
                  f = 0.25d0*(((2.0*j+m+ip-4.0)*(2.0*j+m+ip-3.0)-Cv+cs) &
                    & *f2-cs*f1)/((j-1.0)*(j+m-1.0))
                  if ( j<=kb ) Ck(j) = f
                  f1 = f2
                  f2 = f
               enddo
               fs = f
            endif
            goto 100
         endif
      enddo
 100  su1 = 0.0d0
      do k = 1 , kb
         su1 = su1 + Ck(k)
      enddo
      su2 = 0.0d0
      do k = kb + 1 , nm
         su2 = su2 + Ck(k)
      enddo
      r1 = 1.0d0
      do j = 1 , (n+m+ip)/2
         r1 = r1*(j+0.5d0*(n+m+ip))
      enddo
      r2 = 1.0d0
      do j = 1 , (n-m-ip)/2
         r2 = -r2*j
      enddo
      if ( kb==0 ) then
         s0 = r1/(2.0d0**n*r2*su2)
      else
         s0 = r1/(2.0d0**n*r2*(fl/fs*su1+su2))
      endif
      do k = 1 , kb
         Ck(k) = fl/fs*s0*Ck(k)
      enddo
      do k = kb + 1 , nm
         Ck(k) = s0*Ck(k)
      enddo
      end



!       **********************************

      subroutine sckb(m,n,c,Df,Ck)
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
      implicit none
      real(wp) c , Ck , d1 , d2 , d3 , Df , fac , r , r1 , reg ,&
                     & sum , sw
      integer i , i1 , i2 , ip , k , m , n , nm
      dimension Df(200) , Ck(200)
      if ( c<=1.0d-10 ) c = 1.0d-10
      nm = 25 + int(0.5*(n-m)+c)
      ip = 1
      if ( n-m==2*int((n-m)/2) ) ip = 0
      reg = 1.0d0
      if ( m+nm>80 ) reg = 1.0d-200
      fac = -0.5d0**m
      sw = 0.0d0
      do k = 0 , nm - 1
         fac = -fac
         i1 = 2*k + ip + 1
         r = reg
         do i = i1 , i1 + 2*m - 1
            r = r*i
         enddo
         i2 = k + m + ip
         do i = i2 , i2 + k - 1
            r = r*(i+0.5d0)
         enddo
         sum = r*Df(k+1)
         do i = k + 1 , nm
            d1 = 2.0d0*i + ip
            d2 = 2.0d0*m + d1
            d3 = i + m + ip - 0.5d0
            r = r*d2*(d2-1.0d0)*i*(d3+k)/(d1*(d1-1.0d0)*(i-k)*d3)
            sum = sum + r*Df(i+1)
            if ( abs(sw-sum)<abs(sum)*1.0d-14 ) exit
            sw = sum
         enddo
         r1 = reg
         do i = 2 , m + k
            r1 = r1*i
         enddo
         Ck(k+1) = fac*sum/r1
      enddo
      end



!       **********************************

      subroutine cpdla(n,z,Cdn)
!
!       ===========================================================
!       Purpose: Compute complex parabolic cylinder function Dn(z)
!                for large argument
!       Input:   z   --- Complex argument of Dn(z)
!                n   --- Order of Dn(z) (n = 0,±1,±2,…)
!       Output:  CDN --- Dn(z)
!       ===========================================================
!
      implicit none
      complex(wp) cb0 , Cdn , cr , z
      integer k , n
      cb0 = z**n*exp(-.25d0*z*z)
      cr = (1.0d0,0.0d0)
      Cdn = (1.0d0,0.0d0)
      do k = 1 , 16
         cr = -0.5d0*cr*(2.0*k-n-1.0)*(2.0*k-n-2.0)/(k*z*z)
         Cdn = Cdn + cr
         if ( abs(cr)<abs(Cdn)*1.0d-12 ) exit
      enddo
      Cdn = cb0*Cdn
      end



!       **********************************

      subroutine fcszo(Kf,Nt,Zo)
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
      implicit none
      integer i , it , j , Kf , nr , Nt
      real(wp) pi , psq , px , py , w , w0
      complex(wp) z , zd , zf , zfd , zgd , Zo , zp , zq , zw
      dimension Zo(Nt)
      pi = 3.141592653589793d0
      psq = 0.0d0
      w = 0.0d0
      do nr = 1 , Nt
         if ( Kf==1 ) psq = sqrt(4.0d0*nr-1.0d0)
         if ( Kf==2 ) psq = 2.0d0*nr**(0.5)
         px = psq - log(pi*psq)/(pi*pi*psq**3.0)
         py = log(pi*psq)/(pi*psq)
         z = dcmplx(px,py)
         if ( Kf==2 ) then
            if ( nr==2 ) z = (2.8334,0.2443)
            if ( nr==3 ) z = (3.4674,0.2185)
            if ( nr==4 ) z = (4.0025,0.2008)
         endif
         it = 0
 50      it = it + 1
         if ( Kf==1 ) call cfc(z,zf,zd)
         if ( Kf==2 ) call cfs(z,zf,zd)
         zp = (1.0d0,0.0d0)
         do i = 1 , nr - 1
            zp = zp*(z-Zo(i))
         enddo
         zfd = zf/zp
         zq = (0.0d0,0.0d0)
         do i = 1 , nr - 1
            zw = (1.0d0,0.0d0)
            do j = 1 , nr - 1
               if ( j/=i ) zw = zw*(z-Zo(j))
            enddo
            zq = zq + zw
         enddo
         zgd = (zd-zq*zfd)/zp
         z = z - zfd/zgd
         w0 = w
         w = abs(z)
         if ( it<=50 .and. abs((w-w0)/w)>1.0d-12 ) goto 50
         Zo(nr) = z
      enddo
      end



!       **********************************

      subroutine e1xa(x,e1)
!
!       ============================================
!       Purpose: Compute exponential integral E1(x)
!       Input :  x  --- Argument of E1(x)
!       Output:  E1 --- E1(x) ( x > 0 )
!       ============================================
!
      implicit none
      real(wp) e1 , es1 , es2 , x
      if ( x==0.0 ) then
         e1 = 1.0d+300
      elseif ( x<=1.0 ) then
         e1 = -log(x) + ((((1.07857d-3*x-9.76004d-3)*x+5.519968d-2)*x- &
            & 0.24991055d0)*x+0.99999193d0)*x - 0.57721566d0
      else
         es1 = (((x+8.5733287401d0)*x+18.059016973d0)*x+8.6347608925d0) &
             & *x + 0.2677737343d0
         es2 = (((x+9.5733223454d0)*x+25.6329561486d0)                  &
             & *x+21.0996530827d0)*x + 3.9584969228d0
         e1 = exp(-x)/x*es1/es2
      endif
      end

!       **********************************

      subroutine lpmv0(v,m,x,Pmv)
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
      implicit none
      real(wp) c0 , el , eps , pa , pi , Pmv , pss , psv , pv0 ,&
                     & qr , r , r0 , r1 , r2 , rg , s , s0 , s1 , s2 , v
      real(wp) v0 , vs , x , xq
      integer j , k , m , nv
      pi = 3.141592653589793d0
      el = .5772156649015329d0
      eps = 1.0d-14
      nv = int(v)
      v0 = v - nv
      if ( x==-1.0d0 .and. v/=nv ) then
         if ( m==0 ) Pmv = -1.0d+300
         if ( m/=0 ) Pmv = 1.0d+300
         return
      endif
      c0 = 1.0d0
      if ( m/=0 ) then
         rg = v*(v+m)
         do j = 1 , m - 1
            rg = rg*(v*v-j*j)
         enddo
         xq = sqrt(1.0d0-x*x)
         r0 = 1.0d0
         do j = 1 , m
            r0 = .5d0*r0*xq/j
         enddo
         c0 = r0*rg
      endif
      if ( v0==0.0d0 ) then
!          DLMF 14.3.4, 14.7.17, 15.2.4
         Pmv = 1.0d0
         r = 1.0d0
         do k = 1 , nv - m
            r = 0.5d0*r*(-nv+m+k-1.0d0)*(nv+m+k)/(k*(k+m))*(1.0d0+x)
            Pmv = Pmv + r
         enddo
         Pmv = (-1)**nv*c0*Pmv
      elseif ( x>=-0.35d0 ) then
!             DLMF 14.3.4, 15.2.1
         Pmv = 1.0d0
         r = 1.0d0
         do k = 1 , 100
            r = 0.5d0*r*(-v+m+k-1.0d0)*(v+m+k)/(k*(m+k))*(1.0d0-x)
            Pmv = Pmv + r
            if ( k>12 .and. abs(r/Pmv)<eps ) exit
         enddo
         Pmv = (-1)**m*c0*Pmv
      else
!             DLMF 14.3.5, 15.8.10
         vs = sin(v*pi)/pi
         pv0 = 0.0d0
         if ( m/=0 ) then
            qr = sqrt((1.0d0-x)/(1.0d0+x))
            r2 = 1.0d0
            do j = 1 , m
               r2 = r2*qr*j
            enddo
            s0 = 1.0d0
            r1 = 1.0d0
            do k = 1 , m - 1
               r1 = 0.5d0*r1*(-v+k-1)*(v+k)/(k*(k-m))*(1.0d0+x)
               s0 = s0 + r1
            enddo
            pv0 = -vs*r2/m*s0
         endif
         call psi_spec(v,psv)
         pa = 2.0d0*(psv+el) + pi/tan(pi*v) + 1.0d0/v
         s1 = 0.0d0
         do j = 1 , m
            s1 = s1 + (j*j+v*v)/(j*(j*j-v*v))
         enddo
         Pmv = pa + s1 - 1.0d0/(m-v) + log(0.5d0*(1.0d0+x))
         r = 1.0d0
         do k = 1 , 100
            r = 0.5d0*r*(-v+m+k-1.0d0)*(v+m+k)/(k*(k+m))*(1.0d0+x)
            s = 0.0d0
            do j = 1 , m
               s = s + ((k+j)**2+v*v)/((k+j)*((k+j)**2-v*v))
            enddo
            s2 = 0.0d0
            do j = 1 , k
               s2 = s2 + 1.0d0/(j*(j*j-v*v))
            enddo
            pss = pa + s + 2.0d0*v*v*s2 - 1.0d0/(m+k-v)                 &
                & + log(0.5d0*(1.0d0+x))
            r2 = pss*r
            Pmv = Pmv + r2
            if ( abs(r2/Pmv)<eps ) exit
         enddo
         Pmv = pv0 + Pmv*vs*c0
      endif
      end

!       **********************************

      subroutine lpmv(v,m,x,Pmv)
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
      implicit none
      real(wp) g1 , g2 , p0 , p1 , Pmv , v , v0 , &
                     & vx , x
      integer j , m , mx , neg_m , nv
      if ( x==-1.0d0 .and. v/=int(v) ) then
         if ( m==0 ) Pmv = -dinf()
         if ( m/=0 ) Pmv = dinf()
         return
      endif
      vx = v
      mx = m
!       DLMF 14.9.5
      if ( v<0 ) vx = -vx - 1
      neg_m = 0
      if ( m<0 ) then
         if ( (vx+m+1)>0d0 .or. vx/=int(vx) ) then
            neg_m = 1
            mx = -m
         else
!             We don't handle cases where DLMF 14.9.3 doesn't help
            Pmv = dnan()
            return
         endif
      endif
      nv = int(vx)
      v0 = vx - nv
      if ( nv>2 .and. nv>mx ) then
!          Up-recursion on degree, AMS 8.5.3 / DLMF 14.10.3
         call lpmv0(v0+mx,mx,x,p0)
         call lpmv0(v0+mx+1,mx,x,p1)
         Pmv = p1
         do j = mx + 2 , nv
            Pmv = ((2*(v0+j)-1)*x*p1-(v0+j-1+mx)*p0)/(v0+j-mx)
            p0 = p1
            p1 = Pmv
         enddo
      else
         call lpmv0(vx,mx,x,Pmv)
      endif
      if ( neg_m/=0 .and. abs(Pmv)<1.0d+300 ) then
!          DLMF 14.9.3
         call gamma2(vx-mx+1,g1)
         call gamma2(vx+mx+1,g2)
         Pmv = Pmv*g1/g2*(-1)**mx
      endif
      end


!       **********************************

      subroutine cgama(x,y,Kf,Gr,Gi)
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
      implicit none
      real(wp) a , g0 , Gi , gi1 , Gr , gr1 , pi , si , sr , t ,&
                     & th , th1 , th2 , x , x0 , x1 , y , y1 , z1 , z2
      integer j , k , Kf , na
      dimension a(10)
      pi = 3.141592653589793d0
      data a/8.333333333333333d-02 , -2.777777777777778d-03 ,           &
         & 7.936507936507937d-04 , -5.952380952380952d-04 ,             &
         & 8.417508417508418d-04 , -1.917526917526918d-03 ,             &
         & 6.410256410256410d-03 , -2.955065359477124d-02 ,             &
         & 1.796443723688307d-01 , -1.39243221690590d+00/
      if ( y==0.0d0 .and. x==int(x) .and. x<=0.0d0 ) then
         Gr = 1.0d+300
         Gi = 0.0d0
         return
      elseif ( x<0.0d0 ) then
         x1 = x
         y1 = y
         x = -x
         y = -y
      else
         y1 = 0.0d0
         x1 = x
      endif
      x0 = x
      na = 0
      if ( x<=7.0 ) then
         na = int(7-x)
         x0 = x + na
      endif
      z1 = sqrt(x0*x0+y*y)
      th = atan(y/x0)
      Gr = (x0-.5d0)*log(z1) - th*y - x0 + 0.5d0*log(2.0d0*pi)
      Gi = th*(x0-0.5d0) + y*log(z1) - y
      do k = 1 , 10
         t = z1**(1-2*k)
         Gr = Gr + a(k)*t*cos((2.0d0*k-1.0d0)*th)
         Gi = Gi - a(k)*t*sin((2.0d0*k-1.0d0)*th)
      enddo
      if ( x<=7.0 ) then
         gr1 = 0.0d0
         gi1 = 0.0d0
         do j = 0 , na - 1
            gr1 = gr1 + .5d0*log((x+j)**2+y*y)
            gi1 = gi1 + atan(y/(x+j))
         enddo
         Gr = Gr - gr1
         Gi = Gi - gi1
      endif
      if ( x1<0.0d0 ) then
         z1 = sqrt(x*x+y*y)
         th1 = atan(y/x)
         sr = -sin(pi*x)*dcosh(pi*y)
         si = -cos(pi*x)*dsinh(pi*y)
         z2 = sqrt(sr*sr+si*si)
         th2 = atan(si/sr)
         if ( sr<0.0d0 ) th2 = pi + th2
         Gr = log(pi/(z1*z2)) - Gr
         Gi = -th1 - th2 - Gi
         x = x1
         y = y1
      endif
      if ( Kf==1 ) then
         g0 = exp(Gr)
         Gr = g0*cos(Gi)
         Gi = g0*sin(Gi)
      endif
      end

!       **********************************

      subroutine aswfb(m,n,c,x,Kd,Cv,S1f,S1d)
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
      implicit none
      real(wp) c , Cv , df , eps , pd , pm , S1d , S1f , su1 ,  &
                     & sw , x
      integer ip , k , Kd , m , mk , n , nm , nm2
      dimension df(200) , pm(0:251) , pd(0:251)
      eps = 1.0d-14
      ip = 1
      if ( n-m==2*int((n-m)/2) ) ip = 0
      nm = 25 + int((n-m)/2+c)
      nm2 = 2*nm + m
      call sdmn(m,n,c,Cv,Kd,df)
      call lpmns(m,nm2,x,pm,pd)
      sw = 0.0d0
      su1 = 0.0d0
      do k = 1 , nm
         mk = m + 2*(k-1) + ip
         su1 = su1 + df(k)*pm(mk)
         if ( abs(sw-su1)<abs(su1)*eps ) exit
         sw = su1
      enddo
      S1f = (-1)**m*su1
      su1 = 0.0d0
      do k = 1 , nm
         mk = m + 2*(k-1) + ip
         su1 = su1 + df(k)*pd(mk)
         if ( abs(sw-su1)<abs(su1)*eps ) exit
         sw = su1
      enddo
      S1d = (-1)**m*su1
      end



!       **********************************

      subroutine chgus(a,b,x,Hu,Id)
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
      implicit none
      real(wp) a , b , d1 , d2 , ga , gab , gb , gb2 , h0 ,     &
                     & hmax , hmin , Hu , hu0 , hua , pi , r1 , r2 , x ,&
                     & xg1 , xg2
      integer Id , j
      Id = -100
      pi = 3.141592653589793d0
      call gamma2(a,ga)
      call gamma2(b,gb)
      xg1 = 1.0d0 + a - b
      call gamma2(xg1,gab)
      xg2 = 2.0d0 - b
      call gamma2(xg2,gb2)
      hu0 = pi/sin(pi*b)
      r1 = hu0/(gab*gb)
      r2 = hu0*x**(1.0d0-b)/(ga*gb2)
      Hu = r1 - r2
      hmax = 0.0d0
      hmin = 1.0d+300
      h0 = 0.0d0
      do j = 1 , 150
         r1 = r1*(a+j-1.0d0)/(j*(b+j-1.0d0))*x
         r2 = r2*(a-b+j)/(j*(1.0d0-b+j))*x
         Hu = Hu + r1 - r2
         hua = abs(Hu)
         if ( hua>hmax ) hmax = hua
         if ( hua<hmin ) hmin = hua
         if ( abs(Hu-h0)<abs(Hu)*1.0d-15 ) exit
         h0 = Hu
      enddo
      d1 = log10(hmax)
      d2 = 0.0d0
      if ( hmin/=0.0 ) d2 = log10(hmin)
      Id = 15 - abs(d1-d2)
      end



!       **********************************

      subroutine itth0(x,Tth)
!
!       ===========================================================
!       Purpose: Evaluate the integral H0(t)/t with respect to t
!                from x to infinity
!       Input :  x   --- Lower limit  ( x ≥ 0 )
!       Output:  TTH --- Integration of H0(t)/t from x to infinity
!       ===========================================================
!
      implicit none
      real(wp) f0 , g0 , pi , r , s , t , Tth , tty , x , xt
      integer k
      pi = 3.141592653589793d0
      s = 1.0d0
      r = 1.0d0
      if ( x<24.5d0 ) then
         do k = 1 , 60
            r = -r*x*x*(2.0*k-1.0d0)/(2.0*k+1.0d0)**3
            s = s + r
            if ( abs(r)<abs(s)*1.0d-12 ) exit
         enddo
         Tth = pi/2.0d0 - 2.0d0/pi*x*s
      else
         do k = 1 , 10
            r = -r*(2.0*k-1.0d0)**3/((2.0*k+1.0d0)*x*x)
            s = s + r
            if ( abs(r)<abs(s)*1.0d-12 ) exit
         enddo
         Tth = 2.0d0/(pi*x)*s
         t = 8.0d0/x
         xt = x + .25d0*pi
         f0 = (((((.18118d-2*t-.91909d-2)*t+.017033d0)*t-.9394d-3)      &
            & *t-.051445d0)*t-.11d-5)*t + .7978846d0
         g0 = (((((-.23731d-2*t+.59842d-2)*t+.24437d-2)*t-.0233178d0)   &
            & *t+.595d-4)*t+.1620695d0)*t
         tty = (f0*sin(xt)-g0*cos(xt))/(sqrt(x)*x)
         Tth = Tth + tty
      endif
      end

!       **********************************

      subroutine lgama(Kf,x,Gl)
!
!       ==================================================
!       Purpose: Compute gamma function Г(x) or ln[Г(x)]
!       Input:   x  --- Argument of Г(x) ( x > 0 )
!                KF --- Function code
!                       KF=1 for Г(x); KF=0 for ln[Г(x)]
!       Output:  GL --- Г(x) or ln[Г(x)]
!       ==================================================
!
      implicit none
      real(wp) a , Gl , gl0 , x , x0 , x2 , xp
      integer k , Kf , n
      dimension a(10)
      data a/8.333333333333333d-02 , -2.777777777777778d-03 ,           &
         & 7.936507936507937d-04 , -5.952380952380952d-04 ,             &
         & 8.417508417508418d-04 , -1.917526917526918d-03 ,             &
         & 6.410256410256410d-03 , -2.955065359477124d-02 ,             &
         & 1.796443723688307d-01 , -1.39243221690590d+00/
      x0 = x
      n = 0
      if ( x==1.0 .or. x==2.0 ) then
         Gl = 0.0d0
         goto 100
      elseif ( x<=7.0 ) then
         n = int(7-x)
         x0 = x + n
      endif
      x2 = 1.0d0/(x0*x0)
      xp = 6.283185307179586477d0
      gl0 = a(10)
      do k = 9 , 1 , -1
         gl0 = gl0*x2 + a(k)
      enddo
      Gl = gl0/x0 + 0.5d0*log(xp) + (x0-.5d0)*log(x0) - x0
      if ( x<=7.0 ) then
         do k = 1 , n
            Gl = Gl - log(x0-1.0d0)
            x0 = x0 - 1.0d0
         enddo
      endif
 100  if ( Kf==1 ) Gl = exp(Gl)
      end

!       **********************************

      subroutine lqna(n,x,Qn,Qd)
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
      implicit none
      integer k , n
      real(wp) q0 , q1 , Qd , qf , Qn , x
      dimension Qn(0:n) , Qd(0:n)
      if ( abs(x)==1.0d0 ) then
         do k = 0 , n
            Qn(k) = 1.0d+300
            Qd(k) = -1.0d+300
         enddo
      elseif ( abs(x)<1.0d0 ) then
         q0 = 0.5d0*log((1.0d0+x)/(1.0d0-x))
         q1 = x*q0 - 1.0d0
         Qn(0) = q0
         Qn(1) = q1
         Qd(0) = 1.0d0/(1.0d0-x*x)
         Qd(1) = Qn(0) + x*Qd(0)
         do k = 2 , n
            qf = ((2*k-1)*x*q1-(k-1)*q0)/k
            Qn(k) = qf
            Qd(k) = (Qn(k-1)-x*qf)*k/(1.0d0-x*x)
            q0 = q1
            q1 = qf
         enddo
      endif
      end

!       **********************************

      subroutine dvla(Va,x,Pd)
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
      implicit none
      real(wp) a0 , ep , eps , gl , Pd , pi , r , Va , vl , x , &
                     & x1
      integer k
      pi = 3.141592653589793d0
      eps = 1.0d-12
      ep = exp(-.25*x*x)
      a0 = abs(x)**Va*ep
      r = 1.0d0
      Pd = 1.0d0
      do k = 1 , 16
         r = -0.5d0*r*(2.0*k-Va-1.0)*(2.0*k-Va-2.0)/(k*x*x)
         Pd = Pd + r
         if ( abs(r/Pd)<eps ) exit
      enddo
      Pd = a0*Pd
      if ( x<0.0d0 ) then
         x1 = -x
         call vvla(Va,x1,vl)
         call gamma2(-Va,gl)
         Pd = pi*vl/gl + cos(pi*Va)*Pd
      endif
      end



!       **********************************

      subroutine ik01a(x,Bi0,Di0,Bi1,Di1,Bk0,Dk0,Bk1,Dk1)
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
      implicit none
      real(wp) a , a1 , b , Bi0 , Bi1 , Bk0 , Bk1 , ca , cb ,   &
                     & ct , Di0 , Di1 , Dk0 , Dk1 , el , pi , r , w0 ,  &
                     & ww , x
      real(wp) x2 , xr , xr2
      integer k , k0
      dimension a(12) , b(12) , a1(8)
      pi = 3.141592653589793d0
      el = 0.5772156649015329d0
      x2 = x*x
      if ( x==0.0d0 ) then
         Bi0 = 1.0d0
         Bi1 = 0.0d0
         Bk0 = 1.0d+300
         Bk1 = 1.0d+300
         Di0 = 0.0d0
         Di1 = 0.5d0
         Dk0 = -1.0d+300
         Dk1 = -1.0d+300
         return
      elseif ( x<=18.0d0 ) then
         Bi0 = 1.0d0
         r = 1.0d0
         do k = 1 , 50
            r = 0.25d0*r*x2/(k*k)
            Bi0 = Bi0 + r
            if ( abs(r/Bi0)<1.0d-15 ) exit
         enddo
         Bi1 = 1.0d0
         r = 1.0d0
         do k = 1 , 50
            r = 0.25d0*r*x2/(k*(k+1))
            Bi1 = Bi1 + r
            if ( abs(r/Bi1)<1.0d-15 ) exit
         enddo
         Bi1 = 0.5d0*x*Bi1
      else
         data a/0.125d0 , 7.03125d-2 , 7.32421875d-2 ,                  &
            & 1.1215209960938d-1 , 2.2710800170898d-1 ,                 &
            & 5.7250142097473d-1 , 1.7277275025845d0 ,                  &
            & 6.0740420012735d0 , 2.4380529699556d01 ,                  &
            & 1.1001714026925d02 , 5.5133589612202d02 ,                 &
            & 3.0380905109224d03/
         data b/ - 0.375d0 , -1.171875d-1 , -1.025390625d-1 ,           &
            & -1.4419555664063d-1 , -2.7757644653320d-1 ,               &
            & -6.7659258842468d-1 , -1.9935317337513d0 ,                &
            & -6.8839142681099d0 , -2.7248827311269d01 ,                &
            & -1.2159789187654d02 , -6.0384407670507d02 ,               &
            & -3.3022722944809d03/
         k0 = 12
         if ( x>=35.0 ) k0 = 9
         if ( x>=50.0 ) k0 = 7
         ca = exp(x)/sqrt(2.0d0*pi*x)
         Bi0 = 1.0d0
         xr = 1.0d0/x
         do k = 1 , k0
            Bi0 = Bi0 + a(k)*xr**k
         enddo
         Bi0 = ca*Bi0
         Bi1 = 1.0d0
         do k = 1 , k0
            Bi1 = Bi1 + b(k)*xr**k
         enddo
         Bi1 = ca*Bi1
      endif
      ww = 0.0d0
      if ( x<=9.0d0 ) then
         ct = -(log(x/2.0d0)+el)
         Bk0 = 0.0d0
         w0 = 0.0d0
         r = 1.0d0
         do k = 1 , 50
            w0 = w0 + 1.0d0/k
            r = 0.25d0*r/(k*k)*x2
            Bk0 = Bk0 + r*(w0+ct)
            if ( abs((Bk0-ww)/Bk0)<1.0d-15 ) exit
            ww = Bk0
         enddo
         Bk0 = Bk0 + ct
      else
         data a1/0.125d0 , 0.2109375d0 , 1.0986328125d0 ,               &
            & 1.1775970458984d01 , 2.1461706161499d02 ,                 &
            & 5.9511522710323d03 , 2.3347645606175d05 ,                 &
            & 1.2312234987631d07/
         cb = 0.5d0/x
         xr2 = 1.0d0/x2
         Bk0 = 1.0d0
         do k = 1 , 8
            Bk0 = Bk0 + a1(k)*xr2**k
         enddo
         Bk0 = cb*Bk0/Bi0
      endif
      Bk1 = (1.0d0/x-Bi1*Bk0)/Bi0
      Di0 = Bi1
      Di1 = Bi0 - Bi1/x
      Dk0 = -Bk1
      Dk1 = -Bk0 - Bk1/x
      end

!       **********************************

      subroutine cpbdn(n,z,Cpb,Cpd)
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
      implicit none
      real(wp) a0 , pi , x
      complex(wp) c0 , ca0 , cf , cf0 , cf1 , cfa , cfb , Cpb , Cpd ,    &
               & cs0 , z , z1
      integer k , m , n , n0 , n1 , nm1
      dimension Cpb(0:*) , Cpd(0:*)
      pi = 3.141592653589793d0
      x = dble(z)
      a0 = abs(z)
      c0 = (0.0d0,0.0d0)
      ca0 = exp(-0.25d0*z*z)
      n0 = 0
      if ( n>=0 ) then
         cf0 = ca0
         cf1 = z*ca0
         Cpb(0) = cf0
         Cpb(1) = cf1
         do k = 2 , n
            cf = z*cf1 - (k-1.0d0)*cf0
            Cpb(k) = cf
            cf0 = cf1
            cf1 = cf
         enddo
      else
         n0 = -n
         if ( x<=0.0 .or. abs(z)==0.0 ) then
            cf0 = ca0
            Cpb(0) = cf0
            z1 = -z
            if ( a0<=7.0 ) then
               call cpdsa(-1,z1,cf1)
            else
               call cpdla(-1,z1,cf1)
            endif
            cf1 = sqrt(2.0d0*pi)/ca0 - cf1
            Cpb(1) = cf1
            do k = 2 , n0
               cf = (-z*cf1+cf0)/(k-1.0d0)
               Cpb(k) = cf
               cf0 = cf1
               cf1 = cf
            enddo
         elseif ( a0<=3.0 ) then
            call cpdsa(-n0,z,cfa)
            Cpb(n0) = cfa
            n1 = n0 + 1
            call cpdsa(-n1,z,cfb)
            Cpb(n1) = cfb
            nm1 = n0 - 1
            do k = nm1 , 0 , -1
               cf = z*cfa + (k+1.0d0)*cfb
               Cpb(k) = cf
               cfb = cfa
               cfa = cf
            enddo
         else
            m = 100 + abs(n)
            cfa = c0
            cfb = (1.0d-30,0.0d0)
            do k = m , 0 , -1
               cf = z*cfb + (k+1.0d0)*cfa
               if ( k<=n0 ) Cpb(k) = cf
               cfa = cfb
               cfb = cf
            enddo
            cs0 = ca0/cf
            do k = 0 , n0
               Cpb(k) = cs0*Cpb(k)
            enddo
         endif
      endif
      Cpd(0) = -0.5d0*z*Cpb(0)
      if ( n>=0 ) then
         do k = 1 , n
            Cpd(k) = -0.5d0*z*Cpb(k) + k*Cpb(k-1)
         enddo
      else
         do k = 1 , n0
            Cpd(k) = 0.5d0*z*Cpb(k) - Cpb(k-1)
         enddo
      endif
      end



!       **********************************

      subroutine ik01b(x,Bi0,Di0,Bi1,Di1,Bk0,Dk0,Bk1,Dk1)
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
      implicit none
      real(wp) Bi0 , Bi1 , Bk0 , Bk1 , Di0 , Di1 , Dk0 , Dk1 ,  &
                     & t , t2 , x
      if ( x==0.0d0 ) then
         Bi0 = 1.0d0
         Bi1 = 0.0d0
         Bk0 = 1.0d+300
         Bk1 = 1.0d+300
         Di0 = 0.0d0
         Di1 = 0.5d0
         Dk0 = -1.0d+300
         Dk1 = -1.0d+300
         return
      elseif ( x<=3.75d0 ) then
         t = x/3.75d0
         t2 = t*t
         Bi0 = (((((.0045813d0*t2+.0360768d0)*t2+.2659732d0)*t2+        &
             & 1.2067492d0)*t2+3.0899424d0)*t2+3.5156229d0)*t2 + 1.0d0
         Bi1 = x*((((((.00032411d0*t2+.00301532d0)*t2+.02658733d0)*t2+  &
             & .15084934d0)*t2+.51498869d0)*t2+.87890594d0)*t2+.5d0)
      else
         t = 3.75d0/x
         Bi0 = ((((((((.00392377d0*t-.01647633d0)*t+.02635537d0)*t-     &
             & .02057706d0)*t+.916281d-2)*t-.157565d-2)*t+.225319d-2)   &
             & *t+.01328592d0)*t+.39894228d0)*exp(x)/sqrt(x)
         Bi1 = ((((((((-.420059d-2*t+.01787654d0)*t-.02895312d0)*t+     &
             & .02282967d0)*t-.01031555d0)*t+.163801d-2)*t-.00362018d0) &
             & *t-.03988024d0)*t+.39894228d0)*exp(x)/sqrt(x)
      endif
      if ( x<=2.0d0 ) then
         t = x/2.0d0
         t2 = t*t
         Bk0 = (((((.0000074d0*t2+.0001075d0)*t2+.00262698d0)*t2+       &
             & .0348859d0)*t2+.23069756d0)*t2+.4227842d0)               &
             & *t2 - .57721566d0 - Bi0*log(t)
         Bk1 = ((((((-.00004686d0*t2-.00110404d0)*t2-.01919402d0)*t2-   &
             & .18156897d0)*t2-.67278579d0)*t2+.15443144d0)*t2+1.0d0)   &
             & /x + Bi1*log(t)
      else
         t = 2.0d0/x
         t2 = t*t
         Bk0 = ((((((.00053208d0*t-.0025154d0)*t+.00587872d0)*t-        &
             & .01062446d0)*t+.02189568d0)*t-.07832358d0)               &
             & *t+1.25331414d0)*exp(-x)/sqrt(x)
         Bk1 = ((((((-.00068245d0*t+.00325614d0)*t-.00780353d0)*t+      &
             & .01504268d0)*t-.0365562d0)*t+.23498619d0)*t+1.25331414d0)&
             & *exp(-x)/sqrt(x)
      endif
      Di0 = Bi1
      Di1 = Bi0 - Bi1/x
      Dk0 = -Bk1
      Dk1 = -Bk0 - Bk1/x
      end

!       **********************************

      subroutine beta(p,q,Bt)
!
!       ==========================================
!       Purpose: Compute the beta function B(p,q)
!       Input :  p  --- Parameter  ( p > 0 )
!                q  --- Parameter  ( q > 0 )
!       Output:  BT --- B(p,q)
!       Routine called: GAMMA2 for computing Г(x)
!       ==========================================
!
      implicit none
      real(wp) Bt , gp , gpq , gq , p , ppq , q
      call gamma2(p,gp)
      call gamma2(q,gq)
      ppq = p + q
      call gamma2(ppq,gpq)
      Bt = gp*gq/gpq
      end



!       **********************************

      subroutine lpn(n,x,Pn,Pd)
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
      implicit none
      integer k , n
      real(wp) p0 , p1 , Pd , pf , Pn , x
      dimension Pn(0:n) , Pd(0:n)
      Pn(0) = 1.0d0
      Pn(1) = x
      Pd(0) = 0.0d0
      Pd(1) = 1.0d0
      p0 = 1.0d0
      p1 = x
      do k = 2 , n
         pf = (2.0d0*k-1.0d0)/k*x*p1 - (k-1.0d0)/k*p0
         Pn(k) = pf
         if ( abs(x)==1.0d0 ) then
            Pd(k) = 0.5d0*x**(k+1)*k*(k+1.0d0)
         else
            Pd(k) = k*(p1-x*pf)/(1.0d0-x*x)
         endif
         p0 = p1
         p1 = pf
      enddo
      end

!       **********************************

      subroutine fcoef(Kd,m,q,a,Fc)
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
      implicit none
      real(wp) a , f , f1 , f2 , f3 , Fc , fnan , q ,    &
                     & qm , s , s0 , sp , ss , u , v
      integer i , j , jm , k , kb , Kd , km , m
      dimension Fc(251)
      do i = 1 , 251
         Fc(i) = 0.0d0
      enddo
      if ( abs(q)<=1.0d-7 ) then
!          Expansion up to order Q^1 (Abramowitz & Stegun 20.2.27-28)
         if ( Kd==1 ) then
            jm = m/2 + 1
         elseif ( Kd==2 .or. Kd==3 ) then
            jm = (m-1)/2 + 1
         elseif ( Kd==4 ) then
            jm = m/2
         endif
!          Check for overflow
         if ( jm+1>251 ) then
            fnan = dnan()
            do i = 1 , 251
               Fc(i) = fnan
            enddo
            return
         endif
!          Proceed using the simplest expansion
         if ( Kd==1 .or. Kd==2 ) then
            if ( m==0 ) then
               Fc(1) = 1/sqrt(2.0d0)
               Fc(2) = -q/2.0d0/sqrt(2.0d0)
            elseif ( m==1 ) then
               Fc(1) = 1.0d0
               Fc(2) = -q/8.0d0
            elseif ( m==2 ) then
               Fc(1) = q/4.0d0
               Fc(2) = 1.0d0
               Fc(3) = -q/12.0d0
            else
               Fc(jm) = 1.0d0
               Fc(jm+1) = -q/(4.0d0*(m+1))
               Fc(jm-1) = q/(4.0d0*(m-1))
            endif
         elseif ( Kd==3 .or. Kd==4 ) then
            if ( m==1 ) then
               Fc(1) = 1.0d0
               Fc(2) = -q/8.0d0
            elseif ( m==2 ) then
               Fc(1) = 1.0d0
               Fc(2) = -q/12.0d0
            else
               Fc(jm) = 1.0d0
               Fc(jm+1) = -q/(4.0d0*(m+1))
               Fc(jm-1) = q/(4.0d0*(m-1))
            endif
         endif
         return
      elseif ( q<=1.0d0 ) then
         qm = 7.5 + 56.1*sqrt(q) - 134.7*q + 90.7*sqrt(q)*q
      else
         qm = 17.0 + 3.1*sqrt(q) - .126*q + .0037*sqrt(q)*q
      endif
      km = int(qm+0.5*m)
      if ( km>251 ) then
!          Overflow, generate NaNs
         fnan = dnan()
         do i = 1 , 251
            Fc(i) = fnan
         enddo
         return
      endif
      kb = 0
      s = 0.0d0
      f = 1.0d-100
      u = 0.0d0
      Fc(km) = 0.0d0
      f2 = 0.0d0
      if ( Kd==1 ) then
         do k = km , 3 , -1
            v = u
            u = f
            f = (a-4.0d0*k*k)*u/q - v
            if ( abs(f)<abs(Fc(k+1)) ) then
               kb = k
               Fc(1) = 1.0d-100
               sp = 0.0d0
               f3 = Fc(k+1)
               Fc(2) = a/q*Fc(1)
               Fc(3) = (a-4.0d0)*Fc(2)/q - 2.0d0*Fc(1)
               u = Fc(2)
               f1 = Fc(3)
               do i = 3 , kb
                  v = u
                  u = f1
                  f1 = (a-4.0d0*(i-1.0d0)**2)*u/q - v
                  Fc(i+1) = f1
                  if ( i==kb ) f2 = f1
                  if ( i/=kb ) sp = sp + f1*f1
               enddo
               sp = sp + 2.0d0*Fc(1)**2 + Fc(2)**2 + Fc(3)**2
               ss = s + sp*(f3/f2)**2
               s0 = sqrt(1.0d0/ss)
               do j = 1 , km
                  if ( j<=kb+1 ) then
                     Fc(j) = s0*Fc(j)*f3/f2
                  else
                     Fc(j) = s0*Fc(j)
                  endif
               enddo
               goto 200
            else
               Fc(k) = f
               s = s + f*f
            endif
         enddo
         Fc(2) = q*Fc(3)/(a-4.0d0-2.0d0*q*q/a)
         Fc(1) = q/a*Fc(2)
         s = s + 2.0d0*Fc(1)**2 + Fc(2)**2
         s0 = sqrt(1.0d0/s)
         do k = 1 , km
            Fc(k) = s0*Fc(k)
         enddo
      elseif ( Kd==2 .or. Kd==3 ) then
         do k = km , 3 , -1
            v = u
            u = f
            f = (a-(2.0d0*k-1)**2)*u/q - v
            if ( abs(f)>=abs(Fc(k)) ) then
               Fc(k-1) = f
               s = s + f*f
            else
               kb = k
               f3 = Fc(k)
               goto 50
            endif
         enddo
         Fc(1) = q/(a-1.0d0-(-1)**Kd*q)*Fc(2)
         s = s + Fc(1)*Fc(1)
         s0 = sqrt(1.0d0/s)
         do k = 1 , km
            Fc(k) = s0*Fc(k)
         enddo
         goto 200
 50      Fc(1) = 1.0d-100
         Fc(2) = (a-1.0d0-(-1)**Kd*q)/q*Fc(1)
         sp = 0.0d0
         u = Fc(1)
         f1 = Fc(2)
         do i = 2 , kb - 1
            v = u
            u = f1
            f1 = (a-(2.0d0*i-1.0d0)**2)*u/q - v
            if ( i/=kb-1 ) then
               Fc(i+1) = f1
               sp = sp + f1*f1
            else
               f2 = f1
            endif
         enddo
         sp = sp + Fc(1)**2 + Fc(2)**2
         ss = s + sp*(f3/f2)**2
         s0 = 1.0d0/sqrt(ss)
         do j = 1 , km
            if ( j<kb ) Fc(j) = s0*Fc(j)*f3/f2
            if ( j>=kb ) Fc(j) = s0*Fc(j)
         enddo
      elseif ( Kd==4 ) then
         do k = km , 3 , -1
            v = u
            u = f
            f = (a-4.0d0*k*k)*u/q - v
            if ( abs(f)>=abs(Fc(k)) ) then
               Fc(k-1) = f
               s = s + f*f
            else
               kb = k
               f3 = Fc(k)
               goto 100
            endif
         enddo
         Fc(1) = q/(a-4.0d0)*Fc(2)
         s = s + Fc(1)*Fc(1)
         s0 = sqrt(1.0d0/s)
         do k = 1 , km
            Fc(k) = s0*Fc(k)
         enddo
         goto 200
 100     Fc(1) = 1.0d-100
         Fc(2) = (a-4.0d0)/q*Fc(1)
         sp = 0.0d0
         u = Fc(1)
         f1 = Fc(2)
         do i = 2 , kb - 1
            v = u
            u = f1
            f1 = (a-4.0d0*i*i)*u/q - v
            if ( i/=kb-1 ) then
               Fc(i+1) = f1
               sp = sp + f1*f1
            else
               f2 = f1
            endif
         enddo
         sp = sp + Fc(1)**2 + Fc(2)**2
         ss = s + sp*(f3/f2)**2
         s0 = 1.0d0/sqrt(ss)
         do j = 1 , km
            if ( j<kb ) Fc(j) = s0*Fc(j)*f3/f2
            if ( j>=kb ) Fc(j) = s0*Fc(j)
         enddo
      endif
 200  if ( Fc(1)<0.0d0 ) then
         do j = 1 , km
            Fc(j) = -Fc(j)
         enddo
      endif
      end



!       **********************************

      subroutine sphi(n,x,Nm,Si,Di)
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
      implicit none
      real(wp) cs , Di , f , f0 , f1 , Si , si0 , x
      integer k , m , n , Nm
      dimension Si(0:n) , Di(0:n)
      Nm = n
      if ( abs(x)<1.0d-100 ) then
         do k = 0 , n
            Si(k) = 0.0d0
            Di(k) = 0.0d0
         enddo
         Si(0) = 1.0d0
         Di(1) = 0.333333333333333d0
         return
      endif
      Si(0) = dsinh(x)/x
      Si(1) = -(dsinh(x)/x-dcosh(x))/x
      si0 = Si(0)
      if ( n>=2 ) then
         m = msta1(x,200)
         if ( m<n ) then
            Nm = m
         else
            m = msta2(x,n,15)
         endif
         f = 0.0d0
         f0 = 0.0d0
         f1 = 1.0d0 - 100
         do k = m , 0 , -1
            f = (2.0d0*k+3.0d0)*f1/x + f0
            if ( k<=Nm ) Si(k) = f
            f0 = f1
            f1 = f
         enddo
         cs = si0/f
         do k = 0 , Nm
            Si(k) = cs*Si(k)
         enddo
      endif
      Di(0) = Si(1)
      do k = 1 , Nm
         Di(k) = Si(k-1) - (k+1.0d0)/x*Si(k)
      enddo
      end



!       **********************************

      subroutine pbwa(a,x,W1f,W1d,W2f,W2d)
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
      implicit none
      real(wp) a , d , d1 , d2 , dl , eps , f1 , f2 , g1 , g2 , &
                     & h , h0 , h1 , hl , p0 , r , r1 , ugi , ugr , vgi
      real(wp) vgr , W1d , W1f , W2d , W2f , x , x1 , x2 , y1 , &
                     & y1d , y1f , y2d , y2f
      integer k , l1 , l2 , m
      dimension h(100) , d(80)
      eps = 1.0d-15
      p0 = 0.59460355750136d0
      if ( a==0.0d0 ) then
         g1 = 3.625609908222d0
         g2 = 1.225416702465d0
      else
         x1 = 0.25d0
         y1 = 0.5d0*a
         call cgama(x1,y1,1,ugr,ugi)
         g1 = sqrt(ugr*ugr+ugi*ugi)
         x2 = 0.75d0
         call cgama(x2,y1,1,vgr,vgi)
         g2 = sqrt(vgr*vgr+vgi*vgi)
      endif
      f1 = sqrt(g1/g2)
      f2 = sqrt(2.0d0*g2/g1)
      h0 = 1.0d0
      h1 = a
      h(1) = a
      do l1 = 4 , 200 , 2
         m = l1/2
         hl = a*h1 - 0.25d0*(l1-2.0d0)*(l1-3.0d0)*h0
         h(m) = hl
         h0 = h1
         h1 = hl
      enddo
      y1f = 1.0d0
      r = 1.0d0
      do k = 1 , 100
         r = 0.5d0*r*x*x/(k*(2.0d0*k-1.0d0))
         r1 = h(k)*r
         y1f = y1f + r1
         if ( abs(r1)<=eps*abs(y1f) .and. k>30 ) exit
      enddo
      y1d = a
      r = 1.0d0
      do k = 1 , 99
         r = 0.5d0*r*x*x/(k*(2.0d0*k+1.0d0))
         r1 = h(k+1)*r
         y1d = y1d + r1
         if ( abs(r1)<=eps*abs(y1d) .and. k>30 ) exit
      enddo
      y1d = x*y1d
      d1 = 1.0d0
      d2 = a
      d(1) = 1.0d0
      d(2) = a
      do l2 = 5 , 160 , 2
         m = (l2+1)/2
         dl = a*d2 - 0.25d0*(l2-2.0d0)*(l2-3.0d0)*d1
         d(m) = dl
         d1 = d2
         d2 = dl
      enddo
      y2f = 1.0d0
      r = 1.0d0
      do k = 1 , 79
         r = 0.5d0*r*x*x/(k*(2.0d0*k+1.0d0))
         r1 = d(k+1)*r
         y2f = y2f + r1
         if ( abs(r1)<=eps*abs(y2f) .and. k>30 ) exit
      enddo
      y2f = x*y2f
      y2d = 1.0d0
      r = 1.0d0
      do k = 1 , 79
         r = 0.5d0*r*x*x/(k*(2.0d0*k-1.0d0))
         r1 = d(k+1)*r
         y2d = y2d + r1
         if ( abs(r1)<=eps*abs(y2f) .and. k>30 ) exit
      enddo
      W1f = p0*(f1*y1f-f2*y2f)
      W2f = p0*(f1*y1f+f2*y2f)
      W1d = p0*(f1*y1d-f2*y2d)
      W2d = p0*(f1*y1d+f2*y2d)
      end



!       **********************************

      subroutine rmn1(m,n,c,x,Df,Kd,R1f,R1d)
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
      implicit none
      real(wp) a0 , b0 , c , ck , cx , Df , dj , eps , r , r0 , &
                     & r1 , R1d , R1f , r2 , r3 , reg , sa0 , sj , suc ,&
                     & sud
      real(wp) sum , sw , sw1 , x
      integer ip , j , k , Kd , l , lg , m , n , nm , nm1 , nm2 , np
      dimension ck(200) , Df(200) , sj(0:251) , dj(0:251)
      eps = 1.0d-14
      ip = 1
      nm1 = int((n-m)/2)
      if ( n-m==2*nm1 ) ip = 0
      nm = 25 + nm1 + int(c)
      reg = 1.0d0
      if ( m+nm>80 ) reg = 1.0d-200
      r0 = reg
      do j = 1 , 2*m + ip
         r0 = r0*j
      enddo
      r = r0
      suc = r*Df(1)
      sw = 0.0d0
      do k = 2 , nm
         r = r*(m+k-1.0)*(m+k+ip-1.5d0)/(k-1.0d0)/(k+ip-1.5d0)
         suc = suc + r*Df(k)
         if ( k>nm1 .and. abs(suc-sw)<abs(suc)*eps ) exit
         sw = suc
      enddo
      if ( x==0.0 ) then
         call sckb(m,n,c,Df,ck)
         sum = 0.0d0
         sw1 = 0.0d0
         do j = 1 , nm
            sum = sum + ck(j)
            if ( abs(sum-sw1)<abs(sum)*eps ) exit
            sw1 = sum
         enddo
         r1 = 1.0d0
         do j = 1 , (n+m+ip)/2
            r1 = r1*(j+0.5d0*(n+m+ip))
         enddo
         r2 = 1.0d0
         do j = 1 , m
            r2 = 2.0d0*c*r2*j
         enddo
         r3 = 1.0d0
         do j = 1 , (n-m-ip)/2
            r3 = r3*j
         enddo
         sa0 = (2.0*(m+ip)+1.0)*r1/(2.0**n*c**ip*r2*r3)
         if ( ip==0 ) then
            R1f = sum/(sa0*suc)*Df(1)*reg
            R1d = 0.0d0
         elseif ( ip==1 ) then
            R1f = 0.0d0
            R1d = sum/(sa0*suc)*Df(1)*reg
         endif
         return
      endif
      cx = c*x
      nm2 = 2*nm + m
      call sphj(nm2,cx,nm2,sj,dj)
      a0 = (1.0d0-Kd/(x*x))**(0.5d0*m)/suc
      R1f = 0.0d0
      sw = 0.0d0
      lg = 0
      do k = 1 , nm
         l = 2*k + m - n - 2 + ip
         if ( l==4*int(l/4) ) lg = 1
         if ( l/=4*int(l/4) ) lg = -1
         if ( k==1 ) then
            r = r0
         else
            r = r*(m+k-1.0)*(m+k+ip-1.5d0)/(k-1.0d0)/(k+ip-1.5d0)
         endif
         np = m + 2*k - 2 + ip
         R1f = R1f + lg*r*Df(k)*sj(np)
         if ( k>nm1 .and. abs(R1f-sw)<abs(R1f)*eps ) exit
         sw = R1f
      enddo
      R1f = R1f*a0
      b0 = Kd*m/x**3.0d0/(1.0-Kd/(x*x))*R1f
      sud = 0.0d0
      sw = 0.0d0
      do k = 1 , nm
         l = 2*k + m - n - 2 + ip
         if ( l==4*int(l/4) ) lg = 1
         if ( l/=4*int(l/4) ) lg = -1
         if ( k==1 ) then
            r = r0
         else
            r = r*(m+k-1.0)*(m+k+ip-1.5d0)/(k-1.0d0)/(k+ip-1.5d0)
         endif
         np = m + 2*k - 2 + ip
         sud = sud + lg*r*Df(k)*dj(np)
         if ( k>nm1 .and. abs(sud-sw)<abs(sud)*eps ) exit
         sw = sud
      enddo
      R1d = b0 + a0*c*sud
      end



!       **********************************

      subroutine dvsa(Va,x,Pd)
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
      implicit none
      real(wp) a0 , ep , eps , g0 , g1 , ga0 , gm , Pd , pi ,   &
                     & r , r1 , sq2 , Va , va0 , vm , vt , x
      integer m
      eps = 1.0d-15
      pi = 3.141592653589793d0
      sq2 = sqrt(2.0d0)
      ep = exp(-.25d0*x*x)
      va0 = 0.5d0*(1.0d0-Va)
      if ( Va==0.0 ) then
         Pd = ep
      elseif ( x==0.0 ) then
         if ( va0<=0.0 .and. va0==int(va0) ) then
            Pd = 0.0d0
         else
            call gamma2(va0,ga0)
            Pd = sqrt(pi)/(2.0d0**(-.5d0*Va)*ga0)
         endif
      else
         call gamma2(-Va,g1)
         a0 = 2.0d0**(-0.5d0*Va-1.0d0)*ep/g1
         vt = -.5d0*Va
         call gamma2(vt,g0)
         Pd = g0
         r = 1.0d0
         do m = 1 , 250
            vm = .5d0*(m-Va)
            call gamma2(vm,gm)
            r = -r*sq2*x/m
            r1 = gm*r
            Pd = Pd + r1
            if ( abs(r1)<abs(Pd)*eps ) exit
         enddo
         Pd = a0*Pd
      endif
      end



!       **********************************

      subroutine e1z(z,Ce1)
!
!       ====================================================
!       Purpose: Compute complex exponential integral E1(z)
!       Input :  z   --- Argument of E1(z)
!       Output:  CE1 --- E1(z)
!       ====================================================
!
      implicit none
      real(wp) a0 , el , pi , x , xt
      complex(wp) Ce1 , cr , z , zc , zd , zdc
      integer k
      pi = 3.141592653589793d0
      el = 0.5772156649015328d0
      x = dble(z)
      a0 = abs(z)
!       Continued fraction converges slowly near negative real axis,
!       so use power series in a wedge around it until radius 40.0
      xt = -2*abs(dimag(z))
      if ( a0==0.0d0 ) then
         Ce1 = (1.0d+300,0.0d0)
      elseif ( a0<=5.0 .or. x<xt .and. a0<40.0 ) then
!          Power series
         Ce1 = (1.0d0,0.0d0)
         cr = (1.0d0,0.0d0)
         do k = 1 , 500
            cr = -cr*k*z/(k+1.0d0)**2
            Ce1 = Ce1 + cr
            if ( abs(cr)<=abs(Ce1)*1.0d-15 ) exit
         enddo
         if ( x<=0.0 .and. dimag(z)==0.0 ) then
!     Careful on the branch cut -- use the sign of the imaginary part
!     to get the right sign on the factor if pi.
            Ce1 = -el - log(-z) + z*Ce1 - sign(pi,dimag(z))            &
                & *(0.0d0,1.0d0)
         else
            Ce1 = -el - log(z) + z*Ce1
         endif
      else
!          Continued fraction https://dlmf.nist.gov/6.9
!
!                           1     1     1     2     2     3     3
!          E1 = exp(-z) * ----- ----- ----- ----- ----- ----- ----- ...
!                         Z +   1 +   Z +   1 +   Z +   1 +   Z +
         zc = 0d0
         zd = 1/z
         zdc = 1*zd
         zc = zc + zdc
         do k = 1 , 500
            zd = 1/(zd*k+1)
            zdc = (1*zd-1)*zdc
            zc = zc + zdc

            zd = 1/(zd*k+z)
            zdc = (z*zd-1)*zdc
            zc = zc + zdc

            if ( abs(zdc)<=abs(zc)*1.0d-15 .and. k>20 ) exit
         enddo
         Ce1 = exp(-z)*zc
         if ( x<=0.0 .and. dimag(z)==0.0 ) Ce1 = Ce1 - pi*(0.0d0,1.0d0)
      endif
      end

!       **********************************

      subroutine itjyb(x,Tj,Ty)
!
!       =======================================================
!       Purpose: Integrate Bessel functions J0(t) and Y0(t)
!                with respect to t from 0 to x ( x ≥ 0 )
!       Input :  x  --- Upper limit of the integral
!       Output:  TJ --- Integration of J0(t) from 0 to x
!                TY --- Integration of Y0(t) from 0 to x
!       =======================================================
!
      implicit none
      real(wp) f0 , g0 , pi , t , Tj , Ty , x , x1 , xt
      pi = 3.141592653589793d0
      if ( x==0.0d0 ) then
         Tj = 0.0d0
         Ty = 0.0d0
      elseif ( x<=4.0d0 ) then
         x1 = x/4.0d0
         t = x1*x1
         Tj = (((((((-.133718d-3*t+.2362211d-2)*t-.025791036d0)*t+      &
            & .197492634d0)*t-1.015860606d0)*t+3.199997842d0)           &
            & *t-5.333333161d0)*t+4.0d0)*x1
         Ty = ((((((((.13351d-4*t-.235002d-3)*t+.3034322d-2)*t-         &
            & .029600855d0)*t+.203380298d0)*t-.904755062d0)             &
            & *t+2.287317974d0)*t-2.567250468d0)*t+1.076611469d0)*x1
         Ty = 2.0d0/pi*log(x/2.0d0)*Tj - Ty
      elseif ( x<=8.0d0 ) then
         xt = x - .25d0*pi
         t = 16.0d0/(x*x)
         f0 = ((((((.1496119d-2*t-.739083d-2)*t+.016236617d0)*t-        &
            & .022007499d0)*t+.023644978d0)*t-.031280848d0)             &
            & *t+.124611058d0)*4.0d0/x
         g0 = (((((.1076103d-2*t-.5434851d-2)*t+.01242264d0)*t-         &
            & .018255209)*t+.023664841d0)*t-.049635633d0)               &
            & *t + .79784879d0
         Tj = 1.0d0 - (f0*cos(xt)-g0*sin(xt))/sqrt(x)
         Ty = -(f0*sin(xt)+g0*cos(xt))/sqrt(x)
      else
         t = 64.0d0/(x*x)
         xt = x - .25d0*pi
         f0 = (((((((-.268482d-4*t+.1270039d-3)*t-.2755037d-3)*t+       &
            & .3992825d-3)*t-.5366169d-3)*t+.10089872d-2)               &
            & *t-.40403539d-2)*t+.0623347304d0)*8.0d0/x
         g0 = ((((((-.226238d-4*t+.1107299d-3)*t-.2543955d-3)*t+        &
            & .4100676d-3)*t-.6740148d-3)*t+.17870944d-2)               &
            & *t-.01256424405d0)*t + .79788456d0
         Tj = 1.0d0 - (f0*cos(xt)-g0*sin(xt))/sqrt(x)
         Ty = -(f0*sin(xt)+g0*cos(xt))/sqrt(x)
      endif
      end


!       **********************************

      subroutine chgul(a,b,x,Hu,Id)
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
      implicit none
      real(wp) a , aa , b , Hu , r , r0 , ra , x
      integer Id , k , nm
      logical il1 , il2
      Id = -100
      aa = a - b + 1.0d0
      il1 = a==int(a) .and. a<=0.0
      il2 = aa==int(aa) .and. aa<=0.0
      nm = 0
      if ( il1 ) nm = abs(a)
      if ( il2 ) nm = abs(aa)
!       IL1: DLMF 13.2.7 with k=-s-a
!       IL2: DLMF 13.2.8
      if ( il1 .or. il2 ) then
         Hu = 1.0d0
         r = 1.0d0
         do k = 1 , nm
            r = -r*(a+k-1.0d0)*(a-b+k)/(k*x)
            Hu = Hu + r
         enddo
         Hu = x**(-a)*Hu
         Id = 10
      else
!       DLMF 13.7.3
         Hu = 1.0d0
         r = 1.0d0
         do k = 1 , 25
            r = -r*(a+k-1.0d0)*(a-b+k)/(k*x)
            ra = abs(r)
            if ( k>5 .and. ra>=r0 .or. ra<1.0d-15 ) exit
            r0 = ra
            Hu = Hu + r
         enddo
         Id = abs(log10(ra))
         Hu = x**(-a)*Hu
      endif
      end



!       **********************************

      subroutine gmn(m,n,c,x,Bk,Gf,Gd)
!
!       ===========================================================
!       Purpose: Compute gmn(-ic,ix) and its derivative for oblate
!                radial functions with a small argument
!       ===========================================================
!
      implicit none
      real(wp) Bk , c , eps , Gd , gd0 , gd1 , Gf , gf0 , gw ,  &
                     & x , xm
      integer ip , k , m , n , nm
      dimension Bk(200)
      eps = 1.0d-14
      ip = 1
      if ( n-m==2*int((n-m)/2) ) ip = 0
      nm = 25 + int(0.5*(n-m)+c)
      xm = (1.0d0+x*x)**(-0.5d0*m)
      gf0 = 0.0d0
      gw = 0.0d0
      do k = 1 , nm
         gf0 = gf0 + Bk(k)*x**(2.0*k-2.0)
         if ( abs((gf0-gw)/gf0)<eps .and. k>=10 ) exit
         gw = gf0
      enddo
      Gf = xm*gf0*x**(1-ip)
      gd1 = -m*x/(1.0d0+x*x)*Gf
      gd0 = 0.0d0
      do k = 1 , nm
         if ( ip==0 ) then
            gd0 = gd0 + (2.0d0*k-1.0)*Bk(k)*x**(2.0*k-2.0)
         else
            gd0 = gd0 + 2.0d0*k*Bk(k+1)*x**(2.0*k-1.0)
         endif
         if ( abs((gd0-gw)/gd0)<eps .and. k>=10 ) exit
         gw = gd0
      enddo
      Gd = gd1 + xm*gd0
      end



!       **********************************

      subroutine itjya(x,Tj,Ty)
!
!       ==========================================================
!       Purpose: Integrate Bessel functions J0(t) & Y0(t) with
!                respect to t from 0 to x
!       Input :  x  --- Upper limit of the integral ( x >= 0 )
!       Output:  TJ --- Integration of J0(t) from 0 to x
!                TY --- Integration of Y0(t) from 0 to x
!       =======================================================
!
      implicit none
      real(wp) a , a0 , a1 , af , bf , bg , el , eps , pi , r , &
                     & r2 , rc , rs , Tj , Ty , ty1 , ty2 , x , x2 , xp
      integer k
      dimension a(18)
      pi = 3.141592653589793d0
      el = .5772156649015329d0
      eps = 1.0d-12
      if ( x==0.0d0 ) then
         Tj = 0.0d0
         Ty = 0.0d0
      elseif ( x<=20.0d0 ) then
         x2 = x*x
         Tj = x
         r = x
         do k = 1 , 60
            r = -.25d0*r*(2*k-1.0d0)/(2*k+1.0d0)/(k*k)*x2
            Tj = Tj + r
            if ( abs(r)<abs(Tj)*eps ) exit
         enddo
         ty1 = (el+log(x/2.0d0))*Tj
         rs = 0.0d0
         ty2 = 1.0d0
         r = 1.0d0
         do k = 1 , 60
            r = -.25d0*r*(2*k-1.0d0)/(2*k+1.0d0)/(k*k)*x2
            rs = rs + 1.0d0/k
            r2 = r*(rs+1.0d0/(2.0d0*k+1.0d0))
            ty2 = ty2 + r2
            if ( abs(r2)<abs(ty2)*eps ) exit
         enddo
         Ty = (ty1-x*ty2)*2.0d0/pi
      else
         a0 = 1.0d0
         a1 = 5.0d0/8.0d0
         a(1) = a1
         do k = 1 , 16
            af = ((1.5d0*(k+.5d0)*(k+5.0d0/6.0d0)*a1-.5d0*(k+.5d0)      &
               & *(k+.5d0)*(k-.5d0)*a0))/(k+1.0d0)
            a(k+1) = af
            a0 = a1
            a1 = af
         enddo
         bf = 1.0d0
         r = 1.0d0
         do k = 1 , 8
            r = -r/(x*x)
            bf = bf + a(2*k)*r
         enddo
         bg = a(1)/x
         r = 1.0d0/x
         do k = 1 , 8
            r = -r/(x*x)
            bg = bg + a(2*k+1)*r
         enddo
         xp = x + .25d0*pi
         rc = sqrt(2.0d0/(pi*x))
         Tj = 1.0d0 - rc*(bf*cos(xp)+bg*sin(xp))
         Ty = rc*(bg*cos(xp)-bf*sin(xp))
      endif
      end

!       **********************************

      subroutine rcty(n,x,Nm,Ry,Dy)
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
      implicit none
      real(wp) Dy , rf0 , rf1 , rf2 , Ry , x
      integer k , n , Nm
      dimension Ry(0:n) , Dy(0:n)
      Nm = n
      if ( x<1.0d-60 ) then
         do k = 0 , n
            Ry(k) = -1.0d+300
            Dy(k) = 1.0d+300
         enddo
         Ry(0) = -1.0d0
         Dy(0) = 0.0d0
         return
      endif
      Ry(0) = -cos(x)
      Ry(1) = Ry(0)/x - sin(x)
      rf0 = Ry(0)
      rf1 = Ry(1)
      do k = 2 , n
         rf2 = (2.0d0*k-1.0d0)*rf1/x - rf0
         if ( abs(rf2)>1.0d+300 ) exit
         Ry(k) = rf2
         rf0 = rf1
         rf1 = rf2
      enddo
      Nm = k - 1
      Dy(0) = sin(x)
      do k = 1 , Nm
         Dy(k) = -k*Ry(k)/x + Ry(k-1)
      enddo
      end

!       **********************************

      subroutine lpni(n,x,Pn,Pd,Pl)
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
      implicit none
      integer j , k , n , n1
      real(wp) p0 , p1 , Pd , pf , Pl , Pn , r , x
      dimension Pn(0:n) , Pd(0:n) , Pl(0:n)
      Pn(0) = 1.0d0
      Pn(1) = x
      Pd(0) = 0.0d0
      Pd(1) = 1.0d0
      Pl(0) = x
      Pl(1) = 0.5d0*x*x
      p0 = 1.0d0
      p1 = x
      do k = 2 , n
         pf = (2.0d0*k-1.0d0)/k*x*p1 - (k-1.0d0)/k*p0
         Pn(k) = pf
         if ( abs(x)==1.0d0 ) then
            Pd(k) = 0.5d0*x**(k+1)*k*(k+1.0d0)
         else
            Pd(k) = k*(p1-x*pf)/(1.0d0-x*x)
         endif
         Pl(k) = (x*Pn(k)-Pn(k-1))/(k+1.0d0)
         p0 = p1
         p1 = pf
         if ( k/=2*int(k/2) ) then
            r = 1.0d0/(k+1.0d0)
            n1 = (k-1)/2
            do j = 1 , n1
               r = (0.5d0/j-1.0d0)*r
            enddo
            Pl(k) = Pl(k) + r
         endif
      enddo
      end

!       **********************************

      subroutine klvna(x,Ber,Bei,Ger,Gei,Der,Dei,Her,Hei)
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
      implicit none
      real(wp) Bei , Ber , cn0 , cp0 , cs , Dei , Der , el ,    &
                     & eps , fac , Gei , Ger , gs , Hei , Her , pi ,    &
                     & pn0 , pn1 , pp0 , pp1
      real(wp) qn0 , qn1 , qp0 , qp1 , r , r0 , r1 , rc , rs ,  &
                     & sn0 , sp0 , ss , x , x2 , x4 , xc1 , xc2 , xd ,  &
                     & xe1 , xe2
      real(wp) xt
      integer k , km , m
      pi = 3.141592653589793d0
      el = .5772156649015329d0
      eps = 1.0d-15
      if ( x==0.0d0 ) then
         Ber = 1.0d0
         Bei = 0.0d0
         Ger = 1.0d+300
         Gei = -0.25d0*pi
         Der = 0.0d0
         Dei = 0.0d0
         Her = -1.0d+300
         Hei = 0.0d0
         return
      endif
      x2 = 0.25d0*x*x
      x4 = x2*x2
      if ( abs(x)<10.0d0 ) then
         Ber = 1.0d0
         r = 1.0d0
         do m = 1 , 60
            r = -0.25d0*r/(m*m)/(2.0d0*m-1.0d0)**2*x4
            Ber = Ber + r
            if ( abs(r)<abs(Ber)*eps ) exit
         enddo
         Bei = x2
         r = x2
         do m = 1 , 60
            r = -0.25d0*r/(m*m)/(2.0d0*m+1.0d0)**2*x4
            Bei = Bei + r
            if ( abs(r)<abs(Bei)*eps ) exit
         enddo
         Ger = -(log(x/2.0d0)+el)*Ber + 0.25d0*pi*Bei
         r = 1.0d0
         gs = 0.0d0
         do m = 1 , 60
            r = -0.25d0*r/(m*m)/(2.0d0*m-1.0d0)**2*x4
            gs = gs + 1.0d0/(2.0d0*m-1.0d0) + 1.0d0/(2.0d0*m)
            Ger = Ger + r*gs
            if ( abs(r*gs)<abs(Ger)*eps ) exit
         enddo
         Gei = x2 - (log(x/2.0d0)+el)*Bei - 0.25d0*pi*Ber
         r = x2
         gs = 1.0d0
         do m = 1 , 60
            r = -0.25d0*r/(m*m)/(2.0d0*m+1.0d0)**2*x4
            gs = gs + 1.0d0/(2.0d0*m) + 1.0d0/(2.0d0*m+1.0d0)
            Gei = Gei + r*gs
            if ( abs(r*gs)<abs(Gei)*eps ) exit
         enddo
         Der = -0.25d0*x*x2
         r = Der
         do m = 1 , 60
            r = -0.25d0*r/m/(m+1.0d0)/(2.0d0*m+1.0d0)**2*x4
            Der = Der + r
            if ( abs(r)<abs(Der)*eps ) exit
         enddo
         Dei = 0.5d0*x
         r = Dei
         do m = 1 , 60
            r = -0.25d0*r/(m*m)/(2.d0*m-1.d0)/(2.d0*m+1.d0)*x4
            Dei = Dei + r
            if ( abs(r)<abs(Dei)*eps ) exit
         enddo
         r = -0.25d0*x*x2
         gs = 1.5d0
         Her = 1.5d0*r - Ber/x - (log(x/2.d0)+el)*Der + 0.25*pi*Dei
         do m = 1 , 60
            r = -0.25d0*r/m/(m+1.0d0)/(2.0d0*m+1.0d0)**2*x4
            gs = gs + 1.0d0/(2*m+1.0d0) + 1.0d0/(2*m+2.0d0)
            Her = Her + r*gs
            if ( abs(r*gs)<abs(Her)*eps ) exit
         enddo
         r = 0.5d0*x
         gs = 1.0d0
         Hei = 0.5d0*x - Bei/x - (log(x/2.d0)+el)*Dei - 0.25*pi*Der
         do m = 1 , 60
            r = -0.25d0*r/(m*m)/(2*m-1.0d0)/(2*m+1.0d0)*x4
            gs = gs + 1.0d0/(2.0d0*m) + 1.0d0/(2*m+1.0d0)
            Hei = Hei + r*gs
            if ( abs(r*gs)<abs(Hei)*eps ) return
         enddo
      else
         pp0 = 1.0d0
         pn0 = 1.0d0
         qp0 = 0.0d0
         qn0 = 0.0d0
         r0 = 1.0d0
         km = 18
         if ( abs(x)>=40.0 ) km = 10
         fac = 1.0d0
         do k = 1 , km
            fac = -fac
            xt = 0.25d0*k*pi - int(0.125d0*k)*2.0d0*pi
            cs = cos(xt)
            ss = sin(xt)
            r0 = 0.125d0*r0*(2.0d0*k-1.0d0)**2/k/x
            rc = r0*cs
            rs = r0*ss
            pp0 = pp0 + rc
            pn0 = pn0 + fac*rc
            qp0 = qp0 + rs
            qn0 = qn0 + fac*rs
         enddo
         xd = x/sqrt(2.0d0)
         xe1 = exp(xd)
         xe2 = exp(-xd)
         xc1 = 1.d0/sqrt(2.0d0*pi*x)
         xc2 = sqrt(.5d0*pi/x)
         cp0 = cos(xd+0.125d0*pi)
         cn0 = cos(xd-0.125d0*pi)
         sp0 = sin(xd+0.125d0*pi)
         sn0 = sin(xd-0.125d0*pi)
         Ger = xc2*xe2*(pn0*cp0-qn0*sp0)
         Gei = xc2*xe2*(-pn0*sp0-qn0*cp0)
         Ber = xc1*xe1*(pp0*cn0+qp0*sn0) - Gei/pi
         Bei = xc1*xe1*(pp0*sn0-qp0*cn0) + Ger/pi
         pp1 = 1.0d0
         pn1 = 1.0d0
         qp1 = 0.0d0
         qn1 = 0.0d0
         r1 = 1.0d0
         fac = 1.0d0
         do k = 1 , km
            fac = -fac
            xt = 0.25d0*k*pi - int(0.125d0*k)*2.0d0*pi
            cs = cos(xt)
            ss = sin(xt)
            r1 = 0.125d0*r1*(4.d0-(2.0d0*k-1.0d0)**2)/k/x
            rc = r1*cs
            rs = r1*ss
            pp1 = pp1 + fac*rc
            pn1 = pn1 + rc
            qp1 = qp1 + fac*rs
            qn1 = qn1 + rs
         enddo
         Her = xc2*xe2*(-pn1*cn0+qn1*sn0)
         Hei = xc2*xe2*(pn1*sn0+qn1*cn0)
         Der = xc1*xe1*(pp1*cp0+qp1*sp0) - Hei/pi
         Dei = xc1*xe1*(pp1*sp0-qp1*cp0) + Her/pi
      endif
      end

!       **********************************

      subroutine chgubi(a,b,x,Hu,Id)
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
      implicit none
      real(wp) a , a0 , a1 , a2 , b , da1 , da2 , db1 , db2 ,   &
                     & el , ga , ga1 , h0 , hm1 , hm2 , hm3 , hmax ,    &
                     & hmin , Hu , hu1
      real(wp) hu2 , hw , ps , r , rn , rn1 , s0 , s1 , s2 ,    &
                     & sa , sb , ua , ub , x
      integer Id , id1 , id2 , j , k , m , n
      Id = -100
      el = 0.5772156649015329d0
      n = abs(b-1)
      rn1 = 1.0d0
      rn = 1.0d0
      do j = 1 , n
         rn = rn*j
         if ( j==n-1 ) rn1 = rn
      enddo
      call psi_spec(a,ps)
      call gamma2(a,ga)
      if ( b>0.0 ) then
         a0 = a
         a1 = a - n
         a2 = a1
         call gamma2(a1,ga1)
         ua = (-1)**(n-1)/(rn*ga1)
         ub = rn1/ga*x**(-n)
      else
         a0 = a + n
         a1 = a0
         a2 = a
         call gamma2(a1,ga1)
         ua = (-1)**(n-1)/(rn*ga)*x**n
         ub = rn1/ga1
      endif
      hm1 = 1.0d0
      r = 1.0d0
      hmax = 0.0d0
      hmin = 1.0d+300
      h0 = 0d0
      do k = 1 , 150
         r = r*(a0+k-1.0d0)*x/((n+k)*k)
         hm1 = hm1 + r
         hu1 = abs(hm1)
         if ( hu1>hmax ) hmax = hu1
         if ( hu1<hmin ) hmin = hu1
         if ( abs(hm1-h0)<abs(hm1)*1.0d-15 ) exit
         h0 = hm1
      enddo
      da1 = log10(hmax)
      da2 = 0.0d0
      if ( hmin/=0.0 ) da2 = log10(hmin)
      Id = 15 - abs(da1-da2)
      hm1 = hm1*log(x)
      s0 = 0.0d0
      do m = 1 , n
         if ( b>=0.0 ) s0 = s0 - 1.0d0/m
         if ( b<0.0 ) s0 = s0 + (1.0d0-a)/(m*(a+m-1.0d0))
      enddo
      hm2 = ps + 2.0d0*el + s0
      r = 1.0d0
      hmax = 0.0d0
      hmin = 1.0d+300
      do k = 1 , 150
         s1 = 0.0d0
         s2 = 0.0d0
         if ( b>0.0 ) then
            do m = 1 , k
               s1 = s1 - (m+2.0d0*a-2.0d0)/(m*(m+a-1.0d0))
            enddo
            do m = 1 , n
               s2 = s2 + 1.0d0/(k+m)
            enddo
         else
            do m = 1 , k + n
               s1 = s1 + (1.0d0-a)/(m*(m+a-1.0d0))
            enddo
            do m = 1 , k
               s2 = s2 + 1.0d0/m
            enddo
         endif
         hw = 2.0d0*el + ps + s1 - s2
         r = r*(a0+k-1.0d0)*x/((n+k)*k)
         hm2 = hm2 + r*hw
         hu2 = abs(hm2)
         if ( hu2>hmax ) hmax = hu2
         if ( hu2<hmin ) hmin = hu2
         if ( abs((hm2-h0)/hm2)<1.0d-15 ) exit
         h0 = hm2
      enddo
      db1 = log10(hmax)
      db2 = 0.0d0
      if ( hmin/=0.0 ) db2 = log10(hmin)
      id1 = 15 - abs(db1-db2)
      if ( id1<Id ) Id = id1
      hm3 = 1.0d0
      if ( n==0 ) hm3 = 0.0d0
      r = 1.0d0
      do k = 1 , n - 1
         r = r*(a2+k-1.0d0)/((k-n)*k)*x
         hm3 = hm3 + r
      enddo
      sa = ua*(hm1+hm2)
      sb = ub*hm3
      Hu = sa + sb
      id2 = 0.0d0
      if ( sa/=0.0 ) id1 = int(log10(abs(sa)))
      if ( Hu/=0.0 ) id2 = int(log10(abs(Hu)))
      if ( sa*sb<0.0 ) Id = Id - abs(id1-id2)
      end



!       **********************************

      subroutine cyzo(Nt,Kf,Kc,Zo,Zv)
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
      implicit none
      real(wp) h , w , w0 , x , y
      integer i , it , j , Kc , Kf , nr , Nt
      complex(wp) z , zd , zf , zfd , zgd , Zo , zp , zq , Zv , zw
      dimension Zo(Nt) , Zv(Nt)
      x = 0.0d0
      y = 0.0d0
      h = 0.0d0
      if ( Kc==0 ) then
         x = -2.4d0
         y = 0.54d0
         h = 3.14d0
      elseif ( Kc==1 ) then
         x = 0.89
         y = 0.0
         h = -3.14
      endif
      if ( Kf==1 ) x = -0.503
      if ( Kf==2 ) x = 0.577
      z = dcmplx(x,y)
      w = 0.0d0
      do nr = 1 , Nt
         if ( nr/=1 ) z = Zo(nr-1) - h
         it = 0
 50      it = it + 1
         call cy01(Kf,z,zf,zd)
         zp = (1.0d0,0.0d0)
         do i = 1 , nr - 1
            zp = zp*(z-Zo(i))
         enddo
         zfd = zf/zp
         zq = (0.0d0,0.0d0)
         do i = 1 , nr - 1
            zw = (1.0d0,0.0d0)
            do j = 1 , nr - 1
               if ( j/=i ) zw = zw*(z-Zo(j))
            enddo
            zq = zq + zw
         enddo
         zgd = (zd-zq*zfd)/zp
         z = z - zfd/zgd
         w0 = w
         w = abs(z)
         if ( it<=50 .and. abs((w-w0)/w)>1.0d-12 ) goto 50
         Zo(nr) = z
      enddo
      do i = 1 , Nt
         z = Zo(i)
         if ( Kf==0 .or. Kf==2 ) then
            call cy01(1,z,zf,zd)
            Zv(i) = zf
         elseif ( Kf==1 ) then
            call cy01(0,z,zf,zd)
            Zv(i) = zf
         endif
      enddo
      end



!       **********************************

      subroutine klvnb(x,Ber,Bei,Ger,Gei,Der,Dei,Her,Hei)
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
      implicit none
      real(wp) Bei , Ber , csn , csp , Dei , Der , fxi , fxr ,  &
                     & Gei , Ger , Hei , Her , pi , pni , pnr , ppi ,   &
                     & ppr , ssn , ssp , t
      real(wp) t2 , tni , tnr , tpi , tpr , u , v , x , yc1 ,   &
                     & yc2 , yd , ye1 , ye2
      integer l
      pi = 3.141592653589793d0
      if ( x==0.0d0 ) then
         Ber = 1.0d0
         Bei = 0.0d0
         Ger = 1.0d+300
         Gei = -.25d0*pi
         Der = 0.0d0
         Dei = 0.0d0
         Her = -1.0d+300
         Hei = 0.0d0
      elseif ( x<8.0d0 ) then
         t = x/8.0d0
         t2 = t*t
         u = t2*t2
         Ber = ((((((-.901d-5*u+.122552d-2)*u-.08349609d0)*u+           &
             & 2.64191397d0)*u-32.36345652d0)*u+113.77777774d0)         &
             & *u-64.0d0)*u + 1.0d0
         Bei = t*t*((((((.11346d-3*u-.01103667d0)*u+.52185615d0)*u-     &
             & 10.56765779d0)*u+72.81777742d0)*u-113.77777774d0)        &
             & *u+16.0d0)
         Ger = ((((((-.2458d-4*u+.309699d-2)*u-.19636347d0)*u+          &
             & 5.65539121d0)*u-60.60977451d0)*u+171.36272133d0)         &
             & *u-59.05819744d0)*u - .57721566d0
         Ger = Ger - log(.5d0*x)*Ber + .25d0*pi*Bei
         Gei = t2*((((((.29532d-3*u-.02695875d0)*u+1.17509064d0)*u-     &
             & 21.30060904d0)*u+124.2356965d0)*u-142.91827687d0)        &
             & *u+6.76454936d0)
         Gei = Gei - log(.5d0*x)*Bei - .25d0*pi*Ber
         Der = x*t2*                                                    &
             & ((((((-.394d-5*u+.45957d-3)*u-.02609253d0)*u+.66047849d0)&
             & *u-6.0681481d0)*u+14.22222222d0)*u-4.0d0)
         Dei = x*((((((.4609d-4*u-.379386d-2)*u+.14677204d0)*u-         &
             & 2.31167514d0)*u+11.37777772d0)*u-10.66666666d0)*u+.5d0)
         Her = x*t2*((((((-.1075d-4*u+.116137d-2)*u-.06136358d0)*u+     &
             & 1.4138478d0)*u-11.36433272d0)*u+21.42034017d0)           &
             & *u-3.69113734d0)
         Her = Her - log(.5d0*x)*Der - Ber/x + .25d0*pi*Dei
         Hei = x*((((((.11997d-3*u-.926707d-2)*u+.33049424d0)*u-        &
             & 4.65950823d0)*u+19.41182758d0)*u-13.39858846d0)          &
             & *u+.21139217d0)
         Hei = Hei - log(.5d0*x)*Dei - Bei/x - .25d0*pi*Der
      else
         t = 8.0d0/x
         tnr = 0.0d0
         tni = 0.0d0
         do l = 1 , 2
            v = (-1)**l*t
            tpr = ((((.6d-6*v-.34d-5)*v-.252d-4)*v-.906d-4)             &
                & *v*v+.0110486d0)*v
            tpi = ((((.19d-5*v+.51d-5)*v*v-.901d-4)*v-.9765d-3)         &
                & *v-.0110485d0)*v - .3926991d0
            if ( l==1 ) then
               tnr = tpr
               tni = tpi
            endif
         enddo
         yd = x/sqrt(2.0d0)
         ye1 = exp(yd+tpr)
         ye2 = exp(-yd+tnr)
         yc1 = 1.0d0/sqrt(2.0d0*pi*x)
         yc2 = sqrt(pi/(2.0d0*x))
         csp = cos(yd+tpi)
         ssp = sin(yd+tpi)
         csn = cos(-yd+tni)
         ssn = sin(-yd+tni)
         Ger = yc2*ye2*csn
         Gei = yc2*ye2*ssn
         fxr = yc1*ye1*csp
         fxi = yc1*ye1*ssp
         Ber = fxr - Gei/pi
         Bei = fxi + Ger/pi
         pnr = 0.0d0
         pni = 0.0d0
         do l = 1 , 2
            v = (-1)**l*t
            ppr = (((((.16d-5*v+.117d-4)*v+.346d-4)*v+.5d-6)*v-.13813d-2&
                & )*v-.0625001d0)*v + .7071068d0
            ppi = (((((-.32d-5*v-.24d-5)*v+.338d-4)*v+.2452d-3)*v+      &
                & .13811d-2)*v-.1d-6)*v + .7071068d0
            if ( l==1 ) then
               pnr = ppr
               pni = ppi
            endif
         enddo
         Her = Gei*pni - Ger*pnr
         Hei = -(Gei*pnr+Ger*pni)
         Der = fxr*ppr - fxi*ppi - Hei/pi
         Dei = fxi*ppr + fxr*ppi + Her/pi
      endif
      end

!       **********************************

      subroutine rmn2so(m,n,c,x,Cv,Df,Kd,R2f,R2d)
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
      implicit none
      real(wp) bk , c , ck , ck1 , ck2 , Cv , Df , dn , eps ,   &
                     & gd , gf , h0 , pi , qs , qt , r1d , r1f , R2d ,  &
                     & R2f , sum
      real(wp) sw , x
      integer ip , j , Kd , m , n , nm
      dimension bk(200) , ck(200) , Df(200) , dn(200)
      if ( abs(Df(1))<=1.0d-280 ) then
         R2f = 1.0d+300
         R2d = 1.0d+300
         return
      endif
      eps = 1.0d-14
      pi = 3.141592653589793d0
      nm = 25 + int((n-m)/2+c)
      ip = 1
      if ( n-m==2*int((n-m)/2) ) ip = 0
      call sckb(m,n,c,Df,ck)
      call kmn(m,n,c,Cv,Kd,Df,dn,ck1,ck2)
      call qstar(m,n,c,ck,ck1,qs,qt)
      call cbk(m,n,c,Cv,qt,ck,bk)
      if ( x==0.0d0 ) then
         sum = 0.0d0
         sw = 0.0d0
         do j = 1 , nm
            sum = sum + ck(j)
            if ( abs(sum-sw)<abs(sum)*eps ) exit
            sw = sum
         enddo
         if ( ip==0 ) then
            r1f = sum/ck1
            R2f = -0.5d0*pi*qs*r1f
            R2d = qs*r1f + bk(1)
         elseif ( ip==1 ) then
            r1d = sum/ck1
            R2f = bk(1)
            R2d = -0.5d0*pi*qs*r1d
         endif
         return
      else
         call gmn(m,n,c,x,bk,gf,gd)
         call rmn1(m,n,c,x,Df,Kd,r1f,r1d)
         h0 = atan(x) - 0.5d0*pi
         R2f = qs*r1f*h0 + gf
         R2d = qs*(r1d*h0+r1f/(1.0d0+x*x)) + gd
      endif
      end



!       **********************************

      subroutine bjndd(n,x,Bj,Dj,Fj)
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
      implicit none
      real(wp) Bj , bs , Dj , f , f0 , f1 , Fj , x
      integer k , m , mt , n , nt
      dimension Bj(101) , Dj(101) , Fj(101)
      do nt = 1 , 900
         mt = int(0.5*log10(6.28*nt)-nt*log10(1.36*abs(x)/nt))
         if ( mt>20 ) exit
      enddo
      m = nt
      bs = 0.0d0
      f = 0.0d0
      f0 = 0.0d0
      f1 = 1.0d-35
      do k = m , 0 , -1
         f = 2.0d0*(k+1.0d0)*f1/x - f0
         if ( k<=n ) Bj(k+1) = f
         if ( k==2*int(k/2) ) bs = bs + 2.0d0*f
         f0 = f1
         f1 = f
      enddo
      do k = 0 , n
         Bj(k+1) = Bj(k+1)/(bs-f)
      enddo
      Dj(1) = -Bj(2)
      Fj(1) = -1.0d0*Bj(1) - Dj(1)/x
      do k = 1 , n
         Dj(k+1) = Bj(k) - k*Bj(k+1)/x
         Fj(k+1) = (k*k/(x*x)-1.0d0)*Bj(k+1) - Dj(k+1)/x
      enddo
      end

!       **********************************


      subroutine sphj(n,x,Nm,Sj,Dj)
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
      implicit none
      real(wp) cs , Dj , f , f0 , f1 , sa , sb , Sj , x
      integer k , m , n , Nm
      dimension Sj(0:n) , Dj(0:n)
      Nm = n
      if ( abs(x)<1.0d-100 ) then
         do k = 0 , n
            Sj(k) = 0.0d0
            Dj(k) = 0.0d0
         enddo
         Sj(0) = 1.0d0
         if ( n>0 ) Dj(1) = .3333333333333333d0
         return
      endif
      Sj(0) = sin(x)/x
      Dj(0) = (cos(x)-sin(x)/x)/x
      if ( n<1 ) return
      Sj(1) = (Sj(0)-cos(x))/x
      if ( n>=2 ) then
         sa = Sj(0)
         sb = Sj(1)
         m = msta1(x,200)
         if ( m<n ) then
            Nm = m
         else
            m = msta2(x,n,15)
         endif
         f = 0.0d0
         f0 = 0.0d0
         f1 = 1.0d0 - 100
         do k = m , 0 , -1
            f = (2.0d0*k+3.0d0)*f1/x - f0
            if ( k<=Nm ) Sj(k) = f
            f0 = f1
            f1 = f
         enddo
         cs = 0.0d0
         if ( abs(sa)>abs(sb) ) cs = sa/f
         if ( abs(sa)<=abs(sb) ) cs = sb/f0
         do k = 0 , Nm
            Sj(k) = cs*Sj(k)
         enddo
      endif
      do k = 1 , Nm
         Dj(k) = Sj(k-1) - (k+1.0d0)*Sj(k)/x
      enddo
      end



!       **********************************

      subroutine othpl(Kf,n,x,Pl,Dpl)
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
      implicit none
      real(wp) a , b , c , Dpl , dy0 , dy1 , dyn , Pl , x , y0 ,&
                     & y1 , yn
      integer k , Kf , n
      dimension Pl(0:n) , Dpl(0:n)
      a = 2.0d0
      b = 0.0d0
      c = 1.0d0
      y0 = 1.0d0
      y1 = 2.0d0*x
      dy0 = 0.0d0
      dy1 = 2.0d0
      Pl(0) = 1.0d0
      Pl(1) = 2.0d0*x
      Dpl(0) = 0.0d0
      Dpl(1) = 2.0d0
      if ( Kf==1 ) then
         y1 = x
         dy1 = 1.0d0
         Pl(1) = x
         Dpl(1) = 1.0d0
      elseif ( Kf==3 ) then
         y1 = 1.0d0 - x
         dy1 = -1.0d0
         Pl(1) = 1.0d0 - x
         Dpl(1) = -1.0d0
      endif
      do k = 2 , n
         if ( Kf==3 ) then
            a = -1.0d0/k
            b = 2.0d0 + a
            c = 1.0d0 + a
         elseif ( Kf==4 ) then
            c = 2.0d0*(k-1.0d0)
         endif
         yn = (a*x+b)*y1 - c*y0
         dyn = a*y1 + (a*x+b)*dy1 - c*dy0
         Pl(k) = yn
         Dpl(k) = dyn
         y0 = y1
         y1 = yn
         dy0 = dy1
         dy1 = dyn
      enddo
      end

!       **********************************

      subroutine klvnzo(Nt,Kd,Zo)
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
      implicit none
      real(wp) bei , ber , ddi , ddr , dei , der , gdi , gdr ,  &
                     & gei , ger , hei , her , rt , rt0 , Zo
      integer Kd , m , Nt
      dimension Zo(Nt) , rt0(8)
      rt0(1) = 2.84891
      rt0(2) = 5.02622
      rt0(3) = 1.71854
      rt0(4) = 3.91467
      rt0(5) = 6.03871
      rt0(6) = 3.77268
      rt0(7) = 2.66584
      rt0(8) = 4.93181
      rt = rt0(Kd)
      do m = 1 , Nt
         call klvna(rt,ber,bei,ger,gei,der,dei,her,hei)
         if ( Kd==1 ) then
            rt = rt - ber/der
         elseif ( Kd==2 ) then
            rt = rt - bei/dei
         elseif ( Kd==3 ) then
            rt = rt - ger/her
         elseif ( Kd==4 ) then
            rt = rt - gei/hei
         elseif ( Kd==5 ) then
            ddr = -bei - der/rt
            rt = rt - der/ddr
         elseif ( Kd==6 ) then
            ddi = ber - dei/rt
            rt = rt - dei/ddi
         elseif ( Kd==7 ) then
            gdr = -gei - her/rt
            rt = rt - her/gdr
         else
            gdi = ger - hei/rt
            rt = rt - hei/gdi
         endif
         if ( abs(rt-rt0(Kd))>5.0d-10 ) then
            rt0(Kd) = rt
            cycle
         endif
         Zo(m) = rt
         rt = rt + 4.44d0
      enddo
      end



!       **********************************

      subroutine rswfo(m,n,c,x,Cv,Kf,R1f,R1d,R2f,R2d)
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
      implicit none
      real(wp) c , Cv , df , R1d , R1f , R2d , R2f , x
      integer id , kd , Kf , m , n
      dimension df(200)
      kd = -1
      call sdmn(m,n,c,Cv,kd,df)
      if ( Kf/=2 ) call rmn1(m,n,c,x,df,kd,R1f,R1d)
      if ( Kf>1 ) then
         id = 10
         if ( x>1.0d-8 ) call rmn2l(m,n,c,x,df,kd,R2f,R2d,id)
         if ( id>-1 ) call rmn2so(m,n,c,x,Cv,df,kd,R2f,R2d)
      endif
      end



!       **********************************

      subroutine ch12n(n,z,Nm,Chf1,Chd1,Chf2,Chd2)
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
      implicit none
      complex(wp) cbi , cbj , cbk , cby , cdi , cdj , cdk , cdy , cf1 ,  &
               & cfac , Chd1 , Chd2 , Chf1 , Chf2 , ci , z , zi
      integer k , n , Nm
      real(wp) pi
      dimension cbj(0:250) , cdj(0:250) , cby(0:250) , cdy(0:250) ,     &
              & cbi(0:250) , cdi(0:250) , cbk(0:250) , cdk(0:250)
      dimension Chf1(0:n) , Chd1(0:n) , Chf2(0:n) , Chd2(0:n)
      ci = (0.0d0,1.0d0)
      pi = 3.141592653589793d0
      if ( dimag(z)<0.0d0 ) then
         call cjynb(n,z,Nm,cbj,cdj,cby,cdy)
         do k = 0 , Nm
            Chf1(k) = cbj(k) + ci*cby(k)
            Chd1(k) = cdj(k) + ci*cdy(k)
         enddo
         zi = ci*z
         call ciknb(n,zi,Nm,cbi,cdi,cbk,cdk)
         cfac = -2.0d0/(pi*ci)
         do k = 0 , Nm
            Chf2(k) = cfac*cbk(k)
            Chd2(k) = cfac*ci*cdk(k)
            cfac = cfac*ci
         enddo
      elseif ( dimag(z)>0.0d0 ) then
         zi = -ci*z
         call ciknb(n,zi,Nm,cbi,cdi,cbk,cdk)
         cf1 = -ci
         cfac = 2.0d0/(pi*ci)
         do k = 0 , Nm
            Chf1(k) = cfac*cbk(k)
            Chd1(k) = -cfac*ci*cdk(k)
            cfac = cfac*cf1
         enddo
         call cjynb(n,z,Nm,cbj,cdj,cby,cdy)
         do k = 0 , Nm
            Chf2(k) = cbj(k) - ci*cby(k)
            Chd2(k) = cdj(k) - ci*cdy(k)
         enddo
      else
         call cjynb(n,z,Nm,cbj,cdj,cby,cdy)
         do k = 0 , Nm
            Chf1(k) = cbj(k) + ci*cby(k)
            Chd1(k) = cdj(k) + ci*cdy(k)
            Chf2(k) = cbj(k) - ci*cby(k)
            Chd2(k) = cdj(k) - ci*cdy(k)
         enddo
      endif
      end



!       **********************************

      subroutine jyzo(n,Nt,Rj0,Rj1,Ry0,Ry1)
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
      implicit none
      real(wp) bjn , byn , djn , dyn , fjn , fyn , pi , Rj0 ,   &
                     & Rj1 , Ry0 , Ry1 , x , x0 , xguess
      integer l , n , Nt
      dimension Rj0(Nt) , Rj1(Nt) , Ry0(Nt) , Ry1(Nt)
      pi = 3.141592653589793d0
!       -- Newton method for j_{N,L}
!       1) initial guess for j_{N,1}
      if ( n<=20 ) then
         x = 2.82141 + 1.15859*n
      else
!          Abr & Stg (9.5.14)
         x = n + 1.85576*n**0.33333 + 1.03315/n**0.33333
      endif
      l = 0
!       2) iterate
      xguess = x
 100  x0 = x
      call jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn)
      x = x - bjn/djn
      if ( x-x0<-1 ) x = x0 - 1
      if ( x-x0>1 ) x = x0 + 1
      if ( abs(x-x0)>1.0d-11 ) goto 100
!       3) initial guess for j_{N,L+1}
      if ( l>=1 ) then
         if ( x<=Rj0(l)+0.5 ) then
            x = xguess + pi
            xguess = x
            goto 100
         endif
      endif
      l = l + 1
      Rj0(l) = x
!       XXX: should have a better initial guess for large N ~> 100 here
      x = x + pi + max((0.0972d0+0.0679*n-0.000354*n**2)/l,0d0)
      if ( l<Nt ) goto 100
!       -- Newton method for j_{N,L}'
      if ( n<=20 ) then
         x = 0.961587 + 1.07703*n
      else
         x = n + 0.80861*n**0.33333 + 0.07249/n**0.33333
      endif
      if ( n==0 ) x = 3.8317
      l = 0
      xguess = x
 200  x0 = x
      call jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn)
      x = x - djn/fjn
      if ( x-x0<-1 ) x = x0 - 1
      if ( x-x0>1 ) x = x0 + 1
      if ( abs(x-x0)>1.0d-11 ) goto 200
      if ( l>=1 ) then
         if ( x<=Rj1(l)+0.5 ) then
            x = xguess + pi
            xguess = x
            goto 200
         endif
      endif
      l = l + 1
      Rj1(l) = x
!       XXX: should have a better initial guess for large N ~> 100 here
      x = x + pi + max((0.4955d0+0.0915*n-0.000435*n**2)/l,0d0)
      if ( l<Nt ) goto 200
!       -- Newton method for y_{N,L}
      if ( n<=20 ) then
         x = 1.19477 + 1.08933*n
      else
         x = n + 0.93158*n**0.33333 + 0.26035/n**0.33333
      endif
      l = 0
      xguess = x
 300  x0 = x
      call jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn)
      x = x - byn/dyn
      if ( x-x0<-1 ) x = x0 - 1
      if ( x-x0>1 ) x = x0 + 1
      if ( abs(x-x0)>1.0d-11 ) goto 300
      if ( l>=1 ) then
         if ( x<=Ry0(l)+0.5 ) then
            x = xguess + pi
            xguess = x
            goto 300
         endif
      endif
      l = l + 1
      Ry0(l) = x
!       XXX: should have a better initial guess for large N ~> 100 here
      x = x + pi + max((0.312d0+0.0852*n-0.000403*n**2)/l,0d0)
      if ( l<Nt ) goto 300
!       -- Newton method for y_{N,L}'
      if ( n<=20 ) then
         x = 2.67257 + 1.16099*n
      else
         x = n + 1.8211*n**0.33333 + 0.94001/n**0.33333
      endif
      l = 0
      xguess = x
 400  x0 = x
      call jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn)
      x = x - dyn/fyn
      if ( abs(x-x0)>1.0d-11 ) goto 400
      if ( l>=1 ) then
         if ( x<=Ry1(l)+0.5 ) then
            x = xguess + pi
            xguess = x
            goto 400
         endif
      endif
      l = l + 1
      Ry1(l) = x
!       XXX: should have a better initial guess for large N ~> 100 here
      x = x + pi + max((0.197d0+0.0643*n-0.000286*n**2)/l,0d0)
      if ( l<Nt ) goto 400
      end



!       **********************************

      subroutine ikv(v,x,Vm,Bi,Di,Bk,Dk)
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
      implicit none
      real(wp) a1 , a2 , Bi , bi0 , Bk , bk0 , bk1 , bk2 , ca , &
                     & cb , cs , ct , Di , Dk , f , f1 , f2 , gan ,     &
                     & gap , pi
      real(wp) piv , r , r1 , r2 , sum , v , v0 , v0n , v0p ,   &
                     & Vm , vt , w0 , wa , ww , x , x2
      integer k , k0 , m , n
      dimension Bi(0:*) , Di(0:*) , Bk(0:*) , Dk(0:*)
      pi = 3.141592653589793d0
      x2 = x*x
      n = int(v)
      v0 = v - n
      if ( n==0 ) n = 1
      if ( x<1.0d-100 ) then
         do k = 0 , n
            Bi(k) = 0.0d0
            Di(k) = 0.0d0
            Bk(k) = -1.0d+300
            Dk(k) = 1.0d+300
         enddo
         if ( v==0.0 ) then
            Bi(0) = 1.0d0
            Di(1) = 0.5d0
         endif
         Vm = v
         return
      endif
      piv = pi*v0
      vt = 4.0d0*v0*v0
      if ( v0==0.0d0 ) then
         a1 = 1.0d0
      else
         v0p = 1.0d0 + v0
         call gamma2(v0p,gap)
         a1 = (0.5d0*x)**v0/gap
      endif
      k0 = 14
      if ( x>=35.0 ) k0 = 10
      if ( x>=50.0 ) k0 = 8
      if ( x<=18.0 ) then
         bi0 = 1.0d0
         r = 1.0d0
         do k = 1 , 30
            r = 0.25d0*r*x2/(k*(k+v0))
            bi0 = bi0 + r
            if ( abs(r/bi0)<1.0d-15 ) exit
         enddo
         bi0 = bi0*a1
      else
         ca = exp(x)/sqrt(2.0d0*pi*x)
         sum = 1.0d0
         r = 1.0d0
         do k = 1 , k0
            r = -0.125d0*r*(vt-(2.0d0*k-1.0d0)**2.0)/(k*x)
            sum = sum + r
         enddo
         bi0 = ca*sum
      endif
      m = msta1(x,200)
      if ( m<n ) then
         n = m
      else
         m = msta2(x,n,15)
      endif
      f = 0.0d0
      f2 = 0.0d0
      f1 = 1.0d-100
      ww = 0.0d0
      do k = m , 0 , -1
         f = 2.0d0*(v0+k+1.0d0)/x*f1 + f2
         if ( k<=n ) Bi(k) = f
         f2 = f1
         f1 = f
      enddo
      cs = bi0/f
      do k = 0 , n
         Bi(k) = cs*Bi(k)
      enddo
      Di(0) = v0/x*Bi(0) + Bi(1)
      do k = 1 , n
         Di(k) = -(k+v0)/x*Bi(k) + Bi(k-1)
      enddo
      if ( x>9.0d0 ) then
         cb = exp(-x)*sqrt(0.5d0*pi/x)
         sum = 1.0d0
         r = 1.0d0
         do k = 1 , k0
            r = 0.125d0*r*(vt-(2.0*k-1.0)**2.0)/(k*x)
            sum = sum + r
         enddo
         bk0 = cb*sum
      elseif ( v0==0.0d0 ) then
         ct = -log(0.5d0*x) - 0.5772156649015329d0
         cs = 0.0d0
         w0 = 0.0d0
         r = 1.0d0
         do k = 1 , 50
            w0 = w0 + 1.0d0/k
            r = 0.25d0*r/(k*k)*x2
            cs = cs + r*(w0+ct)
            wa = abs(cs)
            if ( abs((wa-ww)/wa)<1.0d-15 ) exit
            ww = wa
         enddo
         bk0 = ct + cs
      else
         v0n = 1.0d0 - v0
         call gamma2(v0n,gan)
         a2 = 1.0d0/(gan*(0.5d0*x)**v0)
         a1 = (0.5d0*x)**v0/gap
         sum = a2 - a1
         r1 = 1.0d0
         r2 = 1.0d0
         do k = 1 , 120
            r1 = 0.25d0*r1*x2/(k*(k-v0))
            r2 = 0.25d0*r2*x2/(k*(k+v0))
            sum = sum + a2*r1 - a1*r2
            wa = abs(sum)
            if ( abs((wa-ww)/wa)<1.0d-15 ) exit
            ww = wa
         enddo
         bk0 = 0.5d0*pi*sum/sin(piv)
      endif
      bk1 = (1.0d0/x-Bi(1)*bk0)/Bi(0)
      Bk(0) = bk0
      Bk(1) = bk1
      do k = 2 , n
         bk2 = 2.0d0*(v0+k-1.0d0)/x*bk1 + bk0
         Bk(k) = bk2
         bk0 = bk1
         bk1 = bk2
      enddo
      Dk(0) = v0/x*Bk(0) - Bk(1)
      do k = 1 , n
         Dk(k) = -(k+v0)/x*Bk(k) - Bk(k-1)
      enddo
      Vm = n + v0
      end



!       **********************************

      subroutine sdmn(m,n,c,Cv,Kd,Df)
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
      implicit none
      real(wp) a , c , cs , Cv , d , d2k , Df , dk0 , dk1 ,     &
                     & dk2 , f , f0 , f1 , f2 , fl , fs , g , r1 , r3 , &
                     & r4
      real(wp) s0 , su1 , su2 , sw
      integer i , ip , j , k , k1 , kb , Kd , m , n , nm
      dimension a(200) , d(200) , g(200) , Df(200)
      nm = 25 + int(0.5*(n-m)+c)
      if ( c<1.0d-10 ) then
         do i = 1 , nm
            Df(i) = 0d0
         enddo
         Df((n-m)/2+1) = 1.0d0
         return
      endif
      cs = c*c*Kd
      ip = 1
      k = 0
      if ( n-m==2*int((n-m)/2) ) ip = 0
      do i = 1 , nm + 2
         if ( ip==0 ) k = 2*(i-1)
         if ( ip==1 ) k = 2*i - 1
         dk0 = m + k
         dk1 = m + k + 1
         dk2 = 2*(m+k)
         d2k = 2*m + k
         a(i) = (d2k+2.0)*(d2k+1.0)/((dk2+3.0)*(dk2+5.0))*cs
         d(i) = dk0*dk1 + (2.0*dk0*dk1-2.0*m*m-1.0)                     &
              & /((dk2-1.0)*(dk2+3.0))*cs
         g(i) = k*(k-1.0)/((dk2-3.0)*(dk2-1.0))*cs
      enddo
      fs = 1.0d0
      f1 = 0.0d0
      f0 = 1.0d-100
      kb = 0
      Df(nm+1) = 0.0d0
      fl = 0.0d0
      do k = nm , 1 , -1
         f = -((d(k+1)-Cv)*f0+a(k+1)*f1)/g(k+1)
         if ( abs(f)>abs(Df(k+1)) ) then
            Df(k) = f
            f1 = f0
            f0 = f
            if ( abs(f)>1.0d+100 ) then
               do k1 = k , nm
                  Df(k1) = Df(k1)*1.0d-100
               enddo
               f1 = f1*1.0d-100
               f0 = f0*1.0d-100
            endif
         else
            kb = k
            fl = Df(k+1)
            f1 = 1.0d-100
            f2 = -(d(1)-Cv)/a(1)*f1
            Df(1) = f1
            if ( kb==1 ) then
               fs = f2
            elseif ( kb==2 ) then
               Df(2) = f2
               fs = -((d(2)-Cv)*f2+g(2)*f1)/a(2)
            else
               Df(2) = f2
               do j = 3 , kb + 1
                  f = -((d(j-1)-Cv)*f2+g(j-1)*f1)/a(j-1)
                  if ( j<=kb ) Df(j) = f
                  if ( abs(f)>1.0d+100 ) then
                     do k1 = 1 , j
                        Df(k1) = Df(k1)*1.0d-100
                     enddo
                     f = f*1.0d-100
                     f2 = f2*1.0d-100
                  endif
                  f1 = f2
                  f2 = f
               enddo
               fs = f
            endif
            exit
         endif
      enddo
      su1 = 0.0d0
      r1 = 1.0d0
      do j = m + ip + 1 , 2*(m+ip)
         r1 = r1*j
      enddo
      su1 = Df(1)*r1
      do k = 2 , kb
         r1 = -r1*(k+m+ip-1.5d0)/(k-1.0d0)
         su1 = su1 + r1*Df(k)
      enddo
      su2 = 0.0d0
      sw = 0.0d0
      do k = kb + 1 , nm
         if ( k/=1 ) r1 = -r1*(k+m+ip-1.5d0)/(k-1.0d0)
         su2 = su2 + r1*Df(k)
         if ( abs(sw-su2)<abs(su2)*1.0d-14 ) exit
         sw = su2
      enddo
      r3 = 1.0d0
      do j = 1 , (m+n+ip)/2
         r3 = r3*(j+0.5d0*(n+m+ip))
      enddo
      r4 = 1.0d0
      do j = 1 , (n-m-ip)/2
         r4 = -4.0d0*r4*j
      enddo
      s0 = r3/(fl*(su1/fs)+su2)/r4
      do k = 1 , kb
         Df(k) = fl/fs*s0*Df(k)
      enddo
      do k = kb + 1 , nm
         Df(k) = s0*Df(k)
      enddo
      end




!       **********************************

      subroutine ajyik(x,Vj1,Vj2,Vy1,Vy2,Vi1,Vi2,Vk1,Vk2)
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
      implicit none
      real(wp) a0 , b0 , c0 , ck , gn , gn1 , gn2 , gp1 , gp2 , &
                     & pi , pv1 , pv2 , px , qx , r , rp , rp2 , rq ,   &
                     & sk , sum
      real(wp) uj1 , uj2 , uu0 , Vi1 , Vi2 , vil , Vj1 , Vj2 ,  &
                     & vjl , Vk1 , Vk2 , vl , vsl , vv , vv0 , Vy1 ,    &
                     & Vy2 , x , x2 , xk
      integer k , k0 , l
      if ( x==0.0d0 ) then
         Vj1 = 0.0d0
         Vj2 = 0.0d0
         Vy1 = -1.0d+300
         Vy2 = 1.0d+300
         Vi1 = 0.0d0
         Vi2 = 0.0d0
         Vk1 = -1.0d+300
         Vk2 = -1.0d+300
         return
      endif
      pi = 3.141592653589793d0
      rp2 = .63661977236758d0
      gp1 = .892979511569249d0
      gp2 = .902745292950934d0
      gn1 = 1.3541179394264d0
      gn2 = 2.678938534707747d0
      vv0 = 0.444444444444444d0
      uu0 = 1.1547005383793d0
      x2 = x*x
      k0 = 12
      if ( x>=35.0 ) k0 = 10
      if ( x>=50.0 ) k0 = 8
      if ( x<=12.0 ) then
         do l = 1 , 2
            vl = l/3.0d0
            vjl = 1.0d0
            r = 1.0d0
            do k = 1 , 40
               r = -0.25d0*r*x2/(k*(k+vl))
               vjl = vjl + r
               if ( abs(r)<1.0d-15 ) exit
            enddo
            a0 = (0.5d0*x)**vl
            if ( l==1 ) Vj1 = a0/gp1*vjl
            if ( l==2 ) Vj2 = a0/gp2*vjl
         enddo
      else
         do l = 1 , 2
            vv = vv0*l*l
            px = 1.0d0
            rp = 1.0d0
            do k = 1 , k0
               rp = -0.78125d-2*rp*(vv-(4.0*k-3.0)**2.0)                &
                  & *(vv-(4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*x2)
               px = px + rp
            enddo
            qx = 1.0d0
            rq = 1.0d0
            do k = 1 , k0
               rq = -0.78125d-2*rq*(vv-(4.0*k-1.0)**2.0)                &
                  & *(vv-(4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*x2)
               qx = qx + rq
            enddo
            qx = 0.125d0*(vv-1.0)*qx/x
            xk = x - (0.5d0*l/3.0d0+0.25d0)*pi
            a0 = sqrt(rp2/x)
            ck = cos(xk)
            sk = sin(xk)
            if ( l==1 ) then
               Vj1 = a0*(px*ck-qx*sk)
               Vy1 = a0*(px*sk+qx*ck)
            elseif ( l==2 ) then
               Vj2 = a0*(px*ck-qx*sk)
               Vy2 = a0*(px*sk+qx*ck)
            endif
         enddo
      endif
      if ( x<=12.0d0 ) then
         uj1 = 0.0d0
         uj2 = 0.0d0
         do l = 1 , 2
            vl = l/3.0d0
            vjl = 1.0d0
            r = 1.0d0
            do k = 1 , 40
               r = -0.25d0*r*x2/(k*(k-vl))
               vjl = vjl + r
               if ( abs(r)<1.0d-15 ) exit
            enddo
            b0 = (2.0d0/x)**vl
            if ( l==1 ) uj1 = b0*vjl/gn1
            if ( l==2 ) uj2 = b0*vjl/gn2
         enddo
         pv1 = pi/3.0d0
         pv2 = pi/1.5d0
         Vy1 = uu0*(Vj1*cos(pv1)-uj1)
         Vy2 = uu0*(Vj2*cos(pv2)-uj2)
      endif
      if ( x<=18.0 ) then
         do l = 1 , 2
            vl = l/3.0d0
            vil = 1.0d0
            r = 1.0d0
            do k = 1 , 40
               r = 0.25d0*r*x2/(k*(k+vl))
               vil = vil + r
               if ( abs(r)<1.0d-15 ) exit
            enddo
            a0 = (0.5d0*x)**vl
            if ( l==1 ) Vi1 = a0/gp1*vil
            if ( l==2 ) Vi2 = a0/gp2*vil
         enddo
      else
         c0 = exp(x)/sqrt(2.0d0*pi*x)
         do l = 1 , 2
            vv = vv0*l*l
            vsl = 1.0d0
            r = 1.0d0
            do k = 1 , k0
               r = -0.125d0*r*(vv-(2.0d0*k-1.0d0)**2.0)/(k*x)
               vsl = vsl + r
            enddo
            if ( l==1 ) Vi1 = c0*vsl
            if ( l==2 ) Vi2 = c0*vsl
         enddo
      endif
      if ( x<=9.0d0 ) then
         gn = 0.0d0
         do l = 1 , 2
            vl = l/3.0d0
            if ( l==1 ) gn = gn1
            if ( l==2 ) gn = gn2
            a0 = (2.0d0/x)**vl/gn
            sum = 1.0d0
            r = 1.0d0
            do k = 1 , 60
               r = 0.25d0*r*x2/(k*(k-vl))
               sum = sum + r
               if ( abs(r)<1.0d-15 ) exit
            enddo
            if ( l==1 ) Vk1 = 0.5d0*uu0*pi*(sum*a0-Vi1)
            if ( l==2 ) Vk2 = 0.5d0*uu0*pi*(sum*a0-Vi2)
         enddo
      else
         c0 = exp(-x)*sqrt(0.5d0*pi/x)
         do l = 1 , 2
            vv = vv0*l*l
            sum = 1.0d0
            r = 1.0d0
            do k = 1 , k0
               r = 0.125d0*r*(vv-(2.0*k-1.0)**2.0)/(k*x)
               sum = sum + r
            enddo
            if ( l==1 ) Vk1 = c0*sum
            if ( l==2 ) Vk2 = c0*sum
         enddo
      endif
      end



!       **********************************

      subroutine cikvb(v,z,Vm,Cbi,Cdi,Cbk,Cdk)
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
      implicit none
      real(wp) a0 , gan , gap , pi , piv , v , v0 , v0n , v0p , &
                     & Vm , vt , w0
      complex(wp) ca , ca1 , ca2 , cb , Cbi , cbi0 , Cbk , cbk0 , Cdi ,  &
               & Cdk , cf , cf1 , cf2 , ci , ci0 , ckk , cp , cr , cr1 ,&
               & cr2
      complex(wp) cs , csu , ct , cvk , z , z1 , z2
      integer k , k0 , m , n
      dimension Cbi(0:*) , Cdi(0:*) , Cbk(0:*) , Cdk(0:*)
      z1 = z
      z2 = z*z
      a0 = abs(z)
      pi = 3.141592653589793d0
      ci = (0.0d0,1.0d0)
      n = int(v)
      v0 = v - n
      piv = pi*v0
      vt = 4.0d0*v0*v0
      if ( n==0 ) n = 1
      if ( a0<1.0d-100 ) then
         do k = 0 , n
            Cbi(k) = 0.0d0
            Cdi(k) = 0.0d0
            Cbk(k) = -1.0d+300
            Cdk(k) = 1.0d+300
         enddo
         if ( v0==0.0 ) then
            Cbi(0) = (1.0d0,0.0d0)
            Cdi(1) = (0.5d0,0.0d0)
         endif
         Vm = v
         return
      endif
      k0 = 14
      if ( a0>=35.0 ) k0 = 10
      if ( a0>=50.0 ) k0 = 8
      if ( dble(z)<0.0 ) z1 = -z
      if ( a0<18.0 ) then
         if ( v0==0.0 ) then
            ca1 = (1.0d0,0.0d0)
         else
            v0p = 1.0d0 + v0
            call gamma2(v0p,gap)
            ca1 = (0.5d0*z1)**v0/gap
         endif
         ci0 = (1.0d0,0.0d0)
         cr = (1.0d0,0.0d0)
         do k = 1 , 50
            cr = 0.25d0*cr*z2/(k*(k+v0))
            ci0 = ci0 + cr
            if ( abs(cr/ci0)<1.0d-15 ) exit
         enddo
         cbi0 = ci0*ca1
      else
         ca = exp(z1)/sqrt(2.0d0*pi*z1)
         cs = (1.0d0,0.0d0)
         cr = (1.0d0,0.0d0)
         do k = 1 , k0
            cr = -0.125d0*cr*(vt-(2.0d0*k-1.0d0)**2.0)/(k*z1)
            cs = cs + cr
         enddo
         cbi0 = ca*cs
      endif
      m = msta1(a0,200)
      if ( m<n ) then
         n = m
      else
         m = msta2(a0,n,15)
      endif
      cf2 = (0.0d0,0.0d0)
      cf1 = (1.0d-100,0.0d0)
      do k = m , 0 , -1
         cf = 2.0d0*(v0+k+1.0d0)/z1*cf1 + cf2
         if ( k<=n ) Cbi(k) = cf
         cf2 = cf1
         cf1 = cf
      enddo
      cs = cbi0/cf
      do k = 0 , n
         Cbi(k) = cs*Cbi(k)
      enddo
      if ( a0>9.0 ) then
         cb = exp(-z1)*sqrt(0.5d0*pi/z1)
         cs = (1.0d0,0.0d0)
         cr = (1.0d0,0.0d0)
         do k = 1 , k0
            cr = 0.125d0*cr*(vt-(2.0d0*k-1.0d0)**2.0)/(k*z1)
            cs = cs + cr
         enddo
         cbk0 = cb*cs
      elseif ( v0==0.0 ) then
         ct = -log(0.5d0*z1) - 0.5772156649015329d0
         cs = (0.0d0,0.0d0)
         w0 = 0.0d0
         cr = (1.0d0,0.0d0)
         do k = 1 , 50
            w0 = w0 + 1.0d0/k
            cr = 0.25d0*cr/(k*k)*z2
            cp = cr*(w0+ct)
            cs = cs + cp
            if ( k>=10 .and. abs(cp/cs)<1.0d-15 ) exit
         enddo
         cbk0 = ct + cs
      else
         v0n = 1.0d0 - v0
         call gamma2(v0n,gan)
         ca2 = 1.0d0/(gan*(0.5d0*z1)**v0)
         ca1 = (0.5d0*z1)**v0/gap
         csu = ca2 - ca1
         cr1 = (1.0d0,0.0d0)
         cr2 = (1.0d0,0.0d0)
         do k = 1 , 50
            cr1 = 0.25d0*cr1*z2/(k*(k-v0))
            cr2 = 0.25d0*cr2*z2/(k*(k+v0))
            cp = ca2*cr1 - ca1*cr2
            csu = csu + cp
            if ( k>=10 .and. abs(cp/csu)<1.0d-15 ) exit
         enddo
         cbk0 = 0.5d0*pi*csu/sin(piv)
      endif
      Cbk(0) = cbk0
      if ( dble(z)<0.0 ) then
         do k = 0 , n
            cvk = exp((k+v0)*pi*ci)
            if ( dimag(z)<0.0d0 ) then
               Cbk(k) = cvk*Cbk(k) + pi*ci*Cbi(k)
               Cbi(k) = Cbi(k)/cvk
            elseif ( dimag(z)>0.0 ) then
               Cbk(k) = Cbk(k)/cvk - pi*ci*Cbi(k)
               Cbi(k) = cvk*Cbi(k)
            endif
         enddo
      endif
      do k = 1 , n
         ckk = (1.0d0/z-Cbi(k)*Cbk(k-1))/Cbi(k-1)
         Cbk(k) = ckk
      enddo
      Cdi(0) = v0/z*Cbi(0) + Cbi(1)
      Cdk(0) = v0/z*Cbk(0) - Cbk(1)
      do k = 1 , n
         Cdi(k) = -(k+v0)/z*Cbi(k) + Cbi(k-1)
         Cdk(k) = -(k+v0)/z*Cbk(k) - Cbk(k-1)
      enddo
      Vm = n + v0
      end



!       **********************************

      subroutine cikva(v,z,Vm,Cbi,Cdi,Cbk,Cdk)
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
      implicit none
      real(wp) a0 , gan , gap , pi , piv , v , v0 , v0n , v0p , &
                     & Vm , vt , w0 , ws , ws0
      complex(wp) ca , ca1 , ca2 , cb , Cbi , cbi0 , Cbk , cbk0 , cbk1 , &
               & Cdi , Cdk , cf , cf1 , cf2 , cg0 , cg1 , cgk , ci ,    &
               & ci0 , cp
      complex(wp) cr , cr1 , cr2 , cs , csu , ct , cvk , z , z1 , z2
      integer k , k0 , m , n
      dimension Cbi(0:*) , Cdi(0:*) , Cbk(0:*) , Cdk(0:*)
      pi = 3.141592653589793d0
      ci = (0.0d0,1.0d0)
      a0 = abs(z)
      z1 = z
      z2 = z*z
      n = int(v)
      v0 = v - n
      piv = pi*v0
      vt = 4.0d0*v0*v0
      if ( n==0 ) n = 1
      if ( a0<1.0d-100 ) then
         do k = 0 , n
            Cbi(k) = 0.0d0
            Cdi(k) = 0.0d0
            Cbk(k) = -1.0d+300
            Cdk(k) = 1.0d+300
         enddo
         if ( v0==0.0 ) then
            Cbi(0) = (1.0d0,0.0d0)
            Cdi(1) = (0.5d0,0.0d0)
         endif
         Vm = v
         return
      endif
      k0 = 14
      if ( a0>=35.0 ) k0 = 10
      if ( a0>=50.0 ) k0 = 8
      if ( dble(z)<0.0 ) z1 = -z
      if ( a0<18.0 ) then
         if ( v0==0.0 ) then
            ca1 = (1.0d0,0.0d0)
         else
            v0p = 1.0d0 + v0
            call gamma2(v0p,gap)
            ca1 = (0.5d0*z1)**v0/gap
         endif
         ci0 = (1.0d0,0.0d0)
         cr = (1.0d0,0.0d0)
         do k = 1 , 50
            cr = 0.25d0*cr*z2/(k*(k+v0))
            ci0 = ci0 + cr
            if ( abs(cr)<abs(ci0)*1.0d-15 ) exit
         enddo
         cbi0 = ci0*ca1
      else
         ca = exp(z1)/sqrt(2.0d0*pi*z1)
         cs = (1.0d0,0.0d0)
         cr = (1.0d0,0.0d0)
         do k = 1 , k0
            cr = -0.125d0*cr*(vt-(2.0d0*k-1.0d0)**2.0)/(k*z1)
            cs = cs + cr
         enddo
         cbi0 = ca*cs
      endif
      m = msta1(a0,200)
      if ( m<n ) then
         n = m
      else
         m = msta2(a0,n,15)
      endif
      cf2 = (0.0d0,0.0d0)
      cf1 = (1.0d-100,0.0d0)
      do k = m , 0 , -1
         cf = 2.0d0*(v0+k+1.0d0)/z1*cf1 + cf2
         if ( k<=n ) Cbi(k) = cf
         cf2 = cf1
         cf1 = cf
      enddo
      cs = cbi0/cf
      do k = 0 , n
         Cbi(k) = cs*Cbi(k)
      enddo
      if ( a0>9.0 ) then
         cb = exp(-z1)*sqrt(0.5d0*pi/z1)
         cs = (1.0d0,0.0d0)
         cr = (1.0d0,0.0d0)
         do k = 1 , k0
            cr = 0.125d0*cr*(vt-(2.0d0*k-1.0d0)**2.0)/(k*z1)
            cs = cs + cr
         enddo
         cbk0 = cb*cs
      elseif ( v0==0.0 ) then
         ct = -log(0.5d0*z1) - 0.5772156649015329d0
         cs = (0.0d0,0.0d0)
         w0 = 0.0d0
         cr = (1.0d0,0.0d0)
         do k = 1 , 50
            w0 = w0 + 1.0d0/k
            cr = 0.25d0*cr/(k*k)*z2
            cp = cr*(w0+ct)
            cs = cs + cp
            if ( k>=10 .and. abs(cp/cs)<1.0d-15 ) exit
         enddo
         cbk0 = ct + cs
      else
         v0n = 1.0d0 - v0
         call gamma2(v0n,gan)
         ca2 = 1.0d0/(gan*(0.5d0*z1)**v0)
         ca1 = (0.5d0*z1)**v0/gap
         csu = ca2 - ca1
         cr1 = (1.0d0,0.0d0)
         cr2 = (1.0d0,0.0d0)
         ws0 = 0.0d0
         do k = 1 , 50
            cr1 = 0.25d0*cr1*z2/(k*(k-v0))
            cr2 = 0.25d0*cr2*z2/(k*(k+v0))
            csu = csu + ca2*cr1 - ca1*cr2
            ws = abs(csu)
            if ( k>=10 .and. abs(ws-ws0)/ws<1.0d-15 ) exit
            ws0 = ws
         enddo
         cbk0 = 0.5d0*pi*csu/sin(piv)
      endif
      cbk1 = (1.0d0/z1-Cbi(1)*cbk0)/Cbi(0)
      Cbk(0) = cbk0
      Cbk(1) = cbk1
      cg0 = cbk0
      cg1 = cbk1
      do k = 2 , n
         cgk = 2.0d0*(v0+k-1.0d0)/z1*cg1 + cg0
         Cbk(k) = cgk
         cg0 = cg1
         cg1 = cgk
      enddo
      if ( dble(z)<0.0 ) then
         do k = 0 , n
            cvk = exp((k+v0)*pi*ci)
            if ( dimag(z)<0.0d0 ) then
               Cbk(k) = cvk*Cbk(k) + pi*ci*Cbi(k)
               Cbi(k) = Cbi(k)/cvk
            elseif ( dimag(z)>0.0 ) then
               Cbk(k) = Cbk(k)/cvk - pi*ci*Cbi(k)
               Cbi(k) = cvk*Cbi(k)
            endif
         enddo
      endif
      Cdi(0) = v0/z*Cbi(0) + Cbi(1)
      Cdk(0) = v0/z*Cbk(0) - Cbk(1)
      do k = 1 , n
         Cdi(k) = -(k+v0)/z*Cbi(k) + Cbi(k-1)
         Cdk(k) = -(k+v0)/z*Cbk(k) - Cbk(k-1)
      enddo
      Vm = n + v0
      end



!       **********************************

      subroutine cfc(z,Zf,Zd)
!
!       =========================================================
!       Purpose: Compute complex Fresnel integral C(z) and C'(z)
!       Input :  z --- Argument of C(z)
!       Output:  ZF --- C(z)
!                ZD --- C'(z)
!       =========================================================
!
      implicit none
      complex(wp) c , cf , cf0 , cf1 , cg , cr , d , z , z0 , Zd , Zf ,  &
               & zp , zp2
      real(wp) eps , pi , w0 , wa , wa0
      integer k , m
      eps = 1.0d-14
      pi = 3.141592653589793d0
      w0 = abs(z)
      zp = 0.5d0*pi*z*z
      zp2 = zp*zp
      z0 = (0.0d0,0.0d0)
      if ( z==z0 ) then
         c = z0
      elseif ( w0<=2.5 ) then
         cr = z
         c = cr
         wa0 = 0.0d0
         do k = 1 , 80
            cr = -.5d0*cr*(4.0d0*k-3.0d0)/k/(2.0d0*k-1.0d0)             &
               & /(4.0d0*k+1.0d0)*zp2
            c = c + cr
            wa = abs(c)
            if ( abs((wa-wa0)/wa)<eps .and. k>10 ) goto 100
            wa0 = wa
         enddo
      elseif ( w0>2.5 .and. w0<4.5 ) then
         m = 85
         c = z0
         cf1 = z0
         cf0 = (1.0d-100,0.0d0)
         do k = m , 0 , -1
            cf = (2.0d0*k+3.0d0)*cf0/zp - cf1
            if ( k==int(k/2)*2 ) c = c + cf
            cf1 = cf0
            cf0 = cf
         enddo
         c = 2.0d0/(pi*z)*sin(zp)/cf*c
      else
!          See comment at CFS(), use C(z) = iC(-iz)
         if ( dimag(z)>-dble(z) .and. dimag(z)<=dble(z) ) then
!            right quadrant
            d = dcmplx(.5d0,0.0d0)
         elseif ( dimag(z)>dble(z) .and. dimag(z)>=-dble(z) ) then
!            upper quadrant
            d = dcmplx(0.0d0,.5d0)
         elseif ( dimag(z)<-dble(z) .and. dimag(z)>=dble(z) ) then
!            left quadrant
            d = dcmplx(-.5d0,0.0d0)
         else
!            lower quadrant
            d = dcmplx(0.0d0,-.5d0)
         endif
         cr = (1.0d0,0.0d0)
         cf = (1.0d0,0.0d0)
         do k = 1 , 20
            cr = -.25d0*cr*(4.0d0*k-1.0d0)*(4.0d0*k-3.0d0)/zp2
            cf = cf + cr
         enddo
         cr = 1.0d0/(pi*z*z)
         cg = cr
         do k = 1 , 12
            cr = -.25d0*cr*(4.0d0*k+1.0d0)*(4.0d0*k-1.0d0)/zp2
            cg = cg + cr
         enddo
         c = d + (cf*sin(zp)-cg*cos(zp))/(pi*z)
      endif
 100  Zf = c
      Zd = cos(0.5*pi*z*z)
      end



!       **********************************

      subroutine fcs(x,c,s)
!
!       =================================================
!       Purpose: Compute Fresnel integrals C(x) and S(x)
!       Input :  x --- Argument of C(x) and S(x)
!       Output:  C --- C(x)
!                S --- S(x)
!       =================================================
!
      implicit none
      real(wp) c , eps , f , f0 , f1 , g , pi , px , q , r , s ,&
                     & su , t , t0 , t2 , x , xa
      integer k , m
      eps = 1.0d-15
      pi = 3.141592653589793d0
      xa = abs(x)
      px = pi*xa
      t = .5d0*px*xa
      t2 = t*t
      if ( xa==0.0 ) then
         c = 0.0d0
         s = 0.0d0
      elseif ( xa<2.5d0 ) then
         r = xa
         c = r
         do k = 1 , 50
            r = -.5d0*r*(4.0d0*k-3.0d0)/k/(2.0d0*k-1.0d0)               &
              & /(4.0d0*k+1.0d0)*t2
            c = c + r
            if ( abs(r)<abs(c)*eps ) exit
         enddo
         s = xa*t/3.0d0
         r = s
         do k = 1 , 50
            r = -.5d0*r*(4.0d0*k-1.0d0)/k/(2.0d0*k+1.0d0)               &
              & /(4.0d0*k+3.0d0)*t2
            s = s + r
            if ( abs(r)<abs(s)*eps ) goto 100
         enddo
      elseif ( xa<4.5d0 ) then
         m = int(42.0+1.75*t)
         su = 0.0d0
         c = 0.0d0
         s = 0.0d0
         f1 = 0.0d0
         f0 = 1.0d-100
         do k = m , 0 , -1
            f = (2.0d0*k+3.0d0)*f0/t - f1
            if ( k==int(k/2)*2 ) then
               c = c + f
            else
               s = s + f
            endif
            su = su + (2.0d0*k+1.0d0)*f*f
            f1 = f0
            f0 = f
         enddo
         q = sqrt(su)
         c = c*xa/q
         s = s*xa/q
      else
         r = 1.0d0
         f = 1.0d0
         do k = 1 , 20
            r = -.25d0*r*(4.0d0*k-1.0d0)*(4.0d0*k-3.0d0)/t2
            f = f + r
         enddo
         r = 1.0d0/(px*xa)
         g = r
         do k = 1 , 12
            r = -.25d0*r*(4.0d0*k+1.0d0)*(4.0d0*k-1.0d0)/t2
            g = g + r
         enddo
         t0 = t - int(t/(2.0d0*pi))*2.0d0*pi
         c = .5d0 + (f*sin(t0)-g*cos(t0))/px
         s = .5d0 - (f*cos(t0)+g*sin(t0))/px
      endif
 100  if ( x<0.0d0 ) then
         c = -c
         s = -s
      endif
      end

!       **********************************

      subroutine rctj(n,x,Nm,Rj,Dj)
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
      implicit none
      real(wp) cs , Dj , f , f0 , f1 , Rj , rj0 , rj1 , x
      integer k , m , n , Nm
      dimension Rj(0:n) , Dj(0:n)
      Nm = n
      if ( abs(x)<1.0d-100 ) then
         do k = 0 , n
            Rj(k) = 0.0d0
            Dj(k) = 0.0d0
         enddo
         Dj(0) = 1.0d0
         return
      endif
      Rj(0) = sin(x)
      Rj(1) = Rj(0)/x - cos(x)
      rj0 = Rj(0)
      rj1 = Rj(1)
      cs = 0.0d0
      f = 0.0d0
      if ( n>=2 ) then
         m = msta1(x,200)
         if ( m<n ) then
            Nm = m
         else
            m = msta2(x,n,15)
         endif
         f0 = 0.0d0
         f1 = 1.0d-100
         do k = m , 0 , -1
            f = (2.0d0*k+3.0d0)*f1/x - f0
            if ( k<=Nm ) Rj(k) = f
            f0 = f1
            f1 = f
         enddo
         if ( abs(rj0)>abs(rj1) ) cs = rj0/f
         if ( abs(rj0)<=abs(rj1) ) cs = rj1/f0
         do k = 0 , Nm
            Rj(k) = cs*Rj(k)
         enddo
      endif
      Dj(0) = cos(x)
      do k = 1 , Nm
         Dj(k) = -k*Rj(k)/x + Rj(k-1)
      enddo
      end



!       **********************************

      subroutine herzo(n,x,w)
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
      implicit none
      real(wp) f0 , f1 , fd , gd , hd , hf , hn , p , q , r ,   &
                     & r1 , r2 , w , wp , x , z , z0 , zl
      integer i , it , j , k , n , nr
      dimension x(n) , w(n)
      hn = 1.0d0/n
      zl = -1.1611d0 + 1.46d0*n**0.5
      z = 0.0d0
      hf = 0.0d0
      hd = 0.0d0
      do nr = 1 , n/2
         if ( nr==1 ) z = zl
         if ( nr/=1 ) z = z - hn*(n/2+1-nr)
         it = 0
 50      it = it + 1
         z0 = z
         f0 = 1.0d0
         f1 = 2.0d0*z
         do k = 2 , n
            hf = 2.0d0*z*f1 - 2.0d0*(k-1.0d0)*f0
            hd = 2.0d0*k*f1
            f0 = f1
            f1 = hf
         enddo
         p = 1.0d0
         do i = 1 , nr - 1
            p = p*(z-x(i))
         enddo
         fd = hf/p
         q = 0.0d0
         do i = 1 , nr - 1
            wp = 1.0d0
            do j = 1 , nr - 1
               if ( j/=i ) wp = wp*(z-x(j))
            enddo
            q = q + wp
         enddo
         gd = (hd-q*fd)/p
         z = z - fd/gd
         if ( it<=40 .and. abs((z-z0)/z)>1.0d-15 ) goto 50
         x(nr) = z
         x(n+1-nr) = -z
         r = 1.0d0
         do k = 1 , n
            r = 2.0d0*r*k
         enddo
         w(nr) = 3.544907701811d0*r/(hd*hd)
         w(n+1-nr) = w(nr)
      enddo
      if ( n/=2*int(n/2) ) then
         r1 = 1.0d0
         r2 = 1.0d0
         do j = 1 , n
            r1 = 2.0d0*r1*j
            if ( j>=(n+1)/2 ) r2 = r2*j
         enddo
         w(n/2+1) = 0.88622692545276d0*r1/(r2*r2)
         x(n/2+1) = 0.0d0
      endif
      end

!       **********************************

      subroutine jy01b(x,Bj0,Dj0,Bj1,Dj1,By0,Dy0,By1,Dy1)
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
      implicit none
      real(wp) a0 , Bj0 , Bj1 , By0 , By1 , Dj0 , Dj1 , Dy0 ,   &
                     & Dy1 , p0 , p1 , pi , q0 , q1 , t , t2 , ta0 ,    &
                     & ta1 , x
      pi = 3.141592653589793d0
      if ( x==0.0d0 ) then
         Bj0 = 1.0d0
         Bj1 = 0.0d0
         Dj0 = 0.0d0
         Dj1 = 0.5d0
         By0 = -1.0d+300
         By1 = -1.0d+300
         Dy0 = 1.0d+300
         Dy1 = 1.0d+300
         return
      elseif ( x<=4.0d0 ) then
         t = x/4.0d0
         t2 = t*t
         Bj0 = ((((((-.5014415d-3*t2+.76771853d-2)*t2-.0709253492d0)*t2+&
             & .4443584263d0)*t2-1.7777560599d0)*t2+3.9999973021d0)     &
             & *t2-3.9999998721d0)*t2 + 1.0d0
         Bj1 = t*                                                       &
             & (((((((-.1289769d-3*t2+.22069155d-2)*t2-.0236616773d0)*t2&
             & +.1777582922d0)*t2-.8888839649d0)*t2+2.6666660544d0)     &
             & *t2-3.9999999710d0)*t2+1.9999999998d0)
         By0 = (((((((-.567433d-4*t2+.859977d-3)*t2-.94855882d-2)*t2+   &
             & .0772975809d0)*t2-.4261737419d0)*t2+1.4216421221d0)      &
             & *t2-2.3498519931d0)*t2+1.0766115157d0)*t2 + .3674669052d0
         By0 = 2.0d0/pi*log(x/2.0d0)*Bj0 + By0
         By1 = ((((((((.6535773d-3*t2-.0108175626d0)*t2+.107657606d0)*t2&
             & -.7268945577d0)*t2+3.1261399273d0)*t2-7.3980241381d0)    &
             & *t2+6.8529236342d0)*t2+.3932562018d0)*t2-.6366197726d0)/x
         By1 = 2.0d0/pi*log(x/2.0d0)*Bj1 + By1
      else
         t = 4.0d0/x
         t2 = t*t
         a0 = sqrt(2.0d0/(pi*x))
         p0 = ((((-.9285d-5*t2+.43506d-4)*t2-.122226d-3)*t2+.434725d-3) &
            & *t2-.4394275d-2)*t2 + .999999997d0
         q0 = t*(((((.8099d-5*t2-.35614d-4)*t2+.85844d-4)*t2-.218024d-3)&
            & *t2+.1144106d-2)*t2-.031249995d0)
         ta0 = x - .25d0*pi
         Bj0 = a0*(p0*cos(ta0)-q0*sin(ta0))
         By0 = a0*(p0*sin(ta0)+q0*cos(ta0))
         p1 = ((((.10632d-4*t2-.50363d-4)*t2+.145575d-3)*t2-.559487d-3) &
            & *t2+.7323931d-2)*t2 + 1.000000004d0
         q1 = t*                                                        &
            & (((((-.9173d-5*t2+.40658d-4)*t2-.99941d-4)*t2+.266891d-3) &
            & *t2-.1601836d-2)*t2+.093749994d0)
         ta1 = x - .75d0*pi
         Bj1 = a0*(p1*cos(ta1)-q1*sin(ta1))
         By1 = a0*(p1*sin(ta1)+q1*cos(ta1))
      endif
      Dj0 = -Bj1
      Dj1 = Bj0 - Bj1/x
      Dy0 = -By1
      Dy1 = By0 - By1/x
      end

!       **********************************

      subroutine enxb(n,x,En)
!
!       ===============================================
!       Purpose: Compute exponential integral En(x)
!       Input :  x --- Argument of En(x)
!                n --- Order of En(x)  (n = 0,1,2,...)
!       Output:  EN(n) --- En(x)
!       ===============================================
!
      implicit none
      real(wp) En , ens , ps , r , rp , s , s0 , t , t0 , x
      integer j , k , l , m , n
      dimension En(0:n)
      if ( x==0.0 ) then
         En(0) = 1.0d+300
         En(1) = 1.0d+300
         do k = 2 , n
            En(k) = 1.0d0/(k-1.0)
         enddo
         return
      elseif ( x<=1.0 ) then
         En(0) = exp(-x)/x
         s0 = 0.0d0
         do l = 1 , n
            rp = 1.0d0
            do j = 1 , l - 1
               rp = -rp*x/j
            enddo
            ps = -0.5772156649015328d0
            do m = 1 , l - 1
               ps = ps + 1.0d0/m
            enddo
            ens = rp*(-log(x)+ps)
            s = 0.0d0
            do m = 0 , 20
               if ( m/=l-1 ) then
                  r = 1.0d0
                  do j = 1 , m
                     r = -r*x/j
                  enddo
                  s = s + r/(m-l+1.0d0)
                  if ( abs(s-s0)<abs(s)*1.0d-15 ) exit
                  s0 = s
               endif
            enddo
            En(l) = ens - s
         enddo
      else
         En(0) = exp(-x)/x
         m = 15 + int(100.0/x)
         do l = 1 , n
            t0 = 0.0d0
            do k = m , 1 , -1
               t0 = (l+k-1.0d0)/(1.0d0+k/(x+t0))
            enddo
            t = 1.0d0/(x+t0)
            En(l) = exp(-x)*t
         enddo
      endif
      end

!       **********************************

      subroutine sphk(n,x,Nm,Sk,Dk)
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
      implicit none
      real(wp) Dk , f , f0 , f1 , pi , Sk , x
      integer k , n , Nm
      dimension Sk(0:n) , Dk(0:n)
      pi = 3.141592653589793d0
      Nm = n
      if ( x<1.0d-60 ) then
         do k = 0 , n
            Sk(k) = 1.0d+300
            Dk(k) = -1.0d+300
         enddo
         return
      endif
      Sk(0) = 0.5d0*pi/x*exp(-x)
      Sk(1) = Sk(0)*(1.0d0+1.0d0/x)
      f0 = Sk(0)
      f1 = Sk(1)
      do k = 2 , n
         f = (2.0d0*k-1.0d0)*f1/x + f0
         Sk(k) = f
         if ( abs(f)>1.0d+300 ) exit
         f0 = f1
         f1 = f
      enddo
      Nm = k - 1
      Dk(0) = -Sk(1)
      do k = 1 , Nm
         Dk(k) = -Sk(k-1) - (k+1.0d0)/x*Sk(k)
      enddo
      end

!       **********************************

      subroutine enxa(n,x,En)
!
!       ============================================
!       Purpose: Compute exponential integral En(x)
!       Input :  x --- Argument of En(x) ( x ≤ 20 )
!                n --- Order of En(x)
!       Output:  EN(n) --- En(x)
!       Routine called: E1XB for computing E1(x)
!       ============================================
!
      implicit none
      real(wp) e1 , ek , En , x
      integer k , n
      dimension En(0:n)
      En(0) = exp(-x)/x
      call e1xb(x,e1)
      En(1) = e1
      do k = 2 , n
         ek = (exp(-x)-x*e1)/(k-1.0d0)
         En(k) = ek
         e1 = ek
      enddo
      end



!       **********************************

      subroutine gaih(x,Ga)
!
!       =====================================================
!       Purpose: Compute gamma function Г(x)
!       Input :  x  --- Argument of Г(x), x = n/2, n=1,2,…
!       Output:  GA --- Г(x)
!       =====================================================
!
      implicit none
      real(wp) Ga , pi , x
      integer k , m , m1
      pi = 3.141592653589793d0
      if ( x==int(x) .and. x>0.0 ) then
         Ga = 1.0d0
         m1 = int(x-1.0)
         do k = 2 , m1
            Ga = Ga*k
         enddo
      elseif ( x+.5d0==int(x+.5d0) .and. x>0.0 ) then
         m = int(x)
         Ga = sqrt(pi)
         do k = 1 , m
            Ga = 0.5d0*Ga*(2.0d0*k-1.0d0)
         enddo
      endif
      end

!       **********************************

      subroutine pbvv(v,x,Vv,Vp,Pvf,Pvd)
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
      implicit none
      real(wp) f , f0 , f1 , pi , pv0 , Pvd , Pvf , q2p , qe ,  &
                     & s0 , v , v0 , v1 , v2 , vh , Vp , Vv , x , xa
      integer ja , k , kv , l , m , na , nv
      dimension Vv(0:*) , Vp(0:*)
      pi = 3.141592653589793d0
      xa = abs(x)
      vh = v
      v = v + sign(1.0d0,v)
      nv = int(v)
      v0 = v - nv
      na = abs(nv)
      qe = exp(0.25d0*x*x)
      q2p = sqrt(2.0d0/pi)
      ja = 0
      if ( na>=1 ) ja = 1
      f = 0.0d0
      if ( v<=0.0 ) then
         if ( v0==0.0 ) then
            if ( xa<=7.5 ) call vvsa(v0,x,pv0)
            if ( xa>7.5 ) call vvla(v0,x,pv0)
            f0 = q2p*qe
            f1 = x*f0
            Vv(0) = pv0
            Vv(1) = f0
            Vv(2) = f1
         else
            do l = 0 , ja
               v1 = v0 - l
               if ( xa<=7.5 ) call vvsa(v1,x,f1)
               if ( xa>7.5 ) call vvla(v1,x,f1)
               if ( l==0 ) f0 = f1
            enddo
            Vv(0) = f0
            Vv(1) = f1
         endif
         kv = 2
         if ( v0==0.0 ) kv = 3
         do k = kv , na
            f = x*f1 + (k-v0-2.0d0)*f0
            Vv(k) = f
            f0 = f1
            f1 = f
         enddo
      elseif ( x>=0.0 .and. x<=7.5d0 ) then
         v2 = v
         if ( v2<1.0 ) v2 = v2 + 1.0d0
         call vvsa(v2,x,f1)
         v1 = v2 - 1.0d0
         kv = int(v2)
         call vvsa(v1,x,f0)
         Vv(kv) = f1
         Vv(kv-1) = f0
         do k = kv - 2 , 0 , -1
            f = x*f0 - (k+v0+2.0d0)*f1
            if ( k<=na ) Vv(k) = f
            f1 = f0
            f0 = f
         enddo
      elseif ( x>7.5d0 ) then
         call vvla(v0,x,pv0)
         m = 100 + abs(na)
         Vv(1) = pv0
         f1 = 0.0d0
         f0 = 1.0d-40
         do k = m , 0 , -1
            f = x*f0 - (k+v0+2.0d0)*f1
            if ( k<=na ) Vv(k) = f
            f1 = f0
            f0 = f
         enddo
         s0 = pv0/f
         do k = 0 , na
            Vv(k) = s0*Vv(k)
         enddo
      else
         if ( xa<=7.5d0 ) then
            call vvsa(v0,x,f0)
            v1 = v0 + 1.0
            call vvsa(v1,x,f1)
         else
            call vvla(v0,x,f0)
            v1 = v0 + 1.0d0
            call vvla(v1,x,f1)
         endif
         Vv(0) = f0
         Vv(1) = f1
         do k = 2 , na
            f = (x*f1-f0)/(k+v0)
            Vv(k) = f
            f0 = f1
            f1 = f
         enddo
      endif
      do k = 0 , na - 1
         v1 = v0 + k
         if ( v>=0.0d0 ) then
            Vp(k) = 0.5d0*x*Vv(k) - (v1+1.0d0)*Vv(k+1)
         else
            Vp(k) = -0.5d0*x*Vv(k) + Vv(k+1)
         endif
      enddo
      Pvf = Vv(na-1)
      Pvd = Vp(na-1)
      v = vh
      end



!       **********************************

      subroutine clqmn(Mm,m,n,x,y,Cqm,Cqd)
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
      implicit none
      complex(wp) cq0 , cq1 , cq10 , Cqd , cqf , cqf0 , cqf1 , cqf2 ,    &
               & Cqm , z , zq , zs
      integer i , j , k , km , ls , m , Mm , n
      real(wp) x , xc , y
      dimension Cqm(0:Mm,0:n) , Cqd(0:Mm,0:n)
      z = dcmplx(x,y)
      if ( abs(x)==1.0d0 .and. y==0.0d0 ) then
         do i = 0 , m
            do j = 0 , n
               Cqm(i,j) = (1.0d+300,0.0d0)
               Cqd(i,j) = (1.0d+300,0.0d0)
            enddo
         enddo
         return
      endif
      xc = abs(z)
      ls = 0
      if ( dimag(z)==0.0d0 .or. xc<1.0d0 ) ls = 1
      if ( xc>1.0d0 ) ls = -1
      zq = sqrt(ls*(1.0d0-z*z))
      zs = ls*(1.0d0-z*z)
      cq0 = 0.5d0*log(ls*(1.0d0+z)/(1.0d0-z))
      if ( xc<1.0001d0 ) then
         Cqm(0,0) = cq0
         Cqm(0,1) = z*cq0 - 1.0d0
         Cqm(1,0) = -1.0d0/zq
         Cqm(1,1) = -zq*(cq0+z/(1.0d0-z*z))
         do i = 0 , 1
            do j = 2 , n
               Cqm(i,j) = ((2.0d0*j-1.0d0)*z*Cqm(i,j-1)-(j+i-1.0d0)     &
                        & *Cqm(i,j-2))/(j-i)
            enddo
         enddo
         do j = 0 , n
            do i = 2 , m
               Cqm(i,j) = -2.0d0*(i-1.0d0)*z/zq*Cqm(i-1,j)              &
                        & - ls*(j+i-1.0d0)*(j-i+2.0d0)*Cqm(i-2,j)
            enddo
         enddo
      else
         if ( xc>1.1 ) then
            km = 40 + m + n
         else
            km = (40+m+n)*int(-1.0-1.8*log(xc-1.0))
         endif
         cqf2 = (0.0d0,0.0d0)
         cqf1 = (1.0d0,0.0d0)
         do k = km , 0 , -1
            cqf0 = ((2*k+3.0d0)*z*cqf1-(k+2.0d0)*cqf2)/(k+1.0d0)
            if ( k<=n ) Cqm(0,k) = cqf0
            cqf2 = cqf1
            cqf1 = cqf0
         enddo
         do k = 0 , n
            Cqm(0,k) = cq0*Cqm(0,k)/cqf0
         enddo
         cqf2 = 0.0d0
         cqf1 = 1.0d0
         do k = km , 0 , -1
            cqf0 = ((2*k+3.0d0)*z*cqf1-(k+1.0d0)*cqf2)/(k+2.0d0)
            if ( k<=n ) Cqm(1,k) = cqf0
            cqf2 = cqf1
            cqf1 = cqf0
         enddo
         cq10 = -1.0d0/zq
         do k = 0 , n
            Cqm(1,k) = cq10*Cqm(1,k)/cqf0
         enddo
         do j = 0 , n
            cq0 = Cqm(0,j)
            cq1 = Cqm(1,j)
            do i = 0 , m - 2
               cqf = -2.0d0*(i+1)*z/zq*cq1 + (j-i)*(j+i+1.0d0)*cq0
               Cqm(i+2,j) = cqf
               cq0 = cq1
               cq1 = cqf
            enddo
         enddo
      endif
      Cqd(0,0) = ls/zs
      do j = 1 , n
         Cqd(0,j) = ls*j*(Cqm(0,j-1)-z*Cqm(0,j))/zs
      enddo
      do j = 0 , n
         do i = 1 , m
            Cqd(i,j) = ls*i*z/zs*Cqm(i,j) + (i+j)*(j-i+1.0d0)           &
                     & /zq*Cqm(i-1,j)
         enddo
      enddo
      end


!       **********************************

      subroutine segv(m,n,c,Kd,Cv,Eg)
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
      implicit none
      real(wp) a , b , c , cs , Cv , cv0 , d , d2k , dk0 , dk1 ,&
                     & dk2 , e , Eg , f , g , h , s , t , t1 , x1
      real(wp) xa , xb
      integer i , icm , j , k , k1 , Kd , l , m , n , nm , nm1
      dimension b(100) , h(100) , d(300) , e(300) , f(300) , cv0(100) , &
              & a(300) , g(300) , Eg(200)
      if ( c<1.0d-10 ) then
         do i = 1 , n - m + 1
            Eg(i) = (i+m)*(i+m-1.0d0)
         enddo
         goto 100
      endif
      icm = (n-m+2)/2
      nm = 10 + int(0.5*(n-m)+c)
      cs = c*c*Kd
      k = 0
      do l = 0 , 1
         do i = 1 , nm
            if ( l==0 ) k = 2*(i-1)
            if ( l==1 ) k = 2*i - 1
            dk0 = m + k
            dk1 = m + k + 1
            dk2 = 2*(m+k)
            d2k = 2*m + k
            a(i) = (d2k+2.0)*(d2k+1.0)/((dk2+3.0)*(dk2+5.0))*cs
            d(i) = dk0*dk1 + (2.0*dk0*dk1-2.0*m*m-1.0)                  &
                 & /((dk2-1.0)*(dk2+3.0))*cs
            g(i) = k*(k-1.0)/((dk2-3.0)*(dk2-1.0))*cs
         enddo
         do k = 2 , nm
            e(k) = sqrt(a(k-1)*g(k))
            f(k) = e(k)*e(k)
         enddo
         f(1) = 0.0d0
         e(1) = 0.0d0
         xa = d(nm) + abs(e(nm))
         xb = d(nm) - abs(e(nm))
         nm1 = nm - 1
         do i = 1 , nm1
            t = abs(e(i)) + abs(e(i+1))
            t1 = d(i) + t
            if ( xa<t1 ) xa = t1
            t1 = d(i) - t
            if ( t1<xb ) xb = t1
         enddo
         do i = 1 , icm
            b(i) = xa
            h(i) = xb
         enddo
         do k = 1 , icm
            do k1 = k , icm
               if ( b(k1)<b(k) ) then
                  b(k) = b(k1)
                  exit
               endif
            enddo
            if ( k/=1 ) then
               if ( h(k)<h(k-1) ) h(k) = h(k-1)
            endif
 40         x1 = (b(k)+h(k))/2.0d0
            cv0(k) = x1
            if ( abs((b(k)-h(k))/x1)<1.0d-14 ) then
               cv0(k) = x1
               if ( l==0 ) Eg(2*k-1) = cv0(k)
               if ( l==1 ) Eg(2*k) = cv0(k)
            else
               j = 0
               s = 1.0d0
               do i = 1 , nm
                  if ( s==0.0d0 ) s = s + 1.0d-30
                  t = f(i)/s
                  s = d(i) - t - x1
                  if ( s<0.0d0 ) j = j + 1
               enddo
               if ( j<k ) then
                  h(k) = x1
               else
                  b(k) = x1
                  if ( j>=icm ) then
                     b(icm) = x1
                  else
                     if ( h(j+1)<x1 ) h(j+1) = x1
                     if ( x1<b(j) ) b(j) = x1
                  endif
               endif
               goto 40
            endif
         enddo
      enddo
 100  Cv = Eg(n-m+1)
      end


!       **********************************

      subroutine ciknb(n,z,Nm,Cbi,Cdi,Cbk,Cdk)
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
      implicit none
      real(wp) a0 , el , fac , pi , vt
      complex(wp) ca0 , Cbi , Cbk , cbkl , cbs , Cdi , Cdk , cf , cf0 ,  &
               & cf1 , cg , cg0 , cg1 , ci , cr , cs0 , csk0 , z , z1
      integer k , k0 , l , m , n , Nm
      dimension Cbi(0:n) , Cdi(0:n) , Cbk(0:n) , Cdk(0:n)
      pi = 3.141592653589793d0
      el = 0.57721566490153d0
      a0 = abs(z)
      Nm = n
      if ( a0<1.0d-100 ) then
         do k = 0 , n
            Cbi(k) = (0.0d0,0.0d0)
            Cbk(k) = (1.0d+300,0.0d0)
            Cdi(k) = (0.0d0,0.0d0)
            Cdk(k) = -(1.0d+300,0.0d0)
         enddo
         Cbi(0) = (1.0d0,0.0d0)
         Cdi(1) = (0.5d0,0.0d0)
         return
      endif
      z1 = z
      ci = (0.0d0,1.0d0)
      if ( dble(z)<0.0 ) z1 = -z
      if ( n==0 ) Nm = 1
      m = msta1(a0,200)
      if ( m<Nm ) then
         Nm = m
      else
         m = msta2(a0,Nm,15)
      endif
      cbs = 0.0d0
      csk0 = 0.0d0
      cf0 = 0.0d0
      cf1 = 1.0d-100
      do k = m , 0 , -1
         cf = 2.0d0*(k+1.0d0)*cf1/z1 + cf0
         if ( k<=Nm ) Cbi(k) = cf
         if ( k/=0 .and. k==2*int(k/2) ) csk0 = csk0 + 4.0d0*cf/k
         cbs = cbs + 2.0d0*cf
         cf0 = cf1
         cf1 = cf
      enddo
      cs0 = exp(z1)/(cbs-cf)
      do k = 0 , Nm
         Cbi(k) = cs0*Cbi(k)
      enddo
      if ( a0<=9.0 ) then
         Cbk(0) = -(log(0.5d0*z1)+el)*Cbi(0) + cs0*csk0
         Cbk(1) = (1.0d0/z1-Cbi(1)*Cbk(0))/Cbi(0)
      else
         ca0 = sqrt(pi/(2.0d0*z1))*exp(-z1)
         k0 = 16
         if ( a0>=25.0 ) k0 = 10
         if ( a0>=80.0 ) k0 = 8
         if ( a0>=200.0 ) k0 = 6
         do l = 0 , 1
            cbkl = 1.0d0
            vt = 4.0d0*l
            cr = (1.0d0,0.0d0)
            do k = 1 , k0
               cr = 0.125d0*cr*(vt-(2.0*k-1.0)**2)/(k*z1)
               cbkl = cbkl + cr
            enddo
            Cbk(l) = ca0*cbkl
         enddo
      endif
      cg0 = Cbk(0)
      cg1 = Cbk(1)
      do k = 2 , Nm
         cg = 2.0d0*(k-1.0d0)/z1*cg1 + cg0
         Cbk(k) = cg
         cg0 = cg1
         cg1 = cg
      enddo
      if ( dble(z)<0.0 ) then
         fac = 1.0d0
         do k = 0 , Nm
            if ( dimag(z)<0.0 ) then
               Cbk(k) = fac*Cbk(k) + ci*pi*Cbi(k)
            else
               Cbk(k) = fac*Cbk(k) - ci*pi*Cbi(k)
            endif
            Cbi(k) = fac*Cbi(k)
            fac = -fac
         enddo
      endif
      Cdi(0) = Cbi(1)
      Cdk(0) = -Cbk(1)
      do k = 1 , Nm
         Cdi(k) = Cbi(k-1) - k/z*Cbi(k)
         Cdk(k) = -Cbk(k-1) - k/z*Cbk(k)
      enddo
      end


!       **********************************

      subroutine cikna(n,z,Nm,Cbi,Cdi,Cbk,Cdk)
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
      implicit none
      real(wp) a0
      complex(wp) Cbi , cbi0 , cbi1 , Cbk , cbk0 , cbk1 , Cdi , cdi0 ,   &
               & cdi1 , Cdk , cdk0 , cdk1 , cf , cf1 , cf2 , ckk , cs , &
               & z
      integer k , m , n , Nm
      dimension Cbi(0:n) , Cdi(0:n) , Cbk(0:n) , Cdk(0:n)
      a0 = abs(z)
      Nm = n
      if ( a0<1.0d-100 ) then
         do k = 0 , n
            Cbi(k) = (0.0d0,0.0d0)
            Cdi(k) = (0.0d0,0.0d0)
            Cbk(k) = -(1.0d+300,0.0d0)
            Cdk(k) = (1.0d+300,0.0d0)
         enddo
         Cbi(0) = (1.0d0,0.0d0)
         Cdi(1) = (0.5d0,0.0d0)
         return
      endif
      call cik01(z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1)
      Cbi(0) = cbi0
      Cbi(1) = cbi1
      Cbk(0) = cbk0
      Cbk(1) = cbk1
      Cdi(0) = cdi0
      Cdi(1) = cdi1
      Cdk(0) = cdk0
      Cdk(1) = cdk1
      if ( n<=1 ) return
      m = msta1(a0,200)
      if ( m<n ) then
         Nm = m
      else
         m = msta2(a0,n,15)
      endif
      cf2 = (0.0d0,0.0d0)
      cf1 = (1.0d-100,0.0d0)
      do k = m , 0 , -1
         cf = 2.0d0*(k+1.0d0)/z*cf1 + cf2
         if ( k<=Nm ) Cbi(k) = cf
         cf2 = cf1
         cf1 = cf
      enddo
      cs = cbi0/cf
      do k = 0 , Nm
         Cbi(k) = cs*Cbi(k)
      enddo
      do k = 2 , Nm
         if ( abs(Cbi(k-1))>abs(Cbi(k-2)) ) then
            ckk = (1.0d0/z-Cbi(k)*Cbk(k-1))/Cbi(k-1)
         else
            ckk = (Cbi(k)*Cbk(k-2)+2.0d0*(k-1.0d0)/(z*z))/Cbi(k-2)
         endif
         Cbk(k) = ckk
      enddo
      do k = 2 , Nm
         Cdi(k) = Cbi(k-1) - k/z*Cbi(k)
         Cdk(k) = -Cbk(k-1) - k/z*Cbk(k)
      enddo
      end



!       **********************************

      subroutine mtu12(Kf,Kc,m,q,x,F1r,D1r,F2r,D2r)
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
      implicit none
      real(wp) a , bj1 , bj2 , by1 , by2 , c1 , c2 , D1r , D2r ,&
                     & dj1 , dj2 , dy1 , dy2 , eps , F1r , F2r , &
                     & fg , q , qm
      real(wp) u1 , u2 , w1 , w2 , x
      integer ic , k , Kc , kd , Kf , km , m , nm
      dimension fg(251) , bj1(0:251) , dj1(0:251) , bj2(0:251) ,        &
              & dj2(0:251) , by1(0:251) , dy1(0:251) , by2(0:251) ,     &
              & dy2(0:251)
      eps = 1.0d-14
      if ( Kf==1 .and. m==2*int(m/2) ) kd = 1
      if ( Kf==1 .and. m/=2*int(m/2) ) kd = 2
      if ( Kf==2 .and. m/=2*int(m/2) ) kd = 3
      if ( Kf==2 .and. m==2*int(m/2) ) kd = 4
      call cva2(kd,m,q,a)
      if ( q<=1.0d0 ) then
         qm = 7.5 + 56.1*sqrt(q) - 134.7*q + 90.7*sqrt(q)*q
      else
         qm = 17.0 + 3.1*sqrt(q) - .126*q + .0037*sqrt(q)*q
      endif
      km = int(qm+0.5*m)
      if ( km>=251 ) then
         F1r = dnan()
         D1r = dnan()
         F2r = dnan()
         D2r = dnan()
         return
      endif
      call fcoef(kd,m,q,a,fg)
      ic = int(m/2) + 1
      if ( kd==4 ) ic = m/2
      c1 = exp(-x)
      c2 = exp(x)
      u1 = sqrt(q)*c1
      u2 = sqrt(q)*c2
      call jynb(km+1,u1,nm,bj1,dj1,by1,dy1)
      call jynb(km+1,u2,nm,bj2,dj2,by2,dy2)
      w1 = 0.0d0
      w2 = 0.0d0
      if ( Kc/=2 ) then
         F1r = 0.0d0
         do k = 1 , km
            if ( kd==1 ) then
               F1r = F1r + (-1)**(ic+k)*fg(k)*bj1(k-1)*bj2(k-1)
            elseif ( kd==2 .or. kd==3 ) then
               F1r = F1r + (-1)**(ic+k)*fg(k)                           &
                   & *(bj1(k-1)*bj2(k)+(-1)**kd*bj1(k)*bj2(k-1))
            else
               F1r = F1r + (-1)**(ic+k)*fg(k)                           &
                   & *(bj1(k-1)*bj2(k+1)-bj1(k+1)*bj2(k-1))
            endif
            if ( k>=5 .and. abs(F1r-w1)<abs(F1r)*eps ) exit
            w1 = F1r
         enddo
         F1r = F1r/fg(1)
         D1r = 0.0d0
         do k = 1 , km
            if ( kd==1 ) then
               D1r = D1r + (-1)**(ic+k)*fg(k)                           &
                   & *(c2*bj1(k-1)*dj2(k-1)-c1*dj1(k-1)*bj2(k-1))
            elseif ( kd==2 .or. kd==3 ) then
               D1r = D1r + (-1)**(ic+k)*fg(k)                           &
                   & *(c2*(bj1(k-1)*dj2(k)+(-1)**kd*bj1(k)*dj2(k-1))    &
                   & -c1*(dj1(k-1)*bj2(k)+(-1)**kd*dj1(k)*bj2(k-1)))
            else
               D1r = D1r + (-1)**(ic+k)*fg(k)                           &
                   & *(c2*(bj1(k-1)*dj2(k+1)-bj1(k+1)*dj2(k-1))         &
                   & -c1*(dj1(k-1)*bj2(k+1)-dj1(k+1)*bj2(k-1)))
            endif
            if ( k>=5 .and. abs(D1r-w2)<abs(D1r)*eps ) exit
            w2 = D1r
         enddo
         D1r = D1r*sqrt(q)/fg(1)
         if ( Kc==1 ) return
      endif
      F2r = 0.0d0
      do k = 1 , km
         if ( kd==1 ) then
            F2r = F2r + (-1)**(ic+k)*fg(k)*bj1(k-1)*by2(k-1)
         elseif ( kd==2 .or. kd==3 ) then
            F2r = F2r + (-1)**(ic+k)*fg(k)                              &
                & *(bj1(k-1)*by2(k)+(-1)**kd*bj1(k)*by2(k-1))
         else
            F2r = F2r + (-1)**(ic+k)*fg(k)                              &
                & *(bj1(k-1)*by2(k+1)-bj1(k+1)*by2(k-1))
         endif
         if ( k>=5 .and. abs(F2r-w1)<abs(F2r)*eps ) exit
         w1 = F2r
      enddo
      F2r = F2r/fg(1)
      D2r = 0.0d0
      do k = 1 , km
         if ( kd==1 ) then
            D2r = D2r + (-1)**(ic+k)*fg(k)                              &
                & *(c2*bj1(k-1)*dy2(k-1)-c1*dj1(k-1)*by2(k-1))
         elseif ( kd==2 .or. kd==3 ) then
            D2r = D2r + (-1)**(ic+k)*fg(k)                              &
                & *(c2*(bj1(k-1)*dy2(k)+(-1)**kd*bj1(k)*dy2(k-1))       &
                & -c1*(dj1(k-1)*by2(k)+(-1)**kd*dj1(k)*by2(k-1)))
         else
            D2r = D2r + (-1)**(ic+k)*fg(k)                              &
                & *(c2*(bj1(k-1)*dy2(k+1)-bj1(k+1)*dy2(k-1))            &
                & -c1*(dj1(k-1)*by2(k+1)-dj1(k+1)*by2(k-1)))
         endif
         if ( k>=5 .and. abs(D2r-w2)<abs(D2r)*eps ) exit
         w2 = D2r
      enddo
      D2r = D2r*sqrt(q)/fg(1)
      end



!       **********************************

      subroutine cik01(z,Cbi0,Cdi0,Cbi1,Cdi1,Cbk0,Cdk0,Cbk1,Cdk1)
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
      implicit none
      real(wp) a , a0 , a1 , b , pi , w0
      complex(wp) ca , cb , Cbi0 , Cbi1 , Cbk0 , Cbk1 , Cdi0 , Cdi1 ,    &
               & Cdk0 , Cdk1 , ci , cr , cs , ct , cw , z , z1 , z2 ,   &
               & zr , zr2
      integer k , k0
      dimension a(12) , b(12) , a1(10)
      pi = 3.141592653589793d0
      ci = (0.0d0,1.0d0)
      a0 = abs(z)
      z2 = z*z
      z1 = z
      if ( a0==0.0d0 ) then
         Cbi0 = (1.0d0,0.0d0)
         Cbi1 = (0.0d0,0.0d0)
         Cdi0 = (0.0d0,0.0d0)
         Cdi1 = (0.5d0,0.0d0)
         Cbk0 = (1.0d+300,0.0d0)
         Cbk1 = (1.0d+300,0.0d0)
         Cdk0 = -(1.0d+300,0.0d0)
         Cdk1 = -(1.0d+300,0.0d0)
         return
      endif
      if ( dble(z)<0.0 ) z1 = -z
      if ( a0<=18.0 ) then
         Cbi0 = (1.0d0,0.0d0)
         cr = (1.0d0,0.0d0)
         do k = 1 , 50
            cr = 0.25d0*cr*z2/(k*k)
            Cbi0 = Cbi0 + cr
            if ( abs(cr/Cbi0)<1.0d-15 ) exit
         enddo
         Cbi1 = (1.0d0,0.0d0)
         cr = (1.0d0,0.0d0)
         do k = 1 , 50
            cr = 0.25d0*cr*z2/(k*(k+1))
            Cbi1 = Cbi1 + cr
            if ( abs(cr/Cbi1)<1.0d-15 ) exit
         enddo
         Cbi1 = 0.5d0*z1*Cbi1
      else
         data a/0.125d0 , 7.03125d-2 , 7.32421875d-2 ,                  &
            & 1.1215209960938d-1 , 2.2710800170898d-1 ,                 &
            & 5.7250142097473d-1 , 1.7277275025845d0 ,                  &
            & 6.0740420012735d0 , 2.4380529699556d01 ,                  &
            & 1.1001714026925d02 , 5.5133589612202d02 ,                 &
            & 3.0380905109224d03/
         data b/ - 0.375d0 , -1.171875d-1 , -1.025390625d-1 ,           &
            & -1.4419555664063d-1 , -2.7757644653320d-1 ,               &
            & -6.7659258842468d-1 , -1.9935317337513d0 ,                &
            & -6.8839142681099d0 , -2.7248827311269d01 ,                &
            & -1.2159789187654d02 , -6.0384407670507d02 ,               &
            & -3.3022722944809d03/
         k0 = 12
         if ( a0>=35.0 ) k0 = 9
         if ( a0>=50.0 ) k0 = 7
         ca = exp(z1)/sqrt(2.0d0*pi*z1)
         Cbi0 = (1.0d0,0.0d0)
         zr = 1.0d0/z1
         do k = 1 , k0
            Cbi0 = Cbi0 + a(k)*zr**k
         enddo
         Cbi0 = ca*Cbi0
         Cbi1 = (1.0d0,0.0d0)
         do k = 1 , k0
            Cbi1 = Cbi1 + b(k)*zr**k
         enddo
         Cbi1 = ca*Cbi1
      endif
      if ( a0<=9.0 ) then
         cs = (0.0d0,0.0d0)
         ct = -log(0.5d0*z1) - 0.5772156649015329d0
         w0 = 0.0d0
         cr = (1.0d0,0.0d0)
         do k = 1 , 50
            w0 = w0 + 1.0d0/k
            cr = 0.25d0*cr/(k*k)*z2
            cs = cs + cr*(w0+ct)
            if ( abs((cs-cw)/cs)<1.0d-15 ) exit
            cw = cs
         enddo
         Cbk0 = ct + cs
      else
         data a1/0.125d0 , 0.2109375d0 , 1.0986328125d0 ,               &
            & 1.1775970458984d01 , 2.1461706161499d02 ,                 &
            & 5.9511522710323d03 , 2.3347645606175d05 ,                 &
            & 1.2312234987631d07 , 8.401390346421d08 ,                  &
            & 7.2031420482627d10/
         cb = 0.5d0/z1
         zr2 = 1.0d0/z2
         Cbk0 = (1.0d0,0.0d0)
         do k = 1 , 10
            Cbk0 = Cbk0 + a1(k)*zr2**k
         enddo
         Cbk0 = cb*Cbk0/Cbi0
      endif
      Cbk1 = (1.0d0/z1-Cbi1*Cbk0)/Cbi0
      if ( dble(z)<0.0 ) then
         if ( dimag(z)<0.0 ) Cbk0 = Cbk0 + ci*pi*Cbi0
         if ( dimag(z)>0.0 ) Cbk0 = Cbk0 - ci*pi*Cbi0
         if ( dimag(z)<0.0 ) Cbk1 = -Cbk1 + ci*pi*Cbi1
         if ( dimag(z)>0.0 ) Cbk1 = -Cbk1 - ci*pi*Cbi1
         Cbi1 = -Cbi1
      endif
      Cdi0 = Cbi1
      Cdi1 = Cbi0 - 1.0d0/z*Cbi1
      Cdk0 = -Cbk1
      Cdk1 = -Cbk0 - 1.0d0/z*Cbk1
      end

!       **********************************

      subroutine cpsi(x,y,Psr,Psi)
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
      implicit none
      real(wp) a , ct2 , pi , Psi , Psr , ri , rr , th , tm ,   &
                     & tn , x , x0 , x1 , y , y1 , z0 , z2
      integer k , n
      dimension a(8)
      data a/ - .8333333333333d-01 , .83333333333333333d-02 ,           &
         & -.39682539682539683d-02 , .41666666666666667d-02 ,           &
         & -.75757575757575758d-02 , .21092796092796093d-01 ,           &
         & -.83333333333333333d-01 , .4432598039215686d0/
      pi = 3.141592653589793d0
      if ( y==0.0d0 .and. x==int(x) .and. x<=0.0d0 ) then
         Psr = 1.0d+300
         Psi = 0.0d0
      else
         x1 = x
         y1 = y
         if ( x<0.0d0 ) then
            x = -x
            y = -y
         endif
         x0 = x
         n = 0
         if ( x<8.0d0 ) then
            n = 8 - int(x)
            x0 = x + n
         endif
         th = 0.0d0
         if ( x0==0.0d0 .and. y/=0.0d0 ) th = 0.5d0*pi
         if ( x0/=0.0d0 ) th = atan(y/x0)
         z2 = x0*x0 + y*y
         z0 = sqrt(z2)
         Psr = log(z0) - 0.5d0*x0/z2
         Psi = th + 0.5d0*y/z2
         do k = 1 , 8
            Psr = Psr + a(k)*z2**(-k)*cos(2.0d0*k*th)
            Psi = Psi - a(k)*z2**(-k)*sin(2.0d0*k*th)
         enddo
         if ( x<8.0d0 ) then
            rr = 0.0d0
            ri = 0.0d0
            do k = 1 , n
               rr = rr + (x0-k)/((x0-k)**2.0d0+y*y)
               ri = ri + y/((x0-k)**2.0d0+y*y)
            enddo
            Psr = Psr - rr
            Psi = Psi + ri
         endif
         if ( x1<0.0d0 ) then
            tn = tan(pi*x)
            tm = dtanh(pi*y)
            ct2 = tn*tn + tm*tm
            Psr = Psr + x/(x*x+y*y) + pi*(tn-tn*tm*tm)/ct2
            Psi = Psi - y/(x*x+y*y) - pi*tm*(1.0d0+tn*tn)/ct2
            x = x1
            y = y1
         endif
      endif
      end

!       **********************************

      subroutine sphy(n,x,Nm,Sy,Dy)
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
      implicit none
      real(wp) Dy , f , f0 , f1 , Sy , x
      integer k , n , Nm
      dimension Sy(0:n) , Dy(0:n)
      Nm = n
      if ( x<1.0d-60 ) then
         do k = 0 , n
            Sy(k) = -1.0d+300
            Dy(k) = 1.0d+300
         enddo
         return
      endif
      Sy(0) = -cos(x)/x
      f0 = Sy(0)
      Dy(0) = (sin(x)+cos(x)/x)/x
      if ( n<1 ) return
      Sy(1) = (Sy(0)-sin(x))/x
      f1 = Sy(1)
      do k = 2 , n
         f = (2.0d0*k-1.0d0)*f1/x - f0
         Sy(k) = f
         if ( abs(f)>=1.0d+300 ) exit
         f0 = f1
         f1 = f
      enddo
      Nm = k - 1
      do k = 1 , Nm
         Dy(k) = Sy(k-1) - (k+1.0d0)*Sy(k)/x
      enddo
      end

!       **********************************

      subroutine jelp(u,Hk,Esn,Ecn,Edn,Eph)
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
      implicit none
      real(wp) a , a0 , b , b0 , c , d , dn , Ecn , Edn , Eph , &
                     & Esn , Hk , pi , r , sa , t , u
      integer j , n
      dimension r(40)
      pi = 3.14159265358979d0
      a0 = 1.0d0
      b0 = sqrt(1.0d0-Hk*Hk)
      do n = 1 , 40
         a = (a0+b0)/2.0d0
         b = sqrt(a0*b0)
         c = (a0-b0)/2.0d0
         r(n) = c/a
         if ( c<1.0d-7 ) exit
         a0 = a
         b0 = b
      enddo
      dn = 2.0d0**n*a*u
      d = 0.0d0
      do j = n , 1 , -1
         t = r(j)*sin(dn)
         sa = atan(t/sqrt(abs(1.0d0-t*t)))
         d = .5d0*(dn+sa)
         dn = d
      enddo
      Eph = d*180.0d0/pi
      Esn = sin(d)
      Ecn = cos(d)
      Edn = sqrt(1.0d0-Hk*Hk*Esn*Esn)
      end

    end module specfun_module