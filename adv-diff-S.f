c----------------------------------------------------------------------2
      program advdif
c
ccccccccccccccccccccccccccccccccccccccccc
c                                       c
c     Author: C. S. Ng                  c
c                                       c
c     Last change: 2021-07-14           c
c                                       c
ccccccccccccccccccccccccccccccccccccccccc
c
c  This program solve the temperature based on a heating model
c  with both advection and diffusion, as a two-point
c  boundary value problem.
c  For Saturn parameters
c
      integer amax
      PARAMETER (amax=1000000)
      real eta,xi0,beta,p,k0,kf,Lm,Lw,L0,Lf,T0,Tf,vs,kd
      real L,dL,dL2,be3,L1,dk
      real Ta(0:amax),k(0:amax),v(0:amax)
      real a(amax),b(amax),c(amax),f(amax),t(amax)
      real d2T,cc,dd,ss
      integer nmax,i
c
      namelist/in/eta,xi0,beta,p,k0,kf,L0,Lf,Lm,Lw,T0,Tf,vs,nmax
      read(*,in)
      write(*,in)
      write(*,*)eta
      write(*,*)xi0
      write(*,*)beta
      write(*,*)p
      write(*,*)k0
      write(*,*)kf
      write(*,*)L0
      write(*,*)Lf
      write(*,*)Lm
      write(*,*)Lw
      write(*,*)T0
      write(*,*)Tf
      write(*,*)vs
      write(*,*)nmax
      dL = (Lf-L0)/nmax
      dk = (kf-k0)/nmax
      write(*,*) dL
      write(*,*)

c
c eta is a factor proportional to mdot
c xi0 is a factor for the q term
c beta is the index for the density
c p is the index for the q term
c k0 is a constant for the k(L) function as a test
c L0 is the smaller value of the L range
c Lf is the larger value of the L range
c T0 is temperature at L0
c Tf is temperature at Lf
c dT0 is dT/dL at L0
c nmax is the number of points along L (< 1000000)
c L is the normalized radial position
c dL is the increment in the integration
c
c      open(unit=12,file='bigger-2.txt',STATUS='old')
c
c      read(12,10)
c      do i = 0,190
c        read(12,10) rx(i),jx(i),dx(i)
c      enddo
c10    format(e14.8,x,2(e15.8,x))
c
      dL2 = dL*dL
      be3 = 2.*beta/3.
      k(0) = k0
      kd = k0 - kf
      do i = 1,nmax
        L = L0 + i*dL
        k(i) = kf + 0.5*kd*erfc((L-Lm)/Lw)
      enddo
c
      do i = 1,nmax-1
        L = L0 + i*dL
        a(i) = 1. + 0.5*dL*(eta*L**(beta-k(i)-1.)+beta-k(i)-6.)/L
        b(i) = -2.-dL2*(be3*(eta*L**(beta-k(i)-1.)-beta+4.)
     .              -(6.-beta)*(k(i)-beta+3.))/L**2
        c(i) = 1. - 0.5*dL*(eta*L**(beta-k(i)-1.)+beta-k(i)-6.)/L
        f(i) = -xi0*L**(beta-k(i)+p)*dL2
      enddo
      f(1) = f(1) - a(1)*T0
      f(nmax-1) = f(nmax-1) - c(nmax-1)*Tf
c
      call tridia(a,b,c,f,nmax-1,t)
c
      t(nmax) = Tf
      Ta(0) = T0
      Ta(1) = Ta(0) + (xi0*L0**(p+2.)/eta - be3*Ta(0)/L0)*dL
      do i = 2, nmax
        L = L0 + i*dL
        L1 = L - dL
        Ta(i) = Ta(i-2) + (xi0*L1**(p+2.)/eta
     .          - be3*Ta(i-1)/L1)*dL*2.
      enddo
c
      v(0) = 0.
      do i = 1, nmax
        L = L0 + i*dL
        v(i) = vs*(L0**(k0-beta+1.) - L**(k(i)-beta+1.))*L**(beta-2.)
      enddo
c
      write(*,90) "L                        ",
     .            "T                        ",
     .            "Ta                       ",
     .            "v                        "
c
      write(*,100) L0,T0,Ta(0),v(0),k(0)
      do i = 1,nmax
        write(*,100) L0+i*dL,t(i),Ta(i),v(i),k(i)
      enddo
c
90    format(10a25)
100   format(10(e24.14,x))
c
1000  continue
      STOP
      END     
c----------------------------------------------------------------------2
c tridiagonal matrix algorithm
c
c copied from Y. Jaluria K. Torrance ,computational heat transfer pp336
c
c a,b and c are the three elements in each row. With b at the diagonal
c f is the constant on the right-hand side of each equation. n is the
c number of equations and t is the variable to be computed.
c
      subroutine tridia(a,b,c,f,n,t)
c
      integer nmax
      parameter (nmax = 10000)
      integer n,i,j
c
      real a(nmax),b(nmax),c(nmax),f(nmax),t(nmax),d
c
c Elimiante the a's by Gaussian elimination and determine
c the new coefficients.
c
      do i = 2,n
        d = a(i)/b(i-1)
        b(i) = b(i) - c(i-1)*d
        f(i) = f(i) - f(i-1)*d
      enddo
c
c back subsititution
c
      t(n) = f(n)/b(n)
      do i = 1,n-1
        j = n-i
        t(j) = (f(j) - c(j)*t(j+1))/b(j)
      enddo
      return
      end
c----------------------------------------------------------------------2
      subroutine rhs(d2t,t,dt,L,eta,beta,k,xi0,p)
c
      real d2t,c,d,t,dt,L
      real eta,beta,k,xi0,p,be3
c
      be3 = 2.*beta/3.
      d2t = dt*(eta*L**(beta-k)+beta-k-5.)/L
     . + t*(be3*(eta*L**(beta-k)-beta+3.)-(5.-beta)*(k
     . - beta+2.))/L**2 - xi0*L**(beta-k+p)
      return
      end
c----------------------------------------------------------------------2
      subroutine rhs1(c,d,s,L,eta,beta,k,xi0,p)
c
      real c,d,s,L
      real eta,beta,k,xi0,p,be3
c
      be3 = 2.*beta/3.
      c = (eta*L**(beta-k)+beta-k-5.)/L
      d = (be3*(eta*L**(beta-k)-beta+3.)-(5.-beta)*(k-beta+2.))/L**2
      s = - xi0*L**(beta-k+p)
      return
      end
c----------------------------------------------------------------------2
      FUNCTION ERFC(X)
c
c From Numerical Recipes
c
      real erfc,x,gammp,gammq
c
      IF(X.LT.0.)THEN
        ERFC=1.+GAMMP(.5,X**2)
      ELSE
        ERFC=GAMMQ(.5,X**2)
      ENDIF
      RETURN
      END
c----------------------------------------------------------------------2
      FUNCTION GAMMP(A,X)
c
c From Numerical Recipes
c
c      real gammp,a,x,gser,gln,gcf,gammcf
      real gammp,a,x,gln,gammcf
c
      IF (X.LT.0..OR.A.LE.0.) then
        write(*,*) "error in GAMMP"
        return
      endif
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMMP,A,X,GLN)
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.-GAMMCF
      ENDIF
      RETURN
      END
c----------------------------------------------------------------------2
      FUNCTION GAMMQ(A,X)
c
c From Numerical Recipes
c
c      real gammq,a,x,gser,gamser,gln,gcf
      real gammq,a,x,gamser,gln
c
      IF(X.LT.0..OR.A.LE.0.) then
        write(*,*) "error in GAMMQ"
        return
      endif
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.-GAMSER
      ELSE
        CALL GCF(GAMMQ,A,X,GLN)
      ENDIF
      RETURN
      END
c----------------------------------------------------------------------2
      SUBROUTINE GSER(GAMSER,A,X,GLN)
c
c From Numerical Recipes
c
      real gamser,a,x,gln,eps,gammln,ap,sum,del
      integer itmax,n
c
      PARAMETER (ITMAX=1000,EPS=3.E-12)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.) then
          write(*,*)"error in GSER (x < 0)"
          gamser = 0.
          return
        endif
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      write(*,*) "error in GSER, A too large, ITMAX too small"
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END
c----------------------------------------------------------------------2
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
c
c From Numerical Recipes
c
      real gammcf,a,x,gln,eps,gammln,gold,a0,a1,b0,b1,fac,an,ana,
     .     anf,g
      integer itmax,n
c
      PARAMETER (ITMAX=1000,EPS=3.E-12)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=N*1.0
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      write(*,*) "error in GCF, A too large, ITMAX too small"
1     GAMMCF=EXP(-X+A*ALOG(X)-GLN)*G
      RETURN
      END
c----------------------------------------------------------------------2
      FUNCTION GAMMLN(XX)
c
c From Numerical Recipes
c
      real gammln,xx
c
      REAL COF(6),STP,X,TMP,SER,y
      integer j
      DATA COF,STP/76.18009172947146e0,-86.50532032941677e0,
     *    24.01409824083091e0,
     *    -1.231739572450155e0,.1208650973866179e-2,-.5395239384953e-5,
     *     2.5066282746310005e0/
      X=XX
      y=x
      TMP=X+5.5
      TMP=(X+0.5)*LOG(TMP)-TMP
      SER=1.000000000190015
      DO 11 J=1,6
        y=y+1.
        SER=SER+COF(J)/y
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER/x)
      RETURN
      END
c----------------------------------------------------------------------2







