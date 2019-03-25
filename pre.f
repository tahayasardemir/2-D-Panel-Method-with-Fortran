C----------------------------------------------------------------------
C
C.. SPOT; Steady POTential flow solver for airfoils
C.. using a low-order (zeroth) panel method  
C
      program MAIN
      parameter( mx = 301, nx=91 )
      real x(mx),y(mx),xm(mx),ym(mx),costhe(mx),sinthe(mx)
      real xg(mx,nx),yg(mx,nx)
      real Sinf(mx,mx),Ginf(mx,mx),A(mx+1,mx+2)
      real Sigma(mx),Gamma

c..Get input data and Airfoil/panel coordinates
      call INPUT(mx,npanel,x,y,xm,ym,costhe,sinthe)

c..Setup the influence coefficients
      call INFCOEF(mx,npanel,x,y,xm,ym,costhe,sinthe,Sinf,Ginf,A)

c..Solve the system of equations
      call GAUSS(mx,npanel,A)

c..Retrieve the solution
      neqn = npanel+1
      do i = 1, npanel
         Sigma(i) = A(i,neqn+1)
      enddo
      Gamma       = A(neqn,neqn+1)

c..Compute Cp and the aerodynamic load coeffs. Cl,Cd, Cm
      call LOADS(mx,npanel,x,y,xm,ym,costhe,sinthe,Gamma,Sigma,Ginf,Sinf
     >)

      call OGRID(npanel,x,y,xg,yg)

      call QOUT(npanel,x,y,costhe,sinthe,xg,yg,Sigma,Gamma,xm,ym,Sinf,
     >Ginf)

      stop
      end
C----------------------------------------------------------------------
      subroutine INPUT(mx,npanel,x,y,xm,ym,costhe,sinthe)
      real x(mx),y(mx),xm(mx),ym(mx),costhe(mx),sinthe(mx)
      character*10 argv

c---- Retrieve parameters from the command line in Linux.
c     call GETARG(1,argv)
c     read(argv,*) naca
c     call GETARG(2,argv)
c     read(argv,*) npanel
c     call GETARG(3,argv)
c     read(argv,*) alpha

c---  Or set them manually
c     naca   = 12
c     npanel = 100
c     alpha  = 0

c---  Or read from the standart input
      print*,' '
      print*,'   Enter NACA Profile (0012), # of Panels (even), AoA : '
      read(5,*)  naca, npanel, alpha

c---- Sanitize input parameters
      alpha = min(90., alpha)
      alpha = max(-90., alpha)
      alpha = alpha*acos(-1.)/180
      if (mod(npanel,2) .ne. 0 ) npanel = npanel - 1
      npanel = min(mx,npanel) 
      npanel = max(20,npanel)

c---- Generate/read the airfoil coordinates (x,y arrays)
      if( naca .ne. 0) then
        call NACAFOIL(mx,npanel,naca,x,y)
      else
        call READFOIL(mx,npanel,x,y)
      endif

c----Rotate airfoil  if alpha != 0
      if ( alpha .ne. 0) then
        cosa = cos(alpha)
        sina = sin(alpha)
        do i = 1, npanel+1
           xt   = x(i)
           x(i) = xt*cosa + y(i)*sina
           y(i) =-xt*sina + y(i)*cosa
        enddo
       endif

c---write the the airfoil/panel x and y coordinates
      open(2,file='airfoil.dat')
      write(2,'(f8.6,2x,e12.4)') (x(n),y(n),n=1,npanel+1)
      close(2)

c---- Compute the panel properties
      do i = 1, npanel
c-------Control points at the mid-panel locations.
        xm(i) = 0.5*(x(i+1) + x(i))
        ym(i) = 0.5*(y(i+1) + y(i))
c----Panel size and slopes/angles
        dx = x(i+1) - x(i)
        dy = y(i+1) - y(i)
        plen = sqrt(dx**2 + dy**2)
        sinthe(i) = dy/plen
        costhe(i) = dx/plen
      enddo

      return
      end

C----------------------------------------------------------------------
      subroutine READFOIL(mx,npanel,x,y)
      real x(mx),y(mx)
      logical ok
      character fname*36
      data fname /' '/
 
      do while (fname .eq. ' ')
         print*, ''
         print*, '   Enter the filename for the AIRFOIL coordinates: '
         read(*,'(a)') fname
      enddo
      inquire(FILE=fname,EXIST=ok)
      if( .not. ok ) then
         print*, fname,' does NOT EXIST!'
         stop
      endif
      open(1,file=fname,form='formatted')
      read(1,*,end=10,err=10)(x(n),y(n),n=1,mx)
   10 npanel = n-2

      return
      end
c------------------------------------------------------------------------------

      subroutine NACAFOIL(mx,npanel,naca,x,y)
      real x(mx),y(mx)
      pi = acos(-1.)

c..Compute coordinates of NACA 4/5 digit airfoils

c..Decompose the NACA number to determine airfoil properties
       if (( naca .gt. 25999 ).or.( naca .lt. 1 )) naca = 12
       ieps = naca / 1000
       iptmax = naca / 100 - 10 * ieps
       itau = naca - 1000 * ieps - 100 * iptmax
 
c..Set the coefficients.
       epsmax = ieps   * 0.01
       ptmax  = iptmax * 0.1
       tau    = itau   * 0.01
 
c..Error correction for bogus NACA numbers.
       if (( naca .le. 9999 ) .and. ( epsmax .gt. 0 ) .and.
     >     ( ptmax .eq. 0 )) ptmax = 0.1
 
c..If NACA 5 digit coding is used, make neccessary changes.
       if ( ieps .ge. 10 ) then
         if ( ieps .eq. 21 ) then
           ptmax = 0.0580
           ak1 = 361.4
         elseif ( ieps .eq. 22 ) then
           ptmax = 0.1260
           ak1 = 51.64
         elseif ( ieps .eq. 23 ) then
           ptmax = 0.2025
           ak1 = 15.957
         elseif ( ieps .eq. 24 ) then
           ptmax = 0.2900
           ak1 = 6.643
         elseif ( ieps .eq. 25 ) then
           ptmax = 0.3910
           ak1 = 3.230
         endif
         epsmax = ak1 * ptmax**3 / 6
       endif
 
c..initialize indexing for lower surface.
       nlower = npanel / 2
       nupper = npanel - nlower
 
c..Loop over the lower surface.
       index = 0
       do n = 1, nlower
         index = index + 1
         xc = 0.5*( 1 + cos( pi*(n-1)/nlower) )
         call NACA45(naca,tau,epsmax,ptmax,xc,thick,camber,beta)
         x(index) = xc     + thick * sin( beta )
         y(index) = camber - thick * cos( beta )
       enddo
 
c..Loop over the upper surface.
       do n = 1, nupper
         index = index + 1
         xc = 0.5*( 1 - cos(pi*(n-1)/nupper) )
         call NACA45(naca,tau,epsmax,ptmax,xc,thick,camber,beta)
         x(index) = xc     - thick * sin( beta )
         y(index) = camber + thick * cos( beta )
       enddo

c..Set the last point.
       x(index+1) = x(1)
       y(index+1) = y(1)
 
       return
       end

c------------------------------------------------------------------------------

       subroutine NACA45(naca,tau,epsmax,ptmax,xc,thick,camber,beta)
c..Compute the thickness, camber, and angular location of an airfoil point.

c..Compute the thickness
       if ( xc .lt. 1.0E-10 ) then
         thick = 0.0    !..Thickness is corrected when xc is very small.
       else
         thick = tau * 5 * ( 0.2969 * SQRT(xc)
     >               - xc * ( 0.1260
     >               + xc * ( 0.3537
     >               - xc * ( 0.2843
     >               - xc * 0.1015))) )
       endif
 
c..Compute the camber
       if ( epsmax .eq. 0.0 ) then
c..For NACA 4-digit symmetrical arfoils.
         camber = 0.0
         beta = 0.0
       else
         if ( naca .gt. 9999 ) then
c..For NACA 5 digit airfoils
c..Ptmax = m and epsmax = (k_1*m^3)/6 from Abbott and Doenhoff.
           if ( xc .gt. ptmax ) then
             camber = epsmax * ( 1.0 - xc )
             dcamdx = - epsmax
           else
             w = xc / ptmax
             camber = epsmax * ( w**3 - 3*w**2 +(3.-ptmax)*w)
             dcamdx = epsmax/ptmax*(3*w**2 - 6*w + ( 3.0-ptmax))
           endif
         else

c..For NACA 4 digit airfoils.
           if ( xc .gt. ptmax ) then
             camber = epsmax / ( 1.0 - ptmax )**2
     >              * ( 1. + xc - ptmax * 2 ) * ( 1. - xc )
             dcamdx = epsmax * 2 / ( 1.0 - ptmax )**2 * ( ptmax - xc )
           else
             camber = epsmax / ptmax**2 * ( ptmax*2 - xc ) * xc
             dcamdx = epsmax * 2 / ptmax**2  * ( ptmax - xc )
           endif
         endif
 
         beta = atan(dcamdx)
       endif
 
       return
       end
c-----------------------------------------------------------------------

      SUBROUTINE OGRID(np,x,y,xg,yg)
      parameter( mx=301,nx=61 )
      dimension x(mx),y(mx),R(mx),THETA(mx),xg(mx,nx),yg(mx,nx)
      complex Z1,Z2,Z3,Z4,Z5,Z6,Z7
      data jmax/61/, xn,yn/0.008,0./
C
C     ALGEBRAIC O-GRID GENERATOR BASED ON KARMAN-TREFFZ MAPPING
C     ******************************************************
C     Originally written by L.N. Sankar @ Georgia Tech
C
      N     = np+1
      PI    = ACOS(-1.0)
      SLOP1 =  ATAN2((Y(1)-Y(2)),(X(1)-X(2)))
      SLOP2 = -ATAN2((Y(N)-Y(N-1)),(X(N)-X(N-1)))             
      EPS   = ABS(SLOP1 + SLOP2)
      POWER = 2.0 - EPS / PI
      POWER = 1.0 / POWER 
      XT    = 0.5 * (X(1) + X(N))
      YT    = 0.5 * (Y(1) + Y(N)) + 0.000001
      Z1    = CMPLX(XN,YN)
      Z2    = CMPLX(XT,YT)
      Z7    = 0.5 * (Z1+Z2)
      U     = 1.0
      V     = 0.0
      ANGL  = 2.*PI

      DO 30 I = 1, N
      Z3    = CMPLX(X(I),Y(I))
      Z4    = (Z3 - Z2) / (Z3-Z1)
      Z4    = CLOG(Z4) * CMPLX(POWER,0.0)
      Z4    = CEXP(Z4)
      Z5    = CMPLX(1.0,0.0)
      Z4    = (Z5+Z4) * Z7 / (Z5-Z4)
      X4    = REAL(Z4)
      Y4    = AIMAG(Z4)
      R(I)  = SQRT(X4*X4+Y4*Y4)
      A     = X4
      B     = Y4
      ANGL  = ANGL + ATAN2((U*B-V*A),(U*A+V*B))
      U     = A
      V     = B 
      THETA(I)  = ANGL
   30 CONTINUE

C
C     GENERATE THE GRID
C
      DR      = 0.96 / (JMAX-2)
      DO 40 J = 2, JMAX
      ETA     = (J-2) * DR
      ETA1    = 1. / (1. - ETA)
      DO 40 I = 1, N-1
      X4      = R(I) * ETA1 * COS(THETA(I))
      Y4      = R(I) * ETA1 * SIN(THETA(I))
      Z4      = CMPLX(X4,Y4)
      Z4      = (Z4-Z7) / (Z4+Z7)
      Z4      = CMPLX(1.0/POWER,0.0) * CLOG(Z4)
      Z4      = CEXP(Z4)
      Z4      = (Z2 - Z4 * Z1) / (CMPLX(1.0,0.0) - Z4)
      XG(I,J) = REAL(Z4)
      YG(I,J) = AIMAG(Z4)
   40 CONTINUE

C..set grid coordinates at j = 1 and i = n
      do i = 1, n
        xg(i,1) = x(i)
        yg(i,1) = y(i)
      enddo
      do j = 1, jmax
        xg(n,j) = xg(1,j)
        yg(n,j) = yg(1,j)
      enddo

      open(unit=9,file='grid.plt')
      write(9,*) 'VARIABLES = "X", "Y"'
      write(9,*) 'ZONE F=POINT, I=',n,', J=',jmax
      write(9,'(2f12.7)') ((xg(i,j),yg(i,j),i=1,n),j=1,jmax)
      close(9)

      RETURN
      END
c-----------------------------------------------------------------------

      subroutine LOADS(mx,npanel,x,y,xm,ym,costhe,sinthe,
     +                 Gamma,Sigma,Ginf,Sinf)
      real x(mx),y(mx),xm(mx),ym(mx),costhe(mx),sinthe(mx)
      real Sigma(mx),Sinf(mx,mx),Ginf(mx,mx),Cp(mx)

c..Find V_tangential and Cp at the mid-point of i-th panel.
c..use influence coefficients already stored in Sinf and Ginf arrays
      do i = 1, npanel
c..Free Stream contribution
         Vtan = costhe(i) 
c..Add the contributions from source and vortex panels 
         do j = 1, npanel
           Vtan = Vtan - Sigma(j)*Ginf(i,j) + Gamma*Sinf(i,j)
         enddo
         cp(i) = 1.0 - Vtan**2
      enddo

c..write the xmid and cp values for all the panels
      print*, ''
      print*, '  Cp distribution is written in   cp.dat '
      open(1,file='cp.dat')
      write(1,'(f7.4,1x,e12.4)') (xm(n),cp(n),n=1,npanel)
      close(1)

      print*, ''
      print*, '*** Compute Cl, Cd and Cm@x=0.25c in subroutine LOADS'
      print*, ''

c..integrate Cp  to obtain Cl, Cd and Cm about x=0.25
      Cl = 0
      Cd = 0
      Cm = 0.
      do n = 1,npanel
            Cl = Cl + cp(n)*(x(n+1)-x(n))
            Cd = Cd + cp(n)*(y(n+1)-y(n))
            Cm = Cm + cp(n)*((x(n+1)-x(n))*(xm(n)-0.25)+(y(n+1)-y(n))*
     &(ym(n)))
      enddo
c    ...
c..write Cl, Cd and Cm values..
      print*, ' Aerodynamic load Coefs.; Cl, Cd, Cm :',Cl,Cd,Cm

      return
      end
C-----------------------------------------------------------------------
      subroutine INFCOEF(mx,npanel,x,y,xm,ym,costhe,sinthe,Sinf,Ginf,A)
      real x(mx),y(mx),xm(mx),ym(mx),costhe(mx),sinthe(mx)
      real Sinf(mx,mx),Ginf(mx,mx),A(mx+1,mx+2)
 
       neqn   = npanel + 1
       pi     = acos(-1.0)
       pi2inv = 1./(2*pi)
 
c..Initialize coefs.
       do i = 1, neqn
       do j = 1, neqn
         A(i,j) = 0.0
       enddo
       enddo

       DO i = 1,npanel
c..Add the influence of j^th panel on the Vn of i^th.
       DO j = 1,npanel
          if ( j .eq. i ) then      !..influence on itself
             flog = 0.0
             ftan = pi
          else
             dxj  = xm(i) - x(j)
             dxjp = xm(i) - x(j+1)
             dyj  = ym(i) - y(j)
             dyjp = ym(i) - y(j+1)
             flog = 0.5*alog((dxjp**2 + dyjp**2)/(dxj**2  + dyj**2))
             ftan = atan2( dyjp*dxj - dxjp*dyj , dxjp*dxj + dyjp*dyj )
          endif
          ctimtj = costhe(i) * costhe(j) + sinthe(i) * sinthe(j)
          stimtj = sinthe(i) * costhe(j) - costhe(i) * sinthe(j)
c..Store Ginf and Sinf values  to be used in Cp computations
          Sinf(i,j) = pi2inv * ( ftan * ctimtj + flog * stimtj )
          Ginf(i,j) = pi2inv * ( flog * ctimtj - ftan * stimtj )

          A(i,j)    = Sinf(i,j)                !..Source influence
          A(i,neqn) = A(i,neqn) + Ginf(i,j)    !..Vortex influence

c..Find the influence of j^th panel on the Vtan of 1^st and the last panels
          if (( i .eq. 1 ) .or. ( i .eq. npanel )) then
             A(neqn,j)    = A(neqn,j)    - Ginf(i,j)
             A(neqn,neqn) = A(neqn,neqn) + Sinf(i,j)
          endif
       ENDDO
c..Fill in the RHS with the free stream contribution
       A(i,neqn+1) = sinthe(i) 
       ENDDO

c..Fill in the RHS  of 1st and the last panels with free stream contribution
       A(neqn,neqn+1) = -(costhe(1) + costhe(npanel))
 
       return
       end
c-------------------------------------------------------------------
      subroutine GAUSS(mx,npanel,A)
c..Performs Gaussian elimination on matrix A
      real A(mx+1, mx+2)
      
      neqn = npanel+1
      krhs = neqn + 1
c..Do full matrix elimination sequence.
      do i = 2, neqn
         im = i - 1
c..Eliminate the (i-1)th unknown from i-th through neqn-th equations.
         do j = i, neqn
            r = A(j,im) / A(im,im)
            do k = i, krhs
               A(j,k) = A(j,k) - r * A(im,k)
            enddo
         enddo
      enddo
c..Back subtitution.
      k = krhs
      A(neqn,k) = A(neqn,k) / A(neqn,neqn)
      do l = 2, neqn
         i = neqn + 1 - l
         do j = i+1, neqn
            A(i,k) = A(i,k) - A(i,j) * A(j,k)
         enddo
         A(i,k) = A(i,k) / A(i,i)
      enddo

      return
      end
c-------------------------------------------------------------------
      subroutine QOUT(npanel,x,y,costhe,sinthe,xg,yg,Sigma,Gamma,xm,
     +                ym,Sinf,Ginf)
      parameter( mx = 301, jmax=61 )
      real x(mx),y(mx),costhe(mx),sinthe(mx),xg(mx,jmax),yg(mx,jmax)
      real Sigma(mx),u(mx,jmax),v(mx,jmax),cp(mx,jmax),xm(mx),ym(mx)
      real Sinf(mx,mx),Ginf(mx,mx),Sinfgrid(mx,mx),Ginfgrid(mx,mx)

      imax    = npanel+1
      pi      = acos(-1.0)
      pi2inv  = 1./(2*pi)
      yg(1,1) = 0.
      yg(npanel,npanel) = 0.
c..Evaluate the velocity and pressure field
      do j=2,jmax ! j=1 and J=2 are overlapping so we are starting from 2
      do i=1,imax
         xc = xg(i,j) 
         yc = yg(i,j)
         if (i.eq.imax) then ! Consider grid as saperate panels and find their theta
         plengrid = sqrt((xg(1,j)-xg(i,j))**2+(yg(1,j)-yg(i,j))**2)
         costhegrid = (xg(1,j)-xg(i,j))/plengrid
         sinthegrid = (yg(1,j)-yg(i,j))/plengrid
         else
         plengrid = sqrt((xg(i+1,j)-xg(i,j))**2+(yg(i+1,j)-yg(i,j))**2)
         costhegrid = (xg(i+1,j)-xg(i,j))/plengrid
         sinthegrid = (yg(i+1,j)-yg(i,j))/plengrid
         endif    
         u(i,j) = 1.    !..free stream contribution
         v(i,j) = 0.
         u_l = 0. ! initiate local velocities
         v_l = 0. 

         if (i.le.npanel+1.and.j.eq.2) then ! First loop in on the airfoil. Calculate global velocities
c..use influence coefficients already stored in Sinf and Ginf arrays
            Vtan = 0.
c..Add the contributions from source and vortex panels 
            do n = 1, npanel
            Vtan = Vtan - Sigma(n)*Ginf(i,n) + Gamma*Sinf(i,n)
            enddo
            u(i,j) = u(i,j) + Vtan * costhe(i)
            v(i,j) = v(i,j) + Vtan * sinthe(i)
         else ! Other loops are on the grids
            do k=1,npanel  !..add source and vortex contr. from each panel in X&Y
            dxj  = xc - x(k)
            dxjp = xc - x(k+1)
            dyj  = yc - y(k)
            dyjp = yc - y(k+1) 
            flog = 0.5*alog((dxjp**2 + dyjp**2)/(dxj**2  + dyj**2))
            ftan = atan2( dyjp*dxj - dxjp*dyj , dxjp*dxj + dyjp*dyj )
c..Find ctimtj and stimtj again for grid panels from each panel             
            ctimtj = costhegrid * costhe(k) + sinthegrid * sinthe(k)
            stimtj = sinthegrid * costhe(k) - costhegrid * sinthe(k)
c..Get influence coefficents for each panel in the grid
            Sinfgrid(i,j) = pi2inv * ( ftan * ctimtj + flog * stimtj )
            Ginfgrid(i,j) = pi2inv * ( flog * ctimtj - ftan * stimtj )
c..Sum up local tangential(u_l) and normal(v_l) velocities
            u_l = u_l - Sigma(k)*Ginfgrid(i,j) + Gamma*Sinfgrid(i,j)
            v_l = v_l + Sigma(k)*Sinfgrid(i,j) + Gamma*Ginfgrid(i,j)
            enddo
c..Evaluate global velocity vectors using angle of each grid
            u(i,j)  = u(i,j) + (u_l*costhegrid - v_l*sinthegrid)
            v(i,j)  = v(i,j) + (u_l*sinthegrid + v_l*costhegrid)
         endif 
         cp(i,j) = 1. - (u(i,j)**2 + v(i,j)**2) ! get cp distrubution over domain
      enddo
      enddo
      open(unit=10,file='qout.plt')
      write(10,*) 'VARIABLES= "x","y","u","v","Cp"'
      write(10,*) 'ZONE F=POINT, I=',imax,', J=',jmax-1
      write(10,'(5f12.7)') (( xg(i,j),yg(i,j),u(i,j),v(i,j),cp(i,j),
     > i=1,imax),j=2,jmax)
      close(10)

      
      return
      end
