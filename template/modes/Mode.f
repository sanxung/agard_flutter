       program bend_mode
c
c----------------------------------------------------------------------
c
c      AUTHOR: Wendy B. Lessard
c      Version: 2.0
c
c      PURPOSE: The code interpolates the mode shape from a simple FEA beam
c               model onto the ARES rocket stuctured surface grid.
c
c      INPUT: 
c            Unit=1: Surface grid in a MASSOUD (TECPLOT) file format. (modified)
c            Unit=2: File containing the mode shape a simple beam.
c
c      OUTPUT:
c             Files: {project}.bodyN_modeM - New mode shape files for FUN3D.
c
c      VARIABLE DEFINITIONS:
c             PARAMETERS:
c             maxpoints  = Maximum number of surface points.
c             maxelem    = Maximum number of surface elements (triangles)
c             maxmodepts = Maximum number of points to define mode shape.
c
c             i,j,k,l,m,n = do-loop counter, or array indices
c             length      = Reference Length to nondimensionalize surface coordinates.
c             
c             connect     = Array to store surface triangles connectivity data. (data is read in)
c             node        = Array to store corresponding grid node point number for each surface point.
c             coord       = Array to store original surface points. (data is read in)
c             mode_shape  = Array to store displaced surface points.
c             disp_shape_x  = Array to store X-displacement of mode (data is read in.)
c             disp_shape_y  = Array to store Y-displacement of mode (data is read in.)
c             disp_shape_z  = Array to store Z-displacement of mode (data is read in.)
c


c ---------------------------------------------------------------------------
c
       PARAMETER(maxpoints=200000,maxelem=500000,maxmodepts=200
     .           ,nmodes=4)
c
       INTEGER i,j,k,l,m,n,mc
       INTEGER npoints,nelem,nmiddle,nsrf,modepts
       REAL length,displ_fact
       REAL freq,period,dtime,maxdispl


       INTEGER connect(3,maxelem),node(maxpoints)
       REAL coord(3,maxpoints),mode_shape(3,maxpoints,nmodes)
       REAL disp_shape_x(maxmodepts,nmodes)
       REAL disp_shape_y(maxmodepts,nmodes)
       REAL disp_shape_z(maxmodepts,nmodes)
       REAL x(maxmodepts),y(maxmodepts),z(maxmodepts)


       CHARACTER*80 filein_1,filein_2,fileout,project


c---------------------------------------------------------------------------
c---------------------------------------------------------------------------
c      Variables set for case.
c      length & ztransl parameters to adjust data such that the grid and beam coordinates are aligned.
       length   = 3755.36   ! Reference length to scale grid coordinate to match mode coordinates.
c      ztransl  = 334.6     ! Orignal Full Scale CASE : Z - Translation of surface data to match mode coordinates
       ztransl  = 0.0       ! DAC1 CASE: Z - Translation of surface data to match mode coordinates

c      Input files.
       filein_1 = 'wing-445.6_massoud_body1.dat'      !  Input surface file (TECPLOT format)
       filein_2 = '445.6-mode.dat'  !  Mode shape file
c      Output filenames based on project name.
       project = 'wing-445.6'    ! Project Name
       signy   = -1.
       scalecfd = 12.000     
c       nord = 6
c       nlsq = 35
        nord = 10
        nlsq = 45
c---------------------------------------------------------------------------
c---------------------------------------------------------------------------

c      MAIN SECTION

c      READ surface grid from TECPLOT file format
       call read_surface(maxpoints,maxelem,coord,connect,node,
     &                   npoints,nelem,filein_1)
       write(6,*)'READ IN SURFACE POINTS'
       write(6,*)'Surface Points =',npoints,' Surface Triangles =',nelem

c      READ mode shape in from standard curve X,Y file.
       call read_mode(maxmodepts,disp_shape_x,disp_shape_y,
     &                disp_shape_z,x,y,z,modepts,nmodes,filein_2)
       write(6,*)'READING IN MODE DISPLACEMENTS'
       write(6,*)'Number of Structural Node Points =',modepts

c      Compute mode shape for each surface point based on the mode shape.
c
       call bend_surface(maxpoints,maxelem,coord,mode_shape,connect,
     &                   npoints,nelem,maxmodepts,disp_shape_x,
     &                   disp_shape_y,disp_shape_z,x,y,z,signy,scalecfd, 
     &                   modepts,nlsq,nord,nmodes)

       call bend_surface1(maxpoints,maxelem,coord,mode_shape,connect,
     &                    npoints,nelem,maxmodepts,disp_shape_x,
     &                    disp_shape_y,disp_shape_z,x,y,z,signy,
     &                    scalecfd,modepts,nlsq,nord,nmodes)


       call output_surface(maxpoints,maxelem,coord,
     &                     mode_shape,nmodes,connect,npoints,
     &                     nelem,node,project,dtime)


      stop
      end 


c----------------------------------------------------------------------
      subroutine read_surface(maxpoints,maxelem,coord,connect,node,
     &                        npoints,nelem,filein_1)


c     PURPOSE: Read in surface grid from a FUN3D generated MASSOUD 
c              File.
c              NOTE! The MASSOUD header is manually modified such that
c                    it is easier to read the number of points and elemensts.
c                    The strings are remved from the line.
c


      INTEGER n
      INTEGER connect(3,maxelem),ndum
      INTEGER node(maxpoints)
      REAL coord(3,maxpoints)
       
      CHARACTER*80 filein_1,comment,dummy


      OPEN(unit=1,file=filein_1,form='formatted',status='old')

      read (1,1000)dummy
      write(6,1000)dummy
      read (1,1000)dummy
      write(6,1000)dummy


      read(1,10000)npoints,nelem
      write(6,*)'Surface Points =',npoints,'  Triangles=',nelem
!---rtb
!10000 format(23x,i8,6x,i8)

10000 format(23x,i5,4x,i6)
 
!12345678901234567890123
!zone t="mdo body 1", i=50827, j=101359, f=fepoint,  solutiontime= 0.8000000E+03, strandid=0
!10000 format(17x,i6,4x,i6)
!---rtb


      do n = 1, npoints
          read(1,*)coord(1,n),coord(2,n),coord(3,n),node(n)
      end do


      do n = 1, nelem
           read(1,*)connect(1,n),connect(2,n),connect(3,n),ndum
      end do


      CLOSE(1)
c     CLOSE(11)


1000  format(a80)


      return
      end


c----------------------------------------------------------------------
c----------------------------------------------------------------------


      subroutine read_mode(maxmodepts,disp_shape_x,disp_shape_y,
     &                     disp_shape_z,x,y,z,modepts,nmodes,filein_2)


c     PURPOSE: READ in mode shape.  Z-direction is taken along the beam.
c              The mode shapes in the X & Y directions are read in.


      INTEGER maxmodepts,modepts,n
      real max_disp
      REAL disp_shape_x(maxmodepts,nmodes)
      REAL disp_shape_y(maxmodepts,nmodes)
      REAL disp_shape_z(maxmodepts,nmodes)
      REAL x(maxmodepts),y(maxmodepts),z(maxmodepts)


      CHARACTER*80 filein_2,comment


      OPEN(unit=2,file=filein_2,form='formatted',status='old')


c     Read HEADER LINES AS DUMMNY LINES


      read(2,1000)comment
      read(2,*)modepts
      write(6,1000) comment
c     

100   continue

      eps = 1.e-07
      do nm = 1,nmodes
        read(2,1000)comment
        do n = 1, modepts
         read(2,*)x(n),y(n),z(n),disp_shape_x(n,nm),disp_shape_y(n,nm)
     &          ,disp_shape_z(n,nm)
c
c        Hard wire data fix for the 445.6 wing:
c
c        if(abs(y(n)).lt.eps) then
c          disp_shape_x(n,nm) = 0.
c          disp_shape_y(n,nm) = 0.
c          disp_shape_z(n,nm) = 0.
c        end if
        enddo
      enddo

      CLOSE(2)
1000  format(a80)


      return
      end
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
       subroutine bend_surface(maxpoints,maxelem,coord,mode_shape,
     &                         connect,npoints,nelem,maxmodepts,
     &                         disp_shape_x,disp_shape_y,disp_shape_z, 
     &                         x,y,z,signy,scale,modepts,nlsq,nord, 
     &                         nmodes)
c
c     PURPOSE:  Create interpolated surface shape based on mode shape.
c
      INTEGER n
      INTEGER maxpoints,maxelem,npoints,nelem,maxmodepts,modpts
      INTEGER connect(3,maxelem)
      REAL coord(3,maxpoints),mode_shape(3,maxpoints,nmodes)
      REAL disp_shape_x(maxmodepts,nmodes)
      REAL disp_shape_y(maxmodepts,nmodes)
      REAL disp_shape_z(maxmodepts,nmodes)
      REAL x(maxmodepts),y(maxmodepts),z(maxmodepts) 
      dimension iabs(nlsq)
      dimension ul(nord,nord),dist(nlsq),a(nord,nord),c(10),coef(10) 
      dimension x1(10)
c
c
c
      idebug = 0
      do nm = 1,nmodes 
        do n = 1, npoints
          xtp = scale      *coord(1,n)
          ytp = signy*scale*coord(2,n)
          ztp = scale      *coord(3,n)
          iabs = 0
          do l = 1,nlsq
            dist(l) = 1.e+20
            do i = 1,modepts
             do k = 1,l-1
              if(i.eq.iabs(k)) goto 20
             enddo
             test   =
     .            sqrt((x(i)-xtp)*(x(i)-xtp)
     .               + (y(i)-ytp)*(y(i)-ytp))
             if(test.lt.dist(l)) then
               iabs(l) = i
               dist(l) = test 
             end if
20           continue
            enddo
          enddo
          eps = 1.e-2
          do l = 1,nlsq
            if(dist(l).lt.eps) dist(l) = eps
          enddo
          test = 0.
          do l = 1,nlsq
            if(dist(l).gt.test) then
              test = dist(l)
            end if
          enddo
          do l = 1,nlsq
            dist(l) = dist(l)/test
          enddo
          do l = 1,nlsq
            dist(l) = 1.-.95*dist(l)
c           dist(l) = 1.
          enddo
          a = 0.
          c = 0.
          if(nord.eq.6) then
            do l = 1,nlsq 
              x1(1) = 1.
              x1(2) = x(iabs(l))
              x1(3) = x(iabs(l))**2
              x1(4) = y(iabs(l))*x(iabs(l))
              x1(5) = y(iabs(l))
              x1(6) = y(iabs(l))**2
              do jj = 1,nord
                do kk = 1,nord
                  a(kk,jj) = a(kk,jj) + dist(l)*x1(kk)*x1(jj)
                enddo
              enddo
              do kk = 1,nord
                c(kk)   = c(kk) + x1(kk)*dist(l)*
     .                           disp_shape_z(iabs(l),nm) 
              enddo
            enddo
          else if(nord.eq.10) then
            do l = 1,nlsq 
              x1(1) = 1.
              x1(2) = x(iabs(l))
              x1(3) = x(iabs(l))**2
              x1(4) = x(iabs(l))**3
              x1(5) = y(iabs(l))*x(iabs(l))
              x1(6) = y(iabs(l))*x(iabs(l))**2
              x1(7) = x(iabs(l))*y(iabs(l))**2 
              x1(8) = y(iabs(l))
              x1(9) = y(iabs(l))**2
              x1(10)= y(iabs(l))**3
              do jj = 1,nord
                do kk = 1,nord
                  a(kk,jj) = a(kk,jj) + dist(l)*x1(kk)*x1(jj)
                enddo
              enddo
              do kk = 1,nord
                c(kk)   = c(kk) + x1(kk)*dist(l)*
     .                           disp_shape_z(iabs(l),nm) 
              enddo
            enddo
          end if
          call INVDET1(nord,ul,a,dtnrm,detm)
          coef = 0.
          do m = 1,nord
            do l = 1,nord
              coef(l) = coef(l) + ul(l,m)*c(m) 
            enddo
          enddo
          if(idebug.ne.0) then
          if(nord.eq.6) then
            do l = 1,nlsq
              x1(1) = 1.
              x1(2) = x(iabs(l))
              x1(3) = x(iabs(l))**2
              x1(4) = y(iabs(l))*x(iabs(l))
              x1(5) = y(iabs(l))
              x1(6) = y(iabs(l))**2
              dispz = 0.
              do jj = 1,nord
                dispz = dispz + x1(jj)*coef(jj)
              enddo
              write(6,21042) nm,n,l,x1(2),x1(5),dispz,
     .                 disp_shape_z(iabs(l),nm) 
            enddo
          else if(nord.eq.10) then 
            do l = 1,nlsq
              x1(1) = 1.
              x1(2) = x(iabs(l))
              x1(3) = x(iabs(l))**2
              x1(4) = x(iabs(l))**3
              x1(5) = y(iabs(l))*x(iabs(l))
              x1(6) = y(iabs(l))*x(iabs(l))**2
              x1(7) = x(iabs(l))*y(iabs(l))**2 
              x1(8) = y(iabs(l))
              x1(9) = y(iabs(l))**2
              x1(10)= y(iabs(l))**3
              dispz = 0.
              do jj = 1,nord
                dispz = dispz + x1(jj)*coef(jj)
              enddo
              if(xtp.gt.13.420.and.xtp.lt.13.430.
     .          and.ytp.gt.12.622.and.ytp.lt.12.635.
     .          and.nm.eq.3) 
     .          then 
                write(6,21042) n,l,iabs(l),x1(2),x1(8),xtp,ytp
     .                        ,dist(l)
              end if
              if(xtp.gt.13.411.and.xtp.lt.13.420.
     .          and.ytp.gt.12.6148.and.ytp.lt.12.622. 
     .          and.nm.eq.3) 
     .          then 
                write(6,21043) n,l,iabs(l),x1(2),x1(8),xtp,ytp
     .                        ,dist(l)
              end if
            enddo
          end if
          end if
21042     format(' 1st ',3i6,5(1x,e16.8))
21043     format(' 2nd ',3i6,5(1x,e16.8))
          mode_shape(1,n,nm) = 0.
          mode_shape(2,n,nm) = 0.
          mode_shape(3,n,nm) = 0.
          if(nord.eq.10) then
            x1(1) = 1.
            x1(2) = xtp
            x1(3) = xtp**2
            x1(4) = xtp**3
            x1(5) = ytp*xtp
            x1(6) = ytp*xtp**2
            x1(7) = xtp*ytp**2 
            x1(8) = ytp
            x1(9) = ytp**2
            x1(10)= ytp**3
          else if(nord.eq.6) then 
            x1(1) = 1.
            x1(2) = xtp
            x1(3) = xtp**2
            x1(4) = ytp*xtp
            x1(5) = ytp
            x1(6) = ytp**2
          end if
          do jj = 1,nord
            mode_shape(3,n,nm) = mode_shape(3,n,nm) + x1(jj)*coef(jj) 
          enddo
        enddo
      enddo
    


      return
      end
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
       subroutine bend_surface1(maxpoints,maxelem,coord,mode_shape,
     &                          connect,npoints,nelem,maxmodepts,
     &                          disp_shape_x,disp_shape_y,disp_shape_z, 
     &                          x,y,z,signy,scale,modepts,nlsq,nord, 
     &                          nmodes)
c
c     PURPOSE:  Create interpolated surface shape based on mode shape.
c
      INTEGER n
      INTEGER maxpoints,maxelem,npoints,nelem,maxmodepts,modpts
      INTEGER connect(3,maxelem)
      REAL coord(3,maxpoints),mode_shape(3,maxpoints,nmodes)
      REAL disp_shape_x(maxmodepts,nmodes)
      REAL disp_shape_y(maxmodepts,nmodes)
      REAL disp_shape_z(maxmodepts,nmodes)
      REAL x(maxmodepts),y(maxmodepts),z(maxmodepts) 
      dimension iabs(nlsq)
      dimension ul(nord,nord),dist(nlsq),a(nord,nord),c(10),coef(10) 
      dimension x1(10)
c
c
c
      idebug = 1
      do nm = 1,nmodes 
        write(12,31120) modepts,nm 
31120   format(' ZONE I = ',i6,'T="Mode ',i5,'"')
        do n = 1, modepts
          iabs = 0
          do l = 1,nlsq
            dist(l) = 1.e+20
            do i = 1,npoints
             do k = 1,l-1
              if(i.eq.iabs(k)) goto 20
             enddo
             xtp = scale      *coord(1,i)
             ytp = signy*scale*coord(2,i)
             ztp = scale      *coord(3,i)
             test   =
     .            sqrt((x(n)-xtp)*(x(n)-xtp)
     .               + (y(n)-ytp)*(y(n)-ytp)
     .               + (z(n)-ztp)*(z(n)-ztp))
             if(test.lt.dist(l)) then
               iabs(l) = i
               dist(l) = test 
             end if
20           continue
            enddo
          enddo
          do l = 1,nlsq
             xtp = scale      *coord(1,iabs(l))
             ytp = signy*scale*coord(2,iabs(l))
             ztp = scale      *coord(3,iabs(l))
             write(6,21290) x(n),y(n),z(n),xtp,ytp,ztp 
21290        format(6(1x,e16.8))
          enddo
          a = 0.
          c = 0.
          if(nord.eq.6) then
            do l = 1,nlsq 
              xtp = scale      *coord(1,iabs(l))
              ytp = signy*scale*coord(2,iabs(l))
              x1(1) = 1.
              x1(2) = xtp
              x1(3) = xtp**2
              x1(4) = ytp*xtp
              x1(5) = ytp
              x1(6) = ytp**2
              do jj = 1,nord
                do kk = 1,nord
                  a(kk,jj) = a(kk,jj) + x1(kk)*x1(jj)
                enddo
              enddo
              do kk = 1,nord
                c(kk)   = c(kk) + x1(kk)*mode_shape(3,iabs(l),nm) 
              enddo
            enddo
          else if(nord.eq.10) then
            do l = 1,nlsq 
              xtp = scale      *coord(1,iabs(l))
              ytp = signy*scale*coord(2,iabs(l))
              x1(1) = 1.
              x1(2) = xtp
              x1(3) = xtp**2
              x1(4) = xtp**3
              x1(5) = ytp*xtp
              x1(6) = ytp*xtp**2
              x1(7) = xtp*ytp**2 
              x1(8) = ytp
              x1(9) = ytp**2
              x1(10)= ytp**3
              do jj = 1,nord
                do kk = 1,nord
                  a(kk,jj) = a(kk,jj) + x1(kk)*x1(jj)
                enddo
              enddo
              do kk = 1,nord
                c(kk)   = c(kk) + x1(kk)*mode_shape(3,iabs(l),nm) 
              enddo
            enddo
          end if
          call INVDET1(nord,ul,a,dtnrm,detm)
          coef = 0.
          do m = 1,nord
            do l = 1,nord
              coef(l) = coef(l) + ul(l,m)*c(m) 
            enddo
          enddo
          if(debug.ne.0) then
            do l = 1,nlsq
              xtp = scale      *coord(1,iabs(l))
              ytp = signy*scale*coord(2,iabs(l))
              x1(1) = 1.
              x1(2) = xtp
              x1(3) = xtp**2
              x1(4) = xtp**3
              x1(5) = ytp*xtp
              x1(6) = ytp*xtp**2
              x1(7) = xtp*ytp**2 
              x1(8) = ytp
              x1(9) = ytp**2
              x1(10)= ytp**3
              dispz = 0.
              do jj = 1,nord
                dispz = dispz + x1(jj)*coef(jj)
              enddo
              write(6,21042) nm,n,l,xtp,ytp,dispz,
     .                mode_shape(3,iabs(l),nm) 
            enddo
21042     format(3i6,4(1x,e16.8))
          end if
          disp_shape_x1 = 0.
          disp_shape_y1 = 0.
          disp_shape_z1 = 0.
          if(nord.eq.10) then
            x1(1) = 1.
            x1(2) = x(n)
            x1(3) = x(n)**2
            x1(4) = x(n)**3
            x1(5) = y(n)*x(n)
            x1(6) = y(n)*x(n)**2
            x1(7) = x(n)*y(n)**2 
            x1(8) = y(n)
            x1(9) = y(n)**2
            x1(10)= y(n)**3
          else if(nord.eq.6) then 
            x1(1) = 1.
            x1(2) = x(n)
            x1(3) = x(n)**2
            x1(4) = y(n)*x(n)
            x1(5) = y(n)
            x1(6) = y(n)**2
          end if
          do jj = 1,nord
            disp_shape_z1 = disp_shape_z1 + x1(jj)*coef(jj) 
          enddo
          write(12,21020) x(n),y(n),z(n),disp_shape_z1   
     .                 ,disp_shape_z(n,nm) 
21020     format(5(1x,e16.8))
        enddo
      enddo
    


      return
      end

************************************************************************
*
*         SUBROUTINE FOR THE INVERSION OF COMIN USING GAUSS-JORDAN 
*         WITH DOUBLE PIVOTING.
*
************************************************************************
      SUBROUTINE INVDET1(N,UL,COMIN,DTNRM,DETM)
      DIMENSION  UL(N,N),COMIN(N,N)
      DIMENSION C(N,N),J(50),J1(50)
      DO 10 K = 1,n
        DO 12 I = 1,n
          C(I,K) = COMIN(I,K)
12      CONTINUE
10    CONTINUE
      PD = 1.D0
      DO 124 L = 1,N
        DD = 0.D0
        DO 123 K = 1,N
        DD = DD + C(L,K)*C(L,K)
123     END DO
!---rtb
!       DD = DSQRT(DD)
        DD = SQRT(DD)
!---rtb
      PD = PD*DD
124   END DO
      DETM = 1.D0
      DO 125 L = 1,N
      J1(L+20)= L
      J(L+20) = L
125   END DO
      DO 144 L = 1,N
        CC = 0.D0
        M  = L
        M1 = L
        DO 135 K = L,N
          DO 134 K1 = L,N
            IF((ABS(CC)-ABS(C(K1,K))).GE.0.D0) GOTO 134
126         M = K
            M1= K1
            CC = C(K1,K)
134       CONTINUE
135     CONTINUE
127     IF(L.EQ.M) GOTO 170
128     K = J(M+20)
        J(M+20) = J(L+20)
        J(L+20) = K
        DO 137 K = 1,N
          S = C(K,L)
          C(K,L) = C(K,M)
        C(K,M) = S
137     END DO
170     IF(L.EQ.M1) GOTO 140
        K1= J1(M1+20)
        J1(M1+20) = J1(L+20)
        J1(L+20)  = K1
        DO 138 K = 1,N
          S = C(L,K)
          C(L,K) = C(M1,K)
        C(M1,K) = S
138     END DO
140     C(L,L) = 1.D0
        DETM = DETM*CC
        DO 139 M = 1,N
        C(L,M) = C(L,M)/CC
139     END DO
        DO 142 M = 1,N
          IF(L.EQ.M) GOTO 142
129       CC = C(M,L)
          IF(CC.EQ.0.D0) GOTO 142
130       C(M,L) = 0.D0
          DO 141 K = 1,N
          C(M,K) = C(M,K) - CC*C(L,K)
141       END DO
142     CONTINUE
144   CONTINUE
      DO 150 L = 1,N
        IF(J(L+20).EQ.L) GOTO 143
131     M = L
132     M = M + 1
        IF(J(M+20).EQ.L) GOTO 133
136     IF(N.GT.M) GOTO 132
        GOTO 143
133     J(M+20) = J(L+20)
        DO 163 K = 1,N
          CC = C(L,K)
          C(L,K) = C(M,K)
        C(M,K) = CC
163     END DO
        J(L+20) = L
143     IF(J1(L+20).EQ.L) GOTO 150
145     M = L
146     M = M + 1
        IF(J1(M+20).EQ.L) GOTO 148
147     IF(N.GT.M) GOTO 146
        GOTO 150
148     J1(M+20) = J1(L+20)
        DO 149 K = 1,N
          CC = C(K,L)
          C(K,L) = C(K,M)
        C(K,M) = CC
149     END DO
        J1(L+20) = L
150   CONTINUE
      DETM = ABS(DETM)
      DTNRM = DETM/PD
      DO 210 k = 1,N
        DO 212 I = 1,N
          UL(I,K) = C(I,K)
212     CONTINUE
210   CONTINUE
      RETURN
      END


c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine output_surface(maxpoints,maxelem,coord,
     &                          mode_shape,nmodes,connect,
     &                          npoints,nelem,node,project,dtime)


c     PURPOSE:  Output deforming surface data in TECPLOT format
c               compatible with FUN3D moving_grid implementation.
c
      INTEGER n,maxpoints,maxelem,npoints,nelem,l,mc
      INTEGER connect(3,maxelem),node(maxpoints)
      REAL dtime
      REAL coord(3,maxpoints),mode_shape(3,maxpoints,nmodes)
     
      CHARACTER*80 project,fileout
      CHARACTER*3  stat(10)
!---rtb
      character(len=10) ::  string_mode
!---rtb

      DATA stat/ '1','2','3','4','5','6','7','8','9','10'/



c---------------------------------------------------------------------
c     SETTING UP NAMING CONVENTION FOR FUN3D MODE OUTPUT FILE
c     write(6,*)project
      do n = 1, 80
c        write(6,*)n,project(n:n)
         if (project(n:n) .eq. '' .or. project(n:n) .eq. ' ') then
             l = n -1
             goto 1
         end if
      end do

1     continue
      fileout(1:l) = project(1:l)
      fileout(l+1:l+11) = ".body1_mode"

c     MODE 1
      nm = 1
      mc = 1
      nblnk = 0
      do n = 1,3
         if ( stat(mc)(n:n) .eq. ' ' .or. stat(mc)(n:n) .eq. '' ) then
             nblnk = nblnk + n
         end if 
      end do

!---rtb
!     fileout(l+12:l+14) = stat(mc)

!     OPEN(unit=3,file=fileout(1:l+14),form='formatted',
!    &     status='unknown')

!     write(6,21100) fileout(1:l+14),nm 

      write(string_mode,'(i0)') nm
      fileout = trim(project) // '_body1_mode' // trim(string_mode)
     &                        // '.dat'

      open(unit=3,file=fileout,form='formatted',status='unknown')

      write(6,21100) trim(fileout), nm
!---rtb

      write(3,*)'TITLE="',project(1:l),' Mode 1"'
      write(3,*)'VARIABLES= "x" "y" "z" "id" "xmd" "ymd" "zmd"'
      write(3,*)'ZONE  I=',npoints,',  J=',nelem,', F=FEPOINT'


      do n = 1, npoints
         write(3,2000)coord(1,n),coord(2,n),coord(3,n),
     &         node(n),mode_shape(1,n,nm)/12.0,mode_shape(2,n,nm)/12.0,
     &         mode_shape(3,n,nm)/12.0
      end do
      do n = 1, nelem
         write(3,*)connect(1,n),connect(2,n),connect(3,n),connect(3,n)
      end do


      CLOSE(3)

c     Write Second Mode
c     xmd,ymd,zmd --> ymd,xmd,zmd

c     MODE 2
c     SETTING UP NAMING CONVENTION FOR FUN3D MODE OUTPUT FILE
      nm = 2
      mc = 2
      nblnk = 0
      do n = 1,3
         if ( stat(mc)(n:n) .eq. ' ' .or. stat(mc)(n:n) .eq. '' ) then
             nblnk = nblnk + n
         end if
      end do

!---rtb
!     fileout(l+12:l+14) = stat(mc)
 
!     OPEN(unit=3,file=fileout(1:l+14),form='formatted',
!    &     status='unknown')
 
!     write(6,21100) fileout(1:l+14),nm 

      write(string_mode,'(i0)') nm
      fileout = trim(project) // '_body1_mode' // trim(string_mode)
     &                        // '.dat'

      open(unit=3,file=fileout,form='formatted',status='unknown')

      write(6,21100) trim(fileout), nm
!---rtb

      write(3,*)'TITLE="',project(1:l),' Mode 2"'
      write(3,*)'VARIABLES= "x" "y" "z" "id" "xmd" "ymd" "zmd"'
      write(3,*)'ZONE  I=',npoints,',  J=',nelem,', F=FEPOINT'


      do n = 1, npoints
         mode_shape(2,n,nm) = -1. * mode_shape(2,n,nm)
         write(3,2000)coord(1,n),coord(2,n),coord(3,n),
     &         node(n),mode_shape(2,n,nm)/12.0,mode_shape(1,n,nm)/12.0,
     &         mode_shape(3,n,nm)/12.0
      end do
      do n = 1, nelem
         write(3,*)connect(1,n),connect(2,n),connect(3,n),connect(3,n)
      end do


      CLOSE(3)

c     Write Third Mode
c     xmd,ymd,zmd --> ymd,xmd,zmd

c     MODE 3
c     SETTING UP NAMING CONVENTION FOR FUN3D MODE OUTPUT FILE
      nm = 3
      mc = 3
      nblnk = 0
      do n = 1,3
         if ( stat(mc)(n:n) .eq. ' ' .or. stat(mc)(n:n) .eq. '' ) then
             nblnk = nblnk + n
         end if
      end do

!---rtb
!     fileout(l+12:l+14) = stat(mc)

!     OPEN(unit=3,file=fileout(1:l+14),form='formatted',
!    &     status='unknown')

!     write(6,21100) fileout(1:l+14),nm 

      write(string_mode,'(i0)') nm
      fileout = trim(project) // '_body1_mode' // trim(string_mode)
     &                        // '.dat'

      open(unit=3,file=fileout,form='formatted',status='unknown')

      write(6,21100) trim(fileout), nm
!---rtb

      write(3,*)'TITLE="',project(1:l),' Mode 3"'
      write(3,*)'VARIABLES= "x" "y" "z" "id" "xmd" "ymd" "zmd"'
      write(3,*)'ZONE  I=',npoints,',  J=',nelem,', F=FEPOINT'


      do n = 1, npoints
         mode_shape(2,n,nm) = -1. * mode_shape(2,n,nm)
         write(3,2000)coord(1,n),coord(2,n),coord(3,n),
     &         node(n),mode_shape(2,n,nm)/12.0,mode_shape(1,n,nm)/12.0,
     &         mode_shape(3,n,nm)/12.0
      end do
      do n = 1, nelem
         write(3,*)connect(1,n),connect(2,n),connect(3,n),connect(3,n)
      end do


      CLOSE(3)

c     Write Forth Mode
c     xmd,ymd,zmd --> ymd,xmd,zmd

c     MODE 4
c     SETTING UP NAMING CONVENTION FOR FUN3D MODE OUTPUT FILE
      nm = 4
      mc = 4
      nblnk = 0
      do n = 1,3
         if ( stat(mc)(n:n) .eq. ' ' .or. stat(mc)(n:n) .eq. '' ) then
             nblnk = nblnk + n
         end if
      end do

!---rtb
!     fileout(l+12:l+14) = stat(mc)

!     OPEN(unit=3,file=fileout(1:l+14),form='formatted',
!    &     status='unknown')

!     write(6,21100) fileout(1:l+14),nm 

      write(string_mode,'(i0)') nm
      fileout = trim(project) // '_body1_mode' // trim(string_mode)
     &                        // '.dat'

      open(unit=3,file=fileout,form='formatted',status='unknown')

      write(6,21100) trim(fileout), nm
!---rtb

21100 format(' Writing file ',a80,/,' for mode ',i4)

      write(3,*)'TITLE="',project(1:l),' Mode 4"'
      write(3,*)'VARIABLES= "x" "y" "z" "id" "xmd" "ymd" "zmd"'
      write(3,*)'ZONE  I=',npoints,',  J=',nelem,', F=FEPOINT'


      do n = 1, npoints
         mode_shape(2,n,nm) = -1. * mode_shape(2,n,nm)
         write(3,2000)coord(1,n),coord(2,n),coord(3,n),
     &         node(n),mode_shape(2,n,nm)/12.0,mode_shape(1,n,nm)/12.0,
     &         mode_shape(3,n,nm)/12.0
      end do
      do n = 1, nelem
         write(3,*)connect(1,n),connect(2,n),connect(3,n),connect(3,n)
      end do


      CLOSE(3)
1000  format(a3)

!---rtb
!2000  format(3(1x,e13.6),1x,i13,3(1x,e13.6))
2000  format(3(1x,e22.15),1x,i13,3(1x,e22.15))
!---rtb


      return
      end
