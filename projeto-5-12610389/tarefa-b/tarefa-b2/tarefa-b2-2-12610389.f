        PROGRAM orbitas_tres_corpos

        implicit real*8(a-h,o-z)
        parameter(pi = dacos(-1d0))

        parameter(aMT = 6d0*10d24)
        parameter(aMJ = 1000*1.9d0*10d27) 
        parameter(aMS = 2d0*10d30)

        parameter(aMTS = aMT/aMS)
        parameter(aMJS = aMJ/aMS)

        parameter (gms = 4d0*pi**2)
        parameter(gmt = gms*aMTS)
        parameter(gmj = gms*aMJS)
       
        parameter (dt = 0.0001d0)

        dimension r(1:3)
        dimension x(1:3), y(1:3)
        dimension vx(1:3), vy(1:3)
        dimension xm1(1:3), ym1(1:3)
        dimension xp1(1:3), yp1(1:3)
        dimension d(1:3) !TS, TJ, JS

        parameter(r = (/1d0, 5.2d0, 0d0/))

        t_f = 15d0
        t = 0d0
        N = t_f/dt

        aM = aMT + aMS + aMJ
        raM = aM/aMS

        theta = 0d0

        open(10, FILE='saida-b2-2-terra-12610389.dat')
        open(20, FILE='saida-b2-2-jupiter-12610389.dat')
        open(30, FILE='saida-b2-2-sol-12610389.dat')

        !Terra
        x(1) = r(1)
        y(1) = 0d0

        !Jupiter
        x(2) = r(2)*dcos(theta)
        y(2) = r(2)*dsin(theta)
        
        !Sol
        x(3) = 0d0
        y(3) = 0d0

        d(1) = dsqrt((x(1)-x(3))**2+(y(1)-x(3))**2)
        d(3) = dsqrt((x(2)-x(3))**2+(y(2)-x(3))**2)

        tt = dsqrt((gms+gmt)/d(1))
        tj = dsqrt((gms+gmj)/d(3))

        vx(1) = aMJS*tj*dsin(theta)/raM
        vx(2) = -(1+aMTS)*tj*dsin(theta)/raM
        vx(3) = tj*dsin(theta)/raM

        vy(1) = ((1+aMJS)*tt - aMJS*tj*dcos(theta))/raM
        vy(2) = -(aMTS*tt +(1+aMTS)*tj*dcos(theta))/raM
        vy(3) = (-(aMTS*tt)+aMJS*tj*dcos(theta))/raM

        xcm = (x(1) + x(2) + x(3))/aM
        ycm = (y(1) + y(2) + y(3))/aM

        
        do i = 1, 3
              x(i) = x(i) - xcm
              y(i) = y(i) - ycm

              xm1(i) = x(i)
              ym1(i) = y(i)

              x(i) = xm1(i) + vx(i)*dt
              y(i) = ym1(i) + vy(i)*dt

              write(i*10,*) x(i), y(i) 
        end do

        t_i= 0d0
        t_ant = 0d0
        icount = 0

        do j = 1, N
            t = t + dt

            d(1) = dsqrt((x(1)-x(3))**2+(y(1)-x(3))**2)
            d(2) = dsqrt((x(1)-x(2))**2+(y(1)-y(2))**2)
            
            xp1(1) = 2*x(1) - xm1(1) - (dt**2)*
     $      ((gms/d(1)**3)*(x(1)-x(3))+(gmj/d(2)**3)*(x(1)-x(2)))
            yp1(1) = 2*y(1) - ym1(1) - (dt**2)*
     $      ((gms/d(1)**3)*(y(1)-y(3))+(gmj/d(2)**3)*(y(1)-y(2)))

            d(3) = dsqrt((x(2)-x(3))**2+(y(2)-y(3))**2)

            xp1(2) = 2*x(2) - xm1(2) - (dt**2)*
     $      ((gms/d(3)**3)*(x(2)-x(3))+(gmt/d(2)**3)*(x(2)-x(1)))
            yp1(2) = 2*y(2) - ym1(2) - (dt**2)*
     $      ((gms/d(3)**3)*(y(2)-y(3))+(gmt/d(2)**3)*(y(2)-y(1)))

            xp1(3) = 2*x(3) - xm1(3) - (dt**2)*
     $      ((gmt/d(1)**3)*(x(3)-x(1))+(gmj/d(3)**3)*(x(3)-x(2)))
            yp1(3) = 2*y(3) - ym1(3) - (dt**2)*
     $      ((gmt/d(1)**3)*(y(3)-y(1))+(gmj/d(3)**3)*(y(3)-y(2)))

            xcm = (xp1(1) + xp1(2) + xp1(3))/aM
            ycm = (yp1(1) + yp1(2) + yp1(3))/aM

            do i = 1, 3
                     xm1(i) = x(i)
                     ym1(i) = y(i)

                     x(i) = xp1(i)
                     y(i) = yp1(i)

                     x(i) = x(i) - xcm
                     y(i) = y(i) - ycm

                     write(i*10,*) x(i), y(i) 
           end do
        end do

        close(10)
        close(20)
        close(30)
        end program
