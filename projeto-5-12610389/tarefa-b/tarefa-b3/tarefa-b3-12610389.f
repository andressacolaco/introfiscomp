        PROGRAM orbitas_asteroides

        implicit real*8(a-h,o-z)
        
        parameter(pi = dacos(-1d0))
        parameter(aMJ = 1.9d0*10d27) 
        parameter(aMS = 2d0*10d30)
        parameter (gms = 4d0*pi**2)
        parameter(gmj = gms*(aMJ/aMS))

        parameter (dt = 0.0001d0)

        dimension r(1:4)
        dimension x(1:4), y(1:4)
        dimension vx(1:4), vy(1:4)
        dimension xm1(1:4), ym1(1:4)
        dimension xp1(1:4), yp1(1:4)
        dimension rs(1:4), rj(1:3)

        parameter(r = (/3d0, 3.276d0, 3.7d0, 5.2d0/))
        parameter(vx = (/0d0, 0d0, 0d0, 0d0/))
        parameter(vy = (/3.628d0, 3.471d0, 3.267d0, 2.755d0/))

        open(10, FILE='saida-b3-ast-1-12610389.dat')
        open(20, FILE='saida-b3-ast-2-12610389.dat')
        open(30, FILE='saida-b3-ast-3-12610389.dat')
        open(40, FILE='saida-b3-jupiter-12610389.dat')

        t_f = 200d0
        t = 0d0
        N = t_f/dt

        do i = 1, 4
              x(i) = r(i)
              y(i) = 0d0

              write(i*10, *) x(i), y(i)

              xm1(i) = x(i)
              ym1(i) = y(i)

              x(i) = xm1(i) + vx(i)*dt
              y(i) = ym1(i) + vy(i)*dt

              write(i*10, *) x(i), y(i)
        end do

        do j = 1, N
            t = t + dt

            do i = 1, 3
              rs(i) = (x(i)**2+y(i)**2)**0.5d0
              rj(i) = ((x(i)-x(4))**2+(y(i)-y(4))**2)**0.5d0

              xp1(i) = 2*x(i) - xm1(i) - (dt**2)*
     $        ((gms/rs(i)**3)*x(i)+(gmj/rj(i)**3)*(x(i)-x(4)))
              yp1(i) = 2*y(i) - ym1(i) - (dt**2)*
     $        ((gms/rs(i)**3)*y(i)+(gmj/rj(i)**3)*(y(i)-y(4)))

            end do

            !Jupiter
            rs(4) = (x(4)**2+y(4)**2)**0.5d0

            xp1(4) = 2d0*x(4) - xm1(4) - (dt**2)*((gms/rs(4)**3)*x(4))
            yp1(4) = 2d0*y(4) - ym1(4) - (dt**2)*((gms/rs(4)**3)*y(4))


            do i = 1, 4 
              xm1(i) = x(i)
              ym1(i) = y(i)

              x(i) = xp1(i)
              y(i) = yp1(i)

              write(i*10,*) x(i), y(i) 
            end do
       
        end do

        close(10)
        close(20)
        close(30)
        close(40)
        end program
