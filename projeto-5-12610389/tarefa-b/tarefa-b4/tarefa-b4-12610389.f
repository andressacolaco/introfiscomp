        PROGRAM sistema_solar

        implicit real*8(a-h,o-z)
        
        parameter(nb = 9) !Number of bodies
        parameter(pi = dacos(-1d0))
        parameter(aMS = 2d0*10d30)
        parameter (gms = 4d0*pi**2)

        parameter (dt = 0.0005d0)

        dimension aM(1:nb)
        dimension a(1:nb), e(1:nb)
        dimension x(1:nb), y(1:nb)
        dimension vx(1:nb), vy(1:nb)
        dimension xm1(1:nb), ym1(1:nb)
        dimension xp1(1:nb), yp1(1:nb)
        dimension r(nb, nb+1)

        parameter(aM = (/2.4d0*10d23, 4.9d0*10d24, 6d0*10d24, 
     $   6.6d0*10d23, 1.9d0*10d27, 5.7d0*10d26, 8.8d0*10d25, 
     $   1.03d0*10d26, 1.3d0*10d22/))
        parameter(a = (/ 0.39d0, 0.72d0, 1d0, 1.52d0, 5.2d0, 9.24d0, 
     $   19.19d0, 30.06d0, 39.53d0/)) 
        parameter(e = (/ 0.206d0, 0.007d0, 0.017d0, 0.093d0, 0.048d0, 
     $   0.056d0, 0.046d0, 0.010d0, 0.248d0/)) 

        open(10, FILE='saida-b4-mer-12610389.dat')
        open(20, FILE='saida-b4-ven-12610389.dat')
        open(30, FILE='saida-b4-ter-12610389.dat')
        open(40, FILE='saida-b4-mar-12610389.dat')
        open(50, FILE='saida-b4-jup-12610389.dat')
        open(60, FILE='saida-b4-sat-12610389.dat')
        open(70, FILE='saida-b4-ura-12610389.dat')
        open(80, FILE='saida-b4-net-12610389.dat')
        open(90, FILE='saida-b4-plu-12610389.dat')
        
        t_f = 500d0
        t = 0d0
        N = t_f/dt

        do i = 1, nb
              x(i) = a(i)*(1+e(i))
              y(i) = 0d0

              vx(i) = 0d0
              vy(i) = dsqrt(gms*(1-e(i))/a(i)*(1+e(i)))

              write(i*10, *) x(i), y(i)

              xm1(i) = x(i)
              ym1(i) = y(i)

              x(i) = xm1(i) + vx(i)*dt
              y(i) = ym1(i) + vy(i)*dt

              write(i*10, *) x(i), y(i)
        end do

        do j = 1, N
            t = t + dt

            do i = 1, nb
              r(i, 1) = dsqrt(x(i)**2+y(i)**2)

              xp1(i) = 2*x(i) - xm1(i) - (dt**2)*((gms/r(i, 1)**3)*x(i))
              yp1(i) = 2*y(i) - ym1(i) - (dt**2)*((gms/r(i, 1)**3)*y(i))

              do k = 1, nb

                if(k .ne. i) then
                    r(i, k+1) = dsqrt((x(i)-x(k))**2+(y(i)-y(k))**2)
                    xp1(i) = xp1(i) - (dt**2)*
     $               ((((aM(k)/aMs)*gms)/r(i, k+1)**3)*(x(i)-x(k)))
                    yp1(i) = yp1(i) - (dt**2)*
     $               ((((aM(k)/aMs)*gms)/r(i, k+1)**3)*(y(i)-y(k)))
                end if

              end do

              xm1(i) = x(i)
              ym1(i) = y(i)

              x(i) = xp1(i)
              y(i) = yp1(i)

              write(i*10,*) x(i), y(i) 
            end do
       
        end do

        do i = 1, nb
            close(i*10)
        end do

        end program
