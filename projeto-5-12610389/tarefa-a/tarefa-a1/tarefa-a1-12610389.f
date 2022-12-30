        PROGRAM orbitas

        implicit real*8(a-h,o-z)
        parameter(pi = dacos(-1d0))
        parameter (gms = 4*pi**2)
        parameter (dt = 0.0001d0)

        r = 0.72d0

 11     format(A5, f6.4)
        open(10, FILE='saida-a1-12610389.dat')

        t_f = 20d0
        t = 0d0
        N = t_f/dt

        x = r !Raio do planeta
        y = 0d0

        vx = 0d0
        vy = dsqrt(gms/r)

        write(10,*) x, y

        xm1 = x
        ym1 = y

        x = xm1 + vx*dt
        y = ym1 + vy*dt
    
        time_sum = 0d0
        t_ant = 0d0
        icount = 0

        sm1 = 0d0
        sm2 = 0d0

        icount1 = 0
        icount2 = 0

        do i = 1, N
            t = t + dt
            r = (x**2+y**2)**0.5d0
            
            xp1 = 2*x-xm1-(gms/r**3)*x*dt**2
            yp1 = 2*y-ym1-(gms/r**3)*y*dt**2

            xm1 = x
            ym1 = y

            x = xp1
            y = yp1

            !Medida do período
            if (y*ym1 .lt. 0) then
                time_sum = time_sum + (t-t_ant)
                t_ant = t
                icount = icount+1

                if (ym1 .lt. 0) then
                    sm1 = sm1 + abs(x)
                    icount1 = icount1 + 1
                else
                    sm2 = sm2 + abs(x)
                    icount2 = icount2 + 1
                end if
            end if

            write(10, *) x, y
        end do

        T = 2*(time_sum/icount*1d0)
        sm1 = sm1/icount1*1d0
        sm2 = sm2/icount2*1d0
        a = (sm1+sm2)/2d0
        c = r-a
        e = c/a
        razaotr = T**2/r**3

        write(*,*) 'a =', a
        write(*,*) 'v0 =', vy
        write(*, 11) 'e = ', e
        write(*,*) 'T =', T 
        write(*,*) 'T**2/R**3 = ', razaotr

        close(10)
        end program
