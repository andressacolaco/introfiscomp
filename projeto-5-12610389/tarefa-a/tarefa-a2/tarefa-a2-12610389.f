        PROGRAM orbitas_elipticas

        implicit real*8(a-h,o-z)
        parameter(pi = dacos(-1d0))
        parameter (gms = 4*pi**2)
        parameter (dt = 0.0001d0)
 
        a = 0.72d0
        e = 0.007
        rm = a*(1+e)

        open(10, FILE='saida-a2-1-12610389.dat')
        open(20, FILE='saida-a2-2-12610389.dat')
 11     format(A5, f6.4)

        t_f = 10d0
        t = 0d0
        N = t_f/dt

        x = rm !Maior distância ao Sol
        y = 0d0

        vx = 0d0
        vy = dsqrt(gms*(1-e)/rm) 

        write(10,*) x, y

        xm1 = x
        ym1 = y

        x = xm1 + vx*dt
        y = ym1 + vy*dt
    
        time_sum = 0d0
        t_ant = 0d0
        
        sm1 = 0d0
        sm2 = 0d0

        icount1 = 0
        icount2 = 0
        
        theta_ant = 0d0

        do i = 1, N
            t = t + dt
            r = (x**2+y**2)**0.5d0
            theta = datan(y/x) 
            
            xp1 = 2*x-xm1-(gms/r**3)*x*dt**2
            yp1 = 2*y-ym1-(gms/r**3)*y*dt**2

            xm1 = x
            ym1 = y

            x = xp1
            y = yp1

            if (mod(i, 100) .eq. 0) then

                if (theta_ant*theta .gt. 0) then
                    area = 0.5d0*(theta - theta_ant)*r**2
                    write(20,*) area
                end if
                
                theta_ant = theta
            end if

            !Medida do período
            if (y*ym1 .lt. 0) then

                if (ym1 .lt. 0) then
                    time_sum = time_sum + (t-t_ant)
                    t_ant = t

                    sm1 = sm1 + abs(x)
                    icount1 = icount1 + 1
                else
                    sm2 = sm2+abs(x)
                    icount2 = icount2+1
                end if
            end if

            write(10, *) x, y
        end do

        T = (time_sum/icount1*1d0)
        sm1 = sm1/icount1*1d0
        sm2 = sm2/icount2*1d0
        am = (sm1+sm2)/2d0
        c = rm-am
        e = c/am
        razaota = T**2/a**3

        write(*,*) 'a =', am
        write(*, 11) 'e = ', e
        write(*,*) 'T =', T 
        write(*,*) 'T**2/R**3 = ', razaota

        close(10)
        close(20)
        end program
