        PROGRAM orbitas_tres_corpos

        implicit real*8(a-h,o-z)
        parameter(pi = dacos(-1d0))
        parameter(aMT = 6d0*10d24)
        parameter(aMJ = 1.9d0*10d27) 
        parameter(aMS = 2d0*10d30)
        parameter (gms = 4d0*pi**2)
        parameter(gmj = gms*(aMJ/aMS))

        parameter (delta_t = 0.0001d0)

        rt = 1d0
        rj = 5.2d0

        open(10, FILE='saida-b1-terra-12610389.dat')
        open(20, FILE='saida-b1-jupiter-12610389.dat')

        t_f = 30d0
        t = 0d0
        N = t_f/delta_t

        !Terra
        xt = rt
        yt = 0d0

        vxt = 0d0
        vyt = dsqrt(gms/rt)

        !Jupiter
        xj = rj*(1+ej)
        yj = 0d0

        vxj = 0d0
        vyj = dsqrt(gms/rj)

        write(10,*) xt, yt
        write(20,*) xj, yj

        !Terra
        xm1t = xt
        ym1t = yt

        xt = xm1t + vxt*delta_t
        yt = ym1t + vyt*delta_t

        !Jupiter
        xm1j = xj
        ym1j = yj

        xj = xm1j + vxj*delta_t
        yj = ym1j + vyj*delta_t

        write(10,*) xt, yt
        write(20,*) xj, yj

        t_i= 0d0
        t_ant = 0d0
        icount = 0

        do i = 1, N
            t = t + delta_t

            rts = (xt**2+yt**2)**0.5d0
            rtj = ((xt-xj)**2+(yt-yj)**2)**0.5d0
            
            xp1t = 2*xt - xm1t - (delta_t**2)*
     $         ((gms/rts**3)*xt+(gmj/rtj**3)*(xt-xj))
            yp1t = 2*yt - ym1t - (delta_t**2)*
     $         ((gms/rts**3)*yt+(gmj/rtj**3)*(yt-yj))

            rjs = (xj**2+yj**2)**0.5d0

            xp1j = 2*xj - xm1j - (delta_t**2)*
     $         ((gms/rjs**3)*xj+(gmj/rtj**3)*(xj-xt))
            yp1j = 2*yj - ym1j - (delta_t**2)*
     $         ((gms/rjs**3)*yj+(gmj/rtj**3)*(yj-yt))


            xm1t = xt
            ym1t = yt

            xt = xp1t
            yt = yp1t

            xm1j = xj
            ym1j = yj

            xj = xp1j
            yj = yp1j

            write(10, *) xt, yt
            write(20, *) xj, yj
        end do

        close(10)
        close(20)
        end program