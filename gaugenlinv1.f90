program gaugenlin
! Spectral code based in Cholesky decomposition (phase separation included!)
! Same non-linear potential as in standardpert.f90
real, parameter :: dt = 0.02
real, parameter ::  g = 32.0, mu= 16.0, lambda= 16.0
real, parameter :: tpi  = 6.2831853071795864769252867665590057684Q0
real wv
real b(16), H0, a(16)


! Initial conditions from the map in background.f90
       b(1) =  30.0
       b(2) =  0.0
       b(3) =  10.0
       b(4) =  0.0 
       b(5) =  sqrt((0.25*lambda*b(1)*b(1)*b(1)*b(1)+0.5*g*b(1)*b(1)*b(3)*b(3)+0.25*mu*b(3)*b(3)*b(3)*b(3))/(3-0.5*b(2)*b(2)-0.5*b(4)*b(4)))
       H0 = b(5)
       b(6) = 0.0
       wv = tpi*H0*1000.0

 call spectral(b)	
! call tracer(b)

contains
! Equations of motion gauged (considering Cholesky decomposition)
subroutine evalf(y, dydx)
        real y(16), dydx(16)
! Background evolution
! Field 1
        dydx(1) = y(2)
        dydx(2) = (-3+0.5*y(2)*y(2)+0.5*y(4)*y(4))*y(2)-y(1)*(lambda*y(1)*y(1)+g*y(3)*y(3))/(y(5)*y(5))
! Field 2
	dydx(3) = y(4)
	dydx(4) = (-3+0.5*y(2)*y(2)+0.5*y(4)*y(4))*y(4)-y(3)*(mu*y(3)*y(3)+g*y(1)*y(1))/(y(5)*y(5))
! H and log a
	dydx(5) = -0.5*(y(2)*y(2)+y(4)*y(4))*y(5)
	dydx(6) = 1.0
! Phase derivatives
! Phase derivative 1
	dydx(7) = (-3-2*y(10)/y(9))*y(7)+2*(-y(12)*y(13)+y(11)*y(14))/(y(9)*y(11))*y(8) 
! Phase derivative 2
	dydx(8) = (-3-2*y(12)/y(11))*y(8)
! Phase 1
! Evolution for the 1-1 field component
! Mass 1 +3*lambda*y(1)**2+g*y(3)**2
	dydx(9) = y(10)
	dydx(10) = (-3+0.5*y(2)*y(2)+0.5*y(4)*y(4))*y(10)+y(9)*y(7)*y(7)/(y(5)*y(5))-(wv*wv*exp(-2*y(6))+3*lambda*y(1)*y(1)+g*y(3)*y(3))*y(9)/(y(5)*y(5))+2*g*y(3)*y(1)*y(13)*y(9)/(y(11)*y(5)*y(5))
! Evolution for the 2-2 field component
! Mass 2 +3*mu*y(3)**2+g*y(1)**2
	dydx(11) = y(12)
	dydx(12) = (-3+0.5*y(2)*y(2)+0.5*y(4)*y(4))*y(12)+y(11)*y(8)*y(8)/(y(5)*y(5))-(wv*wv*exp(-2*y(6))+3*mu*y(3)*y(3)+g*y(1)*y(1))*y(11)/(y(5)*y(5))-2*g*y(1)*y(3)*y(13)/(y(5)*y(5))
! Evolution of the 1-2 field component
	dydx(13) = y(14)
	dydx(14) = (-3+0.5*y(2)*y(2)+0.5*y(4)*y(4))*y(14)+y(13)*y(8)*y(8)/(y(5)*y(5))-(wv*wv*exp(-2*y(6))+3*lambda*y(1)*y(1)+g*y(3)*y(3))*y(13)/(y(5)*y(5))-2*g*y(1)*y(3)*y(11)/(y(5)*y(5))-2*g*y(1)*y(3)*y(9)*y(9)/(y(11)*y(5)*y(5))
! Phase values
	dydx(15) = y(7)/y(5)
	dydx(16) = y(8)/y(5)
end subroutine evalf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! implicit Gauss-Legendre methods; symplectic with arbitrary Hamiltonian, A-stable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 8th order implicit Gauss-Legendre integrator
subroutine gl8(y, dt)
        integer, parameter :: s = 4, n = 16
        real y(n), g(n,s), dt; integer i, k
        
        ! Butcher tableau for 8th order Gauss-Legendre method
        real, parameter :: a(s,s) = reshape((/ &
                 0.869637112843634643432659873054998518Q-1, -0.266041800849987933133851304769531093Q-1, &
                 0.126274626894047245150568805746180936Q-1, -0.355514968579568315691098184956958860Q-2, &
                 0.188118117499868071650685545087171160Q0,   0.163036288715636535656734012694500148Q0,  &
                -0.278804286024708952241511064189974107Q-1,  0.673550059453815551539866908570375889Q-2, &
                 0.167191921974188773171133305525295945Q0,   0.353953006033743966537619131807997707Q0,  &
                 0.163036288715636535656734012694500148Q0,  -0.141906949311411429641535704761714564Q-1, &
                 0.177482572254522611843442956460569292Q0,   0.313445114741868346798411144814382203Q0,  &
                 0.352676757516271864626853155865953406Q0,   0.869637112843634643432659873054998518Q-1 /), (/s,s/))
        real, parameter ::   b(s) = (/ &
                 0.173927422568726928686531974610999704Q0,   0.326072577431273071313468025389000296Q0,  &
                 0.326072577431273071313468025389000296Q0,   0.173927422568726928686531974610999704Q0  /)
        
        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i))
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl8
! Spectral routine
subroutine spectral(a)
real:: y(16), a(16), uu, vv, cx, cxsq, sx, sxsq, s2x, c2x, omega1, omega2 
integer k, i


! Mode launching every 5 e-folds
	do i=0,60,5
! frequencies (diagonal basis)
		omega1 = sqrt(wv*wv/exp(2*a(6))+0.5*(g+3*lambda)*(a(1)*a(1))+0.5*(g+3*mu)*(a(3)*a(3))+0.5*sqrt((g-3*lambda)**2*a(1)*a(1)*a(1)*a(1)+(g-3*mu)**2*a(3)*a(3)*a(3)*a(3)+a(1)*a(1)*a(3)*a(3)*(14*g**2+6*g*(mu+lambda)-18*lambda*mu)))
		omega2 = sqrt(wv*wv/exp(2*a(6))+0.5*(g+3*lambda)*(a(1)*a(1))+0.5*(g+3*mu)*(a(3)*a(3))-0.5*sqrt((g-3*lambda)**2*a(1)*a(1)*a(1)*a(1)+(g-3*mu)**2*a(3)*a(3)*a(3)*a(3)+a(1)*a(1)*a(3)*a(3)*(14*g**2+6*g*(mu+lambda)-18*lambda*mu)))
! Diagonal correlators 
		uu = exp(-3*a(6))/(2*omega1)
		vv = exp(-3*a(6))/(2*omega2)
! Rotation elements (back to the non-diagonal basis)
		s2x = 4*g*a(1)*a(3)/sqrt(16*g**2*a(1)*a(1)*a(3)*a(3)+((3*lambda-g)*a(1)*a(1)+(g-3*mu)*a(3)*a(3))**2)
		c2x = ((3*lambda-g)*a(1)*a(1)+(g-3*mu)*a(3)*a(3))/sqrt(16*g**2*a(1)*a(1)*a(3)*a(3)+((3*lambda-g)*a(1)*a(1)+(g-3*mu)*a(3)*a(3))**2)
		sxsq = 0.5*(1-c2x)
		cxsq = 0.5*(1+c2x)
		cx = cos(0.5*acos(c2x))
		sx = sin(0.5*acos(c2x))
! Initial conditions (written in Cholesky representation) 
		a(7) = omega1*cx - omega2*sx 
		a(8) = omega1*sx + omega2*cx
		a(10) = 0.0
		a(11) = sqrt(uu*sxsq+vv*cxsq+sqrt(uu*vv)*s2x)
		a(12) = 0.0
		a(13) = (0.5*(uu-vv)*s2x+sqrt(uu*vv)*c2x)/a(11)
		a(9) = sqrt(abs(uu*cxsq+vv*sxsq-sqrt(uu*vv)*s2x-a(13)*a(13)))
		a(14) = 0.0
		a(15) = 0.0
		a(16) = 0.0
      		y = a
		do k =i/dt,80/dt
			call gl8(y,dt)
			if((k-i/dt)*dt==5.0) then
			a(1) = y(1)
			a(2) = y(2)
			a(3) = y(3)
			a(4) = y(4)
			a(5) = y(5)
			end if
		end do
		if(i>10.0) then
		write (*,'(4e32.20e3)') (wv*exp(a(6))), ((y(9)*y(9)+y(13)*y(13))*(wv*wv*wv*exp(3*a(6)))), ((y(11)*y(11))*(wv*wv*wv*exp(3*a(6)))), (abs(y(11)*y(13))*(wv*wv*wv*exp(3*a(6))))
		! Plot in logscale
		end if
		a(6) = i+5.0
	end do
end subroutine spectral 

! Amplitude of perturbations (single mode tracer)
subroutine tracer(y)
real:: y(16), uu, vv, cxsq, sxsq, cx, sx, s2x, c2x, omega1, omega2 
real,parameter :: dt1 = 0.0001
integer l
! Minkowski vacuum set as initial condition for both fields
! Only change the time scale from dt to dt1 when full phase tracing (polar plot) is required, otherwise dt is enough. dt1 is considered by default.
       do l = 0,125/dt1
		call gl8(y, dt1)
		if (l*dt1 == 1.0) then
! frequencies (diagonal basis)
		omega1 = sqrt(wv*wv/exp(2*y(6))+0.5*(g+3*lambda)*(y(1)*y(1))+0.5*(g+3*mu)*(y(3)*y(3))+0.5*sqrt((g-3*lambda)**2*y(1)**4+(g-3*mu)**2*y(3)**4+y(1)**2*y(3)**2*(14*g**2+6*g*(mu+lambda)-18*lambda*mu)))
		omega2 = sqrt(wv*wv/exp(2*y(6))+0.5*(g+3*lambda)*(y(1)*y(1))+0.5*(g+3*mu)*(y(3)*y(3))-0.5*sqrt((g-3*lambda)**2*y(1)**4+(g-3*mu)**2*y(3)**4+y(1)**2*y(3)**2*(14*g**2+6*g*(mu+lambda)-18*lambda*mu)))
! Diagonal correlators 
		uu = exp(-3*y(6))/(2*omega1)
		vv = exp(-3*y(6))/(2*omega2)
! Rotation elements
		s2x = 4*g*y(1)*y(3)/sqrt(16*g**2*y(1)**2*y(3)**2+((3*lambda-g)*y(1)**2+(g-3*mu)*y(3)**2)**2)
		c2x = ((3*lambda-g)*y(1)**2+(g-3*mu)*y(3)**2)/sqrt(16*g**2*y(1)**2*y(3)**2+((3*lambda-g)*y(1)**2+(g-3*mu)*y(3)**2)**2)
		sxsq = 0.5*(1-c2x)
		cxsq = 0.5*(1+c2x)
		cx = cos(0.5*acos(c2x))
		sx = sin(0.5*acos(c2x))
! Initial conditions 
		y(7) = omega1*cx - omega2*sx 
		y(8) = omega2*sx + omega2*cx
		y(10) = 0.0
		y(11) = sqrt(uu*sxsq+vv*cxsq+sqrt(uu*vv)*s2x)
		y(12) = 0.0
		y(13) = (0.5*(uu-vv)*s2x+sqrt(uu*vv)*c2x)/y(11)
		y(9) = sqrt(abs(uu*cxsq+vv*sxsq-sqrt(uu*vv)*s2x-y(13)*y(13)))
		y(14) = 0.0
		y(15) = 0.0
		y(16) = 0.0
		end if
		if (l*dt1>=1.0) then
! Amplitude plots (delete log and plot in logscale, if you require) 
! Does not need dt1 as a time step
!		write (*,'(4e28.17e3)') l*dt, log(y(9)), log(y(11)), log(y(13))
! Full mode evolution (polar plot)
		write (*,'(4e28.17e3)')  y(9)*cos(y(15))+y(13)*cos(y(16)), y(11)*sin(y(16)), l*dt1
		end if
 	end do
end subroutine tracer
end

