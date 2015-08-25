program standardpert; implicit none
! Stendard perturbation code for V(phi,sigma)=1/4lambda*phi^4+1/2*g*phi^2*sigma^2+1/4*mu*sigma^4
real, parameter :: dt = 0.0001
real, parameter :: lambda = 16.0, g = 32.0, mu = 16.0
real, parameter :: tpi  = 6.2831853071795864769252867665590057684Q0
real wv
integer i, l, k; real y(10), a(10), H0

   		a(1) = 30.0
		a(2) = 0.0
		a(3) = 10.0
		a(4) = 0.0
 		a(5) = sqrt((lambda*a(1)**4+2*g*a(1)**2*a(3)**2+mu*a(3)**4)/(12.0-2*a(2)**2-2*a(4)**2))
		H0 = a(5)
		a(6) = 0.0
                wv = tpi*H0*100.0
! call tracer(a)
 call spectral(a)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! implicit Gauss-Legendre methods; symplectic with arbitrary Hamiltonian, A-stable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! evaluate derivatives of equations of motion
subroutine evalf(y, dydx)
        real y(10), dydx(10)
        dydx(1) = y(2)
        dydx(2) = -3*y(2)+0.5*(y(2)*y(2)+y(4)*y(4))*y(2)-lambda*(y(1)*y(1)*y(1))/(y(5)*y(5))-g*y(1)*y(3)*y(3)/(y(5)*y(5))
	dydx(3) = y(4)
	dydx(4)	= -3*y(4)+0.5*(y(2)*y(2)+y(4)*y(4))*y(4)-mu*(y(3)*y(3)*y(3))/(y(5)*y(5))-g*y(3)*y(1)*y(1)/(y(5)*y(5))
	dydx(5) = -0.5*(y(2)*y(2)+y(4)*y(4))*y(5)
	dydx(6) = 1.0
	dydx(7) = y(8)
	dydx(8) = (0.5*y(2)*y(2)+0.5*y(4)*y(4)-3)*y(8)-y(7)*(3*lambda*y(1)*y(1)+wv**2*exp(-2*y(6))+g*y(3)*y(3))/(y(5)*y(5))-2*g*y(1)*y(3)*y(9)/(y(5)*y(5))
	dydx(9) = y(10)
	dydx(10) = (0.5*y(2)*y(2)+0.5*y(4)*y(4)-3)*y(10)-y(9)*(3*mu*y(3)*y(3)+wv**2*exp(-2*y(6))+g*y(1)*y(1))/(y(5)*y(5))-2*g*y(1)*y(3)*y(7)/(y(5)*y(5))	   

end subroutine evalf


! 8th order implicit Gauss-Legendre integrator
subroutine gl8(y, dt)
        integer, parameter :: s = 4, n = 10
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

! Spectral code	
subroutine spectral(a)
real a(10), y(10), x(10), u, v, up, vp, sx, cx,  s2x, c2x, sxsq, cxsq, omega1, omega2
integer k,i 
	do i=0,50,2
! Eigenfrequencies in the diagonal basis at initial conditions surface
		omega1 = sqrt(wv*wv/exp(2*a(6))+0.5*(g+3*lambda)*(a(1)*a(1))+0.5*(g+3*mu)*(a(3)*a(3))+0.5*sqrt((g-3*lambda)**2*a(1)**4+(g-3*mu)**2*a(3)**4+a(1)**2*a(3)**2*(14*g**2+6*g*(mu+lambda)-18*lambda*mu)))
		omega2 = sqrt(wv*wv/exp(2*a(6))+0.5*(g+3*lambda)*(a(1)*a(1))+0.5*(g+3*mu)*(a(3)*a(3))-0.5*sqrt((g-3*lambda)**2*a(1)**4+(g-3*mu)**2*a(3)**4+a(1)**2*a(3)**2*(14*g**2+6*g*(mu+lambda)-18*lambda*mu)))
! Diagonal amplitudes
		u = exp(-1.5*a(6))/sqrt(2*omega1)
		v = exp(-1.5*a(6))/sqrt(2*omega2)
		up = exp(-3*a(6))/(2*a(5)*u)
		vp = exp(-3*a(6))/(2*a(5)*v)
! Rotation elements (change of basis)
		s2x = 4*g*a(1)*a(3)/sqrt(16*g**2*a(1)**2*a(3)**2+((3*lambda-g)*a(1)**2+(g-3*mu)*a(3)**2)**2)
		c2x = ((3*lambda-g)*a(1)**2+(g-3*mu)*a(3)**2)/sqrt(16*g**2*a(1)**2*a(3)**2+((3*lambda-g)*a(1)**2+(g-3*mu)*a(3)**2)**2)
		sxsq = 0.5*(1-c2x)
		cxsq = 0.5*(1+c2x)
		sx = sin(0.5*acos(c2x))
		cx = cos(0.5*acos(c2x))
		y = a
		x = a
! Initial conditions
		y(7) = u*cx - v*sx
		y(8) = 0.0
		y(9) = u*sx + v*cx
		y(10) = 0.0
		x(7) = 0.0
		x(8) = up*cx - vp*sx
		x(9) = 0.0
		x(10) = up*sx + vp*cx
		do k =(i/dt),60/dt
			call gl8(y,dt)	
			call gl8(x,dt)
			if((k-(i/dt))*dt==2.0) then
			a(1) = y(1)
			a(2) = y(2)
			a(3) = y(3)
			a(4) = y(4)
			a(5) = y(5)
			end if
		end do
		if(i>10.0) then
		write (*,'(4e28.17e3)') log(wv*exp(a(6))), log((x(7)*x(7)+y(7)*y(7))*(wv**3*exp(3*a(6)))), log((x(9)*x(9)+y(9)*y(9))*(wv**3*exp(3*a(6)))), log(sqrt(x(7)**2+y(7)**2)*sqrt(x(9)**2+y(9)**2)*(wv**3*exp(3*a(6))))
		end if
		a(6) = 1.0*(i+5.0)
	end do
end subroutine spectral

! Amplitude of perturbations (single mode tracer)
subroutine tracer(a)
real a(10), y(10), x(10), u, v, up, vp, s2x, c2x, cx, sx, sxsq, cxsq, omega1, omega2 
integer l
y=a
	do l = 0,10000
	call gl8(y, dt)
	end do
x=y
! Eigenfrequencies in the diagonal basis at initial conditions surface
		omega1 = sqrt(wv*wv/exp(2*a(6))+0.5*(g+3*lambda)*(a(1)*a(1))+0.5*(g+3*mu)*(a(3)*a(3))+0.5*sqrt((g-3*lambda)**2*a(1)**4+(g-3*mu)**2*a(3)**4+a(1)**2*a(3)**2*(14*g**2+6*g*(mu+lambda)-18*lambda*mu)))
		omega2 = sqrt(wv*wv/exp(2*a(6))+0.5*(g+3*lambda)*(a(1)*a(1))+0.5*(g+3*mu)*(a(3)*a(3))-0.5*sqrt((g-3*lambda)**2*a(1)**4+(g-3*mu)**2*a(3)**4+a(1)**2*a(3)**2*(14*g**2+6*g*(mu+lambda)-18*lambda*mu)))
! Diagonal amplitudes
		u = exp(-1.5*a(6))/sqrt(2*omega1)
		v = exp(-1.5*a(6))/sqrt(2*omega2)
		up = exp(-3*a(6))/(2*a(5)*u)
		vp = exp(-3*a(6))/(2*a(5)*v)
! Rotation elements (change of basis)
		s2x = 4*g*a(1)*a(3)/sqrt(16*g**2*a(1)**2*a(3)**2+((3*lambda-g)*a(1)**2+(g-3*mu)*a(3)**2)**2)
		c2x = ((3*lambda-g)*a(1)**2+(g-3*mu)*a(3)**2)/sqrt(16*g**2*a(1)**2*a(3)**2+((3*lambda-g)*a(1)**2+(g-3*mu)*a(3)**2)**2)
		sxsq = 0.5*(1-c2x)
		cxsq = 0.5*(1+c2x)
		y = a
		x = a
! Initial conditions
		y(7) = u*cx - v*sx
		y(8) = 0.0
		y(9) = u*sx + v*cx
		y(10) = 0.0
		x(7) = 0.0
		x(8) = up*cx - vp*sx
		x(9) = 0.0
		x(10) = up*sx + vp*cx	
       do l = 1/dt,125/dt
		call gl8(y, dt)	
		call gl8(x, dt)	
!		write (*,'(4e28.17e3)') l*dt, 0.5*log(x(7)**2+y(7)**2), 0.5*log(x(9)**2+y(9)**2)
! Polar plot
		write (*,'(4e20.10e3)') y(7), x(9), l*dt
       end do
end subroutine tracer 

end


