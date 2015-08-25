program background; implicit none
! Background evolution for 1/4*lambda*phi^4 + 1/4*mu*sigma^4 + 1/2*g*phi^2*sigma^2
! t = ln(a)
real, parameter :: dt = 0.01
real, parameter :: lambda = 16.0, g = 32.0, mu = 16.0
real, parameter :: unc  = 10**(-6)
real eps
logical endinf
integer i, j, l, k, tend; real a(5)
real x(5)

! Phase space
 do i =1,4
   		a(1) = -10.0*i
		a(2) = -0.5
		a(3) = -30.0*i
		a(4) = -0.5
 		a(5) = sqrt((lambda*a(1)**4+2*g*a(1)**2*a(3)**2+mu*a(3)**4)/(12.0-2*a(2)**2-2*a(4)**2))
		tend = 300000*i
 call tracer(a,tend)
write (*,*) '','','',''
end do

do i =1,4
   		a(1) = 10.0*i
		a(2) = 0.5
		a(3) = 30.0*i
		a(4) = 0.5
 		a(5) = sqrt((lambda*a(1)**4+2*g*a(1)**2*a(3)**2+mu*a(3)**4)/(12.0-2*a(2)**2-2*a(4)**2))
		tend = 300000*i
 call tracer(a,tend)
write (*,*) '','','','',''
end do

! Map of initial conditions 
!do i =0,240,5
!	do j =0,240,5
! Initial conditions for phi
!	a(1) = -120.0 +1.0*i
!	a(2) = 0.0
! Initial conditions for sigma
!	a(3) = -120.0 +1.0*j
!	a(4) = 0.0
!	a(5) = sqrt((lambda*a(1)**4+2*g*a(1)**2*a(3)**2+mu*a(3)**4)/(12.0-2*a(2)**2-2*a(4)**2))
!	endinf = .false.
!	x = a
!	l = 0
!		if (a(1) .ne. 0.0 .or. a(3) .ne. 0.0) then 
!	        do while (endinf == .false.)
!		call gl8(x,dt)
! Epsilon = 1 represents the end of inflation, with "unc" as uncertainty
!		eps = 0.5*(x(2)**2+x(4)**2)
!				if ((1-eps)<= unc) then
!				write (*,'(8g28.20)') a(1) , a(3), l*dt
!				endinf = .true.
!				end if	
!		l = l+1	
!		end do
!		else if (a(1) == 0.0 .or. a(3) == 0.0) then
!		write (*,'(8g28.20)') a(1) , a(3), 0.0
!		end if
!	end do
!	write (*,*) ""
!end do 



contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! implicit Gauss-Legendre methods; symplectic with arbitrary Hamiltonian, A-stable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! evaluate derivatives
subroutine evalf(y, dydx)
        real y(5), dydx(5)
        dydx(1) = y(2)
        dydx(2) = -3*y(2)+0.5*(y(2)*y(2)+y(4)*y(4))*y(2)-lambda*(y(1)*y(1)*y(1))/(y(5)*y(5))-g*y(1)*y(3)*y(3)/(y(5)*y(5))
	dydx(3) = y(4)
	dydx(4)	= -3*y(4)+0.5*(y(2)*y(2)+y(4)*y(4))*y(4)-mu*(y(3)*y(3)*y(3))/(y(5)*y(5))-g*y(3)*y(1)*y(1)/(y(5)*y(5))
	dydx(5) = -0.5*(y(2)**2+y(4)**2)*y(5)
	
end subroutine evalf


! 8th order implicit Gauss-Legendre integrator
subroutine gl8(y, dt)
        integer, parameter :: s = 4, n = 5
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



! Amplitude of perturbations (single mode tracer)
subroutine tracer(a,b)
real a(5)
integer l, b
do l = 0,b
		call gl8(a, dt)	
		write (*,'(8g28.20)') a(1), a(2), a(3), a(4)
       end do
end subroutine tracer 

end






