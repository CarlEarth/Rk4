program testnow
implicit none
!----------
!Example: Damped and driven oscillation
!F0cos(wt)=m*x''+cx'+kx
!Consider a mass-spring system.
real (kind=8) :: m,k,c,F0,J1,J2,J3
real (kind=8) :: ini_t ,final_t,ini_y,ini_ydot
m=0.5d0  !m:mass (unit:kg)
k=8.d0  !spring constant
c=0.03d0 !damping term (for friction force)
F0=1.0d0   !Motor can give the periodic force
!---------------------
J1=c/m
J2=k/m
J3=F0/m !J3*coswt

!print 200, "X''+", J1,"X'+",J2,"X-",J3,"coswt = 0"
WRITE (*, FMT= 200) "X''+", J1,"X'+",J2,"X-",J3,"coswt = 0"
200 format(1x, a, 1x,f5.2 , 1x, a, 1x,f5.2 ,1x, a, 1x,f5.2 ,1x,a)
!100 format (7x,'Name:',f12.3,J3, 7x, 'Id:', 1x, 'Weight:')
ini_t = 0.d0
final_t = 50.d0
ini_y = 0.0d0
ini_ydot = 0.0d0
!Input the ini_t,final_t,ini_y and ini_ydot.
!call evaluate(0.d0,1.d0,-0.4d0,-0.6d0) !###Modify###
call evaluate(ini_t,final_t,ini_y,ini_ydot) !###Modify###
end

!Define the forms of coefficients.
!--------------------------------
real (kind=8) function jae1(t,y)
implicit none
 real (kind=8) :: t,y
!jae1=-2.d0 !###Modify###
jae1=0.06d0 !###Modify###
return
end function jae1
!--------------------------------
real (kind=8) function jae2(t,y)
implicit none
real (kind=8) :: t,y
!jae2=2.d0 !###Modify###
jae2=16.d0 !###Modify###
return
end function jae2
!--------------------------------
real (kind=8) function jae3(t,y)
implicit none
real (kind=8) :: t,y,pi,w
pi=3.1415926d0
w=1.5 !resonance w =jae2^0.5
!jae3=-DEXP(2*t)*DSIN(t) !###Modify###
!jae3=-2.d0*DCOS(2.d0*pi*t)
jae3=-2.0d0*DSIN(3.5d0*t)
return
end function jae3
!==================================
real (kind=8) function f1(u1,u2,j1,j2,j3)
implicit none
real (kind=8) :: t,u1,u2,u3,j1,j2,j3
f1=u2
return
end function f1
!==================================
real (kind=8) function f2(u1,u2,j1,j2,j3)
implicit none
real (kind=8) :: t,u1,u2,u3,j1,j2,j3
f2=(-j1*u2-j2*u1-j3)
return
end function f2
!====================================

subroutine evaluate(ini_t,final_t,ini_y,ini_ydot)
!Evaluate the differential equation : y''+J1y'+J2y+J3=0
!In Runge Kutta 4 ~
implicit none
real (kind=8) :: J1,J2,J3 !(Coeficient of differential equation)
real (kind=8) :: ini_y,ini_ydot !(initial condition y and y')
real (kind=8) :: final_t,ini_t,h,t !(itime )
real (kind=8) :: f1,f2,jae1,jae2,jae3
integer :: i,j
real (kind=8) :: k(4,3)
real (kind=8) :: w1,w2
integer :: Nf=5000000, testit=1
!real(dl) rhode_a(Ncut),w_eff_a(Ncut),a_array(Ncut)

h=(final_t-ini_t)/Nf
t=ini_t
w1=ini_y !(initial y)
w2=ini_ydot !(initial y')
J1=jae1(t,w1)
J2=jae2(t,w1)
J3=jae3(t,w1)
OPEN(UNIT=11,FILE='time_position',STATUS='UNKNOWN')
OPEN(UNIT=10,FILE='J1J2J3',STATUS='UNKNOWN')
OPEN(UNIT=9,FILE='df_result',STATUS='UNKNOWN')


!start RK4

do i=1,Nf

!-------------------------
J1=jae1(t,w1)
J2=jae2(t,w1)
J3=jae3(t,w1)
k(1,1)=h*f1(w1,w2,J1,J2,J3)
k(1,2)=h*f2(w1,w2,J1,J2,J3)

!-------------------------

J1=jae1(t+5.d-1*h,w1+5.d-1*k(1,1))
J2=jae2(t+5.d-1*h,w1+5.d-1*k(1,1))
J3=jae3(t+5.d-1*h,w1+5.d-1*k(1,1))
k(2,1)=h*f1(w1+5.d-1*k(1,1),w2+5.d-1*k(1,2),J1,J2,J3)
k(2,2)=h*f2(w1+5.d-1*k(1,1),w2+5.d-1*k(1,2),J1,J2,J3)

!----------------------------

J1=jae1(t+5.d-1*h,w1+5.d-1*k(2,1))
J2=jae2(t+5.d-1*h,w1+5.d-1*k(2,1))
J3=jae3(t+5.d-1*h,w1+5.d-1*k(2,1))
k(3,1)=h*f1(w1+5.d-1*k(2,1),w2+5.d-1*k(2,2),J1,J2,J3)
k(3,2)=h*f2(w1+5.d-1*k(2,1),w2+5.d-1*k(2,2),J1,J2,J3)

!----------------------------
J1=jae1(t+h,w1+k(3,1))
J2=jae2(t+h,w1+k(3,1))
J3=jae3(t+h,w1+k(3,1))
k(4,1)=h*f1(w1+k(3,1),w2+k(3,2),J1,J2,J3)
k(4,2)=h*f2(w1+k(3,1),w2+k(3,2),J1,J2,J3)

!----------------------------------------------
t=ini_t+i*h
w1=w1+(k(1,1)+2.d0*k(2,1)+2.d0*k(3,1)+k(4,1))/6.d0
w2=w2+(k(1,2)+2.d0*k(2,2)+2.d0*k(3,2)+k(4,2))/6.d0

! print the calulation result
if(MOD(i,1000) .eq. 0) then
WRITE (10, FMT= 109) "t=",t,",J1=", J1,",J2=",J2,",J3=",J3
109 format(1x, a, 1x,f5.2 , 1x, a, 1x,f5.2 ,1x, a, 1x,f5.2 ,1x,a, 1x,f5.2)
WRITE (9, FMT= 110) "t=",t,",x=", w1,",v=",w2
110 format(1x, a, 1x,f5.2 , 1x, a, 1x,f5.2 ,1x, a, 1x,f5.2)
WRITE (11, *) t, w1
endif

!The program will stop when w1 is NAN. 
if (w1 .ne. w1) then
exit
endif

end do
CLOSE(11)
CLOSE(10)
CLOSE(9)

end subroutine evaluate


!end module MBvariables

