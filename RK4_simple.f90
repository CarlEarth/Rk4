program testnow
implicit none
!Input the ini_t,final_t,ini_y and ini_ydot.
call evaluate(0.d0,1.d0,-0.4d0,-0.6d0) !###Modify###
end

!Define the forms of coefficients.
!--------------------------------
real (kind=8) function jae1(t,y)
implicit none
 real (kind=8) :: t,y
jae1=-2.d0 !###Modify###
return
end function jae1
!--------------------------------
real (kind=8) function jae2(t,y)
implicit none
real (kind=8) :: t,y
jae2=2.d0 !###Modify###
return
end function jae2
!--------------------------------
real (kind=8) function jae3(t,y)
implicit none
real (kind=8) :: t,y
jae3=-DEXP(2*t)*DSIN(t) !###Modify###
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
!Evaluate the differential equation : y''+J1*y'+J2y+J3=0
!In Runge Kutta 4 ~
implicit none
real (kind=8) :: J1,J2,J3 !(Coeficient of differential equation)
real (kind=8) :: ini_y,ini_ydot !(initial condition y and y')
real (kind=8) :: final_t,ini_t,h,t !(itime )
real (kind=8) :: f1,f2,jae1,jae2,jae3
integer :: i,j
real (kind=8) :: k(4,3)
real (kind=8) :: w1,w2
integer :: Nf=10000000, testit=1
!real(dl) rhode_a(Ncut),w_eff_a(Ncut),a_array(Ncut)

h=(final_t-ini_t)/Nf
t=ini_t
w1=ini_y !(initial y)
w2=ini_ydot !(initial y')
J1=jae1(t,w1)
J2=jae2(t,w1)
J3=jae3(t,w1)

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
if(MOD(i,1000000) .eq. 0) then
!write(*,*)"n=",i
write(10,*)J1,	J2, J3
write(9,*)t,	w1,	w2
endif

!The program will stop when w1 is NAN. 
if (w1 .ne. w1) then
exit
endif

end do

CLOSE(10)
CLOSE(9)

end subroutine evaluate


!end module MBvariables

