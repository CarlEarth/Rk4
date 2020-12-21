# Rk4
This is a simple Runge Kutta 4 example written in fortran 90.
You can use it to solve second order differential equation.
y''+J1y'+J2y+J3=0
The form of the equation can be set in the fucntion jae(1-3).

Step: (in Linux)
1. gfortran RK4_simple.f90
2. ./a.out

Then it will output the <J1(t),J2(t),J3(t)> and < y'(t),y(t)>