# Rk4
This is a simple Runge Kutta 4 example written in fortran 90.
You can use it to solve second order differential equation.
y''+J1y'+J2y+J3=0
The form of the equation can be set in the fucntion jae(1-3).

Step: (in Linux)
1. gfortran RK4_simple.f90
2. ./a.out
3. python3 plot2.py

or type the command:

sh script.sh

Then it will output the <J1(t),J2(t),J3(t)> , < y'(t),y(t)>, < t, y(t)> 
and the x-t plot animation.

In the code, the physical example is written in the code to help the user 
to understand how to use the code. (Example is Damped and driven oscillation.)

Modify the J1(t),J2(t),J3(t) function to solve your second order differential
problem.

plot1.py: Output the rest x-t plot.
plot3.py: Output the oscillator animation.
