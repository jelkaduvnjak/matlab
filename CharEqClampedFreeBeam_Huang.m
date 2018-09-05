function [ y ] = CharEqClampedFreeBeam_Huang(r,s,b)

dummy1 = r^2+s^2;
dummy2 = (r^2 - s^2)^2+4/b^2;

alfa = 1/sqrt(2)*sqrt(-dummy1+sqrt(dummy2));

beta = 1/sqrt(2)*sqrt(dummy1+sqrt(dummy2));
lambda = alfa/beta;
zeta=(beta^2-s^2)/(alfa^2+s^2);
delta = 1/lambda*(sinh(b*alfa)-sin(b*beta))/(zeta*cosh(b*alfa)+cos(b*beta));
delta_xpos = 0.05;
xpos = (0.0:delta_xpos:3)';
l = 3;
ksi=xpos/l;
D=1;
%D = 1/sqrt((1-beta^2)^2+(2*ksi*beta)^2);

y = 2 + (b^2*(r^2-s^2)^2+2)*cosh(b*alfa)*cos(b*beta)-b*dummy1/sqrt(1-b^2*r^2*s^2)*sinh(b*alfa)*sin(b*beta);
