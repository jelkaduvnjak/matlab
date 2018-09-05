function [ y ] = CharEqSuportedSuportedBeam_Huang(r,s,b)
dummy1 = r^2 + s^2;
dummy2 = (r^2 - s^2)^2+4/b^2;
b = zeros(1,10);
alfa = 1/sqrt(2)*sqrt(-dummy1+sqrt(dummy2));

beta = 1/sqrt(2)*sqrt(dummy1+sqrt(dummy2));
lambda = alfa/beta;
zeta=(beta^2-s^2)/(alfa^2+s^2);


y = sin(b*beta);


