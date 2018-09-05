function [ y ] = CharEqFreeFreeBeam_Huang(r,s,b)

dummy1 = r^2 + s^2;
dummy2 = (r^2 - s^2)^2+4/b^2;

alfa = 1/sqrt(2)*sqrt(-dummy1+sqrt(dummy2));

beta = 1/sqrt(2)*sqrt(dummy1+sqrt(dummy2));
lam = alfa/beta;
zeta=(beta^2-s^2)/(alfa^2+s^2);
%ksi=x/l;
%D = 1/sqrt((1-beta^2)^2+(2*ksi*beta)^2);

y = 2 - 2*cosh(b*alfa)*cos(b*beta)+b/sqrt(1-b^2*r^2*s^2)*(b^2*r^2*(r^2-s^2)^2+(3*r^2-s^2))*sinh(b*alfa)*sin(b*beta);
%y = 2 - 2*cos(b*alfa)*cos(b*beta)+b/sqrt(b^2*r^2*s^2-1)*(b^2*r^2*(r^2-s^2)^2+(3*r^2-s^2))*sin(b*alfa)*sin(b*beta);
%ksi=x/l;
%delta = (cosh(b*alfa)-cos(b*beta))/(lam*sinh(b*alfa-zeta*sin(b*beta)));

%%%%   Mode shapes for Timoshenko free free beam  %%%%
% fn_Huang = D*(cosh(b*alfa*ksi)+lam*delta*sinh(b*alfa*ksi)+cos(b*beta*ksi)/zeta+delta*sin(b*beta*ksi));
% D= 1/sqrt((1-beta^2)^2+(2*ksi*beta)^2);

