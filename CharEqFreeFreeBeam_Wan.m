function [ y ] = CharEq_FreeFreeBeam_Wan(E,I,rho,A,G,mu,l,om)

a = sqrt(E*I/(rho*A));

b = sqrt(I/A*(1+E/(G*mu)));


g1 = sqrt((-b^2*om^2+sqrt(b^4*om^4+4*a^2*om^2))/(2*a^2));

g2 = sqrt((b^2*om^2+sqrt(b^4*om^4+4*a^2*om^2))/(2*a^2));

g3 = g1*(1+rho*b^2*om^2/(G*mu))+E*I/(A*G*mu)*g1^3;

g4 = g2*(1+rho*b^2*om^2/(G*mu))-E*I/(A*G*mu)*g2^3;

y = 2*(g1-g3)*(g2-g4)*(g1^2+rho*om^2/(G*mu))*(g2^2-rho*om^2/(G*mu))*(-1+cos(g2*l)*cosh(g1*l))+((g1-g3)^2*(g2^2-rho*om^2/(G*mu))^2-(g2-g4)^2*(g1^2+rho*om^2/(G*mu))^2)*sin(g2*l)*sinh(g1*l)
