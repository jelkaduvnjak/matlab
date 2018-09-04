# matlab
wan theory
% Characteristics HE100M
rho = 7850;
E = 2.1e11;
nu = 0.3;
G = E/2/(1+nu);
l = 3;
I = 1.143e-05;
A = 5.320e-03;
kappa_ansys = 0.259912;
kappa_SCIA = 1.5785e-03/A;

kappa = kappa_SCIA;
mu = kappa;

omstart = 0;
om = zeros(1,10);
fn_Wan = zeros(1,10);

for j = 2:10
    while om(j)-om(j-1)< 1e-4
        omstart = omstart + 1;
        om(j) = fzero(@(om) CharEqClampedFreeBeam_Wan(E,I,rho,A,G,mu,l,om),omstart);
    end
    fn_Wan(j) = om(j)/2/pi;
end
delta_xpos=0.05;
xpos=(0.0:delta_xpos:3)';
phin_Wan = zeros(length(xpos),10);
phinorm_Wan = zeros(length(xpos),10);
C_norm_Wan = zeros(1,10);

a = sqrt(E*I/(rho*A));
b = sqrt(I/A*(1+E/(G*mu)));
l = 3;
    
for j = 2:10
    g1 = sqrt((-b^2*om(j)^2+sqrt(b^4*om(j)^4+4*a^2*om(j)^2))/(2*a^2));
    g2 = sqrt((b^2*om(j)^2+sqrt(b^4*om(j)^4+4*a^2*om(j)^2))/(2*a^2));
    g3 = g1*(1+rho*b^2*om(j)^2/(G*mu))+E*I/(A*G*mu)*g1^3;
    g4 = g2*(1+rho*b^2*om(j)^2/(G*mu))-E*I/(A*G*mu)*g2^3;
    
    % Solve for c1-4
    Coeff = [g1*g3*cosh(g1*0)    g3*g1*sinh(g1*0)    -g4*g2*cos(g2*0)    -g4*g2*sin(g2*0); ...
         g1*g3*cosh(g1*l)    g3*g1*sinh(g1*l)    -g4*g2*cos(g2*l)    -g4*g2*sin(g2*l); ...
        -(g3-g1)*sinh(g1*0) -(g3-g1)*cosh(g1*0) -(-g4+g2)*sin(g2*0) -(g4-g2)*cos(g2*0); ...
        -(g3-g1)*sinh(g1*l) -(g3-g1)*cosh(g1*l) -(-g4+g2)*sin(g2*l) -(g4-g2)*cos(g2*l)];
    % Since determinant(Coeff) = 0 (otherwise a solution for c1-4 different
    % than 0 cannot be found - > Reduced row echelon form (Gauss-Jordan elimination)
    R = rref(Coeff);
    c1 = -R(1,4);
    c2 = -R(2,4);
    c3 = -R(3,4);
    c4 = 1;
      
    phin_Wan(:,j) = c1*cosh(g1*xpos)+c2*sinh(g1*xpos)+c3*cos(g2*xpos)+c4*sin(g2*xpos);
    C_norm_Wan(j) = sqrt(1/rho/A/trapz(phin_Wan(:,j).^2)/delta_xpos);
    phinorm_Wan(:,j) = C_norm_Wan(j)*phin_Wan(:,j);
end

figure
for mode=2:6
subplot(5,1,mode-1), plot(xpos,-phinorm_Wan(:,mode))
xlabel('x [m]')
ylabel('y')
end
)
function [ y ] = CharEqClampedFreeBeam_Wan(E,I,rho,A,G,mu,l,om)
a = sqrt(E*I/(rho*A));

b = sqrt(I/A*(1+E/(G*mu)));
l=3;


g1 = sqrt((-b^2*om^2+sqrt(b^4*om^4+4*a^2*om^2))/(2*a^2));

g2 = sqrt((b^2*om^2+sqrt(b^4*om^4+4*a^2*om^2))/(2*a^2));

g3 = g1*(1+rho*b^2*om^2/(G*mu))+E*I/(A*G*mu)*g1^3;

g4 = g2*(1+rho*b^2*om^2/(G*mu))-E*I/(A*G*mu)*g2^3;

y =[ g1*g3 0 -g2*g4 0;0 (g3-g1) 0 (g4-g2);cosh(g1*l) sinh(g1*l) cos(g2*l) sin(g2*l);g3*sinh(g1*l) g3*cosh(g1*l) -g4*sin(g2*l) g4*cos(g2*l)]
