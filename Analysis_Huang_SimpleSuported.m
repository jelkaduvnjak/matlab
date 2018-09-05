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

delta_xpos = 0.05;
xpos = (0.0:delta_xpos:3)';
ksi=xpos/l;

kappa = kappa_SCIA

r = sqrt(I/(A*l^2));
s = sqrt(E*I/(kappa*A*G*l^2));

load('fn_EB.mat')

%--------------------------------------------------------------------------
% Simple supported beam
%--------------------------------------------------------------------------


bstart = sqrt(rho*A/(E*I))*l^2*2*pi*fn;
b = zeros(1,10);
fn_Huang = zeros(1,10);

for j = 2:10
    b(j) = fzero(@(b) CharEqSuportedSuportedBeam_Huang(r,s,b),bstart(j))
    omega = sqrt(E*I/(rho*A))*b(j)/l^2;
    fn_Huang(j) = omega/2/pi;
end

% Calculation of modes
delta_xpos=0.05;
xpos=(0.0:delta_xpos:3)';
phin_Huang = zeros(length(xpos),10);
phinorm_Huang = zeros(length(xpos),10);
C_norm_Huang = zeros(1,10);

for j = 2:10
    dummy1 = r^2+s^2;
    dummy2 = (r^2 - s^2)^2+4/b(j)^2;
    alfa = 1/sqrt(2)*sqrt(-dummy1+sqrt(dummy2));
    beta = 1/sqrt(2)*sqrt(dummy1+sqrt(dummy2));
    lambda = alfa/beta;
    zeta = (alfa^2+r^2)/(alfa^2+s^2);
    delta = (1/lambda*sinh(b(j)*alfa)-sin(b(j)*beta))/(zeta*cosh(b(j)*alfa)+cos(b(j)*beta));
    ksi = xpos/l;
    %phin_Huang(:,j) = cosh(b(j)*alfa*ksi)+lambda*delta*sinh(b(j)*alfa*ksi)+cos(b(j)*beta*ksi)/zeta+delta*sin(b(j)*beta*ksi);
    phin_Huang(:,j) = sinh(b(j)*beta*ksi);
    C_norm_Huang(j) = sqrt(1/rho/A/trapz(phin_Huang(:,j).^2)/delta_xpos);
    phinorm_Huang(:,j) = C_norm_Huang(j)*phin_Huang(:,j);
   
end

figure
for mode=2:6
subplot(5,1,mode-1), plot(xpos,phinorm_Huang(:,mode))
xlabel('x [m]')
ylabel('y')
end
