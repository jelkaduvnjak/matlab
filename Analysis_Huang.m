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
% Free free
%--------------------------------------------------------------------------


bstart = sqrt(rho*A/(E*I))*l^2*2*pi*fn;
b = zeros(1,10);
fn_Huang = zeros(1,10);

for j = 2:10
    b(j) = fzero(@(b) CharEqFreeFreeBeam_Huang(r,s,b),bstart(j))
    omega = sqrt(E*I/(rho*A))*b(j)/l^2;
    fn_Huang(j) = omega/2/pi;
end

% Calculation of modes

delta = (cosh(b*alfa)-cos(b*beta))/(lam*sinh(b*alfa-zeta*sin(b*beta)));
fn_Huang = D*(cosh(b*alfa*ksi)+lam*delta*sinh(b*alfa*ksi)+cos(b*beta*ksi)/zeta+delta*sin(b*beta*ksi));
delta_xpos=0.05;
xpos=(0.0:delta_xpos:3)';
phin = zeros(length(xpos),10);
phinorm = zeros(length(xpos),10);

C_norm = zeros(1,10);

for j = 1:10
    phin(:,j) = cos((betaxl(j)/l)*xpos)+cosh((betaxl(j)/l)*xpos)+C13(j)*(sin((betaxl(j)/l)*xpos)+sinh((betaxl(j)/l)*xpos));
    C_norm(j) = sqrt(1/rho/A/trapz(phin(:,j).^2)/delta_xpos);
    phinorm(:,j) = C_norm(j)*phin(:,j);
end

figure
for mode=1:5
subplot(5,1,mode), plot(xpos,phin(:,mode))
xlabel('x [m]')
ylabel('y')
end

%--------------------------------------------------------------------------
% Clamped free
%--------------------------------------------------------------------------
bstart = sqrt(rho*A/(E*I))*l^2*2*pi*fn;
b = zeros(1,10);
fn_Huang = zeros(1,10);

for j = 2:10
    b(j) = fzero(@(b) CharEqClampedFreeBeam_Huang(r,s,b),bstart(j))
    omega = sqrt(E*I/(rho*A))*b(j)/l^2;
    fn_Huang(j) = omega/2/pi;
end

% Calculation of modes for Timoshenko clamped free beam %%%
delta = (cosh(b*alfa)-cos(b*beta))/(lam*sinh(b*alfa-zeta*sin(b*beta)));
fn_Huang = D*(cosh(b*alfa*ksi)-lam*zeta*delta*sinh(b*alfa*ksi)-cos(b*beta*ksi)+delta*sin(b*beta*ksi));
delta_xpos=0.05;
xpos=(0.0:delta_xpos:3)';
phin = zeros(length(xpos),10);
phinorm = zeros(length(xpos),10);

C_norm = zeros(1,10);

for j = 1:10
    phin(:,j) = cos((betaxl(j)/l)*xpos)+cosh((betaxl(j)/l)*xpos)+C13(j)*(sin((betaxl(j)/l)*xpos)+sinh((betaxl(j)/l)*xpos));
    C_norm(j) = sqrt(1/rho/A/trapz(phin(:,j).^2)/delta_xpos);
    phinorm(:,j) = C_norm(j)*phin(:,j);
end

figure
for mode=1:5
subplot(5,1,mode), plot(xpos,phin(:,mode))
xlabel('x [m]')
ylabel('y')
end

