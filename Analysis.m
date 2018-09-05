% Characteristics HE100M
rho = 7850;
E = 2.1e11;
l = 3;
I = 1.143e-05;
A = 5.320e-03;


betaxl = zeros(1,10);

x0 = 0;
betaxl(1) = 0;

for j = 2:10
    betaxl(j) = betaxl(j-1);
    while betaxl(j) == betaxl(j-1)
        x0 = x0 + 2;
        betaxl(j) = fzero('CharEqFreeFreeBeam',[x0]);
    end
end


C13 = zeros(1,10);
C13(1) = 0;

for j = 2:10
    C13(j) = (cos(betaxl(j))-cosh(betaxl(j)))/(sinh(betaxl(j))-sin(betaxl(j)));
end



fn = (1/2/pi)*sqrt(E*I/rho/A)*(betaxl/l).^2;

delta_xpos = 0.05;
xpos = (0.0:delta_xpos:3)';
phin = zeros(length(xpos),10);
phinorm = zeros(length(xpos),10);

C_norm = zeros(1,10);


for j = 1:10
    phin(:,j) = cos((betaxl(j)/l)*xpos)+cosh((betaxl(j)/l)*xpos)+C13(j)*(sin((betaxl(j)/l)*xpos)+sinh((betaxl(j)/l)*xpos));
    C_norm(j) = sqrt(1/rho/A/trapz(phin(:,j).^2)/delta_xpos);
    phinorm(:,j) = C_norm(j)*phin(:,j);
end

figure
for mode=2:6
subplot(5,1,mode-1), plot(xpos,phinorm(:,mode))
xlabel('x [m]')
ylabel('y')
end

print -dpdf -r600 Modes_EB



