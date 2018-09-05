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

% Characteristics HE100M
rho = 2440;
E = 41938e6;
l = 10.16;
I = 0.246e-02;
A = 196035e-06;

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

