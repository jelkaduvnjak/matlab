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

load('fn_EB.mat')

omstart = 2*pi*fn;
om = zeros(1,10);
fn_Wan = zeros(1,10);

for j = 2:10
    om(j) = fzero(@(om) CharEqFreeFreeBeam_Wan(E,I,rho,A,G,mu,l,om),omstart(j));
    fn_Wan(j) = om(j)/2/pi;
end

om = 2*pi*(1:1:2000);
f_om = zeros(1,length(om));

for j = 1: length(om)
    f_om(j) = CharEqFreeFreeBeam_Wan(E,I,rho,A,G,mu,l,om(j));
end
