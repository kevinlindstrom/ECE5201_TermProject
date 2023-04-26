n_e = 1.55338;
n_o = 1.54425;
eps_e = n_e^2;
eps_o = n_o^2;
c = 2.99792E8;
mu0 = pi*4E-7;
eps0 = 8.8541878E-12;
eta0 = sqrt(mu0/eps0);
eta2_e = sqrt(mu0/(eps_e*eps0));
eta2_o = sqrt(mu0/(eps_o*eps0));

lambda = 590E-9;
w = 2*pi*c/lambda;

beta = w*sqrt(mu0*eps0);
beta_e = w*sqrt(mu0*eps_e*eps0);
beta_o = w*sqrt(mu0*eps_o*eps0);

T12_e = 2*eta2_e/(eta2_e + eta0);
T12_o = 2*eta2_o/(eta2_o + eta0);
T23_e = 2*eta0/(eta2_e + eta0);
T23_o = 2*eta0/(eta2_o + eta0);
G12_e = (eta2_e - eta0)/(eta2_e + eta0);
G12_o = (eta2_o - eta0)/(eta2_o + eta0);
G23_e = (eta0 - eta2_e)/(eta2_e + eta0);
G23_o = (eta0 - eta2_o)/(eta2_o + eta0);

E0 = 1;
m = 1;
L = (2*m + 1)*pi/(beta_e - beta_o);
N = 0;

dz = 1E-10;
z = (L:dz:(L+1E-6));

Ex1 = zeros(1,length(z));
Ex2 = zeros(1,length(z));
Ex3 = zeros(1,length(z));
for n=0:N
    Ec = E0/sqrt(2)*exp(-1j*beta*(z-L))*T12_e*T23_e*G23_e^(2*n)...
        *exp(-1j*(2*n+1)* (beta_o*L + (2*m+1)*pi/2 ));
    Ex1 = Ex1 + Ec;
end

N=1;
for n=0:N
    Ec = E0/sqrt(2)*exp(-1j*beta*(z-L))*T12_e*T23_e*G23_e^(2*n)...
        *exp(-1j*(2*n+1)*(beta_o*L + (2*m+1)*pi/2));
    Ex2 = Ex2 + Ec;
end

N=2;
for n=0:N
    Ec = E0/sqrt(2)*exp(-1j*beta*(z-L))*T12_e*T23_e*G23_e^(2*n)...
        *exp(-1j*(2*n+1)*(beta_o*L + (2*m+1)*pi/2));
    Ex3 = Ex3 + Ec;
end

plot(z,normalize(real(Ex1),'range',[-1 1]),"LineWidth",7)
hold on
plot(z,normalize(real(Ex2),'range',[-1 1]),"LineWidth",7)
plot(z,normalize(real(Ex3),'range',[-1 1]),"LineWidth",7)
ylabel("norm$(\vec{\mathcal{E}}_{3,x}(z;t=0))$", 'Interpreter','latex', 'FontSize', 100)
xlabel("$z$ (m)", 'Interpreter','latex', 'FontSize', 100)
title("Normalized Electric Field $\vec{\mathcal{E}}_{3,x}(z;t=0)$ After QWP", 'Interpreter', ...
    'latex', 'FontSize', 100)
grid on
legend("$n=0$", "$n=1$", '$n=2$','interpreter', 'latex', 'FontSize', 25)
ax = gca;
ax.FontSize = 20;
hold off




