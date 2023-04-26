format long

% Look into RCWA methods, could potentially write this as a matrix or in
% series form

%% parameters
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



%% calculations
m = 3;
L = (2*m + 1)*pi/(2*(beta_e - beta_o));

dz = 0.00001;
z = (-0.1:dz:(L/(1E-4)+.1))*1E-4;

E0 = 1;
E1 = [1; 1; 0]*E0/sqrt(2)*exp(-1j*beta*z);
E2 = [T12_e*exp(-1j*beta_e*z); T12_o*exp(-1j*beta_o*z); ...  
    zeros(1,length(z))]*E0/sqrt(2);
% E3 = [exp(-1j*(2*m+1)*pi/2)*T12_e*T23_e;
%     T12_o*T23_o;
%     0]*...
%     E0/sqrt(2)*exp(-1j*beta_o*L)*exp(-1j*(beta*(z-L)));
E3 = [exp(-1j*(2*m+1)*pi/2)*T12_e*T23_e; T12_o*T23_o; 0]*...
    E0/sqrt(2)*exp(-1j*(beta*(z-L)))*exp(-1j*beta_o*L);

E = zeros(3, length(z));
E(:,1:find(z==0)) = E1(:,1:find(z==0));
E(:,(find(z==0)+1):find(z<L, 1, 'last')) = E2(:,(find(z==0)+1):find(z<L, 1, 'last'));
E(:,find(z>L, 1):end) = E3(:,find(z>L, 1):end);
E(3,:) = z;

E_second_trip_mag_e = T12_e*G23_e^2*T23_e;

% Field discrepacy
DeltaE_mag = T12_o*T23_o - T12_e*T23_e;

RP = (T12_e^2*G23_e^4*T23_e^2 + T12_o^2*G23_o^4*T23_o^2)/...
    (T12_e^2*T23_e^2 + T12_o^2*T23_o^2)*100;
    
%% plotting
% plot(z, real(E(1,:)))
% hold on
% plot(z, real(E(2,:)))
% grid on
% xlabel("Distance $z$ (m)", 'Interpreter', 'latex')
% ylabel("Amplitude $\mathbf{E}(z)$ (V/m)", 'Interpreter', 'latex')
% hold off
E_re = real(E);

figure
colormap("turbo")
xp = E_re(1,:);
yp = E_re(2,:);
zp = E_re(3,:)*1E6;
% plot3(xp, yp, zp)
patch([xp nan], [yp nan], [zp nan], [zp nan] ,'FaceColor','none', ...
    'EdgeColor','interp')
xlabel("$x$ (V/m)", 'Interpreter','latex', 'FontSize', 200)
ylabel("$y$ (V/m)", 'Interpreter','latex', 'FontSize', 200)
zlabel("$z$ ($\mu$ m)", 'Interpreter','latex', 'FontSize', 200)
title("Electric Field \textbf{E}$(z; t = 0)$ in QWP", 'Interpreter', ...
    'latex', 'FontSize', 100)
grid on
ax = gca;
ax.FontSize = 20;
ax.GridLineWidth = 2;

disp("Param")
vpa(round([L*1E6, m, DeltaE_mag, RP] ,3))





 



