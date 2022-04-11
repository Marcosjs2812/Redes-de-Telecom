clear all;
clf;

% % % Simulación del modelo de haz de luz
c = 3e8;
eps0 = 8.85* (10^-12); %%Epsilon al vacío
epsr = 4.7; %Permitividad relativa
eps = eps0*epsr;

miu0 = 4*pi*(10^-7);
miur = 1;
miu = miu0*miur;

lambda = 850e-9; %Ventana de operación 1600nm
f0 = c/lambda; %5GHz
omega = 2*pi*f0;
k = omega * sqrt(miu*eps);
%lambda = c/fc;
beta = 2*pi/lambda; %En el espacio libre


%%Para el núcleo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1; %Orden de Bessel
%Obtener las gráficas sobrepuestas de la sol R, con n = 0,1,2,3(hold on)
z = 3;
t0 = 3;

a = 50e-6; %radio del núcleo
r = 0 : 0.5e-6 : a; %radio

h1 = sqrt(k^2-beta^2);
X = r.*h1;

C5 = 1; %20mW laser power
R = C5 * besselj(n,X);

figure(1);
plot(r,(R));
xlabel('r[m]');
ylabel('R');
title("R(r) = C_5J_n(h1*r) function n = "+n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = 0: 0.188: 6*pi;
C3 = 1;
C4 = 1;
PHI = C3.*cos(n*phi) + C4.*sin(n*phi);

%Obtener las gráficas sobrepuestas de la sol PHI, con n = 0,1,2,3(hold on)
figure(2);
plot(phi,(PHI));
xlabel('\phi[rad]');
ylabel('\Phi');
title("\Phi(r) = C_3cos(n\phi)+C_4sen(n\phi) function with n=" + n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C1 = 1;
alpha = 0;
% z = 0: 0.1: 10;
gamma = alpha + i*beta;
Z = C1* exp(-gamma*z);
Z = real(Z);

figure(3);
plot(z,(Z));
xlabel('z[m]');
ylabel('Z');
title("Z(z) = C_1exp(-gamma*Z) ");

%Con parámetros r = 0 -> 50micras, phi = 0->360°, z = 0->10m
[PHI1,R1] = meshgrid(phi,r);

for i = 1:length(R)
   for j = 1:length(PHI)
      Ez(i,j) = R(i)*PHI(j)*Z; 
   end
end

[x1,y1,z1] = pol2cart(PHI1,R1,Ez);

figure(4)
mesh(x1,y1,z1);
xlabel('r[m]');
ylabel('\phi[°]');
zlabel('E_z');
title("Campo eléctrico E_z, z= "+z+" m, t= "+t0+" s, n= "+n);


%%Aplicando las condiciones en frontera para el modo TMnm (Hz=0,fc=)
f0 = 2000e9;
n = 0;
v = 1/sqrt(miu*eps);
p01 = 2.4;
fc_TM01 = (p01*v)/2*pi*a;
lambdac_TM01 = v/fc_TM01;
BetaLambda_TM01 = omega*sqrt(eps*miu)*sqrt(1-((fc_TM01/f0)^2));
omegac = 2*pi*fc_TM01;
kc = omegac*sqrt(eps*miu);
h_TM01 = kc;
X_TM01 = h_TM01.*r;

for ind = 1:length(X_TM01)
    for jnd = 1:length(phi)
        Ez_TM01(ind,jnd) = besselj(n,X_TM01(ind))*(cos(n*phi(jnd)) + sin(n*phi(jnd))) * exp(-BetaLambda_TM01*z); 
    end
end


[x3,y3,z3] = pol2cart(PHI1,R1,Ez_TM01);

figure(5);
mesh(x3,y3,z3);


