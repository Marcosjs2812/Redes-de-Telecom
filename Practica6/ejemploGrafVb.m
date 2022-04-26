clf;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=3e8;
eps0=8.85*(10^-12);
epsr=4.7;
eps=eps0*epsr;
miu0=4*pi*(10^-7);
miur=1;
miu=miu0*miur;
a=50e-6;
n1=1.5;
n2=1.47;




n=0; %%TM01
m=1;
p01=2.4; %m=1 TM01
p02=5.5; %m=1 TM02

%%%%%%%%%%%%%%%%%%%%%
v=1/sqrt(miu*eps);%%velocidad de grupo de onda

fc_TM01=(p01*c)/(2*pi*a);%% fc_TMnm
fc_TM02=(p02*c)/(2*pi*a);%% fc_TMnm

frec_TM01=fc_TM01:0.1e12:10e12;
frec_TM02=fc_TM02:0.1e12:10e12;
%frec_TM01 =1.2e12;

lambda_TM01=v./frec_TM01;
lambda_TM02=v./frec_TM02;

V_TM01=2*pi*a./lambda_TM01.*sqrt(n1^2 - n2^2);%FREC NORM
V_TM02=2*pi*a./lambda_TM02.*sqrt(n1^2 - n2^2);%FREC NORM

omega_TM01=2*pi.*frec_TM01;
omega_TM02=2*pi.*frec_TM02;

Beta_TM01=omega_TM01.*sqrt(eps*miu).*sqrt(1-((fc_TM01./frec_TM01).^2));
k_TM01=omega_TM01.*sqrt(eps*miu)*n1;

Beta_TM02=omega_TM02.*sqrt(eps*miu).*sqrt(1-((fc_TM02./frec_TM02).^2));
k_TM02=omega_TM02.*sqrt(eps*miu)*n1;


figure(1);
plot(V_TM01,Beta_TM01./k_TM01);
grid on;
xlabel("v");
ylabel("\beta/kdis");
title("Grafica v-b para modo TM_{01}");
figure(2);
plot(V_TM02,Beta_TM02./k_TM02,'r');
grid on;
xlabel("v");
ylabel("\beta/k");
title("Grafica v-b para modo TM_{02}");
