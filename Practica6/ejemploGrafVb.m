clf;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=0; %%TM01
m=1;
% p01=2.4; %m=1 TM01
% p02=5.5; %m=1 TM02

TM_pnm = [
    2.4 3.8 5.1 6.18;
    5.52 7.01 8.4 9.2;
    8.65 10.17 11.6 13.01;
    11.7 13.3 14.79 16.2
    ];

TE_pnm = [
    3.8 1.8 3.05 4.2;
    7 3.3 6.7 8;
    10.17 8.5 4.96 11.34;
    13.32 11.7 13.17 14.58
];

calculoModosLP2(TM_pnm(1,2),TM_pnm(1,2),"TM_{11}","TM_{11}"); %Modo LP01
calculoModosLP4(TE_pnm(1,1),TM_pnm(1,1),TM_pnm(1,3),TM_pnm(1,3),"TE_{11}","TM_{11}","TM_{21}","TM_{21}"); %Modo LP11
calculoModosLP4(TE_pnm(1,2),TE_pnm(1,2),TM_pnm(1,4),TM_pnm(1,4),"TE_{11}","TE_{11}","TM_{31}","TM_{31}"); %Modo LP21
calculoModosLP2(TM_pnm(2,2),TM_pnm(2,2),"TM_{12}","TM_{12}"); %Modo LP02
% calculoModosLP4(TE_pnm(1,3),TE_pnm(1,3),TM_pnm(4,1),TM_pnm(4,1),"TE_{11}","TE_{11}","TM_{31}","TM_{31}"); %Modo LP31
calculoModosLP4(TE_pnm(2,1),TM_pnm(2,1),TM_pnm(2,3),TM_pnm(2,3),"TE_{02}","TM_{02}","TM_{22}","TM_{22}"); %Modo LP12




function calculoModosLP1(p1,np1)
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
    
    v=1/sqrt(miu*eps);%%velocidad de grupo de onda

    fc=(p1*c)/(2*pi*a);%% fc_TMnm

    frec_Op=fc:0.1e12:12e12;
    %frec_TM01 =1.2e12;

    lambda_TM01=v./frec_Op;

    V=2*pi*a./lambda_TM01.*sqrt(n1^2 - n2^2);%FREC NORM
    
    omega=2*pi.*frec_Op;

    Beta=omega.*sqrt(eps*miu).*sqrt(1-((fc./frec_Op).^2));
    k=omega.*sqrt(eps*miu)*n1;

    figure();
    plot(V,Beta./k);
    grid on;
    xlabel("v");
    ylabel("\beta/k");
    title("Grafica v-b para modo "+ np1);
    hold on;
end

function calculoModosLP2(p1,p2,np1,np2)
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
    
    v=1/sqrt(miu*eps);%%velocidad de grupo de onda

    fc_TM01=(p1*c)/(2*pi*a);%% fc_TMnm
    fc_TM02=(p2*c)/(2*pi*a);%% fc_TMnm

    frec_TM01=fc_TM01:0.1e12:12e12;
    frec_TM02=fc_TM02:0.1e12:12e12;
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


    figure();
    plot(V_TM01,Beta_TM01./k_TM01);
    
    hold on;
    plot(V_TM02,Beta_TM02./k_TM02,'r');
    grid on;
    xlabel("v");
    ylabel("\beta/k");
    title("Grafica v-b para modo " + np1 + " y para " + np2);
    hold on;
end

function calculoModosLP4(p1,p2,p3,p4,np1,np2,np3,np4)
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
    
    v=1/sqrt(miu*eps);%%velocidad de grupo de onda

    fc_TM01=(p1*c)/(2*pi*a);%% fc_TMnm
    fc_TM02=(p2*c)/(2*pi*a);%% fc_TMnm
    fc_TM03=(p3*c)/(2*pi*a);%% fc_TMnm
    fc_TM04=(p4*c)/(2*pi*a);%% fc_TMnm


    frec_TM01=fc_TM01:0.1e12:12e12;
    frec_TM02=fc_TM02:0.1e12:12e12;
    frec_TM03=fc_TM03:0.1e12:12e12;
    frec_TM04=fc_TM04:0.1e12:12e12;

    %frec_TM01 =1.2e12;

    lambda_TM01=v./frec_TM01;
    lambda_TM02=v./frec_TM02;
    lambda_TM03=v./frec_TM03;
    lambda_TM04=v./frec_TM04;

    V_TM01=2*pi*a./lambda_TM01.*sqrt(n1^2 - n2^2);%FREC NORM
    V_TM02=2*pi*a./lambda_TM02.*sqrt(n1^2 - n2^2);%FREC NORM
    V_TM03=2*pi*a./lambda_TM03.*sqrt(n1^2 - n2^2);%FREC NORM
    V_TM04=2*pi*a./lambda_TM04.*sqrt(n1^2 - n2^2);%FREC NORM


    omega_TM01=2*pi.*frec_TM01;
    omega_TM02=2*pi.*frec_TM02;
    omega_TM03=2*pi.*frec_TM03;
    omega_TM04=2*pi.*frec_TM04;


    Beta_TM01=omega_TM01.*sqrt(eps*miu).*sqrt(1-((fc_TM01./frec_TM01).^2));
    k_TM01=omega_TM01.*sqrt(eps*miu)*n1;

    Beta_TM02=omega_TM02.*sqrt(eps*miu).*sqrt(1-((fc_TM02./frec_TM02).^2));
    k_TM02=omega_TM02.*sqrt(eps*miu)*n1;
    
    Beta_TM03=omega_TM03.*sqrt(eps*miu).*sqrt(1-((fc_TM03./frec_TM03).^2));
    k_TM03=omega_TM03.*sqrt(eps*miu)*n1;
    
    Beta_TM04=omega_TM04.*sqrt(eps*miu).*sqrt(1-((fc_TM04./frec_TM04).^2));
    k_TM04=omega_TM04.*sqrt(eps*miu)*n1;

    figure();
    plot(V_TM01,Beta_TM01./k_TM01);
    
    hold on;
    plot(V_TM02,Beta_TM02./k_TM02);
    
    hold on;
    plot(V_TM03,Beta_TM03./k_TM03);
    
    hold on;
    plot(V_TM04,Beta_TM04./k_TM04);

    grid on;
    xlabel("v");
    ylabel("\beta/k");
    title("Grafica v-b para los modos "+ np1 + "," + np2 + "," + np3 + "," + np4);
    hold on;
end

