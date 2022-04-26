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

%%Para el núcleo
a = 50e-6; %radio del núcleo
r = 0 : 0.5e-6 : a; %radio
phi = 0: 0.188: 6*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[PHI1,R1] = meshgrid(phi,r);
v = 1/sqrt(miu*eps);
%%Aplicando las condiciones en frontera para el modo TMnm (Hz=0,fc=)
pnm = [
    2.4 3.8 5.1 6.18;
    5.52 7.01 8.4 9.2;
    8.65 10.17 11.6 13.01;
    11.7 13.3 14.79 16.2
    ];

n = 0;
m = 1;

for n = 1: 4
   for m = 1: 3
      
        fc_TMnm = (pnm(n,m)*v)/(2*pi*a);
        f0 = fc_TMnm + 200e9; %Para f0>fc_TMnm pero f0 < al siguiente modo fcTMnm
        omega = 2*pi*f0;

        lambdac_TMnm = v/(fc_TMnm);
        BetaLambda_TMnm = omega*sqrt(eps*miu)*sqrt(1-((fc_TMnm/f0)^2));
        omegac = 2*pi*fc_TMnm; %Frecuencia de corte (Pasa alto) del modo TM01
        kc = pnm(n,m)/a;
        h_TMnm = kc;
        X_TMnm = h_TMnm.*r;
        z_TMnm = lambdac_TMnm*2; %Multiplos del doble de la longitud de onda, z = 0-5m,0-10m,0-100m (dependiendo)

        t_TMnm = 1;
        gama = i*BetaLambda_TMnm;

        coeficiente = -i*(f0/fc_TMnm)*(sqrt(1-(fc_TMnm/f0)^2));
        for t_TMnm = 1 : 1: 10
            z_TMnm = lambdac_TMnm*t_TMnm; %Multiplos del doble de la longitud de onda, z = 0-5m,0-10m,0-100m (dependiendo)
            for ind = 1:length(X_TMnm)
                for jnd = 1:length(phi)
                    Hr_TMnm(ind,jnd) = coeficiente*(besselj((n-1),X_TMnm(ind)))*(cos((n-1)*phi(jnd)) - sin((n-1)*phi(jnd))) * exp(i*BetaLambda_TMnm*z_TMnm); 
                end
            end


            [x3,y3,z3] = pol2cart(PHI1,R1,real(Hr_TMnm));

            mesh(x3,y3,z3);
            view(45,45);
            xlabel('r[m]');
            ylabel('\phi [°]');
            zlabel("H_r{TM"+(n-1)+m+"}");
            title("H_r{TM"+(n-1)+m+"}, z ="+z_TMnm+"m, t = "+t_TMnm+"s");
            drawnow limitrate
            pause(0.5)
       end
       figure()
   end
end




