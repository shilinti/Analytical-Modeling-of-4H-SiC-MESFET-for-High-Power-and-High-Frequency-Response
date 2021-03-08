% This code is for simulating Ron-resistance of 4H-SIC
% Revised date:3/7

clc; 
clear all; 
close all; 

% Fixed constants 
K=1.38066e-23; 
q = 1.60219e-19; % electron charge
Na=5e15; %cm^-3
Nb=5e15;%cm^-3
Nd= 1.38e18; %cm^-3
u=800; % electron mobility
Z= 1000e-4; %cm gate width
L= 1e-4; %cm gate length
a= 0.5322e-4; %active channel thickness
Phib = 1.47;   
q = 1.602e-19;
Eo= 9.6; 
E= Eo*8.854e-14; %epsilon
D=4.18e-19; % cm^2/s diffusion coefficient at 2200 c
t= 1800; % s diffusion time 1800s, 60 mins.
Ea= 40e-3;  % Ea is activation energy 
Rp= 0.0123e-4; %cm
Vbs= 0;

%Alpha= Rp/(2* sigma)*sqrt(pai/2);
%R= E*sqrt(pai/2)/(q*Q*sigma);
loop= 50; % loop repeating
loop1 =0.000000001; % loop repeating
Nd = [0.1e16:0.01e16:0.6e16]; % Vds value ranging from -6 to 6
for j=1:length(Nd)
    for i = 1:loop
	Wpp= 1.82e11*(Nd(j)^(-7/8));
	Ron(i+1)= Wpp/(q*u*Nd(j));
    result = abs((Ron(i+1)-Ron(i))/Ron(i+1));
    if(result<loop1)
    res(j) = Ron(i+1);
    break;
    end
    end
    %res(j) = Cgs(i+1);
end
plot (Nd,res,'r')
grid on;
xlabel(' N-drift layer concentration, Nd (/cm^-3)'),
ylabel(' Spercific on-Resistance, Ron-sp (ohm.cm^2) '),
title('Specific Ron versus Nd  characteristics of 4H-Silicon carbide(SIC) MESFETs')
hold on;
hleg1=legend('Z= 1000um L = 1um Na= 5e15 /cm^3');
hold off;

