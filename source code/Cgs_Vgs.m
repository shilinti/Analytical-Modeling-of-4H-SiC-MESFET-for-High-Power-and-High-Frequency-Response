% This code is for simulating Cgs-Vgs of 4H-SIC
% Revised date: 3/7

clc; 
clear all; 
close all; 

% Fixed constants 
K=1.38066e-23; 
q = 1.60219e-19; % electron charge
Na=5e15; %cm^-3
Nb=5e15;%cm^-3
Nd= 1.38e18; %cm^-3
T=2200; 
Eg= 3.26; % band-gap eV
u=800; % electron mobility
Nc=3.25e14*sqrt(T*T*T); 
Nv= 4.8e15*sqrt(T*T*T);
Ni= sqrt(Nc*Nv)*exp(-(Eg*q)/(2*K*T)); 
Z= 1000e-4; %cm gate width
L= 1e-4; %cm gate length
pai= 3.14159; 
Phib = 1.47;   
q = 1.602e-19;
sigma= 0.0032e-4; %cm iron implantationat 140 Kev 
Eo= 9.6; 
E= Eo*8.854e-14; %epsilon
D=4.18e-19; % cm^2/s diffusion coefficient at 2200 c
t= 1800; % s diffusion time 1800s, 60 mins.
Ea= 40e-3;  % Ea is activation energy 
Rp= 0.0123e-4; %cm
Rp1 = Rp + sqrt(sigma*sigma+2*D*t); 
Vbs= 0;
Vbi= 0.0253*(log((Nd*Na)/(Ni*Ni))); %builtin voltage 3
Vt= (K*T)/q; 
loop= 60; % loop repeating
loop1 =0.000000001; % loop repeating
%Vp = ((q*Nd*a*a)/(2*8.8542*10^(-14))) % pinch off voltage (constant)
Vds = 4; % Drain to source voltage
Vgs = [-6:0.1:0]; % Vgs value ranging from 1 to 80
Q=3.75e13;
%Cgs(1)= 5.484e-11;
for j=1:length(Vgs)
    for i = 1:loop
	S= Rp1/(2*(sqrt(sigma^2+2*D*t)));
	Alpha= S*(sqrt(pai/2));
	N= Vbi-Vgs(j)+Vds;
	R= Vbi-Vgs(j);
	M= (4*Alpha*E)/(q*Q*Rp1);
	Alpha= S*(sqrt(pai/2));
	A3= 1- exp(-(S*S))- 2*Alpha*erf(S);
	A1= Alpha^2+A3;
	A2= erf(S)-Alpha;	
	C1= (M*R+A1)^(1/2); %1
	C2= sqrt(N)-sqrt(R);
	C3= sqrt(N)/(2*(C2^2)); %1
	C4= (M*N+A1)^(1/2)/(sqrt(R)); %1
	C5= (M*C2)/C1;
	C6= C1/sqrt(R); %1
	C7= (-1/Q)*(((2*Nb*E)/q)^(1/2))*(2*C2/sqrt(R)+Vds/sqrt(R*N)); 
	B1= (q*Q*Z*L)/2;
	B2= M/(2*C1);
	B3= C3*(C4-C5-C6+C7);
	B4= (pai/2)*E*Z;
    Cgs(i+1)=B1*(B2+B3)+B4; 
    result = abs((Cgs(i+1)-Cgs(i))/Cgs(i+1));
    if(result<loop1)
    res(j) = Cgs(i+1);
    break;
    end
    end
    %res(j) = Cgs(i+1);
end
plot (Vgs,res)
xlabel(' Vgs, Gate-to-source voltage Vgs(V)'),
ylabel(' Cgs, Gate-to-source Capacitance Cgs(F) '),
title('Vgs-Cgs characteristics of 4H-Silicon carbide(SIC) MESFETs')
hold on;
Q=2e13;
for j=1:length(Vgs)
    for i = 1:loop
	S= Rp1/(2*(sqrt(sigma^2+2*D*t)));
	Alpha= S*(sqrt(pai/2));
	N= Vbi-Vgs(j)+Vds;
	R= Vbi-Vgs(j);
	M= (4*Alpha*E)/(q*Q*Rp1);
	Alpha= S*(sqrt(pai/2));
	A3= 1- exp(-(S*S))- 2*Alpha*erf(S);
	A1= Alpha^2+A3;
	A2= erf(S)-Alpha;	
	C1= (M*R+A1)^(1/2); %1
	C2= sqrt(N)-sqrt(R);
	C3= sqrt(N)/(2*(C2^2)); %1
	C4= (M*N+A1)^(1/2)/(sqrt(R)); %1
	C5= (M*C2)/C1;
	C6= C1/sqrt(R); %1
	C7= (-1/Q)*(((2*Nb*E)/q)^(1/2))*(2*C2/sqrt(R)+Vds/sqrt(R*N)); 
	B1= (q*Q*Z*L)/2;
	B2= M/(2*C1);
	B3= C3*(C4-C5-C6+C7);
	B4= (pai/2)*E*Z;
    Cgs(i+1)=B1*(B2+B3)+B4; 
    result = abs((Cgs(i+1)-Cgs(i))/Cgs(i+1));
    if(result<loop1)
    res(j) = Cgs(i+1);
    break;
    end
    end
%     res(j) = Cgs(i+1);
end
plot (Vgs,res,'g')
hold on;
Q=1.5e13;
for j=1:length(Vgs)
    for i = 1:loop
	S= Rp1/(2*(sqrt(sigma^2+2*D*t)));
	Alpha= S*(sqrt(pai/2));
	N= Vbi-Vgs(j)+Vds;
	R= Vbi-Vgs(j);
	M= (4*Alpha*E)/(q*Q*Rp1);
	Alpha= S*(sqrt(pai/2));
	A3= 1- exp(-(S*S))- 2*Alpha*erf(S);
	A1= Alpha^2+A3;
	A2= erf(S)-Alpha;	
	C1= (M*R+A1)^(1/2); %1
	C2= sqrt(N)-sqrt(R);
	C3= sqrt(N)/(2*(C2^2)); %1
	C4= (M*N+A1)^(1/2)/(sqrt(R)); %1
	C5= (M*C2)/C1;
	C6= C1/sqrt(R); %1
	C7= (-1/Q)*(((2*Nb*E)/q)^(1/2))*(2*C2/sqrt(R)+Vds/sqrt(R*N)); 
	B1= (q*Q*Z*L)/2;
	B2= M/(2*C1);
	B3= C3*(C4-C5-C6+C7);
	B4= (pai/2)*E*Z;
    Cgs(i+1)=B1*(B2+B3)+B4; 
    result = abs((Cgs(i+1)-Cgs(i))/Cgs(i+1));
    if(result<loop1)
    res(j) = Cgs(i+1);
    break;
    end
    end
    %res(j) = Cgs(i+1);
end
plot (Vgs,res,'r')
hold on;
hleg1=legend('Q=3.75e13','Q=2e13','Q=1.5e13');
hold off;
