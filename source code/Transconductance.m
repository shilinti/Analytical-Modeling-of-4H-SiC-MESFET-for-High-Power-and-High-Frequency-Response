% This code is for simulating Transconductance-Vgs of 4H-SIC
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
a= 0.5322e-4; %active channel thickness
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
Delta= Vt* log(Nc/Nd);
Vp = ((q*Nd*a*a)/(2*8.8542*10^(-14))) % pinch off voltage (constant)
%Alpha= Rp/(2* sigma)*sqrt(pai/2);
%R= E*sqrt(pai/2)/(q*Q*sigma);
S= Rp1/(2*(sqrt(sigma^2+2*D*t)));
Alpha= S*(sqrt(pai/2));
loop= 180; % loop repeating
loop1 =0.000000001; % loop repeating
Vds = 4; % gate to source voltage (varying)
Vgs = [-15:0.1:3]; % Vds value ranging from -6 to 6
gm(1) = 0; % intially Ids is zero
Q=3.75e13;
for j=1:length(Vgs)
    for i = 1:loop
	V1= q*(Q^2)/(8*Na*E);
	V2= q*Q*sqrt(sigma^2+2*D*t)/sqrt(2*pai*E);
	Ap= (2*Na/Q)*sqrt((2*E/(q*Na))*(Vbi-Vbs+Vp));
	% below is formular of Vth
	B1= Phib-Delta;
	C1= q*Q*Rp1/(2*E);
	C2= erf(Rp1/(sqrt(2)*(sigma^2+2*D*t)));
	D1= 2*E/(q*Na);
	C3= (2*Na/Q)*(D1*(Vbi-Vbs)^(1/2));
	B2= C1*(C2+1-C3);
	C4= (q*Q*sqrt(sigma^2+2*D*t))/(E*sqrt(2*pai));
	C5= exp(-(Rp1)^2/(2*sigma));
	B3= C4*C5;
	Vth= B1-B2-B3;
	% Beloe is formular for gm
	B4= q*u*Z*Q/(4*L);
	B5= 1/(V1*Ap);
	B6= 1/(V2*(Alpha+1-Ap));
	B7= Vgs(j)-Vth;
	gm(i+1)= B4*(B5+B6)*B7;
    result = abs((gm(i+1)-gm(i))/gm(i+1));
    if(result<loop1)
    res(j) = gm(i+1);
    break;
    end
    end
    res(j) = (gm(i+1));
end
plot (Vgs, res)
xlabel(' Drain-to-source Voltage Vds, (V)'),
ylabel(' Transconductance, gm, (siemens)'),
 ylim([0 inf])
title('Vds-gm characteristics of 4H-Silicon carbide(SIC) MESFETs')
hold on;
Q=2e13;
for j=1:length(Vgs)
    for i = 1:loop
	V1= q*(Q^2)/(8*Na*E);
	V2= q*Q*sqrt(sigma^2+2*D*t)/sqrt(2*pai*E);
	Ap= (2*Na/Q)*sqrt((2*E/(q*Na))*(Vbi-Vbs+Vp));
	% below is formular of Vth
	B1= Phib-Delta;
	C1= q*Q*Rp1/(2*E);
	C2= erf(Rp1/(sqrt(2)*(sigma^2+2*D*t)));
	D1= 2*E/(q*Na);
	C3= (2*Na/Q)*(D1*(Vbi-Vbs)^(1/2));
	B2= C1*(C2+1-C3);
	C4= (q*Q*sqrt(sigma^2+2*D*t))/(E*sqrt(2*pai));
	C5= exp(-(Rp1)^2/(2*sigma));
	B3= C4*C5;
	Vth= B1-B2-B3;
	% Beloe is formular for gm
	B4= q*u*Z*Q/(4*L);
	B5= 1/(V1*Ap);
	B6= 1/(V2*(Alpha+1-Ap));
	B7= Vgs(j)-Vth;
	gm(i+1)= B4*(B5+B6)*B7;
    result = abs((gm(i+1)-gm(i))/gm(i+1));
    if(result<loop1)
    res(j) = gm(i+1);
    break;
    end
    end
    res(j) = (gm(i+1));
end
plot (Vgs, res, 'g')
hold on;
Q=1.5e13;
for j=1:length(Vgs)
    for i = 1:loop
	V1= q*(Q^2)/(8*Na*E);
	V2= q*Q*sqrt(sigma^2+2*D*t)/sqrt(2*pai*E);
	Ap= (2*Na/Q)*sqrt((2*E/(q*Na))*(Vbi-Vbs+Vp));
	% below is formular of Vth
	B1= Phib-Delta;
	C1= q*Q*Rp1/(2*E);
	C2= erf(Rp1/(sqrt(2)*(sigma^2+2*D*t)));
	D1= 2*E/(q*Na);
	C3= (2*Na/Q)*(D1*(Vbi-Vbs)^(1/2));
	B2= C1*(C2+1-C3);
	C4= (q*Q*sqrt(sigma^2+2*D*t))/(E*sqrt(2*pai));
	C5= exp(-(Rp1)^2/(2*sigma));
	B3= C4*C5;
	Vth= B1-B2-B3;
	% Beloe is formular for gm
	B4= q*u*Z*Q/(4*L);
	B5= 1/(V1*Ap);
	B6= 1/(V2*(Alpha+1-Ap));
	B7= Vgs(j)-Vth;
	gm(i+1)= B4*(B5+B6)*B7;
    result = abs((gm(i+1)-gm(i))/gm(i+1));
    if(result<loop1)
    res(j) = gm(i+1);
    break;
    end
    end
    res(j) = (gm(i+1));
end
plot (Vgs, res, 'r');
hold on;
hleg1=legend('Q=3.75e13','Q=2e13','Q=1.5e13');
hold off;
