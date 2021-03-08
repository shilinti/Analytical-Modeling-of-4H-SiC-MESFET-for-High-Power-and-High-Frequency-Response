% This code is for simulating Ids-Vds of 4H-SIC
% Revised date: 3/7

clc; 
clear all; 
close all; 

% Fixed constants 
K=1.38066e-23; 
q = 1.60219e-19; % electron charge
Na=5e15; %cm^-3
Nd= 1.38e18; %cm^-3
T=300; 
Eg= 3.26; % band-gap eV
u=800; % electron mobility
Nc=3.25e14*sqrt(T*T*T); 
Nv= 4.8e15*sqrt(T*T*T);
Ni= sqrt(Nc*Nv)*exp(-(Eg*q)/(2*K*T)); 
Z= 1000e-4; %cm gate width
L= 1e-4; %cm gate length
pai= 3.14159; 
Phib = 1.47;  
Q=3.75e13; 
q = 1.602e-19; 
Eo= 9.6; 
E= Eo*8.854e-14; %epsilon
T2= 273; % T1 is annealing temperature
T1= 1900+ T2; % Ea is activation energy 
Ea= 40e-3; 
Rp= 0.0123e-4; %cm  
sigma= 0.0032e-4; %cm iron implantationat 140 Kev
Vbs= 0;
Vbi=  0.0253*(log((Nd*Na)/(Ni*Ni))); %builtin voltage
Vt= (K*T)/q; 
loop= 100; % loop repeating
loop1 =0.000000001; % loop repeating
%Vp = ((q*Nd*a*a)/(2*8.8542*10^(-14))) % pinch off voltage (constant)
Alpha= Rp/(2* sigma)*sqrt(pai/2);
R= E*sqrt(pai/2)/(q*Q*sigma);
C1= erf(Rp/(sigma*sqrt(2)));
A= (q*Q*u*Z)/(2*L);
Vgs = 0; % gate to source voltage (varying)
Vds = [0:1:80]; % Vds value ranging from 1 to 10
Ids(1) = 0; % intially Ids is zero
for j=1:length(Vds)
    for i = 1:loop
	Delta= Vt* log(Nc/Nd);
	A1= Vds(j)*(1+Alpha);
	B1= (2*q*Na*Na)/(3*Q*E);
	B2= (2*E)/(q*Na);
    A2= B1*((B2*(Vbi-Vbs+Vds(j)))^(3/2));
	A3= B1*((B2*(Vbi-Vbs))^(3/2));
	B3= (q*Q*sigma)/(3*E);
	B4= sqrt(2/pai);
	B5= ((Alpha*Alpha)-(2*Alpha*C1)+(C1*C1)+(R*(Vds(j)-Vgs-Phib-Delta)));
	B6= ((Alpha*Alpha)-(2*Alpha*C1)+(C1*C1)+(R*(-Vgs+Phib-Delta)));
	A4= B3*B4*(B5^(3/2));
	A5= B3*B4*(B6^(3/2));
    Ids (i+1)= A*(A1-A2+A3-A4+A5); 
    result = abs((Ids(i+1)-Ids(i))/Ids(i+1));
    if(result<loop1)
    res(j) = Ids(i+1);
    break;
    end
    end
    res(j) = (Ids(i+1))*10000;
end
plot (Vds, res)
xlabel(' Vds, Drain-to-source Voltage Vds(V)'),
ylabel(' Ids, Drain-to-source Current Ids(mA) '),
title('Vds-Ids characteristics of 4H-Silicon carbide(SIC) MESFETs')
hold on;
Vgs = 2;
Vds = [0:1:80]; % Vds value ranging from 1 to 10
Ids(1) = 0; % intially Ids is zero
for j=1:length(Vds)
    for i = 1:loop
	Delta= Vt* log(Nc/Nd);
	A1= Vds(j)*(1+Alpha);
	B1= (2*q*Na*Na)/(3*Q*E);
	B2= (2*E)/(q*Na);
    A2= B1*((B2*(Vbi-Vbs+Vds(j)))^(3/2));
	A3= B1*((B2*(Vbi-Vbs))^(3/2));
	B3= (q*Q*sigma)/(3*E);
	B4= sqrt(2/pai);
	B5= ((Alpha*Alpha)-(2*Alpha*C1)+(C1*C1)+(R*(Vds(j)-Vgs-Phib-Delta)));
	B6= ((Alpha*Alpha)-(2*Alpha*C1)+(C1*C1)+(R*(-Vgs+Phib-Delta)));
	A4= B3*B4*(B5^(3/2));
	A5= B3*B4*(B6^(3/2));
    Ids (i+1)= A*(A1-A2+A3-A4+A5);
    result = abs((Ids(i+1)-Ids(i))/Ids(i+1));
if(result<loop1)
res(j) = Ids(i+1);
break;
end
    end
    res(j) = (Ids(i+1))*10000;
end
plot (Vds, res,'g')
hold on;
Vgs = 4;
Vds = [0:1:80]; % Vds value ranging from 1 to 10
Ids(1) = 0; % intially Ids is zero
for j=1:length(Vds)
    for i = 1:loop
	Delta= Vt* log(Nc/Nd);
	A1= Vds(j)*(1+Alpha);
	B1= (2*q*Na*Na)/(3*Q*E);
	B2= (2*E)/(q*Na);
    A2= B1*((B2*(Vbi-Vbs+Vds(j)))^(3/2));
	A3= B1*((B2*(Vbi-Vbs))^(3/2));
	B3= (q*Q*sigma)/(3*E);
	B4= sqrt(2/pai);
	B5= ((Alpha*Alpha)-(2*Alpha*C1)+(C1*C1)+(R*(Vds(j)-Vgs-Phib-Delta)));
	B6= ((Alpha*Alpha)-(2*Alpha*C1)+(C1*C1)+(R*(-Vgs+Phib-Delta)));
	A4= B3*B4*(B5^(3/2));
	A5= B3*B4*(B6^(3/2));
    Ids (i+1)= A*(A1-A2+A3-A4+A5);
    result = abs((Ids(i+1)-Ids(i))/Ids(i+1));
if(result<loop1)
res(j) = Ids(i+1);
break;
end
    end
    res(j) = (Ids(i+1))*10000;
end
plot (Vds,res,'r')
hold on;
Vgs = 6;
Vds = [0:1:80]; % Vds value ranging from 1 to 10
Ids(1) = 0; % intially Ids is zero
for j=1:length(Vds)
    for i = 1:loop
	Delta= Vt* log(Nc/Nd);
	A1= Vds(j)*(1+Alpha);
	B1= (2*q*Na*Na)/(3*Q*E);
	B2= (2*E)/(q*Na);
    A2= B1*((B2*(Vbi-Vbs+Vds(j)))^(3/2));
	A3= B1*((B2*(Vbi-Vbs))^(3/2));
	B3= (q*Q*sigma)/(3*E);
	B4= sqrt(2/pai);
	B5= ((Alpha*Alpha)-(2*Alpha*C1)+(C1*C1)+(R*(Vds(j)-Vgs-Phib-Delta)));
	B6= ((Alpha*Alpha)-(2*Alpha*C1)+(C1*C1)+(R*(-Vgs+Phib-Delta)));
	A4= B3*B4*(B5^(3/2));
	A5= B3*B4*(B6^(3/2));
    Ids (i+1)= A*(A1-A2+A3-A4+A5);
    result = abs((Ids(i+1)-Ids(i))/Ids(i+1));
if(result<loop1)
res(j) = Ids(i+1);
break;
end
    end
    res(j) = (Ids(i+1))*10000;
end
plot (Vds, res,'k')
hold on;
hleg1=legend('Vgs=0v','Vgs=2v','Vgs=4v','Vgs=6v');
