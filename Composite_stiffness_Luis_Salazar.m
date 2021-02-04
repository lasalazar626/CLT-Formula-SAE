%%%Written: Luis Salazar
clc
clear
%% put the carbon properties
correction = .6%
E_1= 8700000; %lb/in^2
E1=  E_1-(correction*E_1) 
E_2= 8400000; %lb/in^2
E2 = E_1-(correction*E_2);
nu_12= 0.4; %lb/in^2
nu_21= 0.38;
G_1= 4198842.5; %lb/in^2
G =  G_1 -(correction*G_1);
t= (1/3)/25.4; %in
Tension_0 = 105000;% lb/in^2
Tenion_90 = 97000; % lb/in^2
%Shear_s = 
%% Q matrix
Q_11= E1/(1-(nu_12*nu_21));
Q_22= E2/(1-(nu_12*nu_21));
Q_12= (nu_12*E2)/(1-(nu_12*nu_21));
Q_66= G;
%% Enter the ply schedule 
theta = [0,45,0];
%% Q bar Matrix
for i = 1:length(theta);
m=cos((180\pi)*theta(i));
n=sin((180\pi)*theta(i));
Q_bar_11(i)=Q_11*(m^4) + (2*(Q_12 + 2*Q_66)*(m^2)*(n^2)) + (Q_22*n^4);
Q_bar_12(i)=((Q_11+Q_22-(4*Q_66))*(m^2)*(n^2))+ (Q_12*((m^4) + (n^4)));
Q_bar_22(i)=((Q_11)*n^4)+ (2*(Q_12+2*Q_66))*m^2 * n^2 + Q_22*m^4;
Q_bar_16(i)=((Q_11-Q_12-(2*Q_66))*(m^3)*n) + ((Q_12-Q_22+(2*Q_66))*m*(n^3));
Q_bar_26(i)= (Q_11-Q_12-2*Q_66)*n^3*m + (Q_12-Q_22+2*Q_66)*n*m^3;
Q_bar_66(i)= ((Q_11+Q_22-(2*Q_12)-(2*Q_66))*(m^2)*(n^2)) + Q_66*(m^4 + n^4);
end
%% A Matrix
A_11=sum(Q_bar_11)*t;
A_12=sum(Q_bar_12)*t;
A_22=sum(Q_bar_22)*t;
A_16=sum(Q_bar_16)*t;
A_26=sum(Q_bar_26)*t;
A_66=sum(Q_bar_66)*t;
A=[A_11 A_12 A_16; A_12 A_22 A_26; A_16 A_26 A_66];
%% B Matrix
hk=0;
for k=length(theta):-1:1
    if rem(length(theta),2)==0
        hk(length(theta)+1-k)= ((t/2+ (t*(k-1))) - (t*0.5*length(theta)));
    else 
        hk(length(theta)+1-k)= ((t/2)+(t*(k-1))-(t*0.5*length(theta)));
end
end
B_11=0;
B_12=0;
B_22=0;
B_16=0;
B_26=0;
B_66=0;
for i = 1:length(theta)
B_11(i)=Q_bar_11(i)*t*hk(i);
B_12(i)=Q_bar_12(i)*t*hk(i);
B_22(i)=Q_bar_22(i)*t*hk(i);
B_16(i)=Q_bar_16(i)*t*hk(i);
B_26(i)=Q_bar_26(i)*t*hk(i);
B_66(i)=Q_bar_66(i)*t*hk(i);
% B_11=Q_bar_11(i)*t*hk(i)+B_11;
% B_12=Q_bar_12(i)*t*hk(i)+B_12;
% B_22=Q_bar_22(i)*t*hk(i)+B_22;
% B_16=Q_bar_16(i)*t*hk(i)+B_16;
% B_26=Q_bar_26(i)*t*hk(i)+B_26;
% B_66=Q_bar_66(i)*t*hk(i)+B_66;
end
B_11=round(sum(B_11));
B_12=round(sum(B_12));
B_22=round(sum(B_22));
B_16=round(sum(B_16));
B_26=round(sum(B_26));
B_66=round(sum(B_66));
B=[B_11 B_12 B_16; B_12 B_22 B_26; B_16 B_26 B_66];
%% D Matrix
for i=1:length(theta)
D_11(i)=Q_bar_11(i)*(((t^3)/12)+(t*(hk(i)^2)));
D_12(i)=Q_bar_12(i)*(((t^3)/12)+(t*(hk(i)^2)));
D_22(i)=Q_bar_22(i)*(((t^3)/12)+(t*(hk(i)^2)));
D_16(i)=Q_bar_16(i)*(((t^3)/12)+(t*(hk(i)^2)));
D_26(i)=Q_bar_26(i)*(((t^3)/12)+(t*(hk(i)^2)));
D_66(i)=Q_bar_66(i)*(((t^3)/12)+(t*(hk(i)^2)));
end
D_11=(sum(D_11));
D_12=(sum(D_12));
D_22=(sum(D_22));
D_16=(sum(D_16));
D_26=(sum(D_26));
D_66=(sum(D_66));
D=[D_11 D_12 D_16; D_12 D_22 D_26; D_16 D_26 D_66];
%% ABD Matrix
ABD= [A,B;B,D];
%% Modulus of Composite
ABD_inv=[A(2,2),A(2,3),B(1,2),B(2,2),B(2,3);A(2,3),A(3,3),B(1,3),B(2,3),B(3,3);B(1,2),B(1,3),D(1,1),D(1,2),D(1,3);B(2,2),B(2,3),D(1,2),D(2,2),D(2,3);B(2,3),B(3,3),D(1,3),D(2,3),D(3,3)];
mod_1=det(ABD);
mod_2=det(ABD_inv);
E_real=(mod_1/mod_2)/((t*(length(theta)+1)));
%% Enter Core stuff
Core= 3/4; %in
base = 275/25.4; %%in
Core_mod = 4700; % psi
I= (base/12)*((((2*t*length(theta))+Core)^3)-(Core^3)); %%Based on FSAE Structural Equivalency 
EI = E_real*I;
%% Normal Load

%Enter a load
W = [linspace(0,1750)];  %lbs make sure you change this around for your analysis
Length_supports = 400/25.4; %in
h = 2*(t*(length(theta)+1)) + Core; %in
load = [];
for i = 1:length(W)
    load(i) = W(i)*Length_supports/(4*h*base); 
end
Load = max(load) %psi
%% Deflection
distance_centroid = (t*length(theta))+Core;
deflection = [];
for i = 1:length(W)
    deflection(i) = (W(i)*(Length_supports^3))/(24*E_real*base*t*(length(theta)+1)*(distance_centroid^2))...
        + W(i)*(Length_supports)/(4*Core_mod*base*Core);
end
%% Failure (Tsai Hill and Maximum failure Critera)
Strain_matrix=[];
for i = 1:length(load)
Strain_matrix(:,i) = (inv(ABD)) * [load(i);0;0;0;0;0];
end
for i = 1:length(load)
    Strain_matrix(1:3,i) = Strain_matrix(1:3,i)*2;
end
P_strain_matrix= [];
for ss= 1:length(theta)
    for i = 1:length(load)
        m=cos((180\pi)*theta(ss));
        n=sin((180\pi)*theta(ss));
        P_strain_matrix(:,i,ss) = [m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n (m^2 - n^2)]*Strain_matrix(1:3,i);
          P_strain_matrix(:,i,ss)=P_strain_matrix(:,i,ss)*2;
    end
end
Stress=[];
for ss = 1:length(theta)
for i = 1:length(load)
    Stress(:,i,ss)=[Q_11 Q_12 0; Q_12 Q_22 0; 0 0 Q_66]*P_strain_matrix(:,i,ss);
end
end


%% Plots

plot(deflection,W)
