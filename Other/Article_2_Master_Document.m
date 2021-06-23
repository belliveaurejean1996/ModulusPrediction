clear all
clc

%% Inputs

% Mechanical Propreties
Ef=228e9;
Em=3.7e9;

%Image Size
At=10000*6000;

%% Import Data

Data = xlsread('P4_1T.xlsx');

%% Cleaning bad detections: remove 25% outliers from either side of the minor axis

Keep=median(Data(:,5))-0.25*median(Data(:,5))<Data(:,5) & Data(:,5)<median(Data(:,5))+0.25*median(Data(:,5));
Data=Data(Keep,:);

%% Axis Transformations

% calcul de alpha
alpha=wrapTo360(Data(:,7)-90);
corection = alpha > 180;
alpha(corection) = alpha(corection)-180;
corection2 = alpha > 90;
alpha(corection2) = alpha(corection2)-180;

% calcul de beta
beta=acosd(Data(:,5)./Data(:,4));

% vecteur unitaire dans le système d'axe x'y'z'
v=transpose([sind(beta).*sind(alpha),sind(beta).*cosd(alpha),cosd(beta)]);

% transformation au système d'axe xyz
% rotation de l'axe x par -45 deg
R=[1 0 0;
   0 cosd(45) -sind(45);
   0 sind(45) cosd(45)];
w=R*v;

% calcul de Phi
Phi=transpose(acosd(w(3,:)));

% calcul de theta
Theta=transpose(acosd(w(2,:)));

%% Modulus Calculation

%Preperatory calcultations
Af=(Data(:,4)/2).*(Data(:,5)/2)*pi;
Aft=sum(Af);
Am=At-Aft;
Efx=Ef*(cosd(Phi)).^2;

% Weighted Average
E=(sum(Efx.*Af)+Em*Am)/At;

