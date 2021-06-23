%% initialisation
clear all %#ok<CLALL>
clc
%#ok<*NOPTS>

%% Inputs

% Data from Figi
file = 'D:\Université\Matrise\Article\DataOutput\symetrique\SymetricLong.xlsx';
Data = xlsread(file);
Volf = 0.59;    % Fiber volume fraction

% %colors for phi
div_max=90;
div=-div_max:(5):div_max;
n_div=length(div)-1;
color=jet(n_div);

%% Fiji operation
% nettoyage de 25% de la médianne
Keep=median(Data(:,5))-0.25*median(Data(:,5))<Data(:,5) & Data(:,5)<median(Data(:,5))+0.25*median(Data(:,5));
Data=Data(Keep,:);

% calcul de alpha
alpha=wrapTo360(Data(:,7)-90);
corection = alpha > 180;
alpha(corection) = alpha(corection)-180;

%for later
directionindex=ones(length(alpha),1);
signenegatif=alpha < 90;
directionindex(signenegatif)=-1;

%Resume alpha
corection2 = alpha > 90;
alpha(corection2) = alpha(corection2)-180;

% calcul de beta
beta=acosd(Data(:,5)./Data(:,4));

% vecteur unitaire dans le système d'axe x'y'z'
v=transpose([sind(beta).*sind(alpha),sind(beta).*cosd(alpha),cosd(beta)]);

% transformation au système d'axe xyz
% rotation de l'axe x par -45 deg
Gamma = 51;
R=[1    0           0;
   0    cosd(Gamma)    -sind(Gamma);
   0    sind(Gamma)    cosd(Gamma)];
w=R*v;

% calcul de Phi
Phi=transpose(acosd(w(3,:)));
Phi0=Phi;
Phi=Phi.*directionindex;

% calcul de theta
Theta=transpose(acosd(w(2,:)));

%% trace ellipse on image



















