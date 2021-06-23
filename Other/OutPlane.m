%% initialisation
clear all %#ok<CLALL>
clc
%#ok<*NOPTS>

%% Notes

% -------------------------------------------------------------------------
% Basé sur les résultas de ton analyse experimental (Histograme) le module
% de 40 GPa est presque exactement a l'endroit ou on veux.
% -------------------------------------------------------------------------

%% Inputs

% Data from Figi
file = 'D:\Université\Matrise\Article\DataOutput\symetrique\SymetricLong.xlsx';
Data = xlsread(file);

% Material properties
E1 = 135e9;     % Pa
E2 = 10e9;      % Pa
G12 = 5.2e9;    % Pa
Volf = 0.59;    % Fiber volume fraction

Nu12 = 0.22; % Estimation from composite book
Nu21 = (Nu12 * E2)/E1;

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
R=[1    0           0;
   0    cosd(45)    -sind(45);
   0    sind(45)    cosd(45)];
w=R*v;

% calcul de Phi
Phi=transpose(acosd(w(3,:)));
Phi0=Phi;
Phi=Phi.*directionindex;

% calcul de theta
Theta=transpose(acosd(w(2,:)));

%% tensor

L = 12.7; %length of fibre (dosn't car because unit vector)

x = L .* sind(Phi) .* sind(Theta);
y = L .* sind(Phi) .* cosd(Theta);
z = L .* cosd(Phi);

table(Theta, Phi, x, y, z);

Vector = [x,y,z];
n = length(x);

for i = 1:n
    VN(i) = norm(Vector(i,:)); %#ok<SAGROW>
end

Xu = x./VN';
Yu = y./VN';
Zu = z./VN';

UnitV = abs([Xu, Yu, Zu]);

a = eye(3);
a(1,1) = mean(abs(Xu));
a(2,2) = mean(abs(Yu));
a(3,3) = mean(abs(Zu))


