%% Initialisation
clear all %#ok<CLALL>
clc
%#ok<*NOPTS>

%%  Data input
Alpha = 45; % Degres

% parameters
a = [1;1];
b = [0.5;0.5];
RA = [30;150];

% % Material properties
% E1 = 135e9;     % Pa
% E2 = 10e9;      % Pa
% G12 = 5.2e9;    % Pa
% Volf = 0.59;    % Fiber volume fraction
% Nu12 = 0.22;    % Estimation from composite book
% Nu21 = (Nu12 * E2)/E1;
% 
% % colors for Theta
% div_max=90;
% div=-div_max:(10):div_max;
% n_div=length(div)-1;
% color=jet(n_div);

%% Parametres calculation
% Calcul de ThP
PhiP = wrapTo360(180-RA);
corection = PhiP > 270; %#ok<*CHAIN>
PhiP(corection) = PhiP(corection) - 360;
corection2 = PhiP > 180;
PhiP(corection2) = PhiP(corection2)-180;
corection3 = PhiP > 90;
PhiP(corection3) = PhiP(corection3)-180;

% Th directionnal index
directionindex = ones(length(PhiP),1);
signenegatif = PhiP > 0;
directionindex(signenegatif) = -1;

% Calcul de ThP
ThP=acosd(b./a);

% Vecteur unitaire dans le système d'axe x'y'z'
UP=transpose([cosd(ThP),...
              sind(ThP).*cosd(PhiP),...
              sind(ThP).*sind(PhiP)]);

% Transformation au système d'axe xyz
R=[cosd(Alpha)  0  -sind(Alpha);
   0            1  0;
   sind(Alpha) 0  cosd(Alpha)];

% vecteur unitaire dans système d'axe xyz
U=R*UP;

% Calcul de Th
Th = acosd(U(1,:)');
Th = Th.*directionindex;

% Calcul de theta
Phi = (atan2(U(3,:),U(2,:))) * 180/pi;


th = acosd(cosd(Alpha).*cosd(ThP) - sind(Alpha).*sind(ThP).*sind(PhiP))
phi = atand((cosd(Alpha).*sind(ThP).*sind(PhiP) + sind(Alpha).*cosd(ThP))/(sind(ThP).*cosd(PhiP)))
