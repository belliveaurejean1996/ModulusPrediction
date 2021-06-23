%% Initialisation
clear all %#ok<CLALL>
clc
%#ok<*NOPTS>

%%  Data input
% Data from Figi
% Symetric composite
% file = 'D:\Université\Matrise\Article - Modulus Prediction\DataOutput\symetrique\SymetricLong.xlsx';
% Gamma = 51;
% E1 = 130e9;% Pa

% DLF composite
file = 'D:\Université\Matrise\Article - Modulus Prediction\DataOutput\R1p_200x';
Gamma = 45; 
E1 = 135e9;% Pa

% Read excel file
Data = xlsread(file);

% Néttoyage de 25% de la médianne
Keep=median(Data(:,5))-0.25*median(Data(:,5))<Data(:,5) & Data(:,5)<median(Data(:,5))+0.25*median(Data(:,5));
Data=Data(Keep,:);

% parameters
x = Data(:,2);
y = Data(:,3);
b = Data(:,4);
a = Data(:,5);
RA = Data(:,7);

% Material properties
E2 = 10e9;      % Pa
G12 = 5.2e9;    % Pa
Volf = 0.59;    % Fiber volume fraction
Nu12 = 0.22;    % Estimation from composite book
Nu21 = (Nu12 * E2)/E1;

% colors for phi
div_max=90;
div=-div_max:(10):div_max;
n_div=length(div)-1;
color=jet(n_div);

%% Parametres calculation
% Calcul de alpha
alpha=wrapTo360(RA-90);
corection = alpha > 180;
alpha(corection) = alpha(corection)-180;

% Phi directionnal index
directionindex=ones(length(alpha),1);
signenegatif=alpha < 90;
directionindex(signenegatif)=-1;

% Resume alpha
corection2 = alpha > 90;
alpha(corection2) = alpha(corection2)-180;

% Calcul de beta
beta=acosd(a./b);

% Vecteur unitaire dans le système d'axe x'y'z'
v=transpose([sind(beta).*sind(alpha),sind(beta).*cosd(alpha),cosd(beta)]);

% Transformation au système d'axe xyz
% Rotation de l'axe x par -45 deg
R=[1    0               0;
   0    cosd(Gamma)     -sind(Gamma);
   0    sind(Gamma)     cosd(Gamma)];
w=R*v;

% Calcul de Phi
Phi=acosd(w(3,:)');
Phi=Phi.*directionindex;

% Calcul de theta
Theta=acosd(w(2,:)');

% Delete angles larger than thresshold
TH = 65; %degres
rowsToDelete = Phi < -TH | Phi > TH;
Phi(rowsToDelete) = [];
Theta(rowsToDelete) = [];
x(rowsToDelete) = []; %#ok<*NASGU>
y(rowsToDelete) = [];
a(rowsToDelete) = [];
b(rowsToDelete) = [];


%% figures 
figure(1)
histogram(a)
xlabel('Minor axis (Pixels)')
ylabel('Number of fibres')

figure(2)
histogram(b)
xlabel('Major axis (Pixels)')
ylabel('Number of fibres')


figure(3)
histogram(Theta-90)
xlabel('Phi (Degrees)')
ylabel('Number of fibres')

figure(4)
histogram(Phi)
xlabel('Theta (Degrees)')
ylabel('Number of fibres')

Theta = Theta - 90;
save('D2_variable.mat','a','b','Theta','Phi')








