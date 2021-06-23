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
file = 'D:\Université\Matrise\Article - Modulus Prediction\DataOutput\R1p_200x';
Data = xlsread(file);
Gamma = 45; 
E1 = 135e9;% Pa

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

% Keep=Data(:,2)<1e5;
% Data=Data(Keep,:);

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

%% Modified laminate Theorie

% Matrix Q (Doesn't change)
Q = [E1/(1-Nu12*Nu21)           (Nu21*E1)/(1-Nu12*Nu21) 0;
     (Nu12*E2)/(1-Nu12*Nu21)    E2/(1-Nu12*Nu21)        0;
     0                          0                       G12];
 
% Prealocating space for matrices
n = length(Phi);
Sf = zeros(n,1);
St = zeros(n,1);
Sm = zeros(n,1);
Efx = zeros(n,1);
x0 = zeros(n,1);
y0 = zeros(n,1);
A = zeros(3);

% Effective modulus for each fiber
for i = 1 : n
    
    % Tensor matrix (rotation)
    T = [(cosd(Phi(i)))^2               (sind(Phi(i)))^2                2*(sind(Phi(i)))*(cosd(Phi(i)));
         (sind(Phi(i)))^2               (cosd(Phi(i)))^2                -2*(sind(Phi(i)))*(cosd(Phi(i)));
         -(cosd(Phi(i)))*(sind(Phi(i))) (cosd(Phi(i)))*(sind(Phi(i)))   ((cosd(Phi(i)))^2)-((sind(Phi(i)))^2)];  
     
    % S in function of main axis
    L = eye(3).*[1, 1, 2];
    R = eye(3).*[1, 1, 0.5];
    Qbar = inv(T) * Q * L * T * R; %#ok<MINV>
    Sbar = inv(Qbar);
    
    % Modulus in principal coordinate (Each fiber)
%     Efx(i) = 1/Sbar(1,1);
     
    % Area for each fibre
    Sf(i) = (Data(i,4)/2) * (Data(i,5)/2) * pi; % Fiber surface area
    St(i) = Sf(i) / (Volf);                     % Total surface area
    Sm(i) = St(i) - Sf(i);                      % Matrix surface area
    
    % fiber location (Center of elipse)    
    x0(i) = Data(i,2);
    y0(i) = -Data(i,3);
    
    % Matrice A
    Af = Qbar * St(i);  % Each fiber
    A = A + Af;         % All fibers
end

% Total surface area
St_tot = sum(St);
Sf_tot = sum(Sf);
Sm_tot = sum(Sm);  
S = max(x0) * max(abs(y0)); % Area from picture

% Effective modulus for undetected zones

% Tensor matrix (rotation)
phi = 90;
T = [(cosd(phi))^2              (sind(phi))^2             2*(sind(phi))*(cosd(phi));
     (sind(phi))^2              (cosd(phi))^2             -2*(sind(phi))*(cosd(phi));
     -(cosd(phi))*(sind(phi))   (cosd(phi))*(sind(phi))   ((cosd(phi))^2)-((sind(phi))^2)];
     
Qm = inv(T) * Q * L * T * R;      %#ok<MINV>

Am = Qm * (S - St_tot);
A = A + Am;

Ap=inv(A);

% Modulus calculation
%Ex1 = (A(1,1)*A(2,2)-(A(1,2)^2))/(S*A(2,2));
Ex1 = 1/(S*Ap(1,1));
Ex1_GPa = Ex1/1e9

%% Color mapping angle of fibers

% Prealocation array space
z1 = zeros(n,3);
z2 = zeros(n,3);

% Colors division Phi
Phi_max = 90;
div = [0, 25, 60, 90];
colorPhi = [0 1 0; 1 1 1; 1 1 1];
% colorPhi = [0 1 0; 1 1 0; 1 0 0];

% Looping though all fiber location
for i = 1 : n  
    % Phi color
    z1(i,:) = colorPhi(find(div>Phi0(i),1)-1,:);
end

% Plot color map of phi and theta
figure(1)
scatter(x0,y0,10,z1,'filled')
axis([0 max(x0) min(y0) 0])
title('Phi');

GA = length(find(Phi0 < 25.5));
P = GA/n


