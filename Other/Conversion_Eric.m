%% initialisation
% clear all
clc

%% Inputs

% Data from Figi
% Data = xlsread('P4_1T.xlsx');

% Material properties
E1 = 135e9; % Pa
E2 = 10e9;  % Pa
G12 = 5.2e9;% Pa
VolM = 0.34; % Matrix volume fraction

Nu12 = E1/(2*G12) - 1;
Nu21 = (Nu12 * E2)/E1;

%% Figi operation

% nettoyage de 25% de la médianne
Keep=median(Data(:,5))-0.25*median(Data(:,5))<Data(:,5) & Data(:,5)<median(Data(:,5))+0.25*median(Data(:,5));
Data=Data(Keep,:);

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

%% Modified laminate

% Pli matrix Q
Q = [E1/(1-Nu12*Nu21)           (Nu21*E1)/(1-Nu12*Nu21) 0;
     (Nu21*E1)/(1-Nu12*Nu21)    E2/(1-Nu12*Nu21)        0;
     0                          0                       G12];

% Rotate fiber properties to major axis
n = length(Phi);
A11 = zeros(n,1);
Sf = zeros(n,1);
St = zeros(n,1);
Sm = zeros(n,1);

% for loop to go trough all fibers
for i = 1 : n
    
    % Tensor matrix
    T = [(cosd(Phi(i)))^2               (sind(Phi(i)))^2                2*(sind(Phi(i)))*(cosd(Phi(i)));
         (sind(Phi(i)))^2               (cosd(Phi(i)))^2                -2*(sind(Phi(i)))*(cosd(Phi(i)));
         -(cosd(Phi(i)))*(sind(Phi(i))) (cosd(Phi(i)))*(sind(Phi(i)))   (cosd(Phi(i)))^2-(sind(Phi(i)))^2];
    
    % Q in function of main axis
    Qbar = inv(T) * Q * T;  %#ok<MINV>
     
    % Area for each fibre
    Sf(i) = Data(i,4) * Data(i,5) * pi;
    St(i) = Sf(i) / (1-VolM);
    Sm(i) = St(i) - Sf(i);
    
    % Rigidity matrix A
    A = Qbar * St(i);
    A11(i) = A(1,1);
end

Atot = sum(A11) %#ok<NOPTS>
Sftot = sum(Sf);
Sttot = sum(St);
Smtot = sum(Sm);

Ex = 1/(Sttot*Atot)

