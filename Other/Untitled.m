%% initialisation
clear all %#ok<CLALL>
clc

%% Inputs

% Material properties
E1 = 9e9;     % Pa
E2 = 9e9;      % Pa
G12 = 5e9;    % Pa
Volf = 0.59;    % Matrix volume fraction

Nu12 = 0.2;
Nu21 = (Nu12 * E2)/E1;

Phi = [0, 45, -45].';

%% Modified laminate

% Pli matrix Q
Q = [E1/(1-Nu12*Nu21)           (Nu21*E1)/(1-Nu12*Nu21) 0;
     (Nu12*E2)/(1-Nu12*Nu21)    E2/(1-Nu12*Nu21)        0;
     0                          0                       G12];
 
% Rotate fiber properties to major axis
n = length(Phi);
Sf = zeros(n,1);
St = zeros(n,1);
Sm = zeros(n,1);
Efx = zeros(n,1);
EFx = zeros(n,1);

% for loop to go trough all fibers
for i = 1 : n
    
    % Tensor matrix
    T = [(cosd(Phi(i)))^2               (sind(Phi(i)))^2                2*(sind(Phi(i)))*(cosd(Phi(i)));
         (sind(Phi(i)))^2               (cosd(Phi(i)))^2                -2*(sind(Phi(i)))*(cosd(Phi(i)));
         -(cosd(Phi(i)))*(sind(Phi(i))) (cosd(Phi(i)))*(sind(Phi(i)))   ((cosd(Phi(i)))^2)-((sind(Phi(i)))^2)];
     
    % Q in function of main axis
    L = eye(3); L(3,3) = 2;
    R = eye(3); R(3,3) = 0.5;
    Qbar = inv(T) * Q * L * T * R; %#ok<MINV>
    
    Sbar = inv(Qbar);
    EFx(i) = 1/Sbar(1,1);
end


