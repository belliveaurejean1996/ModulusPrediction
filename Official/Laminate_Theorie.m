%% initialisation
clear all %#ok<CLALL>
clc

% error/warning delete
%#ok<*SAGROW>
%#ok<*MINV>
%#ok<*NOPTS>

%% Inputs

% Material properties
E1 = 135e9;     % Pa
E2 = 10e9;      % Pa
G12 = 5.2e9;    % Pa
Nu12 = 0.22;    % Estimation from composite book
Nu21 = (Nu12 * E2)/E1;

% bottom-up configuration
Phi = [-45, 45, 90, 0, -45, 45, 90, 0, -45, 45, 90, 0, -45, 45, 90, 0, 0, 90, 45, -45 ,0, 90, 45, -45, 0, 90, 45, -45, 0, 90, 45, -45].';
% Phi = [-50, 42, 90, 0, -45, 46, 90, 0, -47, 44, 90, 0, -52, 45, 90, 0, 0, 90, 43, -49 ,0, 90, 42, -50, 0, 90, 44, -46, 0, 90, 48, -43].';
thick = 0.125; %mm

% Phi = [45, -45, -45, 45];
% thick = 0.25; %mm


%% Z positions
n = length(Phi);

z(1) = 0;
for i = 2 : n+1
    z(i) = z(i-1) + thick; 
end

mid = max(z)/2;
for i = 1 : n+1
	Dis(i) = z(i) - mid;
end

%% laminate theorie

% Pli matrix Q
Q = [E1/(1-Nu12*Nu21)           (Nu21*E1)/(1-Nu12*Nu21) 0;
     (Nu12*E2)/(1-Nu12*Nu21)    E2/(1-Nu12*Nu21)        0;
     0                          0                       G12];
 
% Prealocate matrix space 
Sf = zeros(n,1);
St = zeros(n,1);
Sm = zeros(n,1);
Efx = zeros(n,1);
EFx = zeros(n,1);
A = zeros(3);
B = zeros(3);
D = zeros(3);

for i = 1 : n
    
    % Tensor matrix
    T = [(cosd(Phi(i)))^2               (sind(Phi(i)))^2                2*(sind(Phi(i)))*(cosd(Phi(i)));
         (sind(Phi(i)))^2               (cosd(Phi(i)))^2                -2*(sind(Phi(i)))*(cosd(Phi(i)));
         -(cosd(Phi(i)))*(sind(Phi(i))) (cosd(Phi(i)))*(sind(Phi(i)))   ((cosd(Phi(i)))^2)-((sind(Phi(i)))^2)];
     
    % Q in function of main axis
    L = eye(3); L(3,3) = 2;
    R = eye(3); R(3,3) = 0.5;
    Qbar = inv(T) * Q * L * T * R; 
    
    % Matrix A
    Ap = Qbar * (Dis(i+1) - Dis(i));
    A = A + Ap;
    
    % Matrix B
    Bp = 1/2 * Qbar * (Dis(i+1)^2 - Dis(i)^2);
    B = B + Bp;    
    
    % Matrix D
    Dp = 1/3 * Qbar * (Dis(i+1)^3 - Dis(i)^3);
    D = D + Dp;
    
end

% Complet stiffness matrix
M = [A, B; B, D];
Minv = inv(M);

% Complet thickness of composite
Thick = length(Phi) * thick; 

% Tensile Modulus (GPa)
Ex = (1/(Thick*Minv(1,1))) / 1e9  



