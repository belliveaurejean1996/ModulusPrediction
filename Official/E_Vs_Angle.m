clear all %#ok<CLALL>
clc

%% Material properties
E1 = 135e9;     % Pa
E2 = 10e9;      % Pa
G12 = 5.2e9;    % Pa
Volf = 0.59;    % Fiber volume fraction

Nu12 = 0.22; % Estimation from composite book
Nu21 = (Nu12 * E2)/E1;

% Pli angle
Theta = 0 : 0.01 : 90;

% Pli matrix Q (Doesn't change)
Q = [E1/(1-Nu12*Nu21)           (Nu21*E1)/(1-Nu12*Nu21) 0;
     (Nu12*E2)/(1-Nu12*Nu21)    E2/(1-Nu12*Nu21)        0;
     0                          0                       G12];

%% Loop through all angles
for i = 1 : length(Theta)
       
    % Tensor matrix (rotation)
    T = [(cosd(Theta(i)))^2                 (sind(Theta(i)))^2                  2*(sind(Theta(i)))*(cosd(Theta(i)));
         (sind(Theta(i)))^2                 (cosd(Theta(i)))^2                  -2*(sind(Theta(i)))*(cosd(Theta(i)));
         -(cosd(Theta(i)))*(sind(Theta(i))) (cosd(Theta(i)))*(sind(Theta(i)))   ((cosd(Theta(i)))^2)-((sind(Theta(i)))^2)];  
     
    % S in function of main axis
    L = eye(3).*[1, 1, 2];
    R = eye(3).*[1, 1, 0.5];
    Qbar = inv(T) * Q * L * T * R; %#ok<MINV>
    Sbar = inv(Qbar);
    
    % Modulus in principal coordinate
    Ex(i) = 1/Sbar(1,1);        %#ok<SAGROW>
    Ex_GPa(i) = Ex(i) / 1e9;    %#ok<SAGROW>
end

%% plot graph
figure(1)

plot(Theta, Ex_GPa, 'k')
axis([0 90 0 140])
xlabel('Ply angle (degrees)','FontSize',16)
ylabel('Effective modulus (GPa)','FontSize',16)
% 
% % Vertical line position
% N = E1*0.10;
% Pos = Theta(find(Ex > N, 1, 'last' ));
% line = xline(Pos,'--r');
% 
% % write excel data
% Output = [Ex_GPa; Theta]';
% filename = 'ModuleVSAngle.xlsx';
% writematrix(Output,filename)



