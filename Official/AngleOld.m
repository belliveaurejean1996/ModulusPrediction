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
file = 'D:\Université\Matrise\Article - Modulus Prediction\DataOutput\BR1_200x';
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
Phi=transpose(acosd(w(3,:)));
Phi=Phi.*directionindex;
figure(1)
histogram(Phi)

% % Calcul de theta
% Theta=transpose(acosd(w(2,:)));
% % Theta=real(asind(w(2,:)'./sind(beta)));
% 
% Delete angles larger than thresshold
TH = 65; %degres
rowsToDelete = Phi < -TH | Phi > TH;
Phi(rowsToDelete) = [];
Theta(rowsToDelete) = [];
x(rowsToDelete) = []; %#ok<*NASGU>
y(rowsToDelete) = [];

% ajustement de l'axe y
y = -y + max(y);

% %% Prealocating space for matrices
% n = length(Phi);
% Sf = zeros(n,1);
% St = zeros(n,1);
% Sm = zeros(n,1);
% c1 = zeros(n,3);
% A = zeros(3);
% B = zeros(3);
% D = zeros(3);
% 
% %% Modified laminate Theorie
% 
% % Matrix Q
% Q = [E1/(1-Nu12*Nu21)           (Nu21*E1)/(1-Nu12*Nu21) 0;
%      (Nu12*E2)/(1-Nu12*Nu21)    E2/(1-Nu12*Nu21)        0;
%      0                          0                       G12];
% 
% % Micrographie dimentions
% Width = max(x) - min(x);
% Higth = max(y) - min(y);
% Mid = Higth / 2;
% 
% % Effective modulus for each fiber
% for i = 1 : n
%     
%     % Tensor matrix (rotation)
%     T = [(cosd(Phi(i)))^2               (sind(Phi(i)))^2                2*(sind(Phi(i)))*(cosd(Phi(i)));
%          (sind(Phi(i)))^2               (cosd(Phi(i)))^2                -2*(sind(Phi(i)))*(cosd(Phi(i)));
%          -(cosd(Phi(i)))*(sind(Phi(i))) (cosd(Phi(i)))*(sind(Phi(i)))   ((cosd(Phi(i)))^2)-((sind(Phi(i)))^2)];  
%      
%     % Q bar
%     L = eye(3).*[1, 1, 2];
%     R = eye(3).*[1, 1, 0.5];
%     Qbar = inv(T) * Q * L * T * R; %#ok<MINV>
%      
%     % Area for each fibre
%     Sf(i) = (b(i)/2) * (a(i)/2) * pi; % Fiber surface area
%     St(i) = Sf(i) / (Volf);  % Total surface area
%     Sm(i) = St(i) - Sf(i);   % Matrix surface area
%     
%     % Equivalent ply distance in Y
%     Tickness = St(i) / Width;
%     Distance = Mid - y(i);
%     Z1 = Distance - (Tickness/2);
%     Z2 = Distance + (Tickness/2);
%     
%     % Matrice A
%     Af = Qbar * (Z2 - Z1);  
%     A = A + Af;
%     
%     % Matrice B
%     Bf = (1/2) * Qbar * (Z2^2 - Z1^2);
%     B = B + Bf;
%     
%     % Matrice D
%     Df = (1/3) * Qbar * (Z2^3 - Z1^3);
%     D = D + Df;
%     
%     % color mapping for Phi angle
%     c1(i,:)=color(find(div>Phi(i),1)-1,:);
% end
% 
% % Total surface area
% St_tot = sum(St);
% Sf_tot = sum(Sf);
% Sm_tot = sum(Sm);  
% S = max(x) * max(abs(y)); % Area from picture
% US = (S - St_tot);
% 
% % Effective modulus for undetected zones
% phi = 90;
% thickness = (US/Width)/5;
% 
% z1(1) = -Higth/2;
% z2(1) = -(Higth/2) + thickness;
% 
% z1(2) = (Mid - 3*Higth/4) - thickness/2;
% z2(2) = (Mid - 3*Higth/4) + thickness/2;
% 
% z1(3) = -thickness/2;
% z2(3) = thickness/2;
% 
% z1(4) = (Mid - Higth/4) - thickness/2;
% z2(4) = (Mid - Higth/4) + thickness/2;
% 
% z1(5) = Mid - thickness;
% z2(5) = 0;
% 
% T = [(cosd(phi))^2              (sind(phi))^2             2*(sind(phi))*(cosd(phi));
%      (sind(phi))^2              (cosd(phi))^2             -2*(sind(phi))*(cosd(phi));
%      -(cosd(phi))*(sind(phi))   (cosd(phi))*(sind(phi))   ((cosd(phi))^2)-((sind(phi))^2)];
%      
% Qm = inv(T) * Q * L * T * R;  %#ok<MINV>
% 
% for i = 1 : 5
%     Am = Qm * (z2(i) - z1(i));
%     A = A + Am;
% 
%     Bm = (1/2) * Qm * (z2(i)^2 - z1(i)^2);
%     B = B + Bm;
% 
%     Dm = (1/3) * Qm * (z2(i)^3 - z1(i)^3);
%     D = D + Dm;
% end
% 
% % Inverse de matrice de régiditer
% M = [A B; B D];
% Minv = inv(M);
% 
% % Modulus calculation
% disp('Predicted tensile modulus')
% Ex1 = 1/(Higth*Minv(1,1));
% Ex_GPa = Ex1/1e9
% 
% %% Orientation index by ply
% % PositivePly = [sum(Phi > 0); sum(Phi > 0)/n];
% % NegativePly = [sum(Phi < 0); sum(Phi < 0)/n];
% % Rows = {'Fibre count';'Proportion'};
% % Orientation = table(PositivePly, NegativePly,'RowNames',Rows)
% 
% %% Color mapping angle of fibers
% % % Plot Fibre orientation of Phi
% % figure(2)
% % scatter(x,y,5,c1,'filled')
% % colorbar
% % colormap(jet(n_div))
% % hcb=colorbar;
% % set(hcb,'YTick',div)
% % axis([0 max(x) 0 round(max(y))])
% % caxis([min(div) max(div)])
% % axis equal
% % title('Phi');
% % yline(Mid);
% 
% % Plot histogram of fibre angles
% % figure(3)
% % h = histogram(Theta);
% % title('Theta')
% % xlabel('Fibre angles')
% % ylabel('Quantity of fibres detected')
% %  
% % figure(4)
% % H = histogram(Phi,100);
% % title('Phi')
% % xlabel('Fibre angles')
% % ylabel('Quantity of fibres detected')
%  
% % figure(10)
% % h = histogram(y);
