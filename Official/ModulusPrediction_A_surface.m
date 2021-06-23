%% initialisation
clear all %#ok<CLALL>
clc
%#ok<*NOPTS>

%% Inputs
% Data from Figi
% Symetric composite
file = 'D:\Université\Matrise\Article - Modulus Prediction\DataOutput\symetrique\SymetricLong.xlsx';
Gamma = 51;
E1 = 130e9;% Pa

% DLF composite
% file = 'D:\Université\Matrise\Article - Modulus Prediction\DataOutput\New echantillon\L4-200x';
% Gamma = 45; 
% E1 = 135e9;% Pa

Data = xlsread(file);

% nettoyage de 25% de la médianne
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

% %colors for phi
div_max=90;
div=-div_max:(5):div_max;
n_div=length(div)-1;
color=jet(n_div);

%% Parametres calculation
% calcul de alpha
alpha=wrapTo360(RA-90);
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
beta=acosd(a./b);

% vecteur unitaire dans le système d'axe x'y'z'
v=transpose([sind(beta).*sind(alpha),sind(beta).*cosd(alpha),cosd(beta)]);

% transformation au système d'axe xyz
% rotation de l'axe x par -45 deg
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

% delete angles larger than thresshold
TH = 65; %degres
rowsToDelete = Phi < -TH | Phi > TH;
Phi(rowsToDelete) = [];
Theta(rowsToDelete) = [];
x(rowsToDelete) = [];
y(rowsToDelete) = [];

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
Am45 = 0;
A0 = 0;
A45 = 0;

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
     
    % Area for each fibre
    Sf(i) = (b(i)/2) * (a(i)/2) * pi; % Fiber surface area
    St(i) = Sf(i) / (Volf);  % Total surface area
    Sm(i) = St(i) - Sf(i);   % Matrix surface area
    
    % fiber location (Center of elipse)    
    x0(i) = x(i);
    y0(i) = -y(i);
    
    % Matrice A
    Af = Qbar * St(i);  % Each fiber
    A = A + Af;         % All fibers
    
    % color mapping for Phi angle
    c1(i,:)=color(find(div>Phi(i),1)-1,:);   %#ok<SAGROW>
end

% Total surface area
St_tot = sum(St);
Sf_tot = sum(Sf);
Sm_tot = sum(Sm);  
S = max(x0) * max(abs(y0)); % Area from picture
US = (S - St_tot);

% Effective modulus for undetected zones
% Tensor matrix (rotation)
phi = 90;
T = [(cosd(phi))^2              (sind(phi))^2             2*(sind(phi))*(cosd(phi));
     (sind(phi))^2              (cosd(phi))^2             -2*(sind(phi))*(cosd(phi));
     -(cosd(phi))*(sind(phi))   (cosd(phi))*(sind(phi))   ((cosd(phi))^2)-((sind(phi))^2)];
     
Qm = inv(T) * Q * L * T * R;      %#ok<MINV>

Am = Qm * US;
A = A + Am;

% Ajustment for width
Width = max(x) - min(x);
Higth = max(y) - min(y);

% inverse of matrix A
Ap=inv(A);

% Modulus calculation
disp('Predicted tensil modulus')
Ex1 = 1/(S*Ap(1,1));
Ex_GPa = Ex1/1e9

%% Color mapping angle of fibers
% Plot Fibre orientation of Phi
% figure(1)
% scatter(x0,y0,10,c1,'filled')
% colorbar
% colormap(jet(n_div))
% hcb=colorbar;
% set(hcb,'YTick',div)
% axis equal
% axis([0 max(x0) round(min(y0)) 0])
% caxis([min(div) max(div)])
% title('Phi');
% 
% % Plot histogram of fibre angles
% figure(3)
% h = histogram(Theta);
% title('Theta')
% xlabel('Angle')
% ylabel('Number of fibres detected')
% 
% figure(4)
% H = histogram(Phi,100);
% title('Phi')
% xlabel('Angle')
% ylabel('Number of fibres detected')
