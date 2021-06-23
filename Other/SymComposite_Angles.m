%% initialisation
clear all %#ok<CLALL>
clc
%#ok<*NOPTS>

%% Inputs

% Data from Figi
file = 'D:\Université\Matrise\Article\DataOutput\symetrique\SymetricLong.xlsx';
Data = xlsread(file);
Volf = 0.59;    % Fiber volume fraction

% nettoyage de 25% de la médianne
Keep=median(Data(:,5))-0.25*median(Data(:,5))<Data(:,5) & Data(:,5)<median(Data(:,5))+0.25*median(Data(:,5));
Data=Data(Keep,:);

% parameters
x = Data(:,2);
y = Data(:,3);
b = Data(:,4);
a = Data(:,5);
RA = Data(:,7);

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
Gamma = 51;
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

%% delete angles larger than thresshold
TH = 65; %degres
rowsToDelete = Phi < -TH | Phi > TH;
Phi(rowsToDelete) = [];
Theta(rowsToDelete) = [];
x(rowsToDelete) = [];
y(rowsToDelete) = [];

%% proportion of fibre angles that are detected
n = length(Phi);

F90 = 0;
F45 = 0;
F0 = 0;
Fm90 = 0;
Fm45 = 0;
UF = 0;

for i = 1:n
    if Phi(i) > -54 && Phi(i) < -40
        Fm45 = Fm45 + 1;
    elseif Phi(i) > -5 && Phi(i) < 5
        F0 = F0 + 1;
    elseif Phi(i) > 40 && Phi(i) < 50
        F45 = F45 + 1;
    else
        UF = UF + 1;
        UFlocX(i,1) = x(i); %#ok<*SAGROW>
        UFlocY(i,1) = -y(i);
        UF_Phi(i,1) = Phi(i);
        z2(i,:)=color(find(div>Phi(i),1)-1,:);
    end
end 

rowsToDelete = UF_Phi == 0;
UF_Phi(rowsToDelete) = [];

disp('la quantité de fibre selon les angles, seulement les fibres detecter')
FibreCount = table(F45, F0, Fm45, UF, n)

P45 = F45/n;
P0 = F0/n;
Pm45 = Fm45/n;
PUF = UF/n;
Counted = (F45 + F0 + Fm45 + UF)/n;

disp('la proportion des angle selon les numbre total de fibre, seulement les fibres detecter')
FibreProportion = table(P45, P0, Pm45, PUF, Counted)

%% proportion of fibre angles that are not detected
 
% Prealocating space for matrices
n = length(Phi);
Sf = zeros(n,1);
St = zeros(n,1);
Sm = zeros(n,1);
x0 = zeros(n,1);
y0 = zeros(n,1);

% Effective modulus for each fiber
for i = 1 : n
    
    % Area for each fibre
    Sf(i) = (b(i)/2) * (a(i)/2) * pi; % Fiber surface area
    St(i) = Sf(i) / (Volf);                     % Total surface area
    Sm(i) = St(i) - Sf(i);                      % Matrix surface area
    
    % fiber location (Center of elipse)    
    x0(i) = x(i);
    y0(i) = -y(i);
 
    z1(i,:)=color(find(div>Phi(i),1)-1,:); 
end

% Total surface area
St_tot = sum(St);
Sf_tot = sum(Sf);
Sm_tot = sum(Sm);  
S = max(x0) * max(abs(y0)); % Area from picture

disp('Proportion de la surface qui a été detecter selon les airs des fibres')
PSurface = St_tot/S

%% Color mapping angle of fibers

figure(1)
scatter(x0,y0,10,z1,'filled')
colorbar
colormap(jet(n_div))
hcb=colorbar;
set(hcb,'YTick',div)
axis equal
axis([0 max(x0) round(min(y0)) 0])
caxis([min(div) max(div)])
title('Phi');

% plot undifined location
figure(2)
scatter(UFlocX,UFlocY,10,z2,'filled')
colorbar
colormap(jet(n_div))
hcb=colorbar;
set(hcb,'YTick',div)
axis equal
axis([0 max(x0) round(min(y0)) 0])
caxis([min(div) max(div)])
title('Undefined fiber orientation');

% histograme of fibre angles
figure(3)
h = histogram(Theta);

figure(4)
H = histogram(Phi);

% Plot histograme of Phi undefined
figure (5)
hist = histogram(UF_Phi);


