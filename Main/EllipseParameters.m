function[a,b,x,y,Phi,Theta,n] = EllipseParameters(Data)

Gamma = 45;

% Néttoyage de 25% de la médianne
Keep=median(Data(:,5))-0.25*median(Data(:,5))<Data(:,5) & Data(:,5)<median(Data(:,5))+0.25*median(Data(:,5));
Data=Data(Keep,:);

% parameters
x = Data(:,2);
y = Data(:,3);
b = Data(:,4);
a = Data(:,5);
RA = Data(:,7);

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
% Rotation de l'axe x par 45 deg
R=[1    0               0;
   0    cosd(Gamma)    -sind(Gamma);
   0    sind(Gamma)    cosd(Gamma)];
w=R*v;

% Calcul de Phi
Phi=transpose(acosd(w(3,:)));
Phi0=Phi;
Phi=Phi.*directionindex;

% Calcul de theta
% Theta=transpose(asind(w(2,:)'/sind(beta)));
Theta=transpose(asind(w(2,:)));

% Delete angles larger than thresshold
TH = 65; %degres
rowsToDelete = Phi < -TH | Phi > TH;
Phi(rowsToDelete) = [];
Theta(rowsToDelete) = [];
x(rowsToDelete) = []; %#ok<*NASGU>
y(rowsToDelete) = [];
a(rowsToDelete) = [];
b(rowsToDelete) = [];

% ajustement de l'axe y
y = -y + max(y);

PT = 3e-3; % plate tickness is mm
FD = 6.2e-6; % Fibre diametre m
FD = (max(y)/PT) * FD; % fibre diametre in Pixels

y = y*cosd(Gamma);
a = ones(length(Phi))*FD;
b = FD./cosd(Phi);

% n
n = length(Phi);

end
