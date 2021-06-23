clear all
clc

D = 0.1; % diametre de la fibre
Th = 45; % angle de la fibre

u = linspace(0,2*pi,50);
v = linspace(0,2*pi,50);
[u,v] = meshgrid(u,v);

% x = v;
% y = D*cos(u);
% z = D*sin(u);

x = cosd(Th)*v - sind(Th)*cos(u);
y = sind(Th)*v + cosd(Th)*cos(u);
z = sin(u);

surf(x,y,z)
xlabel('x')
ylabel('y')
zlabel('z')

