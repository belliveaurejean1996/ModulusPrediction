function[Ex1_ABD] = ModifiedLaminate_ABD(E1,E2,G12,Volf,Nu12,Nu21,Phi,x,y,a,b)

% Prealocating space for matrices
n = length(Phi);
Sf = zeros(n,1);
St = zeros(n,1);
Sm = zeros(n,1);
c1 = zeros(n,3);
A = zeros(3);
B = zeros(3);
D = zeros(3);


% Matrix Q
Q = [E1/(1-Nu12*Nu21)           (Nu21*E1)/(1-Nu12*Nu21) 0;
     (Nu12*E2)/(1-Nu12*Nu21)    E2/(1-Nu12*Nu21)        0;
     0                          0                       G12];

% Micrographie dimentions
Width = max(x) - min(x);
Higth = max(y) - min(y);
Mid = Higth / 2;

% Effective modulus for each fiber
for i = 1 : n
    
    % Tensor matrix (rotation)
    T = [(cosd(Phi(i)))^2               (sind(Phi(i)))^2                2*(sind(Phi(i)))*(cosd(Phi(i)));
         (sind(Phi(i)))^2               (cosd(Phi(i)))^2                -2*(sind(Phi(i)))*(cosd(Phi(i)));
         -(cosd(Phi(i)))*(sind(Phi(i))) (cosd(Phi(i)))*(sind(Phi(i)))   ((cosd(Phi(i)))^2)-((sind(Phi(i)))^2)];  
     
    % Q bar
    L = eye(3).*[1, 1, 2];
    R = eye(3).*[1, 1, 0.5];
    Qbar = inv(T) * Q * L * T * R; %#ok<MINV>
     
    % Area for each fibre
    Sf(i) = (b(i)/2) * (a(i)/2) * pi; % Fiber surface area
    St(i) = Sf(i) / (Volf);  % Total surface area
    Sm(i) = St(i) - Sf(i);   % Matrix surface area
    
    % Equivalent ply distance in Y
    Tickness = St(i) / Width;
    Distance = Mid - y(i);
    Z1 = Distance - (Tickness/2);
    Z2 = Distance + (Tickness/2);
    
    % Matrice A
    Af = Qbar * (Z2 - Z1);  
    A = A + Af;
    
    % Matrice B
    Bf = (1/2) * Qbar * (Z2^2 - Z1^2);
    B = B + Bf;
    
    % Matrice D
    Df = (1/3) * Qbar * (Z2^3 - Z1^3);
    D = D + Df;
    
%     % color mapping for Phi angle
%     c1(i,:)=color(find(div>Phi(i),1)-1,:);
end

% Total surface area
St_tot = sum(St);
Sf_tot = sum(Sf);
Sm_tot = sum(Sm);  
S = max(x) * max(abs(y)); % Area from picture
US = (S - St_tot);

% Effective modulus for undetected zones
phi = 90;
thickness = (US/Width)/5;

z1(1) = -Higth/2;
z2(1) = -(Higth/2) + thickness;

z1(2) = (Mid - 3*Higth/4) - thickness/2;
z2(2) = (Mid - 3*Higth/4) + thickness/2;

z1(3) = -thickness/2;
z2(3) = thickness/2;

z1(4) = (Mid - Higth/4) - thickness/2;
z2(4) = (Mid - Higth/4) + thickness/2;

z1(5) = Mid - thickness;
z2(5) = 0;

T = [(cosd(phi))^2              (sind(phi))^2             2*(sind(phi))*(cosd(phi));
     (sind(phi))^2              (cosd(phi))^2             -2*(sind(phi))*(cosd(phi));
     -(cosd(phi))*(sind(phi))   (cosd(phi))*(sind(phi))   ((cosd(phi))^2)-((sind(phi))^2)];
     
Qm = inv(T) * Q * L * T * R;  %#ok<MINV>

for i = 1 : 5
    Am = Qm * (z2(i) - z1(i));
    A = A + Am;

    Bm = (1/2) * Qm * (z2(i)^2 - z1(i)^2);
    B = B + Bm;

    Dm = (1/3) * Qm * (z2(i)^3 - z1(i)^3);
    D = D + Dm;
end

% Inverse de matrice de régiditer
M = [A B; B D];
Minv = inv(M);

% Modulus calculation
Ex1_ABD = (1/(Higth*Minv(1,1)))*1e-9;

end