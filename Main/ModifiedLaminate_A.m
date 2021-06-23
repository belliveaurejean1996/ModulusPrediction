function[Ex1_A] = ModifiedLaminate_A(E1,E2,G12,Volf,Nu12,Nu21,Phi,x,y,a,b)

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
distance = Mid - (Mid/2);
tickness = US/Width;
z1 = distance - (tickness/2);
z2 = distance + (tickness/2);

T = [(cosd(phi))^2              (sind(phi))^2             2*(sind(phi))*(cosd(phi));
     (sind(phi))^2              (cosd(phi))^2             -2*(sind(phi))*(cosd(phi));
     -(cosd(phi))*(sind(phi))   (cosd(phi))*(sind(phi))   ((cosd(phi))^2)-((sind(phi))^2)];
     
Qm = inv(T) * Q * L * T * R; %#ok<MINV>

Am = Qm * (z2 - z1);
A = A + Am;

% inverse of matrix A
Ainv=inv(A);

% Modulus calculation
Ex1_A = (1/(Higth*Ainv(1,1)))*1e-9;

end