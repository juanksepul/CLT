% Classic laminate theory - Orthotropic material
% This script calculates the mechanical properties for a laminate of
% orthotropic plies

% Load mechanical properties data structure
PlyProperties

E1   = ply.E1*10^9;
E2   = ply.E2*10^9;
mu12 = ply.mu12;
G12  = ply.G12*10^9;
th   = ply.thickness*10^-3;

mu21 = mu12*(E2/E1);

% Assembly flexibility matrix
S = [ 1/E1    , -mu12/E1, 0    ;
      -mu12/E1, 1/E2    , 0    ;
      0       , 0       , 1/G12];
  
Q = inv(S);

theta = 0 : 2 : 90;

E1vector   = zeros(length(theta), 1);
E2vector   = zeros(length(theta), 1);
mu12vector = zeros(length(theta), 1);
G12vector  = zeros(length(theta), 1);

for i = 1 : length(theta)
    T = TMatrix(deg2rad(theta(i)));
    Sbar = T'*S*T;
    E1vector(i) = (1/Sbar(1, 1))*10^-9;
    E2vector(i) = (1/Sbar(2, 2))*10^-9;
    mu12vector(i) = - Sbar(1, 2)*E1vector(i)*10^9;
    G12vector (i) = (1/Sbar(3,3))*10^-9;
end

figure(1)
% Variation of ply properties with theta
plot(theta, E1vector)
hold on
plot(theta, E2vector)
plot(theta, G12vector)
legend ('E_1', 'E_2', 'G_12', 'Location', 'north')
xlabel ('\theta (°)')
ylabel ('GPa')

figure(2)
plot(theta, mu12vector)
legend ('\mu_{12}')
xlabel ('\theta (°)')

% T = TMatrix(theta);
% Sbar = T'*S*T;

% Laminate stacking sequence
LSS = [-45, 45, 0, 90, 45, -45, 0, 0, 90, 90]; % EVEN NUMBER!
nPlies = length(LSS);

% Rotate laminae properties
QbarVector = zeros(3, 3);
R = [1, 0, 0;
     0, 1, 0;
     0, 0, 2];
Rinv = inv(R);
for i = 1 : nPlies
    T = TMatrix(deg2rad(LSS(i)));
    Tinv = inv(T);
    QbarVector(:, :, i) = Tinv*Q*R*T*Rinv;
end

% Laminate thickness
Th = th*nPlies;

z = zeros(1, nPlies);

z0 = -Th/2;
% Z cordinates for upper side (negative)
for i = 1 : nPlies/2
    z(1, i) = z0 + i*th;
end

zN = Th/2;
% Z coordinates for lower side (positive)
for i = 1 : nPlies/2
    z(1, nPlies/2 + i) = zN - (nPlies/2 - i)*th;
end

z = cat(2, z0, z);

A = zeros(3, 3);
B = zeros(3, 3);
D = zeros(3, 3);
% Calculate A, B, D matrices



for i = 1 : nPlies
    A = A + QbarVector(:, :, i).*(z(i + 1)    - z(i)   );
    B = B + QbarVector(:, :, i).*(z(i + 1).^2 - z(i).^2);
    D = D + QbarVector(:, :, i).*(z(i + 1).^3 - z(i).^3);
end


B = 0.5*B;
D = (1/3)*D;
