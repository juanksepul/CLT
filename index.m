% Classic laminate theory - Orthotropic material
% This script calculates the mechanical properties for a laminate of
% orthotropic plies

% Load mechanical properties data structure
PlyProperties

E1 = ply.E1;
E2 = ply.E2;
mu12 = ply.mu12;
G12 = ply.G12;

mu21 = mu12*(E2/E1);

% Assembly flexibility matrix
S = [ 1/E1    , -mu12/E1, 0    ;
      -mu12/E1, 1/E2    , 0    ;
      0       , 0       , 1/G12];

theta = 0 : 1 : 90;
E1vector = zeros(length(theta), 1);
E2vector = zeros(length(theta), 1);
mu12vector = zeros(length(theta), 1);
G12vector = zeros(length(theta), 1);

for i = 1 : length(theta)
    T = TMatrix(deg2rad(theta(i)));
    Sbar = T'*S*T;
    E1vector(i) = 1/Sbar(1, 1);
    E2vector(i) = 1/Sbar(2, 2);
    mu12vector(i) = - Sbar(1, 2) * E1vector(i);
    G12vector (i) = 1/Sbar(3,3);
end

plot(theta, E1vector)
hold on
plot(theta, E2vector)
plot(theta, mu12vector)
plot(theta, G12vector)
legend('E_1', 'E_2', '\mu_12', 'G_12', 'Location', 'north')
xlabel ('\theta (�)')

% T = TMatrix(theta);
% Sbar = T'*S*T;

% Laminate stacking sequence
LSS = [-45, 45, 0, 90, 45, -45, 0, 0, 90, 90];


