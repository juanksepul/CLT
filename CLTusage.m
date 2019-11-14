% CLT Usage
% Example of usage of function CLT

% PlyProperties = [E1, E2, mu12, G12, th]
PlyProperties = [90, 10, 0.3, 5, 0.1];

% LSS = [theta1, theta2, theta3... thetaN]
LSS = [-45, 45, 0, 90, 45, -45, 0, 0, 90, 90];

% Loads and moments [Nx, Ny, Nxy, Mx, My, Mxy] ([0, 0, 0, 0, 0, 0] when is the object of calculation)
Loads = [1000, 2000, 0, 10, 0, 0];

% Deformations [ep_x, ep_y, gam_xy, kap_x, kap_y, kap_xy] ([0, 0, 0, 0, 0, 0] when is the object of calculation)
% Deformations = [-0.0032e-2, 0.0282e-2, 0.0069e-2, 4.1098, -1.0368, 0.5765];
Deformations = [0, 0, 0, 0, 0, 0];

[A, B, D, Loads, Deformations] = CLT(PlyProperties, LSS, Loads, Deformations)