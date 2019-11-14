% CLT Usage
% Example of usage of function CLT

% PlyProperties = [E1, E2, mu12, G12, th]
PlyProperties = [90, 10, 0.3, 5, 0.1];

% LSS = [theta1, theta2, theta3... thetaN]
LSS = [-45, 45, 0, 90, 45, -45, 0, 0, 90, 90];

[A, B, D, Loads, Deformations] = CLT(PlyProperties, LSS)