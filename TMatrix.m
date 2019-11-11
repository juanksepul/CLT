% T matrix
% Function to calculate Transformation matrix
% theta: rotation angle in radians

function T = TMatrix(theta)
t = theta;

T = [(cos(t))^2    , (sin(t))^2   , 2*sin(t)*cos(t)        ;
     (sin(t))^2    , (cos(t))^2   , -2*sin(t)*cos(t)       ;
     -sin(t)*cos(t), sin(t)*cos(t), (cos(t))^2 - (sin(t))^2];