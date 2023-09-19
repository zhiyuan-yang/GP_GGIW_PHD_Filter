function [rmm_rotation, rmm_l, fai_x,fai_y,fai_z] = ...
    get_random_matrix_ellipse1(rmm_extent)
% notations and formmulas are from "Ellipse fitting based approach for extended object tracking"
% Equation (9)
rho_x = abs((rmm_extent(2,2) - rmm_extent(3,3))/(2*rmm_extent(2,3)));
rho_y = abs((rmm_extent(1,1) - rmm_extent(3,3))/(2*rmm_extent(1,3)));
rho_z = abs((rmm_extent(1,1) - rmm_extent(2,2))/(2*rmm_extent(1,2)));

fai_x = atan(-rho_x + sqrt(1 + rho_x^2));
fai_y = atan(-rho_y + sqrt(1 + rho_y^2));
fai_z = atan(-rho_z + sqrt(1 + rho_z^2));

rmm_rotation_x = [1,0,0;0,cos(fai_x),-sin(fai_x);0,sin(fai_x),cos(fai_x)];
rmm_rotation_y = [cos(fai_y),0,sin(fai_y);0,1,0;-sin(fai_y),0,cos(fai_y)];
rmm_rotation_z = [cos(fai_z),-sin(fai_z),0;sin(fai_z),cos(fai_z),0;0,0,1];
rmm_rotation = rmm_rotation_x*rmm_rotation_y*rmm_rotation_z;
% rmm_rotation = [cos(fai_z), -sin(fai_z); sin(fai_z), cos(fai_z)];
rmm_l_sq = diag(rmm_rotation'*rmm_extent*rmm_rotation);
rmm_l = rmm_l_sq.^(1/2);
end