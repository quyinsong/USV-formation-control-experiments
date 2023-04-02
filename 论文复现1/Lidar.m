function [ rho, lambda ] = Lidar( rx, ry )
%LIDAR [ rho, lambda ] = Lidar( rx, ry )
%   obtain LOS range and angle
rho = sqrt(rx^2+ry^2);
lambda = atan2(rx,ry);

end

