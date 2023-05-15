function pi_usv = APF_usv(USVp,usvd,Rusv_up,Rusv_dn)
%APF_O 此处显示有关此函数的摘要
%   此处显示详细说明

% 输出为梯度方向
n = length(usvd(1,:)); % 检测到障碍物的个数
pi_usv = [0 0]';
for i=1:n
    pi_usv = pi_usv+(4*(Rusv_up^2-Rusv_dn^2)*(norm(USVp-usvd(:,i))^2-Rusv_up^2))*...
             (USVp-usvd(:,i))/(norm(USVp-usvd(:,i))^2-Rusv_dn^2)^3;
end