function [nu_c,u_hat,v_hat] = Guidance22( hi, R_psi, zi1, zi2, pj_dot, ri, psi, ts)
%GUIDANCE1 此处显示有关此函数的摘要
%   此处显示详细说明

persistent Ki p1 p2 k1 k2
if isempty(Ki)
    Ki = diag([0.1 0.05]);
    p1 = 0;
    p2 = 0;
    k1 = 25;
    k2 = 25;
end

% 速度估计
theta1 = -pj_dot(1)*cos(psi)-pj_dot(2)*sin(psi)+ri*zi1(2);
theta2 = pj_dot(1)*sin(psi)-pj_dot(2)*cos(psi)-ri*zi1(1);
u_hat = p1+k1*zi1(1);
p1_dot = -k1*p1-k1^2*zi1(1)-k1*theta1;
v_hat = p2+k2*zi1(2);
p2_dot = -k1*p2-k2^2*zi1(2)-k2*theta2;

p1 = p1_dot*ts+p1;
p2 = p2_dot*ts+p2;

nu_c = hi\(-Ki*zi2+R_psi'*pj_dot-[0 v_hat]');
if nu_c(1)<=0.01
    nu_c(1)=0.01;
end

end

