function nu_c = Guidance3( hi, R_psi, zi, pj_dot, nu_bar, ts)
%GUIDANCE1 此处显示有关此函数的摘要
%   此处显示详细说明

persistent Ki S sigma_hat zi_hat Ko1 Ko2
if isempty(Ki)
    Ki = diag([0.05 0.05]);
    S = [0 -1; 1 0];
    sigma_hat = [0 0]';
    zi_hat = zi;
    Ko1 = diag([5 5]);
    Ko2 = diag([15 15]);
end

% 未知估计
ri = nu_bar(2);
vartheta = hi*nu_bar-ri*S*zi-R_psi'*pj_dot;
zi_hat_dot = vartheta+sigma_hat-Ko1*(zi_hat-zi);
sigma_hat_dot = -Ko2*(zi_hat-zi);
zi_hat = zi_hat_dot*ts+zi_hat;
sigma_hat = sigma_hat_dot*ts+sigma_hat;

nu_c = hi\(-Ki*zi+R_psi'*pj_dot-sigma_hat);
if nu_c(1)<=0.01
    nu_c(1)=0.01;
end

end

