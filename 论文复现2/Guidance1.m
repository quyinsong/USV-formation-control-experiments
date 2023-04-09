function [nu_c,sigma_hato] = Guidance1( hi, ai0, zi_bar, nu_bar, vs, R_psi, theta, ri, p0, ts)
%GUIDANCE1 此处显示有关此函数的摘要
%   此处显示详细说明

persistent Ki emi S sigma_hat zi_bar_hat Ko1 Ko2
if isempty(Ki)
    Ki = diag([0.05 0.05]);
    emi = 0.01;
    S = [0 -1; 1 0];
    sigma_hat = [0 0]';
    zi_bar_hat = zi_bar;
    Ko1 = diag([5 5]);
    Ko2 = diag([15 15]);
end

% 未知估计
vartheta = hi*nu_bar-ri*S*zi_bar;
zi_bar_hat_dot = vartheta+sigma_hat-Ko1*(zi_bar_hat-zi_bar);
sigma_hat_dot = -Ko2*(zi_bar_hat-zi_bar);
zi_bar_hat = zi_bar_hat_dot*ts+zi_bar_hat;
sigma_hat = sigma_hat_dot*ts+sigma_hat;

% nu_c = hi\(-Ki*zi_bar/sqrt(norm(zi_bar)^2+emi^2)+ai0*vs*R_psi'*p0-sigma_hat);

nu_c = hi\(-Ki*zi_bar-sigma_hat);

sigma_hato = sigma_hat;


end

