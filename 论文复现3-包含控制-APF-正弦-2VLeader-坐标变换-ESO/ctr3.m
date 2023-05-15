function [tau_c,e_bar] = ctr3( nu_tf_hat, zeta_hat, nu_c, tau, tauc, Gamma, ts )
%CTR 此处显示有关此函数的摘要
%   此处显示详细说明
persistent  T nu_cf K alpha_u alpha_r Tu Tr
if isempty(T)
    T = diag([0.1 0.1]);
    nu_cf = nu_c;
    K = diag([2 2]);
    alpha_u = 0;
    alpha_r = 0;
    Tu = 1;
    Tr = 1;
end

m11 = 25.8;  m22 = 33.8;
M = diag([m11 Gamma/m22]);

% filter
nu_cf_dot = -T\(nu_cf-nu_c);
nu_cf = nu_cf_dot*ts+nu_cf;

delta_u = tauc(1)-tau(1);
delta_r = tauc(2)-tau(2);
alpha_u_dot = -Tu*alpha_u-delta_u/m11;
alpha_r_dot = -Tr*alpha_r-delta_r/M(2,2);
alpha_u = ts*alpha_u_dot+alpha_u;
alpha_r = ts*alpha_r_dot+alpha_r;

nu_tf_hat_bar = [nu_tf_hat(1),nu_tf_hat(3)]';
e_bar = nu_tf_hat_bar-nu_cf-[alpha_u,alpha_r]';
nu_cf_dot = [tanh(nu_cf_dot(1)/5)*5,tanh(nu_cf_dot(2)/5)*5]';
tau_c = M*(-zeta_hat-K*e_bar+nu_cf_dot-[Tu*alpha_u,Tr*alpha_r]');



end

