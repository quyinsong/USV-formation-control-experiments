function [tau_c,e_bar,f_hat] = nnctr1( nu_tf, nu_c, tau, tauc, Gamma, ts )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

u = nu_tf(1); v = nu_tf(2); r = nu_tf(3);
X = [u v r]';  % 神经网络输入
nu_tf_bar = [u,r]';

% 初始化
persistent alpha_u alpha_r Tu Tr T nu_cf W_u W_r c b K
if isempty(alpha_u)
    alpha_u = 0;
    alpha_r = 0;
    Tu = 1;
    Tr = 1;
    T = 0.1;
    nu_cf = nu_c;
    c = 0.1*ones(3,11);
    b = 1*ones(1,11);
    W_u = zeros(11,1);
    W_r = zeros(11,1);
    K = diag([5,5]);
end

m11 = 25.8; m22 = 33.8;
M = diag([m11 Gamma/m22]);


% filter
nu_cf_dot = -T\(nu_cf-nu_c);
nu_cf = nu_cf_dot*ts+nu_cf;

% 辅助系统
delta_u = tauc(1)-tau(1);
delta_r = tauc(2)-tau(2);
alpha_u_dot = -Tu*alpha_u-delta_u/m11;
alpha_r_dot = -Tr*alpha_r-delta_r/M(2,2);
alpha_u = ts*alpha_u_dot+alpha_u;
alpha_r = ts*alpha_r_dot+alpha_r;

% 神经网络
Phi_u = zeros(11,1);
Phi_r = zeros(11,1);
for i = 1:11
    Phi_u(i,1) = exp(-norm(X-c(:,i))^2/(2*b(i)^2));
    Phi_r(i,1) = exp(-norm(X-c(:,i))^2/(2*b(i)^2));
end
fu_hat = W_u'*Phi_u;
fr_hat = W_r'*Phi_r;

f_hat = [fu_hat,fr_hat]';

% 控制率
e_bar = nu_tf_bar-nu_cf-[alpha_u,alpha_r]';
nu_cf_dot = [tanh(nu_cf_dot(1)/5)*5,tanh(nu_cf_dot(2)/5)*5]';
tau_c = M*(-f_hat-K*e_bar+nu_cf_dot-[Tu*alpha_u,Tr*alpha_r]');

% 权值更新
W_u_dot = 2*(e_bar(1)*Phi_u-0.5*W_u);
W_u = W_u_dot*ts+W_u;
W_r_dot = 1*(e_bar(2)*Phi_r-0.5*W_r);
W_r = W_r_dot*ts+W_r;

end

