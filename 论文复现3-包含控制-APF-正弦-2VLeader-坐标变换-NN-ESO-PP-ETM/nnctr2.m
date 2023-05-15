function [tau_c,e_bar,f_hat] = nnctr2( nu_tf, nu_c, tau, tauc, Gamma, q, vartheta,zeta, ts )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

u = nu_tf(1); v = nu_tf(2); r = nu_tf(3);
X = [u v r]';  % 神经网络输入

% 初始化
persistent alpha_u alpha_r Tu Tr T nu_cf W_u W_r c b
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
end

m11 = 25.8; m22 = 33.8;
M = diag([m11 Gamma/m22]);


% filter
nu_cf_dot = -T\(nu_cf-nu_c);
nu_cf = nu_cf_dot*ts+nu_cf;
ucf = nu_cf(1);
rcf = nu_cf(2);
ucf_dot = nu_cf_dot(1);
rcf_dot = nu_cf_dot(2);

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
su = u-ucf-alpha_u;
sr = r-rcf-alpha_r;

e_bar = [su,sr]';

tu1 = 5*su-ucf_dot+Tu*alpha_u+fu_hat+q(1)*vartheta(1);
tr1 = 5*sr-rcf_dot+Tr*alpha_r+fr_hat+q(2)*vartheta(2);

% tu1 = 5*su-ucf_dot+Tu*alpha_u+fu_hat;
% tr1 = 5*sr-rcf_dot+Tr*alpha_r+fr_hat;

% tuc = -(1+0.05)*(su*tu1^2*M(1,1))/sqrt(su^2*tu1^2+0.05);
% trc = -(1+0.05)*(sr*tr1^2*M(2,2))/sqrt(sr^2*tr1^2+0.05);

tuc = -tu1*M(1,1);
trc = -tr1*M(2,2);

tau_c = [tuc,trc]';


% 权值更新
W_u_dot = 2*(e_bar(1)*Phi_u-0.5*W_u);
W_u = W_u_dot*ts+W_u;
W_r_dot = 1*(e_bar(2)*Phi_r-0.5*W_r);
W_r = W_r_dot*ts+W_r;

end

