function [tau_c,fhat] = ctr2( nu_bar, nu_c, tau, tauc, ts )
%CTR 此处显示有关此函数的摘要
%   此处显示详细说明

m11 = 25.8; m33 = 2.76;
M = diag([m11 m33]);

%--------------------------ESO-----------------------------
% persistent Ko1 Ko2 nu_bar_hat f_hat
% if isempty(f_hat)
%     Ko1 = diag([5 5]);
%     Ko2 = diag([10 10]);
%     f_hat = [0 0]';
%     nu_bar_hat = nu_bar;
% end
% nu_bar_hat_dot = f_hat-Ko1*(nu_bar_hat-nu_bar)+M\tau;
% f_hat_dot = -Ko2*(nu_bar_hat-nu_bar);
% nu_bar_hat = nu_bar_hat_dot*ts+nu_bar_hat;
% f_hat = f_hat_dot*ts+f_hat;

%-----------------------------DO-----------------------
persistent Ko varrho
if isempty(varrho)
    Ko = diag([1 1]);
    varrho = [0 0]';
end
f_hat = varrho+Ko*nu_bar;
varrho_dot = -Ko*varrho-Ko*(M\tau+Ko*nu_bar);
varrho = ts*varrho_dot+varrho;

%-------------------1st filter-----------------------
persistent T nu_cf
if isempty(nu_cf)
    T = diag([0.1 0.1]);
    nu_cf = nu_c;
end
nu_cf_dot = -T\(nu_cf-nu_c);
nu_cf = nu_cf_dot*ts+nu_cf;

% --------------------Auxiliary system 1-----------------
% persistent  alpha_u alpha_r gamma_u gamma_r Tu Tr
% if isempty(alpha_u)
%     alpha_u = 0;
%     alpha_r = 0;
%     gamma_u = 1;
%     gamma_r = 1;
%     Tu = 0.1;
%     Tr = 0.1;
% end
% delta_u = tauc(1)-tau(1);
% delta_r = tauc(2)-tau(2);
% alpha_u_dot = cosh(alpha_u)^2*(-Tu*alpha_u-delta_u/m11)/gamma_u;
% alpha_r_dot = cosh(alpha_r)^2*(-Tr*alpha_r-delta_r/m33)/gamma_r;
% alpha_u = ts*alpha_u_dot+alpha_u;
% alpha_r = ts*alpha_r_dot+alpha_r;
% 
% e_bar = nu_bar-nu_cf-[gamma_u*tanh(alpha_u),gamma_r*tanh(alpha_r)]';
% -------------------control law 1--------------------------
% persistent K
% if isempty(K)
%     K = diag([2 2]);
% end
% tau_c = M*(-f_hat-K*e_bar+nu_cf_dot-[Tu*alpha_u,Tr*alpha_r]');

% ----------------------Auxiliary system 2-------------------------------
persistent alpha_u alpha_r Tu Tr
if isempty(alpha_u)
    alpha_u = 0;
    alpha_r = 0;
    Tu = 0.1;
    Tr = 0.1;
end
delta_u = tauc(1)-tau(1);
delta_r = tauc(2)-tau(2);
alpha_u_dot = -Tu*alpha_u-delta_u/m11;
alpha_r_dot = -Tr*alpha_r-delta_r/m33;
alpha_u = ts*alpha_u_dot+alpha_u;
alpha_r = ts*alpha_r_dot+alpha_r;

e_bar = nu_bar-nu_cf-[alpha_u,alpha_r]';

%-------------------------control law 2---------------------
persistent K
if isempty(K)
    K = diag([2 2]);
end
tau_c = M*(-f_hat-K*e_bar+nu_cf_dot-[Tu*alpha_u,Tr*alpha_r]');

%----------------------------------------------------------

fhat = f_hat;

end

