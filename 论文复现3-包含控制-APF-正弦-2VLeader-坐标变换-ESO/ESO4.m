function [eta_hat_o,nu_hat_o,zeta_hat_o] = ESO4(eta_tf,tau,Gamma,w,ts)
%ESO 此处显示有关此函数的摘要
%   此处显示详细说明

psi = eta_tf(3);
persistent eta_hat nu_hat zeta_hat KO1 KO2 KO3
if isempty(eta_hat)
   eta_hat = eta_tf; 
   nu_hat = [0 0 0]';
   zeta_hat = [0 0 0]';
   KO1 = 3*w*diag([1 1 1]);
   KO2 = 3*w^2*diag([1 1 1]);
   KO3 = w^3*diag([1 1 1]);
end
R_psi = [cos(psi) -sin(psi) 0;
         sin(psi) cos(psi) 0;
         0        0        1];
m11 = 25.8; m22 = 33.8;
M = [1/m11,0,m22/Gamma];
tau_bar = [tau(1),0,tau(2)]';
eta_hat_dot = R_psi*nu_hat+KO1*(eta_tf-eta_hat);
nu_hat_dot = R_psi'*KO2*(eta_tf-eta_hat)+zeta_hat+M*tau_bar;
zeta_hat_dot = R_psi'*KO3*(eta_tf-eta_hat);

eta_hat = eta_hat_dot*ts+eta_hat;
nu_hat = nu_hat_dot*ts+nu_hat;
zeta_hat = zeta_hat_dot*ts+zeta_hat;

eta_hat_o = eta_hat;
nu_hat_o = nu_hat;
zeta_hat_o = zeta_hat;


end

