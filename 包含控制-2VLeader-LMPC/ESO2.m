function nu_hat = ESO2(eta_tf,ts)
%ESO 此处显示有关此函数的摘要
%   此处显示详细说明

psi = eta_tf(3);
persistent P K
if isempty(P)
   K = 10*diag([1 1 2]);
   P = -K*eta_tf;
end
R_psi = [cos(psi) -sin(psi) 0;
         sin(psi) cos(psi) 0;
         0        0        1];

nu_hat = R_psi'*(K*eta_tf+P);
P_dot = -K*P-K*K*eta_tf;
P = P_dot*ts+P;


end

