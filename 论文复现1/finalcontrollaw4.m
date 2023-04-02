function tau = finalcontrollaw4( nu, nu_bar, nu_bard,nu_bard_d, ts )
%FINALCONTROLLAW tau = finalcontrollaw1( nu, nu_bar, nu_bard,nu_bard_d, ts )
% get the final control law
% Author : Yinsong Qu
% Date : 7 14 2022
% reference: Adaptive Dynamic Surface Control for Formations of Autonomous Surface Vehicles With Uncertain Dynamics
% inputs: nu = [u v r]'
%         nu_bar = [u r]'
%         nu_bard = [ud rd]'
%         nu_bar_d = [ud_d rd_d]'
%         ts : sample time
% states: W_hat 8x2  NN weight matrix for outputs
%         V_hat 8x8  NN weight matrix for inputs

% %% states Initializing
% persistent W_hat V_hat
% if isempty(W_hat)
%     W_hat = zeros(8,2);
%     V_hat = zeros(8,8);
% end
% 
% %% calculate final control law
% Kv = diag([1 1]);
% nu_bar_e = nu_bar-nu_bard;


%% case 2 所有参数已知不存在干扰
Kv = diag([4 4]);
u = nu(1); v = nu(2); r = nu(3);
m11 = 25.8; m22 = 33.8; m33 = 2.76; m23 = 6.2; m32 = m23;
D11 = 12+2.5*abs(u); D22 = 17+4.5*abs(v); D23 = 0.2; D32 = 0.5;
D33 = 0.5+0.1*abs(r);
C13 = -m22*v-m23*r; C23 = m11*u; C31 = -C13; C32 = -C23;
M = [m11 0 0;0 m22 m23;0 m32 m33];
D = [D11 0 0;0 D22 D23;0 D32 D33];
C = [0 0 C13;0 0 C23;C31 C32 0];
f = -C*nu-D*nu;
f_bar = [f(1) f(3)]';
nu_bar_e = nu_bar-nu_bard;
M_bar = diag([m11 m33]);
tau = -M_bar*Kv*nu_bar_e+M_bar*nu_bard_d-f_bar;



end

