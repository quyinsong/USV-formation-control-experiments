function tau = pathfollowingcontrollaw( ud,rd,nu,tauj,ts )
%PATHFOLLOWINGCONTROLLAW tau = pathfollowingcontrollaw( ud,rd,ts )
%   path following control law
% Inputs: ud : desired velocity
%         rd : desired yaw angular rate
%         ts : sample time
%% 所有参数已知不存在干扰
persistent x2 u_hat r_hat fu_hat fr_hat
if isempty(x2)
    x2 = [0 0]';
    u_hat = 0;
    r_hat = 0;
    fu_hat = 0;
    fr_hat = 0;
end
% first order filter
nu_bard = [ud,rd]';
T = diag([0.1,0.1]);

x2_dot = T\(nu_bard-x2);
x2 = euler2(x2_dot,x2,ts);

Kv = diag([1 1]);
u = nu(1); v = nu(2); r = nu(3);
m11 = 25.8; m22 = 33.8; m33 = 2.76; m23 = 6.2; m32 = m23;
D11 = 12+2.5*abs(u); D22 = 17+4.5*abs(v); D23 = 0.2; D32 = 0.5;
D33 = 0.5+0.1*abs(r);
C13 = -m22*v-m23*r; C23 = m11*u; C31 = -C13; C32 = -C23;
M = [m11 0 0;0 m22 m23;0 m32 m33];
D = [D11 0 0;0 D22 D23;0 D32 D33];
C = [0 0 C13;0 0 C23;C31 C32 0];

% ESO
tau_u = tauj(1);
tau_r = tauj(2);
u_hat_dot = fu_hat+tau_u/m11-5*(u_hat-u);
r_hat_dot = fr_hat+tau_r/m33-5*(r_hat-r);
fu_hat_dot = -20*(u_hat-u);
fr_hat_dot = -20*(r_hat-r);
u_hat = euler2(u_hat_dot,u_hat,ts);
r_hat = euler2(r_hat_dot,r_hat,ts);
fu_hat = euler2(fu_hat_dot,fu_hat,ts);
fr_hat = euler2(fr_hat_dot,fr_hat,ts);

% controller
ue = u-ud;
re = r-rd;
tu = -2*ue*m11-fu_hat*m11+x2_dot(1)*m11;
tr = -1*re*m33-fr_hat*m33+x2_dot(2)*m33;

tau = [tu tr]';




end

