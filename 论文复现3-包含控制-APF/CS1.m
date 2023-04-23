function [y1, y2, f] = CS1( x0, tau, tau_w, ts )
% USV01 xdot = USV( x,tao,taod ) returns the time derivative of 
% the state vector: x = [ u v r x y psi]'  for USV, where
% INPUT: 
% u=x(1): surge velocity (m/s)
% v=x(2): sway velocity (m/s)
% r=x(3): yaw velocity (rad/s)
% x=x(4): position in {n} x-direction (m)
% y=x(5): position in {n} y-direction (m)
% psai=x(6): yaw angle (rad)
% tao=[tx ty tn]':
% wind=[Vw betaw]'
% current=[Vc betac]'
% OUTPUT: 
% xdot=[udot vdot rdot xdot ydot psaidot Vcdot]':time derivative of state vector
% CS2: mass = 15kg, L = 1.255m, maximum surge force = 2N, maximum yaw force
% = 1.5N, MCLab: L = 40m, B = 6.5m
% 
% Author: Quyinsong
% Data: 13rd Jan 2022
% Reference:[1] Adaptive maneuvering, with experiments, for a model ship in a marine control laboratory

% update in 11st April 2023
% Revision: delete wind,current

% check input and state dimentions
if nargin ~=4,error('input number must be 4!');end
if length(x0) ~=6,error('state x number must be 6!');end
if length(tau) ~=2,error('ctr input tao number must be 3!');end
if length(tau_w) ~=3,error('diturbance taod number must be 3!');end

%% USV state:
persistent x
if isempty(x)
    x = x0;
end
u=x(1);
v=x(2);
r=x(3);
nu = [u v r]';
psi=x(6);
%% control input
tu = tau(1); tr = tau(2);
% first order process with time delay
persistent tu1 tr1
if isempty(tu1)
    tu1 = 0;
    tr1 = 0;
end
tu1_dot = (tu-tu1)/0.2;
tr1_dot = (tr-tr1)/0.2;

tu_dot_max = 2; tr_dot_max = 1;

if abs(tu1_dot)>=tu_dot_max
    tu1_dot = tu_dot_max*sign(tu1_dot);
end

if abs(tr1_dot)>=tr_dot_max
    tr1_dot = tr_dot_max*sign(tr1_dot);
end

tu1 = euler2(tu1_dot,tu,ts);
tr1 = euler2(tr1_dot,tr,ts);
% limit
tu_max = 2; tr_max = 1.5;

if abs(tu1)>=tu_max
    tu1 = tu_max*sign(tu1);
end

if abs(tr1)>=tr_max
    tr1 = tr_max*sign(tr1);
end

tu = tu1; tr = tr1;
%% USV parameters 
% m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
% Nrdot = -1; xg = 0.046; Yrdot = 0;
m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
Nrdot = -1; xg = 0; Yrdot = 0;
% ------------------------------------------------------
Xu=-0.72253;         Yv=-0.88965;          Nv=0.1052;
Xuu=-1.32742;        Yr=0.1079;            Nr=-1.900;
Xuuu = -5.8664;      Yvv=-36.47287;        Nvv=5.0437;
                     Yrv=-0.805;           Nrv=0.130;
                     Yvr=-0.845;           Nvr=0.080;
                     Yrr=-3.45;            Nrr=-0.75;               
% ----------------------------------------------------
m11 = m-Xudot; 
m22 = m-Yvdot;
m23 = m*xg-Yrdot;
m32 = m*xg-Nvdot;
m33 = Iz-Nrdot;
% -----------------------------------------------------
c13=-m23*r-m22*v; c23 = m11*u;
c31 = -c13; c32 = -c23;
% -----------------------------------------------------
d11=-Xu-Xuu*abs(u)-Xuuu*u^3;
d22=-Yv-Yvv*abs(v)-Yrv*abs(r);
d23=-Yr-Yvr*abs(v)-Yrr*abs(r);
d32=-Nv-Nvv*abs(v)-Nrv*abs(r);
d33=-Nr-Nvr*abs(v)-Nrr*abs(r);


%% matrix expression

Rbn=[ cos(psi) -sin(psi) 0;
      sin(psi) cos(psi)  0;
      0          0         1];
M=[m11  0    0;
   0   m22 m23;
   0   m32 m33];
Crb=[0     0     c13;
     0     0     c23;
    c31   c32     0 ];
Dv=[d11  0     0;
     0    d22  d23;
     0    d32  d33];
% time derivatives
tau = [tu,0,tr]';
nu_dot=M\(-Crb*nu-Dv*nu+tau+tau_w) ;
eta_dot=Rbn*nu;
xdot=[nu_dot ;eta_dot];

% state update
x = euler2( xdot,x,ts );

%% output
F = M\(-Crb*nu-Dv*nu+tau_w);
y1 = x;
y2 = [tu,tr]';
f = [F(1),F(3)]';

%% components expression
% diturbance
% fdu = tau_w(1);
% fdv = tau_w(2);
% fdr = tau_w(3);
% detM2=m22*m33-m23*m32;
% m0=detM2;
% 
% fu = (-c13*r-d11*u)/m11;
% fv = (m23*c31*u+m23*c32*v-m33*c23*r-(m33*d22-m23*d32)*v-(m33*d23-m23*d33)*r-m23*tr)/m0;
% fr = (-m22*c31*u-m22*c32*v+m32*c23*r-(-m32*d22+m22*d32)*v-(-m32*d23+m22*d33)*r)/m0;
% 
% % time derivatives
% 
% xdot = [fu+tu/m11+fdu/m11;
%         fv+fdv/m22;
%         fr+tr*m22/m0+fdr/m0;
%         u*cos(psai)-v*sin(psai)+Vx;
%         u*sin(psai)+v*cos(psai)+Vy;
%         r];

end



