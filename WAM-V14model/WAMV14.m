function [y1, Thrust, f] = WAMV14( x0, Thrustc, tau_w, ts )
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
% Reference: [1] Adaptive maneuvering, with experiments, for a model ship in a marine control laboratory

% update in 11st April 2023
% Revision: delete wind,current

% check input and state dimentions
if nargin ~=4,error('input number must be 4!');end
if length(x0) ~=6,error('state x number must be 6!');end
if length(Thrustc) ~=2,error('ctr input tao number must be 2!');end
if length(tau_w) ~=3,error('diturbance taod number must be 3!');end
%% USV parameters
L = 4.29;
LWL = 3.21;
T = 0.127;
LGG = 1.27;
BOA = 2.2;
B = 1.83;
Bhull = BOA-B;
rho = 10^3;
m = 150;
Cd = 1.1;

%% USV state:
persistent x
if isempty(x)
    x = x0;
end
u=x(1);
v=x(2);
r=x(3);
U = sqrt(u^2+v^2);
nu = [u v r]';
psi=x(6);
%% control input
Tl = Thrustc(1); Tr = Thrustc(2);
% first order process with time delay
persistent Tl1 Tr1
if isempty(Tl1)
    Tl1 = 0;
    Tr1 = 0;
end
Tl1_dot = (Tl-Tl1)/0.5;
Tr1_dot = (Tr-Tr1)/0.5;

Tl1_dot_max = 100; Tr1_dot_max = 100;

if abs(Tl1_dot)>=Tl1_dot_max
    Tl1_dot = Tl1_dot_max*sign(Tl1_dot);
end

if abs(Tr1_dot)>=Tr1_dot_max
    Tr1_dot = Tr1_dot_max*sign(Tr1_dot);
end

Tl1 = euler2(Tl1_dot,Tl1,ts);
Tr1 = euler2(Tr1_dot,Tr1,ts);
% limit
Tl1_max = 204; Tr1_max = 204;

if abs(Tl1)>=Tl1_max
    Tl1 = Tl1_max*sign(Tl1);
end

if abs(Tr1)>=Tr1_max
    Tr1 = Tr1_max*sign(Tr1);
end

Tl = Tl1; Tr = Tr1;
%% USV parameters 
% m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
% Nrdot = -1; xg = 0.046; Yrdot = 0;
m = 150; Xudot = -11.25; Nvdot = -679.8304; Iz = 600; Yvdot = -195.6398; 
Nrdot = -600.2102; xg = 0; Yrdot = -54.3864;
% ------------------------------------------------------
Xu=-50.5870;         Yvv=-599.3130;        Nvv=-524.3989;
Xuu=-5.8722;         Yrv=-524.3989;        Nrv=-1378;
                     Yvr=-524.3989;        Nvr=-1378;
                     Yrr=-1378;            Nrr=-2996;         
                  
                       
Yv = -20*rho*abs(v)*(1.1+0.045*L/T-0.1*Bhull/T+0.016*(Bhull/T)^2)*(pi*T*L/2);
Nr = -0.02*pi*rho*U*T^2*L^2;
Nv = -0.06*pi*rho*U*T^2*L;
Yr = -6*pi*rho*U*T^2*L;
% ----------------------------------------------------
m11 = m-Xudot; 
m22 = m-Yvdot;
m23 = m*xg-Yrdot;
m32 = m*xg-Nvdot;
m33 = Iz-Nrdot;
% -----------------------------------------------------
c13=-m*v+(Yvdot*v+(Yrdot+Nvdot)*r/2)/200; c23 = m*u-Xudot*u;
c31 = -c13; c32 = -c23;
% -----------------------------------------------------
d11=-Xu-Xuu*abs(u);
d22=-Yv-Yvv*abs(v)-Yvr*abs(r);
d23=-Yr-Yrv*abs(v)-Yrr*abs(r);
d32=-Nv-Nvv*abs(v)-Nvr*abs(r);
d33=-Nr-Nrv*abs(v)-Nrr*abs(r);


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
Thrust = [Tl,Tr]';
BT = [1 1;0 0; B/2 -B/2];
tau = BT*Thrust;
% 对精度进行限制，否则直行会有误差
if abs(tau(3))<0.0000001
    tau(3) = 0;
end
nu_dot=M\(-Crb*nu-Dv*nu+tau+tau_w) ;
% tau = [200;0;30];
% nu_dot=M\(-Crb*nu-Dv*nu+tau+tau_w) ;
eta_dot=Rbn*nu;
xdot=[nu_dot ;eta_dot];

% state update
x = euler2( xdot,x,ts );

%% output
F = M\(-Crb*nu-Dv*nu+tau_w);
y1 = x;
f = [F(1),F(3)]';

end















