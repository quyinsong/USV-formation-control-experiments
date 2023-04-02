function [y1,y2] = ASV2( x0, tau, ts, t )
%USV x = USV( x, tau, ts )
%   3DOF USV model 
% Author : Yinsong Qu
% Date: 7 13 2022
% Reference: Global robust adaptive path following of underactuated ships 2006
% input: x0 = [u0 v0 r0 xn0 yn0 psin0]' : USV states initial value
%        tau = [tau_u tau_r] : control inputs
%        ts : sample time
%        t : simulation time
% state:       x = [u v r xn yn psin] : USV states
% output: y = x
%% states Inintializing 
persistent x
if isempty (x)
    x = x0;
end
%% control input
tu_max = 15; tr_max = 1.5;
tu = tau(1); tr = tau(2);
if tu>=tu_max
    tu=tu_max;
elseif tu<=0
    tu=0;
end
if abs(tr)>=tr_max
    tr = tr_max*sign(tr);
end
% first order process with time delay
persistent tu1 tr1
if isempty(tu1)
    tu1 = 0;
    tr1 = 0;
end
tu1_dot = (tu-tu1)/0.1;
tr1_dot = (tr-tr1)/0.1;
tu1 = euler2(tu1_dot,tu,ts);
tr1 = euler2(tr1_dot,tr,ts);
tu = tu1; tr = tr1;
tau = [tu; 0; tr];
%% parameters
u = x(1); v = x(2); r = x(3); xn = x(4); yn = x(5); psin = x(6);
m11 = 25.8; m22 = 33.8; m33 = 2.76; m23 = 6.2; m32 = m23;
D11 = 12+2.5*abs(u); D22 = 17+4.5*abs(v); D23 = 0.2; D32 = 0.5;
D33 = 0.5+0.1*abs(r);
C13 = -m22*v-m23*r; C23 = m11*u; C31 = -C13; C32 = -C23;
M = [m11 0 0;0 m22 m23;0 m32 m33];
D = [D11 0 0;0 D22 D23;0 D32 D33];
C = [0 0 C13;0 0 C23;C31 C32 0];
g = [0.0122*u*v+0.0142*v^2*r; 0.0257*u*r+0.0193*r^2*v];
tau_w = [-3*cos(0.5*t)*cos(t)+0.3*sin(0.3*t)*cos(0.8*t)-3;
        0.1*cos(0.1*t);0.6*sin(t)*cos(t)];
%% states update
R = [cos(psin) -sin(psin) 0; sin(psin) cos(psin) 0; 0 0 1];
x1 = [u v r]'; x2 = [xn yn psin]';

x_dot = [M\(-C*x1-D*x1+tau)
         R*x1];
x = euler2(x_dot,x,ts);

y1 = x;
y2=[tu,tr]';
end

