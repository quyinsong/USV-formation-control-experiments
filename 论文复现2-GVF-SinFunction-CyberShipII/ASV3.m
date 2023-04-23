function [y1, y2, f] = ASV1( x0, tau, tau_w, ts)
%USV x = USV( x, tau, ts )
%   3DOF USV model 
% Author : Yinsong Qu
% Date: 7 13 2022
% Reference: Distributed Path Following of Multiple Under-Actuated Autonomous Surface Vehicles
% Based on Data-Driven Neural Predictors via Integral Concurrent Learning
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
tu_w = tau_w(1);
tv_w = tau_w(2);
tr_w = tau_w(3);
% if abs(x(3))>=0.8
%    x(3) = sign(x(3))*0.8; 
% end
%% control input
% limit
tu_max = 2; tr_max = 1.5;
tu = tau(1); tr = tau(2);
if tu>=tu_max
    tu=tu_max;
elseif tu<=0
    tu=0;
end
% if abs(tu)>=tu_max
%     tu = tu_max*sign(tu);
% end
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

%% parameters
m11 = 25.8; m22 = 33.8; m33 = 2.76;
u = x(1); v = x(2); r = x(3); xn = x(4); yn = x(5); psin = x(6);
fu = -5.87*u^3-1.33*abs(u)*u-0.72*u+m22*v*r+1.0948*r^2;
fv = -36.5*abs(v)*v-0.8896*v-0.805*v*abs(r)-m11*u*r;
fr = -0.75*abs(r)*r-1.90*r+0.08*abs(v)*r+(m11-m22)*u*v-1.0948*u*r;

%% states update
R = [cos(psin) -sin(psin) 0; sin(psin) cos(psin) 0; 0 0 1];
x1 = [u v r]'; x2 = [xn yn psin]';

x_dot = [fu/m11+(tu+tu_w)/m11
         (fv+tv_w)/m22
         fr/m33+(tr+tr_w)/m33
         R*x1];
x = euler2(x_dot,x,ts);

y1= x;
y2 =[tu tr]'; 
f = [(fu+tu_w)/m11,(fr+tr_w)/m33]';

end

