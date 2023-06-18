function [xo,f,tauo] = follower3( x0,tauc,ts )
%FOLLOWER1 此处显示有关此函数的摘要
%   此处显示详细说明
persistent x
if isempty(x)
    x = x0;
end
%% control input limit
tu = tauc(1); tr = tauc(2);
% first order process with time delay
persistent tu1 tr1
if isempty(tu1)
    tu1 = 0;
    tr1 = 0;
end
tu1_dot = (tu-tu1)/0.5;
tr1_dot = (tr-tr1)/0.5;

tu_dot_max = 10; tr_dot_max = 5;

if abs(tu1_dot)>=tu_dot_max
    tu1_dot = tu_dot_max*sign(tu1_dot);
end

if abs(tr1_dot)>=tr_dot_max
    tr1_dot = tr_dot_max*sign(tr1_dot);
end

tu1 = euler2(tu1_dot,tu,ts);
tr1 = euler2(tr1_dot,tr,ts);
% limit
tu_max = 1.5; tr_max = 1;

if abs(tu1)>=tu_max
    tu1 = tu_max*sign(tu1);
end

if abs(tr1)>=tr_max
    tr1 = tr_max*sign(tr1);
end

tu = tu1; tr = tr1;

R = [cos(x(6)) -sin(x(6)) 0;
     sin(x(6)) cos(x(6)) 0;
     0 0 1];

tau = [tu;0;tr];
%% states update
Xu=0.72253;  Yv = 0.88965; Nr=1.900;
m11 = 25.8;  m22 = 33.8; m33 = 2.76; 
du = Xu; dv = Yv; dr = Nr;
Xuu=1.32742; Nrr=0.75;

M = diag([m11,m22,m33]);
D = diag([du+Xuu*abs(x(1)), dv, dr+Nrr*abs(x(3))]);
C = [0 0 -m22*x(2);
     0 0 m11*x(1);
     m22*x(2) -m11*x(1) 0];
f = -M\D*x(1:3)-M\C*x(1:3);
f = [f(1),f(3)]';
x_dot = [-M\D*x(1:3)-M\C*x(1:3)+M\tau;
         R*x(1:3)];
x = x_dot*ts+x;

xo = x;
tauo = [tu,tr]';
end

