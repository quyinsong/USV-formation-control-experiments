function [ud, rd, posd] = LOS( pc, eta, nu, ud, ts )
%LOS [ud, rd, posd] = LOS( pc, eta, nu, ts )
% curved path LOS guidace law
% Author : Yinsong Qu
% Date : 7 15 2022
%[xd, yd, psid, psif, xe, ye, w] = LOS2( pc, k1, delta, x, y, psi, U, beta, w, ts )
% Author : Quyinsong
% Date: 2022 6 1
% reference : Path-Following Algorithms and Experiments for an Unmanned Surface Vehicle
%%%% INPUT:
% pc : chose curved path
% -------------------pc-----------------------------------
% 1: stright line path; 2: sinusoidal path; 3: circle path
% -------------------pc-----------------------------------
% eta = [x y psi]' : the position and heading angle of USV
% nu = [u v r]' : the speed and yaw angular rate
% ts : sample time

% check input and state dimentions
if nargin ~=5,error('input number must be 5!');end

% states initializing
persistent w
if isempty(w)
    w = 0;
end

% curved path
switch pc
    case 1
        %---------------直线路径-----------------
        xd = 10;
        yd = 2*w;
        xd_dw = 0;
        yd_dw = 2;
        xd_ddw = 0; yd_ddw = 0;
    case 2
        %---------------正弦曲线路径--------------
        xd = 8*cos(0.2*w) + 2*w;
        yd = 4*w;
        xd_dw = -1.6*sin(0.2*w)+2;
        yd_dw = 4;
        xd_ddw = -0.32*cos(0.2*w); yd_ddw = 0;
    case 3
        %---------------圆形路径------------------
        R = 15;
        xd = R*cos(w)+0; 
        yd = R*sin(w)+100; 
        xd_dw = -R*sin(w); 
        yd_dw = R*cos(w);
        xd_ddw = -R*cos(w); yd_ddw = -R*sin(w);
    case 4
        %---------------混合路径------------------
        xd = 10;
        yd = w;
        xd_dw = 0; xd_ddw = 0;
        yd_dw = 1; yd_ddw = 0;
        if w >= 80 
            R = 10;
            xd = R*cos(w-80);
            yd = R*sin(w-80)+80;
            xd_dw = -R*sin(w-80); xd_ddw = -R*cos(w-80);
            yd_dw = R*cos(w-80); yd_ddw = -R*sin(w-80);
        end
        if w >= 80+pi
            xd = -10;
            yd = -(w-80-pi)+80;
            xd_dw = 0; xd_ddw = 0;
            yd_dw = -1; yd_ddw = 0;
        end
        if w >= 160+pi
            xd = R*cos(w-160-pi+pi);
            yd = R*sin(w-160-pi+pi);
            xd_dw = -R*sin(w-160-pi+pi); xd_ddw = -R*cos(w-160-pi+pi);
            yd_dw = R*cos(w-160-pi+pi); yd_ddw = -R*sin(w-160-pi+pi);
        end
        if w >= 160+2*pi
            xd = 10;
            yd = w-160-2*pi;
            xd_dw = 0; xd_ddw = 0;
            yd_dw = 1; yd_ddw = 0;
        end 
    case 5
        xd = 4*cos(w);
        yd = 3*w;
        xd_dw = -4*sin(w);
        yd_dw = 3;
        xd_ddw = -4*cos(w); yd_ddw = 0;
    case 6
        xd = w;
        yd = w;
        xd_dw = 1;
        yd_dw = 1;
        xd_ddw = 0; yd_ddw = 0;
        if w>=15
            xd = 30-w;
            yd = w;
            xd_dw = -1;
            yd_dw = 1;
            xd_ddw = 0; yd_ddw = 0;
        end
    case 7
        xd = w;
        yd = w;
        xd_dw = 1;
        yd_dw = 1;
        xd_ddw = 0; yd_ddw = 0;
end

%---------------------------------------------------------

x = eta(1); y = eta(2); psi = eta(3); 
u = nu(1); v = nu(2); r = nu(3);

kc = abs(xd_dw*yd_ddw-yd_dw*xd_ddw)/sqrt(xd_dw^2+yd_dw^2)^3; % 曲线路径的曲率

psip = atan2(yd_dw,xd_dw);  % 路径虚拟参考点切线与x轴夹角
persistent psipf
if isempty(psipf)
    psipf = psip;
end
psipf_dot = -(psipf-psip)/0.1;
psipf = ts*psipf_dot+psipf;
kcs = kc*sign(psipf_dot);
  
xe = cos(psip) * (x - xd) + sin(psip) * (y - yd);  % USV与虚拟参考点的纵向误差
ye = -sin(psip) * (x - xd) + cos(psip) * (y - yd);  % USV与虚拟参考点的横向误差

betad = atan2(v,ud);
chi_sf = psi + betad - psip; 
k1 = 2;
U = sqrt(u^2+v^2);
up = U*cos(chi_sf) + k1*xe;   % 路径虚拟参考点的切向速度

wdot =  up / sqrt(xd_dw^2 + yd_dw^2);  % 路径参数的更新
w = euler2(wdot, w, ts);
delta = 2;

psiLos = atan2(ye,delta)+betad;
psid = psip-psiLos;
persistent psiLosf
if isempty(psiLosf)
    psiLosf = psiLos;
end
psiLosf_dot = -(psiLosf-psiLos)/0.1;
psiLosf = ts*psiLosf_dot+psiLosf;
psid_dot = kcs*up-psiLosf;

k2 = 2;
rd = -psid_dot-k2*ssa(psi-psid);
posd = [xd,yd]';

end





















