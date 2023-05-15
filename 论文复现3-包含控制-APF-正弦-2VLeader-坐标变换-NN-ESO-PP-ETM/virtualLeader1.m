function [pos,pos_dw,thetao,wo] = virtualLeader1( pc, eF, eL, wl, vs, ts )
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
% if nargin ~=5,error('input number must be 4!');end

% states initializing
% states initializing
persistent theta
if isempty(theta)
   theta = 0; 
end

% curved path
switch pc
    case 1
        %---------------直线路径-----------------
        xd = theta;
        yd = theta;
        xd_dw = 1;
        yd_dw = 1;
        xd_ddw = 0; yd_ddw = 0;
    case 2
        %---------------Wang 正弦曲线路径--------------
        xd = 8*sin(0.5*theta) + 4*theta+60;
        yd = 4*theta;
        xd_dw = 4*cos(0.5*theta)+4;
        yd_dw = 4;
        xd_ddw = -2*sin(0.5*theta); yd_ddw = 0;
    case 3
        %---------------圆形路径------------------
        R = 30;
        xd = R*cos(theta)+30; 
        yd = R*sin(theta)+30; 
        xd_dw = -R*sin(theta); 
        yd_dw = R*cos(theta);
        xd_ddw = -R*cos(theta); yd_ddw = -R*sin(theta);
    case 4
        %---------------混合路径------------------
        xd = 10;
        yd = theta;
        xd_dw = 0; xd_ddw = 0;
        yd_dw = 1; yd_ddw = 0;
        if theta >= 80 
            R = 10;
            xd = R*cos(theta-80);
            yd = R*sin(theta-80)+80;
            xd_dw = -R*sin(theta-80); xd_ddw = -R*cos(theta-80);
            yd_dw = R*cos(theta-80); yd_ddw = -R*sin(theta-80);
        end
        if theta >= 80+pi
            xd = -10;
            yd = -(theta-80-pi)+80;
            xd_dw = 0; xd_ddw = 0;
            yd_dw = -1; yd_ddw = 0;
        end
        if theta >= 160+pi
            xd = R*cos(theta-160-pi+pi);
            yd = R*sin(theta-160-pi+pi);
            xd_dw = -R*sin(theta-160-pi+pi); xd_ddw = -R*cos(theta-160-pi+pi);
            yd_dw = R*cos(theta-160-pi+pi); yd_ddw = -R*sin(theta-160-pi+pi);
        end
        if theta >= 160+2*pi
            xd = 10;
            yd = theta-160-2*pi;
            xd_dw = 0; xd_ddw = 0;
            yd_dw = 1; yd_ddw = 0;
        end 
    case 5
        xd = 4*cos(theta);
        yd = 3*theta;
        xd_dw = -4*sin(theta);
        yd_dw = 3;
        xd_ddw = -4*cos(theta); yd_ddw = 0;
    case 6
        xd = theta;
        yd = theta;
        xd_dw = 1;
        yd_dw = 1;
        xd_ddw = 0; yd_ddw = 0;
        if theta>=15
            xd = 30-theta;
            yd = theta;
            xd_dw = -1;
            yd_dw = 1;
            xd_ddw = 0; yd_ddw = 0;
        end
    case 7
        xd = 0;
        yd = theta;
        xd_dw = 0;
        yd_dw = 1;
        xd_ddw = 0; yd_ddw = 0;
end

%---------------------------------------------------------
kc = abs(xd_dw*yd_ddw-yd_dw*xd_ddw)/sqrt(xd_dw^2+yd_dw^2)^3; % 曲线路径的曲率

psip = atan2(yd_dw,xd_dw);  % 路径虚拟参考点切线与x轴夹角
persistent psipf
if isempty(psipf)
    psipf = psip;
end
psipf_dot = -(psipf-psip)/0.1;
psipf = ts*psipf_dot+psipf;
kcs = kc*sign(psipf_dot);

% 路径参数的更新
% 和参考文献中变量略有出入，w->theta, theta->w
persistent kg1 kg2 wk
if isempty(kg1)
    kg1 = 0.05;
    kg2 = 0.5;
    wk = 0;
end
% wk = -kg*(eF-eL);

wk = -kg1*eF+kg2*eL;
% wk_dot = 0.001*(-eF+3*eL-0.5*wk);
% wk = wk_dot*ts+wk;

theta_dot = vs-wk;
theta = theta_dot*ts+theta;

pos = [xd,yd]';
pos_dw = [xd_dw,yd_dw]';

thetao = theta;
wo = wk;

end





















