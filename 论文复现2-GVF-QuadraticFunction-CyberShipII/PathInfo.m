function [wo,thetao,vs, posd, pd_dw] = PathInfo( pc, z_theta, ud, ts )
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
if nargin ~=4,error('input number must be 4!');end

% states initializing
% states initializing
persistent theta w miu lambda
if isempty(theta)
    theta = 0;
    w = 0;
    miu = 0;
    lambda = 0.1;
end

% curved path
switch pc
    case 1
        %---------------直线路径-----------------
        xd = 2*w;
        yd = 2*w;
        xd_dw = 2;
        yd_dw = 2;
        xd_ddw = 0; yd_ddw = 0;
    case 2
        %---------------Wang 正弦曲线路径--------------
        xd = 8*sin(0.5*w) + 4*w;
        yd = 4*w;
        xd_dw = 4*cos(0.5*w)+4;
        yd_dw = 4;
        xd_ddw = -2*sin(0.5*w); yd_ddw = 0;
    case 3
        %---------------圆形路径------------------
        R = 30;
        xd = R*cos(w)+30; 
        yd = R*sin(w)+30; 
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
        xd = 30;
        yd = w-80;
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
theta_dot = -lambda*(theta+miu*z_theta);
theta = theta_dot*ts+theta;

vs = ud/sqrt(xd_dw^2+yd_dw^2);
w_dot = vs-theta;
w = w_dot*ts+w;

wo = w;
thetao = theta;

posd = [xd,yd]';
pd_dw = [xd_dw,yd_dw]';

end





















