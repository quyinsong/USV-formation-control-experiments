%  ref: L. Liu, D. Wang, Z. Peng and Q. -L. Han, "Distributed Path Following of Multiple Under-Actuated Autonomous Surface Vehicles 
%       Based on Data-Driven Neural Predictors via Integral Concurrent Learning," 
%       in IEEE Transactions on Neural Networks and Learning Systems, vol. 32, no. 12, pp. 5334-5344, Dec. 2021

% by: Yinsong Qu
% date: 4, 7, 2023 

% GVF: realize formation with constrainted boundaries
% date: 4, 11, 2023
% by: Yinsong Qu



clc
clear all
close all

ts = 0.05;
tfinal = 600;
Ns = tfinal/ts;

%% parameters initializition
% USV states
USV1.x0 = [0 0 0 5 -3 90*pi/180]';
USV2.x0 = [0 0 0 8 -2 90*pi/180]';
USV3.x0 = [0 0 0 0 -1 90*pi/180]';
USV4.x0 = [0 0 0 -2 -2 90*pi/180]';
USV5.x0 = [0 0 0 -8 -4 90*pi/180]';
% USV inputs
USV1.tauc=[0 0]';
USV2.tauc=[0 0]';
USV3.tauc=[0 0]';
USV4.tauc= [0 0]';
USV5.tauc=[0 0]';
% desired formation
%--------------------------
%   USV2            USV4
%           USV1
%   USV3            USV5
%------------------------
%----------
% Pijd = [p11d,p12d,p13d,p14d,p15d;
%         p21d,p22d,p23d,p24d,p25d;
%         p31d,p32d,p33d,p34d,p35d;
%         p41d,p42d,p43d,p44d,p45d;
%         p51d,p52d,p53d,p54d,p55d;];
% Pi0d = diag{[p10d,p20d,p30d,p40d,p50d]};
p11d = [0;0]; p12d = [0;0]; p13d = [0;0]; p14d = [0;0]; p15d = [0;0];
p21d = [10;-10]; p22d = [0;0]; p23d = [0;0]; p24d = [0;0]; p25d = [0;0];
p31d = [-10;-10]; p32d = [0;0]; p33d = [0;0]; p34d = [0;0]; p35d = [0;0];
p41d = [-10;10]; p42d = [0;0]; p43d = [0;0]; p44d = [0;0]; p45d = [0;0];
p51d = [10;10]; p52d = [0;0]; p53d = [0;0]; p54d = [0;0]; p55d = [0;0];
Pijd = [p11d,p12d,p13d,p14d,p15d;
        p21d,p22d,p23d,p24d,p25d;
        p31d,p32d,p33d,p34d,p35d;
        p41d,p42d,p43d,p44d,p45d;
        p51d,p52d,p53d,p54d,p55d;];
p10d = [0;0]; p20d = [0;0]; p30d = [0;0]; p40d = [0;0]; p50d = [0;0];
Pi0d = [p10d;p20d;p30d;p40d;p50d];

% communication networks based on Laplacian matrix
I21 = [1 1]';
I2 = diag([1 1]);

A = [0 0 0 0 0;
     1 0 0 0 0;
     1 0 0 0 0;
     1 0 0 0 0;
     1 0 0 0 0];
B = diag([1 0 0 0 0]);
A_bar = Kerector(A,I2);
B_bar = Kerector(B,I2);
D = diag(sum(A,2));
L = D-A;
H = L+B;
H_bar = Kerector(H,I2);

I51 = [1;1;1;1;1];

for k=1:Ns
   % time series
   t = (k-1)*ts;
   
   % path info
   % pc = 1: straight line
   % pc = 2: Wang 2021
   % pc = 3: circle
   % pc = 4: straight line + circle
   % pc = 4: cos
   if t == 0
      pc = 2;
      ud = 0.3; 
      delta0 = 0.3;
      z_w = 0;
   end
   tau_w = [0.3*sin(0.3*t)*cos(0.2*t)+0.2+0.1*randn(),...
            0.1*sin(0.3*t)*cos(0.2*t)+0.1*randn(),...
            0.1*sin(0.2*t)*cos(0.3*t)+0.1*randn()]';
   % USV states update
   [USV1.x,USV1.tau,USV1.f] = CS1( USV1.x0, USV1.tauc, tau_w, ts);
   [USV2.x,USV2.tau,USV2.f] = CS2( USV2.x0, USV2.tauc, tau_w, ts);
   [USV3.x,USV3.tau,USV3.f] = CS3( USV3.x0, USV3.tauc, tau_w, ts);
   [USV4.x,USV4.tau,USV4.f] = CS4( USV4.x0, USV4.tauc, tau_w, ts);
   [USV5.x,USV5.tau,USV5.f] = CS5( USV5.x0, USV5.tauc, tau_w, ts);

   % USV states
   USV1.eta = [USV1.x(4) USV1.x(5) USV1.x(6)]';
   USV2.eta = [USV2.x(4) USV2.x(5) USV2.x(6)]';
   USV3.eta = [USV3.x(4) USV3.x(5) USV3.x(6)]';
   USV4.eta = [USV4.x(4) USV4.x(5) USV4.x(6)]';
   USV5.eta = [USV5.x(4) USV5.x(5) USV5.x(6)]';
    
   USV1.nu = [USV1.x(1) USV1.x(2) USV1.x(3)]';USV1.nu_bar = [USV1.x(1) USV1.x(3)]';
   USV2.nu = [USV2.x(1) USV2.x(2) USV2.x(3)]';USV2.nu_bar = [USV2.x(1) USV2.x(3)]';
   USV3.nu = [USV3.x(1) USV3.x(2) USV3.x(3)]';USV3.nu_bar = [USV3.x(1) USV3.x(3)]';
   USV4.nu = [USV4.x(1) USV4.x(2) USV4.x(3)]';USV4.nu_bar = [USV4.x(1) USV4.x(3)]';
   USV5.nu = [USV5.x(1) USV5.x(2) USV5.x(3)]';USV5.nu_bar = [USV5.x(1) USV5.x(3)]';
   
   % Path Info
   [w,theta,vs,p0,p0_dw] = PathInfo( pc, z_w, ud, ts );
   p0_dw_bar = [p0_dw',p0_dw',p0_dw',p0_dw',p0_dw']';
   
   % P = [p1;p2;p3;p4;p5];
   p1 = [USV1.eta(1);USV1.eta(2)];
   p2 = [USV2.eta(1);USV2.eta(2)];
   p3 = [USV3.eta(1);USV3.eta(2)];
   p4 = [USV4.eta(1);USV4.eta(2)];
   p5 = [USV5.eta(1);USV5.eta(2)];
   P = [p1;p2;p3;p4;p5];
   p0_bar = Kerector(I51,p0);
   
   % calculate formation errors
   % z1 = [z_usv1, z_usv2, z_usv3, z_usv4, z_usv5 ];
%    z1 = H_bar*P-B_bar*p0_bar-B_bar*Pi0d-sum(Kerector(A,I21).*Pijd,2);
   z1 = H_bar*(P-p0_bar)-B_bar*Pi0d-sum(Kerector(A,I21).*Pijd,2);
   R_psi1 = [cos(USV1.eta(3)) -sin(USV1.eta(3));
             sin(USV1.eta(3)) cos(USV1.eta(3))];
   R_psi2 = [cos(USV2.eta(3)) -sin(USV2.eta(3));
             sin(USV2.eta(3)) cos(USV2.eta(3))];
   R_psi3 = [cos(USV3.eta(3)) -sin(USV3.eta(3));
             sin(USV3.eta(3)) cos(USV3.eta(3))];
   R_psi4 = [cos(USV4.eta(3)) -sin(USV4.eta(3));
             sin(USV4.eta(3)) cos(USV4.eta(3))];
   R_psi5 = [cos(USV5.eta(3)) -sin(USV5.eta(3));
             sin(USV5.eta(3)) cos(USV5.eta(3))];
   R_psi = [R_psi1',zeros(2,8);
            zeros(2,2),R_psi2',zeros(2,6);
            zeros(2,4),R_psi3',zeros(2,4);
            zeros(2,6),R_psi4',zeros(2,2);
            zeros(2,8),R_psi5'];
   z2 = R_psi*z1;
   d = Kerector(I51,[delta0,0]');
   z_bar = z2+d;
   z_w = sum(B*Kerector(I51,p0_dw_bar'*R_psi*z_bar),1);
   
   % Guidance
   di_bar = sum(A+B,2);
   h1 = diag([di_bar(1),delta0]);
   h2 = diag([di_bar(2),delta0]);
   h3 = diag([di_bar(3),delta0]);
   h4 = diag([di_bar(4),delta0]);
   h5 = diag([di_bar(5),delta0]);
   
   [USV1.nu_c,USV1.sigmahat] = Guidance1( h1, B(1,1), z_bar(1:2,1), USV1.nu_bar, vs, R_psi1, theta, USV1.nu(3), p0, ts);
   [USV2.nu_c,USV2.sigmahat] = Guidance2( h2, B(2,2), z_bar(3:4,1), USV2.nu_bar, vs, R_psi2, theta, USV2.nu(3), p0, ts);
   [USV3.nu_c,USV3.sigmahat] = Guidance3( h3, B(3,3), z_bar(5:6,1), USV3.nu_bar, vs, R_psi3, theta, USV3.nu(3), p0, ts);
   [USV4.nu_c,USV4.sigmahat] = Guidance4( h4, B(4,4), z_bar(7:8,1), USV4.nu_bar, vs, R_psi4, theta, USV4.nu(3), p0, ts);
   [USV5.nu_c,USV5.sigmahat] = Guidance5( h5, B(5,5), z_bar(9:10,1), USV5.nu_bar, vs, R_psi5, theta, USV5.nu(3), p0, ts);

   % constraint based on GVF
   USV1.beta = atan2(USV1.nu(2),USV1.nu(1));
   USV2.beta = atan2(USV2.nu(2),USV2.nu(1));
   USV3.beta = atan2(USV3.nu(2),USV3.nu(1));
   USV4.beta = atan2(USV4.nu(2),USV4.nu(1));
   USV5.beta = atan2(USV5.nu(2),USV5.nu(1));
   
   if t>=150 && t<=600
   USV1.nu_c(2) = GVF1( USV1.eta(1),USV1.eta(2),USV1.eta(3),USV1.beta,USV1.nu_c(2) );
   USV2.nu_c(2) = GVF2( USV2.eta(1),USV2.eta(2),USV2.eta(3),USV2.beta,USV2.nu_c(2) );
   USV3.nu_c(2) = GVF3( USV3.eta(1),USV3.eta(2),USV3.eta(3),USV3.beta,USV3.nu_c(2) );
   USV4.nu_c(2) = GVF4( USV4.eta(1),USV4.eta(2),USV4.eta(3),USV4.beta,USV4.nu_c(2) );
   USV5.nu_c(2) = GVF5( USV5.eta(1),USV5.eta(2),USV5.eta(3),USV5.beta,USV5.nu_c(2) );
   end
   % calculate guidance disturbances
   uj = [USV1.nu(1),USV2.nu(1),USV3.nu(1),USV4.nu(1),USV5.nu(1)]';
   vj = [USV1.nu(2),USV2.nu(2),USV3.nu(2),USV4.nu(2),USV5.nu(2)]';
   if t==0
       Pijdf = Pijd;
       Pi0df = Pi0d;
   end
   Pijdf_dot = -(Pijdf-Pijd)/0.1;
   Pijdf = Pijdf_dot*ts+Pijdf;
   p1jd_dot = Pijdf_dot(1:2,:); p2jd_dot = Pijdf_dot(3:4,:); p3jd_dot = Pijdf_dot(5:6,:);
   p4jd_dot = Pijdf_dot(7:8,:); p5jd_dot = Pijdf_dot(9:10,:);
   
   Pi0df_dot = -(Pi0df-Pi0d)/0.1;
   Pi0df = Pi0df_dot*ts+Pi0df;
   
   USV1.sigma_1 = -R_psi1'*(A(1,1)*R_psi1*[uj(1),vj(1)]'+A(1,2)*R_psi2*[uj(2),vj(2)]'+A(1,3)*R_psi3*[uj(3),vj(3)]'...
                  +A(1,4)*R_psi4*[uj(4),vj(4)]'+A(1,5)*R_psi5*[uj(5),vj(5)]')-R_psi1'*p2jd_dot(1:2,1);
               
   USV2.sigma_1 = -R_psi2'*(A(2,1)*R_psi1*[uj(1),vj(1)]'+A(2,2)*R_psi2*[uj(2),vj(2)]'+A(2,3)*R_psi3*[uj(3),vj(3)]'...
                  +A(2,4)*R_psi4*[uj(4),vj(4)]'+A(2,5)*R_psi5*[uj(5),vj(5)]')-R_psi2'*p2jd_dot(1:2,1);
               
   USV3.sigma_1 = -R_psi3'*(A(3,1)*R_psi1*[uj(1),vj(1)]'+A(3,2)*R_psi2*[uj(2),vj(2)]'+A(3,3)*R_psi3*[uj(3),vj(3)]'...
                  +A(3,4)*R_psi4*[uj(4),vj(4)]'+A(3,5)*R_psi5*[uj(5),vj(5)]')-R_psi3'*p3jd_dot(1:2,1);
               
   USV4.sigma_1 = -R_psi4'*(A(4,1)*R_psi1*[uj(1),vj(1)]'+A(4,2)*R_psi2*[uj(2),vj(2)]'+A(4,3)*R_psi3*[uj(3),vj(3)]'...
                  +A(4,4)*R_psi4*[uj(4),vj(4)]'+A(4,5)*R_psi5*[uj(5),vj(5)]')-R_psi4'*p4jd_dot(1:2,1);
               
   USV5.sigma_1 = -R_psi5'*(A(5,1)*R_psi1*[uj(1),vj(1)]'+A(5,2)*R_psi2*[uj(2),vj(2)]'+A(5,3)*R_psi3*[uj(3),vj(3)]'...
                  +A(5,4)*R_psi4*[uj(4),vj(4)]'+A(5,5)*R_psi5*[uj(5),vj(5)]')-R_psi5'*p5jd_dot(1:2,1);
               
   USV1.sigma = USV1.sigma_1+[0,di_bar(1)*vj(1)]'-B(1,1)*(vs-theta)*R_psi1'*p0_dw-B(1,1)*Pi0df_dot(1:2,1);
   USV2.sigma = USV2.sigma_1+[0,di_bar(2)*vj(2)]'-B(2,2)*(vs-theta)*R_psi2'*p0_dw-B(2,2)*Pi0df_dot(3:4,1);
   USV3.sigma = USV3.sigma_1+[0,di_bar(3)*vj(3)]'-B(3,3)*(vs-theta)*R_psi3'*p0_dw-B(3,3)*Pi0df_dot(5:6,1);
   USV4.sigma = USV4.sigma_1+[0,di_bar(4)*vj(4)]'-B(4,4)*(vs-theta)*R_psi4'*p0_dw-B(4,4)*Pi0df_dot(7:8,1);
   USV5.sigma = USV5.sigma_1+[0,di_bar(5)*vj(5)]'-B(5,5)*(vs-theta)*R_psi5'*p0_dw-B(5,5)*Pi0df_dot(9:10,1);
   
   % Control
   [USV1.tauc,USV1.fhat] = ctr1( USV1.nu_bar, USV1.nu_c, USV1.tau, USV1.tauc, ts );
   [USV2.tauc,USV2.fhat] = ctr2( USV2.nu_bar, USV2.nu_c, USV2.tau, USV2.tauc, ts );
   [USV3.tauc,USV3.fhat] = ctr3( USV3.nu_bar, USV3.nu_c, USV3.tau, USV3.tauc, ts );
   [USV4.tauc,USV4.fhat] = ctr4( USV4.nu_bar, USV4.nu_c, USV4.tau, USV4.tauc, ts );
   [USV5.tauc,USV5.fhat] = ctr5( USV5.nu_bar, USV5.nu_c, USV5.tau, USV5.tauc, ts );
   
%    % change formation
%    if t>=500
%        p11d = [0;0]; p12d = [0;0]; p13d = [0;0]; p14d = [0;0]; p15d = [0;0];
%        p21d = [0;5]; p22d = [0;0]; p23d = [0;0]; p24d = [0;0]; p25d = [0;0];
%        p31d = [0;10]; p32d = [0;5]; p33d = [0;0]; p34d = [0;0]; p35d = [0;0];
%        p41d = [0;15]; p42d = [0;0]; p43d = [0;5]; p44d = [0;0]; p45d = [0;0];
%        p51d = [0;20]; p52d = [0;0]; p53d = [0;0]; p54d = [0;5]; p55d = [0;0];
%        Pijd = [p11d,p12d,p13d,p14d,p15d;
%                p21d,p22d,p23d,p24d,p25d;
%                p31d,p32d,p33d,p34d,p35d;
%                p41d,p42d,p43d,p44d,p45d;
%                p51d,p52d,p53d,p54d,p55d;];
%    end
%    
%    % change communication networks
%    if t>=500
%        A = [0 0 0 0 0;
%             1 0 0 0 0;
%             0 1 0 0 0;
%             0 0 1 0 0;
%             0 0 0 1 0];
%        B = diag([1 0 0 0 0]);
%        A_bar = Kerector(A,I2);
%        B_bar = Kerector(B,I2);
%        D = diag(sum(A,2));
%        L = D-A;
%        H = L+B;
%        H_bar = Kerector(H,I2);
%    end
   z_bar_norm = [norm(z_bar(1:2)),norm(z_bar(3:4)),norm(z_bar(5:6)),norm(z_bar(7:8)),norm(z_bar(9:10))]';
   
   xout(k,:) = [t,USV1.x', USV2.x', USV3.x', USV4.x', USV5.x', USV1.tau' USV2.tau' USV3.tau' USV4.tau', USV5.tau',...
                USV1.tauc' USV2.tauc' USV3.tauc' USV4.tauc', USV5.tauc', p0',USV1.nu_c',USV2.nu_c',USV3.nu_c',...
                USV4.nu_c',USV5.nu_c',USV1.sigmahat',USV2.sigmahat',USV3.sigmahat',USV4.sigmahat',USV5.sigmahat',USV1.sigma',...
                USV2.sigma',USV3.sigma',USV4.sigma',USV5.sigma',USV1.fhat',USV2.fhat',USV3.fhat',USV4.fhat',USV5.fhat',...
                USV1.f',USV2.f',USV3.f',USV4.f',USV5.f',z_bar_norm'];
    
    
end
%% simulation data
t = xout(:,1);
USV1.x = xout(:,2:7);
USV2.x = xout(:,8:13);
USV3.x = xout(:,14:19);
USV4.x = xout(:,20:25);
USV5.x = xout(:,26:31);
USV1.tau = xout(:,32:33);
USV2.tau = xout(:,34:35);
USV3.tau = xout(:,36:37);
USV4.tau = xout(:,38:39);
USV5.tau = xout(:,40:41);
USV1.tauc = xout(:,42:43);
USV2.tauc = xout(:,44:45);
USV3.tauc = xout(:,46:47);
USV4.tauc = xout(:,48:49);
USV5.tauc = xout(:,50:51);
p0 = xout(:,52:53);
USV1.nu_c = xout(:,54:55);
USV2.nu_c = xout(:,56:57);
USV3.nu_c = xout(:,58:59);
USV4.nu_c = xout(:,60:61);
USV5.nu_c = xout(:,62:63);
USV1.sigmahat = xout(:,64:65);
USV2.sigmahat = xout(:,66:67);
USV3.sigmahat = xout(:,68:69);
USV4.sigmahat = xout(:,70:71);
USV5.sigmahat = xout(:,72:73);
USV1.sigma = xout(:,74:75);
USV2.sigma = xout(:,76:77);
USV3.sigma = xout(:,78:79);
USV4.sigma = xout(:,80:81);
USV5.sigma = xout(:,82:83);
USV1.fhat = xout(:,84:85);
USV2.fhat = xout(:,86:87);
USV3.fhat = xout(:,88:89);
USV4.fhat = xout(:,90:91);
USV5.fhat = xout(:,92:93);
USV1.f = xout(:,94:95);
USV2.f = xout(:,96:97);
USV3.f = xout(:,98:99);
USV4.f = xout(:,100:101);
USV5.f = xout(:,102:103);
z1_bar_norm = xout(:,104);
z2_bar_norm = xout(:,105);
z3_bar_norm = xout(:,106);
z4_bar_norm = xout(:,107);
z5_bar_norm = xout(:,108);


%% PLOTS
% formation
figure(1); hold on

fontsize = 10; fontname = 'Times New Roman';
xrange=[-10 50 200]; yrange = [-10 50 200];
linewid = 1;

Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);
plot(USV1.x(:,5),USV1.x(:,4),'r--'); 
plot(USV2.x(:,5),USV2.x(:,4),'g--'); 
plot(USV3.x(:,5),USV3.x(:,4),'b--'); 
plot(USV4.x(:,5),USV4.x(:,4),'c--'); 
plot(USV5.x(:,5),USV5.x(:,4),'m--'); 
plot(p0(:,2),p0(:,1),'k-');

% USV1-4
for k=1:2000:Ns
    pos1 = [USV1.x(k,4) USV1.x(k,5)]'; 
    pos2 = [USV2.x(k,4) USV2.x(k,5)]'; 
    pos3 = [USV3.x(k,4) USV3.x(k,5)]'; 
    pos4 = [USV4.x(k,4) USV4.x(k,5)]'; 
    pos5 = [USV5.x(k,4) USV5.x(k,5)]';
    modelplot(pos1,USV1.x(k,6),'r-',linewid);
    modelplot(pos2,USV2.x(k,6),'g-',linewid);
    modelplot(pos3,USV3.x(k,6),'b-',linewid);
    modelplot(pos4,USV4.x(k,6),'c-',linewid);
    modelplot(pos5,USV5.x(k,6),'m-',linewid);
end

% ���Ʊ߽�
xx1 = 20:1:100;
xx2 = 50:1:130;
yy1 = xx1+25;
yy2 = xx2-25;
plot(xx1,yy1,'r-',xx2,yy2,'r-','linewid',linewid);

yy1 = xx1+23;
yy2 = xx2-23;
plot(xx1,yy1,'r-.',xx2,yy2,'r-.','linewid',linewid);

yy1 = xx1+21;
yy2 = xx2-21;
plot(xx1,yy1,'b-.',xx2,yy2,'b-.','linewid',linewid);

set(gca,'xTick',Xmin:Xinterval:Xmax);
set(gca,'yTick',Ymin:Yinterval:Ymax);
axis([Xmin Xmax,Ymin Ymax]);

xlabel('y (m)','FontSize',fontsize,'FontName',fontname);
ylabel('x (m)','FontSize',fontsize,'FontName',fontname);

box on

hold off


figure(2)
subplot(511);plot(t,USV1.tau(:,1),'r-',t,USV1.tau(:,2),'b-'); ylim([-2 2]);
subplot(512);plot(t,USV2.tau(:,1),'r-',t,USV2.tau(:,2),'b-'); ylim([-2 2]);
subplot(513);plot(t,USV3.tau(:,1),'r-',t,USV3.tau(:,2),'b-'); ylim([-2 2]);
subplot(514);plot(t,USV4.tau(:,1),'r-',t,USV4.tau(:,2),'b-'); ylim([-2 2]);
subplot(515);plot(t,USV5.tau(:,1),'r-',t,USV5.tau(:,2),'b-'); ylim([-2 2]);
title('force real');

figure(3)
subplot(511);plot(t,USV1.tauc(:,1),'r-',t,USV1.tauc(:,2),'b-'); 
subplot(512);plot(t,USV2.tauc(:,1),'r-',t,USV2.tauc(:,2),'b-');
subplot(513);plot(t,USV3.tauc(:,1),'r-',t,USV3.tauc(:,2),'b-');
subplot(514);plot(t,USV4.tauc(:,1),'r-',t,USV4.tauc(:,2),'b-');
subplot(515);plot(t,USV5.tauc(:,1),'r-',t,USV5.tauc(:,2),'b-');
title('force command');

figure(4)
subplot(511);plot(t,USV1.nu_c(:,1),'r-',t,USV1.nu_c(:,2),'b-'); 
subplot(512);plot(t,USV2.nu_c(:,1),'r-',t,USV2.nu_c(:,2),'b-');
subplot(513);plot(t,USV3.nu_c(:,1),'r-',t,USV3.nu_c(:,2),'b-');
subplot(514);plot(t,USV4.nu_c(:,1),'r-',t,USV4.nu_c(:,2),'b-');
subplot(515);plot(t,USV5.nu_c(:,1),'r-',t,USV5.nu_c(:,2),'b-');
title('speed command');

figure(5)
subplot(511);plot(t,USV1.sigmahat(:,1),'r-',t,USV1.sigmahat(:,2),'b-'); hold on
plot(t,USV1.sigma(:,1),'g--',t,USV1.sigma(:,2),'c--'); hold off

subplot(512);plot(t,USV2.sigmahat(:,1),'r-',t,USV2.sigmahat(:,2),'b-'); hold on
plot(t,USV2.sigma(:,1),'g--',t,USV2.sigma(:,2),'c--'); hold off

subplot(513);plot(t,USV3.sigmahat(:,1),'r-',t,USV3.sigmahat(:,2),'b-'); hold on
plot(t,USV3.sigma(:,1),'g--',t,USV3.sigma(:,2),'c--'); hold off

subplot(514);plot(t,USV4.sigmahat(:,1),'r-',t,USV4.sigmahat(:,2),'b-'); hold on
plot(t,USV4.sigma(:,1),'g--',t,USV4.sigma(:,2),'c--'); hold off

subplot(515);plot(t,USV5.sigmahat(:,1),'r-',t,USV5.sigmahat(:,2),'b-'); hold on
plot(t,USV5.sigma(:,1),'g--',t,USV5.sigma(:,2),'c--'); hold off
title('sigma');

figure(6)
subplot(511);plot(t,USV1.fhat(:,1),'r-',t,USV1.fhat(:,2),'b-'); hold on
plot(t,USV1.f(:,1),'g--',t,USV1.f(:,2),'c--'); hold off
subplot(512);plot(t,USV2.fhat(:,1),'r-',t,USV2.fhat(:,2),'b-'); hold on
plot(t,USV2.f(:,1),'g--',t,USV2.f(:,2),'c--'); hold off
subplot(513);plot(t,USV3.fhat(:,1),'r-',t,USV3.fhat(:,2),'b-'); hold on
plot(t,USV3.f(:,1),'g--',t,USV3.f(:,2),'c--'); hold off
subplot(514);plot(t,USV4.fhat(:,1),'r-',t,USV4.fhat(:,2),'b-'); hold on
plot(t,USV4.f(:,1),'g--',t,USV4.f(:,2),'c--'); hold off
subplot(515);plot(t,USV5.fhat(:,1),'r-',t,USV5.fhat(:,2),'b-'); hold on
plot(t,USV5.f(:,1),'g--',t,USV5.f(:,2),'c--'); hold off
title('f');

figure(7)
subplot(311);plot(t,USV1.x(:,1),'r-');
subplot(312);plot(t,USV2.x(:,2),'r-');
subplot(313);plot(t,USV3.x(:,3),'r-');

title('speed real');

figure(8)
subplot(511);plot(t,USV1.x(:,1)-USV1.nu_c(:,1),'r-',t,USV1.x(:,3)-USV1.nu_c(:,2),'b-');
subplot(512);plot(t,USV2.x(:,1)-USV1.nu_c(:,1),'r-',t,USV2.x(:,3)-USV1.nu_c(:,2),'b-'); 
subplot(513);plot(t,USV3.x(:,1)-USV2.nu_c(:,1),'r-',t,USV3.x(:,3)-USV2.nu_c(:,2),'b-'); 
subplot(514);plot(t,USV4.x(:,1)-USV3.nu_c(:,1),'r-',t,USV4.x(:,3)-USV3.nu_c(:,2),'b-');
subplot(515);plot(t,USV5.x(:,1)-USV4.nu_c(:,1),'r-',t,USV5.x(:,3)-USV4.nu_c(:,2),'b-');

title('speed tracking errors');

figure(9)
plot(t,z1_bar_norm,'r-',t,z2_bar_norm,'b-',t,z3_bar_norm,'g-',t,z4_bar_norm,'c-',t,z5_bar_norm,'m-');

title('formation errors');

% animation
pathX = p0(:,1); pathY = p0(:,2);
area = [Xmin Xmax Ymin Ymax];
X = [USV1.x(:,4),USV2.x(:,4),USV3.x(:,4),USV4.x(:,4),USV5.x(:,4)];
Y = [USV1.x(:,5),USV2.x(:,5),USV3.x(:,5),USV4.x(:,5),USV5.x(:,5)];
Psi = [USV1.x(:,6),USV2.x(:,6),USV3.x(:,6),USV4.x(:,6),USV5.x(:,6),];
myAnimation( X, Y, Psi, pathX, pathY, area, 10 );












