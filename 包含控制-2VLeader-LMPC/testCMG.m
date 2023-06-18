%  ref: Gu, N., Wang, D., Peng, Z., Li, T., & Tong, S. (2021). 
%       Model-Free Containment Control of Underactuated Surface Vessels Under Switching Topologies
%       Based on Guiding Vector Fields and Data-Driven Neural Predictors. IEEE Transactions on Cybernetics, 1–12.
% by: Yinsong Qu
% date: 4, 22, 2023 

% 只考虑两个虚拟领航员

clc
clear all
close all

ts = 0.05;
tfinal = 500;
Ns = tfinal/ts;

%% parameters initializition
% desired formation
%--------------------------
%   USV2            USV4
%           USV1
%   USV3            USV5
%------------------------
%----------

% theta11d = 0; theta12d = 0; theta13d = 5; theta14d = 0;
% theta21d = 0; theta22d = 0; theta23d = 5; theta24d = 5;
% theta31d = -5; theta32d = -5; theta33d = 0; theta34d = 0;
% theta41d = 0; theta42d = -5; theta43d = 0; theta44d = 0;

theta11d = 0; theta12d = 0;
theta21d = 0; theta22d = 0;

Pijd = [theta11d,theta12d;
        theta21d,theta22d];
    
p10d = 0; p20d = 0;
Pi0d = [p10d;p20d];

Pijd = [Pijd,Pi0d];

% communication networks based on Laplacian matrix
I41 = [1 1 1 1]';
I2 = diag([1 1]);

AFF = [0 1 0 0;
       1 0 1 0;
       0 1 0 1;
       0 0 1 0];
AFL = [1 0;
       0 0;
       0 0;
       0 1];

AFS = [0 0 0 0]';
%-------------------
ALF = [1 0 0 0;
       0 0 0 1];
ALL = [0 1;
       1 0];
ALS = [1 0]';
%-------------------
ASF = [0 0 0 0];
ASL = [0 0];
ASS = 0;

A = [AFF,AFL,AFS;
     ALF,ALL,ALS;
     ASF,ASL,ASS];

%-------------------
D = diag(sum(A,2));
L = D-A;
%-----------------
AL = [ALL,ALS];
DL = diag(sum(AL,2));
DL = [DL,zeros(2,1)];
LL = DL-AL;

LF = L(1:4,1:6);
LF_bar = Kerector( LF, I2 );

for k=1:Ns
   % time series
   t = (k-1)*ts;
   
   % path info
   % pc = 1: straight line
   % pc = 2: Wang 2021
   % pc = 3: circle
   % pc = 4: straight line + circle
   if t == 0
      pc = 1;
      us = 0.3; 
      delta0 = 0.5;
      tauc1 = [0 0]';tauc2 = [0 0]';tauc3 = [0 0]';tauc4 = [0 0]';
      x10 = [0 0 0 25 0 0]'; x20 = [0 0 0 18 0 0]';
      x30 = [0 0 0 10 0 0]'; x40 = [0 0 0 0 0 0]';
      M = diag([25.8,33.8,2.76]);
      M_bar = diag([M(1,1),M(3,3)]);
   end
   
   % follower
   [xd1,diu1] = follower1( x10,tauc1,ts );
   [xd2,diu2] = follower2( x20,tauc2,ts );
   [xd3,diu3] = follower3( x30,tauc3,ts );
   [xd4,diu4] = follower4( x40,tauc4,ts );
   pd1 = [xd1(4),xd1(5)]';
   pd2 = [xd2(4),xd2(5)]';
   pd3 = [xd3(4),xd3(5)]';
   pd4 = [xd4(4),xd4(5)]';
   nu1_bar = [xd1(1),xd1(3)]';
   nu2_bar = [xd2(1),xd2(3)]';
   nu3_bar = [xd3(1),xd3(3)]';
   nu4_bar = [xd4(1),xd4(3)]';
   Rd1 = [cos(xd1(6)) -sin(xd1(6));
          sin(xd1(6)) cos(xd1(6))];
   Rd2 = [cos(xd2(6)) -sin(xd2(6));
          sin(xd2(6)) cos(xd2(6))];
   Rd3 = [cos(xd3(6)) -sin(xd3(6));
          sin(xd3(6)) cos(xd3(6))];
   Rd4 = [cos(xd4(6)) -sin(xd4(6));
          sin(xd4(6)) cos(xd4(6))];
   Rbar = [Rd1,zeros(2,6);
           zeros(2,2),Rd2,zeros(2,4);
           zeros(2,4),Rd3,zeros(2,2);
           zeros(2,6),Rd4];
   % super leader and virtual leader
   [pr0,pr0_dw,theta0,vs] = superLeader0( pc, us, ts );

   if t==0
       eTheta = [0 0]';
       z_s = [0 0]';
       Z = zeros(8,1);
   end
   [pr5,pr1_dw,theta1,w1] = virtualLeader1( pc, z_s(1), eTheta(1), vs, ts );
   [pr6,pr2_dw,theta2,w2] = virtualLeader2( pc, z_s(2), eTheta(2), vs, ts );
      
   z_s = [ALF(1,1)*pr1_dw'*Rd1*Z(1:2);
          ALF(2,4)*pr2_dw'*Rd4*Z(7:8)];
   
   Theta = [theta1,theta2,theta0]';
   eTheta = LL*Theta;

   % 计算包含误差
   delta = [delta0,0]';
   delta_bar = repmat(delta,4,1);

   P = [pd1',pd2',pd3',pd4',pr5',pr6']';
   Z = Rbar'*(LF_bar*P)+delta_bar;
   
   d1 = L(1,1); d2 = L(2,2); d3 = L(3,3); d4 = L(4,4);
   G1 = diag([d1,delta0]); G2 = diag([d2,delta0]); 
   G3 = diag([d3,delta0]); G4 = diag([d4,delta0]);
   
   
   f11 = AFF(1,2)*Rd1'*Rd2*xd2(1:2)+AFF(1,3)*Rd1'*Rd3*xd3(1:2)+...
         AFF(1,4)*Rd1'*Rd4*xd4(1:2);
   f22 = AFF(2,1)*Rd2'*Rd1*xd1(1:2)+AFF(2,3)*Rd2'*Rd3*xd3(1:2)+...
         AFF(2,4)*Rd2'*Rd4*xd4(1:2);
   f33 = AFF(3,1)*Rd3'*Rd1*xd1(1:2)+AFF(3,2)*Rd3'*Rd2*xd2(1:2)+...
         AFF(3,4)*Rd3'*Rd4*xd4(1:2);
   f44 = AFF(4,1)*Rd4'*Rd1*xd1(1:2)+AFF(4,2)*Rd4'*Rd2*xd2(1:2)+...
         AFF(4,3)*Rd4'*Rd3*xd3(1:2);
   g11 = AFL(1,1)*Rd1'*pr1_dw*vs+AFL(1,2)*Rd1'*pr2_dw*vs;
   g22 = AFL(2,1)*Rd2'*pr1_dw*vs+AFL(2,2)*Rd2'*pr2_dw*vs;
   g33 = AFL(3,1)*Rd3'*pr1_dw*vs+AFL(3,2)*Rd3'*pr2_dw*vs;
   g44 = AFL(4,1)*Rd4'*pr1_dw*vs+AFL(4,2)*Rd4'*pr2_dw*vs;
   
   if t ==0
       K1 = diag([0.1 0.1]); K2 = K1; K3 = K1; K4 = K1;
   end
   alpha1 = G1\(-K1*Z(1:2)+f11+g11-[0,d1*xd1(2)]');
   alpha2 = G2\(-K2*Z(3:4)+f22+g22-[0,d2*xd2(2)]');
   alpha3 = G3\(-K3*Z(5:6)+f33+g33-[0,d3*xd3(2)]');
   alpha4 = G4\(-K4*Z(7:8)+f44+g44-[0,d4*xd4(2)]');
   
   tauc1 = M_bar*(-5*(nu1_bar-alpha1)-diu1);
   tauc2 = M_bar*(-5*(nu2_bar-alpha2)-diu2);
   tauc3 = M_bar*(-5*(nu3_bar-alpha3)-diu3);
   tauc4 = M_bar*(-5*(nu4_bar-alpha4)-diu4);
   
%    tauc1 = [2 0.5]';
%    tauc2 = [2 0.5]';
%    tauc3 = [2 0.5]';
%    tauc4 = [2 0.5]';
   
   Zn = [norm(Z(1:2)),norm(Z(3:4)),norm(Z(5:6)),norm(Z(7:8))]';
   xout(k,:) = [t,xd1',xd2',xd3',xd4',pr0',pr5',pr6',eTheta',Zn'];
      
end
%% simulation data
t = xout(:,1);
xd1 = xout(:,2:7);
xd2 = xout(:,8:13);
xd3 = xout(:,14:19);
xd4 = xout(:,20:25);
pr0 = xout(:,26:27);
pr5 = xout(:,28:29);
pr6 = xout(:,30:31);
eTheta = xout(:,32:33);
Zn = xout(:,34:37);

figure(1)
plot(xd1(:,5),xd1(:,4),'r-',xd2(:,5),xd2(:,4),'b-',...
     xd3(:,5),xd3(:,4),'g-', xd4(:,5),xd4(:,4),'c-',...
     pr0(:,2),pr0(:,1),'k-',pr5(:,2),pr5(:,1),'m--',pr6(:,2),pr6(:,1),'m--');
 
figure(2)
plot(t,xd1(:,1),'r-',t,xd2(:,1),'b-',...
     t,xd3(:,1),'g-',t,xd4(:,1),'c-');
title('速度u');

figure(3)
plot(t,xd1(:,2),'r-',t,xd2(:,2),'b-',...
     t,xd3(:,2),'g-',t,xd4(:,2),'c-');
title('速度v');

figure(4)
plot(t,xd1(:,3),'r-',t,xd2(:,3),'b-',...
     t,xd3(:,3),'g-',t,xd4(:,3),'c-');
title('角速度r');

figure(5)
plot(t,Zn(:,1),'r-',t,Zn(:,2),'g-',t,Zn(:,3),'b-',t,Zn(:,4),'c-');
title('包含误差');

figure(6)
plot(t,eTheta(:,1),'r-',t,eTheta(:,2),'g-');
title('编队误差');














