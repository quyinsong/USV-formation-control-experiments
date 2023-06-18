%  ref: Gu, N., Wang, D., Peng, Z., & Liu, L. (2019). Distributed containment
%  maneuvering of uncertain under-actuated unmanned surface vehicles guided 
% by multiple virtual leaders with a formation. Ocean Engineering 
% by: Yinsong Qu
% date: 4, 22, 2023 

%% CMG 
% 只考虑两个虚拟领航员

clc
clear all
close all

ts = 0.1; %仿真步长
tfinal = 80;
Ns = tfinal/ts;

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
      us = 0.5; 
      delta0 = 0.5;
      taucr1 = [0 0]';taucr2 = [0 0]';taucr3 = [0 0]';taucr4 = [0 0]';
      USV1.x0 = [0 0 0 25 0 0]'; USV2.x0 = [0 0 0 18 0 0]';
      USV3.x0 = [0 0 0 10 0 0]'; USV4.x0 = [0 0 0 0 0 0]';
      x10 = USV1.x0; x20 = USV2.x0;
      x30 = USV3.x0; x40 = USV4.x0;
      M = diag([25.8,33.8,2.76]);
      M_bar = diag([M(1,1),M(3,3)]);
   end
   
   % follower
   [xd1,diu1,taur1] = follower1( x10,taucr1,ts );
   [xd2,diu2,taur2] = follower2( x20,taucr2,ts );
   [xd3,diu3,taur3] = follower3( x30,taucr3,ts );
   [xd4,diu4,taur4] = follower4( x40,taucr4,ts );
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
       K1 = diag([0.05 0.05]); K2 = K1; K3 = K1; K4 = K1;
   end
   alpha1 = G1\(-K1*Z(1:2)+f11+g11-[0,d1*xd1(2)]');
   alpha2 = G2\(-K2*Z(3:4)+f22+g22-[0,d2*xd2(2)]');
   alpha3 = G3\(-K3*Z(5:6)+f33+g33-[0,d3*xd3(2)]');
   alpha4 = G4\(-K4*Z(7:8)+f44+g44-[0,d4*xd4(2)]');
   
   taucr1 = M_bar*(-5*(nu1_bar-alpha1)-diu1);
   taucr2 = M_bar*(-5*(nu2_bar-alpha2)-diu2);
   taucr3 = M_bar*(-5*(nu3_bar-alpha3)-diu3);
   taucr4 = M_bar*(-5*(nu4_bar-alpha4)-diu4);
   
%    tauc1 = [2 0.5]';
%    tauc2 = [2 0.5]';
%    tauc3 = [2 0.5]';
%    tauc4 = [2 0.5]';
   
   Zn = [norm(Z(1:2)),norm(Z(3:4)),norm(Z(5:6)),norm(Z(7:8))]';
   xrout(k,:) = [t,xd1',xd2',xd3',xd4',pr0',pr5',pr6',eTheta',Zn',taur1',taur2',taur3',taur4'];
      
end

% data
t = xrout(:,1);
xd1 = xrout(:,2:7);
xd2 = xrout(:,8:13);
xd3 = xrout(:,14:19);
xd4 = xrout(:,20:25);
pr0 = xrout(:,26:27);
pr5 = xrout(:,28:29);
pr6 = xrout(:,30:31);
eTheta = xrout(:,32:33);
Zn = xrout(:,34:37);
taur1 = xrout(:,38:39);
taur2 = xrout(:,40:41);
taur3 = xrout(:,42:43);
taur4 = xrout(:,44:45);

tm = 0.5; % 实际采样时间
kk = 1;
for k = 1:tm/ts:length(xd1(:,1))
    Xrm1(6*(kk-1)+1:6*(kk-1)+6,1) = xd1(k,:)';
    taucrm1(2*(kk-1)+1:2*(kk-1)+2,1) = taur1(k,:)';
    xrm(kk,1) = xd1(k,4);
    yrm(kk,1) = xd1(k,5);
    kk = kk+1;
end
Nrm1 = length(Xrm1(:,1));

kk = 1;
for k = 1:tm/ts:length(xd2(:,1))
    Xrm2(6*(kk-1)+1:6*(kk-1)+6,1) = xd2(k,:)';
    taucrm2(2*(kk-1)+1:2*(kk-1)+2,1) = taur2(k,:)';
    kk = kk+1;
end
Nrm2 = length(Xrm2(:,1));

kk = 1;
for k = 1:tm/ts:length(xd3(:,1))
    Xrm3(6*(kk-1)+1:6*(kk-1)+6,1) = xd3(k,:)';
    taucrm3(2*(kk-1)+1:2*(kk-1)+2,1) = taur3(k,:)';
    kk = kk+1;
end
Nrm3 = length(Xrm3(:,1));

kk = 1;
for k = 1:tm/ts:length(xd4(:,1))
    Xrm4(6*(kk-1)+1:6*(kk-1)+6,1) = xd4(k,:)';
    taucrm4(2*(kk-1)+1:2*(kk-1)+2,1) = taur4(k,:)';
    kk = kk+1;
end
Nrm4 = length(Xrm4(:,1));

%% MPC轨迹跟踪
for k = 1:Ns
   t = (k-1)*ts;
   if t==0
       USV1.tauc = [0 0]'; USV2.tauc = [0 0]';
       USV3.tauc = [0 0]'; USV4.tauc = [0 0]';
       N = 3;
   end
   t
   % plant
   tauwl = [0.3*sin(0.08*t)*cos(0.05*t);0.2*sin(0.08*t)*cos(0.05*t);0.1*sin(0.08*t)*cos(0.05*t)];
   if t==0
       wk_u = 0;
       wk_v = 0;
       wk_r = 0;
   end
   miu_u = 5; miu_v = 3; miu_r = 4;
   wu = normrnd(0,0.1); wv = normrnd(0,0.1); wr = normrnd(0,0.1); 
   wk_u_dot = -wk_u+miu_u*wu;
   wk_v_dot = -wk_v+miu_v*wv;
   wk_r_dot = -wk_r+miu_r*wr;
   wk_u = wk_u_dot*ts+wk_u;
   wk_v = wk_v_dot*ts+wk_v;
   wk_r = wk_r_dot*ts+wk_r;  
   tauwh = [wk_u,wk_v,wk_r]';
   
%    tauw = tauwh+tauwl;
   tauw = [0 0 0]';
   
   [USV1.X, USV1.tau, USV1.f] = CS1( USV1.x0, USV1.tauc, tauw, ts );
   [USV2.X, USV2.tau, USV2.f] = CS2( USV2.x0, USV2.tauc, tauw, ts );
   [USV3.X, USV3.tau, USV3.f] = CS3( USV3.x0, USV3.tauc, tauw, ts );
   [USV4.X, USV4.tau, USV4.f] = CS4( USV4.x0, USV4.tauc, tauw, ts );
   [USV1.Xm,USV1.Nm]  = Navigation1( USV1.X, tm, ts );
   [USV2.Xm,USV2.Nm]  = Navigation2( USV2.X, tm, ts );
   [USV3.Xm,USV3.Nm]  = Navigation3( USV3.X, tm, ts );
   [USV4.Xm,USV4.Nm]  = Navigation4( USV4.X, tm, ts );
   % LMPC
   Q = diag([50 0.01 10 500 500 10]);
   R = diag([5 5]);
   F = diag([50 0.01 10 500 500 10]);
   if t == 0
       USV1.Uk_1 = zeros(2*N,1);
       USV2.Uk_1 = zeros(2*N,1);
       USV3.Uk_1 = zeros(2*N,1);
       USV4.Uk_1 = zeros(2*N,1);
       USV1.Xek = zeros(6,1);
       USV2.Xek = zeros(6,1);
       USV3.Xek = zeros(6,1);
       USV4.Xek = zeros(6,1);
       [ Q_bar, R_bar ] = calWeight( Q,R,F,N ); % 计算权值矩阵，只计算一次
   end
   [USV1.Uc_k,J1,USV1.Xhat]=LMPCnavi1(USV1.X,Xrm1,Q_bar,R_bar,taucrm1,N,tm,USV1.Nm);%Xhat为预测的轨迹
   [USV2.Uc_k,J2,USV2.Xhat]=LMPCnavi2(USV2.X,Xrm2,Q_bar,R_bar,taucrm2,N,tm,USV2.Nm);%Xhat为预测的轨迹
   [USV3.Uc_k,J3,USV3.Xhat]=LMPCnavi3(USV3.X,Xrm3,Q_bar,R_bar,taucrm3,N,tm,USV3.Nm);%Xhat为预测的轨迹
   [USV4.Uc_k,J4,USV4.Xhat]=LMPCnavi4(USV4.X,Xrm4,Q_bar,R_bar,taucrm4,N,tm,USV4.Nm);%Xhat为预测的轨迹
   USV1.Uk_1 = USV1.Uc_k;
   USV1.tauc = USV1.Uc_k(1:2);
   USV2.Uk_1 = USV2.Uc_k;
   USV2.tauc = USV2.Uc_k(1:2);
   USV3.Uk_1 = USV3.Uc_k;
   USV3.tauc = USV3.Uc_k(1:2);
   USV4.Uk_1 = USV4.Uc_k;
   USV4.tauc = USV4.Uc_k(1:2);
   
   % animation
   for i = 1:N+1
       USV1.xhat(i,1) = USV1.Xhat(6*(i-1)+4);
       USV1.yhat(i,1) = USV1.Xhat(6*(i-1)+5);
       USV2.xhat(i,1) = USV2.Xhat(6*(i-1)+4);
       USV2.yhat(i,1) = USV2.Xhat(6*(i-1)+5);
       USV3.xhat(i,1) = USV3.Xhat(6*(i-1)+4);
       USV3.yhat(i,1) = USV3.Xhat(6*(i-1)+5);
       USV4.xhat(i,1) = USV4.Xhat(6*(i-1)+4);
       USV4.yhat(i,1) = USV4.Xhat(6*(i-1)+5);
   end
   hn = figure(1);
   hold off
   pos1 = [USV1.X(4),USV1.X(5)]'; psi1 = USV1.X(6);
   pos2 = [USV2.X(4),USV2.X(5)]'; psi2 = USV2.X(6);
   pos3 = [USV3.X(4),USV3.X(5)]'; psi3 = USV3.X(6);
   pos4 = [USV4.X(4),USV4.X(5)]'; psi4 = USV4.X(6);
   linewid = 1;
   % USVs
   modelplot( pos1, psi1, 'r',linewid ); hold on
   modelplot( pos2, psi2, 'g',linewid );
   modelplot( pos3, psi3, 'b',linewid );
   modelplot( pos4, psi4, 'c',linewid );
   if t==0
       pos1Data = []; pos2Data = []; pos3Data = []; pos4Data = [];
   end
   pos1Data(k,:) = [USV1.X(4),USV1.X(5)];
   pos2Data(k,:) = [USV2.X(4),USV2.X(5)];
   pos3Data(k,:) = [USV3.X(4),USV3.X(5)];
   pos4Data(k,:) = [USV4.X(4),USV4.X(5)];
   % trajectories of USVs
   plot(pos1Data(:,2),pos1Data(:,1),'r-');
   plot(pos2Data(:,2),pos2Data(:,1),'g-'); 
   plot(pos3Data(:,2),pos3Data(:,1),'b-');
   plot(pos4Data(:,2),pos4Data(:,1),'c-');
   % desired paths and trajectories of reference points
   plot(xd1(:,5),xd1(:,4),'r--',xd2(:,5),xd2(:,4),'g--',xd3(:,5),xd3(:,4),'b--',...
     xd4(:,5),xd4(:,4),'c--',pr5(:,2),pr5(:,1),'k-',pr6(:,2),pr6(:,1),'k-');
   % pridicted trajectories of USVs
   plot(USV1.yhat,USV1.xhat,'m-','linewid',2);
   plot(USV2.yhat,USV2.xhat,'m-','linewid',2);
   plot(USV3.yhat,USV3.xhat,'m-','linewid',2);
   plot(USV4.yhat,USV4.xhat,'m-','linewid',2);
   % reference points
   h=rectangle('Position',[xd1(k,5)-2,xd1(k,4)-2,2*2,2*2],'Curvature',[1,1],'EdgeColor','k');
   set(h,'LineStyle','-','linewid',1);
   h=rectangle('Position',[xd2(k,5)-2,xd2(k,4)-2,2*2,2*2],'Curvature',[1,1],'EdgeColor','k');
   set(h,'LineStyle','-','linewid',1);
   h=rectangle('Position',[xd3(k,5)-2,xd3(k,4)-2,2*2,2*2],'Curvature',[1,1],'EdgeColor','k');
   set(h,'LineStyle','-','linewid',1);
   h=rectangle('Position',[xd4(k,5)-2,xd4(k,4)-2,2*2,2*2],'Curvature',[1,1],'EdgeColor','k');
   set(h,'LineStyle','-','linewid',1);
   grid on;  
   xlabel('y / m');
   ylabel('x / m');
   title('map');
   Xmin = -10; Xmax = 50;
   Ymin = 0; Ymax = 60;
   area = [Xmin Xmax Ymin Ymax];
   axis(area);
   drawnow;
   % save as gif
   F=getframe(gcf);
   I=frame2im(F);
   [I,map]=rgb2ind(I,256);
   if k == 1
       imwrite(I,map,'test.gif','gif', 'Loopcount',inf,'DelayTime',0.2);
   elseif rem(k,50)==0
       imwrite(I,map,'test.gif','gif','WriteMode','append','DelayTime',0.2);
   end
   pause(0.0001);
   
%    figure(2)
%    hold off
%    if t==0
%        tData = []; uu1 = []; uu2 = []; uu3 = []; uu4 = [];
%    end
%    tData(k,1) = t;
%    uu1(k,1) = USV1.X(1); uu2(k,1) = USV2.X(1);
%    uu3(k,1) = USV3.X(1); uu4(k,1) = USV4.X(1);
%    plot(tData,uu1,'r-'); hold on
%    plot(tData,uu2,'g-'); hold on
%    plot(tData,uu3,'b-'); hold on
%    plot(tData,uu4,'c-'); hold on
%    xlabel('t (s)');
%    ylabel('u (m/s)');
%    title('map');
%    Xmin = 0; Xmax = 100;
%    Ymin = 0; Ymax = 1;
%    area = [Xmin Xmax Ymin Ymax];
%    axis(area);
%    drawnow;
%    pause(0.001);
   % sim data
   %--------------------------USV1--------------------------
   if t==0
       KKKK1 = 1;
       Xdk1 = Xrm1(6*(KKKK1-1)+1:6*(KKKK1-1)+6*(N+1),1);
   end
   if USV1.Nm==0
   if 6*(KKKK1-1)+6*(N+1)>Nrm1
      if 6*(KKKK1-1)+6<=Nrm1
          Xdk1(1:Nrm1-6*(KKKK1-1),1) = Xrm1(6*(KKKK1-1)+1:Nrm1);
          Xdk1(Nrm1-6*(KKKK1-1)+1:6*(N+1),1) = repmat(Xrm1(Nrm1-5:Nrm1,1),(N+1)-(Nrm1/6-(KKKK1-1)),1);
      else
          Xdk1(1:6*(N+1),1) = repmat(Xrm1(Nrm1-5:Nrm1,1),N+1,1);
      end
   else
    Xdk1 = Xrm1(6*(KKKK1-1)+1:6*(KKKK1-1)+6*(N+1),1);
   end 
       KKKK1 = KKKK1+1;
   end
   Xrk1 =Xdk1(1:6,1);

   USV1.Xek = USV1.Xm-Xrk1;
   de1 = norm(USV1.Xek(4:5));
   
   %--------------------------USV2--------------------------
   if t==0
       KKKK2 = 1;
       Xdk2 = Xrm2(6*(KKKK2-1)+1:6*(KKKK2-1)+6*(N+1),1);
   end
   if USV2.Nm==0
   if 6*(KKKK2-1)+6*(N+1)>Nrm2
      if 6*(KKKK2-1)+6<=Nrm2
          Xdk2(1:Nrm2-6*(KKKK2-1),1) = Xrm2(6*(KKKK2-1)+1:Nrm2);
          Xdk2(Nrm2-6*(KKKK2-1)+1:6*(N+1),1) = repmat(Xrm2(Nrm2-5:Nrm2,1),(N+1)-(Nrm2/6-(KKKK2-1)),1);
      else
          Xdk2(1:6*(N+1),1) = repmat(Xrm2(Nrm2-5:Nrm2,1),N+1,1);
      end
   else
    Xdk2 = Xrm2(6*(KKKK2-1)+1:6*(KKKK2-1)+6*(N+1),1);
   end 
       KKKK2 = KKKK2+1;
   end
   Xrk2 =Xdk2(1:6,1);

   USV2.Xek = USV2.Xm-Xrk2;
   de2 = norm(USV2.Xek(4:5));
   %--------------------------USV3--------------------------
   if t==0
       KKKK3 = 1;
       Xdk3 = Xrm3(6*(KKKK3-1)+1:6*(KKKK3-1)+6*(N+1),1);
   end
   if USV3.Nm==0
   if 6*(KKKK3-1)+6*(N+1)>Nrm3
      if 6*(KKKK3-1)+6<=Nrm3
          Xdk3(1:Nrm3-6*(KKKK3-1),1) = Xrm3(6*(KKKK3-1)+1:Nrm3);
          Xdk3(Nrm3-6*(KKKK3-1)+1:6*(N+1),1) = repmat(Xrm3(Nrm3-5:Nrm3,1),(N+1)-(Nrm3/6-(KKKK3-1)),1);
      else
          Xdk3(1:6*(N+1),1) = repmat(Xrm3(Nrm3-5:Nrm3,1),N+1,1);
      end
   else
    Xdk3 = Xrm3(6*(KKKK3-1)+1:6*(KKKK3-1)+6*(N+1),1);
   end 
       KKKK3 = KKKK3+1;
   end
   Xrk3 =Xdk3(1:6,1);

   USV3.Xek = USV3.Xm-Xrk3;
   de3 = norm(USV3.Xek(4:5));
   %--------------------------USV4--------------------------
   if t==0
       KKKK4 = 1;
       Xdk4 = Xrm4(6*(KKKK4-1)+1:6*(KKKK4-1)+6*(N+1),1);
   end
   if USV4.Nm==0
   if 6*(KKKK4-1)+6*(N+1)>Nrm4
      if 6*(KKKK4-1)+6<=Nrm4
          Xdk4(1:Nrm4-6*(KKKK4-1),1) = Xrm4(6*(KKKK4-1)+1:Nrm4);
          Xdk4(Nrm4-6*(KKKK4-1)+1:6*(N+1),1) = repmat(Xrm4(Nrm4-5:Nrm4,1),(N+1)-(Nrm4/6-(KKKK4-1)),1);
      else
          Xdk4(1:6*(N+1),1) = repmat(Xrm4(Nrm4-5:Nrm4,1),N+1,1);
      end
   else
    Xdk4 = Xrm4(6*(KKKK4-1)+1:6*(KKKK4-1)+6*(N+1),1);
   end 
       KKKK4 = KKKK4+1;
   end
   Xrk4 =Xdk4(1:6,1);

   USV4.Xek = USV4.Xm-Xrk4;
   de4 = norm(USV4.Xek(4:5));

   xout(k,:)=[t,J1,J2,J3,J4,de1,de2,de3,de4,USV1.X',USV2.X',USV3.X',USV4.X',...
              USV1.tauc',USV2.tauc',USV3.tauc',USV4.tauc']; 
   
   USV1.xhato(k,:) = USV1.xhat';% 每一时刻预测轨迹数据存储
   USV1.yhato(k,:) = USV1.yhat';% 每一时刻预测轨迹数据存储
   USV2.xhato(k,:) = USV2.xhat';% 每一时刻预测轨迹数据存储
   USV2.yhato(k,:) = USV2.yhat';% 每一时刻预测轨迹数据存储
   USV3.xhato(k,:) = USV3.xhat';% 每一时刻预测轨迹数据存储
   USV3.yhato(k,:) = USV3.yhat';% 每一时刻预测轨迹数据存储
   USV4.xhato(k,:) = USV4.xhat';% 每一时刻预测轨迹数据存储
   USV4.yhato(k,:) = USV4.yhat';% 每一时刻预测轨迹数据存储
end
t = xout(:,1);
J1 = xout(:,2);
J2 = xout(:,3);
J3 = xout(:,4);
J4 = xout(:,5);
de1 = xout(:,6);
de2 = xout(:,7);
de3 = xout(:,8);
de4 = xout(:,9);
x1 = xout(:,10:15);
x2 = xout(:,16:21);
x3 = xout(:,22:27);
x4 = xout(:,28:33);
tauc1 = xout(:,34:35);
tauc2 = xout(:,36:37);
tauc3 = xout(:,38:39);
tauc4 = xout(:,40:41);

%% 绘图
figure(3)
plot(x1(:,5),x1(:,4),'r-',x2(:,5),x2(:,4),'g-',x3(:,5),x3(:,4),'b-',x4(:,5),x4(:,4),'c-',...
     xd1(:,5),xd1(:,4),'r--',xd2(:,5),xd2(:,4),'g--',xd3(:,5),xd3(:,4),'b--',...
     xd4(:,5),xd4(:,4),'c--',pr5(:,2),pr5(:,1),'k-',pr6(:,2),pr6(:,1),'k-');hold on
legend('USV1','USV2','USV3','USV4','Td1','Td2','Td3','Td4','pr5','pr6');
% plot(yrm,xrm,'ro');

figure(4)
plot(t,J1,'r-',t,J2,'g-',t,J3,'b-',t,J4','c-');

figure(5)
plot(t,tauc1(:,1),'r-',t,tauc1(:,2),'b-');

figure(6)
plot(t,x1(:,1),'r-',t,x2(:,1),'g-',t,x3(:,1),'b-',t,x4(:,1),'c-');

figure(7)
plot(t,xd1(:,1),'r-',t,xd2(:,1),'g-',t,xd3(:,1),'b-',t,xd4(:,1),'c-');












