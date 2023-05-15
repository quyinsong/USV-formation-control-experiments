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
tfinal = 600;
Ns = tfinal/ts;

%% parameters initializition
% USV states
USV1.x0 = [0 0 0 12 -3 60*pi/180]';
USV2.x0 = [0 0 0 6 -2 90*pi/180]';
USV3.x0 = [0 0 0 0 -5 80*pi/180]';
USV4.x0 = [0 0 0  -6 -5 70*pi/180]';

% USV1.x0 = [0 0 0 3 -1 60*pi/180]';
% USV2.x0 = [0 0 0 2 -2 90*pi/180]';
% USV3.x0 = [0 0 0 1 -5 80*pi/180]';
% USV4.x0 = [0 0 0  5 -1 70*pi/180]';
% USV inputs
USV1.tauc_etm=[0 0]';
USV2.tauc_etm=[0 0]';
USV3.tauc_etm=[0 0]';
USV4.tauc_etm= [0 0]';
USV1.tauc=[0 0]';
USV2.tauc=[0 0]';
USV3.tauc=[0 0]';
USV4.tauc= [0 0]';
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
DL = diag(sum(ALL,2));
LL = DL-ALL;
HL = LL+diag(ALS);
HL = [HL,-ALS];
I51 = [1;1;1;1;1];

for k=1:Ns
   % time series
   t = (k-1)*ts;
   
   % path info
   % pc = 1: straight line
   % pc = 2: Wang 2021
   % pc = 3: circle
   % pc = 4: straight line + circle
   if t == 0
      pc = 2;
      ud = 0.3; 
      delta0 = 0.1;
      eF1 = [0;0]; eF2 = [0;0]; eF3 = [0;0]; eF4 = [0;0]; 
      eL1 = 0; eL2 = 0; 
      pi1 = [USV1.x0(4);USV1.x0(5)]; pi2 = [USV2.x0(4);USV2.x0(5)]; 
      pi3 = [USV3.x0(4);USV3.x0(5)]; pi4 = [USV4.x0(4);USV4.x0(5)]; 
      pi10 = pi1; pi20 = pi2; pi30 = pi3; pi40 = pi4; 
      pi1_dot = [0;0]; pi2_dot = [0;0]; pi3_dot = [0;0]; pi4_dot = [0;0];
      eF1_s = 0; eF2_s = 0; eF3_s = 0; eF4_s = 0;
      
      O1 = [35 30]'; O2 = [60 60]'; O3 = [70 85]'; O4 = [85 85]';
      RO1 = 2; RO2 = 2; RO3 = 2; RO4 = 2;
      R1_up = RO1+3;  R2_up = RO2+3; R3_up = RO3+3; R4_up = RO3+3; 
      R1_dn = RO1+0.1;  R2_dn = RO2+0.1; R3_dn = RO3+0.1; R4_dn = RO3+0.1;
      
      Rusv_up = 3; Rusv_dn = 1;
   end
%    tau_w = [0.08*sin(0.2*t)*cos(0.1*t)+0.2,...
%             0.05*sin(0.2*t)*cos(0.1*t),...
%             0.08*sin(0.2*t)*cos(0.1*t)]';
        
   % 一阶马尔科夫过程
   if t==0
       wk_u = 0;
       wk_v = 0;
       wk_r = 0;
   end
   miu_u = 4; miu_v = 2; miu_r = 3;
   wu = normrnd(0,0.5); wv = normrnd(0,0.5); wr = normrnd(0,0.5); 
   wk_u_dot = -wk_u+miu_u*wu;
   wk_v_dot = -wk_v+miu_v*wv;
   wk_r_dot = -wk_r+miu_r*wr;
   wk_u = wk_u_dot*ts+wk_u;
   wk_v = wk_v_dot*ts+wk_v;
   wk_r = wk_r_dot*ts+wk_r;  
%    tau_w = [wk_u,wk_v,wk_r]';
   
   tau_w = [0.1*sin(0.1*t)+0.3*sin(0.05*t) 0.05*sin(0.05*t) ...
            0.05*sin(0.1*t)+0.05*sin(0.05*t)]';
%    tau_w = [0 0 0]';

   % USV states update
%    [USV1.x,USV1.tau,USV1.f] = CS1( USV1.x0, USV1.tauc_etm, tau_w, ts);
%    [USV2.x,USV2.tau,USV2.f] = CS2( USV2.x0, USV2.tauc_etm, tau_w, ts);
%    [USV3.x,USV3.tau,USV3.f] = CS3( USV3.x0, USV3.tauc_etm, tau_w, ts);
%    [USV4.x,USV4.tau,USV4.f] = CS4( USV4.x0, USV4.tauc_etm, tau_w, ts);
   
   [USV1.x,USV1.tau,USV1.f, USV1.zeta, USV1.emi, USV1.Gamma] = CS1tf( USV1.x0, USV1.tauc_etm, tau_w, ts );
   [USV2.x,USV2.tau,USV2.f, USV2.zeta, USV2.emi, USV2.Gamma] = CS2tf( USV2.x0, USV2.tauc_etm, tau_w, ts );
   [USV3.x,USV3.tau,USV3.f, USV3.zeta, USV3.emi, USV3.Gamma] = CS3tf( USV3.x0, USV3.tauc_etm, tau_w, ts );
   [USV4.x,USV4.tau,USV4.f, USV4.zeta, USV4.emi, USV4.Gamma] = CS4tf( USV4.x0, USV4.tauc_etm, tau_w, ts );
   
   USV1.zeta_bar = [USV1.zeta(1),USV1.zeta(3)]';
   % USV states and Coordinate transformation
   USV1.eta = [USV1.x(4) USV1.x(5) USV1.x(6)]';
   USV2.eta = [USV2.x(4) USV2.x(5) USV2.x(6)]';
   USV3.eta = [USV3.x(4) USV3.x(5) USV3.x(6)]';
   USV4.eta = [USV4.x(4) USV4.x(5) USV4.x(6)]';
   
   USV1.p = [USV1.eta(1);USV1.eta(2)];
   USV2.p = [USV2.eta(1);USV2.eta(2)];
   USV3.p = [USV3.eta(1);USV3.eta(2)];
   USV4.p = [USV4.eta(1);USV4.eta(2)];
    
   USV1.nu = [USV1.x(1) USV1.x(2) USV1.x(3)]';USV1.nu_bar = [USV1.x(1) USV1.x(3)]';
   USV2.nu = [USV2.x(1) USV2.x(2) USV2.x(3)]';USV2.nu_bar = [USV2.x(1) USV2.x(3)]';
   USV3.nu = [USV3.x(1) USV3.x(2) USV3.x(3)]';USV3.nu_bar = [USV3.x(1) USV3.x(3)]';
   USV4.nu = [USV4.x(1) USV4.x(2) USV4.x(3)]';USV4.nu_bar = [USV4.x(1) USV4.x(3)]';
   
   USV1.v_bar = USV1.nu(2)+USV1.emi*USV1.nu(3);
   USV1.nu_tf = [USV1.x(1),USV1.v_bar,USV1.x(3)]';
   USV2.v_bar = USV2.nu(2)+USV2.emi*USV2.nu(3);
   USV2.nu_tf = [USV2.x(1),USV2.v_bar,USV2.x(3)]';
   USV3.v_bar = USV3.nu(2)+USV3.emi*USV3.nu(3);
   USV3.nu_tf = [USV3.x(1),USV3.v_bar,USV3.x(3)]';
   USV4.v_bar = USV4.nu(2)+USV4.emi*USV4.nu(3);
   USV4.nu_tf = [USV4.x(1),USV4.v_bar,USV4.x(3)]';
   
   [ USV1.x_bar, USV1.y_bar] = CoordinateTrans( USV1.eta(1),USV1.eta(2),USV1.eta(3),USV1.emi );
   [ USV2.x_bar, USV2.y_bar] = CoordinateTrans( USV2.eta(1),USV2.eta(2),USV2.eta(3),USV1.emi );
   [ USV3.x_bar, USV3.y_bar] = CoordinateTrans( USV3.eta(1),USV3.eta(2),USV3.eta(3),USV1.emi );
   [ USV4.x_bar, USV4.y_bar] = CoordinateTrans( USV4.eta(1),USV4.eta(2),USV4.eta(3),USV1.emi );
   
   USV1.p_tf = [ USV1.x_bar, USV1.y_bar]';
   USV2.p_tf = [ USV2.x_bar, USV2.y_bar]';
   USV3.p_tf = [ USV3.x_bar, USV3.y_bar]';
   USV4.p_tf = [ USV4.x_bar, USV4.y_bar]';
   
   USV1.eta_tf = [ USV1.x_bar, USV1.y_bar,USV1.x(6)]';
   USV2.eta_tf = [ USV2.x_bar, USV2.y_bar,USV2.x(6)]';
   USV3.eta_tf = [ USV3.x_bar, USV3.y_bar,USV3.x(6)]';
   USV4.eta_tf = [ USV4.x_bar, USV4.y_bar,USV4.x(6)]';
   
   R_psi1 = [cos(USV1.eta(3)) -sin(USV1.eta(3));
             sin(USV1.eta(3)) cos(USV1.eta(3))];
   R_psi2 = [cos(USV2.eta(3)) -sin(USV2.eta(3));
             sin(USV2.eta(3)) cos(USV2.eta(3))];
   R_psi3 = [cos(USV3.eta(3)) -sin(USV3.eta(3));
             sin(USV3.eta(3)) cos(USV3.eta(3))];
   R_psi4 = [cos(USV4.eta(3)) -sin(USV4.eta(3));
             sin(USV4.eta(3)) cos(USV4.eta(3))];
   R_psi = [R_psi1',zeros(2,6);
            zeros(2,2),R_psi2',zeros(2,4);
            zeros(2,4),R_psi3',zeros(2,2);
            zeros(2,6),R_psi4'];
        
   % super leader and virtual leader
   [pr0,pr0_dw,theta0,vs] = superLeader0( pc, ud, ts );

   [pr1,pr1_dw,theta1] = virtualLeader1( pc, eF1_s, eL1, vs, ts );
   [pr2,pr2_dw,theta2] = virtualLeader2( pc, eF2_s, eL2, vs, ts );
   
%    eF1_s = ALF(1,1)*pr1_dw'*eF1+ALF(1,2)*pr1_dw'*eF2+ALF(1,3)*pr1_dw'*eF3+ALF(1,4)*pr1_dw'*eF4;
%    eF2_s = ALF(2,1)*pr2_dw'*eF1+ALF(2,2)*pr2_dw'*eF2+ALF(2,3)*pr2_dw'*eF3+ALF(2,4)*pr2_dw'*eF4;
   
   eF1_s = ALF(1,1)*pr1_dw'*eF1;
   eF2_s = ALF(2,4)*pr2_dw'*eF4;
   % obstacle detection
   OP = [O1,O2,O3,O4];
%    OP = [O1];
   RO_up = [R1_up,R2_up,R3_up,R4_up]';
   RO_dn = [R1_dn,R2_dn,R3_dn,R4_dn]';
   USVP = [pi1,pi2,pi3,pi4];
   [od1,od2,od3,od4] = odetection(USVP,OP,RO_up,RO_dn);
   
   % collision detection
   [usvd1,usvd2,usvd3,usvd4] = usvdetection(USVP,Rusv_up);
   
   % obstacle avoidance
   if ~isempty(od1)
       pi_o1 = APF_O(pi1,od1);
   else
       pi_o1 = [0;0];
   end
   if ~isempty(od2)
       pi_o2 = APF_O(pi2,od2);
   else
       pi_o2 = [0;0];
   end
   if ~isempty(od3)
       pi_o3 = APF_O(pi3,od3);
   else
       pi_o3 = [0;0];
   end
   if ~isempty(od4)
       pi_o4 = APF_O(pi4,od4);
   else
       pi_o4 = [0;0];
   end
   
   % collision avoidance
   if ~isempty(usvd1)
       pi_usv1 = APF_usv(pi1,usvd1,Rusv_up,Rusv_dn);
   else 
       pi_usv1 = [0;0];
   end
   if ~isempty(usvd2)
       pi_usv2 = APF_usv(pi2,usvd2,Rusv_up,Rusv_dn);
   else 
       pi_usv2 = [0;0];
   end
   if ~isempty(usvd3)
       pi_usv3 = APF_usv(pi3,usvd3,Rusv_up,Rusv_dn);
   else 
       pi_usv3 = [0;0];
   end
   if ~isempty(usvd4)
       pi_usv4 = APF_usv(pi4,usvd4,Rusv_up,Rusv_dn);
   else 
       pi_usv4 = [0;0];
   end
   
   % containment Motion Generator (CMG)
   Theta = [theta1,theta2,theta0]';
   EF = [AFF(1,2)*(pi1-pi2)+AFF(1,3)*(pi1-pi3)+AFF(1,4)*(pi1-pi4)+...
        AFL(1,1)*(pi1-pr1)+AFL(1,2)*(pi1-pr2);
        AFF(2,1)*(pi2-pi1)+AFF(2,3)*(pi2-pi3)+AFF(2,4)*(pi2-pi4)+...
        AFL(2,1)*(pi2-pr1)+AFL(2,2)*(pi2-pr2);
        AFF(3,1)*(pi3-pi1)+AFF(3,2)*(pi3-pi2)+AFF(3,4)*(pi3-pi4)+...
        AFL(3,1)*(pi3-pr1)+AFL(3,2)*(pi3-pr2);
        AFF(4,1)*(pi4-pi1)+AFF(4,2)*(pi4-pi2)+AFF(4,3)*(pi4-pi3)+...
        AFL(4,1)*(pi4-pr1)+AFL(4,2)*(pi4-pr2)];
    
   EL = HL*Theta+sum(AL.*Pijd,2);
   pj1_dot = AFF(1,1)*pi1_dot+AFF(1,2)*pi2_dot+AFF(1,3)*pi3_dot+AFF(1,4)*pi4_dot;
   pj2_dot = AFF(2,1)*pi1_dot+AFF(2,2)*pi2_dot+AFF(2,3)*pi3_dot+AFF(2,4)*pi4_dot;
   pj3_dot = AFF(3,1)*pi1_dot+AFF(3,2)*pi2_dot+AFF(3,3)*pi3_dot+AFF(3,4)*pi4_dot;
   pj4_dot = AFF(4,1)*pi1_dot+AFF(4,2)*pi2_dot+AFF(4,3)*pi3_dot+AFF(4,4)*pi4_dot;

   di1 = sum(A(1,:),2);
   di2 = sum(A(2,:),2);
   di3 = sum(A(3,:),2);
   di4 = sum(A(4,:),2);
   
   eF1 = EF(1:2); eF2 = EF(3:4); eF3 = EF(5:6); eF4 = EF(7:8);
   eL1 = EL(1); eL2 = EL(2);
   
   vsi1 = (AFL(1,1)*pr1_dw+AFL(1,2)*pr2_dw)*vs;
   vsi2 = (AFL(2,1)*pr1_dw+AFL(2,2)*pr2_dw)*vs;
   vsi3 = (AFL(3,1)*pr1_dw+AFL(3,2)*pr2_dw)*vs;
   vsi4 = (AFL(4,1)*pr1_dw+AFL(4,2)*pr2_dw)*vs;
   
%    pi_o1 = [0 0]'; pi_o2 = [0 0]'; pi_o3 = [0 0]'; pi_o4 = [0 0]';
%    pi_usv1 = [0 0]'; pi_usv2 = [0 0]'; pi_usv3 = [0 0]'; pi_usv4 = [0 0]';
   
   [ pi1,pi1_dot] = CMG1( eF1,pj1_dot,vsi1,di1, pi10, pi_o1, pi_usv1, ts );
   [ pi2,pi2_dot] = CMG2( eF2,pj2_dot,vsi2,di2, pi20, pi_o2, pi_usv2, ts );
   [ pi3,pi3_dot] = CMG3( eF3,pj3_dot,vsi3,di3, pi30, pi_o3, pi_usv3, ts );
   [ pi4,pi4_dot] = CMG4( eF4,pj4_dot,vsi4,di4, pi40, pi_o4, pi_usv4, ts );
   
   % test Guidance
%    if t==0
%        x00 = 0;
%        y00 = 0;
%    end
%    x00_dot = 0.2;
%    y00_dot = 0.2;
%    x00 = ts*x00_dot+x00;
%    y00 = ts*y00_dot+y00;
%    d = [delta0,0]';
%    Pi = [x00;y00;x00;y00;x00;y00;x00;y00];
%    z1 = R_psi*(P-Pi);
%    dd = Kerector(I41,d);
%    z2 = z1+dd;
%    hi = diag([1,delta0]);
%    pi1_dot = [x00_dot,y00_dot]';
%    pi2_dot = [x00_dot,y00_dot]';
%    pi3_dot = [x00_dot,y00_dot]';
%    pi4_dot = [x00_dot,y00_dot]';

   % tracking errors
   P_tf = [USV1.p_tf;USV2.p_tf;USV3.p_tf;USV4.p_tf];
   d = [delta0,0]';
   Pi = [pi1;pi2;pi3;pi4];
   z1 = R_psi*(P_tf-Pi);
   dd = Kerector(I41,d);
   z2 = z1+dd;
   hi = diag([1,delta0]);

   % ESO
   USV1.nu_hat = ESO1(USV1.eta_tf,ts);
   USV2.nu_hat = ESO2(USV2.eta_tf,ts);
   USV3.nu_hat = ESO3(USV3.eta_tf,ts);
   USV4.nu_hat = ESO4(USV4.eta_tf,ts);
   USV1.u_hat = USV1.nu_hat(1); USV2.u_hat = USV2.nu_hat(1); USV3.u_hat = USV3.nu_hat(1); USV4.u_hat = USV4.nu_hat(1);
   USV1.v_bar_hat = USV1.nu_hat(2); USV2.v_bar_hat = USV2.nu_hat(2); USV3.v_bar_hat = USV3.nu_hat(2); USV4.v_bar_hat = USV4.nu_hat(2);
   
   % Guidance
   
   USV1.nu_c = Guidance1( hi, R_psi1, z2(1:2), pi1_dot, USV1.nu_hat(1:2));
   USV2.nu_c = Guidance2( hi, R_psi2, z2(3:4), pi2_dot, USV2.nu_hat(1:2));
   USV3.nu_c = Guidance3( hi, R_psi3, z2(5:6), pi3_dot, USV3.nu_hat(1:2));
   USV4.nu_c = Guidance4( hi, R_psi4, z2(7:8), pi4_dot, USV4.nu_hat(1:2));
 
   % Control
%    [USV1.tauc,USV1.e_bar,USV1.zetahat] = ctr1( [USV1.u_hat,USV1.nu_bar(2)]', USV1.nu_c, USV1.tau, USV1.tauc, USV1.Gamma, ts );
%    [USV2.tauc,USV2.e_bar,USV2.zetahat] = ctr2( [USV2.u_hat,USV2.nu_bar(2)]', USV2.nu_c, USV2.tau, USV2.tauc, USV2.Gamma, ts );
%    [USV3.tauc,USV3.e_bar,USV3.zetahat] = ctr3( [USV3.u_hat,USV3.nu_bar(2)]', USV3.nu_c, USV3.tau, USV3.tauc, USV3.Gamma, ts );
%    [USV4.tauc,USV4.e_bar,USV4.zetahat] = ctr4( [USV4.u_hat,USV4.nu_bar(2)]', USV4.nu_c, USV4.tau, USV4.tauc, USV4.Gamma, ts );

%    USV1.nu_tf_hat = [USV1.u_hat,USV1.v_bar_hat,USV1.nu_bar(2)]';
%    USV2.nu_tf_hat = [USV2.u_hat,USV2.v_bar_hat,USV2.nu_bar(2)]';
%    USV3.nu_tf_hat = [USV3.u_hat,USV3.v_bar_hat,USV3.nu_bar(2)]';
%    USV4.nu_tf_hat = [USV4.u_hat,USV4.v_bar_hat,USV4.nu_bar(2)]';
   
   [USV1.tauc,USV1.e_bar,USV1.zetahat] = nnctr1( USV1.nu_hat, USV1.nu_c, USV1.tau, USV1.tauc, USV1.Gamma, ts );
   [USV2.tauc,USV2.e_bar,USV2.zetahat] = nnctr2( USV2.nu_hat, USV2.nu_c, USV2.tau, USV2.tauc, USV2.Gamma, ts );
   [USV3.tauc,USV3.e_bar,USV3.zetahat] = nnctr3( USV3.nu_hat, USV3.nu_c, USV3.tau, USV3.tauc, USV3.Gamma, ts );
   [USV4.tauc,USV4.e_bar,USV4.zetahat] = nnctr4( USV4.nu_hat, USV4.nu_c, USV4.tau, USV4.tauc, USV4.Gamma, ts );

   % ETM
   if t==0
        emiu = 0.3;
        emir = 0.3;
        du = 0.05; dr = 0.01;
   end
   [USV1.tauc_etm,Tk1] = ETM1(USV1.tauc,emiu,emir,du,dr);
   [USV2.tauc_etm,Tk2] = ETM2(USV2.tauc,emiu,emir,du,dr);
   [USV3.tauc_etm,Tk3] = ETM3(USV3.tauc,emiu,emir,du,dr);
   [USV4.tauc_etm,Tk4] = ETM4(USV4.tauc,emiu,emir,du,dr);

   z2_norm = [norm(z2(1:2)),norm(z2(3:4)),norm(z2(5:6)),norm(z2(7:8))]';
   eF1_norm = norm(eF1); eF2_norm = norm(eF2); eF3_norm = norm(eF3); eF4_norm = norm(eF4);
   
   xout(k,:) = [t,USV1.x', USV2.x', USV3.x', USV4.x', USV1.tau' USV2.tau' USV3.tau' USV4.tau',...
                USV1.tauc' USV2.tauc' USV3.tauc' USV4.tauc',z2_norm',wk_u,wk_v,wk_r,pi1',pi2',pi3',pi4',pr1',pr2',...
                Tk1',Tk2',Tk3',Tk4',eF1_norm,eF2_norm,eF3_norm,eF4_norm,eL1,eL2,USV1.u_hat,USV1.v_bar_hat,...
                USV2.u_hat,USV2.v_bar_hat,USV3.u_hat,USV3.v_bar_hat,USV4.u_hat,USV4.v_bar_hat,USV1.v_bar,...
                USV2.v_bar,USV3.v_bar,USV4.v_bar,USV1.zeta_bar',USV1.zetahat',USV1.e_bar',USV2.e_bar',USV3.e_bar',USV4.e_bar'];
    
    
end
%% simulation data
t = xout(:,1);
USV1.x = xout(:,2:7);
USV2.x = xout(:,8:13);
USV3.x = xout(:,14:19);
USV4.x = xout(:,20:25);
USV1.tau = xout(:,26:27);
USV2.tau = xout(:,28:29);
USV3.tau = xout(:,30:31);
USV4.tau = xout(:,32:33);
USV1.tauc = xout(:,34:35);
USV2.tauc = xout(:,36:37);
USV3.tauc = xout(:,38:39);
USV4.tauc = xout(:,40:41);
z1_bar_norm = xout(:,42);
z2_bar_norm = xout(:,43);
z3_bar_norm = xout(:,44);
z4_bar_norm = xout(:,45);
wk_u = xout(:,46);
wk_v = xout(:,47);
wk_r = xout(:,48);
pi1 = xout(:,49:50);
pi2 = xout(:,51:52);
pi3 = xout(:,53:54);
pi4 = xout(:,55:56);
pr1 = xout(:,57:58);
pr2 = xout(:,59:60);
Tk1 = xout(:,61:62);
Tk2 = xout(:,63:64);
Tk3 = xout(:,65:66);
Tk4 = xout(:,67:68);
eF1_norm = xout(:,69);
eF2_norm = xout(:,70);
eF3_norm = xout(:,71);
eF4_norm = xout(:,72);
eL1 = xout(:,73);
eL2 = xout(:,74);
USV1.uhat = xout(:,75);
USV1.v_bar_hat = xout(:,76);
USV2.uhat = xout(:,77);
USV2.v_bar_hat = xout(:,78);
USV3.uhat = xout(:,79);
USV3.v_bar_hat = xout(:,80);
USV4.uhat = xout(:,81);
USV4.v_bar_hat = xout(:,82);
USV1.v_bar = xout(:,83);
USV2.v_bar = xout(:,84);
USV3.v_bar = xout(:,85);
USV4.v_bar = xout(:,86);
USV1.zeta_bar = xout(:,87:88);
USV1.zetahat = xout(:,89:90);
USV1.e_bar = xout(:,91:92);
USV2.e_bar = xout(:,93:94);
USV3.e_bar = xout(:,95:96);
USV4.e_bar = xout(:,97:98);


%% PLOTS
close all
% formation
figure(1); hold on

linewid = 1;
fontsize = 10.5;
fontname = 'Time news roman';
fontweight = 'bold';

yrange=[-30 30 150]; xrange = [-30 30 150];

Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);
plt1 = plot(USV1.x(:,5),USV1.x(:,4),'r-'); 
plt2 = plot(USV2.x(:,5),USV2.x(:,4),'g-'); 
plt3 = plot(USV3.x(:,5),USV3.x(:,4),'b-'); 
plt4 = plot(USV4.x(:,5),USV4.x(:,4),'c-'); 

plt5 = plot(pi1(:,2),pi1(:,1),'r--'); 
plt6 = plot(pi2(:,2),pi2(:,1),'g--'); 
plt7 = plot(pi3(:,2),pi3(:,1),'b--'); 
plt8 = plot(pi4(:,2),pi4(:,1),'c--'); 

plt9 = plot(pr1(:,2),pr1(:,1),'r--'); 
plt10 = plot(pr2(:,2),pr2(:,1),'g--'); 

% plot obstacle
h=rectangle('Position',[O1(2)-R1_up,O1(1)-R1_up,2*R1_up,2*R1_up],'Curvature',[1,1],'EdgeColor','k');
set(h,'LineStyle','--','linewid',1);
h=rectangle('Position',[O2(2)-R2_up,O2(1)-R2_up,2*R2_up,2*R2_up],'Curvature',[1,1],'EdgeColor','k');
set(h,'LineStyle','--','linewid',1);
h=rectangle('Position',[O3(2)-R3_up,O3(1)-R3_up,2*R3_up,2*R3_up],'Curvature',[1,1],'EdgeColor','k');
set(h,'LineStyle','--','linewid',1);
h=rectangle('Position',[O4(2)-R4_up,O4(1)-R4_up,2*R4_up,2*R4_up],'Curvature',[1,1],'EdgeColor','k');
set(h,'LineStyle','--','linewid',1);

h=rectangle('Position',[O1(2)-R1_dn,O1(1)-R1_dn,2*R1_dn,2*R1_dn],'Curvature',[1,1],'EdgeColor','k');
set(h,'LineStyle','--','linewid',1);
h=rectangle('Position',[O2(2)-R2_dn,O2(1)-R2_dn,2*R2_dn,2*R2_dn],'Curvature',[1,1],'EdgeColor','k');
set(h,'LineStyle','--','linewid',1);
h=rectangle('Position',[O3(2)-R3_dn,O3(1)-R3_dn,2*R3_dn,2*R3_dn],'Curvature',[1,1],'EdgeColor','k');
set(h,'LineStyle','--','linewid',1);
h=rectangle('Position',[O4(2)-R4_dn,O4(1)-R4_dn,2*R4_dn,2*R4_dn],'Curvature',[1,1],'EdgeColor','k');
set(h,'LineStyle','--','linewid',1);

h=rectangle('Position',[O1(2)-RO1,O1(1)-RO1,2*RO1,2*RO1],'Curvature',[1,1],'EdgeColor','k');
set(h,'LineStyle','--','linewid',1, 'FaceColor',[0 .5 .5]);
h=rectangle('Position',[O2(2)-RO2,O2(1)-RO2,2*RO2,2*RO2],'Curvature',[1,1],'EdgeColor','k');
set(h,'LineStyle','--','linewid',1, 'FaceColor',[0 .5 .5]);
h=rectangle('Position',[O3(2)-RO3,O3(1)-RO3,2*RO3,2*RO3],'Curvature',[1,1],'EdgeColor','k');
set(h,'LineStyle','--','linewid',1, 'FaceColor',[0 .5 .5]);
h=rectangle('Position',[O4(2)-RO4,O4(1)-RO4,2*RO4,2*RO4],'Curvature',[1,1],'EdgeColor','k');
set(h,'LineStyle','--','linewid',1, 'FaceColor',[0 .5 .5]);


% USV1-4
for k=1:2000:Ns
    pos1 = [USV1.x(k,4) USV1.x(k,5)]'; 
    pos2 = [USV2.x(k,4) USV2.x(k,5)]'; 
    pos3 = [USV3.x(k,4) USV3.x(k,5)]'; 
    pos4 = [USV4.x(k,4) USV4.x(k,5)]'; 
    modelplot(pos1,USV1.x(k,6),'r-',linewid);
    modelplot(pos2,USV2.x(k,6),'g-',linewid);
    modelplot(pos3,USV3.x(k,6),'b-',linewid);
    modelplot(pos4,USV4.x(k,6),'c-',linewid);
    plot(pi1(k,2),pi1(k,1),'ro'); 
    plot(pi2(k,2),pi2(k,1),'go'); 
    plot(pi3(k,2),pi3(k,1),'bo'); 
    plot(pi4(k,2),pi4(k,1),'co'); 
    
    plot(pr1(k,2),pr1(k,1),'ro'); 
    plot(pr2(k,2),pr2(k,1),'go'); 
end

h1=legend([plt1(1),plt2(1),plt3(1),plt4(1)],{'USV1','USV2','USV3','USV4'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight); set(h1,'orientation','horizontal'); 

xlabel('y (m)','FontSize',fontsize,'FontName',fontname);
ylabel('x (m)','FontSize',fontsize,'FontName',fontname);
box on
axis([Xmin Xmax,Ymin Ymax]);
set(gca,'xTick',Xmin:Xinterval:Xmax);
set(gca,'yTick',Ymin:Yinterval:Ymax);

hold off


figure(2)
plt1 = plot(t,USV1.tau(:,1),'r-');hold on
plt2 = plot(t,USV2.tau(:,1),'g-');
plt3 = plot(t,USV3.tau(:,1),'b-');
plt4 = plot(t,USV4.tau(:,1),'c-'); hold off
h1=legend([plt1(1),plt2(1),plt3(1),plt4(1)],...
    {'$$\tau_{1,u}$$','$$\tau_{2,u}$$','$$\tau_{3,u}$$','$$\tau_{4,u}$$'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight); set(h1,'orientation','vertical'); 
ylim([-2 2]); 
xlabel('time (s)');
ylabel('thrust (N)');

figure(3)
plt1 = plot(t,USV1.tau(:,2),'r-');hold on
plt2 = plot(t,USV2.tau(:,2),'g-');
plt3 = plot(t,USV3.tau(:,2),'b-');
plt4 = plot(t,USV4.tau(:,2),'c-'); hold off
h1=legend([plt1(1),plt2(1),plt3(1),plt4(1)],...
    {'$$\tau_{1,u}$$','$$\tau_{2,u}$$','$$\tau_{3,u}$$','$$\tau_{4,u}$$'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight); set(h1,'orientation','vertical'); 
ylim([-1.5 1.5]); 
xlabel('time (s)');
ylabel('moment (Nm)');

figure(4)
subplot(311)
plt1 = plot(t,USV1.x(:,1),'r-','linewid',linewid); hold on
plt2 = plot(t,USV2.x(:,1),'g-','linewid',linewid); 
plt3 = plot(t,USV3.x(:,1),'b-','linewid',linewid); 
plt4 = plot(t,USV4.x(:,1),'c-','linewid',linewid); hold off
ylim([-0.2 0.8]);
h1=legend([plt1(1),plt2(1),plt3(1),plt4(1)],{'USV1','USV2','USV3','USV4'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight); set(h1,'orientation','horizontal');
xlabel('time (s)');
h1 = ylabel('$$u_i$$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

subplot(312);
plot(t,USV1.x(:,2),'r-','linewid',linewid); hold on
plot(t,USV2.x(:,2),'g-','linewid',linewid); 
plot(t,USV3.x(:,2),'b-','linewid',linewid); 
plot(t,USV4.x(:,2),'c-','linewid',linewid); hold off
xlabel('time (s)');
h1 = ylabel('$$v_i$$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
ylim([-0.3 0.3]);

subplot(313);
plot(t,USV1.x(:,3),'r-','linewid',linewid); hold on
plot(t,USV2.x(:,3),'g-','linewid',linewid); 
plot(t,USV3.x(:,3),'b-','linewid',linewid); 
plot(t,USV4.x(:,3),'c-','linewid',linewid);  hold off
xlabel('time (s)');
h1 = ylabel('$$r_i$$ (rad/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
ylim([-0.7 0.7]);

figure(5)
plot(t,z1_bar_norm,'r-',t,z2_bar_norm,'b-',t,z3_bar_norm,'g-',t,z4_bar_norm,'c-');
xlabel('time (s)');
ylabel('formation errors (m)');
h1 = legend('$$||\bar{z}_{1,e}||$$','$$||\bar{z}_{2,e}||$$','$$||\bar{z}_{3,e}||$$','$$||\bar{z}_{4,e}||$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight); set(h1,'orientation','vertical'); 

figure(6)
title('enviroment disturbances');
subplot(311)
plot(t,wk_u,'r-','linewid',linewid);
xlabel('time (s)');
h1 = ylabel('$$wk_u$$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

subplot(312);
plot(t,wk_u,'r-','linewid',linewid);
xlabel('time (s)');
h1 = ylabel('$$wk_u$$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

subplot(313);
plot(t,wk_u,'r-','linewid',linewid);
xlabel('time (s)');
h1 = ylabel('$$wk_u$$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

% triggering time
figure(7)
subplot(211), hold on
for i = 1:Ns
    if Tk1(i,1)==1
        plot(t(i),0.5*Tk1(i,1),'rx');
    end
    if Tk2(i,1)==1
        plot(t(i),1*Tk2(i,1),'gx');
    end
    if Tk3(i,1)==1
        plot(t(i),1.5*Tk3(i,1),'bx');
    end
    if Tk4(i,1)==1
        plot(t(i),2*Tk4(i,1),'cx');
    end
end
hold off
ylim([0 2.5]);
subplot(212), hold on
for i = 1:Ns
    if Tk1(i,2)==1
        plot(t(i),0.5*Tk1(i,2),'rx');
    end
    if Tk2(i,2)==1
        plot(t(i),1*Tk2(i,2),'gx');
    end
    if Tk3(i,2)==1
        plot(t(i),1.5*Tk3(i,2),'bx');
    end
    if Tk4(i,2)==1
        plot(t(i),2*Tk4(i,2),'cx');
    end
end
hold off
ylim([0 2.5]);

% CMG error
figure(8)
plot(t,eF1_norm,'r-',t,eF2_norm,'g-',t,eF3_norm,'b-',t,eF4_norm,'c-');
h1 = legend('$$e_{1}^{[F]}$$','$$e_{2}^{[F]}$$','$$e_{3}^{[F]}$$','$$e_{4}^{[F]}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = ylabel('$$e^{F}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = xlabel('time (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

figure(9)
plot(t,eL1,'r-',t,eL2,'b-');
h1 = legend('$$e_{1}^{[L]}$$','$$e_{2}^{[L]}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = ylabel('$$e^{L}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = xlabel('time (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

figure(10)
subplot(211),plot(t,USV1.x(:,1),'r-',t,USV1.uhat,'b--');
h1 = legend('$$u_1$$','$$\hat{u}_1$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = ylabel('$$\hat{u}_1$$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = xlabel('time (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

subplot(212),plot(t,USV1.v_bar,'r-',t,USV1.v_bar_hat,'b--');
h1 = legend('$$\bar{v}_1$$','$$\hat{v}_1$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = ylabel('$$\hat{v}_1$$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = xlabel('time (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

figure(11)
subplot(211),plot(t,USV1.zeta_bar(:,1),'r-',t,USV1.zetahat(:,1),'b--');
h1 = legend('$$\zeta_{1,u}$$','$$\hat{\zeta}_{1,u}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = ylabel('$$\hat{zeta}_{1,u}$$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = xlabel('time (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

subplot(212),plot(t,USV1.zeta_bar(:,2),'r-',t,USV1.zetahat(:,2),'b--');
h1 = legend('$$\zeta_{1,r}$$','$$\hat{\zeta}_{1,r}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = ylabel('$$\hat{\zeta}_{1,r}$$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = xlabel('time (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

figure(12)
subplot(211),plot(t,USV1.e_bar(:,1),'r-',t,USV2.e_bar(:,1),'b-',t,USV3.e_bar(:,1),'c-',t,USV4.e_bar(:,1),'m-');
h1 = legend('$$\bar(e)_{1,u}$$','$$\bar{e}_{2,u}$$','$$\bar(e)_{3,u}$$','$$\bar{e}_{4,u}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = ylabel('$$\bar{e}_u$$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = xlabel('time (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

subplot(212),plot(t,USV1.e_bar(:,2),'r-',t,USV2.e_bar(:,2),'b--',t,USV3.e_bar(:,2),'c-',t,USV4.e_bar(:,2),'m--');
h1 = legend('$$\bar(e)_{1,u}$$','$$\bar{e}_{2,u}$$','$$\bar(e)_{3,u}$$','$$\bar{e}_{4,u}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = ylabel('$$\bar{e}_r$$ rad/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = xlabel('time (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

% animation
% area = [Xmin Xmax Ymin Ymax];
% X = [USV1.x(:,4),USV2.x(:,4),USV3.x(:,4),USV4.x(:,4)];
% Y = [USV1.x(:,5),USV2.x(:,5),USV3.x(:,5),USV4.x(:,5)];
% Psi = [USV1.x(:,6),USV2.x(:,6),USV3.x(:,6),USV4.x(:,6)];
% CMGX = [pi1(:,1),pi2(:,1),pi3(:,1),pi4(:,1)];
% CMGY = [pi1(:,2),pi2(:,2),pi3(:,2),pi4(:,2)];
% VLX = [pr1(:,1),pr2(:,1)];
% VLY = [pr1(:,2),pr2(:,2)];
% dt = 0.1;
% myAnimation( X, Y, Psi, CMGX, CMGY, VLX, VLY, area, dt, 13 );









