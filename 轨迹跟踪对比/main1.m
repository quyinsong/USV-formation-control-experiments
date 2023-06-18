%  ref: Gu, N., Wang, D., Peng, Z., Li, T., & Tong, S. (2021). 
%       Model-Free Containment Control of Underactuated Surface Vessels Under Switching Topologies
%       Based on Guiding Vector Fields and Data-Driven Neural Predictors. IEEE Transactions on Cybernetics, 1–12.
% by: Yinsong Qu
% date: 4, 22, 2023 

% 两种轨迹跟踪方法对比

clc
clear all
close all

ts = 0.05;
tfinal = 200;
Ns = tfinal/ts;

%% parameters initializition
% USV states
USV1.x0 = [0 0 0 12 8 10*pi/180]';
USV2.x0 = [0 0 0 12 8 10*pi/180]';
% USV inputs
USV1.tauc=[0 0]';
USV2.tauc=[0 0]';

for k=1:Ns
   % time series
   t = (k-1)*ts;
   
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
   tau_w = [0.1*sin(0.02*t)+0.2*sin(0.05*t) 0.05*sin(0.05*t) ...
            0.05*sin(0.02*t)+0.05*sin(0.05*t)]';

   % USV states update
   
   [USV1.x,USV1.tau,USV1.f, USV1.zeta, USV1.emi, USV1.Gamma] = CS1tf( USV1.x0, USV1.tauc, tau_w, ts );
   [USV2.x,USV2.tau,USV2.f, USV2.zeta, USV2.emi, USV2.Gamma] = CS2tf( USV2.x0, USV2.tauc, tau_w, ts );
   
   USV1.zeta_bar = [USV1.zeta(1),USV1.zeta(3)]';
   % USV states and Coordinate transformation
   USV1.eta = [USV1.x(4) USV1.x(5) USV1.x(6)]';
   USV2.eta = [USV2.x(4) USV2.x(5) USV2.x(6)]';
   
   USV1.p = [USV1.eta(1);USV1.eta(2)];
   USV2.p = [USV2.eta(1);USV2.eta(2)];
    
   USV1.nu = [USV1.x(1) USV1.x(2) USV1.x(3)]';USV1.nu_bar = [USV1.x(1) USV1.x(3)]';
   USV2.nu = [USV2.x(1) USV2.x(2) USV2.x(3)]';USV2.nu_bar = [USV2.x(1) USV2.x(3)]';
   
   USV1.v_bar = USV1.nu(2)+USV1.emi*USV1.nu(3);
   USV1.nu_tf = [USV1.x(1),USV1.v_bar,USV1.x(3)]';
   USV2.v_bar = USV2.nu(2)+USV2.emi*USV2.nu(3);
   USV2.nu_tf = [USV2.x(1),USV2.v_bar,USV2.x(3)]';
   
   [ USV1.x_bar, USV1.y_bar] = CoordinateTrans( USV1.eta(1),USV1.eta(2),USV1.eta(3),USV1.emi );
   [ USV2.x_bar, USV2.y_bar] = CoordinateTrans( USV2.eta(1),USV2.eta(2),USV2.eta(3),USV1.emi );
   
   USV1.p_tf = [ USV1.x_bar, USV1.y_bar]';
   USV2.p_tf = [ USV2.x_bar, USV2.y_bar]';
   
   USV1.eta_tf = [ USV1.x_bar, USV1.y_bar,USV1.x(6)]';
   USV2.eta_tf = [ USV2.x_bar, USV2.y_bar,USV2.x(6)]';
   
   R_psi1 = [cos(USV1.eta(3)) -sin(USV1.eta(3));
             sin(USV1.eta(3)) cos(USV1.eta(3))];
   R_psi2 = [cos(USV2.eta(3)) -sin(USV2.eta(3));
             sin(USV2.eta(3)) cos(USV2.eta(3))];
         
   % desired trajectory
   if t == 0
       xd = 8;
       yd = 8;
   end
   xd_dot = 0.2;
   yd_dot = 0.2;
   xd = ts*xd_dot+xd;
   yd = ts*yd_dot+yd;
   pd = [xd,yd]';
   pd_dot = [xd_dot,yd_dot]';
   % tracking errors
   if t == 0
       delta0 = 0.1;
   end
   d = [delta0,0]';
   z1 = R_psi1'*(USV1.p_tf-pd);
   z2 = z1+d;
   hi = diag([1,delta0]);

   % ESO
   USV1.nu_hat = ESO1(USV1.eta_tf,ts);
   USV2.nu_hat = ESO2(USV2.eta_tf,ts);
   USV1.u_hat = USV1.nu_hat(1); USV2.u_hat = USV2.nu_hat(1);
   USV1.v_bar_hat = USV1.nu_hat(2); USV2.v_bar_hat = USV2.nu_hat(2);
   USV1.r_hat = USV1.nu_hat(3); USV2.r_hat = USV2.nu_hat(3);
   
   % Guidance
   [USV1.nu_c,boundary1,boundary2,q1,vartheta1,Thetau1,Thetar1] = Guidance1( hi, R_psi1, z2, pd_dot, USV1.nu_hat, t, ts);
   alpha = [1,1]';
   [USV2.nu_c,pe] = Guidance2( USV2.p_tf,USV2.eta(3),pd,pd_dot,alpha,ts);
   z3 = R_psi2'*pe;
   % Control   
   [USV1.tauc,USV1.e_bar,USV1.zetahat] = nnctr1( USV1.nu_hat, USV1.nu_c, USV1.tau, USV1.tauc, USV1.Gamma,q1,vartheta1,USV1.zeta,ts );
   [USV2.tauc,USV2.e_bar,USV2.zetahat] = nnctr2( USV2.nu_hat, USV2.nu_c, USV2.tau, USV2.tauc, USV2.Gamma,USV2.zeta,ts );

   
   xout(k,:) = [t,USV1.x', USV2.x',pd', z2', z3', boundary1, boundary2 ];
    
    
end
%% simulation data
t = xout(:,1);
USV1.x = xout(:,2:7);
USV2.x = xout(:,8:13);
pd = xout(:,14:15);
z2 = xout(:,16:17);
z3 = xout(:,18:19);
boundary1 = xout(:,20);
boundary2 = xout(:,21);

%% PLOTS
close all
% formation
figure(1); hold on

linewid = 1;
fontsize = 10.5;
fontname = 'Time news roman';
fontweight = 'bold';

yrange=[0 10 50]; xrange = [0 10 50];

Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);
plt1 = plot(USV1.x(:,5),USV1.x(:,4),'r-'); 
plt2 = plot(USV2.x(:,5),USV2.x(:,4),'g-'); 

plot(pd(:,2),pd(:,1),'m--'); 

% USV1-4
for k=1:1000:Ns
    pos1 = [USV1.x(k,4) USV1.x(k,5)]'; 
    pos2 = [USV2.x(k,4) USV2.x(k,5)]';  
    modelplot(pos1,USV1.x(k,6),'r-',linewid);
    modelplot(pos2,USV2.x(k,6),'g-',linewid);
    plot(pd(k,2),pd(k,1),'ms','Markersize',5); 
end

set(gca,'linewidth',1)

h1=legend([plt1(1),plt2(1)],{'TTG method proposed in this paper','TTG method in [11]'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight); set(h1,'orientation','vertical'); 

h1 = xlabel('y (m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

h1 = ylabel('x (m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

box on
axis([Xmin Xmax,Ymin Ymax]);
set(gca,'xTick',Xmin:Xinterval:Xmax);
set(gca,'yTick',Ymin:Yinterval:Ymax);

hold off

figure(2)
plot(t,z2(:,1),'r-',t,z3(:,1),'b-'); hold on
plot(t,5*boundary1,'k--',t,-5*boundary2,'k--'); hold off
h1 = xlabel('time (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

h1 = ylabel('$$z(1)$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight); 

h1=legend('$$z_{1}(1)$$','$$z_{2}(1)$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight); set(h1,'orientation','horizontal'); 

figure(3)
plot(t,z2(:,2),'r-',t,z3(:,2),'b-'); hold on
plot(t,5*boundary1,'k--',t,-5*boundary2,'k--'); hold off
h1 = xlabel('time (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

h1 = ylabel('$$z(2)$$ (m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

h1=legend('$$z_{1}(2)$$','$$z_{2}(2)$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight); set(h1,'orientation','horizontal'); 















