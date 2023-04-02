% Formation control 
% Date: 7 13 2022
% Author : Yinsong Qu
% Reference : Adaptive Dynamic Surface Control for Formations of Autonomous Surface Vehicles With Uncertain Dynamics
% 编队方法： leader-follower
% 制导方法： 视线制导 LOS
% 优点： 只需要测得与领航船的相对LOS距离和角度，不需要领航船的速度和航向角信息，采用神经网络和自适应方法处理外部干扰
% 缺点： 对于领航船不同的运动路径，不能保持好的队形
% 可改进的地方： 可将与领航船之间的跟踪误差转换到领航船坐标系中，获得新的制导率
clc
clear all
close all

ts = 0.02;
tfinal = 300;
Ns = tfinal/ts;

% USV0 : leader 
% USV1~4: follower

%% parameters initializition
% USV states
USV0.x0 = [0 0 0 0 4 90*pi/180]';
USV1.x0 = [0 0 0 15 0 90*pi/180]';
USV2.x0 = [0 0 0 14 0 90*pi/180]';
USV3.x0 = [0 0 0 -10 0 90*pi/180]';
USV4.x0 = [0 0 0 -4 0 90*pi/180]';
% USV inputs
USV0.tauc=[0 0];
USV1.tauc=[0 0];
USV2.tauc=[0 0];
USV3.tauc=[0 0];
USV4.tauc= [0 0];
% desired formation
%--------------------------
%   USV1            USV3
%           USV0
%   USV4            USV2
%------------------------
s1c = [8 45*pi/180]';
s2c = [8 135*pi/180]';
s3c = [8 -135*pi/180]';
s4c = [8 -45*pi/180]';

% desired speed
ud = 0.6;

%% simulation
% preallocated memory
rho =zeros(4,1);
lambda = zeros(4,1);
% simulation start
for k = 1:Ns
    
   % time series
   t = (k-1)*ts;
   % USV states update

   [USV0.x,USV0.tau] = ASV0( USV0.x0, USV0.tauc, ts, t );
   [USV1.x,USV1.tau] = ASV1( USV1.x0, USV1.tauc, ts, t );
   [USV2.x,USV2.tau] = ASV2( USV2.x0, USV2.tauc, ts, t );
   [USV3.x,USV3.tau] = ASV3( USV3.x0, USV3.tauc, ts, t );
   [USV4.x,USV4.tau] = ASV4( USV4.x0, USV4.tauc, ts, t );
   % USV states
   USV0.eta = [USV0.x(4) USV0.x(5) USV0.x(6)]';
   USV1.eta = [USV1.x(4) USV1.x(5) USV1.x(6)]';
   USV2.eta = [USV2.x(4) USV2.x(5) USV2.x(6)]';
   USV3.eta = [USV3.x(4) USV3.x(5) USV3.x(6)]';
   USV4.eta = [USV4.x(4) USV4.x(5) USV4.x(6)]';
   USV0.nu = [USV0.x(1) USV0.x(2) USV0.x(3)]'; 
   USV1.nu = [USV1.x(1) USV1.x(2) USV1.x(3)]';USV1.nu_bar = [USV1.x(1) USV1.x(3)]';
   USV2.nu = [USV2.x(1) USV2.x(2) USV2.x(3)]';USV2.nu_bar = [USV2.x(1) USV2.x(3)]';
   USV3.nu = [USV3.x(1) USV3.x(2) USV3.x(3)]';USV3.nu_bar = [USV3.x(1) USV3.x(3)]';
   USV4.nu = [USV4.x(1) USV4.x(2) USV4.x(3)]';USV4.nu_bar = [USV4.x(1) USV4.x(3)]';

   % path following control for leader USV0 
   pc = 4;
   [u0d, r0d, pos0d] = LOS( pc, USV0.eta, USV0.nu, ud, ts );
   USV0.tauc = pathfollowingcontrollaw( u0d,r0d,USV0.nu,USV0.tau,ts );
   
   % calculate LOS range and angle
   rx = [USV1.eta(1);USV2.eta(1);USV3.eta(1);USV4.eta(1)]-repmat(USV0.eta(1),4,1);
   ry = [USV1.eta(2);USV2.eta(2);USV3.eta(2);USV4.eta(2)]-repmat(USV0.eta(2),4,1);
   
   for i = 1:4
       [rho(i),lambda(i)] = Lidar(rx(i),ry(i));
   end
   
   % formation change
   if t>80
        s1c = [4 45*pi/180]';
        s2c = [4 135*pi/180]';
        s3c = [4 -135*pi/180]';
        s4c = [4 -45*pi/180]';
   end
   
   if t>200
        s1c = [8 45*pi/180]';
        s2c = [8 135*pi/180]';
        s3c = [8 -135*pi/180]';
        s4c = [8 -45*pi/180]';
   end
   
   psij_f = TD1( USV0.eta(3),ts );
%    s1c = [4 90*pi/180-(psij_f-(-45*pi/180))]';
%    s2c = [2 90*pi/180-(psij_f-(-45*pi/180))]';
%    s3c = [2 90*pi/180-(psij_f-45*pi/180)]';
%    s4c = [4 90*pi/180-(psij_f-45*pi/180)]';

%    Rj = [cos(USV0.eta(3)) -sin(USV0.eta(3)); sin(USV0.eta(3)) cos(USV0.eta(3))];
%    s1d = Rj'*[-4 -4]';s2d = Rj'*[-8 -4]';s3d = Rj'*[-8 4]';s4d = Rj'*[-4 4]';
%    s1c = [norm(s1d) atan2(USV0.eta(2)-s1d(2),USV0.eta(1)-s1d(1))]';
%    s2c = [norm(s2d) atan2(USV0.eta(2)-s2d(2),USV0.eta(1)-s2d(1))]';
%    s3c = [norm(s3d) atan2(USV0.eta(2)-s3d(2),USV0.eta(1)-s3d(1))]';
%    s4c = [norm(s4d) atan2(USV0.eta(2)-s4d(2),USV0.eta(1)-s4d(1))]';

   % kinematic congtrol law
   s1 = [rho(1) lambda(1)]';
   [ u1d, r1d, u1d_d, r1d_d, se1, f1, fhat1] = kinematiclaw1( s1, s1c, USV1.eta(3), USV0.nu,  USV0.eta(3), USV1.nu, ts);
   s2 = [rho(2) lambda(2)]';
   [ u2d, r2d, u2d_d, r2d_d, se2, f2, fhat2 ] = kinematiclaw2( s2, s2c, USV2.eta(3), USV0.nu, USV0.eta(3), USV2.nu, ts);
   s3 = [rho(3) lambda(3)]';
   [ u3d, r3d, u3d_d, r3d_d, se3, f3, fhat3 ] = kinematiclaw3( s3, s3c, USV3.eta(3), USV0.nu, USV0.eta(3), USV3.nu, ts);
   s4 = [rho(4) lambda(4)]';
   [ u4d, r4d, u4d_d, r4d_d, se4, f4, fhat4 ] = kinematiclaw4( s4, s4c, USV4.eta(3), USV0.nu, USV0.eta(3), USV4.nu, ts);
   % final control law
   USV1.nu_bard = [u1d r1d]';USV1.nu_bard_d = [u1d_d r1d_d]';
   USV2.nu_bard = [u2d r2d]';USV2.nu_bard_d = [u2d_d r2d_d]';
   USV3.nu_bard = [u3d r3d]';USV3.nu_bard_d = [u3d_d r3d_d]';
   USV4.nu_bard = [u4d r4d]';USV4.nu_bard_d = [u4d_d r4d_d]';
   
   USV1.tauc = finalcontrollaw1( USV1.nu, USV1.nu_bar, USV1.nu_bard,USV1.nu_bard_d, ts );
   USV2.tauc = finalcontrollaw2( USV2.nu, USV2.nu_bar, USV2.nu_bard,USV2.nu_bard_d, ts );
   USV3.tauc = finalcontrollaw3( USV3.nu, USV3.nu_bar, USV3.nu_bard,USV3.nu_bard_d, ts );
   USV4.tauc = finalcontrollaw4( USV4.nu, USV4.nu_bar, USV4.nu_bard,USV4.nu_bard_d, ts );

%    USV1.tau = [5 0]';
%    USV2.tau = [5 0]';
%    USV3.tau = [5 0]';
%    USV4.tau = [5 0]';
   % store simulation data 
   xout(k,:) = [t, USV0.x', USV1.x', USV2.x', USV3.x', USV4.x' USV1.tau' USV2.tau' USV3.tau' USV4.tau' se1' se2' se3' se4' pos0d'...
                s1c' s2c' s3c' s4c' psij_f f1' fhat1'];
    
end
%% simulation data
t = xout(:,1);
USV0.x = xout(:,2:7);
USV1.x = xout(:,8:13);
USV2.x = xout(:,14:19);
USV3.x = xout(:,20:25);
USV4.x = xout(:,26:31);
USV1.tau = xout(:,32:33);
USV2.tau = xout(:,34:35);
USV3.tau = xout(:,36:37);
USV4.tau = xout(:,38:39);
se1 = xout(:,40:41); 
se2 = xout(:,42:43); 
se3 = xout(:,44:45); 
se4 = xout(:,46:47); 
pos0d = xout(:,48:49);
s1c = xout(:,50:51); 
s2c = xout(:,52:53); 
s3c = xout(:,54:55); 
s4c = xout(:,56:57); 
psij_f = xout(:,58); 
f1 = xout(:,59:60);
fhat1 = xout(:,61:62);

%% PLOTS
% formation
figure(1); hold on

fontsize = 10; fontname = 'Times New Roman';
xrange=[-30 50 100]; yrange = [-30 10 30];
linewid = 1;

Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);
plot(USV0.x(:,5),USV0.x(:,4),'r--'); 
plot(USV1.x(:,5),USV1.x(:,4),'g--'); 
plot(USV2.x(:,5),USV2.x(:,4),'b--'); 
plot(USV3.x(:,5),USV3.x(:,4),'c--'); 
plot(USV4.x(:,5),USV4.x(:,4),'m--'); 
plot(pos0d(:,2),pos0d(:,1),'k-');
% USV0
rhoc = [s1c(:,1) s2c(:,1) s3c(:,1) s4c(:,1)];
lambdac = [s1c(:,2) s2c(:,2) s3c(:,2) s4c(:,2)];
NU = 2000;
for k=1:1:Ns
    pos =[USV0.x(k,4) USV0.x(k,5)]';
    if k==1
        modelplot(pos,USV0.x(k,6),'r',linewid);
        % calculate desired formation point and plot 
        posc = repmat(pos,1,4)-[rhoc(k,:).*sin(lambdac(k,:));rhoc(k,:).*cos(lambdac(k,:))]; % 2x4 the four follower desired formation point relative to leader
        xc = posc(1,:); yc = posc(2,:);
        for i =1:4
            x = [xc(i), pos(1)];
            y = [yc(i), pos(2)];
            plot(y,x,'r-','linewid',1);
        end
    end
    if rem(k,NU)==0
        modelplot( pos, USV0.x(k,6), 'r',linewid );
        % calculate desired formation point and plot     
        posc = repmat(pos,1,4)-[rhoc(k,:).*sin(lambdac(k,:));rhoc(k,:).*cos(lambdac(k,:))]; % 2x4 the four follower desired formation point relative to leader
        xc = posc(1,:); yc = posc(2,:);
        for i =1:4
            x = [xc(i), pos(1)];
            y = [yc(i), pos(2)];
            plot(y,x,'r-','linewid',1);
        end
    end   
end
% USV1-4
for k=1:2000:Ns
    pos1 = [USV1.x(k,4) USV1.x(k,5)]'; 
    pos2 = [USV2.x(k,4) USV2.x(k,5)]'; 
    pos3 = [USV3.x(k,4) USV3.x(k,5)]'; 
    pos4 = [USV4.x(k,4) USV4.x(k,5)]'; 
    modelplot(pos1,USV1.x(k,6),'g-',linewid);
    modelplot(pos2,USV2.x(k,6),'b-',linewid);
    modelplot(pos3,USV3.x(k,6),'c-',linewid);
    modelplot(pos4,USV4.x(k,6),'m-',linewid);
end


set(gca,'xTick',Xmin:Xinterval:Xmax);
set(gca,'yTick',Ymin:Yinterval:Ymax);
axis([Xmin Xmax,Ymin Ymax]);

xlabel('y (m)','FontSize',fontsize,'FontName',fontname);
ylabel('x (m)','FontSize',fontsize,'FontName',fontname);

hold off

% force and moment
figure(2)
subplot(4,1,1),plot(t,USV1.tau(:,1),'r-',t,USV1.tau(:,2),'b-');
subplot(4,1,2),plot(t,USV2.tau(:,1),'r-',t,USV2.tau(:,2),'b-');
subplot(4,1,3),plot(t,USV3.tau(:,1),'r-',t,USV3.tau(:,2),'b-');
subplot(4,1,4),plot(t,USV4.tau(:,1),'r-',t,USV4.tau(:,2),'b-');

% formation errors
figure(3)
subplot(2,1,1), plot(t,se1(:,1),'r-',t,se2(:,1),'g-',t,se3(:,1),'b-',t,se4(:,1),'c-');
h1=legend('$$\rho_{e1}$$','$$\rho_{e2}$$','$$\rho_{e3}$$','$$\rho_{e4}$$');
set(h1,'Interpreter','latex')
subplot(2,1,2), plot(t,se1(:,2),'r-',t,se2(:,2),'g-',t,se3(:,2),'b-',t,se4(:,2),'c-');
h1=legend('$$\lambda_{e1}$$','$$\lambda_{e2}$$','$$\lambda_{e3}$$','$$\lambda_{e4}$$');
set(h1,'Interpreter','latex')

figure(4)
plot(t,USV0.x(:,1),'r-',t,USV1.x(:,1),'g-',t,USV2.x(:,1),'b-',t,USV3.x(:,1),'c-',t,USV4.x(:,1),'m-');
legend('u0','u1','u2','u3','u4');

figure(5)
plot(t,psij_f,'r-','linewid',2)

% animation
pathX = pos0d(:,1); pathY = pos0d(:,2);
area = [Xmin Xmax Ymin Ymax];
X = [USV0.x(:,4),USV1.x(:,4),USV2.x(:,4),USV3.x(:,4),USV4.x(:,4)];
Y = [USV0.x(:,5),USV1.x(:,5),USV2.x(:,5),USV3.x(:,5),USV4.x(:,5)];
Psi = [USV0.x(:,6),USV1.x(:,6),USV2.x(:,6),USV3.x(:,6),USV4.x(:,6)];
myAnimation( X, Y, Psi, pathX, pathY, area, 6 );

% 制导律干扰估计效果
figure(7)
subplot(211)
plot(t,f1(:,1),'r-',t,fhat1(:,1),'b--','linewid',1);
h1=legend('$$f_1(1)$$','$$\hat{f}_1(1)$$');
set(h1,'Interpreter','latex')
set(h1,'FontSize',15)
subplot(212)
plot(t,f1(:,2),'r-',t,fhat1(:,2),'b--','linewid',1);
h1=legend('$$f_1(2)$$','$$\hat{f}_1(2)$$');
set(h1,'Interpreter','latex')
set(h1,'FontSize',15)



















