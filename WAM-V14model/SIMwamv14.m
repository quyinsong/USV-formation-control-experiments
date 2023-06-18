 % SIMotter    User editable script for simulation of the Otter USV under
%             feedback control
%
% Calls:      otter.m
%
% Author:    Thor I. Fossen
% Date:      2021-04-25
% Revisions: 
% new version: update by quyinsong  date: 2022 5 27
clc
clear all
close all
%% USER INPUTS
ts  = 0.1;        % sampling time [s]
tfinal = 80;
Ns  = tfinal/ts;		  % number of samples

% initial values for x = [ u v r x y psi ]'
x0 = zeros(6,1);	   

% propeller revolutions (rps)
Thrustc = [150;150];

%% MAIN LOOP
for k=1:Ns+1
   t = (k-1) * ts;                          % time (s)             
   % get otter state output
   tau_w = [0 0 0]';
   [x, Thrust, f] = WAMV14( x0, Thrustc, tau_w, ts );
   if t>=20
       Thrustc = [150;120];
   end
   % store simulation data in a table 
   if t== 0
       Thrust = [0 0]';
   end
   simdata(k,:) = [t x' Thrustc' Thrust']; 
end

%% PLOTS
t    = simdata(:,1); 
nu   = simdata(:,2:4); 
eta  = simdata(:,5:7); 
Thrustc = simdata(:,8:9);
Thrust = simdata(:,10:11);

disp('plot ...');

figure(1)
plot(t,nu(:,1),'linewidt',2)
xlabel('time (s)'),title('Surge velocity (m/s)'),grid
figure(2)
plot(t,nu(:,2),'linewidt',2)
xlabel('time (s)'),title('Sway velocity (m/s)'),grid
figure(3)
plot(t,(180/pi)*nu(:,3),'linewidt',2)
xlabel('time (s)'),title('Yaw rate (deg/s)'),grid

figure(4);hold on
xpos = eta(:,1); ypos = eta(:,2); psi = eta(:,3);
xrange=[-50 50 100]; yrange = [-50 50 150];
for k=1:100:Ns
    pos1 = [eta(k,1) eta(k,2)]'; 
    modelplot(pos1,eta(k,3),'r-',1);
end
Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);
plot(ypos,xpos,'r','linewidth',2)
box on
axis([Xmin Xmax,Ymin Ymax]);
set(gca,'xTick',Xmin:Xinterval:Xmax);
set(gca,'yTick',Ymin:Yinterval:Ymax);
grid on
hold off;

figure(5)
plot(t,Thrustc(:,1),'r--',t,Thrustc(:,2),'b--',...
     t,Thrust(:,1),'r-',t,Thrust(:,2),'b-','linewidt',2)
legend('Tc1','Tc2','T1','T2');
xlabel('time (s)'),title('Force (N)'),grid









