% Author: Yinsong Qu
% Date: 7 13 2022
% test USV model
clc
clear all
close all

ts = 0.02;
t_final = 20;
Ns = t_final/ts;
tau = [20;0];
x0 = [0;0;0;0;0;0];
for k = 1:Ns
    t = (k-1)*ts;
    x = USV( x0, tau, ts, t );
    xout(k,:)=[t x'];
end
t = xout(:,1);
u = xout(:,2);
v = xout(:,3);
r = xout(:,4);
x = xout(:,5);
y = xout(:,6);
psi = xout(:,7);


figure(1); hold on

fontsize = 10; fontname = 'Times New Roman';
xrange=[-2 5 50]; yrange = [-2 5 50];

Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);
plot(y,x); 
for k=1:1:Ns
    pos =[x(k) y(k)]';
    if k==1
        modelplot(pos,psi(k),xrange,yrange);
    end
    if rem(k,500)==0
        modelplot(pos,psi(k),xrange,yrange);
    end   
end
set(gca,'xTick',Xmin:Xinterval:Xmax);
set(gca,'yTick',Ymin:Yinterval:Ymax);
axis([Xmin Xmax,Ymin Ymax]);

xlabel('y (m)','FontSize',fontsize,'FontName',fontname);
ylabel('x (m)','FontSize',fontsize,'FontName',fontname);

hold off

figure(2)
plot(t,u);








    