clc
clear all
close all

x=[140, 130, 100, 100, 100, 120, 150, 160, 160];
y=[10, 40, 70, 90, 100, 130,150,160, 180];
p = polyfit(y, x, 4) % 三次多项式拟合
yy = 10: 1 : 180;
xx = polyval(p, yy) ; % 根据系数向量p计算在xx点处的函数值

x1=[20, 50, 40, 40, 40, 70, 90, 110, 110];
y1=[10, 40, 70, 100, 120 ,140, 160, 170, 190];
p1 = polyfit(y1, x1, 4) % 三次多项式拟合
yy1 = 10: 1 : 190;
xx1 = polyval(p1, yy1) ; % 根据系数向量p计算在xx点处的函数值

yrange=[0 30 200]; xrange = [0 30 200];

Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);

figure(1)
plot(yy, xx, 'ro', y, x, 'b-');hold on
plot(yy1, xx1, 'ro', y1, x1, 'b-');hold on

axis([Xmin Xmax,Ymin Ymax]);
set(gca,'xTick',Xmin:Xinterval:Xmax);
set(gca,'yTick',Ymin:Yinterval:Ymax);
% xxx = 20:0.1:120;
% yyy = 207.3326*xxx.^3-268.2537*xxx.^2+92.9394.*xxx-3.0873;
% plot(xxx', yyy', 'ro'); hold off