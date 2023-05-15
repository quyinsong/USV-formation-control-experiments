function myAnimation( X, Y, Psi, CMGX, CMGY, VLX, VLY, area, dt, figNum )
%MYANIMATION 此处显示有关此函数的摘要
%   此处显示详细说明


figure(figNum)

shipNum = length(X(1,:)); % 判断船的个数
N = length(X(:,1)); % 判断数据点个数

color = ['r','g','b','c','m'];
linewid = 1;

X1 = [];
Y1 = [];
pathX1 = [];
pathY1 = [];
CMGX1 = [];
CMGY2 = [];
VLX1 = [];
VLY2 = [];

k1 = 1;
for k=1:100:N
    hold off
    for i=1:shipNum
        pos = [X(k,i),Y(k,i)]'; psi = Psi(k,i);
        % 绘制无人船形状
        modelplot( pos, psi, color(i),linewid ); hold on
        % 绘制CMG形状和虚拟领航员
        h=rectangle('Position',[CMGY(k,i)-1,CMGX(k,i)-1,2*1,2*1],'Curvature',[1,1],'EdgeColor','k');
        set(h,'LineStyle','-','linewid',1);
        h=rectangle('Position',[VLY(k,i)-1,VLX(k,i)-1,2*1,2*1],'Curvature',[1,1],'EdgeColor','k');
        set(h,'LineStyle','-','linewid',1);
        % 绘制无人船和CMG以及虚拟领航员运动轨迹
        X1(k1,i) = X(k,i);
        Y1(k1,i) = Y(k,i);
        CMGX1(k1,i) = CMGX(k,i);
        CMGY1(k1,i) = CMGY(k,i);
        VLX1(k1,i) = VLX(k,i);
        VLY1(k1,i) = VLY(k,i);
        plot(Y1(:,i),X1(:,i),color(i),'linewid',linewid); hold on 
        plot(CMGY1(:,i),CMGX1(:,i),'k-','linewid',linewid); hold on   
        plot(VLY1(:,i),VLX1(:,i),'k-','linewid',linewid); hold on 
    end
    k1=k1+1;
    grid on;  
    xlabel('y / m');
    ylabel('x / m');
    title('CPF');
    axis(area);
    drawnow;
    pause(dt);
end

end



















