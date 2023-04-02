function myAnimation( X, Y, Psi, pathX, pathY, area, figNum )
%MYANIMATION 此处显示有关此函数的摘要
%   此处显示详细说明


figure(figNum)

shipNum = length(X(1,:)); % 判断船的个数
N = length(X(:,1)); % 判断数据点个数

dt = 0.005;
color = ['r','g','b','c','m'];
linewid = 1;

X1 = [];
Y1 = [];

k1 = 1;
for k=1:20:N
    hold off
    for i=1:shipNum
        pos = [X(k,i),Y(k,i)]'; psi = Psi(k,i);
        modelplot( pos, psi, color(i),linewid ); hold on
        X1(k1,i) = X(k,i);
        Y1(k1,i) = Y(k,i);
        plot(Y1(:,i),X1(:,i),color(i),'linewid',linewid); hold on 
    end
    plot(pathY,pathX,'k-','linewid',linewid); hold on 
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



















