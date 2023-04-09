function myAnimation( X, Y, Psi, pathX, pathY, area, figNum )
%MYANIMATION �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��


figure(figNum)

shipNum = length(X(1,:)); % �жϴ��ĸ���
N = length(X(:,1)); % �ж����ݵ����

dt = 0.005;
color = ['r','g','b','c','m'];
linewid = 1;

X1 = [];
Y1 = [];
pathX1 = [];
pathY1 = [];

k1 = 1;
for k=1:100:N
    hold off
    for i=1:shipNum
        pos = [X(k,i),Y(k,i)]'; psi = Psi(k,i);
        modelplot( pos, psi, color(i),linewid ); hold on
        X1(k1,i) = X(k,i);
        Y1(k1,i) = Y(k,i);
        plot(Y1(:,i),X1(:,i),color(i),'linewid',linewid); hold on 
    end
    pathX1(k1,1) = pathX(k,1);
    pathY1(k1,1) = pathY(k,1);
    plot(pathY1,pathX1,'k-','linewid',linewid); hold on 
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



















