function  plotCircle( R, pos, color, linewid, isFill )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

alpha=0:pi/20:2*pi;%�Ƕ�[0,2*pi]
xo = pos(1); yo = pos(2);
x=R*cos(alpha)+xo;
y=R*sin(alpha)+yo;
plot(y,x,color,'linewidth',linewid)
if isFill
    fill(x,y,color);%���
end
 
end

