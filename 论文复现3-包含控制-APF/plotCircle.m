function  plotCircle( R, pos, color, linewid, isFill )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

alpha=0:pi/20:2*pi;%角度[0,2*pi]
xo = pos(1); yo = pos(2);
x=R*cos(alpha)+xo;
y=R*sin(alpha)+yo;
plot(y,x,color,'linewidth',linewid)
if isFill
    fill(x,y,color);%填充
end
 
end

