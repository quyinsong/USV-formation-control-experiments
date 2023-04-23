function myPlotLine( pos_start,pos_end )
%MYPLOTLINE 此处显示有关此函数的摘要
%   此处显示详细说明

x1 = pos_start(1);
y1 = pos_start(2);
x2 = pos_end(1);
y2 = pos_end(2);

if (x2-x1) >0
    x = x1:0.5:x2;
    k = (y2-y1)/(x2-x1);
    y = k*(x-x1)+y1;
elseif (x2-x1)<0
    x = x2:0.5:x1;
    k = (y2-y1)/(x2-x1);
    y = k*(x-x1)+y1;
else
    y = y1:0.5:y2;
    x = x1*ones(1,length(y));
end
plot(y,x,'k.','linewid',0.3,'Markersize',0.3);

end

