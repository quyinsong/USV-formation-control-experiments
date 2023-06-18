function [ x_bar,y_bar] = CoordinateTrans( x,y,psi,emi )
%COORDINATETRANS 此处显示有关此函数的摘要
%   此处显示详细说明

x_bar = x+emi*cos(psi);
y_bar = y+emi*sin(psi);

end

