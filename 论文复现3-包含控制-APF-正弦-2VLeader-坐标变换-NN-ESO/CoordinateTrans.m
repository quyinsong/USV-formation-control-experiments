function [ x_bar,y_bar] = CoordinateTrans( x,y,psi,emi )
%COORDINATETRANS �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

x_bar = x+emi*cos(psi);
y_bar = y+emi*sin(psi);

end

