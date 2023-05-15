function [ pio,pi_dot] = CMG3( eF,pj_dot,vs,di,p_i0,pi_o, pi_usv,ts )
%CMG1 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

persistent p_i Kg
if isempty(p_i)
    p_i = p_i0;
    Kg = diag([0.1 0.1]);
end
pi_dot = (-Kg*eF+pj_dot+vs)/di-pi_o-pi_usv;
p_i = ts*pi_dot+p_i;
pio = p_i;

end

