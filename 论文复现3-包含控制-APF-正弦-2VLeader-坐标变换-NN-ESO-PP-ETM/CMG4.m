function [ pio,pi_dot] = CMG4( eF,pj_dot,vs,di,p_i0,pi_o, pi_usv,ts )
%CMG1 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

persistent p_i Kg
if isempty(p_i)
    p_i = p_i0;
    Kg = diag([0.1 0.1]);
end
pi_dot = (-Kg*(eF+pi_o+pi_usv)+pj_dot+vs)/di;

% if norm(pi_dot) >= 0.5
%     pi_dot = 0.5*pi_dot/norm(pi_dot);
% end

p_i = ts*pi_dot+p_i;
pio = p_i;

end

