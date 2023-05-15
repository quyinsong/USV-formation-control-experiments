function [ pio,pi_dot] = CMG4( eF,pj_dot,vs,di,p_i0,pi_o, pi_usv,ts )
%CMG1 此处显示有关此函数的摘要
%   此处显示详细说明

persistent p_i Kg
if isempty(p_i)
    p_i = p_i0;
    Kg = diag([0.2 0.2]);
end
pi_dot = (-Kg*eF+pj_dot+vs)/di-pi_o-pi_usv;

if norm(pi_dot) >= 0.5
    pi_dot = 0.5*pi_dot/norm(pi_dot);
end

p_i = ts*pi_dot+p_i;
pio = p_i;

end

