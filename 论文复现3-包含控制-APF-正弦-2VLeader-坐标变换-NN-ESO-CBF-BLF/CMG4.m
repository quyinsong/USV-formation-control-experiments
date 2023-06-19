function [ pio,v_ic] = CMG4( eF,pj_dot,vs,di,p_i0,ts )
%CMG1 此处显示有关此函数的摘要
%   此处显示详细说明

persistent p_i v_i Kv
if isempty(p_i)
    p_i = p_i0;
    v_i = [0 0]';
    Kv = diag([0.2 0.2]);
end
v_ic = (-Kv*eF+pj_dot+vs)/di;

if norm(v_ic) >= 0.5
    v_ic = 0.5*v_ic/norm(v_ic);
end

e_vi = v_i-v_ic;

persistent v_if
if isempty(v_if)
    v_if = v_ic;
end
v_if_dot = -(v_if-v_ic)/0.1;
v_if = v_if_dot*ts+v_if;

u_i = -Kv*e_vi+v_if_dot;
p_i_dot = v_i;
v_i_dot = u_i;
p_i = p_i_dot*ts+p_i;
v_i = v_i_dot*ts+v_i;

pio = p_i;

end

