function [ pio,v_ic] = CMG1( eF,pj_dot,vs,di,p_i0,Ob,ts )
%CMG1 此处显示有关此函数的摘要
%   此处显示详细说明
% Ob = [xo1 yo1 vox1 voy1 uox1 uoy1;
%       xo2 yo2 vox2 voy2 uox1 uox2;
%       ...];

% 得到标称控制率
persistent p_i v_i Kv Ku
if isempty(p_i)
    p_i = p_i0;
    v_i = [0 0]';
    Kv = diag([2 2]);
    Ku = diag([5 5]);
end
v_ic = (-Kv*eF+pj_dot+vs)/di;

if norm(v_ic) >= 0.6
    v_ic = 0.6*v_ic/norm(v_ic);
end

e_vi = v_i-v_ic;

persistent v_if
if isempty(v_if)
    v_if = v_ic;
end
v_if_dot = -(v_if-v_ic)/0.1;
v_if = v_if_dot*ts+v_if;

u_i = -Ku*e_vi+v_if_dot;

% if norm(u_i) >= 2
%     u_i = 2*u_i/norm(u_i);
% end

% CBF
um = 2;
Ds = 2;
alpha = 10;
k = 5;
if ~isempty(Ob)
    for j = 1:length(Ob(:,1))
        Delta_ri(:,j) = p_i-Ob(j,1:2)';
        Delta_vi(:,j) = v_i-Ob(j,3:4)';
        hi(j) = sqrt(4*um*(norm(Delta_ri(:,j))-Ds))+Delta_ri(:,j)'*Delta_vi(:,j)/norm(Delta_ri(:,j));
        bi(j,1) = alpha*hi(j)^3*norm(Delta_ri(:,j))-k*Delta_ri(:,j)'*Delta_vi(:,j)-...
            (Delta_vi(:,j)'*Delta_ri(:,j))^2/norm(Delta_ri(:,j))^2+norm(Delta_vi(:,j))^2+...
            2*um*Delta_vi(:,j)'*Delta_ri(:,j)/sqrt(4*um*(norm(Delta_ri(:,j))-Ds));
        Ai(j,:) = -Delta_ri(:,j)';
    end
    lb = 4*[-1 1]'; ub = 4*[1 1]';
    H = diag([1 1]); f = -2*u_i'*H;
    [u_ic,fval,exitflag,output,lambda] = quadprog(H,f',Ai,bi,[],[],lb,ub);
else 
    u_ic = u_i;
end

% 更新状态
p_i_dot = v_i;
v_i_dot = u_ic;
p_i = p_i_dot*ts+p_i;
v_i = v_i_dot*ts+v_i;

pio = p_i;


end

















