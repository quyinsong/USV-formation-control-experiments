function [nu_c,zeta1,zeta2,q,vartheta,Theta1o,Theta2o] = Guidance3( hi, R_psi, zi, pj_dot, nu_bar, t, ts)
%GUIDANCE1 此处显示有关此函数的摘要
%   此处显示详细说明

uhat = nu_bar(1);
vhat = nu_bar(2);
rhat = nu_bar(3);

persistent Ki Theta1 Theta2
if isempty(Ki)
    Ki = diag([0.5 0.5]);
    Theta1 = 0;
    Theta2 = 0;
end

% 预设性能
a1 = 5; b1 = 5; zeta1 = (1-0.4)*exp(-0.1*t)+0.4;
a2 = 5; b2 = 5; zeta2 = (1-0.4)*exp(-0.2*t)+0.4;
zeta1_dot = -0.2*(1-0.4)*exp(-0.2*t);
zeta2_dot = -0.2*(1-0.4)*exp(-0.2*t);

z1 = zi(1); z2 = zi(2);

% 误差变换
q1 = 0.5*log((b1*z1+a1*b1*zeta1)/(a1*b1*zeta1-a1*z1));
q2 = 0.5*log((b2*z2+a2*b2*zeta2)/(a2*b2*zeta2-a2*z2));

% 制导率
vartheta1 = 0.5*(1/(z1+a1*zeta1)-1/(z1-b1*zeta1));
vartheta2 = 0.5*(1/(z2+a2*zeta2)-1/(z2-b2*zeta2));
varrho=R_psi'*pj_dot-[0;vhat];
uc = varrho(1)+zeta1_dot*z1/zeta1-0.5*q1/vartheta1-...
    0.5*vartheta1*q1-Theta1*sign(vartheta1*q1);
rc = (varrho(2)+zeta2_dot*z2/zeta2-0.05*q2/vartheta2-...
    0.05*vartheta2*q2-Theta2*sign(vartheta2*q2))/hi(2,2);

% uc = varrho(1)-0.2*z1;
% rc = (varrho(2)-0.2*z2)/hi(2,2);

q = [q1,q2]';
vartheta = [vartheta1,vartheta2]';

% 自适应参数更新
Theta1_dot = 5*(-0.5*Theta1+abs(vartheta1*q1));
Theta2_dot = 10*(-0.5*Theta2+abs(vartheta2*q2));

Theta1 = Theta1_dot*ts+Theta1;
Theta2 = Theta2_dot*ts+Theta2;

% nu_c = hi\(-Ki*zi+R_psi'*pj_dot-[0;vhat]);
nu_c=[uc rc]';
% if nu_c(1)<=0.01
%     nu_c(1)=0.01;
% end
Theta1o = Theta1;
Theta2o = Theta2;

end

