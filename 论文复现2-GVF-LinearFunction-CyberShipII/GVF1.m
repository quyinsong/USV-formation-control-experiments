function rc = GVF1( x,y,psi,beta,rc )
%GVF 此处显示有关此函数的摘要
%   此处显示详细说明

% constraints
X = psi+beta;
m = [cos(X), sin(X)]';
E = [0 -1;1 0]; I2 = diag([1 1]);
kro = 0.05;
Ud = 0.3;
% boundary 1
R1 = 21; % 最外层边界，进入该边界后碰撞行为和跟踪行为同时存在
d1 = -2; % 内层边界，进入该边界后只存在避碰行为
phir1 = y+R1-x;
no1 = [-1 1]'; taur1 = E*no1;
H1 = [0 0; 0 0];
kod = -1;

mod1 = kod*taur1-kro*phir1*no1;
mod1_hat = mod1/norm(mod1);

mc1 = Ud*m;
mc1_hat = mc1/norm(mc1);
       
mod1_dot = (E-kro*phir1*I2)*H1*mc1-kro*no1'*mc1*no1;
mod1_hat_dot = -E*(mod1_hat*mod1_hat')*E*mod1_dot/norm(mod1);
psio1_dot = -mod1_hat'*E*mod1_hat_dot;
r_od1 = (psio1_dot-0.5*mc1_hat'*E*mod1_hat)*norm(mc1)/(Ud*mc1_hat'*m);
       
if phir1 <= d1
    f1 = 0; f2 = 1;
elseif phir1 < 0
    l1 = 100; l2 = 5;  % l1 越大则避障优先级越高 l2 越大则路径跟踪优先级越高
    f1 = exp(l1/(d1-phir1)); f2 = exp(l2/phir1);
else
    f1 = 1; f2 = 0;
end
si1 = f2/(f1+f2); zi1 = f1/(f1+f2);

% boundary 2
R2 = 18; % 最外层边界，进入该边界后碰撞行为和跟踪行为同时存在
d2 = -3; % 内层边界，进入该边界后只存在避碰行为
phir2 = x-y+R2;
no2 = [1 -1]'; taur2 = E*no1;
H2 = [0 0; 0 0];

mod2 = kod*taur2-kro*phir2*no2;
mod2_hat = mod2/norm(mod2);

mc2 = Ud*m;
mc2_hat = mc2/norm(mc2);
       
mod2_dot = (E-kro*phir2*I2)*H2*mc2-kro*no2'*mc2*no2;
mod2_hat_dot = -E*(mod2_hat*mod2_hat')*E*mod2_dot/norm(mod2);
psio2_dot = -mod2_hat'*E*mod2_hat_dot;
r_od2 = (psio2_dot-0.5*mc2_hat'*E*mod2_hat)*norm(mc2)/(Ud*mc2_hat'*m);
       
if phir2 <= d2
    f1 = 0; f2 = 1;
elseif phir2 < 0
    l1 = 100; l2 = 5;  % l1 越大则避障优先级越高 l2 越大则路径跟踪优先级越高
    f1 = exp(l1/(d2-phir1)); f2 = exp(l2/phir1);
else
    f1 = 1; f2 = 0;
end
si2 = f2/(f1+f2); zi2 = f1/(f1+f2);

rc = si1*r_od1+si2*r_od2+zi1*zi2*rc;








end

