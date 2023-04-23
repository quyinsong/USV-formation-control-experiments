function rc = GVF4( x,y,psi,beta,rc,U )
%GVF �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

% constraints
X = psi+beta;
m = [cos(X), sin(X)]';
E = [0 -1;1 0]; I2 = diag([1 1]);
kro = 0.05;
Ud = 0.4;
% boundary 1
R1 = 41; % �����߽磬����ñ߽����ײ��Ϊ�͸�����Ϊͬʱ����
d1 = -2; % �ڲ�߽磬����ñ߽��ֻ���ڱ�����Ϊ
phir1 = 5*sin(0.05*y)+R1-x;
no1 = [-1 0.25*cos(0.05*y)]'; taur1 = E*no1;
H1 = [0 0; 0 -0.25*0.05*sin(0.05*y)];
kod1 = -1;

mod1 = kod1*taur1-kro*phir1*no1;
mod1_hat = mod1/norm(mod1);

mod1_dot = (kod1*E-kro*phir1*I2)*Ud*H1*m-kro*Ud*no1'*m*no1;
mod1_hat_dot = -E*(mod1_hat*mod1_hat')*E*mod1_dot/norm(mod1);
psio1_dot = -mod1_hat'*E*mod1_hat_dot;
r_od1 = psio1_dot-0.5*m'*E*mod1_hat;
       
if phir1 <= d1
    f1 = 0; f2 = 1;
elseif phir1 < 0
    l1 = 5; l2 = 1;  % l1 Խ����������ȼ�Խ�� l2 Խ����·���������ȼ�Խ��
    f1 = exp(l1/(d1-phir1)); f2 = exp(l2/phir1);
else
    f1 = 1; f2 = 0;
end
si1 = f2/(f1+f2); zi1 = f1/(f1+f2);

% boundary 2
R2 = 19; % �����߽磬����ñ߽����ײ��Ϊ�͸�����Ϊͬʱ����
d2 = -2; % �ڲ�߽磬����ñ߽��ֻ���ڱ�����Ϊ
phir2 = x-5*sin(0.05*y)-R2;
no2 = [1 -0.25*cos(0.05*y)]'; taur2 = E*no2;
H2 = [0 0; 0 0.25*0.05*sin(0.05*y)];

kod2 = 1;
mod2 = kod2*taur2-kro*phir2*no2;
mod2_hat = mod2/norm(mod2);

mod2_dot = (kod2*E-kro*phir2*I2)*Ud*H2*m-kro*Ud*no2'*m*no2;
mod2_hat_dot = -E*(mod2_hat*mod2_hat')*E*mod2_dot/norm(mod2);
psio2_dot = -mod2_hat'*E*mod2_hat_dot;
r_od2 = psio2_dot-0.5*m'*E*mod2_hat;
       
if phir2 <= d2
    f1 = 0; f2 = 1;
elseif phir2 < 0
    l1 = 5; l2 = 1;  % l1 Խ����������ȼ�Խ�� l2 Խ����·���������ȼ�Խ��
    f1 = exp(l1/(d2-phir2)); f2 = exp(l2/phir2);
else
    f1 = 1; f2 = 0;
end
si2 = f2/(f1+f2); zi2 = f1/(f1+f2);

rc = si1*r_od1+si2*r_od2+zi1*zi2*rc;









end

