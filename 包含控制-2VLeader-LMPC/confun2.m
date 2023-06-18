function [gc,hc] = confun2( U )
%CONFUN 此处显示有关此函数的摘要
%   此处显示详细说明
load Xk2_data
load NN
load Tm 
load Obstacle

Xk_p(1:6,1) = X_k;
fxk(1:6,1) = zeros(6,1);
%---------------------------------------------
P = [ 1 0; 0 0; 0 1; 0 0; 0 0; 0 0];
m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
Nrdot = -1; xg = 0.046; Yrdot = 0;          
m11 = m-Xudot; 
m22 = m-Yvdot;
m23 = m*xg-Yrdot;
m32 = m*xg-Nvdot;
m33 = Iz-Nrdot;
M=[m11  0    0;
   0   m22 m23;
   0   m32 m33];
M = [inv(M),zeros(3,3);
     zeros(3,6)];
%-----------------------------------
for k=1:N
    [fxk(6*k+1:6*k+6,1)] = fX( Xk_p(6*(k-1)+1:6*(k-1)+6,1));
    Xk_p(6*k+1:6*k+6,1) = tm*(fxk(6*k+1:6*k+6)+M*P*U(2*(k-1)+1:2*(k-1)+2,1))+Xk_p(6*(k-1)+1:6*(k-1)+6);
end

load Weight

% 判断与障碍物之间的距离
% dis = norm(X_k(4:5)-Ob(1:2));
% if dis<=Ob(3)
%     for k = 1:N
%         gc(k,1) = -norm(Xk_p(6*k+4:6*k+5,1)-Ob(1:2,1))+Ob(4);
%         hc = [];
%     end
% else
%     gc = [];
%     hc = [];
% end
gc = [];
hc = [];
end

