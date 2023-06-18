function [U_k,J,Xhat]=NMPC2(X_k,Xr_k,Q_bar,R_bar,U0,N,tm,Nm,Ob)

% Xe_k 前一时刻误差系统状态值
% Xr_k 前一时刻期望目标值
% Q_bar R_bar 优化权值矩阵
% U0 前一时刻最优控制序列
% N 预测步数
% tm 采样时间
% Nm 采样间隔

save Weight Q_bar R_bar
save Xk2_data X_k
save Xdk2_data Xr_k
save NN N
save Tm tm
save Obstacle Ob
Nrm = length(Xr_k(:,1));
save Nrm2 Nrm;
%基于实际情况，给输入加约束
I6 = diag(ones(2*N,1));
A = [I6;-I6];
b = [repmat([2 1.5]',N,1);-repmat([-2 -1.5]',N,1)];% 执行器幅值约束
persistent J_1 Uk_1
if isempty(J_1)
    J_1 = 0;
    Uk_1 = zeros(2*N,1);
end
% 采样时刻求解最优控制序列
persistent KK
if isempty(KK)
    KK = 1;
end
save KK2 KK
Aeq = [];
beq = [];
lb = [];
ub =[];
if Nm==0
    [Uk_1,J_1]=fmincon('costfun2',U0,A,b,Aeq,beq,lb,ub,'confun2'); %非线性求解器
    KK = KK+1;
end
J = J_1;
U_k = Uk_1;

% 计算预测轨迹
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
persistent Xhat_1
if isempty(Xhat_1)
    Xhat_1 = repmat(X_k,N+1,1);
end
fxk(1:6,1) = zeros(6,1);
Xhat_1(1:6,1) = X_k;
% 预测未来轨迹
if Nm == 0
    for k=1:N
        [fxk(6*k+1:6*k+6,1)] = fX( Xhat_1(6*(k-1)+1:6*(k-1)+6,1));
        Xhat_1(6*k+1:6*k+6,1) = tm*(fxk(6*k+1:6*k+6)+M*P*U_k(2*(k-1)+1:2*(k-1)+2,1))+Xhat_1(6*(k-1)+1:6*(k-1)+6);
    end
end

Xhat = Xhat_1;

end



