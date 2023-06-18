function [Uc_k,J,Xhat]=LMPCnavi4(x_k,Xr_k,Q_bar,R_bar,Urk,N,tm,Nm)
%--------
% Xk = M*Xk_1+C*Uk_1
%%%%%%%%%%%%%%%%%%%%%%%%
n = size(x_k,1); %% A矩阵是n * n矩阵，得到A矩阵的维数
p = 2; %% B矩阵是n * p矩阵，得到B矩阵的维数
Nrxm = length(Xr_k(:,1)); % 获取所有轨迹点向量的总长度
Nrum = length(Urk(:,1)); % 获取所有轨迹点向量的总长度
% 每次采样时刻都向前取6*(N+1)个期望轨迹点,即下一采样时刻预测时域内需要跟踪的轨迹点
persistent KK
if isempty(KK)
    KK = 1;
end
Xd = Xr_k; %获取所有轨迹点数据
% 取当前需要的轨迹点数据和参考控制数据
if 6*(KK-1)+6*(N+1)>Nrxm
    if 6*(KK-1)+6<=Nrxm
        Xdk(1:Nrxm-6*(KK-1),1) = Xd(6*(KK-1)+1:Nrxm);
        Xdk(Nrxm-6*(KK-1)+1:6*(N+1),1) = repmat(Xd(Nrxm-5:Nrxm,1),(N+1)-(Nrxm/6-(KK-1)),1);
    else
        Xdk(1:6*(N+1),1) = repmat(Xd(Nrxm-5:Nrxm,1),N+1,1);
    end
else
    Xdk = Xd(6*(KK-1)+1:6*(KK-1)+6*(N+1),1);
end

if 2*(KK-1)+2*(N)>Nrum
    if 2*(KK-1)+2<=Nrum
        Urdk(1:Nrum-2*(KK-1),1) = Urk(2*(KK-1)+1:Nrum);
        Urdk(Nrum-2*(KK-1)+1:2*(N),1) = repmat(Urk(Nrum-1:Nrum,1),(N)-(Nrum/2-(KK-1)),1);
    else
        Urdk(1:2*(N),1) = repmat(Urk(Nrxm-1:Nrum,1),N,1);
    end
else
    Urdk = Urk(2*(KK-1)+1:2*(KK-1)+2*(N),1);
end
%接下来计算完整的M矩阵与C矩阵，和f矩阵
tmp = eye(n); %定义一个n阶单位阵，工具人
M = [eye(n);zeros(N*n,n)]; %% 初始化M矩阵,第一个分块矩阵置单位阵，其余矩阵置零
C = zeros((N+1)*n,N*p); %% 初始化C矩阵，置零
xd_k = Xdk(1:6,1);
for i = 1:N
    [Ak,Bk] = calAkBk( x_k, tm ); % 计算当前时刻在期望轨迹点线性化的系统矩阵
%     xd_k = Xdk(6*i+1:6*i+6,1); % 获取预测时域内第i时刻的轨迹点
    rows = i*n + (1:n);%行数，因为是分块矩阵所以从1至n;
    C(rows, :) = [tmp*Bk, C(rows-n, 1:end-p)];%用遍历的方法将C矩阵填满;
    tmp = Ak*tmp;%每次都左乘一次A矩阵;
    M(rows,:) = tmp;%写满M矩阵;
end

G = M'*Q_bar*M;%定义M矩阵，事实上在代价函数中，这和输入无关，并没有被用到
E = M'*Q_bar*C;%定义E矩阵
H = C'*Q_bar*C + R_bar;%定义H矩阵
f = x_k'*E-Xdk'*Q_bar*C;%由于quadprog函数的定义，需要把其写成矩阵相乘形式
% 求解当前采样时刻后预测时域内的代价函数
% 计算约束
b1 = 1000*[3 3 2 0.5 0.6 0.6]'; b1rep = repmat(b1,N+1,1); % 跟踪误差幅值约束
b2 = -1000*[3 3 2 0.5 0.6 0.6]'; b2rep = repmat(b2,N+1,1); % 跟踪误差幅值约束
b3 = [(b1rep-M*x_k)',(M*x_k-b2rep)']';
D = [C ;-C];
Aep=[];Bep=[];
lb=[-2 -1.5]';ub=[2 1.5]'; % 执行器幅值约束
lb = repmat(lb,N,1); ub = repmat(ub,N,1);
% 当前采样时刻后预测时域内求解误差系统最优控制序列
persistent Uk_pre J_pre
if isempty(Uk_pre)
    Uk_pre = zeros(2*N,1);
    J_pre = 0;
end
if Nm == 0
    Uk_pre = quadprog(H,f',[],[],Aep,Bep,lb,ub);%求解最优的U_k值
    J_pre = x_k'*G*x_k+2*x_k'*E*Uk_pre+Uk_pre'*H*Uk_pre-2*Xdk'*Q_bar*C*Uk_pre;
    KK = KK+1;
end
J = J_pre;
U_k = Uk_pre;
% 预测误差系统未来轨迹
persistent Xehat_pre
if isempty(Xehat_pre)
    Xehat_pre = M*x_k+C*U_k;
end
if Nm == 0
    Xehat_pre = M*x_k+C*U_k;
end
Xehat = Xehat_pre; % 预测的误差系统状态
% 得到真实的预测轨迹
Xhat = Xehat; 

% 求解真实的最优控制序列
Uc_k = U_k;


end




