function [ Q_bar, R_bar ] = calWeight( Q,R,F,N )
%CALWEIGHT 此处显示有关此函数的摘要
%   此处显示详细说明

% 计算权值矩阵
S_q = size(Q,1);%得到Q矩阵维度
S_r = size(R,1);%得到R矩阵维度
Q_bar = zeros((N+1)*S_q,(N+1)*S_q);%定义Q_bar矩阵维度
R_bar = zeros(N*S_r,N*S_r);%定义R_bar矩阵维度
for i = 0:N-1
    Q_bar(i*S_q+1:(i+1)*S_q,i*S_q+1:(i+1)*S_q) = Q;%把对角线上写满Q
end
Q_bar(N*S_q+1:(N+1)*S_q, N*S_q+1:(N+1)*S_q) = F;%最后一块写上F
for i = 0:N-1
    R_bar(i*S_r+1:(i+1)*S_r, i*S_r+1:(i+1)*S_r) = R;%对角线上写满R
end

end

