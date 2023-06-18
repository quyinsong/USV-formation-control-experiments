function [ Q_bar, R_bar ] = calWeight( Q,R,F,N )
%CALWEIGHT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

% ����Ȩֵ����
S_q = size(Q,1);%�õ�Q����ά��
S_r = size(R,1);%�õ�R����ά��
Q_bar = zeros((N+1)*S_q,(N+1)*S_q);%����Q_bar����ά��
R_bar = zeros(N*S_r,N*S_r);%����R_bar����ά��
for i = 0:N-1
    Q_bar(i*S_q+1:(i+1)*S_q,i*S_q+1:(i+1)*S_q) = Q;%�ѶԽ�����д��Q
end
Q_bar(N*S_q+1:(N+1)*S_q, N*S_q+1:(N+1)*S_q) = F;%���һ��д��F
for i = 0:N-1
    R_bar(i*S_r+1:(i+1)*S_r, i*S_r+1:(i+1)*S_r) = R;%�Խ�����д��R
end

end

