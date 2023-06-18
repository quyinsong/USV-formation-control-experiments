function [Uc_k,J,Xhat]=LMPCnavi4(x_k,Xr_k,Q_bar,R_bar,Urk,N,tm,Nm)
%--------
% Xk = M*Xk_1+C*Uk_1
%%%%%%%%%%%%%%%%%%%%%%%%
n = size(x_k,1); %% A������n * n���󣬵õ�A�����ά��
p = 2; %% B������n * p���󣬵õ�B�����ά��
Nrxm = length(Xr_k(:,1)); % ��ȡ���й켣���������ܳ���
Nrum = length(Urk(:,1)); % ��ȡ���й켣���������ܳ���
% ÿ�β���ʱ�̶���ǰȡ6*(N+1)�������켣��,����һ����ʱ��Ԥ��ʱ������Ҫ���ٵĹ켣��
persistent KK
if isempty(KK)
    KK = 1;
end
Xd = Xr_k; %��ȡ���й켣������
% ȡ��ǰ��Ҫ�Ĺ켣�����ݺͲο���������
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
%����������������M������C���󣬺�f����
tmp = eye(n); %����һ��n�׵�λ�󣬹�����
M = [eye(n);zeros(N*n,n)]; %% ��ʼ��M����,��һ���ֿ�����õ�λ�������������
C = zeros((N+1)*n,N*p); %% ��ʼ��C��������
xd_k = Xdk(1:6,1);
for i = 1:N
    [Ak,Bk] = calAkBk( x_k, tm ); % ���㵱ǰʱ���������켣�����Ի���ϵͳ����
%     xd_k = Xdk(6*i+1:6*i+6,1); % ��ȡԤ��ʱ���ڵ�iʱ�̵Ĺ켣��
    rows = i*n + (1:n);%��������Ϊ�Ƿֿ�������Դ�1��n;
    C(rows, :) = [tmp*Bk, C(rows-n, 1:end-p)];%�ñ����ķ�����C��������;
    tmp = Ak*tmp;%ÿ�ζ����һ��A����;
    M(rows,:) = tmp;%д��M����;
end

G = M'*Q_bar*M;%����M������ʵ���ڴ��ۺ����У���������޹أ���û�б��õ�
E = M'*Q_bar*C;%����E����
H = C'*Q_bar*C + R_bar;%����H����
f = x_k'*E-Xdk'*Q_bar*C;%����quadprog�����Ķ��壬��Ҫ����д�ɾ��������ʽ
% ��⵱ǰ����ʱ�̺�Ԥ��ʱ���ڵĴ��ۺ���
% ����Լ��
b1 = 1000*[3 3 2 0.5 0.6 0.6]'; b1rep = repmat(b1,N+1,1); % ��������ֵԼ��
b2 = -1000*[3 3 2 0.5 0.6 0.6]'; b2rep = repmat(b2,N+1,1); % ��������ֵԼ��
b3 = [(b1rep-M*x_k)',(M*x_k-b2rep)']';
D = [C ;-C];
Aep=[];Bep=[];
lb=[-2 -1.5]';ub=[2 1.5]'; % ִ������ֵԼ��
lb = repmat(lb,N,1); ub = repmat(ub,N,1);
% ��ǰ����ʱ�̺�Ԥ��ʱ����������ϵͳ���ſ�������
persistent Uk_pre J_pre
if isempty(Uk_pre)
    Uk_pre = zeros(2*N,1);
    J_pre = 0;
end
if Nm == 0
    Uk_pre = quadprog(H,f',[],[],Aep,Bep,lb,ub);%������ŵ�U_kֵ
    J_pre = x_k'*G*x_k+2*x_k'*E*Uk_pre+Uk_pre'*H*Uk_pre-2*Xdk'*Q_bar*C*Uk_pre;
    KK = KK+1;
end
J = J_pre;
U_k = Uk_pre;
% Ԥ�����ϵͳδ���켣
persistent Xehat_pre
if isempty(Xehat_pre)
    Xehat_pre = M*x_k+C*U_k;
end
if Nm == 0
    Xehat_pre = M*x_k+C*U_k;
end
Xehat = Xehat_pre; % Ԥ������ϵͳ״̬
% �õ���ʵ��Ԥ��켣
Xhat = Xehat; 

% �����ʵ�����ſ�������
Uc_k = U_k;


end




