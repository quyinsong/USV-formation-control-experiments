function [U_k,J,Xhat]=NMPC2(X_k,Xr_k,Q_bar,R_bar,U0,N,tm,Nm,Ob)

% Xe_k ǰһʱ�����ϵͳ״ֵ̬
% Xr_k ǰһʱ������Ŀ��ֵ
% Q_bar R_bar �Ż�Ȩֵ����
% U0 ǰһʱ�����ſ�������
% N Ԥ�ⲽ��
% tm ����ʱ��
% Nm �������

save Weight Q_bar R_bar
save Xk2_data X_k
save Xdk2_data Xr_k
save NN N
save Tm tm
save Obstacle Ob
Nrm = length(Xr_k(:,1));
save Nrm2 Nrm;
%����ʵ��������������Լ��
I6 = diag(ones(2*N,1));
A = [I6;-I6];
b = [repmat([2 1.5]',N,1);-repmat([-2 -1.5]',N,1)];% ִ������ֵԼ��
persistent J_1 Uk_1
if isempty(J_1)
    J_1 = 0;
    Uk_1 = zeros(2*N,1);
end
% ����ʱ��������ſ�������
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
    [Uk_1,J_1]=fmincon('costfun2',U0,A,b,Aeq,beq,lb,ub,'confun2'); %�����������
    KK = KK+1;
end
J = J_1;
U_k = Uk_1;

% ����Ԥ��켣
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
% Ԥ��δ���켣
if Nm == 0
    for k=1:N
        [fxk(6*k+1:6*k+6,1)] = fX( Xhat_1(6*(k-1)+1:6*(k-1)+6,1));
        Xhat_1(6*k+1:6*k+6,1) = tm*(fxk(6*k+1:6*k+6)+M*P*U_k(2*(k-1)+1:2*(k-1)+2,1))+Xhat_1(6*(k-1)+1:6*(k-1)+6);
    end
end

Xhat = Xhat_1;

end



