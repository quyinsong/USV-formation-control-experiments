function [Xm,Nm]  = Navigation4( X, Tm, Ts )
%STATES �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
persistent N Xm_pre
if isempty(N)
    N = Tm/Ts;
    Xm_pre = X;
end
N = N-1;
Nm = N;
if Nm==0
    Xm_pre = X;
    N = Tm/Ts;
end
Xm = Xm_pre;
end

