function [Xm,Nm]  = Navigation4( X, Tm, Ts )
%STATES 此处显示有关此函数的摘要
%   此处显示详细说明
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

