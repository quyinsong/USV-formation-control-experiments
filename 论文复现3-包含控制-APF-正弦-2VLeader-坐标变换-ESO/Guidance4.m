function nu_c = Guidance4( hi, R_psi, zi2, pj_dot, v_hat)
%GUIDANCE1 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

persistent Ki
if isempty(Ki)
    Ki = diag([0.1 0.05]);
end

nu_c = hi\(-Ki*zi2+R_psi'*pj_dot-[0 v_hat]');
if nu_c(1)<=0.01
    nu_c(1)=0.01;
end

end