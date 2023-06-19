function nu_c = Guidance3( hi, R_psi, zi, pj_dot, nu_bar)
%GUIDANCE1 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

persistent Ki
if isempty(Ki)
    Ki = diag([0.5 0.5]);
end

nu_c = hi\(-Ki*zi+R_psi'*pj_dot-[0;nu_bar(2)]);
if nu_c(1)<=0.01
    nu_c(1)=0.01;
end

end

