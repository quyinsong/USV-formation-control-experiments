function pi_usv = APF_usv(USVp,usvd,Rusv_up,Rusv_dn)
%APF_O �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

% ���Ϊ�ݶȷ���
n = length(usvd(1,:)); % ��⵽�ϰ���ĸ���
pi_usv = [0 0]';
for i=1:n
    pi_usv = pi_usv+(4*(Rusv_up^2-Rusv_dn^2)*(norm(USVp-usvd(:,i))^2-Rusv_up^2))*...
             (USVp-usvd(:,i))/(norm(USVp-usvd(:,i))^2-Rusv_dn^2)^3;
end