function pi_o = APF_O(USVp,od)
%APF_O �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

% ���Ϊ�ݶȷ���
n = length(od(1,:)); % ��⵽�ϰ���ĸ���
O = od(1:2,:);
Ro_up = od(3,:);
Ro_dn = od(4,:);
pi_o = zeros(2,1);
for i=1:n
    pi_o = pi_o+(4*(Ro_up(i)^2-Ro_dn(i)^2)*(norm(USVp-O(:,i))^2-Ro_up(i)^2))*(USVp-O(:,i))/(norm(USVp-O(:,i))^2-Ro_dn(i)^2)^3;
end

