function [nu_c,pe] = Guidance2( p,psi,pd,pd_dot,alpha,ts)
%GUIDANCE1 此处显示有关此函数的摘要
%   此处显示详细说明

pe = p-pd;
k= 2;
F = k*(alpha'*pe)*pe-alpha*(pe'*pe);
psid = atan2(F(2),F(1));
if F(2)==0&&F(1)==0
    psid = psi;
end
persistent psidf
if isempty(psidf)
    psidf = psid;
end
psidf_dot = -(psidf-psid)/0.1;
psidf = ts*psidf_dot+psidf;

ku = 0.1;
kr = 1;
uc = ku*norm(pe)/sqrt(norm(pe)^2+2)+norm(pd_dot);
psie = psi- psid;
rc = -kr*psie/sqrt(psie^2+2)+psidf_dot;
nu_c = [uc,rc]';

end

