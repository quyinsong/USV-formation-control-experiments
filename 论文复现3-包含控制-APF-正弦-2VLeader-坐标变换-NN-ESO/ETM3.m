function [tauc_etm,Tk] = ETM3( tauc,emiu,emir,du,dr )
%ETM1 此处显示有关此函数的摘要
%   此处显示详细说明
% tauc_u_etm tauc_r_etm 为零阶保持器保持的状态，在需要更新的时候才会被更新
persistent tauc_u_etm tauc_r_etm
if isempty(tauc_u_etm)
    tauc_u_etm = tauc(1);
    tauc_r_etm = tauc(2);
end

tuk = 0;
trk = 0;

Xu = tauc_u_etm-tauc(1);
Xr = tauc_r_etm-tauc(2);

if abs(Xu)>=emiu*abs(tauc_u_etm)+du
    tauc_u_etm = tauc(1);
    tuk = 1;
end

if abs(Xr)>=emir*abs(tauc_r_etm)+dr
    tauc_r_etm = tauc(2);
    trk = 1;
end

tauc_etm = [tauc_u_etm,tauc_r_etm]';
Tk = [tuk,trk]';

end

