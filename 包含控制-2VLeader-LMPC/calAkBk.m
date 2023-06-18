function [Ak,Bk] = calAkBk( xr_k, tm )
%CALAK 此处显示有关此函数的摘要
%   此处显示详细说明
% desired trajectory points
urk = xr_k(1);
vrk = xr_k(2);
rrk = xr_k(3);
psirk = xr_k(6);
% model parameters
m11 = 25.8; m22 = 33.8; m33 = 2.76;
Xu=-0.72253; Yv=-0.88965; Nr=-1.900;
afuu = Xu/m11; afuv = m22*rrk/m11; afur = m22*vrk/m11;
afvu = -m11*rrk/m22; afvv = Yv/m22; afvr = -m11*urk/m22;
afru = (m11-m22)*vrk/m33; afrv = (m11-m22)*urk/m33; afrr = Nr/m33;

% calculate system matrix
n = size(xr_k,1);
Al = [afuu afuv afur 0 0 0;
      afvu afvv afvr 0 0 0;
      afru afrv afrr 0 0 0;
      cos(psirk) -sin(psirk) 0 0 0 -urk*sin(psirk)-vrk*cos(psirk);
      sin(psirk) cos(psirk) 0 0 0 urk*cos(psirk)-vrk*sin(psirk);
      0 0 1 0 0 0];
I = diag(ones(n,1));
Ak = I+tm*Al;
  
Bl = [1/m11 0; 0 0; 0 1/m33; 0 0; 0 0; 0 0];
Bk = tm*Bl;
  
  
end

