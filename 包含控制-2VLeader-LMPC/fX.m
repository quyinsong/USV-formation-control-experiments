function fxk = fX( Xk )
%FX 此处显示有关此函数的摘要
%   此处显示详细说明

uk = Xk(1);
vk = Xk(2);
rk = Xk(3);
xk = Xk(4);
yk = Xk(5);
psik = Xk(6);

nuk = [uk,vk,rk]';
%% USV parameters 
% m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
% Nrdot = -1; xg = 0.046; Yrdot = 0;
m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
Nrdot = -1; xg = 0.046; Yrdot = 0;
% ------------------------------------------------------
Xu=-0.72253;         Yv=-0.88965;          Nv=0.1052;
Xuu=-1.32742;        Yr=0.1079;            Nr=-1.900;
Xuuu = -5.8664;      Yvv=-36.47287;        Nvv=5.0437;
                     Yrv=-0.805;           Nrv=0.130;
                     Yvr=-0.845;           Nvr=0.080;
                     Yrr=-3.45;            Nrr=-0.75;               
% ----------------------------------------------------
m11 = m-Xudot; 
m22 = m-Yvdot;
m23 = m*xg-Yrdot;
m32 = m*xg-Nvdot;
m33 = Iz-Nrdot;
% -----------------------------------------------------
c13=-m23*rk-m22*vk; c23 = m11*uk;
c31 = -c13; c32 = -c23;
% -----------------------------------------------------
d11=-Xu-Xuu*abs(uk);
d22=-Yv-Yvv*abs(vk)-Yrv*abs(rk);
d23=-Yr-Yvr*abs(vk)-Yrr*abs(rk);
d32=-Nv-Nvv*abs(vk)-Nrv*abs(rk);
d33=-Nr-Nvr*abs(vk)-Nrr*abs(rk);


%% matrix expression

Rbnk=[ cos(psik) -sin(psik) 0;
      sin(psik) cos(psik)  0;
      0          0         1];
M=[m11  0    0;
   0   m22 m23;
   0   m32 m33];
Crbk=[0     0     c13;
     0     0     c23;
    c31   c32     0 ];
Dvk=[d11  0     0;
     0    d22  d23;
     0    d32  d33];
% next uotputs

fxk = [M\(-Crbk*nuk-Dvk*nuk);
       Rbnk*nuk] ;


end

