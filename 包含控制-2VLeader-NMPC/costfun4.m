function f = costfun4( U )
%COSTFUN 此处显示有关此函数的摘要
%   此处显示详细说明
load Xk4_data
load Xdk4_data
load NN
load Tm 
load KK4
Xk_p(1:6,1) = X_k;
Xd = Xr_k;
fxk(1:6,1) = zeros(6,1);
%---------------------------------------------
P = [ 1 0; 0 0; 0 1; 0 0; 0 0; 0 0];
m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
Nrdot = -1; xg = 0.046; Yrdot = 0;          
m11 = m-Xudot; 
m22 = m-Yvdot;
m23 = m*xg-Yrdot;
m32 = m*xg-Nvdot;
m33 = Iz-Nrdot;
M=[m11  0    0;
   0   m22 m23;
   0   m32 m33];
M = [inv(M),zeros(3,3);
     zeros(3,6)];
%-----------------------------------
for k=1:N
    [fxk(6*k+1:6*k+6,1)] = fX( Xk_p(6*(k-1)+1:6*(k-1)+6,1));
    Xk_p(6*k+1:6*k+6,1) = tm*(fxk(6*k+1:6*k+6)+M*P*U(2*(k-1)+1:2*(k-1)+2,1))+Xk_p(6*(k-1)+1:6*(k-1)+6);
end

load Weight
load Nrm4

if 6*(KK-1)+6*(N+1)>Nrm
    if 6*(KK-1)+6<=Nrm
        Xdk(1:Nrm-6*(KK-1),1) = Xd(6*(KK-1)+1:Nrm);
        Xdk(Nrm-6*(KK-1)+1:6*(N+1),1) = repmat(Xd(Nrm-5:Nrm,1),(N+1)-(Nrm/6-(KK-1)),1);
    else
        Xdk(1:6*(N+1),1) = repmat(Xd(Nrm-5:Nrm,1),N+1,1);
    end
else
    Xdk = Xd(6*(KK-1)+1:6*(KK-1)+6*(N+1),1);
end
    
f = (Xk_p-Xdk)'*Q_bar*(Xk_p-Xdk)+U'*R_bar*U;

end










