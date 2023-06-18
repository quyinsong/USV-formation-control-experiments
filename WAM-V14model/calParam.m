clc
clear all
close all
L = 4.29;
LWL = 3.21;
T = 0.127;
LGG = 1.27;
BOA = 2.2;
B = 1.83;
Bhull = BOA-B;
rho = 10^3;
m = 150;
Cd = 1.1;

Nvdot = -2.5*pi*rho*T^2*((L-LGG)^2+LGG^2)/2;
Nrdot = -1.2*(4.75*pi*rho*Bhull*T^4/4+pi*rho*T^2*((L-LGG)^3+LGG^3)/3);
Xudot = -0.075*m;
Yrdot = -0.2*pi*rho*T^2*((L-LGG)^2+LGG^2)/2;
Yvdot = -0.9*pi*rho*T^2*L;

Xu = -50.897;

Xuu = -5.8722;
Yvv = -rho*T*Cd*L;
Yvr = -rho*T*1.1*((L-LGG)^2-LGG^2)/2;
Yrv = Yvr;
Yrr = -rho*T*Cd*((L-LGG)^3+LGG^3)/3;
Nvv = Yvr;
Nvr = Yrr;
Nrv = Yrr;
Nrr = -rho*T*Cd*((L-LGG)^4+LGG^4)/4;











