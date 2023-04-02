function [ y1, y2, y1_d, y2_d, se, f, fhat ] = kinematiclaw1( s, sc, psi, nu_j, psi_j, nu_i, ts)
%KINEMATICLAW1 [ output_args ] = kinematiclaw1( s, sd, ts, t )
%   calculate desired kinematic control law
% Date: 7 13 2022
% inputs : s= [rho lambda]' : LOS range and angle
%          sd = [rhod lambdad]' : desired LOS range and angle
%          psi : heading angle
%          nu_j = [uj vj rj]' : speed of leader
%          psi_j : heading angle of leader
%          ts : sample time
% states: varrho_hat : adaptive law
% outputs: y1: ud : desire surge speed
%          y2: rd : desired yaw angular rate
%          y1_d: ud_d : the derivative of ud
%          y2_d: rd_d : the derivative of rd
%          se = s-sd : formation errors (LOS range and angle errors)

persistent varrho_hat psid ud rd sd
if isempty(varrho_hat)
   varrho_hat = 0; 
   psid = 0;
   ud = 0;
   rd = 0;
   sd = [0 0]';
end
% parameters
delta1 = 1; delta2 = 1;
Gamma_varrho = 100;
K_varrho = 0.5;
Ks = diag([1 1]);
C = diag([ 1/10 1/10]);% 1st order filter
gamma1 = 0.05; gamma2 = 0.05; gamma3 = 0.05; % 1st order filter
kr = 2;

% 1st order filter
sd_dot = C*(sc-sd);
sd = euler2(sd_dot,sd,ts);

se = s-sd; % error vector
rho = s(1);
lambda = s(2);
B = diag([1 1/rho]);
Y = diag([tanh(se'*B(:,1))/delta1, tanh(se'*B(:,2))/delta2]);
% adaptive law
varrho_hat_d = Gamma_varrho*(Y*B'*se-K_varrho*varrho_hat);
varrho_hat = euler2(varrho_hat_d,varrho_hat,ts);
% kinematic law
%--------假设Delta已知--------------
u_j = nu_j(1); v_j = nu_j(2); v_i = nu_i(2); u_i = nu_i(1);
Delta = [v_i*cos(psi_j+lambda)-u_j*sin(psi_j+lambda)+v_j*cos(psi_j+lambda)
         -v_i*sin(psi_j+lambda)-u_j*cos(psi_j+lambda)+v_j*cos(psi_j+lambda)];
f = B*Delta;
% w_bar = B\(sd_dot-Ks*se-Delta);
%--------------end-- ---------------
%--------假设Delta未知--------------
% w_bar = B\(sd_dot-Ks*se)-Y*varrho_hat;
%--------假设Delta未知  状态观测器估计--------------
w = [u_i*sin(psi+lambda) u_i*cos(psi+lambda)]';
persistent s_hat f_hat ksh1 ksh2
if isempty(s_hat)
   s_hat = s;
   f_hat = zeros(2,1);
   ksh1=  diag([5 5]);
   ksh2 = diag([1 1]);
end
s_hat_dot = B*w+f_hat-ksh1*(s_hat-s);  %f_hat = B*delta_hat
f_hat_dot = -ksh2*(s_hat-s);
s_hat = s_hat_dot*ts+s_hat;
f_hat = ts*f_hat_dot+f_hat;
fhat = f_hat;

w_bar = B\(sd_dot-Ks*se)-f_hat;
%--------------end-- ---------------
% virtual control signal
psi_bar = atan2(w_bar(1),w_bar(2))-lambda;
u_bar = sin(psi_bar+lambda)*w_bar(1)+cos(psi_bar+lambda)*w_bar(2);
% 1st order filter

psid_dot = (psi_bar-psid)/gamma1;
psid = euler2(psid_dot,psid,ts);
ud_dot = (u_bar-ud)/gamma2;
ud = euler2(ud_dot,ud,ts);
% virtual yaw control law
psie = ssa(psi-psid);

r_bar = psid_dot-kr*psie;
% 1st order filter
rd_dot = (r_bar-rd)/gamma3;
rd = euler2(rd_dot,rd,ts);
% outputs
y1 = ud;
y2 = rd;
y1_d = ud_dot;
y2_d = rd_dot;

end

