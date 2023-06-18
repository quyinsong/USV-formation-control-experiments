function  modelplot( pos, psi, color,linewid )
% MODELPLOT modelplot(pos,psai) :plot the geometrical model of ship in real time
%
% Input:  pos; ship real time position in {n}, pos=[x y]'; 
%         psai: real time yaw angle
% Output: none
%
% Auther: Quyinsong
% Data:   12nd Jan 2022
% new version: update by Quyinsong   Date: 2022 5 28
%% check of input dimensions
if nargin~=4 
    error('the number of input must be 6');end
if length(pos)~=2 
    error('the number of pos must be 2');end
if length(psi)~=1 
    error('the number of psai must be 1');end
%% ship geometrical shape is characterized by five point:
% xb1,xb2,xb3,xb4,xb5, their coodinates are described in {b} 
% ship  length: 4.5m   width: 1.5m
L = 2; L1 = L/4; W = 0.8;
xb1=[L/2-L1 -W/2]';xb2=[L/2 0]';xb3=[L/2-L1 W/2]';xb4=[-L/2 W/2]';xb5=[-L/2 -W/2]';
%% rotation matrix from {b} to {n}
Rb_n=[cos(psi) -sin(psi);
      sin(psi) cos(psi)]; 
%% trasfer the five point in {b} to {n}
xn1=Rb_n*xb1+pos;
xn2=Rb_n*xb2+pos;
xn3=Rb_n*xb3+pos;
xn4=Rb_n*xb4+pos;
xn5=Rb_n*xb5+pos;
%% plot
N=pos(1); E=pos(2);
plot([xn1(2) xn2(2)],[xn1(1) xn2(1)],color,[xn2(2) xn3(2)],[xn2(1) xn3(1)],color,...
[xn3(2) xn4(2)],[xn3(1) xn4(1)],color,[xn4(2) xn5(2)],[xn4(1) xn5(1)],color,...
[xn5(2) xn1(2)],[xn5(1) xn1(1)],color,'linewidth',1); % plot ship model
XX = [xn1(1) xn2(1) xn3(1) xn4(1) xn5(1) ]; YY = [xn1(2) xn2(2) xn3(2) xn4(2) xn5(2) ];
fill(YY,XX,color);

end

