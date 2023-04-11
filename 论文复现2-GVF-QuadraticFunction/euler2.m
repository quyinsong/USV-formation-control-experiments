function xnext = euler2( xdot,x,ts )
% EULER  xnext = euler2(xdot,x,ts) Integrate a system of ordinary differential equations using 
%	  Euler's 2nd-order method.
% x(k+1)=x(k)+ts*(f(x(k),u(k)))
% INPUT:
% xdot: f(x(k),u(k))
% x: x(k)
% ts: sample time
% OUTPUT:
% xnext: x(k+1)
% Author:   Quyisnong
% Date:     14th Jan 2022

xnext=x+ts*xdot;

end

