function y = TD1( vd,ts )
%DSC 1-st TD
persistent x
if isempty(x)
    x = 0;
end
T = 1;
x_dot = -(x-vd)/T;
x = euler2(x_dot,x,ts);

y = x;
end

