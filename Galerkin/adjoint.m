%
% u=adjoint(m,u0)
%
% Integrate the heat equation backward in time using the Crank-Nicolson
% implicit Euler method.
%
% The system of equations for each step in time is of the form
%
%    A*u^(n+1)=B*u^(n)
%
% where A and B are tridiagonal matrices.  

function u = adjoint(m,u0)


global A;
global B;

% Do the time steps.

u=u0;
R = B'*inv(A');
%R = (A')\B';
for i=1:m
  u = R*u;
  u(1) = 0 ; u(end) = 0;
end
