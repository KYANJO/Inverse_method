%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ufinal=forward(m,n,deltat,deltax,D,u0)
%
% Integrate the heat equation forward in time using the Crank-Nicolson
% implicit Euler method.
%
% The system of equations for each step in time is of the form
%
%    A*u^(k+1)=B*u^(k)
%
% where A and B are tridiagonal matrices.  
%
% The input parameters are:
%    m         Number of time steps.
%    u0        Initial solution 
%
% Outputs are:
%    u         The solution at time T=m*deltat
%    A,B       The Crank-Nicolson sparse matrices.

function u = forward(m,u0)


global A;
global B;

%
% Do the time steps.
%
u=u0;
R = A\B;
for i=1:m
  u = R*u;
  %u(1) = 0 ; u(end) = 0;
end
